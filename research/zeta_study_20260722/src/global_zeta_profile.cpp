#include <TChain.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

constexpr double C = 0.299792458; // m/ns
constexpr double U = 0.9314941;
struct Event {
  double beta, mass, e0;
  std::array<double,3> dl;
  std::array<double,4> edep;
  std::array<double,4> meas;
  std::array<double,4> truth;
  int qs;
  double active12_ratio, passive23_ratio, approx_zeta;
};

double invBetaP(double p,double m){ return std::sqrt(1.0+(m/p)*(m/p)); }
double kinetic(double p,double m){ return std::sqrt(p*p+m*m)-m; }
struct Pt{double z,p;};

double pAt(const std::vector<Pt>& v,double z){
  if(z>v.front().z || z<v.back().z) return std::numeric_limits<double>::quiet_NaN();
  for(size_t i=0;i+1<v.size();++i){
    if(v[i].z+1e-9>=z && z>=v[i+1].z-1e-9){
      double dz=v[i].z-v[i+1].z;
      if(std::abs(dz)<1e-10) return 0.5*(v[i].p+v[i+1].p);
      double a=(v[i].z-z)/dz;
      return v[i].p*(1-a)+v[i+1].p*a;
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

double avgInvBeta(const std::vector<Pt>& v,double za,double zb,double m){
  if(za<zb) std::swap(za,zb);
  if(za>v.front().z || zb<v.back().z || za<=zb) return std::numeric_limits<double>::quiet_NaN();
  std::vector<double> zs{za,zb};
  for(auto &q:v) if(q.z<za && q.z>zb) zs.push_back(q.z);
  std::sort(zs.begin(),zs.end(),std::greater<double>());
  zs.erase(std::unique(zs.begin(),zs.end(),[](double a,double b){return std::abs(a-b)<1e-6;}),zs.end());
  double integ=0;
  for(size_t i=0;i+1<zs.size();++i){
    double z1=zs[i],z2=zs[i+1],zm=0.5*(z1+z2);
    double p1=pAt(v,z1),p2=pAt(v,z2),pm=pAt(v,zm);
    if(!(p1>0&&p2>0&&pm>0)) return std::numeric_limits<double>::quiet_NaN();
    integ+=(z1-z2)*(invBetaP(p1,m)+4*invBetaP(pm,m)+invBetaP(p2,m))/6.0;
  }
  return integ/(za-zb);
}

double centeredSSE(const std::array<double,4>& a,const std::array<double,4>& b){
  double r[4],mean=0; for(int i=0;i<4;i++){r[i]=a[i]-b[i];mean+=r[i];} mean/=4;
  double s=0; for(double x:r)s+=(x-mean)*(x-mean); return s;
}

bool model(const Event&e,double z,std::array<double,4>&t){
  t={0,0,0,0}; double en=e.e0;
  for(int i=1;i<4;i++){
    en-=z*e.edep[i-1];
    if(!(en>e.mass) || !std::isfinite(en)) return false;
    double ib=1.0/std::sqrt(1.0-(e.mass/en)*(e.mass/en));
    t[i]=t[i-1]+e.dl[i-1]/100.0/C*ib;
  }
  return true;
}

double objective(const std::vector<Event>&ev,double z,bool truth,double bmax,int qmode,long &used){
  long n=0; long double s=0; std::array<double,4> t;
  for(auto&e:ev){
    if(e.beta>=bmax) continue;
    if(qmode==1 && e.qs!=1111) continue;
    if(qmode==2 && e.qs==1111) continue;
    if(!model(e,z,t)) continue;
    s+=centeredSSE(t,truth?e.truth:e.meas); n++;
  }
  used=n; return n?double(s/n):std::numeric_limits<double>::infinity();
}

double golden(const std::function<double(double)>&f,double a,double b){
  const double gr=(std::sqrt(5.0)-1)/2; double c=b-gr*(b-a),d=a+gr*(b-a),fc=f(c),fd=f(d);
  for(int i=0;i<55;i++){
    if(fc<fd){b=d;d=c;fd=fc;c=b-gr*(b-a);fc=f(c);} else {a=c;c=d;fc=fd;d=a+gr*(b-a);fd=f(d);}
  }
  return 0.5*(a+b);
}

void qstats(std::vector<double> v,const std::string& name){
  v.erase(std::remove_if(v.begin(),v.end(),[](double x){return !std::isfinite(x);}),v.end());
  std::sort(v.begin(),v.end()); if(v.empty()){std::cout<<name<<" n=0\n";return;}
  auto q=[&](double p){double u=p*(v.size()-1);size_t i=(size_t)u;double a=u-i;return v[i]*(1-a)+v[std::min(i+1,v.size()-1)]*a;};
  long double sm=0;for(double x:v)sm+=x;
  std::cout<<name<<" n="<<v.size()<<" mean="<<double(sm/v.size())<<" q05="<<q(.05)<<" q16="<<q(.16)<<" median="<<q(.5)<<" q84="<<q(.84)<<" q95="<<q(.95)<<"\n";
}

int main(int argc,char**argv){
  std::string listfile=argc>1?argv[1]:"input_Z2.list";
  int nfiles=argc>2?std::stoi(argv[2]):10;
  std::ifstream in(listfile); std::string fn; TChain ch("amstreea"); int nf=0;
  while(nf<nfiles && std::getline(in,fn)){if(fn.empty())continue;ch.Add(fn.c_str());nf++;}
  int mpar=0,tof_qs=0; float mch=0,pm[21],coo[21][3],ed[4],leng[4],tl[4],pos[4][3];
  ch.SetBranchStatus("*",0);
  for(auto n:{"mpar","mch","mevmom1","mevcoo1","tof_edep","tof_leng","tof_tl","tof_pos","tof_qs"})ch.SetBranchStatus(n,1);
  ch.SetBranchAddress("mpar",&mpar);ch.SetBranchAddress("mch",&mch);ch.SetBranchAddress("mevmom1",pm);ch.SetBranchAddress("mevcoo1",coo);
  ch.SetBranchAddress("tof_edep",ed);ch.SetBranchAddress("tof_leng",leng);ch.SetBranchAddress("tof_tl",tl);ch.SetBranchAddress("tof_pos",pos);ch.SetBranchAddress("tof_qs",&tof_qs);
  std::vector<Event> ev; ev.reserve(ch.GetEntries()/4); long nall=ch.GetEntries(),n4=0,nphys=0,ntruth=0;
  for(long ie=0;ie<nall;ie++){
    ch.GetEntry(ie); double mass=(mpar==69?14.899169:2.0*mch*U); if(!(mass>0&&pm[4]>0))continue;
    double beta=pm[4]/std::sqrt(pm[4]*pm[4]+mass*mass); if(!(beta>0.2&&beta<0.9))continue;
    bool ok=true;for(int i=0;i<4;i++)if(!(tl[i]!=-1&&std::isfinite(tl[i])&&ed[i]>0&&std::isfinite(ed[i])&&std::isfinite(leng[i])&&std::isfinite(pos[i][2])))ok=false;
    if(!ok)continue;n4++;
    std::array<double,3> dl{{leng[1]-leng[0],leng[2]-leng[1],leng[3]-leng[2]}};for(double x:dl)if(!(x>0&&x<300))ok=false;if(!ok)continue;
    std::vector<Pt> pts;for(int i=0;i<21;i++)if(pm[i]>0&&std::isfinite(pm[i])&&std::isfinite(coo[i][2]))pts.push_back({coo[i][2],pm[i]});
    std::sort(pts.begin(),pts.end(),[](auto&a,auto&b){return a.z>b.z;});
    std::vector<Pt> uniq;for(auto&q:pts){if(uniq.empty()||std::abs(q.z-uniq.back().z)>1e-5)uniq.push_back(q);else uniq.back().p=0.5*(uniq.back().p+q.p);}pts.swap(uniq);
    std::array<double,4> truth{{0,0,0,0}};for(int i=1;i<4;i++){double ai=avgInvBeta(pts,pos[i-1][2],pos[i][2],mass);if(!std::isfinite(ai)){ok=false;break;}truth[i]=truth[i-1]+dl[i-1]/100.0/C*ai;}if(!ok)continue;ntruth++;
    Event e; e.beta=beta;e.mass=mass;e.e0=std::sqrt(pm[4]*pm[4]+mass*mass);e.dl=dl;e.qs=tof_qs;e.truth=truth;
    for(int i=0;i<4;i++){e.edep[i]=ed[i]*1e-3;e.meas[i]=tl[i];}
    // True active losses bracketed by checkpoint pairs: S1 4->3, S2 6->5.
    double k4=kinetic(pm[4],mass),k3=kinetic(pm[3],mass),k6=kinetic(pm[6],mass),k5=kinetic(pm[5],mass),k15=kinetic(pm[15],mass);
    double den=e.edep[0]+e.edep[1];
    e.active12_ratio=den>0?((k4-k3)+(k6-k5))/den:NAN;
    e.passive23_ratio=den>0?(k5-k15)/den:NAN;
    double pS2=pAt(pts,pos[1][2]),pS3=pAt(pts,pos[2][2]);
    double kS2=kinetic(pS2,mass),kS3=kinetic(pS3,mass);
    e.approx_zeta=den>0?(k4-0.5*(kS2+kS3))/den:NAN;
    // Stable model up to zeta=3; reject cases that would stop / produce NaN in the region of interest.
    if(e.e0-3.0*(e.edep[0]+e.edep[1]+e.edep[2])<=e.mass)continue;
    ev.push_back(e);nphys++;
  }
  std::cout<<std::setprecision(8)<<"FILES="<<nf<<" raw="<<nall<<" fourhit_lowbeta="<<n4<<" truth_integrable="<<ntruth<<" stable_events="<<ev.size()<<" qs1111="<<std::count_if(ev.begin(),ev.end(),[](auto&e){return e.qs==1111;})<<"\n";
  for(int qmode=0;qmode<=1;qmode++){
    std::cout<<(qmode?"[QS1111]":"[ALL]")<<"\n";
    for(double bm:{0.5,0.6,0.7,0.8,0.9}){
      for(bool tr:{true,false}){
        long used=0;auto ff=[&](double z){long n;return objective(ev,z,tr,bm,qmode,n);};double z=golden(ff,0.0,3.0);double fmin=ff(z);objective(ev,z,tr,bm,qmode,used);
        std::cout<<(tr?"truth":"measured")<<" beta<"<<bm<<" zeta="<<z<<" mse="<<fmin<<" used="<<used<<" f0="<<ff(0)<<" f1="<<ff(1)<<" f1.9="<<ff(1.9)<<" f2="<<ff(2)<<"\n";
      }
    }
  }
  for(int qmode=0;qmode<=1;qmode++){
    std::vector<double>a,p,z;for(auto&e:ev)if(e.beta<0.7&&(!qmode||e.qs==1111)){a.push_back(e.active12_ratio);p.push_back(e.passive23_ratio);z.push_back(e.approx_zeta);}
    std::string pre=qmode?"qs1111_":"all_";qstats(a,pre+"active12_true_over_reco");qstats(p,pre+"central_passive_over_reco12");qstats(z,pre+"approx_timing_zeta");
  }
}
