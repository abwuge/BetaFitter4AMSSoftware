#include <TChain.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
constexpr double U=0.9314941;
struct Acc{long n=0; long double den=0; std::array<long double,8> num{};};
double massOf(int p,float z){return p==69?14.899169:2.0*std::abs(z)*U;}
double kin(double p,double m){return std::sqrt(p*p+m*m)-m;}
double sAtZ(double z,const double zl[4],const double sl[4]){
  int seg=0;if(z<=zl[1]&&z>=zl[2])seg=1;else if(z<zl[2])seg=2;else if(z>zl[0])seg=0;else if(z<zl[3])seg=2;
  double dz=zl[seg]-zl[seg+1]; return sl[seg]+(zl[seg]-z)/dz*(sl[seg+1]-sl[seg]);
}
// Uniform loss over [sa,sb], contribution to arrival hit at sh.
double kernel(double sa,double sb,double sh){if(sb<sa)std::swap(sa,sb);if(sh<=sa)return 0;double u2=std::min(sb,sh),du=u2-sa;if(du<=0)return 0;return (sh*du-0.5*(u2*u2-sa*sa))/(sb-sa);}
int cat(int a,int b){
 if(a==4&&b==3)return 0; if(a==3&&b==6)return 1; if(a==6&&b==5)return 2; if(a==5&&b==7)return 3;
 if((a==7&&b==8)||(a==8&&b==9)||(a==9&&b==10)||(a==10&&b==11)||(a==11&&b==12)||(a==12&&b==13)||(a==13&&b==15))return 4;
 if(a==15&&b==14)return 5; if(a==14&&b==16)return 6; if(a==16&&b==17)return 7; return -1;
}
int main(int argc,char**argv){std::string list=argc>1?argv[1]:"input_Z2.list";int nfmax=argc>2?std::stoi(argv[2]):100;
 std::ifstream in(list);std::string f;TChain ch("amstreea");int nf=0;while(nf<nfmax&&std::getline(in,f)){if(f.empty())continue;ch.Add(f.c_str());nf++;}
 int mpar,qs;float mch,pm[21],coo[21][3],ed[4],leng[4],pos[4][3],tl[4];ch.SetBranchStatus("*",0);for(auto n:{"mpar","mch","mevmom1","mevcoo1","tof_edep","tof_leng","tof_pos","tof_tl","tof_qs"})ch.SetBranchStatus(n,1);
 ch.SetBranchAddress("mpar",&mpar);ch.SetBranchAddress("mch",&mch);ch.SetBranchAddress("mevmom1",pm);ch.SetBranchAddress("mevcoo1",coo);ch.SetBranchAddress("tof_edep",ed);ch.SetBranchAddress("tof_leng",leng);ch.SetBranchAddress("tof_pos",pos);ch.SetBranchAddress("tof_tl",tl);ch.SetBranchAddress("tof_qs",&qs);
 const int order[]={4,3,6,5,7,8,9,10,11,12,13,15,14,16,17};std::array<Acc,10> ac;long raw=ch.GetEntries(),good=0,badloss=0;
 for(long ie=0;ie<raw;ie++){ch.GetEntry(ie);double m=massOf(mpar,mch);if(!(m>0&&pm[4]>0))continue;double beta=pm[4]/std::sqrt(pm[4]*pm[4]+m*m);if(!(beta>0.2&&beta<0.9))continue;bool ok=true;for(int i=0;i<4;i++)if(!(ed[i]>0&&tl[i]!=-1&&std::isfinite(leng[i])&&std::isfinite(pos[i][2])))ok=false;for(int id:order)if(!(pm[id]>0&&std::isfinite(coo[id][2])))ok=false;if(!ok)continue;
 double sl[4]={0,leng[1]-leng[0],leng[2]-leng[0],leng[3]-leng[0]},zl[4]={pos[0][2],pos[1][2],pos[2][2],pos[3][2]};if(!(sl[1]>0&&sl[2]>sl[1]&&sl[3]>sl[2]&&zl[0]>zl[1]&&zl[1]>zl[2]&&zl[2]>zl[3]))continue;
 std::array<double,4> G{};for(int ih=0;ih<4;ih++)for(int j=0;j<3;j++)if(sl[ih]>sl[j])G[ih]+=ed[j]*1e-3*(sl[ih]-sl[j]);
 std::array<std::array<double,4>,8> H{};bool lossok=true;for(size_t k=0;k+1<sizeof(order)/sizeof(order[0]);k++){int a=order[k],b=order[k+1],c=cat(a,b);double d=kin(pm[a],m)-kin(pm[b],m);if(c<0||d<-1e-5){lossok=false;break;}if(d<0)d=0;double sa=sAtZ(coo[a][2],zl,sl),sb=sAtZ(coo[b][2],zl,sl);for(int ih=0;ih<4;ih++)H[c][ih]+=d*kernel(sa,sb,sl[ih]);}if(!lossok){badloss++;continue;}
 double gm=0;for(double x:G)gm+=x;gm/=4;for(double &x:G)x-=gm;std::array<double,8> dot{};double den=0;for(double x:G)den+=x*x;if(!(den>0))continue;
 for(int c=0;c<8;c++){double hm=0;for(double x:H[c])hm+=x;hm/=4;for(int i=0;i<4;i++)dot[c]+=G[i]*(H[c][i]-hm);}
 double E=std::sqrt(pm[4]*pm[4]+m*m),kap=m*m/std::pow(E*E-m*m,1.5),w=kap*kap;good++;
 for(int mode=0;mode<2;mode++){if(mode==1&&qs!=1111)continue;for(int ib=0;ib<5;ib++){double bm=0.5+0.1*ib;if(beta>=bm)continue;int ix=mode*5+ib;ac[ix].n++;ac[ix].den+=w*den;for(int c=0;c<8;c++)ac[ix].num[c]+=w*dot[c];}}
 }
 const char* names[]={"S1_bracket","gap_S1S2","S2_bracket","S2_to_L2","inner_tracker_magnet","S3_bracket","gap_S3S4","S4_bracket"};
 std::cout<<std::setprecision(8)<<"list="<<list<<" files="<<nf<<" raw="<<raw<<" good="<<good<<" badloss="<<badloss<<"\n";
 for(int mode=0;mode<2;mode++){std::cout<<(mode?"[QS1111]":"[ALL]")<<"\n";for(int ib=0;ib<5;ib++){auto&a=ac[mode*5+ib];std::cout<<"beta<"<<0.5+0.1*ib<<" n="<<a.n;long double sum=0;for(int c=0;c<8;c++){long double v=a.den?a.num[c]/a.den:0;sum+=v;std::cout<<" "<<names[c]<<"="<<double(v);}std::cout<<" total="<<double(sum)<<"\n";}}
}
