#include <TChain.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
constexpr double U=.9314941;
double massOf(int p,float z){return p==69?14.899169:2*std::abs(z)*U;} double K(double p,double m){return std::sqrt(p*p+m*m)-m;}
struct S{std::vector<double>v;std::array<long,300>h{};};
void out(std::string n,S&s){std::sort(s.v.begin(),s.v.end());auto q=[&](double p){double u=p*(s.v.size()-1);size_t i=u;double a=u-i;return s.v[i]*(1-a)+s.v[std::min(i+1,s.v.size()-1)]*a;};int im=std::max_element(s.h.begin(),s.h.end())-s.h.begin();long double sm=0;for(double x:s.v)sm+=x;std::cout<<n<<" n="<<s.v.size()<<" mean="<<double(sm/s.v.size())<<" median="<<q(.5)<<" q16="<<q(.16)<<" q84="<<q(.84)<<" mode_bin="<<(im+.5)*.01<<" mode_count="<<s.h[im]<<"\n";}
int main(){std::ifstream in("input_Z2.list");std::string f;TChain c("amstreea");while(std::getline(in,f))if(!f.empty())c.Add(f.c_str());int mp,qs;float z,p[21],e[4],tl[4];c.SetBranchStatus("*",0);for(auto n:{"mpar","mch","mevmom1","tof_edep","tof_tl","tof_qs"})c.SetBranchStatus(n,1);c.SetBranchAddress("mpar",&mp);c.SetBranchAddress("mch",&z);c.SetBranchAddress("mevmom1",p);c.SetBranchAddress("tof_edep",e);c.SetBranchAddress("tof_tl",tl);c.SetBranchAddress("tof_qs",&qs);S all,un,ov,lowun,lowov;long N=c.GetEntries();for(long i=0;i<N;i++){c.GetEntry(i);double m=massOf(mp,z);if(!(p[4]>0&&p[7]>0&&e[0]>0&&e[1]>0&&tl[0]!=-1&&tl[1]!=-1))continue;double b=p[4]/sqrt(p[4]*p[4]+m*m);double r=(K(p[4],m)-K(p[7],m))/((e[0]+e[1])*1e-3);if(!(r>0&&r<3))continue;bool u=qs/100==11;for(S*x:{&all,u?&un:&ov,(b<.7?(u?&lowun:&lowov):nullptr)})if(x){x->v.push_back(r);int ib=std::min(299,std::max(0,int(r/.01)));x->h[ib]++;}}std::cout<<std::setprecision(7);out("all",all);out("S1S2_unoverlap",un);out("S1S2_overlap",ov);out("beta<.7_unoverlap",lowun);out("beta<.7_overlap",lowov);}
