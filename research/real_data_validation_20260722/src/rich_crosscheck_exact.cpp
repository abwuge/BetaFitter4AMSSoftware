#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
struct D{std::vector<double>v;void add(double x){if(std::isfinite(x)&&fabs(x)<3)v.push_back(x);}void p(std::string n){if(v.size()<30)return;sort(v.begin(),v.end());auto q=[&](double p){double a=p*(v.size()-1);size_t i=a;double f=a-i;return v[i]*(1-f)+v[std::min(i+1,v.size()-1)]*f;};long double s=0,s2=0;for(double x:v){s+=x;s2+=x*x;}double m=s/v.size(),r=sqrt(std::max(0.0,(double)(s2/v.size()-m*m))),med=q(.5);long long t1=0,t2=0;for(double x:v){if(fabs(x-med)>.1)t1++;if(fabs(x-med)>.5)t2++;}std::cout<<n<<" n="<<v.size()<<" mean="<<m<<" rms="<<r<<" med="<<med<<" sig68="<<(q(.84)-q(.16))/2<<" tail.1="<<(double)t1/v.size()<<" tail.5="<<(double)t2/v.size()<<"\n";}};
int main(){std::string p="/lustre02/user/yachen/public/amsd69n_B1130_Data_TofCalibN2/1306976244_0.root";TFile fr(p.c_str());auto raw=(TTree*)fr.Get("amstreea");TFile fb2("/tmp/real_center_global_Z2.root"),fb6("/tmp/real_center_global_Z6.root"),fb8("/tmp/real_center_global_Z8.root");TTree*bt[3]={(TTree*)fb2.Get("betaTree"),(TTree*)fb6.Get("betaTree"),(TTree*)fb8.Get("betaTree")};std::cout<<"raw="<<raw->GetEntries()<<" beta="<<bt[0]->GetEntries()<<"\n";
 float tb,rb[3],rpb;int ru;bool rg,rc,rnaf;raw->SetBranchStatus("*",0);for(auto b:{"tof_betah","rich_beta","rich_pb","rich_used","rich_good","rich_clean","rich_NaF"})raw->SetBranchStatus(b,1);raw->SetBranchAddress("tof_betah",&tb);raw->SetBranchAddress("rich_beta",rb);raw->SetBranchAddress("rich_pb",&rpb);raw->SetBranchAddress("rich_used",&ru);raw->SetBranchAddress("rich_good",&rg);raw->SetBranchAddress("rich_clean",&rc);raw->SetBranchAddress("rich_NaF",&rnaf);
 double bl[3],bn[3],R[3];int Zout[3];for(int k=0;k<3;k++){bt[k]->SetBranchAddress("linearBeta",&bl[k]);bt[k]->SetBranchAddress("nonlinearBeta",&bn[k]);bt[k]->SetBranchAddress("innerRigidity",&R[k]);bt[k]->SetBranchAddress("Z",&Zout[k]);}std::map<std::string,D>L,N;long long bad=0,sel=0;for(Long64_t i=0;i<raw->GetEntries();i++){raw->GetEntry(i);for(int k=0;k<3;k++)bt[k]->GetEntry(i);int Z=Zout[0],k=Z==2?0:Z==6?1:Z==8?2:-1;if(k<0)continue;if(Zout[1]!=Z||Zout[2]!=Z||fabs(tb-bl[k])>1e-6)bad++;double br=rb[1];if(!(rg&&rc&&ru>=4&&rpb>.01&&br>.70&&br<1.01&&bl[k]>.4&&bl[k]<1.1&&bn[k]>0&&R[k]>0))continue;sel++;double dl=1/bl[k]-1/br,dn=1/bn[k]-1/br;std::string rad=rnaf?"NaF":"Agl",bb=br<.8?"b.70-.80":br<.9?"b.80-.90":br<.95?"b.90-.95":"b.95-1.01";for(auto key:{std::string("all"),rad,rad+"_"+bb,rad+"_Z"+std::to_string(Z)}){L[key].add(dl);N[key].add(dn);}}
 std::cout<<"alignBad="<<bad<<" selected="<<sel<<"\n";for(auto&x:L){x.second.p("L_"+x.first);N[x.first].p("N_"+x.first);} }
