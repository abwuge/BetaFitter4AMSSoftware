#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>
struct D{std::vector<double>v;void add(double x){if(std::isfinite(x)&&fabs(x)<2)v.push_back(x);}void p(std::string n){if(v.size()<20)return;sort(v.begin(),v.end());auto q=[&](double p){double a=p*(v.size()-1);size_t i=a;double f=a-i;return v[i]*(1-f)+v[std::min(i+1,v.size()-1)]*f;};std::cout<<n<<" n="<<v.size()<<" med="<<q(.5)<<" sig68="<<(q(.84)-q(.16))/2<<" q05="<<q(.05)<<" q95="<<q(.95)<<"\n";}};
int main(){TFile f("results/10.root");auto t=(TTree*)f.Get("betaTree");double br,bl,bn,R;int Z;t->SetBranchAddress("betaRigidity",&br);t->SetBranchAddress("linearBeta",&bl);t->SetBranchAddress("nonlinearBeta",&bn);t->SetBranchAddress("innerRigidity",&R);t->SetBranchAddress("Z",&Z);std::map<std::string,D>L,N;for(Long64_t i=0;i<t->GetEntries();i++){t->GetEntry(i);if(!(br>.35&&br<1&&bl>0&&bn>0&&R>0))continue;std::string b=br<.6?"b.35-.6":br<.7?"b.6-.7":br<.8?"b.7-.8":br<.9?"b.8-.9":"b.9-1";std::string z="Z"+std::to_string(Z);L[z+"_"+b].add(1/bl-1/br);N[z+"_"+b].add(1/bn-1/br);}for(auto&x:L){x.second.p("L_"+x.first);N[x.first].p("N_"+x.first);} }
