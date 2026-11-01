#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TSystem.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Measured flight branches only. No simulation file or truth branch is opened.
void extractMassStudy(const char* manifest="config/data_sources.conf",
                      const char* output="results/mass_study/flight.root",
                      int nTrain=64, int nValidation=32) {
  const std::vector<std::string> names={"run","event","tof_betah","tof_tl","tof_leng","tof_barid","tof_pos","tof_dir","tof_edep","tof_bpatt","tk_edep","tk_pos","tk_q","tk_qin","tk_rigidity","tk_rigidity1","tk_rigidityk","tk_erigidity1","tk_chis1","irich","itrtrack","rich_itrtrack","rich_beta","rich_ebeta1","rich_good","rich_clean","rich_NaF","rich_pb","rich_pmt","rich_used","rich_BetaConsistency","rich_PMTChargeConsistency","rich_npe","rich_pos","rich_tile"};
  struct Source {std::string path; int role; unsigned long long key;};
  auto hash=[](const std::string& s){unsigned long long h=1469598103934665603ULL;for(unsigned char c:s){h^=c;h*=1099511628211ULL;}return h;};
  std::ifstream in(manifest);if(!in)throw std::runtime_error("Missing manifest");
  std::vector<Source> pool[2],sources;std::string line;
  while(std::getline(in,line)){if(line.empty()||line[0]=='#')continue;std::istringstream s(line);std::string role,path;int z;s>>role>>z>>path;if(path.empty())continue;int r=role=="validation";pool[r].push_back({path,r,hash(path)});}
  for(int r=0;r<2;r++){std::sort(pool[r].begin(),pool[r].end(),[](const Source&a,const Source&b){return a.key<b.key;});int lim=r?nValidation:nTrain;if(lim>0&&int(pool[r].size())>lim)pool[r].resize(lim);sources.insert(sources.end(),pool[r].begin(),pool[r].end());}
  std::sort(sources.begin(),sources.end(),[](const Source&a,const Source&b){return a.path<b.path;});
  std::string outpath(output);gSystem->mkdir(gSystem->DirName(output),true);
  TFile out(output,"RECREATE");if(out.IsZombie())throw std::runtime_error("Output open failed");
  TTree tree("flight","Measured low-rigidity flight sample; no MC truth");
  Long64_t run,event,entry;int role,source,z,richOK,naf;double R,u,producerU,richU,richError,features[8],rigidity[63],rigidityK[28],erigidity[21],chis[189];
  tree.Branch("rigidityK",rigidityK,"rigidityK[28]/D");
  tree.Branch("run",&run);tree.Branch("event",&event);tree.Branch("entry",&entry);tree.Branch("role",&role);tree.Branch("source",&source);tree.Branch("Z",&z);tree.Branch("R",&R);tree.Branch("u",&u);tree.Branch("producerU",&producerU);tree.Branch("richU",&richU);tree.Branch("richError",&richError);tree.Branch("richOK",&richOK);tree.Branch("NaF",&naf);tree.Branch("features",features,"features[8]/D");tree.Branch("rigidity",rigidity,"rigidity[63]/D");tree.Branch("erigidity",erigidity,"erigidity[21]/D");tree.Branch("chis",chis,"chis[189]/D");
  std::ofstream audit(outpath+".sources.tsv");audit<<"source\trole\tpath\traw\tkept\n";
  std::set<std::pair<Long64_t,Long64_t>> ids;std::map<Long64_t,int> runroles;Long64_t total=0;int count=0;
  for(const auto& src:sources){
    TFile file(src.path.c_str(),"READ");auto t=dynamic_cast<TTree*>(file.Get("amstreea"));if(!t)throw std::runtime_error("Missing input tree: "+src.path);
    t->SetBranchStatus("*",0);std::map<std::string,TLeaf*> leaf;
    for(const auto& name:names){auto l=t->GetLeaf(name.c_str());if(!l)throw std::runtime_error("Missing measured leaf: "+name);t->SetBranchStatus(name.c_str(),1);leaf[name]=l;}
    auto v=[&](const std::string& name,int k=0){return leaf.at(name)->GetValue(k);};
    role=src.role;source=++count;Long64_t kept=0;
    for(entry=0;entry<t->GetEntries();entry++){
      t->GetEntry(entry);double q=v("tk_qin",2)<2.5?v("tk_q",1):v("tk_qin",2);if(!std::isfinite(q))continue;z=int(q+.5);if(z!=1&&z!=2&&z!=6&&z!=8)continue;
      R=v("tk_rigidity");if(!std::isfinite(R)||R<1||R>12||v("tof_bpatt")!=4444)continue;
      bool good=true;for(int k=1;k<8;k++){good&=std::isfinite(v("tk_edep",k))&&v("tk_edep",k)>0;for(int j=0;j<3;j++)good&=std::isfinite(v("tk_pos",3*k+j));}good&=v("tk_pos",14)>0&&v("tk_pos",17)<0;if(!good)continue;
      double x[4],tt[4],xm=0,tm=0,px=0,py=0;
      for(int k=0;k<4;k++){x[k]=v("tof_leng",k)/29.9792458;tt[k]=v("tof_tl",k);xm+=x[k]/4;tm+=tt[k]/4;px+=v("tof_pos",3*k)/240;py+=v("tof_pos",3*k+1)/240;good&=std::isfinite(x[k])&&std::isfinite(tt[k])&&std::isfinite(v("tof_pos",3*k))&&std::isfinite(v("tof_pos",3*k+1))&&v("tof_edep",k)>0&&std::isfinite(v("tof_edep",k))&&v("tof_barid",k)>=0&&v("tof_barid",k)<10;}
      if(!good)continue;double den=0,num=0;for(int k=0;k<4;k++){x[k]-=xm;tt[k]-=tm;den+=x[k]*x[k];num+=x[k]*tt[k];}if(den<=1e-12)continue;u=num/den;if(!std::isfinite(u)||u<=0)continue;double span=*std::max_element(x,x+4)-*std::min_element(x,x+4);if(span<=0)continue;
      double residual[4];for(int k=0;k<4;k++)residual[k]=tt[k]-u*x[k];
      features[0]=1;features[1]=px;features[2]=py;features[3]=px*px;features[4]=py*py;features[5]=std::log(span/4.0);features[6]=(residual[0]-residual[1])/span;features[7]=(residual[2]-residual[3])/span;
      producerU=1/v("tof_betah");
      for(int k=0;k<63;k++)rigidity[k]=v("tk_rigidity1",k);for(int k=0;k<28;k++)rigidityK[k]=v("tk_rigidityk",k);for(int k=0;k<21;k++)erigidity[k]=v("tk_erigidity1",k);for(int k=0;k<189;k++)chis[k]=v("tk_chis1",k);
      double b=v("rich_beta",0),bc=v("rich_BetaConsistency"),pc=v("rich_PMTChargeConsistency");naf=v("rich_NaF")!=0;richU=std::isfinite(b)&&b!=0?1/b:NAN;richError=v("rich_ebeta1");double xx=v("rich_pos",0),yy=v("rich_pos",1),rad=xx*xx+yy*yy,xy=std::max(std::abs(xx),std::abs(yy));int tile=v("rich_tile");std::set<int>badtiles={3,7,12,20,87,100,108};
      richOK=v("irich")>=0&&v("itrtrack")>=0&&v("rich_itrtrack")==v("itrtrack")&&v("rich_good")&&v("rich_clean")&&v("rich_pb")>=.01&&v("rich_pmt")>=3&&v("rich_used")>=4&&std::isfinite(b)&&b>=(naf?.78:.955)&&b<=1.02&&std::isfinite(richError)&&richError>0&&richError<.005&&std::isfinite(bc)&&bc>=0&&bc<=(naf?.01:.005)&&std::isfinite(pc)&&pc>=0&&pc<=15&&std::isfinite(v("rich_npe",0))&&std::isfinite(v("rich_npe",1))&&std::isfinite(v("rich_npe",2))&&v("rich_npe",2)>0&&v("rich_npe",0)>=.4*v("rich_npe",2)&&v("rich_npe",1)>=(naf?1.:2.)&&std::isfinite(xx)&&std::isfinite(yy)&&std::isfinite(v("rich_pos",2))&&(naf?xy<=17:(xy>=19&&rad<=3500&&badtiles.count(tile)==0));
      run=Long64_t(v("run"));event=Long64_t(v("event"));auto rr=runroles.find(run);if(rr!=runroles.end()&&rr->second!=role)throw std::runtime_error("Run overlap");runroles[run]=role;if(!ids.insert({run,event}).second)throw std::runtime_error("Duplicate selected event");tree.Fill();kept++;
    }
    audit<<source<<'\t'<<role<<'\t'<<src.path<<'\t'<<t->GetEntries()<<'\t'<<kept<<'\n';total+=kept;std::cout<<"source "<<count<<"/"<<sources.size()<<" raw="<<t->GetEntries()<<" kept="<<kept<<" total="<<total<<std::endl;
  }
  out.cd();tree.Write();std::ostringstream provenance;provenance<<"No MC truth accessed. Whole-run train/validation split reused. Hash-selected "<<nTrain<<" train and "<<nValidation<<" validation sources. Selection: measured charge, four TOF layers, tracker hits, 1<R<12 GV. Branches:";for(auto& name:names)provenance<<' '<<name;TNamed metadata("provenance",provenance.str().c_str());metadata.Write();out.Close();std::cout<<"DONE events="<<total<<" output="<<output<<std::endl;
}
