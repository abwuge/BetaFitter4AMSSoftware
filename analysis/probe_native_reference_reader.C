// Read selected parent-event branches without detector reconstruction.
#define main native_producer_main
#include "../../amsd69n/analysis_amsd69n.C"
#undef main
#include <TBranchElement.h>
int main(int argc,char **argv){
 if(argc!=2)return 2;
 TFile *f=TFile::Open(argv[1],"READ");if(!f||f->IsZombie())return 3;
 TTree *t=dynamic_cast<TTree*>(f->Get("AMSRoot"));if(!t)return 4;
 AMSEventR *event=new AMSEventR;
 t->SetBranchStatus("*",0);
 const char* names[]={"ev.fHeader","ev.fMCEventg","ev.fTofMCCluster","ev.fTrMCCluster","ev.fBetaH","ev.fTrTrack"};
 for(int j=0;j<6;++j){UInt_t found=0;t->SetBranchStatus(names[j],1,&found);std::cout<<"BRANCH "<<names[j]<<" found "<<found<<std::endl;}
 std::cout<<"PARENT_BIND "<<t->SetBranchAddress("ev.",&event)<<std::endl;
 for(int i=0;i<3;++i){
  Long64_t bytes=t->GetEntry(i);std::cout<<"ENTRY "<<i<<" bytes "<<bytes<<" pointer "<<event<<" id "<<event->Run()<<"/"<<event->Event()<<" nMC "<<event->nMCEventg()<<std::endl;
  MCEventgR *g=event->nMCEventg()?event->pMCEventg(0):0;
  if(g)std::cout<<"GEN "<<g->Particle<<" "<<g->Momentum<<" "<<g->Mass<<" "<<g->trkID<<std::endl;
 }
 return 0;
}
