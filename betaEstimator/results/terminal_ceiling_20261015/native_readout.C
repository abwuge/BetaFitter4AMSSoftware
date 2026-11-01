// Minimal native reconstructed-counter plane/timing verification.
#include <map>
#include <memory>
#include <limits>
#define main native_producer_main
#include "/nas01/home/hxwu/AMSSoft/amsd69n/analysis_amsd69n.C"
#undef main
int main(int argc,char **argv){
 if(argc!=4)return 2;
 std::unique_ptr<TFile> input(TFile::Open(argv[1],"READ"));if(!input||input->IsZombie())return 3;
 TTree *src=dynamic_cast<TTree*>(input->Get("AMSRoot"));if(!src)return 4;
 AMSEventR *ev=new AMSEventR;src->SetBranchStatus("*",0);
 const char *names[]={"ev.fHeader","ev.fTofClusterH","ev.fBetaH"};
 for(int k=0;k<3;++k){UInt_t found=0;src->SetBranchStatus(names[k],1,&found);if(found!=1)return 5;}
 if(src->SetBranchAddress("ev.",&ev)<0)return 6;
 TBranch *header=src->GetBranch("ev.fHeader");std::map<std::pair<unsigned int,unsigned int>,Long64_t> lookup;
 for(Long64_t i=0;i<src->GetEntries();++i){header->GetEntry(i);if(!lookup.emplace(std::make_pair(ev->Run(),ev->Event()),i).second)return 7;}
 TFile output(argv[3],"CREATE");if(output.IsZombie())return 8;
 unsigned int run,event;int requested_track,requested_beta,bar[4],found,associated;
 float position[4][3],counter_z[4],native_time[4],raw_time[4][2],edep[4];
 TTree tree("readout","Native reconstructed counter planes, corrected times and raw side times");
 tree.Branch("run",&run,"run/i");tree.Branch("event",&event,"event/i");tree.Branch("found",&found,"found/I");tree.Branch("associated",&associated,"associated/I");
 tree.Branch("position",position,"position[4][3]/F");tree.Branch("counter_z",counter_z,"counter_z[4]/F");tree.Branch("native_time",native_time,"native_time[4]/F");tree.Branch("raw_time",raw_time,"raw_time[4][2]/F");tree.Branch("edep",edep,"edep[4]/F");
 std::ifstream targets(argv[2]);if(!targets)return 9;int count=0,missing=0;float nan=std::numeric_limits<float>::quiet_NaN();
 while(targets>>run>>event>>requested_track>>requested_beta>>bar[0]>>bar[1]>>bar[2]>>bar[3]){
  found=associated=0;std::fill(&position[0][0],&position[0][0]+12,nan);std::fill(counter_z,counter_z+4,nan);std::fill(native_time,native_time+4,nan);std::fill(&raw_time[0][0],&raw_time[0][0]+8,nan);std::fill(edep,edep+4,nan);
  auto it=lookup.find(std::make_pair(run,event));
  if(it!=lookup.end()&&src->GetEntry(it->second)>0){
   found=ev->Run()==run&&ev->Event()==event;BetaHR *beta=ev->pBetaH(requested_beta);associated=beta&&beta->iTrTrack()==requested_track;
   if(beta)for(int l=0;l<4;++l){TofClusterHR *cl=0;for(int j=0;j<beta->NTofClusterH();++j){TofClusterHR *candidate=ev->pTofClusterH(beta->iTofClusterH(j));if(candidate&&candidate->Layer==l&&candidate->Bar==bar[l])cl=candidate;}native_time[l]=beta->gTofBetaPar().Time[l];counter_z[l]=TOFGeom::Sci_pz[l][bar[l]];if(cl){for(int j=0;j<3;++j)position[l][j]=cl->Coo[j];for(int j=0;j<2;++j)raw_time[l][j]=cl->Rtime[j];edep[l]=cl->GetEdep();}}
  }
  if(!found)++missing;tree.Fill();++count;
 }
 output.cd();tree.Write();output.Close();std::cout<<"READOUT_COMPLETE "<<count<<" missing "<<missing<<std::endl;return missing?10:0;
}
