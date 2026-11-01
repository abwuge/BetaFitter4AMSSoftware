// Extract primary reference records after measured-only MC models are frozen.
#include <limits>
#include <vector>
#include <algorithm>
#define main native_producer_main
#include "../../amsd69n/analysis_amsd69n.C"
#undef main


#include <map>
#include <memory>
class ParentReferenceReader {
 std::unique_ptr<TFile> file;
 TTree *tree;
 AMSEventR *event;
 TBranch *header;
 std::map<std::pair<unsigned int,unsigned int>,Long64_t> index;
public:
 ParentReferenceReader(const char *path):file(TFile::Open(path,"READ")),tree(0),event(0),header(0) {
  if(!file || file->IsZombie())return;
  tree=dynamic_cast<TTree*>(file->Get("AMSRoot"));if(!tree)return;
  event=new AMSEventR;
  tree->SetBranchStatus("*",0);
  const char *names[]={"ev.fHeader","ev.fMCEventg","ev.fTofMCCluster","ev.fTrMCCluster","ev.fBetaH","ev.fTrTrack"};
  for(int j=0;j<6;++j){UInt_t found=0;tree->SetBranchStatus(names[j],1,&found);if(found!=1){tree=0;return;}}
  if(tree->SetBranchAddress("ev.",&event)<0){tree=0;return;}
  header=tree->GetBranch("ev.fHeader");if(!header){tree=0;return;}
  for(Long64_t entry=0;entry<tree->GetEntries();++entry){
   if(header->GetEntry(entry)<=0){std::cerr<<"HEADER_FAILURE "<<entry<<std::endl;tree=0;return;}
   std::pair<unsigned int,unsigned int> key(event->Run(),event->Event());
   if(index.count(key)){std::cerr<<"DUPLICATE_IDENTITY\n";tree=0;return;}
   index[key]=entry;
  }
 }
 bool valid()const{return tree!=0;}
 AMSEventR *GetEventFast(unsigned int run,unsigned int number){
  auto it=index.find(std::make_pair(run,number));if(it==index.end())return 0;
  if(tree->GetEntry(it->second)<=0)return 0;
  return event;
 }
};

int main(int argc,char **argv) {
    if(argc!=4){std::cerr<<"usage: mc_reference raw_file_or_list targets output.root\n";return 2;}
    ParentReferenceReader chain(argv[1]);if(!chain.valid())return 3;
    TFile output(argv[3],"CREATE");if(output.IsZombie())return 4;
    TTree tree("reference","Native primary MC acceptance and TOF-cluster identities");
    unsigned int run,event;int found,requested_track,requested_beta;
    int gen_particle,gen_track,gen_parent,gen_nskip,nmc,track_exists,beta_exists,beta_track;
    float gen_mass,gen_charge,gen_momentum;
    float gen_pos[3],gen_dir[3];
    int bar[4];
    std::vector<int> plane,particle,track,parent,nskip;
    std::vector<float> momentum,mass,charge,position,direction;
    std::vector<int> tof_layer,tof_bar,tof_particle,tof_track,tof_parent;
    std::vector<float> tof_time,tof_beta,tof_edep,tof_step,tof_position;
    std::vector<int> tk_layer,tk_particle,tk_track,tk_primary;
    std::vector<float> tk_momentum,tk_position;
#define S(NAME,TYPE) tree.Branch(#NAME,&NAME,#NAME TYPE)
#define V(NAME) tree.Branch(#NAME,&NAME)
    S(run,"/i");S(event,"/i");S(found,"/I");S(requested_track,"/I");S(requested_beta,"/I");
    S(gen_particle,"/I");S(gen_track,"/I");S(gen_parent,"/I");S(gen_nskip,"/I");S(nmc,"/I");
    S(track_exists,"/I");S(beta_exists,"/I");S(beta_track,"/I");
    S(gen_mass,"/F");S(gen_charge,"/F");S(gen_momentum,"/F");
    tree.Branch("gen_pos",gen_pos,"gen_pos[3]/F");tree.Branch("gen_dir",gen_dir,"gen_dir[3]/F");
    tree.Branch("measured_bar",bar,"measured_bar[4]/I");
    V(plane);V(particle);V(track);V(parent);V(nskip);V(momentum);V(mass);V(charge);V(position);V(direction);
    V(tof_layer);V(tof_bar);V(tof_particle);V(tof_track);V(tof_parent);V(tof_time);V(tof_beta);V(tof_edep);V(tof_step);V(tof_position);
    V(tk_layer);V(tk_particle);V(tk_track);V(tk_primary);V(tk_momentum);V(tk_position);
#undef S
#undef V
    std::ifstream targets(argv[2]);if(!targets)return 5;
    int total=0,missing=0;const float nan=std::numeric_limits<float>::quiet_NaN();
    while(targets>>run>>event>>requested_track>>requested_beta>>bar[0]>>bar[1]>>bar[2]>>bar[3]) {
        found=0;nmc=0;track_exists=beta_exists=0;beta_track=-1;
        gen_particle=gen_track=gen_parent=gen_nskip=-999;gen_mass=gen_charge=gen_momentum=nan;
        std::fill(gen_pos,gen_pos+3,nan);std::fill(gen_dir,gen_dir+3,nan);
        plane.clear();particle.clear();track.clear();parent.clear();nskip.clear();momentum.clear();mass.clear();charge.clear();position.clear();direction.clear();
        tof_layer.clear();tof_bar.clear();tof_particle.clear();tof_track.clear();tof_parent.clear();tof_time.clear();tof_beta.clear();tof_edep.clear();tof_step.clear();tof_position.clear();
        tk_layer.clear();tk_particle.clear();tk_track.clear();tk_primary.clear();tk_momentum.clear();tk_position.clear();
        AMSEventR *ev=chain.GetEventFast(run,event);found=ev && ev->Run()==run && ev->Event()==event;
        if(found) {
            track_exists=ev->pTrTrack(requested_track)!=0;
            BetaHR *betah=ev->pBetaH(requested_beta);beta_exists=betah!=0;if(betah)beta_track=betah->iTrTrack();
            nmc=ev->nMCEventg();MCEventgR *gen=nmc?ev->pMCEventg(0):0;
            if(gen) {
                gen_particle=gen->Particle;gen_track=gen->trkID;gen_parent=gen->parentID;gen_nskip=gen->Nskip;
                gen_mass=gen->Mass;gen_charge=gen->Charge;gen_momentum=gen->Momentum;
                for(int j=0;j<3;++j){gen_pos[j]=gen->Coo[j];gen_dir[j]=gen->Dir[j];}
            }
            for(int j=1;j<nmc;++j) {
                MCEventgR *v=ev->pMCEventg(j);if(!v || v->Nskip < -1020 || v->Nskip > -1000)continue;
                plane.push_back(-v->Nskip-1000);particle.push_back(v->Particle);track.push_back(v->trkID);parent.push_back(v->parentID);nskip.push_back(v->Nskip);
                momentum.push_back(v->Momentum);mass.push_back(v->Mass);charge.push_back(v->Charge);
                for(int k=0;k<3;++k){position.push_back(v->Coo[k]);direction.push_back(v->Dir[k]);}
            }
            for(int j=0;j<ev->nTofMCCluster();++j) {
                TofMCClusterR *v=ev->pTofMCCluster(j);if(!v)continue;
                int l=v->GetLayer();if(l<0||l>3||v->GetBar()!=bar[l])continue;
                tof_layer.push_back(l);tof_bar.push_back(v->GetBar());tof_particle.push_back(v->Particle);tof_track.push_back(v->GtrkID);tof_parent.push_back(v->ParentNo);
                tof_time.push_back(v->TOF);tof_beta.push_back(v->Beta);tof_edep.push_back(v->EdepR);tof_step.push_back(v->Step);
                for(int k=0;k<3;++k)tof_position.push_back(v->Coo[k]);
            }
            for(int j=0;j<ev->nTrMCCluster();++j) {
                TrMCClusterR *v=ev->pTrMCCluster(j);if(!v)continue;int il=std::abs(v->GetTkId()/100);
                if(il==8)il=1;else if(il<8)il=il+1;--il;if(il!=4)continue;
                tk_layer.push_back(il);tk_particle.push_back(v->GetPart());tk_track.push_back(v->GetGtrkID());tk_primary.push_back(v->IsPrimary()?1:0);tk_momentum.push_back(v->GetMomentum());
                for(int k=0;k<3;++k)tk_position.push_back(nan); // No Tracker geometry initialization in this reader.
            }
        } else ++missing;
        tree.Fill();if(++total%250==0)std::cout<<"REFERENCE "<<total<<" missing "<<missing<<std::endl;
    }
    output.cd();tree.Write();output.Close();std::cout<<"DONE "<<total<<" missing "<<missing<<std::endl;
    return missing?6:0;
}
