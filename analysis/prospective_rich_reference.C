// Read the existing flight RICH association and error after prediction freeze.
// No ring selection by TOF agreement, isotope, or fitted mass is performed.
#include <limits>
#define main native_producer_main
#include "../../amsd69n/analysis_amsd69n.C"
#undef main

int main(int argc,char **argv) {
    if(argc!=4){std::cerr<<"usage: rich_reference raw_files targets output.root\n";return 2;}
    AMSChain chain;std::ifstream input(argv[1]);std::string path;int files=0;
    while(std::getline(input,path)){if(path.empty())continue;if(chain.Add(path.c_str())<=0)return 3;++files;}
    if(!files)return 3;
    AMSAnalysis analysis;analysis.LoadHeader();
    TFile output(argv[3],"CREATE");if(output.IsZombie())return 4;
    TTree tree("reference","Existing native RICH ring and Tracker association");
    unsigned int run=0,event=0;
    int itrtrack=-1,irich=-1,rich_itrtrack=-1,found=0,ring_found=0;
    float beta[3],error;const float nan=std::numeric_limits<float>::quiet_NaN();
    tree.Branch("run",&run,"run/i");tree.Branch("event",&event,"event/i");
    tree.Branch("itrtrack",&itrtrack,"itrtrack/I");tree.Branch("irich",&irich,"irich/I");
    tree.Branch("rich_itrtrack",&rich_itrtrack,"rich_itrtrack/I");tree.Branch("found",&found,"found/I");
    tree.Branch("ring_found",&ring_found,"ring_found/I");tree.Branch("beta",beta,"beta[3]/F");tree.Branch("error",&error,"error/F");
    std::ifstream targets(argv[2]);if(!targets)return 5;
    int total=0,missing=0;
    while(targets>>run>>event>>itrtrack>>irich) {
        found=ring_found=0;rich_itrtrack=-1;beta[0]=beta[1]=beta[2]=error=nan;
        AMSEventR *ev=chain.GetEventFast(run,event);found=ev&&ev->Run()==run&&ev->Event()==event;
        if(found) {
            ev->RecordRTIRun();analysis.unbiasetag0=1;analysis.unbiasetag=0;analysis.istkaddhit=0;
            analysis.InitEvent(ev);analysis.isreal=true;analysis.LoadHeaderEv();
            RichRingR *ring=ev->pRichRing(irich);
            if(ring){ring_found=1;rich_itrtrack=ring->iTrTrack();beta[0]=ring->getBeta();beta[1]=ring->Beta;beta[2]=ring->BetaRefit;error=ring->getBetaError();}
        }else ++missing;
        tree.Fill();if(++total%100==0)std::cout<<"REFERENCE "<<total<<" missing "<<missing<<std::endl;
    }
    output.cd();tree.Write();output.Close();std::cout<<"DONE "<<total<<" missing "<<missing<<std::endl;
    return missing?6:0;
}
