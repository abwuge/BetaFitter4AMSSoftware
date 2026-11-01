// Export the existing betaFitter material geometry without loading its fitted coefficients.
#include <limits>
#include <vector>
#include <algorithm>
#define main native_producer_main
#include "../../amsd69n/analysis_amsd69n.C"
#undef main
#include "../../betaFitter/src/HybridEnergyLoss.C"

int main(int argc,char **argv) {
    if(argc!=3){std::cerr<<"usage: material_geometry measured_geometry.txt output.root\n";return 2;}
    AMSAnalysis analysis;analysis.LoadHeader();
    std::ifstream input(argv[1]);if(!input)return 3;
    TFile output(argv[2],"CREATE");if(output.IsZombie())return 4;
    TTree tree("geometry","Native material geometry, no learned response coefficients");
    unsigned int run,event;double rigidity;int bar[4],valid;float length[4],pos[9][3],dir[9][3],tracker_z[9];
    float center_z,tof_z[4],boundaries[5];
    std::vector<float> high,low,step_length,abundance;
    std::vector<int> segment,at_high,at_low;
    tree.Branch("run",&run,"run/i");tree.Branch("event",&event,"event/i");tree.Branch("valid",&valid,"valid/I");
    tree.Branch("center_z",&center_z,"center_z/F");tree.Branch("tof_z",tof_z,"tof_z[4]/F");tree.Branch("boundaries",boundaries,"boundaries[5]/F");
    tree.Branch("high",&high);tree.Branch("low",&low);tree.Branch("length",&step_length);tree.Branch("abundance",&abundance);
    tree.Branch("segment",&segment);tree.Branch("at_high",&at_high);tree.Branch("at_low",&at_low);
    int total=0,invalid=0;
    while(input>>run>>event>>rigidity) {
        for(int i=0;i<4;++i)input>>bar[i];for(int i=0;i<4;++i)input>>length[i];
        for(int i=0;i<9;++i)for(int j=0;j<3;++j)input>>pos[i][j];
        for(int i=0;i<9;++i)for(int j=0;j<3;++j)input>>dir[i][j];
        if(!input)return 5;for(int i=0;i<9;++i)tracker_z[i]=pos[i][2];
        HybridMaterialProfile p=HybridEnergyLossCalibration::BuildProfile(pos[4],dir[4],rigidity,bar,length,tracker_z,1.0);
        valid=p.valid;center_z=p.centerZ;
        for(int i=0;i<4;++i)tof_z[i]=p.tofCenterZ[i];for(int i=0;i<5;++i)boundaries[i]=p.boundaryZ[i];
        high.clear();low.clear();step_length.clear();abundance.clear();segment.clear();at_high.clear();at_low.clear();
        for(const auto &s:p.steps){high.push_back(s.zHigh);low.push_back(s.zLow);step_length.push_back(s.length);segment.push_back(s.segment);at_high.push_back(s.tofAtHigh);at_low.push_back(s.tofAtLow);for(double v:s.abundance)abundance.push_back(v);}
        if(!valid)++invalid;tree.Fill();if(++total%2000==0)std::cout<<"GEOMETRY "<<total<<" invalid "<<invalid<<std::endl;
    }
    output.cd();tree.Write();output.Close();std::cout<<"DONE "<<total<<" invalid "<<invalid<<std::endl;return 0;
}
