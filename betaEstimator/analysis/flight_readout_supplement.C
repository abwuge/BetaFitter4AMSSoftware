// Extract measured readout with explicit native event/track association.
// Existing calibrated times are not replaced. No RICH or truth value is used.
#include <limits>
#include <vector>
#include <algorithm>
#define main native_producer_main
#include "../../amsd69n/analysis_amsd69n.C"
#undef main

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "usage: readout_supplement raw_files.txt targets.txt output.root\n";
        return 2;
    }
    AMSChain chain;
    std::ifstream inputs(argv[1]);
    std::string path;
    int nfile=0;
    while (std::getline(inputs,path)) {
        if (path.empty()) continue;
        if (chain.Add(path.c_str())<=0) return 3;
        ++nfile;
    }
    if (!nfile) return 3;
    AMSAnalysis analysis;
    analysis.LoadHeader();
    TFile output(argv[3],"CREATE");
    if (output.IsZombie()) return 4;
    TTree tree("readout","Measured TOF and Tracker readout, no fitted labels");
    unsigned int run=0,event=0;
    int found=0,selected=0,track_id=-1,beta_id=-1,requested_beta=-1;
    int tk_id[9],tk_mult[9],tk_nstrip[9][2],tk_offset[9][2],tk_seed[9][2];
    unsigned int tk_status[9],tof_status[4],tof_sidebit[4][2];
    int tof_bar[4];
    float tk_sum[9][2],tk_signedsum[9][2],tk_edep[9][3],tk_slope[9][2][2];
    float tof_adc[4][2],tof_dadc[4][2][3],tof_rawtime[4][2],tof_ht[4][2],tof_sht[4][2];
    float tof_aq2[4][2],tof_dq2[4][2][3],tof_edep[4][3],tof_coo[4][3];
    std::vector<float> strip_signal,strip_sigma;
    std::vector<int> strip_address,strip_status,strip_overflow;
    tree.Branch("run",&run,"run/i"); tree.Branch("event",&event,"event/i");
    tree.Branch("found",&found,"found/I"); tree.Branch("selected",&selected,"selected/I");
    tree.Branch("itrtrack",&track_id,"itrtrack/I");
    tree.Branch("ibetah",&beta_id,"ibetah/I");
    tree.Branch("requested_beta",&requested_beta,"requested_beta/I");
    tree.Branch("tk_id",tk_id,"tk_id[9]/I");
    tree.Branch("tk_mult",tk_mult,"tk_mult[9]/I");
    tree.Branch("tk_status",tk_status,"tk_status[9]/i");
    tree.Branch("tk_nstrip",tk_nstrip,"tk_nstrip[9][2]/I");
    tree.Branch("tk_offset",tk_offset,"tk_offset[9][2]/I");
    tree.Branch("tk_seed",tk_seed,"tk_seed[9][2]/I");
    tree.Branch("tk_sum",tk_sum,"tk_sum[9][2]/F");
    tree.Branch("tk_signedsum",tk_signedsum,"tk_signedsum[9][2]/F");
    tree.Branch("tk_edep",tk_edep,"tk_edep[9][3]/F");
    tree.Branch("tk_slope",tk_slope,"tk_slope[9][2][2]/F");
    tree.Branch("strip_signal",&strip_signal); tree.Branch("strip_sigma",&strip_sigma);
    tree.Branch("strip_address",&strip_address); tree.Branch("strip_status",&strip_status);
    tree.Branch("strip_overflow",&strip_overflow);
    tree.Branch("tof_barid",tof_bar,"tof_barid[4]/I");
    tree.Branch("tof_status",tof_status,"tof_status[4]/i");
    tree.Branch("tof_sidebit",tof_sidebit,"tof_sidebit[4][2]/i");
    tree.Branch("tof_adca",tof_adc,"tof_adca[4][2]/F");
    tree.Branch("tof_dadc",tof_dadc,"tof_dadc[4][2][3]/F");
    tree.Branch("tof_rtime",tof_rawtime,"tof_rtime[4][2]/F");
    tree.Branch("tof_rtime_ht",tof_ht,"tof_rtime_ht[4][2]/F");
    tree.Branch("tof_rtime_sht",tof_sht,"tof_rtime_sht[4][2]/F");
    tree.Branch("tof_aq2",tof_aq2,"tof_aq2[4][2]/F");
    tree.Branch("tof_dq2",tof_dq2,"tof_dq2[4][2][3]/F");
    tree.Branch("tof_edep",tof_edep,"tof_edep[4][3]/F");
    tree.Branch("tof_coo",tof_coo,"tof_coo[4][3]/F");
    std::ifstream targets(argv[2]);
    if (!targets) return 5;
    int total=0,missing=0,bad_association=0;
    const float nan=std::numeric_limits<float>::quiet_NaN();
    while (targets >> run >> event >> track_id >> requested_beta) {
        found=selected=0; beta_id=-1;
        std::fill(tk_id,tk_id+9,0); std::fill(tk_mult,tk_mult+9,-1);
        std::fill(tk_status,tk_status+9,0);
        std::fill(&tk_nstrip[0][0],&tk_nstrip[0][0]+18,0);
        std::fill(&tk_offset[0][0],&tk_offset[0][0]+18,-1);
        std::fill(&tk_seed[0][0],&tk_seed[0][0]+18,-1);
        std::fill(&tk_sum[0][0],&tk_sum[0][0]+18,nan);
        std::fill(&tk_signedsum[0][0],&tk_signedsum[0][0]+18,nan);
        std::fill(&tk_edep[0][0],&tk_edep[0][0]+27,nan);
        std::fill(&tk_slope[0][0][0],&tk_slope[0][0][0]+36,nan);
        std::fill(tof_bar,tof_bar+4,-1); std::fill(tof_status,tof_status+4,0);
        std::fill(&tof_sidebit[0][0],&tof_sidebit[0][0]+8,0);
        std::fill(&tof_adc[0][0],&tof_adc[0][0]+8,nan);
        std::fill(&tof_dadc[0][0][0],&tof_dadc[0][0][0]+24,nan);
        std::fill(&tof_rawtime[0][0],&tof_rawtime[0][0]+8,nan);
        std::fill(&tof_ht[0][0],&tof_ht[0][0]+8,nan);
        std::fill(&tof_sht[0][0],&tof_sht[0][0]+8,nan);
        std::fill(&tof_aq2[0][0],&tof_aq2[0][0]+8,nan);
        std::fill(&tof_dq2[0][0][0],&tof_dq2[0][0][0]+24,nan);
        std::fill(&tof_edep[0][0],&tof_edep[0][0]+12,nan);
        std::fill(&tof_coo[0][0],&tof_coo[0][0]+12,nan);
        strip_signal.clear();strip_sigma.clear();strip_address.clear();strip_status.clear();strip_overflow.clear();
        AMSEventR *ev=chain.GetEventFast(run,event);
        found=ev && ev->Run()==run && ev->Event()==event;
        if (!found) ++missing;
        else {
            ev->RecordRTIRun();
            analysis.unbiasetag0=1;analysis.unbiasetag=0;analysis.istkaddhit=0;
            analysis.InitEvent(ev); analysis.isreal=true; analysis.LoadHeaderEv();
            TrTrackR *track=ev->pTrTrack(track_id);
            if (track) {
                TofRecH::BuildOpt=0;TofRecH::ReBuild();
                for (int i=0;i<ev->nBetaH();++i) {
                    BetaHR *candidate=ev->pBetaH(i);
                    if(candidate && candidate->iTrTrack()==track_id) {beta_id=i;break;}
                }
                if (beta_id==requested_beta && beta_id>=0) {
                    selected=1;
                    BetaHR *beta=ev->pBetaH(beta_id);
                    for(int l=0;l<4;++l) {
                        TofClusterHR *hit=beta->GetClusterHL(l);
                        if(!hit) continue;
                        tof_bar[l]=hit->Bar;tof_status[l]=hit->Status;
                        for(int s=0;s<2;++s) {
                            tof_sidebit[l][s]=hit->SideBitPat[s];tof_adc[l][s]=hit->Aadc[s];
                            tof_rawtime[l][s]=hit->Rtime[s];tof_ht[l][s]=hit->Rhtime[s];
                            tof_sht[l][s]=hit->Rshtime[s];tof_aq2[l][s]=hit->AQ2[s];
                            for(int p=0;p<3;++p){tof_dadc[l][s][p]=hit->Dadc[s][p];tof_dq2[l][s][p]=hit->DQ2[s][p];}
                        }
                        for(int p=0;p<3;++p){tof_edep[l][p]=hit->GetEdep(p);tof_coo[l][p]=hit->Coo[p];}
                    }
                    for(int l=0;l<9;++l) {
                        TrRecHitR *hit=track->GetHitLJ(l+1);
                        if(!hit) continue;
                        tk_id[l]=hit->GetTkId();tk_mult[l]=hit->GetResolvedMultiplicity();tk_status[l]=hit->getstatus();
                        for(int s=0;s<2;++s) {
                            TrClusterR *cl=s==0?hit->GetXCluster():hit->GetYCluster();
                            if(!cl) continue;
                            tk_nstrip[l][s]=cl->GetNelem();tk_offset[l][s]=strip_signal.size();tk_seed[l][s]=cl->GetSeedIndex(0);
                            tk_sum[l][s]=cl->GetTotSignal(0);tk_signedsum[l][s]=0;
                            tk_slope[l][s][0]=cl->GetDxDz();tk_slope[l][s][1]=cl->GetDyDz();
                            for(int j=0;j<cl->GetNelem();++j) {
                                float signal=cl->GetSignal(j,0);tk_signedsum[l][s]+=signal;
                                strip_signal.push_back(signal);strip_sigma.push_back(cl->GetSigma(j));
                                strip_address.push_back(cl->GetAddress(j));strip_status.push_back(cl->GetStatus(j));
                                strip_overflow.push_back(cl->IsOverflowADC(j,0)?1:0);
                            }
                        }
                        // These calibrated diagnostics are kept separate from raw ADC.
                        for(int s=0;s<3;++s) tk_edep[l][s]=hit->GetEdepH(s)*1.e-3;
                    }
                } else ++bad_association;
            } else ++bad_association;
        }
        tree.Fill();
        if (++total%25==0) std::cout<<"Readout "<<total<<" missing "<<missing<<" association_fail "<<bad_association<<std::endl;
    }
    output.cd();tree.Write();output.Close();
    std::cout<<"DONE requested="<<total<<" missing="<<missing<<" association_fail="<<bad_association<<std::endl;
    return (missing||bad_association)?6:0;
}
