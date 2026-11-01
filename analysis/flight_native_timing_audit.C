// Inspect native timing associations and every calibrated discriminator pulse.
// No RICH or simulation label is read; no production object is written back.
#include <limits>
#include <vector>
#include <algorithm>
#define main native_producer_main
#include "../../amsd69n/analysis_amsd69n.C"
#undef main

int main(int argc,char **argv) {
    if(argc!=4) {std::cerr<<"usage: timing_audit raw_files targets output.root\n";return 2;}
    AMSChain chain; std::ifstream inputs(argv[1]);std::string path;
    int nfiles=0;
    while(std::getline(inputs,path)) {if(path.empty())continue;if(chain.Add(path.c_str())<=0)return 3;++nfiles;}
    if(!nfiles)return 3;
    AMSAnalysis analysis;analysis.LoadHeader();
    TFile output(argv[3],"CREATE");if(output.IsZombie())return 4;
    TTree tree("timing","Native timing association audit");
    unsigned int run=0,event=0;int track_id=-1,requested_beta=-1,beta_id=-1,found=0,selected=0;
    int usehit=0,sumhit=0,buildtype=0,beta_status=0;
    float beta_value=0,beta_error=0,t0=0,chi2=0;
    int bar[4],pattern[4],used[4],exists[4],goodtime[4],pulse_code[4];
    unsigned int status[4],sidebit[4][2];
    float time[4],etime[4],beta_time[4],beta_etime[4],pathlen[4],residual[4];
    float corrected[4][2][3],raw[4][2][3],adc[4][2];
    int pulse_offset[4][2][3],pulse_count[4][2][3],rawside_index[4][2];
    std::vector<float> pulse_time;
#define B(NAME,FORMAT) tree.Branch(#NAME,NAME,#NAME FORMAT)
    tree.Branch("run",&run,"run/i");tree.Branch("event",&event,"event/i");
    tree.Branch("itrtrack",&track_id,"itrtrack/I");tree.Branch("ibetah",&beta_id,"ibetah/I");
    tree.Branch("requested_beta",&requested_beta,"requested_beta/I");
    tree.Branch("found",&found,"found/I");tree.Branch("selected",&selected,"selected/I");
    tree.Branch("usehit",&usehit,"usehit/I");tree.Branch("sumhit",&sumhit,"sumhit/I");
    tree.Branch("buildtype",&buildtype,"buildtype/I");tree.Branch("beta_status",&beta_status,"beta_status/I");
    tree.Branch("beta_value",&beta_value,"beta_value/F");tree.Branch("beta_error",&beta_error,"beta_error/F");
    tree.Branch("t0",&t0,"t0/F");tree.Branch("chi2",&chi2,"chi2/F");
    B(bar,"[4]/I");B(pattern,"[4]/I");B(used,"[4]/I");B(exists,"[4]/I");B(goodtime,"[4]/I");B(pulse_code,"[4]/I");
    B(status,"[4]/i");B(sidebit,"[4][2]/i");B(time,"[4]/F");B(etime,"[4]/F");B(beta_time,"[4]/F");B(beta_etime,"[4]/F");B(pathlen,"[4]/F");B(residual,"[4]/F");
    B(raw,"[4][2][3]/F");B(corrected,"[4][2][3]/F");B(adc,"[4][2]/F");
    B(pulse_offset,"[4][2][3]/I");B(pulse_count,"[4][2][3]/I");B(rawside_index,"[4][2]/I");
    tree.Branch("pulse_time",&pulse_time);
#undef B
    std::ifstream targets(argv[2]);if(!targets)return 5;
    const float nan=std::numeric_limits<float>::quiet_NaN();int total=0,failed=0;
    while(targets>>run>>event>>track_id>>requested_beta) {
        found=selected=0;beta_id=-1;usehit=sumhit=buildtype=beta_status=0;beta_value=beta_error=t0=chi2=nan;
        std::fill(bar,bar+4,-1);std::fill(pattern,pattern+4,-1);std::fill(used,used+4,0);std::fill(exists,exists+4,0);std::fill(goodtime,goodtime+4,0);std::fill(pulse_code,pulse_code+4,-999);
        std::fill(status,status+4,0);std::fill(&sidebit[0][0],&sidebit[0][0]+8,0);
        std::fill(time,time+4,nan);std::fill(etime,etime+4,nan);std::fill(beta_time,beta_time+4,nan);std::fill(beta_etime,beta_etime+4,nan);std::fill(pathlen,pathlen+4,nan);std::fill(residual,residual+4,nan);
        std::fill(&raw[0][0][0],&raw[0][0][0]+24,nan);std::fill(&corrected[0][0][0],&corrected[0][0][0]+24,nan);std::fill(&adc[0][0],&adc[0][0]+8,nan);
        std::fill(&pulse_offset[0][0][0],&pulse_offset[0][0][0]+24,-1);std::fill(&pulse_count[0][0][0],&pulse_count[0][0][0]+24,0);std::fill(&rawside_index[0][0],&rawside_index[0][0]+8,-1);pulse_time.clear();
        AMSEventR *ev=chain.GetEventFast(run,event);found=ev&&ev->Run()==run&&ev->Event()==event;
        if(found) {
            ev->RecordRTIRun();analysis.unbiasetag0=1;analysis.unbiasetag=0;analysis.istkaddhit=0;
            analysis.InitEvent(ev);analysis.isreal=true;analysis.LoadHeaderEv();
            TrTrackR *track=ev->pTrTrack(track_id);
            if(track) {
                TofRecH::BuildOpt=0;TofRecH::ReBuild();
                for(int i=0;i<ev->nBetaH();++i){BetaHR *b=ev->pBetaH(i);if(b&&b->iTrTrack()==track_id){beta_id=i;break;}}
                if(beta_id>=0&&beta_id==requested_beta) {
                    selected=1;BetaHR *b=ev->pBetaH(beta_id);
                    usehit=b->GetUseHit();sumhit=b->GetSumHit();buildtype=b->GetBuildType();beta_status=b->GetStatus();
                    beta_value=b->GetBeta();beta_error=b->GetEBetaV();t0=b->GetT0();chi2=b->GetChi2T();
                    for(int l=0;l<4;++l) {
                        pattern[l]=b->GetPattern(l);used[l]=b->IsBetaUseHL(l);exists[l]=b->TestExistHL(l);
                        beta_time[l]=b->GetTime(l);beta_etime[l]=b->GetETime(l);pathlen[l]=b->GetTkTFLen(l);residual[l]=b->GetTResidual(l);
                        TofClusterHR *hit=b->GetClusterHL(l);if(!hit)continue;
                        bar[l]=hit->Bar;time[l]=hit->Time;etime[l]=hit->ETime;goodtime[l]=hit->IsGoodTime();status[l]=hit->Status;
                        std::vector<float> lt[2],ht[2],sht[2];pulse_code[l]=hit->GetSTimeAll(lt,ht,sht);
                        for(int s=0;s<2;++s) {
                            sidebit[l][s]=hit->SideBitPat[s];rawside_index[l][s]=hit->iTofRawSide(s);adc[l][s]=hit->Aadc[s];
                            raw[l][s][0]=hit->Rtime[s];raw[l][s][1]=hit->Rhtime[s];raw[l][s][2]=hit->Rshtime[s];
                            corrected[l][s][0]=hit->Stime[s];corrected[l][s][1]=hit->Shtime[s];corrected[l][s][2]=hit->Sshtime[s];
                            for(int q=0;q<3;++q) {const std::vector<float>&v=q==0?lt[s]:(q==1?ht[s]:sht[s]);pulse_offset[l][s][q]=pulse_time.size();pulse_count[l][s][q]=v.size();pulse_time.insert(pulse_time.end(),v.begin(),v.end());}
                        }
                    }
                }
            }
        }
        if(!selected)++failed;tree.Fill();++total;
    }
    output.cd();tree.Write();output.Close();std::cout<<"DONE "<<total<<" failed "<<failed<<std::endl;return failed?6:0;
}
