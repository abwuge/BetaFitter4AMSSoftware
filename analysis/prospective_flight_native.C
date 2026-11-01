// Native measured supplements for the frozen ADC/time estimators.
// Reuse the verified light refit path. No RICH object or MC label is accessed.
#include <limits>
#include <vector>
#include <algorithm>
#define main native_producer_main
#include "../../amsd69n/analysis_amsd69n.C"
#undef main

int main(int argc, char **argv) {
    if (argc != 4) { std::cerr << "usage: prospective_native raw_files targets output.root\n"; return 2; }
    AMSChain chain; std::ifstream input(argv[1]); std::string path; int files=0;
    while (std::getline(input, path)) { if(path.empty()) continue; if(chain.Add(path.c_str())<=0) return 3; ++files; }
    if(!files) return 3;
    AMSAnalysis analysis; analysis.LoadHeader();
    TFile output(argv[3], "CREATE"); if(output.IsZombie()) return 4;
    TTree tree("measured", "Native association, Kalman rigidity and time/ADC checks");
    unsigned int run=0,event=0;
    int itrtrack=-1,requested_beta=-1,ibetah=-1,found=0,selected=0,fit=-1;
    int used[4],good[4],bar[4],pattern[4],sumhit=0,usehit=0,bpatt=0;
    float R[4],sigma_r,chi[3],beta,time[4],etime[4],length[4],adc[4][2],edep[4];
    double R_L5=0;
    const float nan=std::numeric_limits<float>::quiet_NaN();
#define S(NAME,FMT) tree.Branch(#NAME,&NAME,#NAME FMT)
#define A(NAME,FMT) tree.Branch(#NAME,NAME,#NAME FMT)
    S(run,"/i");S(event,"/i");S(itrtrack,"/I");S(requested_beta,"/I");S(ibetah,"/I");S(found,"/I");S(selected,"/I");S(fit,"/I");
    S(sumhit,"/I");S(usehit,"/I");S(bpatt,"/I");S(sigma_r,"/F");S(beta,"/F");S(R_L5,"/D");
    A(R,"[4]/F");A(chi,"[3]/F");A(used,"[4]/I");A(good,"[4]/I");A(bar,"[4]/I");A(pattern,"[4]/I");
    A(time,"[4]/F");A(etime,"[4]/F");A(length,"[4]/F");A(adc,"[4][2]/F");A(edep,"[4]/F");
#undef S
#undef A
    std::ifstream targets(argv[2]); if(!targets) return 5;
    int total=0,missing=0,rejected=0;float charge_measurement;
    while(targets>>run>>event>>itrtrack>>requested_beta>>charge_measurement) {
        found=selected=0;ibetah=fit=-1;sumhit=usehit=bpatt=0;sigma_r=beta=nan;R_L5=nan;
        std::fill(R,R+4,nan);std::fill(chi,chi+3,nan);std::fill(used,used+4,0);std::fill(good,good+4,0);
        std::fill(bar,bar+4,-1);std::fill(pattern,pattern+4,-1);std::fill(time,time+4,nan);std::fill(etime,etime+4,nan);
        std::fill(length,length+4,nan);std::fill(edep,edep+4,nan);std::fill(&adc[0][0],&adc[0][0]+8,nan);
        AMSEventR *ev=chain.GetEventFast(run,event);found=ev&&ev->Run()==run&&ev->Event()==event;
        if(found) {
            ev->RecordRTIRun();analysis.unbiasetag0=1;analysis.unbiasetag=0;analysis.istkaddhit=0;
            analysis.InitEvent(ev);analysis.isreal=true;analysis.LoadHeaderEv();analysis.itrtrack=itrtrack;
            TrTrackR *track=ev->pTrTrack(itrtrack);
            if(track) {
                TofRecH::BuildOpt=0;TofRecH::ReBuild();
                for(int i=0;i<ev->nBetaH();++i) {BetaHR *candidate=ev->pBetaH(i);if(candidate&&candidate->iTrTrack()==itrtrack){ibetah=i;break;}}
                analysis.ibetah=ibetah;
                BetaHR *b=ev->pBetaH(ibetah);
                if(b&&ibetah==requested_beta) {
                    sumhit=b->GetSumHit();usehit=b->GetUseHit();bpatt=b->GetBetaPattern();beta=b->GetBeta();
                    for(int l=0;l<4;++l) {
                        used[l]=b->IsBetaUseHL(l);pattern[l]=b->GetPattern(l);length[l]=b->GetTkTFLen(l);
                        TofClusterHR *hit=b->GetClusterHL(l);if(!hit)continue;
                        good[l]=hit->IsGoodTime();bar[l]=hit->Bar;time[l]=hit->Time;etime[l]=hit->ETime;edep[l]=hit->GetEdep();
                        for(int side=0;side<2;++side)adc[l][side]=hit->Aadc[side];
                    }
                    int refit=0,charge=0;float mass=0;
                    analysis.InitTrkReFit(charge_measurement,track->GetAdvancedFitCharge(),refit,mass,charge,6);
                    TrFit trfit;
                    for(int version=0;version<3;++version) {
                        int option=(version==0?0:version==1?20:30)+refit%10;
                        fit=track->iTrTrackPar(trfit,6,3,option,mass,charge);
                    }
                    if(fit>=0) {
                        for(int z=0;z<4;++z)R[z]=track->GetRigidity(fit,z);
                        sigma_r=track->GetErrRinv(fit);chi[0]=track->GetNormChisqX(fit);chi[1]=track->GetNormChisqY(fit);chi[2]=track->GetChisq(fit);
                        AMSPoint p5;AMSDir d5;R_L5=R[1];trfit.InterpolateKalman(1.698,p5,d5,R_L5);
                        selected=1;
                    }
                }
            }
        }else ++missing;
        if(!selected)++rejected;tree.Fill();++total;
        if(total%100==0)std::cout<<"NATIVE "<<total<<" missing "<<missing<<" rejected "<<rejected<<std::endl;
    }
    output.cd();tree.Write();output.Close();std::cout<<"DONE "<<total<<" missing "<<missing<<" rejected "<<rejected<<std::endl;
    return missing?6:0;
}
