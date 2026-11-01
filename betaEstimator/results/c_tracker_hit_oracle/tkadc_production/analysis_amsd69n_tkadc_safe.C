///--Author Qi Yan (qyan@cern.ch)
///--Processing DST for AMS data analysis
#include <signal.h>
#include <TH1.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLeaf.h>
#include <TGraphErrors.h>
#include <TApplication.h>
#include <TChain.h>
#include <TLine.h>
#include <TSystem.h>
#include <TStyle.h>
#include <iomanip>
#include "TMath.h"
#include "TMinuit.h"
#include "TTree.h"
#include <sstream>
#include "TString.h"
#include <fstream>
#include <iostream>
#include "root_RVSP.h"
#include "amschain.h"
#include "Tofrec02_ihep.h"
#include "TrdKCluster.h"
#include "bcorr.h"
#include <time.h>
#include "Track.h"
#include "tkdcards.h"
#include "TrInnerDzDB.h"
#include "TrExtAlignDB.h"
#include "AbsKalmanFitter.h"
#include "KalmanFitStatus.h"
#include "KalmanFittedStateOnPlane.h"
#include "KalmanFitterInfo.h"
//#include "/afs/cern.ch/work/q/qyan/AMSVDEV/QTool/readfile3.C"
#include "readfile3.C" //2023.12.12
#include "TofCalib.h" //use my TOF Time Calibration parameter //2024.08.24
//#define _USEFAST_
#ifdef _HEINNERPRESCALE_
#define _HEL1PRESCALE_
#endif
#ifdef _IONINNERPRESCALE_
#define _IONL1PRESCALE_
#endif
#if defined (_USENOLINEARCOR_) || defined (_USECALIB_)
#define _USETKADC_ 
#define _GOODTKTOFPAR_
#define _PRHESAVEUNBIASE_
#endif

#ifndef _USETKADC_
#define _USETKADC_
#endif

//--2024.11.17
#include "Tofsim02.h"

//2024.08.24
//#define MYTOFTIMECALIB
const int NMyTOFTimeCalibVer = 2;

//2024.09.05
//#define COMBINE_Z_TOFTIMECALIB
const int NTOFTimeCalibVer = 5;

//TOF Beta without calibration, used as reference on the performance of TOF Time Calibration //2024.09.05
// #define TOFBETA_WITHOUT_CAL
const int NLTHTSHT = 3;

//2024.10.23
// #define NEWTOFTIMECALIB
//const int NZ_NewCalib = 3;
//const int ZList_NewCalib[] = {3, 6, 8};
//const int NZ_NewCalib = 2; //=0: -1, =1: mc_charge //2024.11.04
//const int NZ_NewCalib = 9; //=0: -1, =1: mc_charge //2024.11.04
const int NZ_NewCalib = 3; //=0: -1, i.e. 2-8, =1: -2, i.e. 3-8, =2: mc_charge, first read the calibration for all the charge and read correspondingly  //2024.11.04
const int ZList_NewCalib[] = {-1, -2, 2, 3, 4, 5, 6, 7, 8};

//const int NMinChis_NewCalib = 3; //0: MinChis, 1: MinChisDer_DerQuadSum, 2: MinChisDer
//const int NMinChis_NewCalib = 1; //0: MinChis, 1: MinChisDer_DerQuadSum, 2: MinChisDer
//2024.12.10
const int NMinChis_NewCalib = 2; //0: MinChis, 1: MinChis, combine Min of Time and Coo, 2: MinChisDer_DerQuadSum, 3: MinChisDer
const int NExcludeEdge_NewCalib = 2;

//const int NFix_NewCalib = 2;
const int NFix_NewCalib = 3;
const double dFix_NewCalib[] = {0.3, 0.4, 0.5};

namespace tofconst{
  const int NECALL=18;
  const int NTKL=9;
  const int NTKS=2;
  const int NTKDIS=4;
  const float TKDIS[NTKDIS]={1,2,4,6};//dis=cm
  const int NTOFL=4;
  const int NTOFBM=10;
  const int NTOFS=2;
  const int NBAR[NTOFL]={8,8,10,8};
  const int TRANS[NTOFL]={0,1,1,0};//trans pos
  const float TOFZ[NTOFL]={65.2,62.1,-62.1,-65.2};//central cm  
  const float Sci_pt[NTOFL][NTOFBM]={
      -45.5,  -28.75, -17.25, -5.75, 5.75, 17.25, 28.75, 45.5,  0., 0.,
      -47.,   -28.75, -17.25, -5.75, 5.75, 17.25, 28.75, 47.,   0., 0.,
-55., -40.25, -28.75, -17.25, -5.75, 5.75, 17.25, 28.75, 40.25, 55.,
      -47.25, -28.75, -17.25, -5.75, 5.75, 17.25, 28.75, 47.25, 0., 0.,
   };
   const float qxlow[]={0.9, 1.9, 2.85, 3.8, 4.7, 5.6, 6.65, 7.6};
   const float qxhig[]={1.15,2.2, 3.2,  4.2, 5.2, 6.2, 7.2, 8.2};
   const float tk_pz[NTKL]={159.04,53.05,29.22,25.24,1.706,-2.292,-25.26,-29.25,-136.03};
#ifdef _USENEWL1L9G_
   const double cirr[]={62+2,62+2,46+2,46+2,46+2,46+2,46+2,46+2,43+2};//New L1L9 Geometry(+-2cm)
   const double ciry[]={47+2,40+2,44+2,44+2,36+2,36+2,44+2,44+2,29+2};
#else
   const double cirr[]={62,62,46,46,46,46,46,46,43};
   const double ciry[]={47,40,44,44,36,36,44,44,29};
#endif
//--Static information
#ifndef _EVENTSCALE_ 
   const int prscale=1;
   const int prunbias=20;//1/20 for unbias
#else
   const int prscale=10;
   const int prunbias=200;//1/20 for unbias
#endif
   int    prunbiasmc=10;
   unsigned int prtime=0;
   int    sstat[100]={0};
   int    sstree[100]={0};
   double cput[100]={0};
   //smear betah //2024.01.03
   #ifdef _DOBETAMCTUNE_
   #ifdef _NEWTOFMCTUNE_ //2024.01.30
   //const int NDT = 5;
   //const float MCTuneDT[NDT] = {-10, -15, -20, -30, -40};
   const int NDT = 10;
   const float MCTuneDT[NDT] = {-30, -20, -15, -10, -5, 5, 10, 15, 20, 30};
   const int NST = 1;
   const float MCTuneST[NST] = {0};
   #else
   const int NDT = 5;
   const int NST = 5;
   //const int MCTuneST=20;
   //const int NST = 2;
   #ifdef _FINETUNE_
   const float MCTuneDT[NDT] = {-25, -27.5, -30, -32.5, -35};
   const float MCTuneST[NST] = {5, 7.5, 10, 12.5, 15};
   #else
   const float MCTuneDT[NDT] = {0, -20, -40, -60, -80};
   const float MCTuneST[NST] = {-40, -20, 0, 20, 40};
   #endif //_FINETUNE_
   #endif //_NEWTOFMCTUNE_
   #endif //_DOBETAMCTUNE_
}
using namespace tofconst;

//TofCalib *tofcalib = NULL;
TofCalib *tofcalib[NMyTOFTimeCalibVer]; //2024.09.16
TofCalib tofcalibCombZ[NTOFTimeCalibVer]; //2024.09.05
TofCalib tofcalibNull; //2024.09.05
//2024.10.23
//TofCalib tofcalibNewCalib[9][NMinChis_NewCalib][NExcludeEdge_NewCalib];
TofCalib tofcalibNewCalib[NZ_NewCalib][NMinChis_NewCalib][NExcludeEdge_NewCalib];

class AMSAnalysis{
  public:
   AMSAnalysis():fout(0),pev(0){tout[0]=tout[1]=tout[2]=0;};
   ~AMSAnalysis(){};
   void BookFile(char *ofile);
   int  LoadHeader();
   int  LoadHeaderEv();
   int  LoadHeaderMC(int mcz=0,int vers=0,int mparid=0);
   int  InitTrkReFit(float tkiq,float tkfz,int &refit,float &mass,int &chrg,int algo=1);
//---
   int  GetTreeSize(TTree *tree);
   int  AddBranch_Header (TTree *tree);
   int  AddBranch_Trigger(TTree *tree);
   int  AddBranch_RTI(TTree *tree);
   int  AddBranch_MC (TTree *tree);
   int  AddBranch_Particle(TTree *tree);
   int  AddBranch_Trk (TTree *tree);
   int  AddBranch_Tof (TTree *tree);
   int  AddBranch_Trd (TTree *tree);
   int  AddBranch_Rich(TTree *tree);
   int  AddBranch_Ecal(TTree *tree);
   int  AddBranch_Anti(TTree *tree);
   int  AddBranch_TrkParticle(TTree *tree);
//---
   int  Save();
   void InitEvent(AMSEventR *ev);   
   int  Select();//1 Write-Pre 2 Write-All
   bool Select_Header();
   bool Select_Trigger();
   bool Select_MC(int opt=0);
   bool Select_RTI();
   bool Select_Particle();
//----*
   bool Select_L1L9Particle();
   bool Select_MCL1L9Particle();
   bool Select_2ndParticle();
   bool Select_LinkParticle();
//----
   bool Select_Trk(int sel=0);//pre selction 1 or 0
   bool Select_Tof();
   bool Select_Trd();
   bool Select_Rich();
   bool Select_Ecal();  
   bool Select_Anti();
   bool Select_TrkParticle();
   int  Fill(int it);
//--
   TFile *fout;
   TTree *tout[3];

   AMSEventR *pev;
   int iparindex;
//--Sum Control
   int          unbiasetag0;
   int          unbiasetag;
   int          issmearext;
   int          istkaddhit;
//--Header
   bool         isreal;
   unsigned int run;
   unsigned int event;
   int          amsdn;
   int          version;
   unsigned int errorb;
   unsigned int time[2];
   float        yaw;
   float        pitch;
   float        roll;
   float        thetam;//magnetic theta
   float        phim;
   bool         isbadrun;

//---Trigger
   int          nlevel1;
   int          ndaqev;
   int          ndaqerr;
   unsigned int daqlen;
   int          daqtyerr;
   unsigned int daqjrmerr;
   unsigned int daqjlen[24];
   unsigned int daqsdlen;
   unsigned int daqsdlenu;
   float        livetime;
   int          physbpatt;
   int          jmembpatt;
   int          physbpatt1;
   int          jmembpatt1;
   int          physbpatt2;
   int          jmembpatt2;
   int          ecalflag;
   int          tofflag[2];
   int          antipatt;

//---RTI
   int          irti;
   float        zenith;
   float        thetas;
   float        phis;
   float        rads;
   float        glat;//galatic latti
   float        glong;//galatic longti
   bool         issaa;
   float        mcutoff[4][2];
   float        mcutoffi[4][2];
   float        rtilf;
   float        rtinev;
   float        rtinerr;
   float        rtintrig;
   float        rtinpar;
   int          rtigood;
   float        rtinexl[2][2];
   float        rticdif[2][3];

//--MC
   int          mpar; // MC particle number
   float        mmom; // MC particle momentum
   float        mch;  // MC particle charge
   float        mmass; // MC particle TOF
   int          mtof_pass; // MC TOF pass (1111: all pass, 0000: all fail, left to right: 1st, 2nd, 3rd, 4th)
   float        mevcoo[3];
   float        mevdir[3];
   float        mfcoo[3];
   float        mfdir[3];
   float        mpare[2];
   int          mparp[2];
   float        mparc[2];
   float        msume;
   float        msumc;
   float        mtrdmom[20];
   int          mtrdpar[20];
   float        mtrdcoo[20][3];
   int          ntrmccl;
   float        mtrmom[9];
   int          mtrpar[9];
   float        mtrcoo[9][3];
   float        mtrcoo1[9][3];
   int          mtrpri[9];
   int          mtrz[9];
   int          ntofmccl;
   float        mtofdep[4];
   int          mtofpri[4];
   int          mtofbar[4];
//---
/**
 * Detector Layers Z Positions (cm):
 * Layer  Position  Label
 *  0     158.92    L1
 *  1     141.20      
 *  2     86.10       
 *  3     64.425      
 *  4     65.975      
 *  5     61.325      
 *  6     62.875      
 *  7     53.06     L2
 *  8     29.228    L3
 *  9     25.212    L4
 * 10     1.698     L5
 * 11    -2.318     L6
 * 12    -25.212    L7
 * 13    -29.228    L8
 * 14    -62.875      
 * 15    -61.325      
 * 16    -64.425      
 * 17    -69.975      
 * 18    -135.882   L9
 * 19    -142.792     
 * 20    -158.457     
 */
   float        mevcoo1[21][3]; // MC particle coordinate in each layer
   float        mevdir1[21][3]; // MC particle direction in each layer
   float        mevmom1[21];    // MC particle momentum in each layer

//---L1L9Particle
   int          l1i;
   int          l9i;
   float        tk_l1q;
   float        tk_l9q;
   float        tk_l1qxy[2];
   float        tk_l9qxy[2];
   int          tk_l1qs;
   int          tk_l9qs;
   int          ibetahs;
   int          itrdtracks;
   int          tof_nhits;
   int          tof_hsumhus;
   float        betahs;
   float        tof_chiscs;
   float        tof_chists;
   float        tof_qls[4];
   float        trd_thetas;
   float        trd_phis;
   float        trd_coos[3];
   float        trd_rqs;
   float        ecal_pos[3];
   float        ecal_ens;
   float        ecal_dis;
   float        tk_l1mds;
   float        tk_l9mds;
   float        tk_pos1s[9][3];
   float        tk_dir1s[3];
   float        tk_l1qvs;
   float        tk_l9qvs;
   float        tk_exqvn[2][2];//YJ+HL
   float        tk_l1qvr;
   float        tk_l9qvr;
   float        cutoffps; 

//---Particle
   int          nparticle;
   int          nbetah;
   int          nbeta;
   int          ncharge;
   int          ntrack;
   int          ntrrawcl;
   int          ntrcl;
   int          ntrhit;
   int          ntrdtrack;
   int          nshow;
   int          necalhit;
   int          ntofclh;
   int          ntrdseg;
   int          ntrdrawh;
   int          ntrdcl;
   int          nrich;
   int          nrichb;
   int          nrichh;
   int          btstat;
   int          ibeta;
   int          ibetah;
   int          itrtrack;
   int          itrdtrack;
   int          icharge;
   int          ishow;
   int          irich;
   int          irichb;
   float        show_dis;
   float        ec_pos[3];
   float        ec_dir[3];
   float        cutoffp[2];
   float        cutoffpi[2];
//---2nd Particle
  int          betah2p;//Pattern
  float        betah2q;//Q
  float        betah2r;//Rigidity
  int          betah2hb[2];//Hit
  float        betah2ch;//Chis
   
//---Track
   bool         tk_l1l9par;
   bool         tk_bgpar;
   bool         tk_bhgpar;
   bool         tk_ubpar;
   bool         tk_ngpar;
   int          tk_rtype;
   int          tk_nhiti[2];
   int          tk_hitb[2];//XY
   float        tk_iq;
   float        tk_1q;
   float        tk_9q;
   float        tk_q[2];
   float        tk_qr;//norig and beta correction
   float        tk_qrms[2];
   int          tk_qhit[2];
   int          tk_rz[2][2];
   float        tk_rpz[2][2];
   float        tk_ql[9];//Lay1+2+..9
   float        tk_edep[9];//energy deposition of track-associated hit (MeV)
   float        tk_ql2[9][2];
   int          tk_qls[9];
   int          tk_qlsn;
   float        tk_qin[2][3];
   float        tk_qrmn[2][3];
   float        tk_qln[2][9][3];//Lay1+2+..9
   float        tk_qlqs[9][2][5];
   int          tk_lid[9];
   int          tk_cad[9];
   int          tk_stripx[9];
   float        tk_adc[9][2][11];//L+XY+ii
   float        tk_cof[9][2];//L+XY
   float        tk_iso[9][2];//Y stripe isolation
//---External Two Layers
   float        tk_exql[2][3];
   float        tk_exqln[2][2][3];//YJ+LU
   int          tk_exqls[2];
   float        tk_exdis[2];
   int          tk_exlid[2];
   int          tk_excad[2];
//---Full+Inner+InnerL0+InnerL8+InnerUP+InnerDo
   float        tk_fzm[2];
   float        tk_rigidity1[3][3][7];//Algo+PGCIEMAT+Span
   float        tk_rigidity;
   float        tk_rigidityk[4][7];//Kalman/Z+Span
   float        tk_erigidity1[3][7];//Algo+Span
   float        tk_chis1[3][3][7][3];////Algo+PGCIEMAT+Span+XYA
   float        tk_dir0[3][3][4];//Algo+PGCIEMAT+Span+X
   float        tk_bcor;
   float        tk_res[9][2];//Residual Layer+XY
   float        tk_res1[9][2];//Residual Layer+XY
   float        tk_rigidityi[9];
   float        tk_chisi[9];
   float        tk_cdif[2][3];//PG-CIEMAT Layer19+XYZ
   float        tk_pos[9][3];
   float        tk_dir[9][3];
   float        tk_hitc[9][3];
   float        tk_oq[2];
   float        tk_oel[9]; 
   int          tk_ohitl[9];
//-----Kalman
   float        tk_rigidityr[4][3][7];//Z+Algo+Span
   float        tk_posr[4][3][7][3];//Z+Algo+Span+P
   float        tk_dirr[4][3][7][3];//Z+Algo+Span+P

//---TOF
   int          tof_btype;
   int          tof_bpatt;
   float        tof_beta;
   float        tof_betah; // BetaHR reconstructed beta
   float        tof_ebetah;
   float        tof_ebetah_n;
   float        tof_betahmc;
   int          tof_hsumh;
   int          tof_hsumhu;
   float        tof_chist;
   float        tof_chisc;
   float        tof_chist_n;
   float        tof_chisc_n;
   float        tof_q[2];
   float        tof_qrms[2];
   int          tof_z;
   float        tof_pz;
   int          tof_zr;
   float        tof_pzr;
   int          tof_zud[2];
   float        tof_pud[2];
   int          tof_pass;
   int          tof_qs;//QStat 
   int          tof_nclhl[4];
   int          tof_barid[4];
   float        tof_pos[4][3];
   float        tof_dir[4][3];
   float        tof_ql[4];
   float        tof_ql1[4];
   float        tof_qlc[6][4];
   float        tof_qlpm[4][2][4];
   float        tof_qlr[3][4];//All+Anode+Dynode From BetaH
   float        tof_tl[4];  // Tof counter hit time
   float        tof_etl[4]; // Tof counter hit time error
   float        tof_tres[4]; //YC //2024.01.03
   float        tof_cres[4];
   float		tof_leng[4]; //YC //2024.06.02
   unsigned int tof_sbit[4][2];
   float        tof_oq[4][2];//Two Max Non-BetaH ClusterH
   int          tof_ob[4][2];
   int          tof_oncl;
   //beta using HT and SHT //2024.06.27
   float        tof_betah_HT;
   float        tof_betah_SHT;
   //time from HT and SHT //2024.07.02
   float        tof_tl_HT[4];
   float        tof_tl_SHT[4];
   //beta using Only S2S3, which should not influenced by energy loss. Store using S1S2S3S4 and only S1S4 as well to compare //2024.06.27
   float        tof_betah_S2S3;
   float        tof_betah_S1S4;
   float        tof_betah_S1S2S3S4;
   //add S1S2 and S3S4 as well, the difference due to energy loss should be largest for them //2024.07.02
   float        tof_betah_S1S2;
   float        tof_betah_S3S4;
   //energy loss in each layer //2024.06.27
   float        tof_edep[4];
   //calculate the 1st iteraion of beta change due to energy loss, and refit the 2nd iteration beta //2024.07.02
   float		tof_deltaBetah[10][4];
   float		tof_betah_iter[10];
   float		tof_deltaBetah_HT[10][4];
   float		tof_betah_iter_HT[10];
   float		tof_deltaBetah_SHT[10][4];
   float		tof_betah_iter_SHT[10];
   //--calculate the beta change from energy loss directly, using average mass, and then refit //2024.07.22
   float		tof_betah_el;
   float		tof_deltaBetah_el[4];
   float		tof_deltaT_el[4]; //time change due to energy loss
   //store the residual as well //2024.07.29
   float        tof_tres_el[4];
   float        tof_cres_el[4];
   //store the residual, chis, for HT/SHT as well //2024.07.04
   float        tof_tres_HT[4];
   float        tof_cres_HT[4];
   float        tof_chist_HT;
   float        tof_chisc_HT;
   float        tof_chist_n_HT;
   float        tof_chisc_n_HT;
   float        tof_tres_SHT[4];
   float        tof_cres_SHT[4];
   float        tof_chist_SHT;
   float        tof_chisc_SHT;
   float        tof_chist_n_SHT;
   float        tof_chisc_n_SHT;
   //combine beta fit using LT, HT and SHT //2024.07.28
   float		tof_betah_LTHTSHT;
   //--store the nHT and nSHT //2024.09.16
   int			nHTime[4][2];
   int			nSHTime[4][2];
   float		tof_rtime[4][2];
   float		tof_rtime_ht[4][2];
   float		tof_rtime_sht[4][2];
   //--store Anode adc //2024.09.19
   float		tof_Aadc[4][2];
   //--tof using my calibration //2024.08.23
   #ifdef MYTOFTIMECALIB
   float        tof_betah_MyCalib[NMyTOFTimeCalibVer][3][2];
   float        tof_tl_MyCalib[NMyTOFTimeCalibVer][3][2][4];
   float        tof_tres_MyCalib[NMyTOFTimeCalibVer][3][2][4];
   float        tof_cres_MyCalib[NMyTOFTimeCalibVer][3][2][4];
   float        tof_chist_MyCalib[NMyTOFTimeCalibVer][3][2];
   float        tof_chist_n_MyCalib[NMyTOFTimeCalibVer][3][2];
   float        tof_chisc_MyCalib[NMyTOFTimeCalibVer][3][2];
   float        tof_chisc_n_MyCalib[NMyTOFTimeCalibVer][3][2];
   //float        tof_betah_MyCalib[3][2];
   //float        tof_tl_MyCalib[3][2][4];
   //float        tof_tres_MyCalib[3][2][4];
   //float        tof_cres_MyCalib[3][2][4];
   //float        tof_chist_MyCalib[3][2];
   //float        tof_chist_n_MyCalib[3][2];
   //float        tof_chisc_MyCalib[3][2];
   //float        tof_chisc_n_MyCalib[3][2];
   #endif //MYTOFTIMECALIB
   //--tof using my calibration for LT, with new parameters using pr, he4, o16, pr+he4//2024.09.05
   #ifdef COMBINE_Z_TOFTIMECALIB
   float        tof_betah_CombZCalib[NTOFTimeCalibVer][2];
   float        tof_tl_CombZCalib[NTOFTimeCalibVer][2][4];
   float        tof_tres_CombZCalib[NTOFTimeCalibVer][2][4];
   float        tof_cres_CombZCalib[NTOFTimeCalibVer][2][4];
   float        tof_chist_CombZCalib[NTOFTimeCalibVer][2];
   float        tof_chist_n_CombZCalib[NTOFTimeCalibVer][2];
   float        tof_chisc_CombZCalib[NTOFTimeCalibVer][2];
   float        tof_chisc_n_CombZCalib[NTOFTimeCalibVer][2];
   #endif //COMBINE_Z_TOFTIMECALIB
   //TOF Beta without calibration, used as reference on the performance of TOF Time Calibration //2024.09.05
   #ifdef TOFBETA_WITHOUT_CAL
   float        tof_betah_NoCalib[NLTHTSHT];
   float        tof_tl_NoCalib[NLTHTSHT][4];
   float        tof_tres_NoCalib[NLTHTSHT][4];
   float        tof_cres_NoCalib[NLTHTSHT][4];
   float        tof_chist_NoCalib[NLTHTSHT];
   float        tof_chist_n_NoCalib[NLTHTSHT];
   float        tof_chisc_NoCalib[NLTHTSHT];
   float        tof_chisc_n_NoCalib[NLTHTSHT];
   #endif //TOFBETA_WITHOUT_CAL
   //2024.10.23
   #ifdef NEWTOFTIMECALIB
   float        tof_betah_NewCalib		[NZ_NewCalib][NMinChis_NewCalib][NExcludeEdge_NewCalib][3][NFix_NewCalib];
   float        tof_tl_NewCalib			[NZ_NewCalib][NMinChis_NewCalib][NExcludeEdge_NewCalib][3][NFix_NewCalib][4];
   float        tof_tres_NewCalib		[NZ_NewCalib][NMinChis_NewCalib][NExcludeEdge_NewCalib][3][NFix_NewCalib][4];
   float        tof_cres_NewCalib		[NZ_NewCalib][NMinChis_NewCalib][NExcludeEdge_NewCalib][3][NFix_NewCalib][4];
   float        tof_chist_NewCalib		[NZ_NewCalib][NMinChis_NewCalib][NExcludeEdge_NewCalib][3][NFix_NewCalib];
   float        tof_chisc_NewCalib		[NZ_NewCalib][NMinChis_NewCalib][NExcludeEdge_NewCalib][3][NFix_NewCalib];
   float        tof_chist_n_NewCalib	[NZ_NewCalib][NMinChis_NewCalib][NExcludeEdge_NewCalib][3][NFix_NewCalib];
   float        tof_chisc_n_NewCalib	[NZ_NewCalib][NMinChis_NewCalib][NExcludeEdge_NewCalib][3][NFix_NewCalib];
   #endif //NEWTOFTIMECALIB
   
   //smear betah //2024.01.03
   #ifdef _DOBETAMCTUNE_
   //float        tof_betah_Smear[NDT][2*NST+1];
   //float        tof_chist_Smear[NDT][2*NST+1];
   //float        tof_chisc_Smear[NDT][2*NST+1];
   //float        tof_chist_n_Smear[NDT][2*NST+1];
   //float        tof_chisc_n_Smear[NDT][2*NST+1];
   //float        tof_tl_Smear[NDT][2*NST+1][4];
   //float        tof_tres_Smear[NDT][2*NST+1][4];
   //float        tof_cres_Smear[NDT][2*NST+1][4];
   float        tof_betah_Smear[NDT][NST];
   float        tof_chist_Smear[NDT][NST];
   float        tof_chisc_Smear[NDT][NST];
   float        tof_chist_n_Smear[NDT][NST];
   float        tof_chisc_n_Smear[NDT][NST];
   float        tof_tl_Smear[NDT][NST][4];
   float        tof_tres_Smear[NDT][NST][4];
   float        tof_cres_Smear[NDT][NST][4];
   #endif //_DOBETAMCTUNE_

//--TrdTrack
   int          nitrdseg;
   float        distrd[3];///0x ,1y, 2theta
   float        trd_rq;
   int          trd_rz;
   float        trd_rpz;
   int          trd_pass;
   int          trd_statk;
   int          trd_nhitk;
   float        trd_lik[3];//e/p e/He p/He
   float        trd_like[3];
   int          trd_onhitk;
   float        trd_oampk;
   float        trd_qk[5];//A+U+D+Non-Dedx
   int          trd_qnhk[3];
   float        trd_ipch;
   float        trd_qrmsk[3];
   float        trd_amplk[20];//raw amps (ADC-chan)
   float        trd_pathlk[20];
   float        trd_ampsk;//tot.adc, path_length corrected
   int          trd_amphitk;//good hits for above
   float        trd_pathsk;//tot path_length

//--Rich
   int          rich_itrtrack;
   float        rich_beta[3];//getbeta+default+refit
   float        rich_ebeta;
   float        rich_ebeta1; 
   float        rich_pb;//Prob To Ring
   float        rich_udis;
   int          rich_hit;
   int          rich_used;
   int          rich_usedm;
   int          rich_stat;
   float        rich_q[2];//q2 from ring+q from ChargeR
   float        rich_width;
   float        rich_npe[3];//Collect and expect
   bool         rich_good;
   bool         rich_clean;
   bool         rich_NaF;
   int          rich_pmt;
   int          rich_pmt1;
   float        rich_PMTChargeConsistency;
   float        rich_BetaConsistency;
   float        rich_pos[3];
   float        rich_theta;
   float        rich_phi;
   float        rich_n[2];
   int          rich_tile;
   float        rich_distb;
   int          rich_cstat;
   float        rich_betap;
   float        rich_likp;
   float        rich_pbp;
   int          rich_usedp;
   float        rich_qp;
   int          rich_statp;
   float        rich_pzp[3];//3most probZ
   int          rich_rz;
   float        rich_rq;
   float        rich_rpz;

//---Ecal
   unsigned int ecal_stat;
   float        ecal_en[4];//EnergyE EnergyA EnergyD EnergyC
   float        ecal_enc[3];//GetCorrect-Energy(vdev-add)
   float        ecal_rearl;//EnergyE RearLeak
   float        ecal_bdt[3];//BDT+BDTA+BDTChis
   float        ecal_q;
   int          ecal_nhit[18];//Hit Number for each Layer
   float        ecal_el[18];//Energy for each Layer
   float        ecal_eh[18];//Max Cell Energy for each Layer

//---Anti
   int          anti_nhit;

//---TrkParticle
   int          itrtracks;
   int          tk_hitbs[2];
   float        tk_qis[2];
   float        tk_qirmss[2];
   float        tk_qlss[2][9];
   float        tk_rigiditys[5];
   float        tk_chiss[5][3];
   float        tk_pos2s[9][3];
   float        tk_dir2s[2][3];
   int          tk_ibetahs;
   float        tk_betahs;
   float        tk_tofchiscs;
   float        tk_tofchists;
   int          tk_tofsumhus;
   float        tk_tofqls[4];
};



void AMSAnalysis::BookFile(char *ofile){
  fout=new TFile(ofile,"RECREATE");

//--book tree
  char trnam[1000];
  for(int it=0;it<3;it++){
    if     (it==0)sprintf(trnam,"amstreeb");
    else if(it==1)sprintf(trnam,"amstreea");
    else if(it==2)sprintf(trnam,"amstreeu");
    tout[it]=new TTree(trnam,trnam);
    sstree[0]=AddBranch_Header(tout[it]);
    sstree[1]=AddBranch_Trigger(tout[it]);
    sstree[2]=AddBranch_RTI(tout[it]);
    sstree[3]=AddBranch_Particle(tout[it]);
    if(it==0)continue;
    sstree[4]=AddBranch_Trk(tout[it]);
    sstree[5]=AddBranch_Tof(tout[it]);
    sstree[7]=AddBranch_Trd(tout[it]);
    sstree[8]=AddBranch_Rich(tout[it]);
    sstree[9]=AddBranch_Ecal(tout[it]);
    sstree[10]=AddBranch_Anti(tout[it]);
    sstree[11]=AddBranch_TrkParticle(tout[it]);
 }

}

int AMSAnalysis::GetTreeSize(TTree *utree){
   int treelen=0;
   TObjArray* leaves=utree->GetListOfLeaves();
   int nleaves=leaves->GetEntriesFast();
   for (int i = 0; i < nleaves; i++) {
     TLeaf* leaf = (TLeaf*) leaves->UncheckedAt(i);
     int len = leaf->GetLen();
     treelen +=len;
   }
   return treelen;
}


int AMSAnalysis::AddBranch_Header(TTree *utree){//12
    int lenb=GetTreeSize(utree);
    utree->Branch("run",      &run,      "run/i");
    utree->Branch("isbadrun", &isbadrun, "isbadrun/O");
    utree->Branch("isreal",   &isreal,   "isreal/O");
    utree->Branch("event",    &event,    "event/i");
    utree->Branch("version",  &version,  "version/I");
    utree->Branch("amsdn",    &amsdn,    "amsdn/I");
    utree->Branch("errorb",   &errorb,   "errorb/i");
    utree->Branch("time",      time,     "time[2]/i");
    utree->Branch("thetas",   &thetas,   "thetas/F");
    utree->Branch("phis",     &phis,     "phis/F");
    utree->Branch("rads",     &rads,     "rads/F");
    utree->Branch("yaw",      &yaw,      "yaw/F");
    utree->Branch("pitch",    &pitch,    "pitch/F");
    utree->Branch("roll",     &roll,     "roll/F");
    int lene=GetTreeSize(utree);
    return lene-lenb;
}

int AMSAnalysis::AddBranch_Trigger(TTree *utree){//14
    int lenb=GetTreeSize(utree);
    utree->Branch("nlevel1",    &nlevel1,    "nlevel1/I");
    utree->Branch("ndaqev",     &ndaqev,     "ndaqev/I");
    utree->Branch("ndaqerr",    &ndaqerr,    "ndaqerr/I");
    utree->Branch("daqlen",     &daqlen,     "daqlen/i");
    utree->Branch("daqtyerr",   &daqtyerr,   "daqtyerr/I");
    utree->Branch("daqjrmerr",  &daqjrmerr,  "daqjrmerr/i");
    utree->Branch("daqjlen",     daqjlen,    "daqjlen[24]/i");
    utree->Branch("daqsdlen",   &daqsdlen,   "daqsdlen/i");
    utree->Branch("daqsdlenu",  &daqsdlenu,  "daqsdlenu/i");
    utree->Branch("livetime",   &livetime,   "livetime/F");
    utree->Branch("physbpatt",  &physbpatt,  "physbpatt/I");
    utree->Branch("jmembpatt",  &jmembpatt,  "jmembpatt/I");
    utree->Branch("physbpatt1", &physbpatt1, "physbpatt1/I");
    utree->Branch("jmembpatt1", &jmembpatt1, "jmembpatt1/I");
    utree->Branch("physbpatt2", &physbpatt2, "physbpatt2/I");
    utree->Branch("jmembpatt2", &jmembpatt2, "jmembpatt2/I");
    utree->Branch("ecalflag",   &ecalflag,   "ecalflag/I");
    utree->Branch("tofflag",     tofflag,    "tofflag[2]/I");
    utree->Branch("antipatt",   &antipatt,"   antipatt/I");
    int lene=GetTreeSize(utree);
    return lene-lenb;
}

int AMSAnalysis::AddBranch_RTI(TTree *utree){//22
    int lenb=GetTreeSize(utree);
    utree->Branch("irti",     &irti,     "irti/I");
    utree->Branch("zenith",   &zenith,   "zenith/F");
    utree->Branch("glong",    &glong,    "glong/F");
    utree->Branch("issaa",    &issaa,    "issaa/O");
    utree->Branch("glat",     &glat,     "glat/F");
    utree->Branch("thetam",   &thetam,   "thetam/F");
    utree->Branch("phim",     &phim,     "phim/F");
    utree->Branch("mcutoff",   mcutoff,  "mcutoff[4][2]/F");//Minus+Plus 
    utree->Branch("mcutoffi",  mcutoffi, "mcutoffi[4][2]/F");//Minus+Plus 
    utree->Branch("rtilf",    &rtilf,    "rtilf/F");
    utree->Branch("rtinev",   &rtinev,   "rtinev/F");
    utree->Branch("rtinerr",  &rtinerr,  "rtinerr/F");
    utree->Branch("rtintrig", &rtintrig, "rtintrig/F"); 
    utree->Branch("rtinpar",  &rtinpar,  "rtinpar/F"); 
    utree->Branch("rtigood",  &rtigood,  "rtigood/I");
    utree->Branch("rtinexl",   rtinexl,  "rtinexl[2][2]/F");
    utree->Branch("rticdif",   rticdif,  "rticdif[2][3]/F");
    int lene=GetTreeSize(utree); 
    return lene-lenb;
}

int AMSAnalysis::AddBranch_MC(TTree *utree){
    int lenb=GetTreeSize(utree);
    utree->Branch("mpar",     &mpar,     "mpar/I");
    utree->Branch("mmom",     &mmom,     "mmom/F");
    utree->Branch("mch",      &mch,      "mch/F");
    utree->Branch("mmass",    &mmass,    "mmass/F");
    utree->Branch("mtof_pass",&mtof_pass,"mtof_pass/I");
    utree->Branch("mevcoo",    mevcoo,   "mevcoo[3]/F");
    utree->Branch("mevdir",    mevdir,   "mevdir[3]/F");
    utree->Branch("mfcoo",     mfcoo,    "mfcoo[3]/F");
    utree->Branch("mfdir",     mfdir,    "mfdir[3]/F");
    utree->Branch("mpare",     mpare,    "mpare[2]/F");
    utree->Branch("mparp",     mparp,    "mparp[2]/I");
    utree->Branch("mparc",     mparc,    "mparc[2]/F");
    utree->Branch("msume",    &msume,    "msume/F");
    utree->Branch("msumc",    &msumc,    "msumc/F");
    utree->Branch("mtrdmom",   mtrdmom,  "mtrdmom[20]/F");
    utree->Branch("mtrdpar",   mtrdpar,  "mtrdpar[20]/I");
    utree->Branch("mtrdcoo",   mtrdcoo,  "mtrdcoo[20][3]/F");
    utree->Branch("ntrmccl",  &ntrmccl,  "ntrmccl/I");
    utree->Branch("mtrmom",    mtrmom,   "mtrmom[9]/F");
    utree->Branch("mtrpar",    mtrpar,   "mtrpar[9]/I");
    utree->Branch("mtrcoo",    mtrcoo,   "mtrcoo[9][3]/F");
    utree->Branch("mtrcoo1",   mtrcoo1,  "mtrcoo1[9][3]/F");
    utree->Branch("mtrpri",    mtrpri,   "mtrpri[9]/I");
    utree->Branch("mtrz",      mtrz,     "mtrz[9]/I");
    utree->Branch("ntofmccl", &ntofmccl, "ntofmccl/I");
    utree->Branch("mtofdep",   mtofdep,  "mtofdep[4]/F");
    utree->Branch("mtofpri",   mtofpri,  "mtofpri[4]/I");
    utree->Branch("mtofbar",   mtofbar,  "mtofbar[4]/I");
//---
    utree->Branch("mevcoo1",    mevcoo1,   "mevcoo1[21][3]/F");
    utree->Branch("mevdir1",    mevdir1,   "mevdir1[21][3]/F");
    utree->Branch("mevmom1",    mevmom1,   "mevmom1[21]/F");
//---
    int lene=GetTreeSize(utree);
    return lene-lenb;
}

int AMSAnalysis::AddBranch_Particle(TTree *utree){
    int lenb=GetTreeSize(utree);
    utree->Branch("tk_l1q",     &tk_l1q,     "tk_l1q/F");
    utree->Branch("tk_l9q",     &tk_l9q,     "tk_l9q/F");
    utree->Branch("tk_l1qxy",   tk_l1qxy,    "tk_l1qxy[2]/F");
    utree->Branch("tk_l9qxy",   tk_l9qxy,    "tk_l9qxy[2]/F");
    utree->Branch("tk_l1qs",    &tk_l1qs,    "tk_l1qs/I");
    utree->Branch("tk_l9qs",    &tk_l9qs,    "tk_l9qs/I");
    utree->Branch("ibetahs",    &ibetahs,    "ibetahs/I");
    utree->Branch("itrdtracks", &itrdtracks, "itrdtracks/I");
    utree->Branch("tof_nhits",  &tof_nhits,  "tof_nhits/I");
    utree->Branch("tof_hsumhus",&tof_hsumhus,"tof_hsumhus/I");
    utree->Branch("betahs",     &betahs,     "betahs/F");
    utree->Branch("tof_chiscs", &tof_chiscs, "tof_chiscs/F");
    utree->Branch("tof_chists", &tof_chists, "tof_chists/F");
    utree->Branch("tof_qls",     tof_qls,    "tof_qls[4]/F");
    utree->Branch("trd_thetas", &trd_thetas, "trd_thetas/F");
    utree->Branch("trd_phis",   &trd_phis,   "trd_phis/F");
    utree->Branch("trd_coos",    trd_coos,   "trd_coos[3]/F");
    utree->Branch("trd_rqs",    &trd_rqs,    "trd_rqs/F");
    utree->Branch("ecal_pos",    ecal_pos,   "ecal_pos[3]/F");
    utree->Branch("ecal_ens",   &ecal_ens,   "ecal_ens/F");
    utree->Branch("ecal_dis",   &ecal_dis,   "ecal_dis/F");
    utree->Branch("tk_l1mds",   &tk_l1mds,   "tk_l1mds/F");
    utree->Branch("tk_l9mds",   &tk_l9mds,   "tk_l9mds/F");
    utree->Branch("tk_pos1s",    tk_pos1s,   "tk_pos1s[9][3]/F");
    utree->Branch("tk_dir1s",    tk_dir1s,   "tk_dir1s[3]/F");
    utree->Branch("tk_l1qvs",   &tk_l1qvs,   "tk_l1qvs/F");
    utree->Branch("tk_l9qvs",   &tk_l9qvs,   "tk_l9qvs/F");
    utree->Branch("tk_exqvn",    tk_exqvn,   "tk_exqvn[2][2]/F");
    utree->Branch("tk_l1qvr",   &tk_l1qvr,   "tk_l1qvr/F");
    utree->Branch("tk_l9qvr",   &tk_l9qvr,   "tk_l9qvr/F");
    utree->Branch("cutoffps",   &cutoffps,   "cutoffps/F");

   //--Particle //22
    utree->Branch("nparticle",  &nparticle,  "nparticle/I");
    utree->Branch("nbetah",     &nbetah,     "nbetah/I");
    utree->Branch("nbeta",      &nbeta,      "nbeta/I");
    utree->Branch("ntrack",     &ntrack,     "ntrack/I");
    utree->Branch("ntrrawcl",   &ntrrawcl,   "ntrrawcl/I");
    utree->Branch("ntrcl",      &ntrcl,      "ntrcl/I");
    utree->Branch("ntrhit",     &ntrhit,     "ntrhit/I");
    utree->Branch("ntrdtrack",  &ntrdtrack,  "ntrdtrack/I");
    utree->Branch("nrich",      &nrich,      "nrich/I");
    utree->Branch("nrichb",     &nrichb,     "nrichb/I");
    utree->Branch("nrichh",     &nrichh,     "nrichh/I");
    utree->Branch("nshow",      &nshow,      "nshow/I");
    utree->Branch("necalhit",   &necalhit,   "necalhit/I");
    utree->Branch("ntofclh",    &ntofclh,    "ntofclh/I");
    utree->Branch("ntrdseg",    &ntrdseg,    "ntrdseg/I");
    utree->Branch("ntrdrawh",   &ntrdrawh,   "ntrdrawh/I");
    utree->Branch("ntrdcl",     &ntrdcl,     "ntrdcl/I");
    utree->Branch("btstat",     &btstat,     "btstat/I"); //Status 1:Over cutoff, 2:Under cutoff, 3:Trapped
    utree->Branch("ibeta",      &ibeta,      "ibeta/I");
    utree->Branch("ibetah",     &ibetah,     "ibetah/I");
    utree->Branch("icharge",    &icharge,    "icharge/I");
    utree->Branch("itrtrack",   &itrtrack,   "itrtrack/I");
    utree->Branch("itrdtrack",  &itrdtrack,  "itrdtrack/I");
    utree->Branch("ishow",      &ishow,      "ishow/I");
    utree->Branch("irich",      &irich,      "irich/I");
    utree->Branch("irichb",     &irichb,     "irichb/I");
    utree->Branch("iparindex",  &iparindex,  "iparindex/I");
    utree->Branch("show_dis",   &show_dis,   "show_dis/F");
    utree->Branch("ec_pos",      ec_pos,     "ec_pos[3]/F");
    utree->Branch("ec_dir",      ec_dir,     "ec_dir[3]/F");
    utree->Branch("cutoffp",     cutoffp,    "cutoffp[2]/F");
    utree->Branch("cutoffpi",    cutoffpi,   "cutoffpi[2]/F");
    //--2ndParticle //6
    utree->Branch("betah2p",    &betah2p,    "betah2p/I");
    utree->Branch("betah2q",    &betah2q,    "betah2q/F");
    utree->Branch("betah2r",    &betah2r,    "betah2r/F");
    utree->Branch("betah2hb",    betah2hb,   "betah2hb[2]/I");
    utree->Branch("betah2ch",   &betah2ch,   "betah2ch/F");
    int lene=GetTreeSize(utree);
    return lene-lenb;
}

int AMSAnalysis::AddBranch_Trk(TTree *utree){//52+380-54
    int lenb=GetTreeSize(utree);
    utree->Branch("tk_rtype",     &tk_rtype,       "tk_rtype/I");
    utree->Branch("tk_hitb",       tk_hitb,        "tk_hitb[2]/I");
    utree->Branch("tk_rz",         tk_rz,          "tk_rz[2][2]/I");//All+Inner
    utree->Branch("tk_rpz",        tk_rpz,         "tk_rpz[2][2]/F");//
    utree->Branch("tk_qr",        &tk_qr,          "tk_qr/F");
    utree->Branch("tk_q",          tk_q,           "tk_q[2]/F");//All+Inner
    utree->Branch("tk_qrms",       tk_qrms,        "tk_qrms[2]/F");
    utree->Branch("tk_qhit",       tk_qhit,        "tk_qhit[2]/I");
    utree->Branch("tk_ql",         tk_ql,          "tk_ql[9]/F");//Temp Only Layer1+2
    utree->Branch("tk_edep",       tk_edep,        "tk_edep[9]/F");
    utree->Branch("tk_ql2",        tk_ql2,         "tk_ql2[9][2]/F");
    utree->Branch("tk_qls",        tk_qls,         "tk_qls[9]/I");//TK-QL Status
    utree->Branch("tk_qlsn",      &tk_qlsn,        "tk_qlsn/I");//TK_QL Status YJ
    utree->Branch("tk_qin",        tk_qin,         "tk_qin[2][3]/F");
    utree->Branch("tk_qrmn",       tk_qrmn,        "tk_qrmn[2][3]/F");
    utree->Branch("tk_qln",        tk_qln,         "tk_qln[2][9][3]/F");
    utree->Branch("tk_lid",        tk_lid,         "tk_lid[9]/I");//TK-Ladder-Id
    utree->Branch("tk_cad",        tk_cad,         "tk_cad[9]/I");
    utree->Branch("tk_cof",        tk_cof,         "tk_cof[9][2]/F");
#ifdef _USETKADC_
    utree->Branch("tk_qlqs",       tk_qlqs,        "tk_qlqs[9][2][5]/F");
    utree->Branch("tk_stripx",     tk_stripx,      "tk_stripx[9]/I");
    utree->Branch("tk_adc",        tk_adc,         "tk_adc[9][2][11]/F");
#endif
    utree->Branch("tk_iso",        tk_iso,         "tk_iso[9][2]/F"); 
///--
    utree->Branch("tk_exql",       tk_exql,        "tk_exql[2][3]/F");
    utree->Branch("tk_exqln",      tk_exqln,       "tk_exqln[2][2][3]/F");
    utree->Branch("tk_exqls",      tk_exqls,       "tk_exqls[2]/I");//TK-QL Status
    utree->Branch("tk_exdis",      tk_exdis,       "tk_exdis[2]/F");
    utree->Branch("tk_exlid",      tk_exlid,       "tk_exlid[2]/I");
    utree->Branch("tk_excad",      tk_excad,       "tk_excad[2]/I");
///-All Span Rigidity
    utree->Branch("tk_res",        tk_res,         "tk_res[9][2]/F");
#ifdef _USELAYERRES_
    utree->Branch("tk_res1",       tk_res1,        "tk_res1[9][2]/F");
#endif
    utree->Branch("tk_fzm",        tk_fzm,         "tk_fzm[2]/F");
    utree->Branch("tk_rigidity1",  tk_rigidity1,   "tk_rigidity1[3][3][7]/F");
    utree->Branch("tk_rigidity",  &tk_rigidity,    "tk_rigidity/F");
    utree->Branch("tk_rigidityk",  tk_rigidityk,   "tk_rigidityk[4][7]/F");
    utree->Branch("tk_erigidity1", tk_erigidity1,  "tk_erigidity1[3][7]/F");
    utree->Branch("tk_chis1",      tk_chis1,       "tk_chis1[3][3][7][3]/F");
    utree->Branch("tk_dir0",       tk_dir0,        "tk_dir0[3][3][4]/F");
    utree->Branch("tk_rigidityi",  tk_rigidityi,   "tk_rigidityi[9]/F");
    utree->Branch("tk_chisi",      tk_chisi,       "tk_chisi[9]/F");
    utree->Branch("tk_bcor",      &tk_bcor,        "tk_bcor/F");
    utree->Branch("tk_cdif",       tk_cdif,        "tk_cdif[2][3]/F");
    utree->Branch("tk_pos",        tk_pos,         "tk_pos[9][3]/F");
    utree->Branch("tk_dir",        tk_dir,         "tk_dir[9][3]/F");
    utree->Branch("tk_hitc",       tk_hitc,        "tk_hitc[9][3]/F");
    utree->Branch("tk_oq",         tk_oq,          "tk_oq[2]/F");
    utree->Branch("tk_oel",        tk_oel,         "tk_oel[9]/F");
    utree->Branch("tk_ohitl",      tk_ohitl,       "tk_ohitl[9]/I");
#ifdef _USEKALMANFIT_
    utree->Branch("tk_rigidityr",  tk_rigidityr,   "tk_rigidityr[4][3][7]/F");
    utree->Branch("tk_posr",       tk_posr,        "tk_posr[4][3][7][3]/F");
    utree->Branch("tk_dirr",       tk_dirr,        "tk_dirr[4][3][7][3]/F");
#endif
    int lene=GetTreeSize(utree);
    return lene-lenb;
}

int AMSAnalysis::AddBranch_Tof(TTree *utree){//45
    int lenb=GetTreeSize(utree);
    utree->Branch("tof_btype",   &tof_btype,   "tof_btype/I");
    utree->Branch("tof_bpatt",   &tof_bpatt,   "tof_bpatt/I");
    utree->Branch("tof_beta",    &tof_beta,    "tof_beta/F");
    utree->Branch("tof_betah",   &tof_betah,   "tof_betah/F");
    utree->Branch("tof_ebetah",  &tof_ebetah,  "tof_ebetah/F");
    utree->Branch("tof_ebetah_n",&tof_ebetah_n,"tof_ebetah_n/F");
    utree->Branch("tof_betahmc", &tof_betahmc, "tof_betahmc/F");
    utree->Branch("tof_hsumh",   &tof_hsumh,   "tof_hsumh/I"); //BetaH
    utree->Branch("tof_hsumhu",  &tof_hsumhu,  "tof_hsumhu/I");//
    utree->Branch("tof_nclhl",    tof_nclhl,   "tof_nclhl[4]/I");
    utree->Branch("tof_barid",    tof_barid,   "tof_barid[4]/I");
    utree->Branch("tof_pos",      tof_pos,     "tof_pos[4][3]/F");
    utree->Branch("tof_dir",      tof_dir,     "tof_dir[4][3]/F");
    utree->Branch("tof_pass",    &tof_pass,    "tof_pass/I");
    utree->Branch("tof_oq",       tof_oq,      "tof_oq[4][2]/F");
    utree->Branch("tof_ob",       tof_ob,      "tof_ob[4][2]/I");
    utree->Branch("tof_oncl",    &tof_oncl,    "tof_oncl/I");
    utree->Branch("tof_chist",   &tof_chist,   "tof_chist/F");
    utree->Branch("tof_chisc",   &tof_chisc,   "tof_chisc/F");
    utree->Branch("tof_chist_n", &tof_chist_n, "tof_chist_n/F");
    utree->Branch("tof_chisc_n", &tof_chisc_n, "tof_chisc_n/F");
    utree->Branch("tof_z",       &tof_z,       "tof_z/I");
    utree->Branch("tof_pz",      &tof_pz,      "tof_pz/F");
    utree->Branch("tof_zr",      &tof_zr,      "tof_zr/I");
    utree->Branch("tof_pzr",     &tof_pzr,     "tof_pzr/F");
    utree->Branch("tof_zud",      tof_zud,     "tof_zud[2]/I");
    utree->Branch("tof_pud",      tof_pud,     "tof_pud[2]/F");
    utree->Branch("tof_q",        tof_q,       "tof_q[2]/F");//Max+Trun TofChargeH+BetaH
    utree->Branch("tof_qrms",     tof_qrms,    "tof_qrms[2]/F");
    utree->Branch("tof_ql",       tof_ql,      "tof_ql[4]/F");
    utree->Branch("tof_ql1",      tof_ql1,     "tof_ql1[4]/F");
#ifdef _USETKADC_
    utree->Branch("tof_qlc",      tof_qlc,     "tof_qlc[6][4]/F");
    utree->Branch("tof_qlpm",     tof_qlpm,    "tof_qlpm[4][2][4]/F");
#endif
    utree->Branch("tof_qlr",      tof_qlr,     "tof_qlr[3][4]/F");
    utree->Branch("tof_tl",       tof_tl,      "tof_tl[4]/F");
    utree->Branch("tof_etl",      tof_etl,     "tof_etl[4]/F");
    utree->Branch("tof_tres",     tof_tres,    "tof_tres[4]/F"); //YC //2024.01.03
    utree->Branch("tof_leng",     tof_leng,    "tof_leng[4]/F"); //YC //2024.06.02
    utree->Branch("tof_cres",     tof_cres,    "tof_cres[4]/F");
    utree->Branch("tof_sbit",     tof_sbit,    "tof_sbit[4][2]/i");
    utree->Branch("tof_qs"  ,    &tof_qs,      "tof_qs/I");//Q-Status
	//2024.06.27
    utree->Branch("tof_betah_HT",	&tof_betah_HT,		"tof_betah_HT/F");
    utree->Branch("tof_betah_SHT",	&tof_betah_SHT,		"tof_betah_SHT/F");
    utree->Branch("tof_tl_HT",       tof_tl_HT,      "tof_tl_HT[4]/F"); //2024.07.02
    utree->Branch("tof_tl_SHT",       tof_tl_SHT,      "tof_tl_SHT[4]/F"); //2024.07.02
    utree->Branch("tof_betah_S2S3",	&tof_betah_S2S3,	"tof_betah_S2S3/F");
    utree->Branch("tof_betah_S1S4",	&tof_betah_S1S4,	"tof_betah_S1S4/F");
    utree->Branch("tof_betah_S1S2S3S4",	&tof_betah_S1S2S3S4,	"tof_betah_S1S2S3S4/F");
	//add S1S2 and S3S4 as well, the difference due to energy loss should be largest for them //2024.07.02
    utree->Branch("tof_betah_S1S2",	&tof_betah_S1S2,	"tof_betah_S1S2/F");
    utree->Branch("tof_betah_S3S4",	&tof_betah_S3S4,	"tof_betah_S3S4/F");
    utree->Branch("tof_edep",		tof_edep,			"tof_edep[4]/F");
	//2024.07.02
    /*utree->Branch("tof_deltaBetah",       tof_deltaBetah,      "tof_deltaBetah[4]/F");
    utree->Branch("tof_betah_iter",      &tof_betah_iter,    "tof_betah_iter/F");
    utree->Branch("tof_deltaBetah_HT",    tof_deltaBetah_HT,   "tof_deltaBetah_HT[4]/F");
    utree->Branch("tof_betah_iter_HT",   &tof_betah_iter_HT, "tof_betah_iter_HT/F");
    utree->Branch("tof_deltaBetah_SHT",   tof_deltaBetah_SHT,  "tof_deltaBetah_SHT[4]/F");
    utree->Branch("tof_betah_iter_SHT",  &tof_betah_iter_SHT, "tof_betah_iter_SHT/F");*/
    utree->Branch("tof_deltaBetah",       tof_deltaBetah,      "tof_deltaBetah[10][4]/F");
    utree->Branch("tof_betah_iter",      tof_betah_iter,     "tof_betah_iter[10]/F");
    utree->Branch("tof_deltaBetah_HT",    tof_deltaBetah_HT,   "tof_deltaBetah_HT[10][4]/F");
    utree->Branch("tof_betah_iter_HT",   tof_betah_iter_HT,  "tof_betah_iter_HT[10]/F");
    utree->Branch("tof_deltaBetah_SHT",   tof_deltaBetah_SHT,  "tof_deltaBetah_SHT[10][4]/F");
    utree->Branch("tof_betah_iter_SHT",  tof_betah_iter_SHT, "tof_betah_iter_SHT[10]/F");
	//--calculate the beta change from energy loss directly, using average mass, and then refit //2024.07.22
    utree->Branch("tof_betah_el",		&tof_betah_el,		"tof_betah_el/F");
    utree->Branch("tof_deltaBetah_el",	tof_deltaBetah_el,	"tof_deltaBetah_el[4]/F");
    utree->Branch("tof_deltaT_el",	tof_deltaT_el,	"tof_deltaT_el[4]/F");
	//store the residual as well //2024.07.29
    utree->Branch("tof_tres_el",  tof_tres_el, "tof_tres_el[4]/F");
    utree->Branch("tof_cres_el",  tof_cres_el, "tof_cres_el[4]/F");
	//store the residual, chis, for HT/SHT as well //2024.07.04
    utree->Branch("tof_tres_HT",  tof_tres_HT, "tof_tres_HT[4]/F");
    utree->Branch("tof_cres_HT",  tof_cres_HT, "tof_cres_HT[4]/F");
    utree->Branch("tof_chist_HT",  &tof_chist_HT, "tof_chist_HT/F");
    utree->Branch("tof_chist_n_HT",  &tof_chist_n_HT, "tof_chist_n_HT/F");
    utree->Branch("tof_chisc_HT",  &tof_chisc_HT, "tof_chisc_HT/F");
    utree->Branch("tof_chisc_n_HT",  &tof_chisc_n_HT, "tof_chisc_n_HT/F");
    utree->Branch("tof_tres_SHT",  tof_tres_SHT, "tof_tres_SHT[4]/F");
    utree->Branch("tof_cres_SHT",  tof_cres_SHT, "tof_cres_SHT[4]/F");
    utree->Branch("tof_chist_SHT",  &tof_chist_SHT, "tof_chist_SHT/F");
    utree->Branch("tof_chist_n_SHT",  &tof_chist_n_SHT, "tof_chist_n_SHT/F");
    utree->Branch("tof_chisc_SHT",  &tof_chisc_SHT, "tof_chisc_SHT/F");
    utree->Branch("tof_chisc_n_SHT",  &tof_chisc_n_SHT, "tof_chisc_n_SHT/F");
	//combine beta fit using LT, HT and SHT //2024.07.28
    utree->Branch("tof_betah_LTHTSHT",	&tof_betah_LTHTSHT,		"tof_betah_LTHTSHT/F");
	//--store the nHT and nSHT //2024.09.16
    utree->Branch("nHTime",		nHTime,		"nHTime[4][2]/I");
    utree->Branch("nSHTime",	nSHTime,	"nSHTime[4][2]/I");
    utree->Branch("tof_rtime",	tof_rtime,	"tof_rtime[4][2]/F");
    utree->Branch("tof_rtime_ht",	tof_rtime_ht,	"tof_rtime_ht[4][2]/F");
    utree->Branch("tof_rtime_sht",	tof_rtime_sht,	"tof_rtime_sht[4][2]/F");
	//--store Anode adc //2024.09.19
    utree->Branch("tof_Aadc",	tof_Aadc,	"tof_Aadc[4][2]/F");
	//--tof using my calibration //2024.08.23
	#ifdef MYTOFTIMECALIB
    utree->Branch("tof_betah_MyCalib",	tof_betah_MyCalib,		Form("tof_betah_MyCalib[%d][3][2]/F", NMyTOFTimeCalibVer));
    utree->Branch("tof_tl_MyCalib",		tof_tl_MyCalib,			Form("tof_tl_MyCalib[%d][3][2][4]/F", NMyTOFTimeCalibVer));
    utree->Branch("tof_tres_MyCalib",	tof_tres_MyCalib,		Form("tof_tres_MyCalib[%d][3][2][4]/F", NMyTOFTimeCalibVer));
    utree->Branch("tof_cres_MyCalib",	tof_cres_MyCalib,		Form("tof_cres_MyCalib[%d][3][2][4]/F", NMyTOFTimeCalibVer));
    utree->Branch("tof_chist_MyCalib",	tof_chist_MyCalib,		Form("tof_chist_MyCalib[%d][3][2]/F", NMyTOFTimeCalibVer));
    utree->Branch("tof_chist_n_MyCalib",tof_chist_n_MyCalib,	Form("tof_chist_n_MyCalib[%d][3][2]/F", NMyTOFTimeCalibVer));
    utree->Branch("tof_chisc_MyCalib",	tof_chisc_MyCalib,		Form("tof_chisc_MyCalib[%d][3][2]/F", NMyTOFTimeCalibVer));
    utree->Branch("tof_chisc_n_MyCalib",tof_chisc_n_MyCalib,	Form("tof_chisc_n_MyCalib[%d][3][2]/F", NMyTOFTimeCalibVer));
	
    //utree->Branch("tof_betah_MyCalib",	tof_betah_MyCalib,		"tof_betah_MyCalib[3][2]/F");
    //utree->Branch("tof_tl_MyCalib",		tof_tl_MyCalib,			"tof_tl_MyCalib[3][2][4]/F");
    //utree->Branch("tof_tres_MyCalib",	tof_tres_MyCalib,		"tof_tres_MyCalib[3][2][4]/F");
    //utree->Branch("tof_cres_MyCalib",	tof_cres_MyCalib,		"tof_cres_MyCalib[3][2][4]/F");
    //utree->Branch("tof_chist_MyCalib",	tof_chist_MyCalib,		"tof_chist_MyCalib[3][2]/F");
    //utree->Branch("tof_chist_n_MyCalib",tof_chist_n_MyCalib,	"tof_chist_n_MyCalib[3][2]/F");
    //utree->Branch("tof_chisc_MyCalib",	tof_chisc_MyCalib,		"tof_chisc_MyCalib[3][2]/F");
    //utree->Branch("tof_chisc_n_MyCalib",tof_chisc_n_MyCalib,	"tof_chisc_n_MyCalib[3][2]/F");
	#endif //MYTOFTIMECALIB
	//--tof using my calibration for LT, with new parameters using pr, he4, o16, pr+he4//2024.09.05
	#ifdef COMBINE_Z_TOFTIMECALIB
    //utree->Branch("tof_betah_CombZCalib",	tof_betah_CombZCalib,	"tof_betah_CombZCalib[4][2]/F");
    //utree->Branch("tof_tl_CombZCalib",		tof_tl_CombZCalib,		"tof_tl_CombZCalib[4][2][4]/F");
    //utree->Branch("tof_tres_CombZCalib",	tof_tres_CombZCalib,	"tof_tres_CombZCalib[4][2][4]/F");
    //utree->Branch("tof_cres_CombZCalib",	tof_cres_CombZCalib,	"tof_cres_CombZCalib[4][2][4]/F");
    //utree->Branch("tof_chist_CombZCalib",	tof_chist_CombZCalib,	"tof_chist_CombZCalib[4][2]/F");
    //utree->Branch("tof_chist_n_CombZCalib",	tof_chist_n_CombZCalib,	"tof_chist_n_CombZCalib[4][2]/F");
    //utree->Branch("tof_chisc_CombZCalib",	tof_chisc_CombZCalib,	"tof_chisc_CombZCalib[4][2]/F");
    //utree->Branch("tof_chisc_n_CombZCalib",	tof_chisc_n_CombZCalib,	"tof_chisc_n_CombZCalib[4][2]/F");
    utree->Branch("tof_betah_CombZCalib",	tof_betah_CombZCalib,	Form("tof_betah_CombZCalib[%d][2]/F", NTOFTimeCalibVer));
    utree->Branch("tof_tl_CombZCalib",		tof_tl_CombZCalib,		Form("tof_tl_CombZCalib[%d][2][4]/F", NTOFTimeCalibVer));
    utree->Branch("tof_tres_CombZCalib",	tof_tres_CombZCalib,	Form("tof_tres_CombZCalib[%d][2][4]/F", NTOFTimeCalibVer));
    utree->Branch("tof_cres_CombZCalib",	tof_cres_CombZCalib,	Form("tof_cres_CombZCalib[%d][2][4]/F", NTOFTimeCalibVer));
    utree->Branch("tof_chist_CombZCalib",	tof_chist_CombZCalib,	Form("tof_chist_CombZCalib[%d][2]/F", NTOFTimeCalibVer));
    utree->Branch("tof_chist_n_CombZCalib",	tof_chist_n_CombZCalib,	Form("tof_chist_n_CombZCalib[%d][2]/F", NTOFTimeCalibVer));
    utree->Branch("tof_chisc_CombZCalib",	tof_chisc_CombZCalib,	Form("tof_chisc_CombZCalib[%d][2]/F", NTOFTimeCalibVer));
    utree->Branch("tof_chisc_n_CombZCalib",	tof_chisc_n_CombZCalib,	Form("tof_chisc_n_CombZCalib[%d][2]/F", NTOFTimeCalibVer));
	#endif //COMBINE_Z_TOFTIMECALIB

	//TOF Beta without calibration, used as reference on the performance of TOF Time Calibration //2024.09.05
	#ifdef TOFBETA_WITHOUT_CAL
    utree->Branch("tof_betah_NoCalib",		tof_betah_NoCalib,		Form("tof_betah_NoCalib[%d]/F", NLTHTSHT));
    utree->Branch("tof_tl_NoCalib",			tof_tl_NoCalib,			Form("tof_tl_NoCalib[%d][4]/F", NLTHTSHT));
    utree->Branch("tof_tres_NoCalib",		tof_tres_NoCalib,		Form("tof_tres_NoCalib[%d][4]/F", NLTHTSHT));
    utree->Branch("tof_cres_NoCalib",		tof_cres_NoCalib,		Form("tof_cres_NoCalib[%d][4]/F", NLTHTSHT));
    utree->Branch("tof_chist_NoCalib",		tof_chist_NoCalib,		Form("tof_chist_NoCalib[%d]/F", NLTHTSHT));
    utree->Branch("tof_chist_n_NoCalib",	tof_chist_n_NoCalib,	Form("tof_chist_n_NoCalib[%d]/F", NLTHTSHT));
    utree->Branch("tof_chisc_NoCalib",		tof_chisc_NoCalib,		Form("tof_chisc_NoCalib[%d]/F", NLTHTSHT));
    utree->Branch("tof_chisc_n_NoCalib",	tof_chisc_n_NoCalib,	Form("tof_chisc_n_NoCalib[%d]/F", NLTHTSHT));
	#endif //TOFBETA_WITHOUT_CAL
	
	//2024.10.23
	#ifdef NEWTOFTIMECALIB
    utree->Branch("tof_betah_NewCalib",		tof_betah_NewCalib,		Form("tof_betah_NewCalib[%d][%d][%d][3][%d]/F", NZ_NewCalib, NMinChis_NewCalib, NExcludeEdge_NewCalib, NFix_NewCalib));
    utree->Branch("tof_tl_NewCalib",		tof_tl_NewCalib,		Form("tof_tl_NewCalib[%d][%d][%d][3][%d][4]/F", NZ_NewCalib, NMinChis_NewCalib, NExcludeEdge_NewCalib, NFix_NewCalib));
    utree->Branch("tof_tres_NewCalib",		tof_tres_NewCalib,		Form("tof_tres_NewCalib[%d][%d][%d][3][%d][4]/F", NZ_NewCalib, NMinChis_NewCalib, NExcludeEdge_NewCalib, NFix_NewCalib));
    utree->Branch("tof_cres_NewCalib",		tof_cres_NewCalib,		Form("tof_cres_NewCalib[%d][%d][%d][3][%d][4]/F", NZ_NewCalib, NMinChis_NewCalib, NExcludeEdge_NewCalib, NFix_NewCalib));
    utree->Branch("tof_chist_NewCalib",		tof_chist_NewCalib,		Form("tof_chist_NewCalib[%d][%d][%d][3][%d]/F", NZ_NewCalib, NMinChis_NewCalib, NExcludeEdge_NewCalib, NFix_NewCalib));
    utree->Branch("tof_chisc_NewCalib",		tof_chisc_NewCalib,		Form("tof_chisc_NewCalib[%d][%d][%d][3][%d]/F", NZ_NewCalib, NMinChis_NewCalib, NExcludeEdge_NewCalib, NFix_NewCalib));
    utree->Branch("tof_chist_n_NewCalib",	tof_chist_n_NewCalib,	Form("tof_chist_n_NewCalib[%d][%d][%d][3][%d]/F", NZ_NewCalib, NMinChis_NewCalib, NExcludeEdge_NewCalib, NFix_NewCalib));
    utree->Branch("tof_chisc_n_NewCalib",	tof_chisc_n_NewCalib,	Form("tof_chisc_n_NewCalib[%d][%d][%d][3][%d]/F", NZ_NewCalib, NMinChis_NewCalib, NExcludeEdge_NewCalib, NFix_NewCalib));
	#endif //NEWTOFTIMECALIB
	
	//smear betah //2024.01.03
    #ifdef _DOBETAMCTUNE_
    //utree->Branch("tof_betah_Smear",   tof_betah_Smear,   Form("tof_betah_Smear[%d][%d]/F", NDT, 2*NST+1));
    //utree->Branch("tof_chist_Smear",   tof_chist_Smear,   Form("tof_chist_Smear[%d][%d]/F", NDT, 2*NST+1));
    //utree->Branch("tof_chisc_Smear",   tof_chisc_Smear,   Form("tof_chisc_Smear[%d][%d]/F", NDT, 2*NST+1));
    //utree->Branch("tof_chist_n_Smear",   tof_chist_n_Smear,   Form("tof_chist_n_Smear[%d][%d]/F", NDT, 2*NST+1));
    //utree->Branch("tof_chisc_n_Smear",   tof_chisc_n_Smear,   Form("tof_chisc_n_Smear[%d][%d]/F", NDT, 2*NST+1));
    //utree->Branch("tof_tl_Smear",   tof_tl_Smear,   Form("tof_tl_Smear[%d][%d][4]/F", NDT, 2*NST+1));
    //utree->Branch("tof_tres_Smear",   tof_tres_Smear,   Form("tof_tres_Smear[%d][%d][4]/F", NDT, 2*NST+1));
    //utree->Branch("tof_cres_Smear",   tof_cres_Smear,   Form("tof_cres_Smear[%d][%d][4]/F", NDT, 2*NST+1));
    utree->Branch("tof_betah_Smear",   tof_betah_Smear,   Form("tof_betah_Smear[%d][%d]/F", NDT, NST));
    utree->Branch("tof_chist_Smear",   tof_chist_Smear,   Form("tof_chist_Smear[%d][%d]/F", NDT, NST));
    utree->Branch("tof_chisc_Smear",   tof_chisc_Smear,   Form("tof_chisc_Smear[%d][%d]/F", NDT, NST));
    utree->Branch("tof_chist_n_Smear",   tof_chist_n_Smear,   Form("tof_chist_n_Smear[%d][%d]/F", NDT, NST));
    utree->Branch("tof_chisc_n_Smear",   tof_chisc_n_Smear,   Form("tof_chisc_n_Smear[%d][%d]/F", NDT, NST));
    utree->Branch("tof_tl_Smear",   tof_tl_Smear,   Form("tof_tl_Smear[%d][%d][4]/F", NDT, NST));
    utree->Branch("tof_tres_Smear",   tof_tres_Smear,   Form("tof_tres_Smear[%d][%d][4]/F", NDT, NST));
    utree->Branch("tof_cres_Smear",   tof_cres_Smear,   Form("tof_cres_Smear[%d][%d][4]/F", NDT, NST));
    #endif //_DOBETAMCTUNE_
    int lene=GetTreeSize(utree);
    return lene-lenb;
}

int AMSAnalysis::AddBranch_Trd(TTree *utree){//25
    int lenb=GetTreeSize(utree);
    utree->Branch("nitrdseg",    &nitrdseg,   "nitrdseg/I");
    utree->Branch("distrd",       distrd,     "distrd[3]/F");
    utree->Branch("trd_pass",    &trd_pass,   "trd_pass/I");
    utree->Branch("trd_statk",   &trd_statk,  "trd_statk/I");
    utree->Branch("trd_nhitk",   &trd_nhitk,  "trd_nhitk/I");
    utree->Branch("trd_lik",      trd_lik,    "trd_lik[3]/F");
    utree->Branch("trd_like",     trd_like,   "trd_like[3]/F");
    utree->Branch("trd_onhitk",  &trd_onhitk, "trd_onhitk/I");
    utree->Branch("trd_oampk",   &trd_oampk,  "trd_oampk/F");
    utree->Branch("trd_qk",       trd_qk,     "trd_qk[5]/F");//All+Up+Down
    utree->Branch("trd_qnhk",     trd_qnhk,   "trd_qnhk[3]/I");
    utree->Branch("trd_ipch",    &trd_ipch,   "trd_ipch/F");
    utree->Branch("trd_qrmsk",    trd_qrmsk,  "trd_qrmsk[3]/F");
    utree->Branch("trd_rq",      &trd_rq,     "trd_rq/F");
    utree->Branch("trd_rz",      &trd_rz,     "trd_rz/I");
    utree->Branch("trd_rpz",     &trd_rpz,    "trd_rpz/F");
    utree->Branch("trd_amplk",    trd_amplk,  "trd_amplk[20]/F");
    utree->Branch("trd_pathlk",   trd_pathlk, "trd_pathlk[20]/F");
    utree->Branch("trd_ampsk",   &trd_ampsk,  "trd_ampsk/F");
    utree->Branch("trd_pathsk",  &trd_pathsk, "trd_pathsk/F");
    utree->Branch("trd_amphitk", &trd_amphitk,"trd_amphitk/I");
    int lene=GetTreeSize(utree);
    return lene-lenb;
}

int AMSAnalysis::AddBranch_Rich(TTree *utree){//25
    int lenb=GetTreeSize(utree);
    utree->Branch("rich_itrtrack",&rich_itrtrack,"rich_itrtrack/I");
    utree->Branch("rich_beta",     rich_beta,  "rich_beta[3]/F");
    utree->Branch("rich_ebeta",   &rich_ebeta, "rich_ebeta/F");
    utree->Branch("rich_ebeta1",   &rich_ebeta1, "rich_ebeta1/F");
    utree->Branch("rich_pb",      &rich_pb,    "rich_pb/F");
    utree->Branch("rich_udis",    &rich_udis,  "rich_udis/F");
    utree->Branch("rich_hit",     &rich_hit,  "rich_hit/I");
    utree->Branch("rich_used",    &rich_used,  "rich_used/I");
    utree->Branch("rich_usedm",   &rich_usedm, "rich_usedm/I");
    utree->Branch("rich_stat",    &rich_stat,  "rich_stat/I");
    utree->Branch("rich_q",        rich_q,     "rich_q[2]/F");
    utree->Branch("rich_width",   &rich_width, "rich_width/F");
    utree->Branch("rich_npe",      rich_npe,   "rich_npe[3]/F");
    utree->Branch("rich_good",    &rich_good,  "rich_good/O");
    utree->Branch("rich_clean",   &rich_clean, "rich_clean/O");
    utree->Branch("rich_NaF",     &rich_NaF,   "rich_NaF/O");
    utree->Branch("rich_pmt",     &rich_pmt,    "rich_pmt/I");
    utree->Branch("rich_pmt1",    &rich_pmt1,   "rich_pmt1/I");
    utree->Branch("rich_PMTChargeConsistency",&rich_PMTChargeConsistency,"rich_PMTChargeConsistency/F");
    utree->Branch("rich_BetaConsistency",&rich_BetaConsistency,"rich_BetaConsistency/F");
    utree->Branch("rich_pos",      rich_pos,    "rich_pos[3]/F");
    utree->Branch("rich_theta",   &rich_theta,  "rich_theta/F");
    utree->Branch("rich_phi",     &rich_phi,    "rich_phi/F");
    utree->Branch("rich_n",        rich_n,      "rich_n[2]/F");
    utree->Branch("rich_tile",    &rich_tile,   "rich_tile/I");
    utree->Branch("rich_distb",   &rich_distb,  "rich_distb/F");
    utree->Branch("rich_cstat",   &rich_cstat,  "rich_cstat/I");
//--Lip
    utree->Branch("rich_betap",   &rich_betap, "rich_betap/F");
    utree->Branch("rich_likp",    &rich_likp,  "rich_likp/F");
    utree->Branch("rich_pbp",     &rich_pbp,   "rich_pbp/F");
    utree->Branch("rich_usedp",   &rich_usedp, "rich_usedp/I");
    utree->Branch("rich_statp",   &rich_statp, "rich_statp/I");
    utree->Branch("rich_qp",      &rich_qp,    "rich_qp/F");
    utree->Branch("rich_pzp",      rich_pzp,   "rich_pzp[3]/F");
//--RICHZ
    utree->Branch("rich_rz",      &rich_rz,    "rich_rz/I");
    utree->Branch("rich_rq",      &rich_rq,    "rich_rq/F");
    utree->Branch("rich_rpz",     &rich_rpz,   "rich_rpz/F");
    int lene=GetTreeSize(utree);
    return lene-lenb;
}

int AMSAnalysis::AddBranch_Ecal(TTree *utree){
   int lenb=GetTreeSize(utree);
   utree->Branch("ecal_stat",   &ecal_stat,    "ecal_stat/i"); 
   utree->Branch("ecal_en",      ecal_en,      "ecal_en[4]/F"); //EnergyE EnergyA EnergyD
   utree->Branch("ecal_enc",     ecal_enc,     "ecal_enc[3]/F");//GetCorrect-Energy(vdev-add)
   utree->Branch("ecal_rearl",  &ecal_rearl,   "ecal_rearl/F");
   utree->Branch("ecal_bdt",     ecal_bdt,     "ecal_bdt[3]/F");
   utree->Branch("ecal_q",      &ecal_q,       "ecal_q/F");
   utree->Branch("ecal_nhit",    ecal_nhit,    "ecal_nhit[18]/I");
   utree->Branch("ecal_el",      ecal_el,      "ecal_el[18]/F");
   utree->Branch("ecal_eh",      ecal_eh,      "ecal_eh[18]/F");
   int lene=GetTreeSize(utree);
   return lene-lenb;
}

int AMSAnalysis::AddBranch_Anti(TTree *utree){
   int lenb=GetTreeSize(utree);
   utree->Branch("anti_nhit",   &anti_nhit,    "anti_nhit/I");
   int lene=GetTreeSize(utree);
   return lene-lenb;
}

int AMSAnalysis::AddBranch_TrkParticle(TTree *utree){
   int lenb=GetTreeSize(utree);
   utree->Branch("itrtracks",   &itrtracks,    "itrtracks/I");
   utree->Branch("tk_hitbs",     tk_hitbs,     "tk_hitbs[2]/I");
   utree->Branch("tk_qis",       tk_qis,       "tk_qis[2]/F");//All+Inner
   utree->Branch("tk_qirmss",    tk_qirmss,    "tk_qirmss[2]/F");
   utree->Branch("tk_qlss",      tk_qlss,      "tk_qlss[2][9]/F");
   utree->Branch("tk_rigiditys", tk_rigiditys, "tk_rigiditys[5]/F");
   utree->Branch("tk_chiss",     tk_chiss,     "tk_chiss[5][3]/F");
   utree->Branch("tk_pos2s",     tk_pos2s,     "tk_pos2s[9][3]/F");
   utree->Branch("tk_dir2s",     tk_dir2s,     "tk_dir2s[2][3]/F");
   utree->Branch("tk_ibetahs",  &tk_ibetahs,   "tk_ibetahs/I");
   utree->Branch("tk_betahs",   &tk_betahs,    "tk_betahs/F");
   utree->Branch("tk_tofchiscs",&tk_tofchiscs, "tk_tofchiscs/F");
   utree->Branch("tk_tofchists",&tk_tofchists, "tk_tofchists/F");
   utree->Branch("tk_tofsumhus",&tk_tofsumhus, "tk_tofsumhus/I");
   utree->Branch("tk_tofqls",    tk_tofqls,    "tk_tofqls[4]/F");
   int lene=GetTreeSize(utree);
   return lene-lenb;
}


int AMSAnalysis::LoadHeader(){

///--RTI
     AMSSetupR::RTI::UseLatest();

///--Track
     TkDBc::UseFinal();
     TRMCFFKEY_DEF::ReadFromFile = 0;
     TRFITFFKEY_DEF::ReadFromFile = 0;
     TrInnerDzDB::ForceFromTDV = 1;
     cout<<"NewN TkVersion Used"<<endl;

    return 0;
}


int AMSAnalysis::LoadHeaderEv(){

//---Center of Gravity
    TRCLFFKEY.ClusterCofGOpt=1;
    if(isreal){TrLinearEtaDB::SetLinearCluster();}
    return 0;
}


int  AMSAnalysis::LoadHeaderMC(int mcz,int vers,int mparid){

     TRFITFFKEY.magtemp = 0;
     TRFITFFKEY.Zshift=-1;
//---Load default tuning parameter
     TRMCFFKEY.MCtuneDmax = 100e-4;
     TRMCFFKEY.MCtuneDs[0]= 0;//InnerX
     TRMCFFKEY.MCtuneDs[1]= -1.000001;//InnerY
     TRMCFFKEY.MCtuneDy9  = 0;
/*     TRMCFFKEY.OuterSmearing[0][1] = -7.5e-4;//L1Y
     TRMCFFKEY.OuterSmearing[1][1] = -8.3e-4;//L9Y
     TRMCFFKEY.OuterSmearing[0][0] = 7.5e-4;//L1X
     TRMCFFKEY.OuterSmearing[1][0] = 8.5e-4;//L9X*/
     TRMCFFKEY.OuterSmearing[0][1] = -11.4e-4;//L1Y
     TRMCFFKEY.OuterSmearing[1][1] = -12.6e-4;//L9Y
     TRMCFFKEY.OuterSmearing[0][0] = 11.4e-4;//L1X
     TRMCFFKEY.OuterSmearing[1][0] = 12.6e-4;//L9X
/*     TRMCFFKEY.OuterSmearing[0][1] = -14.5e-4;//L1Y
     TRMCFFKEY.OuterSmearing[1][1] = -16.0e-4;//L9Y
     TRMCFFKEY.OuterSmearing[0][0] = 14.5e-4;//L1X
     TRMCFFKEY.OuterSmearing[1][0] = 16.0e-4;//L9X*/
//---V6 external layer tunning
     TRMCFFKEY.OuterSmearingN[0][1] = 7.1e-4;//L1Y
     TRMCFFKEY.OuterSmearingN[1][1] = 7.9e-4;//L9Y
     TRMCFFKEY.OuterSmearingN[0][0] = 8.e-4;//L1X
     TRMCFFKEY.OuterSmearingN[1][0] = 8.e-4;//L9X
     TRMCFFKEY.MCscat[0] = 0; 
     TRMCFFKEY.MCscat[1] = 0; 
     TRMCFFKEY.MCscat[2] = 0;
#ifdef _USEMCTKRAW_
     TRMCFFKEY.MCtuneDs[1]= -1.000001;
#else
     if (abs(mcz)==9){//F
        TRMCFFKEY.MCtuneDs[1]= -900.;
     }
     else if (abs(mcz)==10){//Neon
        TRMCFFKEY.MCtuneDs[1]= -1000.;
     }
     else if (abs(mcz)==11){//Sodium
        TRMCFFKEY.MCtuneDs[1]= -1100.;
     }
     else if (abs(mcz)==12){//Magnesium
        TRMCFFKEY.MCtuneDs[1]= -1200.;
     }
     else if (abs(mcz)==13){//Al
        TRMCFFKEY.MCtuneDs[1]= -1300.;
     }
     else if (abs(mcz)==14){//Silicon
        TRMCFFKEY.MCtuneDs[1]= -1400.;
     }
     else if (abs(mcz)==16){//Sulfur
        TRMCFFKEY.MCtuneDs[1]= -1600.;
     }
     else if (abs(mcz)==26){//Ferrum
        TRMCFFKEY.MCtuneDs[1]= -2600.;
     }
     else if (abs(mcz)==8){//Oxygen
//       if(vers>=1106)TRMCFFKEY.MCtuneDs[1]= -1.000001;
       if     (vers>=1200)TRMCFFKEY.MCtuneDs[1]= -820.;
       else               TRMCFFKEY.MCtuneDs[1]= -810.;
     }
     else if (abs(mcz)==7){//Nitrogen
//       if(vers>=1106)TRMCFFKEY.MCtuneDs[1]= -1.000001;
       TRMCFFKEY.MCtuneDs[1]= -710.;
     }
     else if (abs(mcz)==6){//Carbon
//       if     (vers>=1106)TRMCFFKEY.MCtuneDs[1]= -1.000001;
       if     (vers>=1200)TRMCFFKEY.MCtuneDs[1]= -680.;
       else               TRMCFFKEY.MCtuneDs[1]= -670.;
     }
     else if (abs(mcz)==5){//Boron
//       if     (vers>=1106)TRMCFFKEY.MCtuneDs[1]= -1.000001;
       TRMCFFKEY.MCtuneDs[1]= -550.;
     }
     else if (abs(mcz)==4){//Beryllium
       TRMCFFKEY.MCtuneDs[1]= -410.;
     }
     else if (abs(mcz)==3){//Lithium
//       if(vers>=1106)TRMCFFKEY.MCtuneDs[1]= -1.000001;
       TRMCFFKEY.MCtuneDs[1]= -340.; 
     }
     else if (abs(mcz)==2){//Helium
//       if     (vers>=1106)TRMCFFKEY.MCtuneDs[1]= -1.000001;
       TRMCFFKEY.MCtuneDs[1]= -230.;
     }
     else if (abs(mcz)==1){
       bool ismcproton=(abs(mparid)==14||abs(mparid)==15);
       if(vers>=1106){
//         TRMCFFKEY.MCtuneDs[1]= -1.000001;
         TRMCFFKEY.MCtuneDs[1]= -110.;
       }
       else if(vers>=928&&ismcproton){
         TRMCFFKEY.MCtuneDs[0]= -9.0e-4;
         TRMCFFKEY.MCtuneDs[1]=  101.0;
       }
       if(vers<900&&ismcproton){//Proton
         TRMCFFKEY.MCtuneDy9  = 1.0e-4;  //<B900
         TRMCFFKEY.MCscat[0] = -15;      // This is NEW
         TRMCFFKEY.MCscat[1] = -440.06;  // This is NEW
         TRMCFFKEY.MCscat[2] = -440.05;  // This is NEW
       }
    }
#endif
//---TRDK
  TrdKCluster::ForceReadAlignment=0;
  TrdKCluster::ForceReadCalibration=0;
  TrdKCluster::ForceReadXePressure=0;
  TrdKCluster::SetDefaultMCXePressure(900);
  return 0;
}

int AMSAnalysis::Select(){

//---IsReal
   isreal=(pev->nMCEventg()==0);
   static bool hasev=0;
   if(hasev==0){
      for(int it=0;it<3;it++)if(isreal==0){sstree[2]+=AddBranch_MC(tout[it]);}
      hasev=1;
   }

//----
   clock_t tb,te;
   sstat[0]++;

//---Datacard
   LoadHeaderEv();

//---Header 0
   tb=clock();
   bool selhead=Select_Header();
   te=clock();
   cput[0]+=te-tb;
   if(!selhead)return 0;
   sstat[1]++;

//---Trigger 1
   tb=te;
   bool seltrig=Select_Trigger();
   te=clock();
   cput[1]+=(te-tb);
   if(!seltrig)return 0;
   sstat[2]++;

//----RTI or MC
   tb=te;
   bool selrti=isreal?Select_RTI():Select_MC(0);
   if(isreal==0)LoadHeaderMC(int(mch),version,mpar);
   te=clock();
   cput[2]+=(te-tb);
   if(!selrti)return 0;
   sstat[3]++;

//---Particle 2
   tb=te;
   bool selpart=Select_Particle();
   te=clock();
   cput[3]+=(te-tb); 
   if(!selpart)return 0;
   sstat[4]++;
 
//----MC2
   if(isreal==0){
     tb=te;
     selrti=Select_MC(1);
     te=clock();
     cput[2]+=(te-tb);
//     if(!selrti)return 0;
   } 
 
//--Preselection Track 3
   tb=te;
   bool seltrk1=Select_Trk(0);
   te=clock();
   cput[4]+=(te-tb);
//   if(!seltrk1)return 1;
   sstat[5]++;

//--Strict Selection Tof 4
   tb=te;
   bool seltof=Select_Tof();
   te=clock();
   cput[5]+=(te-tb);
//   if(!seltof)return 1;
   sstat[6]++;

//---Track 5
   tb=te;
   bool seltrk0=Select_Trk(1);
   te=clock();
   cput[6]+=(te-tb);
//   if(!seltrk0)return 1;
   sstat[7]++;

//---Trd 6
   tb=te;
   bool seltrd=Select_Trd();
   te=clock();
   cput[7]+=(te-tb);
//   if(!seltrd)return 1;
   sstat[8]++;

//---Rich 7
   tb=te;
   bool selrich=Select_Rich();
   te=clock();
   cput[8]+=(te-tb);
//   if(!selrich)return 1;
   sstat[9]++;   

//---Ecal 8
   tb=te;
   bool selecal=Select_Ecal();
   te=clock();
   cput[9]+=(te-tb);
//   if(!selecal)return 1;
   sstat[10]++;

//--Anti 9
   tb=te;
   bool selanti=Select_Anti();
   te=clock();
   cput[10]+=(te-tb);
//   if(!selanti)return 1;
   sstat[11]++;

//--TrkParticle 10
   tb=te;
   bool seltrkpart=Select_TrkParticle();
   te=clock();
   cput[11]+=(te-tb);
//   if(!seltrkpart)return 1;
   sstat[12]++;

   return true;
}

//---Header
bool  AMSAnalysis::Select_Header(){

   run=pev->Run();
   event=pev->Event();
   version=pev->Version();
//   cout<<"version="<<version<<endl;
   errorb=pev->fHeader.Error;
   time[0]=pev->fHeader.Time[0];
   time[1]=pev->fHeader.Time[1];
//--ISS
   yaw=   pev->fHeader.Yaw;
   pitch= pev->fHeader.Pitch;
   roll=  pev->fHeader.Roll;
//--Magnetic
   thetam=pev->fHeader.ThetaM;
   phim=pev->fHeader.PhiM;
//---only real data
   if(isreal==0)return true;

   isbadrun=pev->isBadRun(run);
   return true;
}

//---Trigger
bool AMSAnalysis::Select_Trigger(){

    nlevel1=pev->nLevel1();
    livetime=pev->LiveTime();
    if(nlevel1>0){
      Level1R *lv=pev->pLevel1(0);
//      livetime=  lv->LiveTime;
      physbpatt= lv->PhysBPatt;
      jmembpatt= lv->JMembPatt;
//Add MC
      if(isreal==0){
         jmembpatt1=physbpatt1=0;
         lv->RebuildTrigPatt(jmembpatt1,physbpatt1);
       }
//Add ISS
      else {
         jmembpatt2=physbpatt2=0; 
         lv->RestorePhysBPat();
         physbpatt2=lv->PhysBPatt;
         jmembpatt2=lv->JMembPatt;
      }
//---
      tofflag[0]=lv->TofFlag1;
      tofflag[1]=lv->TofFlag2;
      ecalflag=lv->EcalFlag;
      antipatt=  lv->AntiPatt;
    }
//--DAQ
    ndaqev=pev->nDaqEvent();
    ndaqerr=0;
    daqlen=0;
    daqtyerr=0;
    if(pev->fStatus&(1<<30))daqtyerr=1;
    if(ndaqev>=1){
       DaqEventR *pdaq=pev->pDaqEvent(0);
       if(pdaq->HasHWError())ndaqerr++;
       daqlen=pdaq->Length;
       if(pdaq->L3RunError())daqtyerr|=(1<<1);
       if(pdaq->L3EventError())daqtyerr|=(1<<2);
       if(pdaq->L3ProcError())daqtyerr|=(1<<3);
       for(int ino=0;ino<=6;ino++){
         if(pdaq->L3NodeError(ino))daqtyerr|=(1<<(4+ino));
       }
       daqjrmerr=pdaq->JRoomError;
       for(int ii=0;ii<24;ii++)daqjlen[ii]=pdaq->JLength[ii];
       daqsdlen=pdaq->SDLength();
       daqsdlenu=pdaq->SDLengthU();
    }
    else {
       daqjrmerr=0;
       for(int ii=0;ii<24;ii++)daqjlen[ii]=0;
       daqsdlen=daqsdlenu=0;
    }
    return true;
}
//---MC
bool AMSAnalysis::Select_MC(int opt){

//---
 if(isreal)return true;
 int mnpart=pev->nMCEventg();

//---MCEveng0
  if(opt==0){

    ntrmccl=pev->nTrMCCluster();
    ntofmccl=pev->nTofMCCluster();
//----
    MCEventgR *mcev0=pev->pMCEventg(0);
    MCEventgR *mcev=pev->GetPrimaryMC();
//  if(mcev0!=mcev)cout<<"Error MCEvent(0)"<<" run="<<run<<" event="<<event<<endl;
//----
    mpar=mcev->Particle;
    mmom=mcev->Momentum;
    mch= mcev->Charge;
    mmass=mcev->Mass;

//----Coo+Dir
    for(int ico=0;ico<3;ico++){
       mevcoo[ico]=mcev->Coo[ico]; mevdir[ico]=mcev->Dir[ico];
    }

  }
//---MCEvengI
  else {
 
    AMSPoint pnt(mevcoo[0],mevcoo[1],mevcoo[2]);
    AMSDir dir(mevdir[0],mevdir[1],mevdir[2]);

//---MCEveng0 TOF Pass
    AMSPoint tofpnt;
    float disedge;
    mtof_pass=0;
    int ntofp=0;
    for(int ilay=0;ilay<4;ilay++){
      int tk_pass=ParticleR::IsPassTOF(ilay,pnt,dir,tofpnt,disedge);
      if(tk_pass>=0){mtof_pass+=int(pow(10.,3-ilay));ntofp++;}
    }

//---MCEvengN Interaction Point Find Max Kine
    for(int ic=0;ic<3;ic++){mfcoo[ic]=-1000;mfdir[ic]=-1000;}
    for(int is=0;is<2;is++){
      mpare[is]=mparp[is]=mparc[is]=0;
    }
    for(int ipar=1;ipar<mnpart;ipar++){
      MCEventgR *mcev=pev->pMCEventg(ipar);
      int parid=mcev->Particle;
//    if(parid<0)cout<<"parid="<<parid<<endl;
      if(fabs(parid)==1||fabs(parid)==2||fabs(parid)==3)continue;//gamma+positron+electron
      if(mcev->Coo[2]<-138||mcev->Nskip==-2||fabs(mcev->Particle)==fabs(mpar))continue;
      if(mcev->parentID!=1)continue;
      if(mcev->Momentum>mpare[0]){
        mpare[0]=mcev->Momentum;
        for(int ic=0;ic<3;ic++){
          mfcoo[ic]=mcev->Coo[ic];
          mfdir[ic]=mcev->Dir[ic];
        }
      }
    }
//Again
    msume=msumc=0;
    for(int iz=0;iz<21;iz++){
      mevmom1[iz]=-1000;
      for(int ic=0;ic<3;ic++){mevcoo1[iz][ic]=-1000;mevdir1[iz][ic]=-1000;}
    }
    for(int ipar=1;ipar<mnpart;ipar++){
      MCEventgR *mcev=pev->pMCEventg(ipar);
      if(mcev->Coo[2]==mfcoo[2]){
         msume+=mcev->Momentum;
         msumc+=mcev->Charge;
         if(mcev->Momentum==mpare[0]){mpare[0]=mcev->Momentum;mparp[0]=mcev->Particle;mparc[0]=mcev->Charge;}//MaxE
         if(fabs(mcev->Charge)>mparc[1]||(fabs(mcev->Charge)==mparc[1]&&mcev->Momentum==mpare[0])){
            mpare[1]=mcev->Momentum;mparp[1]=mcev->Particle;mparc[1]=mcev->Charge;
          }//MaxCh
      }
//--AC
      int iskip=mcev->Nskip;
      if(iskip>=-1020&&iskip<=-1000){
         iskip=fabs(iskip)-1000;
         mevmom1[iskip]=mcev->Momentum;
         for(int ic=0;ic<3;ic++){
           mevcoo1[iskip][ic]=mcev->Coo[ic];
           mevdir1[iskip][ic]=mcev->Dir[ic];
        }
      }
    }

//---TrdMCCluster MaxE
    for(int ilay=0;ilay<20;ilay++){mtrdmom[ilay]=0;mtrdpar[ilay]=0;for(int ic=0;ic<3;ic++)mtrdcoo[ilay][ic]=0;}
    for(int itrd=0;itrd<pev->nTrdMCCluster();itrd++){
      TrdMCClusterR *trdmc=pev->pTrdMCCluster(itrd);
      int ilay=trdmc->Layer;
      if(trdmc->Ekin>mtrdmom[ilay]){
        mtrdmom[ilay]=trdmc->Ekin;
        mtrdpar[ilay]=trdmc->ParticleNo;
        for(int ic=0;ic<3;ic++)mtrdcoo[ilay][ic]=trdmc->Xgl[ic];
      }
    }

//--TrackMCCluster MaxE
    for(int ilay=0;ilay<9;ilay++){
      mtrmom[ilay]=mtrpri[ilay]=mtrpar[ilay]=mtrz[ilay]=0;
      for(int ic=0;ic<3;ic++){mtrcoo[ilay][ic]=0;mtrcoo1[ilay][ic]=0;}
    }
    for(int itr=0;itr<pev->nTrMCCluster();itr++){
      TrMCClusterR *trmc=pev->pTrMCCluster(itr);
      int il=fabs(trmc->GetTkId()/100);
      if     (il==8)il=1;
      else if(il<8)il=il+1;
      il--;
      if(trmc->GetMomentum()>mtrmom[il]){
        mtrmom[il]=trmc->GetMomentum();
        mtrpar[il]=trmc->GetPart();
        for(int ic=0;ic<3;ic++){mtrcoo[il][ic]=trmc->GetStartPoint()[ic];mtrcoo1[il][ic]=trmc->GetXgl()[ic];}
        mtrpri[il]=trmc->IsPrimary()?1:0;
        mtrz[il]=trmc->Status;
      }
    }

//--TofMCCluster
    for(int ilay=0;ilay<4;ilay++){
      mtofdep[ilay]=mtofpri[ilay]=0;
      mtofbar[ilay]=-1;
    }
    double mtofe[4][10]={{0}};
    double mtofepri[4][10]={{0}};
    for(int imccl=0;imccl<pev->nTofMCCluster();imccl++){
      TofMCClusterR *tofmc=pev->pTofMCCluster(imccl);
      int ilay=tofmc->GetLayer();
      int ibar=tofmc->GetBar();
      mtofe[ilay][ibar]+=tofmc->EdepR*1000;//GeV->MeV
      if(tofmc->ParentNo==0){
        mtofpri[ilay]=1;
        mtofepri[ilay][ibar]+=tofmc->EdepR*1000;
      }
    }
    for(int ilay=0;ilay<4;ilay++){
      double maxdep=0;
      for(int ibar=0;ibar<10;ibar++){
        if(mtofpri[ilay]==1&&mtofepri[ilay][ibar]>maxdep){//Primary
          mtofbar[ilay]=ibar;
          maxdep=mtofdep[ilay]=mtofepri[ilay][ibar];
        }
        else if(mtofpri[ilay]==0&&mtofe[ilay][ibar]>maxdep){//No
          mtofbar[ilay]=ibar;
          maxdep=mtofdep[ilay]=mtofe[ilay][ibar];
        }
      }
    } 


  }

  return true;

}


//---RTI
bool AMSAnalysis::Select_RTI(){
    
    irti=-1;
    if(isreal==0)return true;
    AMSSetupR::RTI a;
    if(pev->GetRTI(a)!=0)return true;

//--Time Cut Provide By VC
     bool cut[10]={0};
     cut[0]=(a.ntrig/a.nev>0.98);
     cut[1]=(a.npart/a.ntrig>0.07/1600*a.ntrig&&a.npart/a.ntrig<0.25);
     cut[2]=(a.lf>0.5);
     cut[3]=(a.zenith<40);
     cut[3]=(a.nerr>=0&&a.nerr/a.nev<0.1);
     cut[4]=(a.npart>0&&a.nev<1800);
     bool tcut=(cut[0]&&cut[1]&&cut[2]&&cut[3]&&cut[4]);

//----RTI var
     zenith=a.zenith;
     thetas=a.theta;
     phis=a.phi;
     rads=a.r;
     glat=a.glat;
     glong=a.glong;
     issaa=a.IsInSAA();
//---cutoff
     for(int ifv=0;ifv<4;ifv++){
       for(int ipn=0;ipn<2;ipn++){mcutoff[ifv][ipn]=a.cf[ifv][ipn];mcutoffi[ifv][ipn]=a.cfi[ifv][ipn];}
     }
//---more 
     rtilf=a.lf;
     rtinev=a.nev;
     rtinerr=a.nerr;
     rtintrig=a.ntrig;
     rtinpar=a.npart;
     rtigood=a.good;
     irti=0;
//----
     AMSPoint nexl[2]; 
     pev->GetRTInL1L9(nexl,a.utime,90);
     for(int uexl=0;uexl<2;uexl++){
       for(int ixy=0;ixy<2;ixy++)rtinexl[uexl][ixy]=nexl[uexl][ixy];
       AMSPoint nxyz,dxyz;
       pev->GetRTIdL1L9(uexl,nxyz,dxyz,a.utime,60);//PG-CIEMAT
       for(int ixy=0;ixy<3;ixy++)rticdif[uexl][ixy]=dxyz[ixy];
     }
//     return tcut;
     return true;
}

//---Particle
bool AMSAnalysis::Select_Particle(){

//---Load Status
   nparticle=pev->nParticle();
   nbeta=pev->nBeta();
   ncharge=pev->nCharge();
   ntrack=pev->nTrTrack();
   ntrrawcl=pev->nTrRawCluster();
   ntrcl=pev->nTrCluster();
   ntrhit=pev->nTrRecHit();
   nshow=pev->nEcalShower();
   necalhit=pev->nEcalHit();
   nrich=pev->nRichRing();
   nrichb=pev->nRichRingB();
   nrichh=pev->nRichHit();
//---
   ntrdtrack=pev->nTrdTrack();
   ntrdseg=pev->nTrdSegment();
   ntrdrawh=pev->nTrdRawHit();
   ntrdcl=pev->nTrdCluster();
   ntofclh=-1;
   anti_nhit=pev->nAntiCluster();

//---PreScale 
#ifdef _GOODTKTOFPAR_
  if(pev->nTrTrack()<1)return false;//require Tk
#endif
//----TkIn Particle
  float tkiqm[3]={0};//L1orL9 +L1L9 
  bool l1hy=0;
  bool l9hy=0;
  bool ngpar=0;
  for(int itr=0;itr<pev->nTrTrack();itr++){
    TrTrackR *trk=pev->pTrTrack(itr);
    l1hy=((trk->GetBitPatternJ()&(1<<0))>0);//L1Hit
    l9hy=((trk->GetBitPatternJ()&(1<<8))>0);//L9Hit
    float ntkq=-1;
#if  defined (_HEINNERPRESCALE_) || defined (_IONINNERPRESCALE_)
    ntkq=trk->GetInnerQ();
    if(ntkq>tkiqm[0]){
       tkiqm[0]=ntkq;//Inner
    }
#endif
    if(!l1hy&&!l9hy)continue;//Both
    if(ntkq<0)ntkq=trk->GetInnerQ();
    if(ntkq>tkiqm[1]){//L1 or L9
       tkiqm[1]=ntkq;
     }
    if(l1hy&&l9hy&&ntkq>tkiqm[2]){//L1&&L9
       tkiqm[2]=ntkq;
    }
#ifdef _SAVENEGSCALE_
    int nhiti=0;
    for(int ilay=0+1;ilay<9-1;ilay++){
       if((trk->GetBitPatternJ()&(1<<ilay))>0)nhiti++;
    }
    if(nhiti>=5&&ntkq>0.6){
      float tkrig=trk->GetRigidity();
      float nchis=trk->GetNormChisqY();
      if(nchis==0)nchis=9999;
      float obeta=0;
      for(int ibh=0;ibh<pev->nBetaH();ibh++){
        if(pev->pBetaH(ibh)->pTrTrack()==trk)obeta=pev->pBetaH(ibh)->GetBeta();
      }
      if(fabs(obeta)>0.4&&fabs(tkrig)>0.7&&tkrig/obeta<0&&nchis<20)ngpar=1;
    }
#endif
  }

  bool isintkbgpar=0;
//----Default
#if defined (_PRHEPRESCALE3_)
   isintkbgpar=(tkiqm[1]>1.5||tkiqm[2]>0);//L1Inner He,L1InnerL9 Pr
#elif defined (_HEL1L9PRESCALE_)
   isintkbgpar=(tkiqm[2]>1.5);//L1InnerL9 He+Full
#endif
//---P/He/Ion Inner+L1Inner
#if defined (_PRL1PRESCALE_)
   isintkbgpar=(isintkbgpar||tkiqm[1]>0);//L1Inner Pr+He+Ion 
#endif   
#if defined (_HEINNERPRESCALE_)
   isintkbgpar=(isintkbgpar||tkiqm[0]>1.5);//Inner He+Ion
#elif defined (_HEL1PRESCALE_)
   isintkbgpar=(isintkbgpar||tkiqm[1]>1.5);//L1Inner He+Ion
#endif
#if defined (_IONINNERPRESCALE_)
   isintkbgpar=(isintkbgpar||tkiqm[0]>2.5);//Inner Z>2
#elif defined (_IONL1PRESCALE_)
   isintkbgpar=(isintkbgpar||tkiqm[1]>2.5);//L1Inner Z>2
#endif
//----Unbiase
#ifdef _PRHESAVEUNBIASE_
   if(unbiasetag)isintkbgpar=(isintkbgpar||tkiqm[0]>0);//unbiased P+He
#endif
//-----
#ifdef _GOODTKTOFPAR_
  if(!isintkbgpar&&!unbiasetag0)return false;
#else
  bool isl1l9bgpar=Select_L1L9Particle();
  bool isl1l9mcpar=Select_MCL1L9Particle(); 
  if(!isintkbgpar&&!ngpar&&!isl1l9bgpar&&!isl1l9mcpar&&!unbiasetag0)return false;
#endif

//--TOFClusterH   
  if(ntofclh<0){
     TofRecH::BuildTofClusterH();
     ntofclh=pev->nTofClusterH();
   }
 
//--Paritcle
   TofRecH::BuildOpt=0;
   TofRecH::BuildBetaH();
   nbetah=pev->nBetaH();

//---Efficiency Mis-Match Continue to Search(For Z>=1)
  float tksq=0,tofsq=0;
  for(int ibh=0;ibh<nbetah;ibh++){
//---TOF sq
      BetaHR *betah=pev->pBetaH(ibh);
      float ntofsq=0;float ntksq=0;
      for(int ilay=0;ilay<4;ilay++){ 
         float tofql=betah->GetQL(ilay);
         if(tofql>0)ntofsq+=tofql;
      }
      TrTrackR *track=betah->pTrTrack();
      if(track){ntksq=track->GetInnerQ()*track->GetNhitsI();}//InnerQ*InnerNHits
//---First BetaH
     if(ibetah<0){ibetah=ibh;tofsq=ntofsq;}
//---Fist Track
     if(itrtrack<0&&betah->iTrTrack()>=0){//First Track Order
        ibetah=ibh; itrtrack=betah->iTrTrack();
        tofsq=ntofsq;tksq=ntksq;
     }
//--If Second Tracker Max TkQ
     else if(itrtrack>=0&&betah->iTrTrack()>=0&&ntksq>tksq){
        ibetah=ibh; itrtrack=betah->iTrTrack();
        tofsq=ntofsq;tksq=ntksq;
     }
//---If NoTrack Max TofQ
     else if(itrtrack<0&&betah->iTrTrack()<0&&ntofsq>tofsq){//Then Eneryg> Raw
        ibetah=ibh; itrtrack=betah->iTrTrack();
        tofsq=ntofsq;tksq=ntksq;
     }
//---end try
  }

///---Double Check(Recover) For Wrong-Tk
  if(nbetah>1&&itrtrack>=0&&ibetah>=0){
//---Check Good Candidate
     int fid=0;
     map<pair<int,pair<double,int> >,int>ibhg;
     for(int ibh=0;ibh<nbetah;ibh++){
       BetaHR *betah=pev->pBetaH(ibh);
       if(betah->iTrTrack()<0)continue;
       float obeta=betah->GetBeta();
       if(obeta<0.4)continue;
//---InnerNHit Suppress
       int ntk_hitb[2]={0};
       TrTrackR *track=betah->pTrTrack();
       ntk_hitb[0]=track->GetBitPatternXYJ();
       ntk_hitb[1]=track->GetBitPatternJ();
       int nhiti[2]={0};
       for(int ilay=0+1;ilay<9-1;ilay++){
         for(int ixy=0;ixy<2;ixy++){
           if((ntk_hitb[ixy]&(1<<ilay))>0)nhiti[ixy]++;
         }
       }
       if(nhiti[1]<5)continue;//Inner-Not enough
//--Charge Suppress Wrong-Tk
       double tkiq=track->GetInnerQ(obeta,fid);
       if(tkiq<0.6)continue;
       int npi=track->GetNPlanesI();
//--L1 or L9
       bool lhxy[2]={0};
       for(int iex=0;iex<2;iex++){
         int ilay=(iex==0)?0:8;
         lhxy[iex]=((ntk_hitb[1]&(1<<ilay))>0&&(ntk_hitb[0]&(1<<ilay))>0);//HitXY
         if(!lhxy[iex])continue;
         double tkexq=track->GetLayerJQ(ilay+1,obeta,fid);
         if     (tkiq>=2.5)          lhxy[iex]=(tkexq>2.2&&tkexq>tkiq*0.5);//Z>=3
         else if(tkiq>=1.5&&tkiq<2.5)lhxy[iex]=(tkexq>1.5&&tkexq<3);//Z>=2
         else                        lhxy[iex]=(tkexq>0.5&&tkexq<2);//Z>=1
       }
       int nex=0,iz=0;
       vector<pair<int,int> > mex;
       if(lhxy[0]&&lhxy[1])mex.push_back(make_pair(TrTrackR::kFitLayer8|TrTrackR::kFitLayer9,2));
       if(lhxy[0]){mex.push_back(make_pair(TrTrackR::kFitLayer8,1));nex++;}
       if(lhxy[1]){mex.push_back(make_pair(TrTrackR::kFitLayer9,1));nex++;}
       mex.push_back(make_pair(0,0));
//-----
       double tkrig=track->GetRigidity();
       int iid=(track->Gettrdefaultfit()&(~(TrTrackR::kFitLayer8|TrTrackR::kFitLayer9)));//inner tracker
       for(unsigned int isp=0;isp<mex.size();isp++){//highest priority first
         int mfit=iid|mex[isp].first;
         if(!track->FitDone(mfit))continue;
         if(isp==0)tkrig=track->GetRigidity(mfit);//max-span
         if(track->GetNormChisqY(mfit)<20){//max-span good-track
           tkrig=track->GetRigidity(mfit);
           iz=2+mex[isp].second;
           break;
         }
       }
//---Rigidity Suppress Wrong-Tk
       if(fabs(tkrig)<0.7)continue;
       iz+=1000*int(tkiq+0.5)+nex*100+npi*10;
       int nhitsw=2*(nex+npi)+(track->GetNhitsI()-npi);
       double quality=tkiq*fabs(tkrig)*nhitsw;
       pair<double,int> qval(quality,nhitsw);
       ibhg[make_pair(iz,qval)]=ibh;
     }
     if(ibhg.size()>=1){
       map<pair<int,pair<double,int> >,int>::reverse_iterator sit=ibhg.rbegin();
       map<pair<int,pair<double,int> >,int>::reverse_iterator rit=sit;
       ++rit;
       for(;rit!=ibhg.rend();++rit){
         int iz[2]={0},gz[2]={0},nex[2]={0},npi[2]={0},gchis[2]={0},nhitsw[2]={0};
         double quality[2]={0};
         for(int j=0;j<2;j++){
           iz[j]=(j==0)?(sit->first).first:(rit->first).first;
           gz[j]=iz[j]/1000;//charge
           if(gz[j]>=3)gz[j]=3;//Z=3 is the highest Z group
           nex[j]=(iz[j]/100)%10;
           npi[j]=(iz[j]/10)%10;
           gchis[j]=iz[j]%10;
           quality[j]=(j==0)?(sit->first).second.first:(rit->first).second.first;
           nhitsw[j]=(j==0)?(sit->first).second.second:(rit->first).second.second; 
         }
         int ib=0;
         if(gz[1]==gz[0]){//the same group
           if(iz[0]/1000>iz[1]/1000){//different charge
             if     (nex[0]>=nex[1]&&npi[0]>npi[1])ib=0;
             else if(npi[0]>npi[1])ib=(quality[1]>quality[0]&&(gchis[1]>0&&gchis[0]==0))?1:0;
             else if(nex[0]>nex[1])ib=(quality[1]>quality[0]&&(gchis[1]>=gchis[0]-1&&gchis[1]>0))?1:0;
             else                  ib=(quality[1]>quality[0]&&(gchis[1]>=gchis[0]-1||gchis[1]>0))?1:0;
           }
           else if(nex[0]>nex[1]){//same charge, different external 
             if     (npi[0]>npi[1])ib=0;
             else if(npi[0]==npi[1]&&nex[0]==nex[1]+2)ib=(quality[1]>quality[0]&&(gchis[1]>=gchis[0]-1&&gchis[1]>0))?1:0;
             else                  ib=(quality[1]>quality[0]&&(gchis[1]>=gchis[0]-1||gchis[1]>0))?1:0;
           }
           else if(npi[0]>npi[1]){
             ib=(quality[1]>quality[0]&&(gchis[1]>=gchis[0]+2))?1:0; 
           }
           else if(gchis[0]>0&&gchis[1]==0){
             ib=(quality[1]>quality[0]*2.&&nhitsw[1]>nhitsw[0])?1:0;
           }
           else {
             ib=(quality[1]>quality[0])?1:0; 
           }
         }
         else if(gz[1]<=1){//Z>=3||Z=2 VS Z=1
           if(gz[0]>=3){//Z>=3 VS Z=1
             if(npi[0]<=3&&nex[1]==nex[0]+2){if(quality[1]>quality[0]&&(gchis[1]==4||gchis[0]==0))ib=1;}//Z>=3 3-pla VS Z=1 L1IL9
           }
           else {//Z=2 VS Z=1
             if     ((npi[0]<=3&&nex[0]==0)&&nex[1]>nex[0]){if(quality[1]>quality[0]&&(gchis[1]>=3||gchis[0]==0))ib=1;}//Z=2 3-pla VS Z=1 L1Inner
             else if(npi[1]>=npi[0]&&nex[1]==nex[0]+2){if(quality[1]>quality[0]&&(gchis[1]==4||gchis[0]==0))ib=1;}//Z=2 4-pla VS Z=1 L1IL9
           }
         }
         else {//Z>=3 VS Z=2
           if     ((npi[0]<=3&&nex[0]==0)&&(npi[1]>npi[0])){if(quality[1]>quality[0]&&gchis[1]>=gchis[0])ib=1;}//Z>=3 3-pla VS Z=2 4-pli
           else if((npi[0]<=3&&nex[0]==0)&&(nex[1]>nex[0])){if(quality[1]>quality[0]&&(gchis[1]>=3||gchis[0]==0))ib=1;}//Z>=3 3-pla VS Z=2 L1I
           else if(npi[1]>=npi[0]&&nex[1]==nex[0]+2){if(quality[1]>quality[0]&&(gchis[1]==4||gchis[0]==0))ib=1;}//Z>=3 4-pla VS Z=2 L1IL9
         }
         if(ib==1)sit=rit;
       }//end of iteration
       ibetah=sit->second;
       itrtrack=pev->pBetaH(ibetah)->iTrTrack();
     }//end of size
//----
  }

//  cout<<"ubetah="<<ibetah<<endl;
#ifdef _GOODTKTOFPAR_
  if(ibetah<0||itrtrack<0)return false;
#endif


//---2nd Particle
    Select_2ndParticle();

//--Link Particle
    Select_LinkParticle();


  return true;

}

bool  AMSAnalysis::Select_MCL1L9Particle(){

  if(isreal)return false; 
//--Hit Number
   bool mmlh[9]={0};
   for(int il=0;il<NTKL;il++){
      float coox=(tk_pz[il]-mevcoo[2])*mevdir[0]/mevdir[2]+mevcoo[0];
      float cooy=(tk_pz[il]-mevcoo[2])*mevdir[1]/mevdir[2]+mevcoo[1];
      float disr=coox*coox+cooy*cooy;
      if(il==9-1)mmlh[il]=(fabs(coox)<cirr[il]&&fabs(cooy)<ciry[il]);
      else       mmlh[il]=(sqrt(disr)<cirr[il]&&fabs(cooy)<ciry[il]);
   }
//---
   int mmnhitf=0;
   for(int il=0+1;il<9-1;il++){
      if(mmlh[il]){mmnhitf++;}
   }

   bool mcl1l9tkpar=(mmlh[0]&&mmlh[8]);//L1L9
   float mchcut=1.5;
#ifdef _PRL1PRESCALE_
   mchcut=0.5;
#endif 
#if defined (_PRL1PRESCALE_) || defined (_HEL1PRESCALE_) || defined (_IONL1PRESCALE_)
   if(!mcl1l9tkpar&&(fabs(mch)>mchcut)){//MC Z>2
       mcl1l9tkpar=(mmnhitf>=5&&mmlh[0]); //L1Inner
   }
#endif
  return mcl1l9tkpar;
}


bool AMSAnalysis::Select_L1L9Particle(){
 
  tk_l9qs=tk_l1qs=tk_l1qvs=tk_l9qvs=tk_l1qvr=tk_l9qvr=0;
  tk_l9q=tk_l1q=0;
  tk_l9qxy[0]=tk_l9qxy[1]=tk_l1qxy[0]=tk_l1qxy[1]=0;
  for(int iyh=0;iyh<2;iyh++){tk_exqvn[iyh][0]=tk_exqvn[iyh][1]=0;}
  tk_l1mds=tk_l9mds=1000;
  //---Pre Scale
  if(ntrdtrack<1&&ntrack<1){return false;}


  AMSPoint postr;AMSDir dirtr;
//----InnerTk Interpolate L1 L9 Geometry
  float tkiqm[3]={0};
  for(int itr=0;itr<pev->nTrTrack();itr++){
    bool tkl1g=0,tkl9g=0;
    for(int ilay=0;ilay<NTKL;ilay++){//L1Pos+L9Pos Geometry
        if(ilay!=0&&ilay!=NTKL-1)continue;
        double length=pev->pTrTrack(itr)->InterpolateLayerJ(ilay+1,postr,dirtr);//L1+L9
        if(length<=-1000)continue;
        float disr=sqrt(postr[0]*postr[0]+postr[1]*postr[1]);
        if(ilay==0)tkl1g=(disr<cirr[ilay]&&fabs(postr[1])<ciry[ilay]);//L1
        else       tkl9g=(fabs(postr[0])<cirr[ilay]&&fabs(postr[1])<ciry[ilay]);//L9
     }
     if(!tkl1g&&!tkl9g)continue;
     float ntkq=pev->pTrTrack(itr)->GetInnerQ();
     if(ntkq>tkiqm[1])tkiqm[1]=ntkq;
     if(tkl1g&&tkl9g&&ntkq>tkiqm[2])tkiqm[2]=ntkq;
  }

///---Ecal Events
  float ecalen=0;
  for(int ish=0;ish<pev->nEcalShower();ish++){
     EcalShowerR *show=pev->pEcalShower(ish);
     if(show->EnergyD>ecalen)ecalen=show->EnergyD; 
  }
  bool ecalpar=(ecalen/1000.>0.6);

  bool l1l9tkpar=0;//L1Inner
//------Default
#if defined (_PRHEPRESCALE3_)
  l1l9tkpar=(tkiqm[1]>1.5||tkiqm[2]>0);//L1Inner He,L1InnerL9 Pr
#elif defined (_HEL1L9PRESCALE_)
  l1l9tkpar=(tkiqm[2]>1.5);//L1InnerL9 He
#endif
//------L1Inner
#ifdef _PRL1PRESCALE_
  l1l9tkpar=(l1l9tkpar||tkiqm[1]>0);//L1Inner Pr
#endif
#ifdef _HEL1PRESCALE_
  l1l9tkpar=(l1l9tkpar||tkiqm[1]>1.5);//L1Inner He
#endif
#ifdef _IONL1PRESCALE_
  l1l9tkpar=(l1l9tkpar||tkiqm[1]>2.5);//L1Inner Z>2
#endif
  if(unbiasetag){
    l1l9tkpar=(tkiqm[1]>1.5||tkiqm[2]>0);//L1Inner He,L1InnerL9 Pr
//#ifdef _HEINNERPRESCALE_
    l1l9tkpar=(l1l9tkpar||tkiqm[1]>0);
//#endif
  }

//----TrdInterpolate L1 L9 Geometry
  int l1l9trdpar=0; 
  for(int itr=0;itr<pev->nTrdTrack();itr++){
    bool trdl1g=0,trdl9g=0;
    for(int ilay=0;ilay<NTKL;ilay++){//L1Pos+L9Pos Geometry
      if(ilay!=0&&ilay!=NTKL-1)continue;
      pev->pTrdTrack(itr)->Interpolate(tk_pz[ilay],postr,dirtr);//L1+L9
      float disr=sqrt(postr[0]*postr[0]+postr[1]*postr[1]);
      if(ilay==0)trdl1g=(disr<cirr[ilay]&&fabs(postr[1])<ciry[ilay]);//L1
      else       trdl9g=(fabs(postr[0])<cirr[ilay]&&fabs(postr[1])<ciry[ilay]);//L9
    }
    if(trdl1g&&trdl9g){l1l9trdpar=2;break;}//L1InnerL9
    else if(trdl1g)   {l1l9trdpar=1;}//L1Inner
  }

//---PreScale
   bool tktrdscal=0;
#if defined (_HEL1L9PRESCALE_) || defined (_PRHEPRESCALE3_)
   tktrdscal=(l1l9tkpar||l1l9trdpar>=2);//Track|TRDL1InnerL9G
#endif
#if defined (_PRL1PRESCALE_) || defined (_HEL1PRESCALE_) || defined (_IONL1PRESCALE_)
   tktrdscal=(l1l9tkpar||l1l9trdpar>=1);//Track|TRDL1InnerG
#endif
  if(unbiasetag)tktrdscal=(l1l9tkpar||l1l9trdpar>=1);//Track|TRDL1InnerG
  if(!tktrdscal)return false;

//------
  float ubeta=1,urig=0,umass=0,udxdz=0,udydz=0;
  int umul=0; 
//---L19Hit
  for(int itkh=0;itkh<pev->nTrRecHit();itkh++){//MaxQ
    TrRecHitR* tkhit=pev->pTrRecHit(itkh);
    int layj=tkhit->GetLayerJ();
    if(layj==1||layj==9){//have to exclude noisy strip
      float nlqy=tkhit->GetQ(1,ubeta,urig,umass,umul,udxdz,udydz);//YQ
      float nlqx=tkhit->GetQ(0,ubeta,urig,umass,umul,udxdz,udydz);//XQ
      float nlq =tkhit->GetQ(2,ubeta,urig,umass,umul,udxdz,udydz);
      int   nqst=tkhit->GetQStatus();//have to define multiplicity here
      float &l19qy= (layj==1)? tk_l1qxy[1]: tk_l9qxy[1];//L19Y
      float &l19qx= (layj==1)? tk_l1qxy[0]: tk_l9qxy[0];//L19X
      float &l19q=  (layj==1)? tk_l1q:      tk_l9q;
      int   &l19qs= (layj==1)? tk_l1qs:     tk_l9qs;
      int   &l19i = (layj==1)? l1i:         l9i;//L19index
      if ((nlqy>1.2*l19qy)||(nlqy>0.99*l19qy&&nlqx>1.05*l19qx)){l19i=itkh;l19qy=nlqy;l19qx=nlqx;l19q=nlq;l19qs=nqst;}
    }
  }


  bool tkl1scal=0;
#ifdef _HEL1L9PRESCALE_ 
  tkl1scal=(l1i>=0&&tk_l1q>1.5&&l1l9trdpar>=2);
#endif
#ifdef _PRHEPRESCALE3_
  tkl1scal=(l1i>=0&&tk_l1q>1.5&&l1l9trdpar>=2);
  tkl1scal=(tkl1scal||(l1i>=0&&tk_l1q>0.6&&l1l9trdpar>=2&&ecalpar));//ECAL Events
#endif
//---------
#ifdef _PRL1PRESCALE_
  tkl1scal=(tkl1scal||(l1i>=0&&tk_l1q>0.5&&l1l9trdpar>=1));//PrL1InnerG
#endif
#ifdef _HEL1PRESCALE_
  tkl1scal=(tkl1scal||(l1i>=0&&tk_l1q>1.5&&l1l9trdpar>=1));//HeL1InnerG
#endif
#ifdef _IONL1PRESCALE_
  tkl1scal=(tkl1scal||(l1i>=0&&tk_l1q>2.5&&l1l9trdpar>=1));//Z>2 L1InnerG
#endif
  if(unbiasetag){tkl1scal=(tkl1scal||(l1i>=0&&tk_l1q>1.5&&l1l9trdpar>=1)||(l1i>=0&&tk_l1q>0.6&&l1l9trdpar>=2));}//(HeL1InnerG || PrL1InnerL9G)
  if(!tkl1scal)return l1l9tkpar;

//----Build TrdParticle
//   TofRecH::BuildOpt=31000;
  TofRecH::BuildOpt=1000;//TrdParticle
  TofRecH::ReBuild();

//-----BuildTofClusterH
  ntofclh=pev->nTofClusterH();
  int nbetahs=pev->nBetaH();

//--Select TrdBetaH
  float tofsqs=0;
  for(int ibh=0;ibh<nbetahs;ibh++){
     BetaHR *betah=pev->pBetaH(ibh);
     if(betah->iTrdTrack()<0)continue; 
     float ntofsq=0;
     for(int ilay=0;ilay<4;ilay++){
        float tofql=betah->GetQL(ilay);
        if(tofql>0)ntofsq+=tofql;
      }
      if(ntofsq>tofsqs){
       ibetahs=ibh;
       tofsqs=ntofsq;
     }
  }
  if(ibetahs<0){return l1l9tkpar;}
  BetaHR *betah=pev->pBetaH(ibetahs);
  itrdtracks=betah->iTrdTrack();
  tof_nhits=betah->GetSumHit();
  tof_hsumhus=betah->GetUseHit();
  TofBetaPar normbetapar=betah->GetNormTofBetaPar();
  tof_chists=normbetapar.Chi2T;
  tof_chiscs=normbetapar.Chi2C;
  betahs=betah->GetBeta();
  float tofq=0; int ntofq=0;
  for(int ilay=0;ilay<4;ilay++){
    tof_qls[ilay]=betah->GetQL(ilay,2,TofClusterHR::DefaultQOpt,111111,betahs,0);//beta correction only
    if(tof_qls[ilay]>0){tofq+=tof_qls[ilay]; ntofq++;}
  }
  if(ntofq>=1)tofq/=ntofq;

//---Ecal Match
  double ecalz=-143; double utm;
  betah->TInterpolate(ecalz,postr,dirtr,utm);
  for(int ic=0;ic<3;ic++){ecal_pos[ic]=postr[ic];}
  ecal_ens=-1000; ecal_dis=1000;
  for(int ish=0;ish<pev->nEcalShower();ish++){
    EcalShowerR *show=pev->pEcalShower(ish);
    betah->TInterpolate(show->Entry[2],postr,dirtr,utm);
    double necal_dis=pow(show->Entry[0]-postr[0],2)+pow(show->Entry[1]-postr[1],2);
    necal_dis=sqrt(necal_dis);
    if(necal_dis<ecal_dis){
      ecal_dis=necal_dis;
      ecal_ens=show->EnergyD;
    }
  }

//--L1L9Hit X+Y Match + L1L9Q
  for(int uexl=0;uexl<2;uexl++){
    int ilay=(uexl==0)?0:8; 
    int iexl=(uexl==0)?l1i:l9i;
    float qex=(uexl==0)?tk_l1q:tk_l9q;
    if(iexl<0)continue; 
    TrRecHitR* tkhit=pev->pTrRecHit(iexl);
    TVector3 ptrk,dtrk;
    betah->TInterpolate(tk_pz[ilay],postr,dirtr,utm);
    for(int ii=0;ii<3;ii++){ptrk[ii]=postr[ii]; dtrk[ii]=dirtr[ii];}
    int mult; double dist,xaddress,dadz[2]={0};
    AMSPlaneM plm=tkhit->FindMult(dist,mult,xaddress,ptrk,dtrk);
    //Check Proton+Helium
    float cutexdis=10;//10cm
    if(dist>cutexdis&&qex<2.5&&tofq<2.5){//Low amplitude Refind L1|L9 hit
      float qexn=0,exprobn=0;
      for(int itkh=0;itkh<pev->nTrRecHit();itkh++){
        TrRecHitR* tkhit1=pev->pTrRecHit(itkh);
        if(tkhit1->GetLayerJ()!=ilay+1||tkhit1->GetYCluster()==0||tkhit1->GetXCluster()==0||itkh==iexl)continue;
        int mult1; double dist1;
        AMSPlaneM plm1=tkhit1->FindMult(dist1,mult1,xaddress,ptrk,dtrk); 
        if(dist1<cutexdis){
          float nlq=tkhit1->GetQ(1,ubeta,urig,umass,umul,udxdz,udydz);//MaxY
          float nlp=tkhit1->GetProb();//X
          bool qmatch=(tofq>1.5)? (nlq>1.5): (nlq>0.5);//charge match
          if(qmatch&&(nlq>1.1*qexn||(nlq>0.99*qexn&&nlp>exprobn))){tkhit=tkhit1;plm=plm1;dist=dist1;mult=mult1;qexn=nlq;exprobn=nlp;}
        }
      }
    }
    TVector3 dloc=plm.getDLocal();
    for(int ixy=0;ixy<2;ixy++)dadz[ixy]=(dloc[2]!=0)?dloc[ixy]/dloc[2]:0;
    float &mds=(uexl==0)?tk_l1mds:tk_l9mds;
    float &qvr=(uexl==0)?tk_l1qvr:tk_l9qvr;
    float &qvs=(uexl==0)?tk_l1qvs:tk_l9qvs;
    mds=dist;
    qvr=tkhit->GetQ(2,betahs,urig,umass,umul,udxdz,udydz);
    qvs=tkhit->GetQ(2,betahs,urig,umass,mult,dadz[0],dadz[0]);
    if(qvs>1.5){
      tk_exqvn[0][uexl]=tkhit->GetQYJ(2,betahs,urig,mult,dadz[0],dadz[1]);
      tk_exqvn[1][uexl]=tkhit->GetQH(2,betahs,urig,mult,dadz[0],dadz[1]);
    }
  }

//---Interpolate TkL
  int nhiti=0;
  bool trdl1l9g[2]={0};
  for(int ilay=0;ilay<NTKL;ilay++){
    betah->TInterpolate(tk_pz[ilay],postr,dirtr,utm);
    for(int ic=0;ic<3;ic++){tk_pos1s[ilay][ic]=postr[ic];if(ilay==0)tk_dir1s[ic]=dirtr[ic];}
    float disr=sqrt(postr[0]*postr[0]+postr[1]*postr[1]);
    if(ilay!=0&&ilay!=NTKL-1){  
      if(disr<cirr[ilay]&&fabs(postr[1])<ciry[ilay])nhiti++;
    }
    else { 
      if(ilay==0)trdl1l9g[0]=(disr<cirr[ilay]&&fabs(postr[1])<ciry[ilay]);//L1
      else       trdl1l9g[1]=(fabs(postr[0])<cirr[ilay]&&fabs(postr[1])<ciry[ilay]);//L9
    }
  }

//---TrdTrack
  TrdTrackR *trdtracks=pev->pTrdTrack(itrdtracks);
  if(trdtracks){
    trd_thetas=trdtracks->Theta;
    trd_phis=trdtracks->Phi;
    for(int ic=0;ic<3;ic++)trd_coos[ic]=trdtracks->Coo[ic];
    trdtracks->ComputeCharge(1);
    trd_rqs=trdtracks->Q;
  }

//--TrdParticle Cutoff
  if(isreal){
    double cfp=0;
    AMSDir dirp(tk_dir1s[0],tk_dir1s[1],tk_dir1s[2]);
    if(dirp[2]>0)dirp=dirp*(-1);
    if(betahs<0)dirp=dirp*(-1);
    pev->GetStoermerCutoff(cfp,1,dirp);
    cutoffps=cfp;//positive charge cutoff
  }
  else cutoffps=mmom/mch;

  bool l1trdpar=0;
#ifdef _HEL1L9PRESCALE_
  l1trdpar=(tk_l1qvs>1.5&&tk_l1mds<10&&(tof_qls[0]>1.5||tof_qls[1]>1.5)&&(tof_qls[2]>1.5||tof_qls[3]>1.5)&&trdl1l9g[1]);//L1IL9-He
#endif
#ifdef _PRHEPRESCALE3_
  l1trdpar=(tk_l1qvs>1.5&&tk_l1mds<10&&(tof_qls[0]>1.5||tof_qls[1]>1.5)&&(tof_qls[2]>1.5||tof_qls[3]>1.5)&&trdl1l9g[1]);//L1IL9-He
  l1trdpar=(l1trdpar||(tk_l1qvs>0.6&&tk_l1mds<10&&(tof_qls[0]>0.6||tof_qls[1]>0.6)&&(tof_qls[2]>0.6||tof_qls[3]>0.6)&&trdl1l9g[1]&&ecalpar));//L1IL9E-Pr
#endif
#ifdef _PRL1PRESCALE_
  l1trdpar=(l1trdpar||(tk_l1qvs>0.5&&tk_l1mds<10&&(tof_qls[0]>0.5&&tof_qls[1]>0.5)&&(tof_qls[2]>0.5&&tof_qls[3]>0.5)));//L1I-Pr
#endif
#ifdef _HEL1PRESCALE_
  l1trdpar=(l1trdpar||(tk_l1qvs>1.5&&tk_l1mds<10&&(tof_qls[0]>1.5||tof_qls[1]>1.5)&&(tof_qls[2]>1.5||tof_qls[3]>1.5)));//L1I-He
#endif
#ifdef _IONL1PRESCALE_
  l1trdpar=(l1trdpar||(tk_l1qvs>2.5&&tk_l1mds<10&&(tof_qls[0]>2.5||tof_qls[1]>2.5)&&(tof_qls[2]>2.5||tof_qls[3]>2.5)));//L1I-Ion
#endif
  if(unbiasetag){
    bool l1trdparbz=(tk_l1qvs>1.5&&tk_l1mds<10&&(tof_qls[0]>1.5||tof_qls[1]>1.5)&&(tof_qls[2]>1.5||tof_qls[3]>1.5));//L1I-He
    bool l1trdparpr=(tk_l1qvs>0.6&&tk_l1mds<10&&(tof_qls[0]>0.6||tof_qls[1]>0.6)&&(tof_qls[2]>0.6||tof_qls[3]>0.6)&&trdl1l9g[1]); //L1IL9-Pr
    l1trdpar=(l1trdpar||l1trdparbz||l1trdparpr);
  }
  if(!l1trdpar){return l1l9tkpar;}

  return true;
}


bool AMSAnalysis::Select_2ndParticle(){

//---Secondary TrTrack
      betah2r=0;
      betah2p=0;
      betah2q=0;
      betah2hb[0]=betah2hb[1]=0;
      betah2ch=10000;

//---1nd Not Exist
     if(ibetah<0||itrtrack<0||nbetah<=1)return false;
     BetaHR *betah=pev->pBetaH(ibetah);
     float sbeta=betah->GetBeta();

//---Require Tof+Tracker
     for(int ibh=0;ibh<nbetah;ibh++){
         if(ibh==ibetah)continue;//not selected
         BetaHR *betah2=pev->pBetaH(ibh);
         if(betah2->iTrTrack()<0||betah2->iTrTrack()==itrtrack)continue;//second Tracker
//---Non First BetaH Matched TofClusterH 
         int nclmatch=0;
         for(int ilay=0;ilay<NTOFL;ilay++){
           if(!betah2->TestExistHL(ilay))continue;
           if(betah2->GetClusterHL(ilay)==betah->GetClusterHL(ilay))continue;
           nclmatch++;
         }
//---Rigidity+HitB
         float nrig= betah2->pTrTrack()->GetRigidity();
         float nchis=betah2->pTrTrack()->GetChisq();
         float ntkq= betah2->pTrTrack()->GetInnerQ();
         int ntkhb[2];
         ntkhb[0]=betah2->pTrTrack()->GetBitPatternXYJ();
         ntkhb[1]=betah2->pTrTrack()->GetBitPatternJ();
//----
         int nhit2i[2][2]={0};
         for(int ilay=0+1;ilay<9-1;ilay++){
             for(int ixy=0;ixy<2;ixy++){
               if((betah2hb[ixy]&(1<<ilay))>0)nhit2i[0][ixy]++;
               if((ntkhb[ixy]&(1<<ilay))>0)nhit2i[1][ixy]++;
             }
         }
         bool goodtkpr=(nhit2i[0][0]>=3&&nhit2i[0][1]>=5&&betah2r/sbeta>0);
         bool goodtknw=(nhit2i[1][0]>=3&&nhit2i[1][1]>=5&&nrig/sbeta>0);
//-------
         int fw=0;
//----Good TkNW
         if(goodtknw){
            if     (!goodtkpr)               {fw=1;}//tk-first now found
            else if(fabs(nrig)>fabs(betah2r)){fw=1;}//already found Largest Rig
         }
         else if(!goodtkpr){
            if     (nclmatch>betah2p)               {fw=1;}
            else if(nclmatch==betah2p&&ntkq>betah2q){fw=1;}
         }
         if(fw!=1)continue;
         betah2r=nrig;
         betah2p=nclmatch;
         betah2q=ntkq;
         betah2hb[0]=ntkhb[0]; betah2hb[1]=ntkhb[1];
         betah2ch=nchis;
//---
     }
     
    if(betah2q==0)return false;
    return true;

}


bool AMSAnalysis::Select_LinkParticle(){

//--Particle index    
   if(ibetah>=0&&itrtrack>=0){
     for(int ipar=0;ipar<nparticle;ipar++){
       if(pev->pParticle(ipar)->iTrTrack()==itrtrack){iparindex=ipar;break;}
     }
   }
  else if(ibetah<0&&itrtrack<0){
     float tksq=0;
     for(int itr=0;itr<ntrack;itr++){
        TrTrackR *trk=pev->pTrTrack(itr);
        float ntksq=trk->GetQ();
        if(ntksq>tksq){itrtrack=itr;tksq=ntksq;}
     }
   }

//-----index 
   if(iparindex>=0){
     btstat=pev->pParticle(iparindex)->BT_status;
     ibeta=pev->pParticle(iparindex)->iBeta();
     itrdtrack=pev->pParticle(iparindex)->iTrdTrack();
     irich=pev->pParticle(iparindex)->iRichRing();
     irichb=pev->pParticle(iparindex)->iRichRingB();
     icharge=pev->pParticle(iparindex)->iCharge();
//-----
     trd_pass=(pev->pParticle(iparindex)->IsInsideTRD())?1:0;
  }
  else if(itrtrack>=0){
    for(int ir=0;ir<nrich;ir++){
      if(pev->pRichRing(ir)->iTrTrack()==itrtrack){irich=ir;break;}
    }
    for(int irb=0;irb<nrichb;irb++){
       if(pev->pRichRingB(irb)->iTrTrack()==itrtrack){irichb=irb;break;}
    }
    for(int ib=0;ib<nbeta;ib++){
       if(pev->pBeta(ib)->iTrTrack()==itrtrack){ibeta=ib;break;}
    }
    for(int ich=0;ich<ncharge;ich++){
      if(pev->pCharge(ich)->iBetaH()==ibetah){icharge=ich;break;}
    }
  }

  show_dis=1000;
  if(ibetah>=0){
//----Ecal Match 
    BetaHR *betah=pev->pBetaH(ibetah);
    double ecalz=-143; double utm;
    AMSPoint postr;AMSDir dirtr;
    betah->TInterpolate(ecalz,postr,dirtr,utm);
    for(int ic=0;ic<3;ic++){ec_pos[ic]=postr[ic];ec_dir[ic]=dirtr[ic];}
    for(int ish=0;ish<pev->nEcalShower();ish++){
        EcalShowerR *show=pev->pEcalShower(ish);
        betah->TInterpolate(show->Entry[2],postr,dirtr,utm);
        double necal_dis=pow(show->Entry[0]-postr[0],2)+pow(show->Entry[1]-postr[1],2);
        necal_dis=sqrt(necal_dis);
        if(necal_dis<show_dis){
           show_dis=necal_dis;
           ishow=ish; 
        }
    }
//---Particle Cutoff
    if(isreal){
      double cfp=0;
      float obeta=betah->GetBeta();
      betah->TInterpolate(tk_pz[obeta>0?0:8],postr,dirtr,utm);//L1||L8 direction
      if(dirtr[2]>0)dirtr=dirtr*(-1);
      if(obeta<0)dirtr=dirtr*(-1);
      for(int ipn=0;ipn<2;ipn++){
          int sign=(ipn==0)?-1:1;//negtive(postive) charge
          pev->GetStoermerCutoff(cfp,sign,dirtr);
          cutoffp[ipn]=cfp;
       }
     }
     else cutoffp[0]=cutoffp[1]=mmom/mch; 
 } 
 

  return true; 
 
}

//---Init Trk ReFit
int AMSAnalysis::InitTrkReFit(float tkiq,float tkfz,int &refit,float &mass,int &chrg,int algo){

//--Fiting Charge and Mass 
    chrg=int(tkiq+0.5);
    if(chrg<1)chrg=1;
    mass=TrFit::Mproton;//Proton
    if(chrg>=2)mass=TrFit::Mhelium/2*chrg;//Helium+Ion
    if(!isreal&&version<1106)TRFITFFKEY.ErcHeY=0.65;//OldVersion MC
    else                     TRFITFFKEY.ErcHeY=0;//Newweight (NewMC||ISS)
//--ISS Setting
    if(isreal){
       TRFITFFKEY.Zshift=2;
#if defined (_USENOLINEARCOR_)
       TrClusterR::UnsetLinearityCorrection();//No Linearity Correction 
#endif
//       refit=3;//Force ReFit
//       if(tkfz<0||chrg==int(tkfz+0.5))refit=1;
       if(tkfz<0||chrg==int(tkfz+0.5))refit=0;
       else                           refit=3;
//       if(refit==3)cout<<"refit="<<refit<<" chrg="<<tkiq<<","<<tkfz<<endl;
    }
//--MC Smearing
    else  {
      if(issmearext==0){//once
        TrExtAlignDB::SmearExtAlign();//MC Smear Ext-Layer
        TRCLFFKEY.UseSensorAlign = 0;
        issmearext=1;
      }
      TRFITFFKEY.Zshift=-1;
      refit=3;//Force ReFit
    }

    return 0;
}

//---Trk
bool AMSAnalysis::Select_Trk(int sel){
    int mfit=-1;
    TrTrackR *trk=pev->pTrTrack(itrtrack);
    if(!trk)return false;
    float umass=0; int fid=0;
    BetaHR *betah=pev->pBetaH(ibetah);
    float obeta=betah?betah->GetBeta():1;
    float frig=trk->GetRigidity(fid);

//---Charge
   if(sel==0){

//---New charge All 
     for(int isp=0;isp<2;isp++){
//--Float Z
        mean_t tkqa;
        if(isp==0)tkqa=trk->GetQ_all(obeta,fid);
        else      {tk_qr=trk->GetInnerQ(); tkqa=trk->GetInnerQ_all(obeta,fid);}
        tk_q[isp]   =tkqa.Mean;
        tk_qrms[isp]=tkqa.RMS;
        tk_qhit[isp]=tkqa.NPoints;
//---TkZ
       vector<like_t> likelihood;
       if(isp==0)trk->GetZ(likelihood,obeta,fid);
       else      trk->GetInnerZ(likelihood,obeta,fid);
       int nz=likelihood.size()>=2?2:likelihood.size();
       for(int ip=0;ip<nz;ip++){//Max + Second Max
        tk_rz[ip][isp]=likelihood[ip].Z;
        tk_rpz[ip][isp]=likelihood[ip].GetNormProb();
      }
//----YJ
       if(isp==1){
         for(int ixy=0;ixy<3;ixy++){
           if(tk_q[1]<=1.5)continue;
           TrQYJTrack tkqyj=trk->GetQYJ_all(ixy,obeta,fid);
           tk_qin[0][ixy]=tkqyj.InnerQ;
           tk_qrmn[0][ixy]=tkqyj.InnerQRMS;
           mean_t tkqha=trk->GetInnerQH_all(ixy,obeta,fid);
           tk_qin[1][ixy]=tkqha.Mean;
           tk_qrmn[1][ixy]=tkqha.RMS;
//           cout<<"ixy="<<ixy<<" innerQ="<<tk_q[isp]<<","<<tk_qin[0][ixy]<<","<<tk_qin[1][ixy]<<endl;
         }
       }
     }

//---New Charge Layer     
     for(int il=0;il<9;il++){
       tk_lid[il]=0;
       tk_qls[il]=0;
       tk_edep[il]=0;
       if(il==0)tk_qlsn=0;
       tk_ql2[il][0]=tk_ql2[il][1]=0;
       tk_ql[il]=trk->GetLayerJQ(il+1,obeta,fid);
       for(int ixy=0;ixy<3;ixy++){
         tk_qln[0][il][ixy]=tk_qln[1][il][ixy]=0;
         if(tk_q[1]<=1.5)continue;
         TrQYJHit tkyjhit=trk->GetLayerQYJ_all(il+1,ixy,obeta,fid);
         tk_qln[0][il][ixy]=tkyjhit.Q;
         tk_qln[1][il][ixy]=trk->GetLayerJQH(il+1,ixy,obeta,fid);
         if(ixy>=2)continue;
         for(int is=0;is<5;is++)tk_qlqs[il][ixy][is]=tkyjhit.qcl[ixy]._qs[is];
       }
       TrRecHitR* tkhit =trk->GetHitLJ(il+1);
       if(!tkhit)continue;
       tk_edep[il]=tkhit->GetEdepH(2)*1.e-3;//GetEdepH returns keV
       tk_lid[il]=tkhit->GetTkId();
       tk_qls[il]=tkhit->GetQStatus();
       if(il==0)tk_qlsn=tkhit->GetQStatusYJ();
       tk_ql2[il][0]=tkhit->GetQ(0,obeta,frig);//XQ
       tk_ql2[il][1]=tkhit->GetQ(1,obeta,frig);//YQ
     }
    

//---Tk quality
    int algo=1;// 1:Choutko 2:Alcaraz
    int patt=0;// 3:Inner 5:L1I 6:IL9 7:L1IL9
    int refit=3;
    float mass=TrFit::Mproton;
    int chrg=1;

//---Hit Pattern
    tk_rtype=trk->GetRecType();
    tk_hitb[0]=trk->GetBitPatternXYJ(); tk_hitb[1]=trk->GetBitPatternJ();
    for(int ilay=0+1;ilay<9-1;ilay++){
      for(int ixy=0;ixy<2;ixy++){
        if((tk_hitb[ixy]&(1<<ilay))>0)tk_nhiti[ixy]++;
      }
    }
    tk_fzm[0]=trk->GetAdvancedFitCharge(); tk_fzm[1]=trk->GetAdvancedFitMass();
    float tkfz=tk_fzm[0];
    if(ibetah<0){tkfz=-1;}//ibetah=-1 reassigned tkfz

//---
//    int uyh=-1;
    int uyh=0;
    tk_iq=tk_q[1];
    tk_1q=tk_ql[0];
    tk_9q=tk_ql[8];
/*    if(uyh>=0){
      tk_iq=tk_qin[uyh][2];
      tk_1q=tk_qln[uyh][0][2];
      tk_9q=tk_qln[uyh][8][2];
    }*/
    if(uyh>=0&&tk_qin[uyh][2]>8.5){
      tk_iq=tk_qin[uyh][2];
      tk_1q=tk_qln[uyh][0][2];
      tk_9q=tk_qln[uyh][8][2];
    }
    tk_l1l9par  =(tk_iq>0&&tk_1q>0&&tk_9q>0&&(tk_hitb[0]&(1<<0))>0&&(tk_hitb[0]&(1<<8))>0);
    tk_bgpar    =(tk_iq>1.5);
    tk_bhgpar   =(tk_iq>2.5);
    tk_ubpar    =(tk_iq>0&&unbiasetag>0);
    bool tk_mcz1=(tk_iq>0&&!isreal&&fabs(mch)==1);

//---Algo
    for(int ialgo=0;ialgo<3;ialgo++){
      if     (ialgo==0)algo=1;//Choutko
//      else if(ialgo==1)algo=2;//Alcaraz
      else if(ialgo==1)algo=7;//GBL
//      else if(ialgo==2)algo=13;//ChikanianF //without energy loss
      else if(ialgo==2)algo=6;//Kalman
//---Init Track ReFit
      InitTrkReFit(tk_iq,tkfz,refit,mass,chrg,algo);
      for(int ival=0;ival<3;ival++){
        if     (ival==0)refit=refit%10;//PG
        else if(ival==1)refit=20+refit%10;//PG+CIEMAT
        else if(ival==2)refit=30+refit%10;//V6
        int pal=refit/10;
        for(int isp=0;isp<=6;isp++){
          if     (isp==0)patt=0;//Max-Span
          else if(isp==1)patt=3;//Inner
          else if(isp==2)patt=5;//Inner+Layer1
          else if(isp==3)patt=6;//Inner+Layer9
          else if(isp==4)patt=7;//Force Full Span(change to 4)
          else if(isp==5)patt=1;//Inner Up
          else if(isp==6)patt=2;//Inner Down
//--Rigidity
          TrFit trfit;
          if(refit%10>=1&&(algo%10==3||algo%10==6||algo%10==7)){//Chikanian/Kalman/GBL scale
            mfit=-2;
#ifndef _USEFAST_
            if(tk_bgpar||tk_mcz1||tk_ubpar)mfit=-1;
#ifdef _SAVENEGSCALE_
            if(tk_ngpar)mfit=-1;
#endif
#endif
//            if((algo%10==6||algo%10==7)&&(patt==1||patt==2))mfit=-2;
            if(mfit>=-1)mfit=trk->iTrTrackPar(trfit,algo,patt,refit,mass,chrg);
            else if(refit%10<=1)mfit=trk->iTrTrackPar(trfit,algo,patt,(refit/10)*10,mass,chrg);//try with no refit
          }
          else mfit=trk->iTrTrackPar(trfit,algo,patt,refit,mass,chrg);
//--clear
          float nrig=0,nrigp[5]={0};
          float enrig=9999;
          float nchis[3]={9999,9999,9999};
          float ndir0[3]={0,0,-1};
          float ntkres[9][3];
          for(int il=0;il<9;il++)for(int ixy=0;ixy<3;ixy++){ntkres[il][ixy]=9999;}
          if(mfit>=0){
            for(int iz=-1;iz<4;iz++){
              if(isreal&&pal!=3)nrigp[iz+1]=trk->GetCorrectedRigidity(mfit,3,iz);
              else              nrigp[iz+1]=trk->GetRigidity(mfit,iz);
            }
            int izf=-1;
            nrig=nrigp[izf+1];//default
            float bcor=trk->GetBcorr(mfit);//positive
//---Add Magnetic correction
            if(isreal&&bcor==1){//Not applied
              float bcor1=1,bcor2=1;
              int bret1= MagnetVarp::btempcor(bcor1, 0, 1);
              int bret2= MagnetVarp::btempcor(bcor2, 0, 2);
              if      (bret1==0 && bret2==0)bcor=(bcor1+bcor2)/2;
              else if (bret1!=0 && bret2==0)bcor=bcor2;
              nrig*=bcor;
              bcor=-bcor;//put negtive in this case
            }
            if(tk_bcor==0)tk_bcor=bcor;
//----
            enrig=trk->GetErrRinv(mfit);
            nchis[0]=trk->GetNormChisqX(mfit);
            nchis[1]=trk->GetNormChisqY(mfit);
            nchis[2]=trk->GetChisq(mfit);
            for(int ixy=0;ixy<3;ixy++){ndir0[ixy]=trk->GetDir(mfit)[ixy];}
            if((ialgo==0)&&(ival==2)&&(isp<=4)){if(ibetah>=0&&fabs(obeta)>0.4&&tk_nhiti[1]>=5&&tk_iq>0.6&&fabs(nrig)>0.7&&nrig/obeta<0&&nchis[1]<20)tk_ngpar=1;}
            for(int il=0;il<NTKL;il++){
              AMSPoint pnt=trk->GetResidualJ(il+1,mfit);
              for(int ixy=0;ixy<2;ixy++)ntkres[il][ixy]=pnt[ixy];
              if((ialgo==0)&&(ival==2)&&(isp==0)){//V6
                AMSPlaneM plm=trk->GetHitCooLJN(il+1);
                TVector3 m=plm.getMGlobal();
                for(int ixy=0;ixy<3;ixy++)tk_hitc[il][ixy]=m[ixy];
              }
            }
#ifdef _USEKALMANFIT_
            const int nzl=4;
            double pzl[nzl]={195.00, 158.92, 53.06,-69.975};
//            double pzl[nzl]={195.00, 158.92, 53.06,-135.882};
            for(int iz=0;iz<nzl;iz++){
              if(ival!=2)break;
              AMSPoint postr;AMSDir dirtr;
              double rigtr=nrig;
              double zl=pzl[iz];
              if     (iz==0&&!isreal)zl=mevcoo[2];//TOI
              else if(iz==1&&!isreal)zl=mevcoo1[0][2];//L1
              else if(iz==2&&!isreal)zl=mevcoo1[7][2];//L2
              else if(iz==3&&!isreal)zl=mevcoo1[17][2];//LTOF and RICH
//              else if(iz==3&&!isreal)zl=mevcoo1[18][2];//L9
              if(zl==-1000)zl=pzl[iz];
              if(algo%10==6){
                trfit.InterpolateKalman(zl,postr,dirtr,rigtr);
                if(isreal&&pal!=3)rigtr=pev->GetCorrectedRigidity(rigtr,0,3);//rigscale correction
              }
              else if(algo%10==7){
                trfit.InterpolateGbl(zl,postr,dirtr);
              }
              else trk->Interpolate(zl,postr,dirtr,mfit);
              tk_rigidityr[iz][ialgo][isp]=rigtr;
              for(int ic=0;ic<3;ic++){tk_posr[iz][ialgo][isp][ic]=postr[ic];tk_dirr[iz][ialgo][isp][ic]=dirtr[ic];}
            }
#endif
          }
//--Select Copy To Data
          tk_rigidity1[ialgo][ival][isp]=nrig;
          if((algo%10==6)&&(ival==2)){for(int iz=0;iz<4;iz++)tk_rigidityk[iz][isp]=nrigp[iz+1];}//V6
          if(ival==2)tk_erigidity1[ialgo][isp]=enrig;
          for(int ixy=0;ixy<3;ixy++){
            tk_chis1[ialgo][ival][isp][ixy]=nchis[ixy];
            if((ixy==0)&&(isp<4))tk_dir0[ialgo][ival][isp]=ndir0[ixy];
          }
          if((ialgo==0)&&(ival==2)&&(isp==1)){//V6
            for(int il=0;il<NTKL;il++){for(int ixy=0;ixy<2;ixy++)tk_res[il][ixy]=ntkres[il][ixy];}
          }
          if((ialgo==0)&&pal==2&&(isp==0)){//PG+CIEMAT
            for(int uexl=0;uexl<2;uexl++){
              int lay=(uexl==0)?1:9;
              for(int ixy=0;ixy<3;ixy++)tk_cdif[uexl][ixy]=(trk->GetPG_CIEMAT_diff(lay))[ixy];
            }
          }
#ifdef _USELAYERRES_
//----Exclude Residual using this layer(inner-layer)
          if((ialgo==0)&&(ival==2)&&(isp==0)){
            int intkpatt=0;
            for(int il=0+1;il<9-1;il++){if(trk->TestHitLayerJ(il+1))intkpatt+=9*int(pow(10.,il));}
            for(int il=0;il<9;il++){
//              int lpatt=99999990;//(L9)L8L7L6L5L4L3L2(L1)
              int lpatt=intkpatt;
              lpatt=lpatt-((lpatt/int(pow(10.,il)))%10)*int(pow(10.,il));
              int lmfit=trk->iTrTrackPar(algo,lpatt,(refit/10)*10+3,mass,chrg);//force refit
              float nrigi=0,nchisi=9999;
              float nresi[2]={9999,9999};
              if(lmfit>=0){
                if(isreal&&pal!=3)nrigi=trk->GetCorrectedRigidity(lmfit,3);
                else              nrigi=trk->GetRigidity(lmfit);
                nchisi=trk->GetNormChisqY(lmfit);
                AMSPoint pnt=trk->GetResidualJ(il+1,lmfit);
                for(int ixy=0;ixy<2;ixy++)nresi[ixy]=pnt[ixy];
              }
              tk_rigidityi[il]=nrigi;
              tk_chisi[il]=nchisi;
              for(int ixy=0;ixy<2;ixy++)tk_res1[il][ixy]=nresi[ixy];
            }
         }
#endif
//----
      }//end isp
     }//end ival
    }//end ialgo
//----
    trk->RecalcHitCoordinates(0,0,0,chrg);//using default id to restore coordinate
  }
  
  else {

//--Track Pos
    for(int ilay=0;ilay<NTKL;ilay++){
      double length,rigidity;
      AMSPlaneM pl=trk->InterpolateLayerJN(length,rigidity,ilay+1);
      for(int ic=0;ic<3;ic++){
        tk_pos[ilay][ic]=(length<=-1000)?(ic==2?0:9999):pl.getP()[ic];
        tk_dir[ilay][ic]=(length<=-1000)?(ic==2?-1:0):pl.getD()[ic];
      }
      //---
      if(ilay!=0&&ilay!=NTKL-1)continue;
      int uexl=(ilay==0)?0:1;
      int iexl=(ilay==0)?l1i:l9i;
      tk_exql[uexl][0]=tk_exql[uexl][1]=tk_exql[uexl][2]=0;
      for(int ixy=0;ixy<3;ixy++){tk_exqln[0][uexl][ixy]=tk_exqln[1][uexl][ixy]=0;}
      tk_exqls[uexl]=0;
      tk_exdis[uexl]=1000;
      tk_exlid[uexl]=tk_excad[uexl]=0; 
      if(iexl<0||length<=-1000)continue;//interpolation failed
      TrRecHitR* tkhit=pev->pTrRecHit(iexl);
      int mult; double dist,xaddress,dadz[2]={0};
      AMSPlaneM plm=tkhit->FindMult(dist,mult,xaddress,pl.getP(),pl.getD(),-1,tk_iq);
      tk_exdis[uexl]=dist;
      TVector3 dloc=plm.getDLocal();
      for(int ixy=0;ixy<2;ixy++)dadz[ixy]=(dloc[2]!=0)?dloc[ixy]/dloc[2]:0;
      for(int ixy=0;ixy<3;ixy++){
        tk_exql[uexl][ixy]=tkhit->GetQ(ixy,obeta,frig,umass,mult,dadz[0],dadz[1]);//X+Y+XYQ 
        if(tk_q[1]<=1.5)continue;
        tk_exqln[0][uexl][ixy]=tkhit->GetQYJ(ixy,obeta,frig,mult,dadz[0],dadz[1]);
        tk_exqln[1][uexl][ixy]=tkhit->GetQH(ixy,obeta,frig,mult,dadz[0],dadz[1]);
      }
      tk_exqls[uexl]=tkhit->GetQStatus();
      tk_exlid[uexl]=tkhit->GetTkId();
      for(int ixy=0;ixy<2;ixy++){
        TrClusterR *cl=(ixy==0)?tkhit->GetXCluster():tkhit->GetYCluster();
        if(!cl)continue;
        int opt=(isreal)? (TrClusterR::kAsymEta|TrClusterR::kOverflow): TrClusterR::kAsym;
        int cstrip,sstrip;//seed/second strip index
        int noverflow=cl->GetSeedSecondIndex(cstrip,sstrip,opt);//number of ADC overflow strips nearby the seed strip
        tk_excad[uexl]+=int(pow(10000.,ixy))*cl->GetAddress(cstrip);//1000000000*K7+100000000*Mult+10000*y(address)+x(address)
      }
   }

//---Others/Calibration
    for(int il=0;il<9;il++){
       tk_cad[il]=tk_stripx[il]=tk_iso[il][0]=tk_iso[il][1]=-1;
       tk_cof[il][0]=tk_cof[il][1]=0;
       for(int ii=0;ii<11;ii++){
         tk_adc[il][0][ii]=tk_adc[il][1][ii]=0;//5strips[2] largest
       }
       TrRecHitR* tkhit =trk->GetHitLJ(il+1);
       if(!tkhit)continue;
       int k7=tkhit->GetLadder()->IsK7()?1:0;
       int imul=tkhit->GetResolvedMultiplicity();
       if(imul<0||imul>=9)imul=9;//impossible
       tk_cad[il]=(k7*10+imul)*100000000;
       for(int ixy=0;ixy<2;ixy++){
         TrClusterR *cl=(ixy==0)?tkhit->GetXCluster():tkhit->GetYCluster();
         if(!cl)continue;
         int opt=(isreal)? (TrClusterR::kAsymEta|TrClusterR::kOverflow): TrClusterR::kAsym;//QY p/n new asymmetry for ISS
         int cstrip,sstrip;//seed/second strip index
         int noverflow=cl->GetSeedSecondIndex(cstrip,sstrip,opt);//number of ADC overflow strips nearby the seed strip
         tk_cad[il]+=int(pow(10000.,ixy))*cl->GetAddress(cstrip);//1000000000*K7+100000000*Mult+10000*y(address)+x(address)
         tk_cof[il][ixy]=cl->GetLinearCofG(-1,imul,opt,int(tk_iq+0.5));//X+Y linearized center-gravity with respect to the seed [-0.5,0.5]
         if(ixy==0){//X
           int sensor;
           int mcstrip=cl->GetSeedIndex();//muliplicity defined by old-cof algorithm with old-seed
           int sensor_cad=cl->GetSensorAddress(sensor,cstrip,imul,mcstrip,2);//sesor/local sensor-address for seed
           int striptype=cl->GetSeedSecondPairType(imul,mcstrip,opt,2)%1000;//type of seed-second strips location
           int ik=0;
           if     (striptype==102||striptype==202)ik=1;//niters=1(2)/K5/type=2
           else if(striptype==113)                ik=2;//niters=1   /K7/type=3(sensor-address: 64-159)
           else if(striptype==114)                ik=3;//niters=1   /K7/type=4(sensor-address: 0-63,160-223)
           else if(striptype== 14)                ik=4;//niters=0   /K7/type=4(sensor-address: 0-63,160-223)
           tk_stripx[il]=ik+sensor*10+sensor_cad*1000+noverflow*1000000;
         }
//---adc
         double cadc=0,sadc=0,tadc=0,iadc=0;
         for(int ii=0;ii<cl->GetNelem();ii++){
           float nadc=cl->GetSignal(ii,opt);
           if(nadc<=0)continue;
           if(abs(ii-cstrip)<=5){tk_adc[il][ixy][ii-cstrip+5]+=nadc;}
           else                 {tk_adc[il][ixy][ii-cstrip+5<0?0:10]+=nadc;}
           if     (ii==cstrip)cadc=nadc;
           else if(abs(ii-cstrip)==1){if(ii==sstrip)sadc=nadc;else tadc=nadc;}
           else if(abs(ii-cstrip)>=2)iadc+=nadc;//non-position isolate strip sum-up
         }
         if(cadc!=0){tadc=tadc/(cadc+sadc);iadc=iadc/(iadc+cadc+sadc);}
         if(ixy==1){
           tk_iso[il][0]=tadc;
           if(!cl->IsGoodVA(2))tk_iso[il][0]+=10;
           tk_iso[il][1]=iadc;
         }
      }
   }


//---2nd Trk
   tk_oq[0]=tk_oq[1]=-1;
   for(int itr=0;itr<ntrack;itr++){
      if(itr==itrtrack)continue;
      float otkq=pev->pTrTrack(itr)->GetInnerQ();
      if     (otkq>tk_oq[0]){tk_oq[1]=tk_oq[0];tk_oq[0]=otkq;}
      else if(otkq>tk_oq[1]){tk_oq[1]=otkq;}
   }

//--NonTrk Cluster
   float cdis=1.5;//1.5cm
   float cedep=0.05;//0.05MeV
   if(!tk_bgpar)cedep=0.03;//Proton
   float cdism=0.3;//0.3cm
   mfit=trk->Gettrdefaultfit()&(~(TrTrackR::kFitLayer8|TrTrackR::kFitLayer9));//inner track
   for(int il=0;il<NTKL;il++){
     tk_oel[il]=tk_ohitl[il]=0;
     if(mfit<0)continue;
//--TrkPos
     double length,rigidity;
     AMSPlaneM pl=trk->InterpolateLayerJN(length,rigidity,il+1,mfit);//no hit
     if(length<=-1000)continue;//interpolation failed
//--TrkCl
     TrRecHitR  *tkhit=trk->GetHitLJ(il+1);
     TrClusterR *tkcl[2]={0};
     float eltr[2]={0};
     for(int ixy=0;ixy<2;ixy++){
       if(!tkhit)continue;
       TrClusterR *cl=(ixy==0)?tkhit->GetXCluster():tkhit->GetYCluster();
       if(cl){tkcl[ixy]=cl;eltr[ixy]=cl->GetEdep();}
     } 
//--NonTrkCl
     float oeltr[2]={0}; int ohittr[2]={0};
     for(int icl=0;icl<pev->nTrCluster();icl++){
        TrClusterR *cl=pev->pTrCluster(icl);
        if(!cl||cl->GetLayerJ()!=il+1)continue;
        int ixy=cl->GetSide();
        if(cl==tkcl[ixy])continue;//Not Trk-Cluster
        int mult; double dist,xaddress;
        cl->FindMult(dist,mult,xaddress,pl.getP(),pl.getD());
        if(dist>10)continue;//Not >10cm
        float edep=cl->GetEdep();
        if     (eltr[ixy]>0&&dist<cdis&&edep>cedep)oeltr[ixy]+=edep;//TrkHas-Layer
        else if(eltr[ixy]<=0&&dist<cdism&&edep>cedep)oeltr[ixy]+=edep;//TrkMis-Layer
        if(edep>cedep)ohittr[ixy]++;//<10cm
     }
//--
     int uxy=1;//Y
     if(eltr[uxy]>0){tk_oel[il]=oeltr[uxy]/(eltr[uxy]+oeltr[uxy]);}//TrkHas-Layer Positive
     else           {tk_oel[il]=-oeltr[uxy];} //TrkMis-Layer Negtive
     tk_ohitl[il]=ohittr[uxy];
   }

//----IGRF cutoff(Tracker Extraplation)
   if(isreal){
#ifndef _USEFAST_
     if(ibetah>=0&&tk_bgpar&&(tk_ubpar||tk_ngpar)){//unbiase||rig/beta<0
       double length,rigidity;
       AMSPlaneM pl=trk->InterpolateLayerJN(length,rigidity,obeta>0?1:9);//no hit
       AMSDir dirtr(0,0,-1); double cfp=0; 
       if(length>-1000)dirtr.setd(pl.getD()[0],pl.getD()[1],pl.getD()[2]);
       if(dirtr[2]>0)dirtr=dirtr*(-1);
       if(obeta<0)dirtr=dirtr*(-1);
       for(int ipn=0;ipn<2;ipn++){
          int sign=(ipn==0)?-1:1;//negtive(postive) charge
          int cfstat=pev->GetIGRFCutoff(cfp,sign,dirtr);
          if(cfstat==0)cutoffpi[ipn]=cfp;
       }
      }
#endif
    }
    else cutoffpi[0]=cutoffpi[1]=mmom/mch;
//---
  }//end else
   return true;
}

bool AMSAnalysis::Select_Trd(){

    TrdTrackR *trdtrack=pev->pTrdTrack(itrdtrack);
    if(trdtrack){
       nitrdseg=trdtrack->NTrdSegment();
//---TrdTrack Match
       AMSDir dir;AMSPoint pos;
       TrTrackR *trk=pev->pTrTrack(itrtrack);
       if(trk){
         trk->Interpolate(trdtrack->Coo[2],pos,dir);
         distrd[0]=trdtrack->Coo[0]-pos.x();
         distrd[1]=trdtrack->Coo[1]-pos.y();
         AMSDir trddir(trdtrack->Theta,trdtrack->Phi);
         distrd[2]=acos(fabs(trddir.prod(dir)))*180/3.1415926;
       }
    }

//---TrdTrack-Default-Charge
    ChargeR *charge=pev->pCharge(icharge);
    if(charge){
      ChargeSubDR *zTrd=charge->getSubD("AMSChargeTRD");
      if(zTrd){
        trd_rz = TMath::Max(Int_t(zTrd->ChargeI[0]),1);
        trd_rq = zTrd->Q;
        trd_rpz= zTrd->getProb();
      }
    }

 
//-----TrdK
   for(int i=0;i<20;i++){trd_amplk[i]=trd_pathlk[i]=0;}
   TrTrackR *trk=pev->pTrTrack(itrtrack);
   if(!trk||tk_iq<=0)return false;
//---
   float ecalen=0;
   EcalShowerR *show= pev->pEcalShower(ishow);
   if(show){ecalen=show->EnergyE;}
   bool ecalpar=(ecalen>1);
//---
#ifndef _USEFAST_
   if(tk_bgpar||tk_l1l9par||tk_ubpar||ecalpar)trd_statk=-500;
#ifdef _SAVENEGSCALE_
   if(tk_ngpar)trd_statk=-500;
#endif
#endif
   if(trd_statk<=-1000)return false;
//---Track
   int fid=trk->iTrTrackPar(1,0,20); // Get any prefered fit code
   double trdlikr[3]={-1,-1,-1},trdlikre[3]={-1,-1,-1},trdlik[3]={-1,-1,-1}; //LikelihoodRatio :  e/P, e/H, P/H
   int NHits=0;  //number of hits in Likelihood Calculation
   float threshold=15;  //ADC above will be taken into account in Likelihood Calculation
   static unsigned int prtimek=0;//reduce huge output
   static bool prstatk=1;
   if(pev->UTime()==prtimek&&!prstatk){trd_statk=-200;return false;}
   prtimek=pev->UTime();prstatk=1;
   TrdKCluster trdcluster=TrdKCluster(pev,trk,fid);
//   return false;
   int IsReadAlignmentOK=trdcluster.IsReadAlignmentOK; //0: No Alignment, 1: Static Alignment Layer level,  2: Dynamic
   int IsReadCalibOK=trdcluster.IsReadCalibOK;  // 0: No Gain Calibration, 1: Gain Calibration Succeeded
   if(IsReadAlignmentOK==0||IsReadCalibOK==0){prstatk=0;trd_statk=-200;return false;} 
   trd_statk=trdcluster.GetLikelihoodRatio_TrTrack(threshold,trdlikr,NHits);//Track=1
   trdcluster.GetOffTrackHit_TrTrack(trd_onhitk,trd_oampk);
   if(trd_statk!=1){
      trd_statk=trdcluster.GetLikelihoodRatio_TRDRefit(threshold,trdlikr,NHits);
      if(trd_statk==1)trd_statk+=10;//TRD=11
      trdcluster.GetOffTrackHit_TRDRefit(trd_onhitk,trd_oampk);
   }
   trd_nhitk=NHits;//Number of surrounding fired tubes
   if(trd_statk%10==1){//trd_statk=1 || trd_statk=11
      if(ecalpar){
        if(trd_statk==1)trdcluster.GetLikelihoodRatio_TrTrack(threshold,trdlikre,NHits,ecalen,trdlik);
        else            trdcluster.GetLikelihoodRatio_TRDRefit(threshold,trdlikre,NHits,ecalen,trdlik);
      }
      for(int i=0;i<3;i++){trd_lik[i]=trdlikr[i];trd_like[i]=trdlikre[i];}
      AMSPoint pnt; AMSDir dir;
      int x= trdcluster.GetTrTrackExtrapolation(pnt,dir);
      if(x>=0){
       for(int i=0;i<trdcluster.NHits();i++){
         TrdKHit *hit=trdcluster.GetHit(i);
         if(!hit)continue;
         float tube_length=hit->Tube_Track_3DLength(&pnt,&dir);
         if(tube_length<=0.02)continue;
         float ampk = hit->TRDHit_Amp;
         if(ampk<=0)continue;
         if(ampk>32767)ampk=32767;
         trd_ampsk+=(ampk/tube_length);
         trd_pathsk+=tube_length;
         trd_amphitk++;
         int ilay=hit->TRDHit_Layer;
         if(ampk<=trd_amplk[ilay])continue;//MaxADC
         trd_amplk[ilay]=ampk;
         trd_pathlk[ilay]=tube_length;
       }
      }
    }
   ///TRD-Charge
     trdcluster.CalculateTRDCharge(0,tof_betah);
     trd_qk[0]=   trdcluster.GetTRDCharge();
     trd_qk[1]=   trdcluster.GetTRDChargeUpper();
     trd_qk[2]=   trdcluster.GetTRDChargeLower();
     trd_qrmsk[0]=trdcluster.GetTRDChargeError();
     trd_qnhk[0]= trdcluster.GetQNHit();
     trd_qnhk[1]= trdcluster.GetQNHitUpper();
     trd_qnhk[2]= trdcluster.GetQNHitLower();
     trd_ipch=    trdcluster.GetIPChi2();//impact Chis

///--End Fragmentation delta
     trdcluster.CalculateTRDCharge(1,tof_betah);
     trd_qk[3]=   trdcluster.GetTRDCharge();
     trd_qrmsk[1]=trdcluster.GetTRDChargeError();
     trdcluster.CalculateTRDCharge(2,tof_betah);
     trd_qk[4]=   trdcluster.GetTRDCharge();
     trd_qrmsk[2]=trdcluster.GetTRDChargeError();


    return true;
}

bool AMSAnalysis::Select_Tof(){

//---
     int fid=0;
     float frig=(itrtrack>=0)?pev->pTrTrack(itrtrack)->GetRigidity(fid):0;//default
//--
     BetaR *beta=pev->pBeta(ibeta);
     if(beta)tof_beta=beta->Beta; 
	 
	 //--store HT and SHT, which need to rebuild BetaHR, mainly need to rebuild TofClusterHR to have HT/SHT time and error //2024.06.27
	 //now ReBuild at the beginning, TofClusterHR should have HT/SHT already
	 TofRecH::BuildOpt = 0;
	 TofRecH::ReBuild();
	 BetaHR *betahHTSHT = pev->pBetaH(ibetah);
//---
     BetaHR *betah=pev->pBetaH(ibetah);
      
     if(!betah)return false;
    
     //--do TOF MC Tuning //2023.12.30
     //#ifdef _DOBETAMCTUNE_
     #ifdef _DOFINALTUNE_ //2024.01.14
     TofBetaPar betaparSmear=betah->gTofBetaPar();
     if (!isreal)
     {
	   if (mch==2)
	   {
         TofMCPar::MCtuneDT=-48.68;
         TofMCPar::MCtuneST=18.80;
	   }
	   else if (mch==3)
	   {
         TofMCPar::MCtuneDT=-32.58;
         TofMCPar::MCtuneST=9.55;
	   }
	   else if (mch==4)
	   {
         TofMCPar::MCtuneDT=-23.36;
         TofMCPar::MCtuneST=12.91;
	   }
	   else if (mch==5)
	   {
         TofMCPar::MCtuneDT=-10.46;
         TofMCPar::MCtuneST=11.89;
	   }
	   else if (mch==6)
	   {
         TofMCPar::MCtuneDT=0;
         TofMCPar::MCtuneST=12.57;
	   }
	   else if (mch==7)
	   {
         //TofMCPar::MCtuneDT=-17.08;
         //TofMCPar::MCtuneST=8.50;
		 //based on data after LTofQ cut
         TofMCPar::MCtuneDT=-9.16;
         TofMCPar::MCtuneST=11.51;
	   }
	   else if (mch==8)
	   {
         TofMCPar::MCtuneDT=-11.12;
         TofMCPar::MCtuneST=13.09;
	   }
	   
       //int iBetaMCTune = betah->DoMCtune();
	   int iBetaMCTune = betah->DoMCtune(betaparSmear, 0); //retrieve fit parameter
       if (iBetaMCTune<0) cout << " !!!!!! Select_Particle: fail to tune MC TOF, BetaHR::DoMCtune()=" << iBetaMCTune << endl;
     }
     #endif //_DOFINALTUNE_
     //#endif //_DOBETAMCTUNE_
	 
     tof_btype=betah->GetBuildType();
     tof_bpatt=betah->GetBetaPattern();

     tof_betah=betah->GetBeta();
     tof_betahmc=(!isreal)?betah->GetMCBeta(1):tof_betah;
     tof_ebetah=betah->GetEBetaV();
     tof_chist=betah->GetChi2T();
     tof_chisc=betah->GetChi2C();
//---Normalized
     TofBetaPar normbetapar=betah->GetNormTofBetaPar();
     tof_chist_n=normbetapar.Chi2T; 
     tof_chisc_n=normbetapar.Chi2C;
     tof_ebetah_n=normbetapar.InvErrBeta; 
//----
     tof_hsumh=betah->GetSumHit();
     tof_hsumhu=betah->GetUseHit();

//---TOF Layer
     tof_pass=0;
     AMSPoint pnt; AMSDir dir; AMSPoint tofpnt;
     float disedge;double utm;
	 //pointer to TofCluster in each layer //2024.06.27
     TofClusterHR *tof_clh[4] = {NULL, NULL, NULL, NULL};
     for(int ilay=0;ilay<NTOFL;ilay++){
//---Clear
        tof_barid[ilay]=tof_tl[ilay]=-1;
        tof_oq[ilay][0]=tof_oq[ilay][1]=0;
        tof_ob[ilay][0]=tof_ob[ilay][1]=-1;
        tof_sbit[ilay][0]=tof_sbit[ilay][1]=0;
        tof_tres[ilay]=9999; //YC //2024.01.03
        tof_cres[ilay]=9999;
		tof_leng[ilay]=0; //YC //2024.06.02
		//2024.06.27
		tof_edep[ilay]=0;
		//2024.07.02
		tof_tl_HT[ilay]=tof_tl_SHT[ilay]=-1;
		for (int iter=0; iter<10; iter++)
		{
			tof_deltaBetah[iter][ilay]=tof_deltaBetah_HT[iter][ilay]=tof_deltaBetah_SHT[iter][ilay]=-9999;
		}
		//2024.07.22
		tof_deltaBetah_el[ilay]=-9999;
		//2024.09.16
		for (int iside=0; iside<2; iside++)
		{
			nHTime[ilay][iside] = nSHTime[ilay][iside] = -1;
			tof_rtime[ilay][iside] = -1;
			tof_rtime_ht[ilay][iside] = -1;
			tof_rtime_sht[ilay][iside] = -1;
			tof_Aadc[ilay][iside] = -1;
		} //iside
//---All
        tof_nclhl[ilay]=betah->GetAllFireHL(ilay);
        betah->TInterpolate (TOFGeom::GetMeanZ(ilay),pnt,dir,utm);
        for(int ic=0;ic<3;ic++)
        {
          tof_pos[ilay][ic]=pnt[ic];
          tof_dir[ilay][ic]=dir[ic];
        }
        int tk_pass=ParticleR::IsPassTOF(ilay,pnt,dir,tofpnt,disedge);
        if(tk_pass>=0){tof_pass+=int(pow(10.,3-ilay));}
//----TOFCL
        if(!betah->TestExistHL(ilay))continue;
        TofClusterHR *tofclh=betah->GetClusterHL(ilay);
        tof_barid[ilay]=tofclh->Bar;
        if(tofclh->IsGoodTime()){
           tof_tl[ilay]=tofclh->Time;
           tof_etl[ilay]=tofclh->ETime;
           tof_tres[ilay]=betah->GetTResidual(ilay); //YC //2024.01.03
           tof_cres[ilay]=betah->GetCResidual(ilay,TOFGeom::Proj[ilay]);
           #ifdef _DOFINALTUNE_ //2024.01.17
           tof_tl[ilay]=betaparSmear.Time[ilay];
           tof_tres[ilay]=betaparSmear.TResidual[ilay];
           tof_cres[ilay]=betaparSmear.CResidual[ilay][TOFGeom::Proj[ilay]];
           #endif //_DOFINALTUNE_
		   //deposit energy //2024.06.27
		   tof_edep[ilay]=tofclh->GetEdep();
		   
		   //--store tof_rtime //2024.09.17
		   for (int iside=0; iside<2; iside++)
		   {
			   tof_rtime[ilay][iside] = tofclh->Rtime[iside];
			  //  tof_rtime_ht[ilay][iside] = tofclh->Rhtime[iside];
			  //  tof_rtime_sht[ilay][iside] = tofclh->Rshtime[iside];
			   tof_Aadc[ilay][iside] = tofclh->Aadc[iside];
		   }
        }
		//--tof path length //2024.06.02
		if (itrtrack>=0)
		{
			AMSPoint tfcoo(tofclh->Coo);
			AMSDir tkdir(0,0,1);
			AMSPoint tkcoo;
			double theta, phi, sleng;
			pev->pTrTrack(itrtrack)->interpolate(tfcoo,tkdir,tkcoo,theta,phi,sleng);
			tof_leng[ilay] = sleng;
		}
        for(int is=0;is<2;is++)tof_sbit[ilay][is]=tofclh->SideBitPat[is];
		
		//2024.06.27
		tof_clh[ilay] = tofclh;
     }
	 //beta using only S2S3, which should not influenced by energy loss //2024.06.27
	 tof_betah_S2S3=0;
	 tof_betah_S1S4=0;
	 tof_betah_S1S2S3S4=0;
	 for (int iter=0; iter<10; iter++) tof_betah_iter[iter]=0; //2024.07.02
     //TofClusterHR *tof_clhS2S3[4] = {NULL, tof_clh[1], tof_clh[2], NULL};
	 //S2S3
	 if (tof_clh[1]!=NULL && tof_clh[2]!=NULL)
	 {
		 //TofBetaPar fitparS2S3;
		 //TofRecH::GetHead()->TofTrackFit(tof_clhS2S3, fitparS2S3);
		 //use BetaRefit
		 TofBetaPar fitparS2S3 = betah->gTofBetaPar();
		 betah->BetaReFit(fitparS2S3,110,1,0); //using ilay=1 ilay=2 only
		 tof_betah_S2S3=fitparS2S3.Beta;
	 }
	 //S1S4
	 if (tof_clh[0]!=NULL && tof_clh[3]!=NULL)
	 {
		 //use BetaRefit
		 TofBetaPar fitparS1S4 = betah->gTofBetaPar();
		 betah->BetaReFit(fitparS1S4,1001,1,0); //using ilay=0 ilay=3 only
		 tof_betah_S1S4=fitparS1S4.Beta;
	 }
	 //add S1S2 and S3S4 as well, the difference due to energy loss should be largest for them //2024.07.02
	 tof_betah_S1S2=0;
	 if (tof_clh[0]!=NULL && tof_clh[1]!=NULL)
	 {
		 TofBetaPar fitparS1S2 = betah->gTofBetaPar();
		 /*cout << " S1S2: " << endl; //XXX
		 cout << " UseHit before refit: " << fitparS1S2.UseHit << endl; //XXX*/
		 //betah->BetaReFit(fitparS1S2,11,1,0); //using ilay=1 ilay=2 only
		 betah->BetaReFit(fitparS1S2,1100,1,0); //using ilay=3 ilay=4 only
		 tof_betah_S1S2=fitparS1S2.Beta;
		 //XXX
		 /*cout << " UseHit after refit: " << fitparS1S2.UseHit << endl; //XXX
		 for (int ilay=0; ilay<4; ilay++)
		 {
			 cout << fitparS1S2.Time[ilay] << endl;
		 } //ilay*/
	 }
	 tof_betah_S3S4=0;
	 if (tof_clh[2]!=NULL && tof_clh[3]!=NULL)
	 {
		 TofBetaPar fitparS3S4 = betah->gTofBetaPar();
		 /*cout << " S3S4: " << endl; //XXX
		 cout << " UseHit before refit: " << fitparS3S4.UseHit << endl; //XXX*/
		 //betah->BetaReFit(fitparS3S4,1100,1,0); //using ilay=3 ilay=4 only
		 betah->BetaReFit(fitparS3S4,11,1,0); //using ilay=3 ilay=4 only
		 tof_betah_S3S4=fitparS3S4.Beta;
		 //XXX
		 /*cout << " UseHit after refit: " << fitparS3S4.UseHit << endl; //XXX
		 for (int ilay=0; ilay<4; ilay++)
		 {
			 cout << fitparS3S4.Time[ilay] << endl;
		 } //ilay*/
	 }
	 
	 //--same value as TofRech::cvel, which is a protected number and cannot be used //2024.07.02
	 float cvel(29.979);
	 //TofBetaPar fitparS1S2S3S4;
	 //TofRecH::GetHead()->TofTrackFit(tof_clh, fitparS1S2S3S4);
	 //use BetaReFit
	 TofBetaPar fitparS1S2S3S4 = betah->gTofBetaPar();
	 //betah->BetaReFit(fitparS1S2S3S4,-1,1,0); //using all layer, should be the same as tof_betah
	 betah->BetaReFit(fitparS1S2S3S4,1111,1,0); //using all layer, should be the same as tof_betah
	 tof_betah_S1S2S3S4=fitparS1S2S3S4.Beta;
	 //cout << "tof_betah=" << tof_betah << " betah-S1S2S3S4=" << tof_betah-tof_betah_S1S2S3S4 << " S2S3=" << tof_betah_S2S3 << " S1S4=" << tof_betah_S1S4 << endl; //XXX
	 /*//XXX
	 cout << " tof_betah=" << tof_betah << endl;
	 cout << " S1S2S3S4=" << tof_betah_S1S2S3S4 << endl;
	 cout << " S2S3=" << tof_betah_S2S3 << endl;
	 cout << " S1S2=" << tof_betah_S1S2 << endl;
	 cout << " S3S4=" << tof_betah_S3S4 << endl;
	 cout << tof_clh[0] << " " << tof_clh[1] << " " << tof_clh[2] << " " << tof_clh[3] << endl;
	 //XXX*/
	 //calculate the beta change //2024.07.02 //will put into member function of TofBetaPar //TODO
	 for (int iter=0; iter<10; iter++)
	 {
		 for (int ilay=0; ilay<NTOFL; ilay++)
		 {
			 if (ilay==0)
			 {
				 tof_deltaBetah[iter][0] = -fitparS1S2S3S4.Beta + (fitparS1S2S3S4.Len[0]-fitparS1S2S3S4.Len[1])/(cvel*fitparS1S2S3S4.Time[0]-cvel*fitparS1S2S3S4.T0-fitparS1S2S3S4.Len[1]/fitparS1S2S3S4.Beta);
			 }
			 else if (ilay==3)
			 {
				 tof_deltaBetah[iter][3] =  fitparS1S2S3S4.Beta - (fitparS1S2S3S4.Len[3]-fitparS1S2S3S4.Len[2])/(cvel*fitparS1S2S3S4.Time[3]-cvel*fitparS1S2S3S4.T0-fitparS1S2S3S4.Len[2]/fitparS1S2S3S4.Beta);
			 }
			 else continue;
		 } //ilay
		 //--refit the beta with corrected t1 and t4
		 for (int ilay=0; ilay<NTOFL; ilay++)
		 {
			 if (ilay==0)
			 {
				 fitparS1S2S3S4.Time[0] = fitparS1S2S3S4.Time[0]+(1/fitparS1S2S3S4.Beta -1/(fitparS1S2S3S4.Beta+tof_deltaBetah[iter][0])) * (fitparS1S2S3S4.Len[0] - fitparS1S2S3S4.Len[1])/cvel;
			 }
			 else if (ilay==3)
			 {
				 fitparS1S2S3S4.Time[3] = fitparS1S2S3S4.Time[3]+(1/fitparS1S2S3S4.Beta -1/(fitparS1S2S3S4.Beta-tof_deltaBetah[iter][3])) * (fitparS1S2S3S4.Len[3] - fitparS1S2S3S4.Len[2])/cvel;
			 }
			 else continue;
		 } //ilay
		 betah->BetaReFit(fitparS1S2S3S4,-1,1,0);
		 tof_betah_iter[iter]=fitparS1S2S3S4.Beta;
	 } //iter
	 //cout << " betah=" << tof_betah_S1S2S3S4 << " betah_iter10=" << tof_betah_iter[9] << endl; //XXX
	 //cout << " betah chage: " << tof_deltaBetah[9][0] << " " << tof_deltaBetah[9][3] << endl; //XXX
	 
	 //--store HT and SHT, which need to rebuild BetaHR, mainly need to rebuild TofClusterHR to have HT/SHT time and error //2024.06.27
	 //cout << " rebuild TOFBetah=" << betahHTSHT->GetBeta() << " diff between tof_betah=" << betahHTSHT->GetBeta()-tof_betah << endl; //XXX
	 //use HT
	 tof_betah_HT=0;
	 for (int iter=0; iter<10; iter++) tof_betah_iter_HT[iter]=0; //2024.07.02
	 TofBetaPar fitparHT = betah->gTofBetaPar();
	 for (int ilay=0; ilay<NTOFL; ilay++)
	 {
		 if(!betah->TestExistHL(ilay))continue;
		 if (tof_clh[ilay]==NULL) continue;
		//  fitparHT.Time[ilay] = tof_clh[ilay]->HTime;
		//  fitparHT.ETime[ilay] = tof_clh[ilay]->EHTime;
		 //2024.07.02
		 tof_tl_HT[ilay] = fitparHT.Time[ilay];
		 //2024.09.16
		//  for (int iside=0; iside<2; iside++) nHTime[ilay][iside] = tof_clh[ilay]->nHTime[iside];
	 } //ilay
	 betah->BetaReFit(fitparHT,-1,1,0);
	 tof_betah_HT=fitparHT.Beta;
	 //2024.07.04
	 tof_chist_HT=fitparHT.Chi2T;
	 tof_chisc_HT=fitparHT.Chi2C;
     for(int ilay=0;ilay<NTOFL;ilay++)
	 {
        tof_tres_HT[ilay]=9999;
        tof_cres_HT[ilay]=9999;
		
        if(!betah->TestExistHL(ilay))continue;
        TofClusterHR *tofclh=betah->GetClusterHL(ilay);
        if(tofclh->IsGoodTime())
		{
           tof_tres_HT[ilay]=fitparHT.TResidual[ilay];
           tof_cres_HT[ilay]=fitparHT.CResidual[ilay][TOFGeom::Proj[ilay]];
		}
	 } //ilay
	 //normlized chis
	 betah->BetaReFit(fitparHT,-1,2,0);
	 tof_chist_n_HT=fitparHT.Chi2T;
	 tof_chisc_n_HT=fitparHT.Chi2C;
	 
	 //calculate the beta change //2024.07.02 //will put into member function of TofBetaPar //TODO
	 for (int iter=0; iter<10; iter++)
	 {
		 for (int ilay=0; ilay<NTOFL; ilay++)
		 {
			 if (ilay==0)
			 {
				 //tof_deltaBetah_HT[0] = -fitparHT.Beta + (fitparHT.Len[0]-fitparHT.Len[1])/(TofRecH::cvel*fitparHT.Time[0]-TofRecH::cvel*fitparHT.T0-fitparHT.Len[1]/fitparHT.Beta);
				 tof_deltaBetah_HT[iter][0] = -fitparHT.Beta + (fitparHT.Len[0]-fitparHT.Len[1])/(cvel*fitparHT.Time[0]-cvel*fitparHT.T0-fitparHT.Len[1]/fitparHT.Beta);
			 }
			 else if (ilay==3)
			 {
				 //tof_deltaBetah_HT[3] =  fitparHT.Beta - (fitparHT.Len[3]-fitparHT.Len[2])/(TofRecH::cvel*fitparHT.Time[3]-TofRecH::cvel*fitparHT.T0-fitparHT.Len[2]/fitparHT.Beta);
				 tof_deltaBetah_HT[iter][3] =  fitparHT.Beta - (fitparHT.Len[3]-fitparHT.Len[2])/(cvel*fitparHT.Time[3]-cvel*fitparHT.T0-fitparHT.Len[2]/fitparHT.Beta);
			 }
			 else continue;
		 } //ilay
		 //--refit the beta with corrected t1 and t4
		 for (int ilay=0; ilay<NTOFL; ilay++)
		 {
			 if (ilay==0)
			 {
				 //fitparHT.Time[0] = fitparHT.Time[0]+(1/fitparHT.Beta -1/(fitparHT.Beta+tof_deltaBetah_HT[0])) * (fitparHT.Len[0] - fitparHT.Len[1])/TofRecH::cvel;
				 fitparHT.Time[0] = fitparHT.Time[0]+(1/fitparHT.Beta -1/(fitparHT.Beta+tof_deltaBetah_HT[iter][0])) * (fitparHT.Len[0] - fitparHT.Len[1])/cvel;
			 }
			 else if (ilay==3)
			 {
				 //fitparHT.Time[3] = fitparHT.Time[3]+(1/fitparHT.Beta -1/(fitparHT.Beta-tof_deltaBetah_HT[3])) * (fitparHT.Len[3] - fitparHT.Len[2])/TofRecH::cvel;
				 fitparHT.Time[3] = fitparHT.Time[3]+(1/fitparHT.Beta -1/(fitparHT.Beta-tof_deltaBetah_HT[iter][3])) * (fitparHT.Len[3] - fitparHT.Len[2])/cvel;
			 }
			 else continue;
		 } //ilay
		 betah->BetaReFit(fitparHT,-1,1,0);
		 tof_betah_iter_HT[iter]=fitparHT.Beta;
	 } //iter
	 //cout << " betah_HT=" << tof_betah_HT << " betah_iter10_HT=" << tof_betah_iter_HT[9] << endl; //XXX
	 //cout << " betah_HT chage: " << tof_deltaBetah_HT[9][0] << " " << tof_deltaBetah_HT[9][3] << endl; //XXX
	 
	 //use SHT
	 tof_betah_SHT=0;
	 for (int iter=0; iter<10; iter++) tof_betah_iter_SHT[iter]=0;
	 TofBetaPar fitparSHT = betah->gTofBetaPar();
	 for (int ilay=0; ilay<NTOFL; ilay++)
	 {
		 if(!betah->TestExistHL(ilay))continue;
		 if (tof_clh[ilay]==NULL) continue;
		//  fitparSHT.Time[ilay] = tof_clh[ilay]->SHTime;
		//  fitparSHT.ETime[ilay] = tof_clh[ilay]->ESHTime;
		 //2024.07.02
		 tof_tl_SHT[ilay] = fitparSHT.Time[ilay];
		 //2024.09.16
		//  for (int iside=0; iside<2; iside++) nSHTime[ilay][iside] = tof_clh[ilay]->nSHTime[iside];
	 } //ilay
	 betah->BetaReFit(fitparSHT,-1,1,0);
	 tof_betah_SHT=fitparSHT.Beta;
	 //2024.07.04
	 tof_chist_SHT=fitparSHT.Chi2T;
	 tof_chisc_SHT=fitparSHT.Chi2C;
     for(int ilay=0;ilay<NTOFL;ilay++)
	 {
        tof_tres_SHT[ilay]=9999;
        tof_cres_SHT[ilay]=9999;
		
        if(!betah->TestExistHL(ilay))continue;
        TofClusterHR *tofclh=betah->GetClusterHL(ilay);
        if(tofclh->IsGoodTime())
		{
           tof_tres_SHT[ilay]=fitparSHT.TResidual[ilay];
           tof_cres_SHT[ilay]=fitparSHT.CResidual[ilay][TOFGeom::Proj[ilay]];
		}
	 } //ilay
	 //normlized chis
	 betah->BetaReFit(fitparSHT,-1,2,0);
	 tof_chist_n_SHT=fitparSHT.Chi2T;
	 tof_chisc_n_SHT=fitparSHT.Chi2C;
	 
	 //calculate the beta change //2024.07.02 //will put into member function of TofBetaPar //TODO
	 for (int iter=0; iter<10; iter++)
	 {
		 for (int ilay=0; ilay<NTOFL; ilay++)
		 {
			 if (ilay==0)
			 {
				 //tof_deltaBetah_SHT[0] = -fitparSHT.Beta + (fitparSHT.Len[0]-fitparSHT.Len[1])/(TofRecH::cvel*fitparSHT.Time[0]-TofRecH::cvel*fitparSHT.T0-fitparSHT.Len[1]/fitparSHT.Beta);
				 tof_deltaBetah_SHT[iter][0] = -fitparSHT.Beta + (fitparSHT.Len[0]-fitparSHT.Len[1])/(cvel*fitparSHT.Time[0]-cvel*fitparSHT.T0-fitparSHT.Len[1]/fitparSHT.Beta);
			 }
			 else if (ilay==3)
			 {
				 //tof_deltaBetah_SHT[3] =  fitparSHT.Beta - (fitparSHT.Len[3]-fitparSHT.Len[2])/(TofRecH::cvel*fitparSHT.Time[3]-TofRecH::cvel*fitparSHT.T0-fitparSHT.Len[2]/fitparSHT.Beta);
				 tof_deltaBetah_SHT[iter][3] =  fitparSHT.Beta - (fitparSHT.Len[3]-fitparSHT.Len[2])/(cvel*fitparSHT.Time[3]-cvel*fitparSHT.T0-fitparSHT.Len[2]/fitparSHT.Beta);
			 }
			 else continue;
		 } //ilay
		 //--refit the beta with corrected t1 and t4
		 for (int ilay=0; ilay<NTOFL; ilay++)
		 {
			 if (ilay==0)
			 {
				 //fitparSHT.Time[0] = fitparSHT.Time[0]+(1/fitparSHT.Beta -1/(fitparSHT.Beta+tof_deltaBetah_SHT[0])) * (fitparSHT.Len[0] - fitparSHT.Len[1])/TofRecH::cvel;
				 fitparSHT.Time[0] = fitparSHT.Time[0]+(1/fitparSHT.Beta -1/(fitparSHT.Beta+tof_deltaBetah_SHT[iter][0])) * (fitparSHT.Len[0] - fitparSHT.Len[1])/cvel;
			 }
			 else if (ilay==3)
			 {
				 //fitparSHT.Time[3] = fitparSHT.Time[3]+(1/fitparSHT.Beta -1/(fitparSHT.Beta-tof_deltaBetah_SHT[3])) * (fitparSHT.Len[3] - fitparSHT.Len[2])/TofRecH::cvel;
				 fitparSHT.Time[3] = fitparSHT.Time[3]+(1/fitparSHT.Beta -1/(fitparSHT.Beta-tof_deltaBetah_SHT[iter][3])) * (fitparSHT.Len[3] - fitparSHT.Len[2])/cvel;
			 }
			 else continue;
		 } //ilay
		 betah->BetaReFit(fitparSHT,-1,1,0);
		 tof_betah_iter_SHT[iter]=fitparSHT.Beta;
	 } //iter
	 //cout << " betah_SHT=" << tof_betah_SHT << " betah_iter10_SHT=" << tof_betah_iter_SHT[9] << endl; //XXX
	 //cout << " betah_SHT chage: " << tof_deltaBetah_SHT[9][0] << " " << tof_deltaBetah_SHT[9][3] << endl; //XXX
	 #ifdef TEST
	 cout << " tof_betah_HT=" << tof_betah_HT << " tof_betah_SHT=" << tof_betah_SHT << endl; //XXX
	 for (int il=0; il<4; il++)
	 {
		 if(!betah->TestExistHL(il))continue;
		 if (tof_clh[il]==NULL) continue;
		 cout << "il=" << il+1 << " " << tof_tl_HT[il] << " " << tof_clh[il]->Shtime[0] << " " << tof_clh[il]->Shtime[1] << endl;
	 } //il
	 #endif //TEST
	 
	 //combine beta fit using LT, HT and SHT //2024.07.28
	 tof_betah_LTHTSHT=0;
	 TofBetaPar fitparLTHTSHT = betah->gTofBetaPar();
	 double dtof_tl_ht[4], dtof_tl_sht[4];
	 for (int ilay=0; ilay<4; ilay++) dtof_tl_ht[ilay] = tof_tl_HT[ilay], dtof_tl_sht[ilay]=tof_tl_SHT[ilay];
	//  betah->BetaReFit_HTSHT(fitparLTHTSHT,dtof_tl_ht,dtof_tl_sht,-1,1,0);
	 tof_betah_LTHTSHT=fitparLTHTSHT.Beta;
	 /*cout << " ****** " << endl;
	 cout << " tof beta:	" << tof_betah << endl;
	 cout << " HT:			" << tof_betah_HT << endl;
	 cout << " SHT:			" << tof_betah_SHT << endl;
	 cout << " LT+HT+SHT:	" << tof_betah_LTHTSHT << endl;*/
	 
	 //--calculate the beta using my TOF Time Calibration //2024.08.25
	 #ifdef MYTOFTIMECALIB
	 //--set time and bar
	 for (int imycalibver=0; imycalibver<NMyTOFTimeCalibVer; imycalibver++)
	 {
		 tofcalib[imycalibver]->InitTimeADC();
		 for (int ilay=0; ilay<4; ilay++)
		 {
			 if(!betah->TestExistHL(ilay))continue;
			 if (tof_clh[ilay]==NULL) continue;
			 for (int iside=0; iside<2; iside++)
			 {
				 //raw time
				 tofcalib[imycalibver]->dRawT[0][ilay][iside] = tof_clh[ilay]->Rtime[iside];
				 //tofcalib[imycalibver]->dRawT[1][ilay][iside] = tof_clh[ilay]->Shtime[iside];
				 //tofcalib[imycalibver]->dRawT[2][ilay][iside] = tof_clh[ilay]->Sshtime[iside];
				 //correct raw time
				 tofcalib[imycalibver]->dRawT[1][ilay][iside] = tof_clh[ilay]->Rhtime[iside];
				 tofcalib[imycalibver]->dRawT[2][ilay][iside] = tof_clh[ilay]->Rshtime[iside];
				 
				 //adc
				 tofcalib[imycalibver]->dADC[ilay][iside] = tof_clh[ilay]->Aadc[iside];
			 } //iside
			 
			 //bar id
			 tofcalib[imycalibver]->ibar[ilay] = tof_clh[ilay]->Bar;
		 } //ilay
		 
		 //cout << " ****** run=" << run << " event=" << event << " ****** " << endl; //XXX
		 //cout << " tof_betah=" << tof_betah << endl;
		 for (int ilthtsht=0; ilthtsht<3; ilthtsht++)
		 {
			 tofcalib[imycalibver]->CalTime(ilthtsht);
			 //cout << endl; //XXX
			 //cout << " ilthtsht=" << ilthtsht << endl;
			 
			 /*for (int ilay=0; ilay<4; ilay++)
			 {
				 if(!betah->TestExistHL(ilay))continue;
				 if (tof_clh[ilay]==NULL) continue;
				 cout << " L" << ilay+1 << " RTime: " << tofcalib[imycalibver]->dRawT[ilthtsht][ilay][0] << " " << tofcalib[imycalibver]->dRawT[ilthtsht][ilay][1] << endl; //XXX
			 } //ilay*/
			 
			 for (int ifix=0; ifix<2; ifix++)
			 {
				 //cout << endl; //XXX
				 //cout << " ifix=" << ifix << endl;
				 tof_betah_MyCalib[imycalibver][ilthtsht][ifix] = 0;
				 tof_chist_MyCalib[imycalibver][ilthtsht][ifix] = 0;
				 tof_chisc_MyCalib[imycalibver][ilthtsht][ifix] = 0;
				 //tof_chist_n_MyCalib[imycalibver][ilthtsht][ifix] = 0;
				 //tof_chisc_n_MyCalib[imycalibver][ilthtsht][ifix] = 0;
				 TofBetaPar fitparMyCalib = betah->gTofBetaPar();
				 
				 for (int ilay=0; ilay<4; ilay++)
				 {
					 tof_tl_MyCalib[imycalibver][ilthtsht][ifix][ilay] = 0;
					 
					 if(!betah->TestExistHL(ilay))continue;
					 if (tof_clh[ilay]==NULL) continue;
#ifdef TEST
					 cout << " L" << ilay+1 << " Cal T: " << tofcalib[imycalibver]->dTimeCalib[ilthtsht][ifix][ilay][0] << " " << tofcalib[imycalibver]->dTimeCalib[ilthtsht][ifix][ilay][1] << " Comb T: " << tofcalib[imycalibver]->dTimeCom[ilthtsht][ifix][ilay] << endl; //XXX
					 cout << " Raw T: " << tofcalib[imycalibver]->dRawT[ilthtsht][ilay][0] << " " << tofcalib[imycalibver]->dRawT[ilthtsht][ilay][1] << endl; //XXX
#endif //TEST
					 fitparMyCalib.Time[ilay] = tofcalib[imycalibver]->dTimeCom[ilthtsht][ifix][ilay];
					 tof_tl_MyCalib[imycalibver][ilthtsht][ifix][ilay] = tofcalib[imycalibver]->dTimeCom[ilthtsht][ifix][ilay];
				 } //ilay
				 betah->BetaReFit(fitparMyCalib,-1,1,0);
				 tof_betah_MyCalib[imycalibver][ilthtsht][ifix] = fitparMyCalib.Beta;
				 tof_chist_MyCalib[imycalibver][ilthtsht][ifix] = fitparMyCalib.Chi2T;
				 tof_chisc_MyCalib[imycalibver][ilthtsht][ifix] = fitparMyCalib.Chi2C;
				 for (int ilay=0; ilay<4; ilay++)
				 {
					 tof_cres_MyCalib[imycalibver][ilthtsht][ifix][ilay] = 9999;
					 tof_tres_MyCalib[imycalibver][ilthtsht][ifix][ilay] = 9999;
					 
					 if(!betah->TestExistHL(ilay))continue;
					 if (tof_clh[ilay]==NULL) continue;
					 
					 tof_tres_MyCalib[imycalibver][ilthtsht][ifix][ilay] = fitparMyCalib.TResidual[ilay];
					 tof_cres_MyCalib[imycalibver][ilthtsht][ifix][ilay] = fitparMyCalib.CResidual[ilay][TOFGeom::Proj[ilay]];
				 } //ilay
				 //--normalized chis
				 betah->BetaReFit(fitparMyCalib,-1,2,0);
				 tof_chist_n_MyCalib[imycalibver][ilthtsht][ifix] = fitparMyCalib.Chi2T;
				 tof_chisc_n_MyCalib[imycalibver][ilthtsht][ifix] = fitparMyCalib.Chi2C;
				#ifdef TEST
	    	 	if (ilthtsht==1) cout << " MyTofCalib ilthtsht=" << ilthtsht << " ifix=" << ifix << " " << tof_betah_MyCalib[imycalibver][ilthtsht][ifix] << endl; //XXX
				#endif //TEST
			 } //ifix
		 } //ilthtsht
	 } //imycalibver
	 #endif //MYTOFTIMECALIB
	 
	 //--tof using my calibration for LT, with new parameters using pr, he4, o16, pr+he4//2024.09.05
	 #ifdef COMBINE_Z_TOFTIMECALIB
	 //--calculate the time
	 //cout << " ****** run=" << run << " event=" << event << " ****** " << endl; //XXX
	 //cout << " tof_betah=" << tof_betah << endl;
	 
	 for (int icalibver=0; icalibver<NTOFTimeCalibVer; icalibver++)
	 {
		 #ifdef TEST
		 cout << endl; //XXX
		 cout << " icalibver=" << icalibver << endl;
		 #endif //TEST
		 
		 //--set time and bar
		 tofcalibCombZ[icalibver].InitTimeADC();
		 for (int ilay=0; ilay<4; ilay++)
		 {
			 if(!betah->TestExistHL(ilay))continue;
			 if (tof_clh[ilay]==NULL) continue;
			 for (int iside=0; iside<2; iside++)
			 {
				 //raw time
				 tofcalibCombZ[icalibver].dRawT[0][ilay][iside] = tof_clh[ilay]->Rtime[iside];
				 //tofcalibCombZ[icalibver].dRawT[1][ilay][iside] = tof_clh[ilay]->Shtime[iside];
				 //tofcalibCombZ[icalibver].dRawT[2][ilay][iside] = tof_clh[ilay]->Sshtime[iside];
				 //correct raw time
				 tofcalibCombZ[icalibver].dRawT[1][ilay][iside] = tof_clh[ilay]->Rhtime[iside];
				 tofcalibCombZ[icalibver].dRawT[2][ilay][iside] = tof_clh[ilay]->Rshtime[iside];
				 
				 //adc
				 tofcalibCombZ[icalibver].dADC[ilay][iside] = tof_clh[ilay]->Aadc[iside];
			 } //iside
			 
			 //bar id
			 tofcalibCombZ[icalibver].ibar[ilay] = tof_clh[ilay]->Bar;
		 } //ilay
		 
		 for (int ilthtsht=0; ilthtsht<1; ilthtsht++)
		 {
			 tofcalibCombZ[icalibver].CalTime(ilthtsht);
			 //cout << endl; //XXX
			 //cout << " ilthtsht=" << ilthtsht << endl;
			 
			 /*for (int ilay=0; ilay<4; ilay++)
			 {
				 if(!betah->TestExistHL(ilay))continue;
				 if (tof_clh[ilay]==NULL) continue;
				 cout << " L" << ilay+1 << " RTime: " << tofcalibCombZ[icalibver].dRawT[ilthtsht][ilay][0] << " " << tofcalibCombZ[icalibver].dRawT[ilthtsht][ilay][1] << endl; //XXX
			 } //ilay*/
			 
			 for (int ifix=0; ifix<2; ifix++)
			 {
		 		#ifdef TEST
				 cout << endl; //XXX
				 cout << " ifix=" << ifix << endl;
				#endif //TEST
				 tof_betah_CombZCalib[icalibver][ifix] = 0;
				 tof_chist_CombZCalib[icalibver][ifix] = 0;
				 tof_chisc_CombZCalib[icalibver][ifix] = 0;
				 TofBetaPar fitparCombZCalib = betah->gTofBetaPar();
				 
				 for (int ilay=0; ilay<4; ilay++)
				 {
					 tof_tl_CombZCalib[icalibver][ifix][ilay] = 0;
					 
					 if(!betah->TestExistHL(ilay))continue;
					 if (tof_clh[ilay]==NULL) continue;
		 			#ifdef TEST
					 cout << " L" << ilay+1 << " Cal T: " << tofcalibCombZ[icalibver].dTimeCalib[ilthtsht][ifix][ilay][0] << " " << tofcalibCombZ[icalibver].dTimeCalib[ilthtsht][ifix][ilay][1] << " Comb T: " << tofcalibCombZ[icalibver].dTimeCom[ilthtsht][ifix][ilay] << endl; //XXX
		 			#endif //TEST
					 fitparCombZCalib.Time[ilay] = tofcalibCombZ[icalibver].dTimeCom[ilthtsht][ifix][ilay];
					 tof_tl_CombZCalib[icalibver][ifix][ilay] = tofcalibCombZ[icalibver].dTimeCom[ilthtsht][ifix][ilay];
				 } //ilay
				 betah->BetaReFit(fitparCombZCalib,-1,1,0);
				 tof_betah_CombZCalib[icalibver][ifix] = fitparCombZCalib.Beta;
				 tof_chist_CombZCalib[icalibver][ifix] = fitparCombZCalib.Chi2T;
				 tof_chisc_CombZCalib[icalibver][ifix] = fitparCombZCalib.Chi2C;
				 for (int ilay=0; ilay<4; ilay++)
				 {
					 tof_cres_CombZCalib[icalibver][ifix][ilay] = 9999;
					 tof_tres_CombZCalib[icalibver][ifix][ilay] = 9999;
					 
					 if(!betah->TestExistHL(ilay))continue;
					 if (tof_clh[ilay]==NULL) continue;
					 
					 tof_tres_CombZCalib[icalibver][ifix][ilay] = fitparCombZCalib.TResidual[ilay];
					 tof_cres_CombZCalib[icalibver][ifix][ilay] = fitparCombZCalib.CResidual[ilay][TOFGeom::Proj[ilay]];
				 } //ilay
				 
				 //--normalized chis
				 betah->BetaReFit(fitparCombZCalib,-1,2,0);
				 tof_chist_n_CombZCalib[icalibver][ifix] = fitparCombZCalib.Chi2T;
				 tof_chisc_n_CombZCalib[icalibver][ifix] = fitparCombZCalib.Chi2C;
				#ifdef TEST
				 cout << " icalibver=" << icalibver << " ifix=" << ifix << " " << tof_betah_CombZCalib[icalibver][ifix] << endl; //XXX
				#endif //TEST
			 } //ifix
		 } //ilthtsht
	 } //icalibver
	 #endif //COMBINE_Z_TOFTIMECALIB
	 
	 //--tof using my calibration for LT, with new parameters using pr, he4, o16, pr+he4//2024.09.05
	 #ifdef TOFBETA_WITHOUT_CAL
	 //--calculate the time
	 //cout << " ****** run=" << run << " event=" << event << " ****** " << endl; //XXX
	 //cout << " tof_betah=" << tof_betah << endl;
	 
	//--set time and bar
	tofcalibNull.InitTimeADC();
	for (int ilay=0; ilay<4; ilay++)
	{
	    if(!betah->TestExistHL(ilay))continue;
	    if (tof_clh[ilay]==NULL) continue;
	    for (int iside=0; iside<2; iside++)
	    {
	   	 //raw time
	   	 tofcalibNull.dRawT[0][ilay][iside] = tof_clh[ilay]->Rtime[iside];
	   	 //tofcalibNull.dRawT[1][ilay][iside] = tof_clh[ilay]->Shtime[iside];
	   	 //tofcalibNull.dRawT[2][ilay][iside] = tof_clh[ilay]->Sshtime[iside];
		 //correct raw time
	   	 tofcalibNull.dRawT[1][ilay][iside] = tof_clh[ilay]->Rhtime[iside];
	   	 tofcalibNull.dRawT[2][ilay][iside] = tof_clh[ilay]->Rshtime[iside];
	   	 
	   	 //adc
	   	 tofcalibNull.dADC[ilay][iside] = tof_clh[ilay]->Aadc[iside];
	    } //iside
	    
	    //bar id
	    tofcalibNull.ibar[ilay] = tof_clh[ilay]->Bar;
	} //ilay
	
	for (int ilthtsht=0; ilthtsht<3; ilthtsht++)
	{
	    tofcalibNull.CalTime(ilthtsht);
	    //cout << endl; //XXX
	    //cout << " ilthtsht=" << ilthtsht << endl;
	    
	    /*for (int ilay=0; ilay<4; ilay++)
	    {
	   	 if(!betah->TestExistHL(ilay))continue;
	   	 if (tof_clh[ilay]==NULL) continue;
	   	 cout << " L" << ilay+1 << " RTime: " << tofcalibNull.dRawT[ilthtsht][ilay][0] << " " << tofcalibNull.dRawT[ilthtsht][ilay][1] << endl; //XXX
	    } //ilay*/
	    
	    for (int ifix=0; ifix<1; ifix++)
	    {
		#ifdef TEST
	   	 cout << endl; //XXX
	   	 cout << " ifix=" << ifix << endl;
	   	#endif //TEST
	   	 tof_betah_NoCalib[ilthtsht] = 0;
	   	 tof_chist_NoCalib[ilthtsht] = 0;
	   	 tof_chisc_NoCalib[ilthtsht] = 0;
	   	 TofBetaPar fitparCombZCalib = betah->gTofBetaPar();
	   	 
	   	 for (int ilay=0; ilay<4; ilay++)
	   	 {
	   		 tof_tl_NoCalib[ilthtsht][ilay] = 0;
	   		 
	   		 if(!betah->TestExistHL(ilay))continue;
	   		 if (tof_clh[ilay]==NULL) continue;
				#ifdef TEST
	   		 cout << " L" << ilay+1 << " Cal T: " << tofcalibNull.dTimeCalib[ilthtsht][ifix][ilay][0] << " " << tofcalibNull.dTimeCalib[ilthtsht][ifix][ilay][1] << " Comb T: " << tofcalibNull.dTimeCom[ilthtsht][ifix][ilay] << endl; //XXX
				#endif //TEST
	   		 fitparCombZCalib.Time[ilay] = tofcalibNull.dTimeCom[ilthtsht][ifix][ilay];
	   		 tof_tl_NoCalib[ilthtsht][ilay] = tofcalibNull.dTimeCom[ilthtsht][ifix][ilay];
	   	 } //ilay
	   	 betah->BetaReFit(fitparCombZCalib,-1,1,0);
	   	 tof_betah_NoCalib[ilthtsht] = fitparCombZCalib.Beta;
	   	 tof_chist_NoCalib[ilthtsht] = fitparCombZCalib.Chi2T;
	   	 tof_chisc_NoCalib[ilthtsht] = fitparCombZCalib.Chi2C;
	   	 for (int ilay=0; ilay<4; ilay++)
	   	 {
	   		 tof_cres_NoCalib[ilthtsht][ilay] = 9999;
	   		 tof_tres_NoCalib[ilthtsht][ilay] = 9999;
	   		 
	   		 if(!betah->TestExistHL(ilay))continue;
	   		 if (tof_clh[ilay]==NULL) continue;
	   		 
	   		 tof_tres_NoCalib[ilthtsht][ilay] = fitparCombZCalib.TResidual[ilay];
	   		 tof_cres_NoCalib[ilthtsht][ilay] = fitparCombZCalib.CResidual[ilay][TOFGeom::Proj[ilay]];
	   	 } //ilay
	   	 
	   	 //--normalized chis
	   	 betah->BetaReFit(fitparCombZCalib,-1,2,0);
	   	 tof_chist_n_NoCalib[ilthtsht] = fitparCombZCalib.Chi2T;
	   	 tof_chisc_n_NoCalib[ilthtsht] = fitparCombZCalib.Chi2C;
	   	#ifdef TEST
	   	 cout << " ilthtsht=" << ilthtsht << " ifix=" << ifix << " " << tof_betah_NoCalib[ilthtsht] << endl; //XXX
		 if (ilthtsht==1)
		 {
			 for (int il=0; il<4; il++)
			 {
				 if(!betah->TestExistHL(il))continue;
				 if (tof_clh[il]==NULL) continue;
				 cout << "il=" << il+1 << " " << tof_tl_NoCalib[ilthtsht][il] << " " << -(tof_clh[il]->Shtime[0]+tof_clh[il]->Shtime[1])/2. << " " << tof_clh[il]->Shtime[0] << " " << tof_clh[il]->Shtime[1] << endl;
			 }
		 }
	   	#endif //TEST
	    } //ifix
	} //ilthtsht
	 #endif //TOFBETA_WITHOUT_CAL
	 
	 //-- //2024.10.25
	 #ifdef NEWTOFTIMECALIB
	 //--set time and bar
	 for (int iqcalib=0; iqcalib<NZ_NewCalib; iqcalib++)
	 {
		for (int iminchis=0; iminchis<NMinChis_NewCalib; iminchis++)
		{
			for (int iexcludeedge=0; iexcludeedge<NExcludeEdge_NewCalib; iexcludeedge++)
			{
				 //tofcalibNewCalib[iqcalib][iminchis][iexcludeedge].InitTimeADC();
				 //2024.11.04
				 //TofCalib *tc = NULL;
				 //if (iqcalib==0 || iqcalib==1) tc = &tofcalibNewCalib[iqcalib][iminchis][iexcludeedge];
				 //else tc = &tofcalibNewCalib[((int)mch)-1+1][iminchis][iexcludeedge];
				 //2024.11.16
				 TofCalib *tc = &tofcalibNewCalib[iqcalib][iminchis][iexcludeedge];
				 tc->InitTimeADC();
				 
				 for (int ilay=0; ilay<4; ilay++)
				 {
					 if(!betah->TestExistHL(ilay))continue;
					 if (tof_clh[ilay]==NULL) continue;
					 for (int iside=0; iside<2; iside++)
					 {
						 //raw time
						 tc->dRawT[0][ilay][iside] = tof_clh[ilay]->Rtime[iside];
						 //tc->dRawT[1][ilay][iside] = tof_clh[ilay]->Shtime[iside];
						 //tc->dRawT[2][ilay][iside] = tof_clh[ilay]->Sshtime[iside];
						 //correct raw time
						 tc->dRawT[1][ilay][iside] = tof_clh[ilay]->Rhtime[iside];
						 tc->dRawT[2][ilay][iside] = tof_clh[ilay]->Rshtime[iside];
						 
						 //adc
						 tc->dADC[ilay][iside] = tof_clh[ilay]->Aadc[iside];
					 } //iside
					 
					 //bar id
					 tc->ibar[ilay] = tof_clh[ilay]->Bar;
				 } //ilay
				 
				 //cout << " ****** run=" << run << " event=" << event << " ****** " << endl; //XXX
				 //cout << " tof_betah=" << tof_betah << endl;
				 for (int ilthtsht=0; ilthtsht<3; ilthtsht++)
				 {
					 tc->CalTime(ilthtsht);
					 //cout << endl; //XXX
					 //cout << " ilthtsht=" << ilthtsht << endl;
					 
					 /*for (int ilay=0; ilay<4; ilay++)
					 {
						 if(!betah->TestExistHL(ilay))continue;
						 if (tof_clh[ilay]==NULL) continue;
						 cout << " L" << ilay+1 << " RTime: " << tc->dRawT[ilthtsht][ilay][0] << " " << tc->dRawT[ilthtsht][ilay][1] << endl; //XXX
					 } //ilay*/
					 
					 for (int ifix=0; ifix<NFix_NewCalib; ifix++)
					 {
						 //cout << endl; //XXX
						 //cout << " ifix=" << ifix << endl;
						 tof_betah_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix] = 0;
						 tof_chist_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix] = 0;
						 tof_chisc_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix] = 0;
						 //tof_chist_n_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix] = 0;
						 //tof_chisc_n_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix] = 0;
						 TofBetaPar fitparMyCalib = betah->gTofBetaPar();
						 
						 for (int ilay=0; ilay<4; ilay++)
						 {
							 tof_tl_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix][ilay] = 0;
							 
							 if(!betah->TestExistHL(ilay))continue;
							 if (tof_clh[ilay]==NULL) continue;
#ifdef TEST
							 cout << " L" << ilay+1 << " Cal T: " << tc->dTimeCalib[ilthtsht][ifix][ilay][0] << " " << tc->dTimeCalib[ilthtsht][ifix][ilay][1] << " Comb T: " << tc->dTimeCom[ilthtsht][ifix][ilay] << endl; //XXX
							 cout << " Raw T: " << tc->dRawT[ilthtsht][ilay][0] << " " << tc->dRawT[ilthtsht][ilay][1] << endl; //XXX
#endif //TEST
							 fitparMyCalib.Time[ilay] = tc->dTimeCom[ilthtsht][ifix][ilay];
							 tof_tl_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix][ilay] = tc->dTimeCom[ilthtsht][ifix][ilay];
						 } //ilay
						 betah->BetaReFit(fitparMyCalib,-1,1,0);
						 tof_betah_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix] = fitparMyCalib.Beta;
						 tof_chist_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix] = fitparMyCalib.Chi2T;
						 tof_chisc_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix] = fitparMyCalib.Chi2C;
						 for (int ilay=0; ilay<4; ilay++)
						 {
							 tof_cres_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix][ilay] = 9999;
							 tof_tres_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix][ilay] = 9999;
							 
							 if(!betah->TestExistHL(ilay))continue;
							 if (tof_clh[ilay]==NULL) continue;
							 
							 tof_tres_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix][ilay] = fitparMyCalib.TResidual[ilay];
							 tof_cres_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix][ilay] = fitparMyCalib.CResidual[ilay][TOFGeom::Proj[ilay]];
						 } //ilay
						 //--normalized chis
						 betah->BetaReFit(fitparMyCalib,-1,2,0);
						 tof_chist_n_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix] = fitparMyCalib.Chi2T;
						 tof_chisc_n_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix] = fitparMyCalib.Chi2C;
						#ifdef TEST
						if (ilthtsht==1) cout << " MyTofCalib ilthtsht=" << ilthtsht << " ifix=" << ifix << " " << tof_betah_NewCalib[iqcalib][iminchis][iexcludeedge][ilthtsht][ifix] << endl; //XXX
						#endif //TEST
					 } //ifix
				 } //ilthtsht
			} //iexcludeedge
		} //iminchis
	 } //iqcalib
	 #endif //NEWTOFTIMECALIB
	 
//---TOF ClusterH
    for(int iclh=0;iclh<ntofclh;iclh++){
       TofClusterHR *tofclh=pev->pTofClusterH(iclh);
       int ilay=tofclh->Layer;
       if(tofclh->Bar==tof_barid[ilay])continue;
       if(tofclh->GetQSignal()>tof_oq[ilay][0]){//>Max NonBetaH Cluster
          tof_oq[ilay][1]=tof_oq[ilay][0];
          tof_oq[ilay][0]=tofclh->GetQSignal();
          tof_ob[ilay][1]=tof_ob[ilay][0];
          tof_ob[ilay][0]=tofclh->Bar;
       }
       else if(tofclh->GetQSignal()>tof_oq[ilay][1]){//>Second Max NonBetaH Cluster
          tof_oq[ilay][1]=tofclh->GetQSignal();
          tof_ob[ilay][1]=tofclh->Bar;
       }
    }

//----
   float obetac=tof_betah;
//--TOFCharge All
   tof_qs=0;
   TofChargeHR *tofhc=betah->pTofChargeH();
   int nlay=0;
   tof_zr=tofhc->GetZ(nlay,tof_pzr);

//---TOFQL new
   tofhc->ReFit(obetac,TofChargeHR::DefaultQOptIonW,frig); 
   tof_q[0]=tofhc->GetQ(nlay,tof_qrms[0],1111);//All
   tof_q[1]=tofhc->GetQ(nlay,tof_qrms[1],-2);//Truncate 
///---Z All
   tof_z   =tofhc->GetZ(nlay,tof_pz);
   tof_zud[0]=tofhc->GetZ(nlay,tof_pud[0],0,1100);
   tof_zud[1]=tofhc->GetZ(nlay,tof_pud[1],0,11);
//---TOFQL
   for(int ilay=0;ilay<NTOFL;ilay++){
//---Clear
     for(int ity=0;ity<3;ity++){tof_qlr[ity][ilay]=0;}
     tof_ql[ilay]=tof_ql1[ilay]=0;
//---BetaH Check
     tof_ql1[ilay]=betah->GetQL(ilay,2,TofClusterHR::DefaultQOptIonW,111111,obetac,frig);
     tof_qlr[0][ilay]=betah->GetQL(ilay,2,TofClusterHR::DefaultQOptIonW,111111,1,25);//no-beta(no-rigidity) correction
//---QL
     if(!tofhc->TestExistHL(ilay))continue;
     tof_ql[ilay] =tofhc->GetQL(ilay,22);//New TOFL
     tof_qlr[1][ilay]=tofhc->GetQL(ilay,1);//Anode
     tof_qlr[2][ilay]=tofhc->GetQL(ilay,0);//Dynode
     if(tofhc->IsGoodQPathL(ilay)){tof_qs+=int(pow(10.,3-ilay));}
   }

//---
   int ncls[4];
   tof_oncl=pev->GetNTofClustersInTime(betah, ncls); 

//---Calibration only
#ifdef _USETKADC_
  for(int ic=0;ic<3;ic++){
    if     (ic==0)tofhc->ReFit(obetac,(TofChargeHR::DefaultQOptIonW&(~TofRecH::kBetaCor)&(~TofRecH::kRigidityCor)),frig);
    else if(ic==1)tofhc->ReFit(obetac,(TofChargeHR::DefaultQOptIonW|TofRecH::kBetaRigV1Cor),frig);
    else          tofhc->ReFit(obetac,(TofChargeHR::DefaultQOptIonW|TofRecH::kBetaRigV2Cor),frig);
    for(int ilay=0;ilay<NTOFL;ilay++){
      for(int ity=0;ity<2;ity++){tof_qlc[ic*2+ity][ilay]=0;}
      if(ic==2){
        for(int is=0;is<2;is++){for(int ipm=0;ipm<4;ipm++)tof_qlpm[ilay][is][ipm]=0;}
      }
      if(!tofhc->TestExistHL(ilay))continue;
      tof_qlc[ic*2+0][ilay]=tofhc->GetQL(ilay,1);//Anode
      tof_qlc[ic*2+1][ilay]=tofhc->GetQL(ilay,0);//Dynode
//------
      if(ic==2){
        for(int is=0;is<2;is++){
          tof_qlpm[ilay][is][0]=betah->GetQLPM(ilay,1,is,0,(TofChargeHR::DefaultQOptIonW|TofRecH::kBetaRigV2Cor),frig);//Anode
          for(int ipm=0;ipm<3;ipm++){
             tof_qlpm[ilay][is][ipm+1]=betah->GetQLPM(ilay,0,is,ipm,(TofChargeHR::DefaultQOptIonW|TofRecH::kBetaRigV2Cor),frig);//Dynode
          }
        }
      }
//------
    }
  }
#endif
	 
	 //--calculate the beta change from energy loss directly, using average mass, and then refit //2024.07.22
	 //--a few issue: 1) what energy loss to be used, 2) how to determind the average mass, according to TOF Z or Tracker Z, 3) iteration might still be needed as the delta beta depends on beta as well
	 //--might causing some bias, as cannot be done for beta>1 as gamma will be meaningless. How to keep the gaussian behavior? or it's fine as this will only be important for small beta events? TODO
	 //if some hits are missing, should be OK as at the end we will require 4/4, but if two layers in upper or lower have hits, than the time can be corrected //2024.07.23
	 TofBetaPar fitpar_el= betah->gTofBetaPar();
	 double dGamma = 1./sqrt(1-fitpar_el.Beta*fitpar_el.Beta);
	 const double dMass[30] = {1., 4., 6.5, 8., 10.7, 12, 14.5, 16, 19., 20., 23., 24., 27., 28., 31., 32., 35., 36., 39., 40., 0., 0., 0., 0., 0., 56., 0., 58., 0., 0.};
	 const double dMassUnit = 0.93149;
	 double dTimeCor[4] = {0};
	 bool bBothHitUD[2];
	 bBothHitUD[0] = tof_clh[0]!=NULL && tof_clh[1]!=NULL;
	 bBothHitUD[1] = tof_clh[2]!=NULL && tof_clh[3]!=NULL;
	 
	 if (bBothHitUD[0] || bBothHitUD[1])
	 {
		 for (int ilay=0; ilay<NTOFL; ilay++)
		 {
			 tof_deltaBetah_el[ilay]=0;
			 tof_deltaT_el[ilay]=0;
			 if (ilay==0 || ilay==3)
			 {
				 //should have both measurement for Upper or Lower
				 if (ilay==0 && (!bBothHitUD[0])) continue;
				 else if (ilay==3 && (!bBothHitUD[1])) continue;
				 
				 //tof_deltaBetah_el[ilay] = 1./(fitpar_el.Beta*pow(dGamma,3)*dMass[tof_z-1])*tof_edep[ilay]*1e-3;
				 tof_deltaBetah_el[ilay] = 1./(fitpar_el.Beta*pow(dGamma,3)*dMass[tof_z-1]*dMassUnit)*tof_edep[ilay]*1e-3; //2024.07.29
				 if (ilay==3) tof_deltaBetah_el[ilay] = -tof_deltaBetah_el[ilay];
				 //--time correction due to energy loss //2024.07.29 move back to this order
				 int ilay0, ilay1;
				 if (ilay==0) ilay0=0, ilay1=1;
				 else if (ilay==3) ilay0=3, ilay1=2;
				 /*//should always be later layer - previous layer
				 if (ilay==0) ilay0=1, ilay1=0;
				 else if (ilay==3) ilay0=3, ilay1=2;*/
				 dTimeCor[ilay] = (1/fitpar_el.Beta-1/(fitpar_el.Beta+tof_deltaBetah_el[ilay])) * (fitpar_el.Len[ilay0] - fitpar_el.Len[ilay1])/cvel;
				 //rethink about the path lenght order, as for S1 and S4 are not the same, the old order should be correct //2024.07.29
				 fitpar_el.Time[ilay] += dTimeCor[ilay];
				 tof_deltaT_el[ilay] = dTimeCor[ilay];
			 }
			 else continue;
		 } //ilay
		 betah->BetaReFit(fitpar_el,-1,1,0);
		 tof_betah_el=fitpar_el.Beta;
		 //store time and coo residual as well //2024.07.29
		 for (int ilay=0; ilay<4; ilay++)
		 {
           tof_tres_el[ilay]=fitpar_el.TResidual[ilay];
           tof_cres_el[ilay]=fitpar_el.CResidual[ilay][TOFGeom::Proj[ilay]];
		 } //ilay
		 /*//XXX
		 if (tof_betah<1)
		 {
			 cout << " *** " << endl;
			 cout << " energy loss cor:" << endl;
			 cout << " tofz=" << tof_z << " dMass=" << dMass[tof_z-1] << endl;
			 cout << tof_betah << " " << tof_betah_el << endl;
			 cout << "beta change: " << tof_deltaBetah_el[0] << " " << tof_deltaBetah_el[3] << endl;
			 cout << " Gamma=" << dGamma << endl;
			 cout << " EL [GeV]: " << tof_edep[0]*1e-3 << " " << tof_edep[3]*1e-3 << endl;
			 cout << " Time Cor: " << endl;
			 cout << dTimeCor[0] << " " << dTimeCor[3] << endl;
			 cout << " Path length: " << endl;
			 for (int ilay=0; ilay<4; ilay++)
			 {
				 cout << "ilay=" << ilay << " " << fitpar_el.Len[ilay] << endl;
			 } //ilay
			 cout << " Time: " << endl;
			 for (int ilay=0; ilay<4; ilay++)
			 {
				 cout << "ilay=" << ilay << " " << fitpar_el.Time[ilay] << endl;
			 } //ilay
			 cout << " Time 0: " << fitpar_el.T0 << endl;
		 }
		 //XXX*/
	 }
 
  //--do TOF MC Tuning //2023.12.30
  #ifdef _DOBETAMCTUNE_
  if (!isreal)
  {
    #ifdef _NEWTOFMCTUNE_ //2024.01.30
	TofMCPar::MCtuneDS=4;
    #endif //_NEWTOFMCTUNE_
	
	for (int iDT=0; iDT<NDT; iDT++)
	{
		TofMCPar::MCtuneDT=MCTuneDT[iDT];
		//for (int iST=-NST; iST<NST+1; iST++)
		for (int iST=0; iST<NST; iST++)
		{
			//TofMCPar::MCtuneST=iST*MCTuneST;
			TofMCPar::MCtuneST=MCTuneST[iST];
			
			if (TofMCPar::MCtuneDT==0 && TofMCPar::MCtuneST==0) continue;
			//cout << "TofMCPar::MCtuneDS=" << TofMCPar::MCtuneDS << endl; //XXX
			//cout << "TofMCPar::MCtuneDT=" << TofMCPar::MCtuneDT << " TofMCPar::MCtuneST=" << TofMCPar::MCtuneST << endl; //XXX
			TofBetaPar betaparSmear=betah->gTofBetaPar();
			/*cout << "Before tune" << endl;
			cout << betaparSmear.Time[0] << " " << betaparSmear.Time[1] << " " << betaparSmear.Time[2] << " " << betaparSmear.Time[3] << endl;
			cout << betaparSmear.Beta << endl;
			cout << betaparSmear.TResidual[0] << endl;*/
			int iBetaMCTune = betah->DoMCtune(betaparSmear, 0);
			if (iBetaMCTune<0)
			{
				cout << " !!!!!! Select_Particle: fail to tune MC TOF, BetaHR::DoMCtune()=" << iBetaMCTune << endl;
				continue;
			}
			/*cout << "After tune " << iDT << " " << iST+NST << endl;
			cout << betaparSmear.Time[0] << " " << betaparSmear.Time[1] << " " << betaparSmear.Time[2] << " " << betaparSmear.Time[3] << endl;
			cout << betaparSmear.Beta << endl;
			cout << betaparSmear.TResidual[0] << endl;*/
			//--store result
            //tof_betah_Smear[iDT][iST+NST]=betaparSmear.Beta;
            //tof_chist_Smear[iDT][iST+NST]=betaparSmear.Chi2T;
            //tof_chisc_Smear[iDT][iST+NST]=betaparSmear.Chi2C;
            tof_betah_Smear[iDT][iST]=betaparSmear.Beta;
            tof_chist_Smear[iDT][iST]=betaparSmear.Chi2T;
            tof_chisc_Smear[iDT][iST]=betaparSmear.Chi2C;
//---Normalized
            //TofBetaPar normbetapar=betah->GetNormTofBetaPar();
            TofBetaPar normbetapar=betaparSmear;
			betah->BetaReFit(normbetapar,-1,2);
            //tof_chist_n_Smear[iDT][iST+NST]=normbetapar.Chi2T; 
            //tof_chisc_n_Smear[iDT][iST+NST]=normbetapar.Chi2C;
            tof_chist_n_Smear[iDT][iST]=normbetapar.Chi2T; 
            tof_chisc_n_Smear[iDT][iST]=normbetapar.Chi2C;
//---TOF Layer
            for(int ilay=0;ilay<NTOFL;ilay++)
			{
//---Clear
			  //tof_tl_Smear[iDT][iST+NST][ilay]=-1;
              //tof_tres_Smear[iDT][iST+NST][ilay]=9999;
              //tof_cres_Smear[iDT][iST+NST][ilay]=9999;
			  tof_tl_Smear[iDT][iST][ilay]=-1;
              tof_tres_Smear[iDT][iST][ilay]=9999;
              tof_cres_Smear[iDT][iST][ilay]=9999;
//----TOFCL
              if(!betah->TestExistHL(ilay))continue;
              TofClusterHR *tofclh=betah->GetClusterHL(ilay);
              //tof_barid[ilay]=tofclh->Bar;
              if(tofclh->IsGoodTime())
			  {
                 //tof_tl_Smear[iDT][iST+NST][ilay]=betaparSmear.Time[ilay];
                 //tof_tres_Smear[iDT][iST+NST][ilay]=betaparSmear.TResidual[ilay];
                 //tof_cres_Smear[iDT][iST+NST][ilay]=betaparSmear.CResidual[ilay][TOFGeom::Proj[ilay]];
                 tof_tl_Smear[iDT][iST][ilay]=betaparSmear.Time[ilay];
                 tof_tres_Smear[iDT][iST][ilay]=betaparSmear.TResidual[ilay];
                 tof_cres_Smear[iDT][iST][ilay]=betaparSmear.CResidual[ilay][TOFGeom::Proj[ilay]];
				 /*cout << "iDT=" << iDT << " iST+NST=" << iST+NST << " ilay=" << ilay << endl; //XXX
				 cout << "tof_tl_Smear=" << tof_tl_Smear[iDT][iST+NST][ilay] << endl;
				 cout << "tof_tres_Smear=" << tof_tres_Smear[iDT][iST+NST][ilay] << endl;
				 cout << "tof_cres_Smear=" << tof_cres_Smear[iDT][iST+NST][ilay] << endl;*/
              }
           }
		} //iST
	} //iDT
  }
  #endif //_DOBETAMCTUNE_

  return true; 
}

bool AMSAnalysis::Select_Rich(){
//---
    RichRingR *ring=pev->pRichRing(irich);
    if(ring){
      rich_itrtrack=ring->iTrTrack();
      rich_beta[0]=ring->getBeta();
      rich_beta[1]=ring->Beta; 
      rich_beta[2]=ring->BetaRefit;
      rich_ebeta=ring->ErrorBeta;
      rich_ebeta1=ring->getBetaError();
      rich_pb=ring->getProb();
      rich_udis=ring->UDist;
      rich_hit=ring->getHits();// 
      rich_used=ring->getUsedHits();
      rich_usedm=ring->getReflectedHits();
      rich_q[0]=ring->getCharge2Estimate();
      rich_q[1]=ring->NpExp>0?sqrt(ring->NpCol/ring->NpExp):-1;
      rich_width=ring->getWidth();
      rich_npe[0]=ring->getPhotoelectrons();
      rich_npe[1]=ring->getExpectedPhotoElectrons();
      rich_npe[2]=RichHitR::getCollectedPhotoElectrons();
      rich_stat=ring->Status;
      rich_good=ring->IsGood();
      rich_clean=ring->IsClean();
      rich_NaF=ring->IsNaF();
      rich_pmt=ring->getPMTs();
      rich_pmt1=RichHitR::getPMTs()-RichHitR::getPMTs(false);
      rich_PMTChargeConsistency=ring->getPMTChargeConsistency();
      rich_BetaConsistency=ring->getBetaConsistency();
      rich_pos[0]=ring->getTrackEmissionPoint()[0];
      rich_pos[1]=ring->getTrackEmissionPoint()[1];
      rich_pos[2]=ring->getTrackEmissionPoint()[2];
      rich_theta=ring->getTrackTheta();
      rich_phi=ring->getTrackPhi();
      rich_n[0]=ring->getIndexUsed();
      rich_n[1]=ring->getRawIndexUsed();
      rich_tile=ring->getTileIndex();
      rich_distb=ring->DistanceTileBorder();
      rich_cstat=ring->PmtCorrectionsFailed();
    }
//-----
     RichRingBR *ringb=pev->pRichRingB(irichb);
     if(ringb){
       rich_betap=ringb->Beta;
       rich_likp= ringb->Likelihood;
       rich_pbp=  ringb->ProbKolm;
       rich_usedp=ringb->Used;
       rich_qp=   ringb->ChargeRec; 
       rich_statp=ringb->Status;
       for(int iz=0;iz<3;iz++){rich_pzp[iz]=ringb->ChargeProb[iz];}
     }
   
//---Rich-Default-Charge
   ChargeR *charge=pev->pCharge(icharge);
   if(charge){
     ChargeSubDR *zRing=charge->getSubD("AMSChargeRich");
     if(zRing){
       rich_rz = TMath::Max(Int_t(zRing->ChargeI[0]),1);
       rich_rq = zRing->Q;
       rich_rpz= zRing->getProb();
     }
   }
 
 
  return true;
}

bool AMSAnalysis::Select_Ecal(){
   
   EcalShowerR *show= pev->pEcalShower(ishow);
   ecal_q=0;
   for(int ilay=0;ilay<NECALL;ilay++){ecal_el[ilay]=ecal_eh[ilay]=ecal_nhit[ilay]=0;}
   if(show){
     ecal_stat=show->Status;
     ecal_en[0]=show->EnergyE;
     ecal_en[1]=show->EnergyA;
     ecal_en[2]=show->EnergyD;
     ecal_en[3]=show->EnergyC;
     ecal_rearl=show->RearLeakPI;
     ecal_q=show->EcalChargeEstimator();
     for(int i2dcl=0;i2dcl<show->NEcal2DCluster();i2dcl++){//2dCl
       Ecal2DClusterR *e2dcl=show->pEcal2DCluster(i2dcl);
       for(int icl=0;icl<e2dcl->NEcalCluster();icl++){//1dCl
         EcalClusterR *ecl=e2dcl->pEcalCluster(icl);
         for(int ihit=0;ihit<ecl->NEcalHit();ihit++){//Nhit
           EcalHitR *ehit=ecl->pEcalHit(ihit);
           float edep=ehit->Edep;
           int ilay=ehit->Plane; 
           ecal_nhit[ilay]++;
           ecal_el[ilay]+=edep;
           if(edep>ecal_eh[ilay]){ecal_eh[ilay]=edep;}
         }
       }
     }
#ifndef _USEFAST_
     bool ecalpar=(ecal_en[0]>1);
     bool ecalbepar=(ecal_en[0]>1);//>1GV
     if(ecalpar){
        ecal_bdt[0]=show->GetEcalBDT();
        ecal_bdt[1]=show->GetEcalBDT(7,1);
//       ecal_bdt[1]=show->EcalStandaloneEstimatorV2();
//       if(ecalbepar)ecal_bdt[2]=show->GetEcalBDTCHI2();
     }
#endif
     ecal_enc[0]=show->GetCorrectedEnergy(2,2,0);//e+-/ps/oves(default)
     ecal_enc[1]=show->GetCorrectedEnergy(2,2,1);//e+-/ps/nves
     ecal_enc[2]=show->GetCorrectedEnergy(2,0,0);//e+-/vc/onves
   }

   return true;
}

bool AMSAnalysis::Select_Anti(){
  return true;
}


bool AMSAnalysis::Select_TrkParticle(){

//----Unbiased-Tracker(TOF efficiency study) Max-Momentum
  if(ntrack>=1)itrtracks=0;
  if(ntrack>1){
    float tkmom=0;
    for(int itr=0;itr<ntrack;itr++){
      TrTrackR *trk=pev->pTrTrack(itr);
      if(!trk)continue;
      float ntkmom=fabs(trk->GetRigidity())*trk->GetQ();
      if(ntkmom>tkmom){itrtracks=itr;tkmom=ntkmom;}
    }
  }

//---TrTrackQ+Rigidity(good for Z>=2 particle)
  if(itrtracks<0)return false;
  TrTrackR *trk=pev->pTrTrack(itrtracks);
  if(!trk)return false;
  for(int ibh=0;ibh<nbetah;ibh++){
    BetaHR *betah=pev->pBetaH(ibh);
    if(betah && betah->iTrTrack()==itrtracks){tk_ibetahs=ibh;break;}
  }
  BetaHR *betah=pev->pBetaH(tk_ibetahs);

 for(int iq=0;iq<2;iq++){
    float obeta=1;int fid=0;
    if(iq==1&&betah)obeta=betah->GetBeta(); //Beta Correction
    mean_t tkqa=trk->GetInnerQ_all(obeta,fid);//Inner
    tk_qis[iq]   =tkqa.Mean;
    tk_qirmss[iq]=tkqa.RMS;
    for(int il=0;il<9;il++){//L1L9Q
      tk_qlss[iq][il]=trk->GetLayerJQ(il+1,obeta,fid);
    }
  }

//---Tk quality
  int algo=1;// 1:Choutko 2:Alcaraz
  int patt=0;// 3:Inner 5:L1I 6:IL9 7:L1IL9
  int refit=3;
  float mass=TrFit::Mproton;
  int chrg=1;

//---Init Track ReFit
  float tkfz=trk->GetAdvancedFitCharge();
  tkfz=-1;

  InitTrkReFit(tk_qis[0],tkfz,refit,mass,chrg,algo); 

//---Hit Pattern
  tk_hitbs[0]=trk->GetBitPatternXYJ(); tk_hitbs[1]=trk->GetBitPatternJ();

//---Algo
  for(int isp=0;isp<=4;isp++){
    if     (isp==0)patt=0;//Full Span
    else if(isp==1)patt=3;//Inner
    else if(isp==2)patt=5;//Inner+Layer1
    else if(isp==3)patt=6;//Inner+Layer9
    else if(isp==4)patt=7;//Force Full Span
    refit=30+refit%10;//V6
    int mfit=trk->iTrTrackPar(algo,patt,refit,mass,chrg); 
    tk_rigiditys[isp]=0;
    for(int ichis=0;ichis<3;ichis++){tk_chiss[isp][ichis]=9999;}
    if(mfit>=0){
      tk_rigiditys[isp]=trk->GetRigidity(mfit);
      tk_chiss[isp][0]=trk->GetNormChisqX(mfit);
      tk_chiss[isp][1]=trk->GetNormChisqY(mfit);
      tk_chiss[isp][2]=trk->GetChisq(mfit);
    }
  }
  trk->RecalcHitCoordinates(0,0,0,chrg);//using default id to restore coordinate

//--Track Pos
  for(int ilay=0;ilay<NTKL;ilay++){
    double length,rigidity;
    AMSPlaneM pl=trk->InterpolateLayerJN(length,rigidity,ilay+1);
    int uexl=(ilay==0)?0:1;
    for(int ic=0;ic<3;ic++){
      tk_pos2s[ilay][ic]=(length<=-1000)?(ic==2?0:9999):pl.getP()[ic];
      if(ilay!=0&&ilay!=NTKL-1)continue;
      tk_dir2s[uexl][ic]=(length<=-1000)?(ic==2?-1:0):pl.getD()[ic];
    }
  }

//---TOF
  if(betah){
    tk_betahs=betah->GetBeta();
    TofBetaPar normbetapar=betah->GetNormTofBetaPar();
    tk_tofchists=normbetapar.Chi2T;
    tk_tofchiscs=normbetapar.Chi2C;
    tk_tofsumhus=betah->GetUseHit();
    TofChargeHR *tofhc=betah->pTofChargeH();
    float obetac=tk_betahs;
    tofhc->ReFit(obetac,TofChargeHR::DefaultQOptIonW,tk_rigiditys[0]);
    for(int ilay=0;ilay<NTOFL;ilay++){
      tk_tofqls[ilay]=0;
      if(!tofhc->TestExistHL(ilay))continue;
      tk_tofqls[ilay]=tofhc->GetQL(ilay,22);//New TOFL
    }
  }

  return true;

}


void AMSAnalysis::InitEvent(AMSEventR *ev){
//---amsdn
   amsdn=169; //vers*100+amsd
//---Header
   run=0;event=0;version=0;errorb=0;time[0]=0;time[1]=0;isbadrun=0;
//---Trigger
   nlevel1=-1;
//---RTI
   irti=-1;
//---MC
   mpar=mch=0;
   ntrmccl=ntofmccl=-1;
//--Particle
   nparticle=nbetah=nbeta=ntrack=ntrrawcl=ntrcl=ntrhit=ntrdtrack=nrich=nrichb=nshow=necalhit=ntofclh=ncharge=-1;
   tk_l1q=0;
   trd_pass=-1;
//---
   l1i=l9i=itrtrack=ibetah=ibeta=iparindex=itrdtrack=ishow=irich=irichb=icharge=-1;
   cutoffp[0]=cutoffp[1]=cutoffpi[0]=cutoffpi[1]=0;
//---L1L9Par
   ibetahs=itrdtracks=itrtracks=tk_ibetahs=-1;
   tof_nhits=tof_hsumhus=tk_tofsumhus=-1;
   betahs=tk_betahs=-3;
   tk_qis[0]=tk_qis[1]=0;
   cutoffps=trd_thetas=trd_phis=trd_rqs=0;
//---TOF
   tof_z=0;
   tof_betah=tof_betahmc=tof_beta=-3;
   tof_hsumh=tof_hsumhu=-1;
//---Track
   issmearext=0;
   tk_l1l9par=tk_bgpar=tk_bhgpar=tk_ubpar=tk_ngpar=0;
   tk_iq=tk_1q=tk_9q=-1;
   tk_qr=tk_q[0]=tk_q[1]=tk_rz[0][0]=tk_rz[0][1]=tk_rz[1][0]=tk_rz[1][1]=-1;
   for(int ic=0;ic<2;ic++){
     for(int ixy=0;ixy<3;ixy++){tk_qin[ic][ixy]=tk_qrmn[ic][ixy]=0;}
   }
   tk_rtype=tk_hitb[0]=tk_hitb[1]=tk_nhiti[0]=tk_nhiti[1]=0;
   tk_bcor=tk_rigidity=tk_fzm[0]=tk_fzm[1]=0;
//---Trd
   trd_statk=-1000;
   trd_nhitk=trd_onhitk=trd_oampk=-1;
   trd_lik[0]=trd_lik[1]=trd_lik[2]=trd_like[0]=trd_like[1]=trd_like[2]=-3;
   trd_ampsk=trd_amphitk=trd_pathsk=0;
   trd_qk[0]=0;
   trd_rq=trd_rz=-1;
//--Rich
   rich_itrtrack=-1;
   rich_beta[0]=rich_beta[1]=rich_betap=-3;
   rich_rz=rich_rq=-1;
//--Ecal
   ecal_stat=0;
   ecal_en[0]=ecal_en[1]=ecal_en[2]=ecal_en[3]=ecal_rearl=-3;
   ecal_enc[0]=ecal_enc[1]=ecal_enc[2]=-3;
   ecal_bdt[0]=ecal_bdt[1]=ecal_bdt[2]=-3;
//--Anti
   anti_nhit=-1;

  pev=ev;
}

int AMSAnalysis::Fill(int it){
  if(it!=0){
    TrTrackR *trk=(pev&&itrtrack>=0)?pev->pTrTrack(itrtrack):0;
    tk_rigidity=trk?trk->GetRigidity():0;
  }
  if     (it==0){return tout[0]->Fill();}
  else if(it==1){return tout[1]->Fill();}
  else          {return tout[2]->Fill();}
}

int AMSAnalysis::Save(){
    int fsave=0;
    fout->cd();
    for(int it=0;it<3;it++)tout[it]->Write();
    if(fout->IsZombie())fsave=1;
    fout->Close();
    return fsave;
}

int stop=0;
void shandler(int sig){
  stop=1;
  return;
}

int main(int argc,char ** argv){//in_file ofile num
    char ifile[1000];
    char ofile[1000];
    int rmode=0;//readmode 0 txt 1 root
    int num=-1;
    int nstat=9;
    int nrun=-1;
    cout<<"argc="<<argc<<endl;
    if(argc>=3){
      sprintf(ifile,"%s",argv[1]);
      sprintf(ofile,"%s",argv[2]);
      cout<<"ifile="<<ifile<<endl;
      cout<<"ofile="<<ofile<<endl;
    }
    else {
        printf("Usage: %s  <infile_txt> <ofile.root> <max_events_number>\n",argv[0]);
        exit(-1);
    }

    if (argc>=4){
       rmode=atoi(argv[3]);
       cout<<"rmode="<<rmode<<endl;
    }
    if(argc>=5){
       num=atoi(argv[4]);
       if(num==0)num=-1;
       cout<<"num="<<num<<endl;
    }
//    char *otxt="./bac/";
    char *otxt=0;
    unsigned int srun=0;
    if(argc>=6){
      if(rmode==2){srun=stoul(argv[5]);   cout<<"srun="<<srun<<endl;}
      else        {otxt=argv[5];          cout<<"otxt="<<otxt<<endl;}
    }

    char *ofilecp=0;
    unsigned int sevent=0;
    if(argc>=7){
      if(rmode==2){sevent=stoul(argv[6]); cout<<"sevent="<<sevent<<endl;}
      else        {ofilecp=argv[6];       cout<<"ofilecp="<<ofilecp<<endl;}
    }

    signal(SIGTERM, shandler);
    signal(SIGINT, shandler);

    AMSChain ams;
    int nfile=-1;
    int nread[2]={-1,-1};
    if(rmode==1||rmode==2){
        ams.Add(ifile);
        nread[0]=1;
    }
    else{
       nread[1]=nread[0]=readfile(ifile,ams,nfile,0,otxt,rmode);
       if(nread[0]<0){//read-error
         printf("error files read:%s \n",ifile);
         return 5;
       }
    }//read file

//---
    AMSAnalysis *analysis =new AMSAnalysis();
    analysis->BookFile(ofile);
    analysis->LoadHeader();
///---Baichi
    Long64_t num2= ((num==-1)||(num>ams.GetEntries()))?ams.GetEntries():num;
    cout<<"num2="<<num2<<endl;
    int nevout=300000;
    static unsigned int prrun=0;
    static Long64_t nevtree=0;//Local Tree Events
    static string bf="";//bad Tree
    static string pf="";//Pre Tree
    AMSEventR::fRunList.clear();//Clear AMSEventR:fRunList
/*#ifdef _USEONEEV_
    srun=604311619;
    sevent=751735;
#endif*/
    if(srun!=0&&sevent!=0)num2=1;
	
    //--do TOF MC Tuning //2023.12.30
/*#ifdef _DOBETAMCTUNE_
#ifdef _BETAMCTUNEONLYSMEAR_
	//--only smearing
    TofMCPar::MCtuneDT=-40; //40ps time resolution smearing
    TofMCPar::MCtuneST=0;
#elif defined _BETAMCTUNEONLYSHIFT_
	//--only shift
    TofMCPar::MCtuneDT=0;
    TofMCPar::MCtuneST=-40; //40ps shift
#else
	//--smearing and shift
    TofMCPar::MCtuneDT=-40; //40ps time resolution smearing
    TofMCPar::MCtuneST=-40; //40ps shift
#endif //_BETAMCTUNEONLYSMEAR_ || _BETAMCTUNEONLYSHIFT_
    cout << " Select_Particle doing TOF MC Tuning, TofMCPar::MCtuneDT=" << TofMCPar::MCtuneDT << " TofMCPar::MCtuneST=" << TofMCPar::MCtuneST << endl;
#endif //_DOBETAMCTUNE_*/
	
	//--read my tof time calibration //2024.08.24
	#ifdef MYTOFTIMECALIB
	/*tofcalib = new TofCalib();
	for (int ilthtsht=0; ilthtsht<3; ilthtsht++)
	{
		//string sCalibName = "TofTimeCalib_o16";
		string sCalibName = "TofTimeCalib_he4";
		//tofcalib->ReadCalibParFromRoot(ilthtsht);
		if (ilthtsht==0) tofcalib->ReadCalibParFromRoot(ilthtsht);
		else if (ilthtsht>=1) tofcalib->ReadCalibParFromFile(ilthtsht, sCalibName);
		for (int ifix=0; ifix<2; ifix++) tofcalib->InfoCalibPar(ilthtsht, ifix);
	}*/
	for (int imycalibver=0; imycalibver<NMyTOFTimeCalibVer; imycalibver++)
	{
		cout << " imycalibver=" << imycalibver << endl;
		string sCalibName;
		if (imycalibver==0) sCalibName = "TofTimeCalib_TofCalibN_Z3_n200000";
		else if (imycalibver==1) sCalibName = "TofTimeCalib_TofCalibN_Z8_n200000";
		
		tofcalib[imycalibver] = new TofCalib();
		for (int ilthtsht=0; ilthtsht<3; ilthtsht++)
		{
			tofcalib[imycalibver]->ReadCalibParFromFile(ilthtsht, sCalibName);
			for (int ifix=0; ifix<2; ifix++) tofcalib[imycalibver]->InfoCalibPar(ilthtsht, ifix);
		} //ilthtsht
	} //imycalibver
	#endif //MYTOFTIMECALIB
	
	//--read my calibration for LT, with new parameters using pr, he4, o16, pr+he4//2024.09.05
	#ifdef COMBINE_Z_TOFTIMECALIB
	for (int icalibver=0; icalibver<NTOFTimeCalibVer; icalibver++)
	{
		string sCalibName;
		if (icalibver==0) sCalibName = "TofTimeCalib_pr";
		else if (icalibver==1) sCalibName = "TofTimeCalib_he4";
		else if (icalibver==2) sCalibName = "TofTimeCalib_o16";
		else if (icalibver==3) sCalibName = "TofTimeCalib_pr_he4";
		else if (icalibver==4) sCalibName = "TofTimeCalib_all";
		for (int ilthtsht=0; ilthtsht<1; ilthtsht++) //now only LT
		{
			tofcalibCombZ[icalibver].ReadCalibParFromFile(ilthtsht, sCalibName);
			for (int ifix=0; ifix<2; ifix++) tofcalibCombZ[icalibver].InfoCalibPar(ilthtsht, ifix);
		}
	} //icalibver
	#endif //COMBINE_Z_TOFTIMECALIB
	
	#ifdef TOFBETA_WITHOUT_CAL
	cout << " Tof Time Calibration, Null parameter: " << endl;
	for (int ilthtsht=0; ilthtsht<3; ilthtsht++)
	{
		for (int ifix=0; ifix<2; ifix++) tofcalibNull.InfoCalibPar(ilthtsht, ifix);
	} //ilthtsht
	cout << endl;
	#endif //TOFBETA_WITHOUT_CAL
	
	//--
	ifstream infile;
	infile.open(ifile);
	string sFileName;
	infile >> sFileName;
	infile.close();
	cout << sFileName << endl; //XXX
	sFileName = sFileName.substr(0, sFileName.find_last_of("/"));
	cout << sFileName << endl; //XXX
	sFileName = sFileName.substr(sFileName.find_last_of("/")+1);
	cout << sFileName << endl; //XXX
	sFileName = sFileName.substr(0, sFileName.find_first_of("."));
	cout << sFileName << endl; //XXX
	
	int MCCharge = -1;
	if (sFileName=="pr") MCCharge=1;
	else if (sFileName=="he3" || sFileName=="he4") MCCharge=2;
	else if (sFileName=="li6" || sFileName=="li7") MCCharge=3;
	else if (sFileName=="be7" || sFileName=="be9" || sFileName=="be10") MCCharge=4;
	else if (sFileName=="b10" || sFileName=="b11") MCCharge=5;
	else if (sFileName=="c12" || sFileName=="c13") MCCharge=6;
	else if (sFileName=="n14" || sFileName=="n15") MCCharge=7;
	else if (sFileName=="o16" || sFileName=="o17") MCCharge=8;
	
	cout << " MCCharge=" << MCCharge << endl;
	
	//2024.10.23
	#ifdef NEWTOFTIMECALIB
	cout << " Tof Time Calibration, All sets of parameters: " << endl;
	//--
	for (int iqcalib=0; iqcalib<NZ_NewCalib; iqcalib++)
	//for (int iqcalib=0; iqcalib<9; iqcalib++)
	{
		for (int iminchis=0; iminchis<NMinChis_NewCalib; iminchis++)
		{
			for (int iexcludeedge=0; iexcludeedge<NExcludeEdge_NewCalib; iexcludeedge++)
			{
				//string sCalibName = Form("TofTimeCalib_TofCalibN_Z%d_n200000", ZList_NewCalib[iqcalib]);
				////2024.11.04
				//string sCalibName = Form("TofTimeCalib_TofCalibN_Z%d", ZList_NewCalib[iqcalib]);
				//if (ZList_NewCalib[iqcalib]==-1) sCalibName = Form("TofTimeCalib_TofCalibN_All");
				//else if (ZList_NewCalib[iqcalib]==-2) sCalibName = Form("TofTimeCalib_TofCalibN_Z3_Z8");
				//2024.11.16
				string sCalibName = Form("TofTimeCalib_TofCalibN_Z%d", MCCharge);
				if (iqcalib==0) sCalibName = Form("TofTimeCalib_TofCalibN_All");
				else if (iqcalib==1) sCalibName = Form("TofTimeCalib_TofCalibN_Z3_Z8");
				//2024.12.11
				if (iexcludeedge==1) sCalibName += "_ExcludeEdge";
				
				string sCalibSuffix = "";
				if (iexcludeedge==1) sCalibSuffix += "_ExcludeEdge";
				if (iminchis==1) sCalibSuffix += "_CombMinTimeCoo"; //2024.12.10
				else if (iminchis==2) sCalibSuffix += "_MinChisDer_DerQuadSum";
				else if (iminchis==3) sCalibSuffix += "_MinChisDer";
				
				cout << sCalibName << endl; //XXX
				tofcalibNewCalib[iqcalib][iminchis][iexcludeedge].Charge = ZList_NewCalib[iqcalib];
				for (int ilthtsht=0; ilthtsht<3; ilthtsht++)
				{
					tofcalibNewCalib[iqcalib][iminchis][iexcludeedge].ReadCalibParFromFile(ilthtsht, sCalibName, sCalibSuffix);
					for (int ifix=0; ifix<NFix_NewCalib; ifix++) tofcalibNewCalib[iqcalib][iminchis][iexcludeedge].InfoCalibPar(ilthtsht, ifix);
				} //ilthtsht
			} //iexcludeedge
		} //iminchis
	} //iqcalib
	cout << endl;
	#endif //NEWTOFTIMECALIB
	
    for(Long64_t entry=0; entry<num2; entry++){
        AMSEventR *ev=NULL;
        if(num2==1)ev=(AMSEventR *)ams.GetEventFast(srun,sevent);
#if defined (_USEEVENTORDER_)
        else       ev=(AMSEventR *)ams.GetOrderedEntry(entry);
#else
        else       ev=(AMSEventR *)ams.GetEvent(entry);
#endif
        if(ev==NULL)continue;
        if(entry>=ams.GetEntries())break;
        if(ev->Run()!=prrun){cout<<"ProcessInfo NewRun="<<ev->Run()<<endl;prrun=ev->Run();}
        if(entry%nevout==0)cout<<"ProcessInfo Entry="<<entry<<" Run="<<ev->Run()<<" Event="<<ev->Event()<<endl; 
        int records=ev->RecordRTIRun();//Run to AMSEventR:fRunList
//----Check Run in AMSEventR:fRunList
        TTree *evtree=ev->Tree();
        string nf=evtree->GetCurrentFile()->GetName();
        if (nf==bf)continue;//BadRunList
        if (nf!=pf){//New Tree
          nevtree=0;
          if(records==2){//Exist in RunList
            cout<<"Error ExistFile in RunList="<<nf<<endl;
            pf=bf=nf;//Exist to  BadRunList
            continue;
          }
          else {
            pf=nf;
          }
        }
        else {
          nevtree++;
        }
//---
       if(num2==1)cout<<"run="<<ev->Run()<<" event="<<ev->Event()<<endl;
//       cout<<"run="<<ev->Run()<<" event="<<ev->Event()<<endl;
//---Set Unbiase
       analysis->unbiasetag0=-1;//Pre Unbiase-Tag
       analysis->unbiasetag=0;//ISS
       analysis->istkaddhit=0;//AddLostHit-Tag
//---Set Pre Unbiase
       if(entry%prscale==0)analysis->unbiasetag0=0;
//       cout<<"current event="<<ev->Tree()->GetReadEvent()<<" "<<ev->Tree()->GetEntries()<<" cev"<<nevtree<<endl;
       if((evtree->GetReadEvent()==0)||(evtree->GetReadEvent()==evtree->GetEntries()-1)||(nevtree==0)||(nevtree==evtree->GetEntries()-1)){//tree bein+end |entry begin+end 
         analysis->unbiasetag0=1;
       }
       if(analysis->unbiasetag0>=0){
         analysis->InitEvent(ev);
         int ksel=analysis->Select();
         if(ksel==1){
           int ktree=analysis->Fill(1);
           if(ktree<0){cout<<"Error Write Tree1"<<endl;nstat=6;goto ENDList;}
         }
       }
//----Set ISS|MC Unbiase
        if(ev->nMCEventg()){//MC
           if(fabs(analysis->mch)==1)prunbiasmc=1;//P/e+ prunbiasmc=1
           if(entry%prunbiasmc==0)analysis->unbiasetag=1;
        }
        else {
           if(entry%prunbias==0)analysis->unbiasetag=1;
        } 
        if(analysis->unbiasetag){
          analysis->InitEvent(ev);
          int ksel=analysis->Select();
          if(ksel==1){
            int ktree=analysis->Fill(2);
            if(ktree<0){cout<<"Error Write Tree2"<<endl;nstat=6;goto ENDList;}
          }
        }
        prtime=ev->UTime();
//---
        if(stop==1){cout<<"Force stop"<<endl;nstat=4;goto ENDList;}//force stop
    }

//--Write To List
  {
    map <string,unsigned int >runlistn;
    for(map<unsigned int,AMSSetupR::RunI>::iterator it=AMSEventR::fRunList.begin();it!=AMSEventR::fRunList.end();it++){//Get Run by Run
       AMSSetupR::RunI runi=it->second; //run information
       for(unsigned int k=0;k<runi.fname.size();k++)runlistn[runi.fname[k]]=runi.run;
    }
    nrun=runlistn.size();
    AMSEventR::fRunList.clear();
    int fsave=analysis->Save();
    if(nread[1]>0&&otxt!=0){
      vector<string> sfile;//check again the list before write
      nread[1]=readfile(ifile,sfile,nfile,0,otxt,rmode);
      if(nread[1]==0){nread[0]=nrun=0;}//no-run
    }
    if(fsave){cout<<"Error Saving File="<<ofile<<endl;nstat=6;}
    else if(nrun>0){
//---CopyToEOS
      int cpstat=cpofile(ofile,ofilecp,0);
//---Recording to txt if CopyToEOS success
      if(cpstat==0){WriteList(otxt,runlistn); nstat=0;}//normal
      else         {nstat=3;}//copy-error
    }
    else {nstat=2;}//no-run
    if((nrun!=nread[0])&&(nstat==0||nstat==2))nstat+=10;//nrun!=nread[0]
  }

///--List Status
ENDList:
   const char *sname[12]={
      "Head  ", "Trig  ","RTI|MC","Part    ","Trk0    ","Tof     ","Trk1    ","Trd     ","Rich    ","Ecal    ","Anti    ","TrkPart ",
    };
///---
    cout<<"<<---------Begin Static-------------->>"<<endl;
    cout<<"   Event("<<sstat[0]<<")"<<"   CPU"<<"   TreeSize"<<endl;
    for(int icut=0;icut<12;icut++){
      cout<<sname[icut]<<"   "<<sstat[icut+1]<<"   "<<(double(cput[icut]/CLOCKS_PER_SEC))<<"   "<<sstree[icut]<<endl;
     }
    cout<<"    Process Run="<<nrun<<" Stat="<<nstat<<endl;
    if(nrun!=nread[0])cout<<"    Read Run="<<nread[0]<<endl;
    cout<<"<<---------End Static AMSJOB-------------->>"<<endl;

//-----
   return nstat;//nstat=0 normal, nstat=2 no-run, nstat=3 copy-error, nstat=4 stop-error, nstat=5 read-error,nstat=6 save-error,nstat=9 initialization-error,nstat+10 nrun!=nread[0]
   

}
