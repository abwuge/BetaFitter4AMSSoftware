#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <iostream>
#include <cmath>
#include "root.h"
#include "TrFit.h"
#include "tkdcards.h"
#include "TkDBc.h"
int main(int argc,char**argv){
 if(argc<3){std::cerr<<"input.root output.root [N]";return 1;}
 TFile input(argv[1]);TTree*t=(TTree*)input.Get("events"); if(!t)return 2;
 float hit[9][3],truth[21][3],dir[21][3],mom[21],pos[9][3],rdir[9][3];double R;int event_index,bits[2];
 t->SetBranchAddress("tk_hitc",hit);t->SetBranchAddress("mevcoo1",truth);t->SetBranchAddress("mevdir1",dir);t->SetBranchAddress("mevmom1",mom);t->SetBranchAddress("tk_pos",pos);t->SetBranchAddress("tk_dir",rdir);t->SetBranchAddress("ref_R",&R);t->SetBranchAddress("ref_event_index",&event_index);t->SetBranchAddress("tk_hitb",bits);
 TFile out(argv[2],"RECREATE");TTree o("probe","native pure-field trajectories, not fitted estimates");
 double ideal[9][3],jac[9][2][5],idealrec[9][3],jacrec[9][2][5],field[9][3];
 double fits[4];
 o.Branch("event_index",&event_index,"event_index/I");o.Branch("ideal",ideal,"ideal[9][3]/D");o.Branch("jac",jac,"jac[9][2][5]/D");o.Branch("idealrec",idealrec,"idealrec[9][3]/D");o.Branch("jacrec",jacrec,"jacrec[9][2][5]/D");o.Branch("field",field,"field[9][3]/D");o.Branch("fits",fits,"fits[4]/D");
 long n=t->GetEntries();if(argc>3)n=std::min(n,atol(argv[3]));
 for(long e=0;e<n;e++){
  t->GetEntry(e);
  for(int m=0;m<2;m++){
   double x=m?pos[4][0]:truth[10][0],y=m?pos[4][1]:truth[10][1],z=m?pos[4][2]:truth[10][2];
   double tx=m?rdir[4][0]/rdir[4][2]:dir[10][0]/dir[10][2],ty=m?rdir[4][1]/rdir[4][2]:dir[10][1]/dir[10][2];
   double r=m?R:mom[10]/6.;
   for(int k=-1;k<5;k++){
    for(int s=(k<0?0:-1);s<=1;s+=2){
     double eps=(k<2?0.001:(k<4?1e-5:0.005));
     double xx=x,yy=y,txx=tx,tyy=ty,rr=r;
     if(k==0)xx+=s*eps;if(k==1)yy+=s*eps;if(k==2)txx+=s*eps;if(k==3)tyy+=s*eps;if(k==4)rr*=exp(s*eps);
     TrProp p(AMSPoint(xx,yy,z),AMSDir(txx,tyy,1),rr);
     for(int l=0;l<9;l++){
      TrProp q=p;double zz=hit[l][2];q.Propagate(zz);AMSPoint v=q.GetP0();
      for(int c=0;c<3;c++){
       if(k<0){if(m)idealrec[l][c]=v[c];else ideal[l][c]=v[c];}
       else if(c<2){double& j=m?jacrec[l][c][k]:jac[l][c][k];if(s<0)j=-v[c]/(2*eps);else j+=v[c]/(2*eps);}
      }
     }
     if(k<0)break;
    }
   }
  }
  for(int l=0;l<9;l++){TrProp::GuFld(hit[l][0],hit[l][1],hit[l][2],field[l]);}
  // Controlled Choutko hit replacement at identical geometry and fixed errors.
  int mi[9]={0,7,8,9,10,11,12,13,18};
  for(int k=0;k<4;k++){
   TrFit f;
   for(int l=1;l<8;l++){
    int j=mi[l];if(!(bits[1]&(1<<l))||truth[j][2]<-999)continue;
    double dz=hit[l][2]-truth[j][2];
    double x=hit[l][0],y=hit[l][1];
    if(k>=2){x=truth[j][0]+dir[j][0]/dir[j][2]*dz;y=truth[j][1]+dir[j][1]/dir[j][2]*dz;}
    double ex=(bits[0]&(1<<l))?0.002:0;
    f.Add(x,y,hit[l][2],ex,0.0005,0.00001);
   }
   double stat=f.DoFit(TrFit::CHOUTKO,k%2,0,6,11.1779);
   fits[k]=(stat<0)?0:f.GetRigidity();
  }
  o.Fill();if(e%2000==0)std::cout<<"PROBE "<<e<<"/"<<n<<" R="<<R<<" B="<<field[4][0]<<" dr="<<ideal[4][1]-truth[10][1]<<" fits="<<fits[0]<<","<<fits[1]<<","<<fits[2]<<","<<fits[3]<<std::endl;
 }
 o.Write();out.Close();return 0;
}
