#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <cmath>
#include "root.h"
#include "TrFit.h"

static const double DLOGP[6] = {
 -0.00036466808523982763, -0.0008004530100151896,
 -0.0003674109175335616, -0.0008043602865654975,
 -0.0003658904315670952,  -0.0008044497226364911
};

static void make_profile(double x0,double y0,double z0,double tx0,double ty0,double R0,
                         const float hit[9][3], double yout[7], double tout[7]) {
  // Inner layers L2..L8 correspond hit[1]..hit[7], anchor at L5=hit[4].
  yout[3]=y0; tout[3]=ty0;
  // Downward L5->L8.
  double x=x0,y=y0,z=z0,tx=tx0,ty=ty0,R=R0;
  for(int i=3;i<6;i++){
    double Rend=R*std::exp(DLOGP[i]);
    double Rmid=std::sqrt(R*Rend);
    TrProp p(AMSPoint(x,y,z),AMSDir(tx,ty,1),Rmid);
    p.Propagate(hit[i+2][2]);
    AMSPoint q=p.GetP0();
    x=q[0]; y=q[1]; z=q[2]; tx=p.GetDxDz(); ty=p.GetDyDz(); R=Rend;
    yout[i+1]=y; tout[i+1]=ty;
  }
  // Upward L5->L2, reverse the expected energy loss.
  x=x0;y=y0;z=z0;tx=tx0;ty=ty0;R=R0;
  for(int i=2;i>=0;i--){
    double Rend=R*std::exp(-DLOGP[i]);
    double Rmid=std::sqrt(R*Rend);
    TrProp p(AMSPoint(x,y,z),AMSDir(tx,ty,1),Rmid);
    p.Propagate(hit[i+1][2]);
    AMSPoint q=p.GetP0();
    x=q[0]; y=q[1]; z=q[2]; tx=p.GetDxDz(); ty=p.GetDyDz(); R=Rend;
    yout[i]=y; tout[i]=ty;
  }
}

int main(int argc,char**argv){
  if(argc<3){std::cerr<<"input.root output.root\n";return 1;}
  TFile in(argv[1]); TTree *t=(TTree*)in.Get("events"); if(!t)return 2;
  float hit[9][3],pos[9][3],dir[9][3]; double R; int event_index;
  t->SetBranchAddress("tk_hitc",hit); t->SetBranchAddress("tk_pos",pos);
  t->SetBranchAddress("tk_dir",dir); t->SetBranchAddress("ref_R",&R);
  t->SetBranchAddress("ref_event_index",&event_index);
  TFile out(argv[2],"RECREATE"); TTree o("reverse","R-conditioned no-scatter reference with mean energy loss");
  double yref[7],tref[7],jac[7][3];
  o.Branch("event_index",&event_index,"event_index/I");
  o.Branch("yref",yref,"yref[7]/D"); o.Branch("tref",tref,"tref[7]/D");
  o.Branch("jac",jac,"jac[7][3]/D");
  long n=t->GetEntries();
  for(long e=0;e<n;e++){
    t->GetEntry(e);
    double x=pos[4][0],y=pos[4][1],z=pos[4][2];
    double tx=dir[4][0]/dir[4][2],ty=dir[4][1]/dir[4][2];
    make_profile(x,y,z,tx,ty,R,hit,yref,tref);
    for(int k=0;k<3;k++){
      double eps=(k==0?0.001:(k==1?1e-5:0.005));
      double ym[7],tm[7],yp[7],tp[7];
      double xm=x,ym0=y,txm=tx,tym=ty,Rm=R;
      double xp=x,yp0=y,txp=tx,typ=ty,Rp=R;
      if(k==0){ym0-=eps;yp0+=eps;}
      if(k==1){tym-=eps;typ+=eps;}
      if(k==2){Rm*=std::exp(-eps);Rp*=std::exp(eps);}
      make_profile(xm,ym0,z,txm,tym,Rm,hit,ym,tm);
      make_profile(xp,yp0,z,txp,typ,Rp,hit,yp,tp);
      for(int l=0;l<7;l++) jac[l][k]=(yp[l]-ym[l])/(2*eps);
    }
    o.Fill();
    if(e%10000==0) std::cout<<"REV "<<e<<"/"<<n<<" R="<<R<<" y2="<<yref[0]<<" y8="<<yref[6]<<std::endl;
  }
  o.Write(); out.Close(); return 0;
}
