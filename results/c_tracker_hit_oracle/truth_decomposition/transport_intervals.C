#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <cmath>
#include "root.h"
#include "TrFit.h"
int main(){
 TFile f("aligned.root");TTree*t=(TTree*)f.Get("events");float p[21],x[21][3],d[21][3];int ei;
 t->SetBranchAddress("mevmom1",p);t->SetBranchAddress("mevcoo1",x);t->SetBranchAddress("mevdir1",d);t->SetBranchAddress("ref_event_index",&ei);
 TFile out("intervals.root","RECREATE");TTree o("intervals","truth transport increments relative to native field propagation");
 double kick[8][2],bend[8][2],offset[8][2],dlogp[8];
 o.Branch("event_index",&ei,"event_index/I");o.Branch("kick",kick,"kick[8][2]/D");o.Branch("bend",bend,"bend[8][2]/D");o.Branch("offset",offset,"offset[8][2]/D");o.Branch("dlogp",dlogp,"dlogp[8]/D");
 int idx[9]={0,7,8,9,10,11,12,13,18};
 for(long e=0;e<t->GetEntries();e++){
  t->GetEntry(e);
  for(int j=0;j<8;j++){
   int a=idx[j],b=idx[j+1];
   if(p[a]<=0||p[b]<=0){for(int k=0;k<2;k++)kick[j][k]=bend[j][k]=offset[j][k]=1e10;dlogp[j]=1e10;continue;}
   double tx=d[a][0]/d[a][2],ty=d[a][1]/d[a][2];
   TrProp prop(AMSPoint(x[a][0],x[a][1],x[a][2]),AMSDir(tx,ty,1),(p[a]+p[b])/12.);
   prop.Propagate(x[b][2]);AMSPoint pos=prop.GetP0();double slopes[2]={prop.GetDxDz(),prop.GetDyDz()};
   for(int k=0;k<2;k++){kick[j][k]=atan(d[b][k]/d[b][2])-atan(slopes[k]);bend[j][k]=atan(slopes[k])-atan(d[a][k]/d[a][2]);offset[j][k]=x[b][k]-pos[k];}
   dlogp[j]=log(p[b]/p[a]);
  }
  o.Fill();
 }
 o.Write();out.Close();return 0;
}
