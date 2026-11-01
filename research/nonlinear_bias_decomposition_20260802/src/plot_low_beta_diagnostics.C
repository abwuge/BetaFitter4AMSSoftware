#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
constexpr int NB=14;
constexpr double B0=.25, B1=.95, BW=(B1-B0)/NB;
struct R{int fi,z,b2,b3; double b,dt,raw,corr,la2,la3,e2,e3,tu,tl,au,al;};
struct S{int n=0; double m=NAN,q16=NAN,q84=NAN;};
struct P{std::string n,l; int c,m; std::vector<R> v;};
using F=std::function<double(const R&)>;

double q(std::vector<double> v,double p){
  if(v.empty()) return NAN; std::sort(v.begin(),v.end());
  double x=p*(v.size()-1); size_t i=(size_t)x; double f=x-i;
  return v[i]*(1-f)+v[std::min(i+1,v.size()-1)]*f;
}
S stat(const std::vector<double>&v){S s; s.n=v.size(); if(!s.n)return s; s.m=q(v,.5);s.q16=q(v,.16);s.q84=q(v,.84);return s;}
int bin(double b){if(!(b>=B0&&b<B1))return -1;return std::max(0,std::min(NB-1,(int)((b-B0)/BW)));}
std::vector<std::string> split(const std::string&x){std::vector<std::string>o;std::stringstream ss(x);std::string a;while(std::getline(ss,a,','))o.push_back(a);return o;}
std::vector<R> read(const std::string&p){
  std::ifstream in(p); std::string x; std::getline(in,x); std::vector<R>o;
  while(std::getline(in,x)){auto a=split(x);if(a.size()<16)continue;try{
    R r{std::stoi(a[0]),std::stoi(a[1]),std::stoi(a[3]),std::stoi(a[4]),
      std::stod(a[2]),std::stod(a[5]),std::stod(a[6]),std::stod(a[7]),
      std::stod(a[8]),std::stod(a[9]),std::stod(a[10]),std::stod(a[11]),
      std::stod(a[12]),std::stod(a[13]),std::stod(a[14]),std::stod(a[15])};
    if(std::isfinite(r.b)&&bin(r.b)>=0)o.push_back(r);
  }catch(...){}}
  return o;
}
TGraphAsymmErrors* graph(const P&p,F f,const std::string&n,int color=-1,int marker=-1){
  std::array<std::vector<double>,NB>a;for(auto&r:p.v){int k=bin(r.b);double y=f(r);if(k>=0&&std::isfinite(y))a[k].push_back(y);}
  auto*g=new TGraphAsymmErrors;g->SetName(n.c_str());int ip=0;
  for(int k=0;k<NB;k++){auto s=stat(a[k]);if(s.n<10)continue;double x=B0+(k+.5)*BW;g->SetPoint(ip,x,s.m);g->SetPointError(ip,BW/2,BW/2,s.m-s.q16,s.q84-s.m);ip++;}
  int c=color<0?p.c:color, m=marker<0?p.m:marker;g->SetLineColor(c);g->SetMarkerColor(c);g->SetMarkerStyle(m);g->SetMarkerSize(.9);g->SetLineWidth(2);return g;
}
TH2D* density(const P&p,F fx,F fy,const std::string&n,const std::string&t,int nx,double x0,double x1,int ny,double y0,double y1){
  auto*h=new TH2D(n.c_str(),t.c_str(),nx,x0,x1,ny,y0,y1);h->SetDirectory(nullptr);for(auto&r:p.v){double x=fx(r),y=fy(r);if(std::isfinite(x)&&std::isfinite(y))h->Fill(x,y);}return h;
}
void zero(){TLine l(B0,0,B1,0);l.SetLineStyle(2);l.SetLineColor(kGray+2);l.SetLineWidth(2);l.DrawClone();}
void text(const std::string&s){TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.036);t.DrawLatex(.14,.92,s.c_str());}
void axes(TMultiGraph&m,const char*x,const char*y,double lo,double hi){m.GetXaxis()->SetLimits(B0,B1);m.GetXaxis()->SetTitle(x);m.GetYaxis()->SetTitle(y);m.GetYaxis()->SetTitleOffset(1.28);m.SetMinimum(lo);m.SetMaximum(hi);}
void save(TCanvas*c,const std::string&d,const std::string&n){c->SaveAs((d+"/"+n+".png").c_str());c->SaveAs((d+"/"+n+".pdf").c_str());}

void timing(const P&p,const std::string&d,TDirectory*od){
  auto*gr=graph(p,[](auto&r){return 1000*r.raw;},p.n+"_raw_dt23_ps",kRed+1,20);
  auto*gc=graph(p,[](auto&r){return 1000*r.corr;},p.n+"_applied_correction23_ps",kBlue+1,21);
  auto*gt=graph(p,[](auto&r){return 1000*r.dt;},p.n+"_corrected_dt23_ps",kBlack,20);
  auto*h=density(p,[](auto&r){return r.b;},[](auto&r){return 1000*r.dt;},p.n+"_corrected_dt23_density",";#beta_{center}^{MC};#Delta t_{23}^{corr} [ps]",70,B0,B1,180,-450,450);
  auto*c=new TCanvas((p.n+"_timing_canvas").c_str(),"",1350,950);c->Divide(2,2,.008,.008);
  c->cd(1);gPad->SetGrid();auto*m1=new TMultiGraph;m1->Add(gr,"PZ");m1->Draw("A");axes(*m1,"#beta_{center}^{MC}","raw S2-S3 - checkpoint truth [ps]",-450,450);zero();text(p.l+": before slew/delay correction");
  c->cd(2);gPad->SetGrid();auto*m2=new TMultiGraph;m2->Add(gc,"PZ");m2->Draw("A");axes(*m2,"#beta_{center}^{MC}","applied S2-S3 correction [ps]",-450,450);zero();text("(tof_{tl}-raw)_{S3}-(tof_{tl}-raw)_{S2}");
  c->cd(3);gPad->SetGrid();auto*m3=new TMultiGraph;m3->Add(gt,"PZ");m3->Draw("A");axes(*m3,"#beta_{center}^{MC}","corrected S2-S3 - checkpoint truth [ps]",-140,140);zero();text("standard TofClusterHR::Time after correction");
  c->cd(4);gPad->SetRightMargin(.13);gPad->SetLogz();h->Draw("COLZ");gt->Draw("PZ SAME");zero();text("event density; black = median and 16-84% interval");
  save(c,d,"tof_timewalk_vs_beta_"+p.n);od->cd();gr->Write();gc->Write();gt->Write();h->Write();c->Write();
}

void energy(const P&p,const std::string&d,TDirectory*od){
  F truth=[](auto&r){return r.tu+r.tl;}, active=[](auto&r){return r.au+r.al;};
  F miss=[&](auto&r){return truth(r)-active(r);};
  F rat=[&](auto&r){double a=active(r);return a>0?truth(r)/a:NAN;};
  F ru=[](auto&r){return r.au>0?r.tu/r.au:NAN;}, rl=[](auto&r){return r.al>0?r.tl/r.al:NAN;};
  auto*gt=graph(p,truth,p.n+"_truth_total_dE",kRed+1,20);auto*ga=graph(p,active,p.n+"_active_total_dE",kBlue+1,21);
  auto*gm=graph(p,miss,p.n+"_missing_total_dE",kMagenta+1,22);auto*gr=graph(p,rat,p.n+"_truth_over_active",kBlack,20);
  auto*gu=graph(p,ru,p.n+"_truth_over_active_upper",kRed+1,24);auto*gl=graph(p,rl,p.n+"_truth_over_active_lower",kBlue+1,25);
  double xm=0,ym=0;for(auto&r:p.v){xm=std::max(xm,active(r));ym=std::max(ym,truth(r));}xm=std::max(20.,std::min(600.,1.02*xm));ym=std::max(100.,std::min(5000.,1.02*ym));
  auto*h=density(p,active,truth,p.n+"_truth_vs_active",";tracker active #Sigma tk_edep [MeV];checkpoint total #DeltaE [MeV]",120,0,xm,140,0,ym);
  auto*c=new TCanvas((p.n+"_energy_canvas").c_str(),"",1350,950);c->Divide(2,2,.008,.008);
  c->cd(1);gPad->SetGrid();auto*m1=new TMultiGraph;m1->Add(gt,"PZ");m1->Add(ga,"PZ");m1->Draw("A");m1->GetXaxis()->SetLimits(B0,B1);m1->GetXaxis()->SetTitle("#beta_{center}^{MC}");m1->GetYaxis()->SetTitle("S2-center-S3 energy loss [MeV]");m1->SetMinimum(0);m1->GetYaxis()->SetTitleOffset(1.28);auto*l1=new TLegend(.52,.70,.88,.87);l1->SetBorderSize(0);l1->SetFillStyle(0);l1->AddEntry(gt,"checkpoint total #DeltaE","lp");l1->AddEntry(ga,"tracker active #Sigma tk_edep","lp");l1->Draw();text(p.l+": total central energy loss");
  c->cd(2);gPad->SetGrid();auto*m2=new TMultiGraph;m2->Add(gm,"PZ");m2->Draw("A");m2->GetXaxis()->SetLimits(B0,B1);m2->GetXaxis()->SetTitle("#beta_{center}^{MC}");m2->GetYaxis()->SetTitle("checkpoint #DeltaE - active #Sigma tk_edep [MeV]");m2->SetMinimum(0);m2->GetYaxis()->SetTitleOffset(1.28);text("energy change absent from active-silicon deposits");
  c->cd(3);gPad->SetGrid();auto*m3=new TMultiGraph;m3->Add(gr,"PZ");m3->Add(gu,"PZ");m3->Add(gl,"PZ");m3->Draw("A");m3->GetXaxis()->SetLimits(B0,B1);m3->GetXaxis()->SetTitle("#beta_{center}^{MC}");m3->GetYaxis()->SetTitle("checkpoint #DeltaE / tracker active #Sigma tk_edep");m3->SetMinimum(0);m3->GetYaxis()->SetTitleOffset(1.28);auto*l3=new TLegend(.53,.67,.88,.87);l3->SetBorderSize(0);l3->SetFillStyle(0);l3->AddEntry(gr,"S2-center-S3 total","lp");l3->AddEntry(gu,"S2 to center","lp");l3->AddEntry(gl,"center to S3","lp");l3->Draw();
  c->cd(4);gPad->SetRightMargin(.13);gPad->SetLogz();h->Draw("COLZ");TLine id(0,0,std::min(xm,ym),std::min(xm,ym));id.SetLineStyle(2);id.SetLineColor(kGray+2);id.DrawClone();text("event comparison; dashed = total #DeltaE equals active tk_edep");
  save(c,d,"central_energy_loss_vs_beta_"+p.n);od->cd();gt->Write();ga->Write();gm->Write();gr->Write();gu->Write();gl->Write();h->Write();c->Write();
}

void combined(const std::vector<P>&ps,const std::string&d,TDirectory*od){
  auto*c=new TCanvas("combined_low_beta_summary","",1350,620);c->Divide(2,1,.008,.008);
  c->cd(1);gPad->SetGrid();auto*mt=new TMultiGraph;auto*lt=new TLegend(.58,.68,.88,.87);lt->SetBorderSize(0);lt->SetFillStyle(0);std::vector<TGraphAsymmErrors*> gs;
  for(auto&p:ps){auto*g=graph(p,[](auto&r){return 1000*r.dt;},p.n+"_combined_dt");mt->Add(g,"PZ");lt->AddEntry(g,p.l.c_str(),"lp");gs.push_back(g); }mt->Draw("A");axes(*mt,"#beta_{center}^{MC}","corrected S2-S3 - checkpoint truth [ps]",-140,140);zero();lt->Draw();text("residual after standard slew/delay correction");
  c->cd(2);gPad->SetGrid();auto*mr=new TMultiGraph;auto*lr=new TLegend(.58,.68,.88,.87);lr->SetBorderSize(0);lr->SetFillStyle(0);
  for(auto&p:ps){auto*g=graph(p,[](auto&r){double a=r.au+r.al;return a>0?(r.tu+r.tl)/a:NAN;},p.n+"_combined_ratio");mr->Add(g,"PZ");lr->AddEntry(g,p.l.c_str(),"lp");gs.push_back(g); }mr->Draw("A");mr->GetXaxis()->SetLimits(B0,B1);mr->GetXaxis()->SetTitle("#beta_{center}^{MC}");mr->GetYaxis()->SetTitle("checkpoint total #DeltaE / active #Sigma tk_edep");mr->SetMinimum(0);mr->SetMaximum(18);mr->GetYaxis()->SetTitleOffset(1.28);lr->Draw();text("central S2-center-S3 material effect");
  save(c,d,"combined_low_beta_diagnostics");od->cd();for(auto*g:gs)g->Write();c->Write();
}
}

void plot_low_beta_diagnostics(const char*base="/nas01/home/hxwu/AMSSoft/betaFitter/research/nonlinear_bias_decomposition_20260802"){
  gStyle->SetOptStat(0);gStyle->SetPadTickX(1);gStyle->SetPadTickY(1);gStyle->SetLegendFont(42);
  std::string b(base),e=b+"/results/plot_events",o=b+"/plots/root_diagnostics";gSystem->mkdir(o.c_str(),true);
  std::vector<P>ps{{"He4","^{4}He",kGreen+2,20,read(e+"/He4_events.csv")},{"C12","^{12}C",kBlue+1,21,read(e+"/C12_events.csv")},{"O16","^{16}O",kRed+1,22,read(e+"/O16_events.csv")}};
  TFile f((o+"/low_beta_diagnostics.root").c_str(),"RECREATE");for(auto&p:ps){std::cout<<p.n<<" events="<<p.v.size()<<"\n";auto*d=f.mkdir(p.n.c_str());timing(p,o,d);energy(p,o,d);}auto*d=f.mkdir("combined");combined(ps,o,d);f.Write();f.Close();std::cout<<"wrote "<<o<<"\n";
}
