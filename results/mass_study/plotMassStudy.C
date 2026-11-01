#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TColor.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

namespace MassPlots {
struct Row {int z,rich,tag,bin;double R,u[4],ur,m2[4],score[3],forced,ri,rbefore,rafter;};
double quantile(std::vector<double>v,double p){if(v.empty())return NAN;std::sort(v.begin(),v.end());double j=p*(v.size()-1);size_t i=size_t(j);return v[i]+(j-i)*(v[std::min(i+1,v.size()-1)]-v[i]);}
double w68(const std::vector<double>&v){return .5*(quantile(v,.84)-quantile(v,.16));}
const char* labels[4]={"Raw TOF","Mass offset","Mass calibration","RICH reference"};
int colors[4]={kP6Blue,kP6Blue+1,kP6Blue+2,kP6Blue+3};
void style(TH1*h,int k){h->SetLineColor(colors[k]);h->SetLineWidth(3);h->SetStats(false);h->GetXaxis()->SetTitleSize(.047);h->GetYaxis()->SetTitleSize(.047);h->GetXaxis()->SetLabelSize(.037);h->GetYaxis()->SetLabelSize(.037);h->GetYaxis()->SetTitleOffset(1.0);}
void legend(TLegend&l){l.SetFillStyle(0);l.SetBorderSize(0);l.SetTextSize(.036);}
void text(double x,double y,const std::string&s,int color=kBlack,double size=.033){TLatex a;a.SetNDC();a.SetTextColor(color);a.SetTextSize(size);a.DrawLatex(x,y,s.c_str());}
std::vector<std::string> pages;
void save(TCanvas&c,const std::string&dir,const std::string&name){c.Modified();c.Update();c.SaveAs((dir+"/plots/"+name+".png").c_str());std::string p=dir+"/plots/"+name+".pdf";c.SaveAs(p.c_str());pages.push_back(p);}
}

void plotMassStudy(const char*dest="results/mass_study"){
  using namespace MassPlots;std::string dir(dest);gSystem->mkdir((dir+"/plots").c_str(),true);gROOT->SetStyle("Pub");gROOT->ForceStyle();gStyle->SetOptStat(0);gStyle->SetOptTitle(1);gStyle->SetTitleX(.5);gStyle->SetTitleAlign(23);gStyle->SetTitleY(.985);gStyle->SetTitleH(.055);gStyle->SetTitleBorderSize(0);gStyle->SetTitleTextColor(kBlack);gStyle->SetTitleFont(42,"");gStyle->SetOptFit(0);gStyle->SetTitleFontSize(.045);gStyle->SetPadLeftMargin(.12);gStyle->SetPadRightMargin(.04);gStyle->SetPadBottomMargin(.13);gStyle->SetPadTopMargin(.10);
  TFile in((dir+"/predictions.root").c_str());auto t=dynamic_cast<TTree*>(in.Get("predictions"));if(!t)throw std::runtime_error("Missing predictions");Row r;
  t->SetBranchAddress("Z",&r.z);t->SetBranchAddress("R",&r.R);t->SetBranchAddress("richOK",&r.rich);t->SetBranchAddress("tag",&r.tag);t->SetBranchAddress("bin",&r.bin);t->SetBranchAddress("rawU",&r.u[0]);t->SetBranchAddress("offsetU",&r.u[1]);t->SetBranchAddress("calibratedU",&r.u[2]);t->SetBranchAddress("richU0",&r.ur);t->SetBranchAddress("mass2Raw",&r.m2[0]);t->SetBranchAddress("mass2Offset",&r.m2[1]);t->SetBranchAddress("mass2Calibrated",&r.m2[2]);t->SetBranchAddress("scoreRaw",&r.score[0]);t->SetBranchAddress("scoreCalibrated",&r.score[1]);t->SetBranchAddress("nearestMass",&r.forced);t->SetBranchAddress("referenceU",&r.u[3]);t->SetBranchAddress("mass2Reference",&r.m2[3]);t->SetBranchAddress("scoreReference",&r.score[2]);t->SetBranchAddress("Rinjected",&r.ri);t->SetBranchAddress("mass2RigidityInjectionBefore",&r.rbefore);t->SetBranchAddress("mass2RigidityInjectionAfter",&r.rafter);
  std::vector<Row>data[3];for(Long64_t i=0;i<t->GetEntries();i++){t->GetEntry(i);data[r.z].push_back(r);}
  TFile histout((dir+"/plots.root").c_str(),"RECREATE");std::ofstream fits(dir+"/gaussian_core.tsv");fits<<std::setprecision(10)<<"Z\tmethod\tN\tmean\tmean_error\tsigma\tsigma_error\tW68\tfit_status\n";
  TCanvas c("mass_study_canvas","Mass study",1280,720);std::string pdf=dir+"/mass-study.pdf";pages.clear();
  for(int z=1;z<=2;z++){
    std::string zs="Z"+std::to_string(z);auto&d=data[z];double lo=z==1?-2.:-5.,hi=z==1?7.:28.;
    std::unique_ptr<TH1D> h[4];for(int k=0;k<4;k++){h[k].reset(new TH1D((zs+"_mass2_"+std::to_string(k)).c_str(),("Z = "+std::to_string(z)+", 2.5 < R_{0} < 6 GV;m_{eff}^{2} [GeV^{2}/c^{4}];Counts").c_str(),180,lo,hi));h[k]->SetDirectory(nullptr);style(h[k].get(),k);for(auto&v:d)h[k]->Fill(v.m2[k]);}
    c.Clear();c.SetLogy(true);double mx=0;for(int k=0;k<4;k++)mx=std::max(mx,h[k]->GetMaximum());h[0]->SetMaximum(mx*6);h[0]->SetMinimum(.65);for(int k=0;k<4;k++)h[k]->Draw(k?"hist same":"hist");TLegend lg(.63,.63,.93,.87);legend(lg);for(int k=0;k<4;k++){auto e=lg.AddEntry(h[k].get(),labels[k],"l");e->SetTextColor(colors[k]);}lg.Draw();text(.15,.84,"N = "+std::to_string(d.size()));save(c,dir,zs+"_mass2");histout.cd();for(int k=0;k<4;k++)h[k]->Write();

    std::vector<double>res[4];for(auto&v:d)if(v.rich)for(int k=0;k<4;k++)res[k].push_back(v.u[k]-v.ur);
    std::unique_ptr<TH1D>hr[4];std::unique_ptr<TF1>f[4];for(int k=0;k<4;k++){hr[k].reset(new TH1D((zs+"_rich_residual_"+std::to_string(k)).c_str(),("Z = "+std::to_string(z)+", NaF reference;u_{TOF,eff} - u_{RICH#rightarrow0};Counts").c_str(),180,-.16,.20));hr[k]->SetDirectory(nullptr);style(hr[k].get(),k);for(double v:res[k])hr[k]->Fill(v);double median=quantile(res[k],.5),si=.741301109*(quantile(res[k],.75)-quantile(res[k],.25));f[k].reset(new TF1((zs+"_g_"+std::to_string(k)).c_str(),"gaus",median-3*si,median+3*si));f[k]->SetParameters(hr[k]->GetMaximum(),median,si);int status=hr[k]->Fit(f[k].get(),"QRS0");f[k]->SetLineColor(colors[k]);f[k]->SetLineStyle(2);f[k]->SetLineWidth(3);fits<<z<<'\t'<<labels[k]<<'\t'<<res[k].size()<<'\t'<<f[k]->GetParameter(1)<<'\t'<<f[k]->GetParError(1)<<'\t'<<std::abs(f[k]->GetParameter(2))<<'\t'<<f[k]->GetParError(2)<<'\t'<<w68(res[k])<<'\t'<<status<<'\n';}
    c.Clear();c.SetLogy(false);mx=0;for(int k=0;k<4;k++)mx=std::max(mx,hr[k]->GetMaximum());hr[0]->SetMaximum(mx*1.20);for(int k=0;k<4;k++)hr[k]->Draw(k?"hist same":"hist");for(int k=0;k<4;k++)f[k]->Draw("same");text(.15,.84,"N = "+std::to_string(res[0].size()));text(.15,.79,"0.78 < #beta_{RICH} < 0.97",kBlack,.030);for(int k=0;k<4;k++){double y=.83-.16*k;text(.66,y,labels[k],colors[k]);text(.66,y-.046,Form("#mu = %.4f #pm %.4f",f[k]->GetParameter(1),f[k]->GetParError(1)),colors[k],.030);text(.66,y-.090,Form("#sigma = %.4f #pm %.4f",std::abs(f[k]->GetParameter(2)),f[k]->GetParError(2)),colors[k],.030);}save(c,dir,zs+"_rich_residual");histout.cd();for(int k=0;k<4;k++)hr[k]->Write();

    std::array<TGraphErrors,4>graphs;double centers[4]={2.75,3.5,4.5,5.5};for(int k=0;k<4;k++){graphs[k].SetLineColor(colors[k]);graphs[k].SetMarkerColor(colors[k]);graphs[k].SetLineWidth(3);graphs[k].SetMarkerStyle(20+k);graphs[k].SetMarkerSize(1.25);int p=0;for(int b=0;b<4;b++){std::vector<double>values;for(auto&v:d)if(v.rich&&v.bin==b)values.push_back(v.u[k]-v.ur);if(values.size()<40)continue;graphs[k].SetPoint(p,centers[b],w68(values));graphs[k].SetPointError(p,0,0);p++;}}
    c.Clear();auto frame=c.DrawFrame(2.4,0,6.1,z==1?.10:.075);frame->SetTitle(("Z = "+std::to_string(z)+", NaF reference;R_{0} [GV];W_{68}(u_{TOF,eff} - u_{RICH#rightarrow0})").c_str());style(frame,0);for(int k=0;k<4;k++)graphs[k].Draw("LP same");TLegend gl(.61,.61,.92,.86);legend(gl);for(int k=0;k<4;k++){auto e=gl.AddEntry(&graphs[k],labels[k],"lp");e->SetTextColor(colors[k]);}gl.Draw();save(c,dir,zs+"_width_by_rigidity");

    TGraph roc[3];for(int k=0;k<3;k++){std::vector<std::pair<double,int>>s;int n0=0,n1=0;for(auto&v:d)if(v.tag>=0){s.push_back({v.score[k],v.tag});n0+=v.tag==0;n1+=v.tag==1;}std::sort(s.begin(),s.end(),[](const std::pair<double,int>&a,const std::pair<double,int>&b){return a.first>b.first;});int a=0,b=0,p=0;roc[k].SetPoint(p++,0,0);for(size_t i=0;i<s.size();i++){a+=s[i].second==0;b+=s[i].second==1;if(i+1==s.size()||s[i+1].first!=s[i].first)roc[k].SetPoint(p++,double(a)/std::max(1,n0),double(b)/std::max(1,n1));}roc[k].SetLineColor(colors[k==0?0:k+1]);roc[k].SetLineWidth(3);}
    c.Clear();frame=c.DrawFrame(0,0,1,1.03);frame->SetTitle(("Z = "+std::to_string(z)+", RICH mass tags;Light-tag false positive rate;Heavy-tag efficiency").c_str());style(frame,0);roc[0].Draw("L same");roc[1].Draw("L same");roc[2].Draw("L same");TLegend rl(.56,.19,.93,.39);legend(rl);rl.AddEntry(&roc[0],labels[0],"l")->SetTextColor(colors[0]);rl.AddEntry(&roc[1],labels[2],"l")->SetTextColor(colors[2]);rl.AddEntry(&roc[2],labels[3],"l")->SetTextColor(colors[3]);rl.Draw();text(.53,.45,"Tags share Tracker information",kBlack,.030);save(c,dir,zs+"_rich_tag_roc");

    std::unique_ptr<TH1D>ht[4];const char*name4[4]={"Light tag: raw","Heavy tag: raw","Light tag: calibrated","Heavy tag: calibrated"};for(int k=0;k<4;k++){ht[k].reset(new TH1D((zs+"_tag_"+std::to_string(k)).c_str(),("Z = "+std::to_string(z)+", RICH mass tags;m_{eff}^{2} [GeV^{2}/c^{4}];Counts").c_str(),120,lo,hi));ht[k]->SetDirectory(nullptr);style(ht[k].get(),0);ht[k]->SetLineColor(kP6Blue+k);for(auto&v:d)if(v.tag==k%2)ht[k]->Fill(v.m2[k<2?0:2]);}c.Clear();c.SetLogy(true);mx=0;for(int k=0;k<4;k++)mx=std::max(mx,ht[k]->GetMaximum());ht[0]->SetMaximum(mx*8);ht[0]->SetMinimum(.65);for(int k=0;k<4;k++)ht[k]->Draw(k?"hist same":"hist");TLegend tl(.64,.65,.94,.87);legend(tl);for(int k=0;k<4;k++)tl.AddEntry(ht[k].get(),name4[k],"l")->SetTextColor(kP6Blue+k);tl.Draw();save(c,dir,zs+"_tagged_mass2");histout.cd();for(int k=0;k<4;k++)ht[k]->Write();

    TH1D forced((zs+"_forced_mass2").c_str(),"",180,lo,hi);forced.SetDirectory(nullptr);style(&forced,1);for(auto&v:d)forced.Fill(v.forced*v.forced);c.Clear();c.SetLogy(true);h[0]->SetTitle(("Z = "+std::to_string(z)+", negative control;m_{eff}^{2} [GeV^{2}/c^{4}];Counts").c_str());h[0]->SetMaximum(std::max(forced.GetMaximum(),h[0]->GetMaximum())*5);h[0]->Draw("hist");forced.Draw("hist same");TLegend fl(.59,.73,.94,.87);legend(fl);fl.AddEntry(h[0].get(),"Raw TOF","l")->SetTextColor(colors[0]);fl.AddEntry(&forced,"Forced nearest isotope mass","l")->SetTextColor(colors[1]);fl.Draw();text(.57,.67,"Zero width does not establish separation",kBlack,.029);save(c,dir,zs+"_forced_peak_control");histout.cd();forced.Write();
    TGraph bias[2];for(int k=0;k<2;k++){bias[k].SetLineColor(kP6Blue+k);bias[k].SetMarkerColor(kP6Blue+k);bias[k].SetLineWidth(3);bias[k].SetMarkerStyle(20+k);bias[k].SetMarkerSize(1.3);int point=0;for(int b=0;b<12;b++){double low=-.12+.02*b,high=low+.02;std::vector<double>values;for(auto&v:d){double x=std::log(v.ri/v.R);if(x>=low&&x<high)values.push_back((k?v.rafter:v.rbefore)-v.m2[2]);}if(values.size()<30)continue;bias[k].SetPoint(point++,low+.01,quantile(values,.5));}}
    double xmax=.02,ymax=.02;for(int k=0;k<2;k++)for(int j=0;j<bias[k].GetN();j++){double xx,yy;bias[k].GetPoint(j,xx,yy);xmax=std::max(xmax,std::abs(xx)+.015);ymax=std::max(ymax,1.4*std::abs(yy));}c.Clear();c.SetLogy(false);frame=c.DrawFrame(-xmax,-ymax,xmax,ymax);frame->SetTitle(("Z = "+std::to_string(z)+", controlled rigidity bias;Injected ln(R'/R);Median #Delta m_{eff}^{2} [GeV^{2}/c^{4}]").c_str());style(frame,0);bias[0].Draw("LP same");bias[1].Draw("LP same");TLegend il(.16,.70,.58,.86);legend(il);il.AddEntry(&bias[0],"R bias, fixed calibration","lp")->SetTextColor(kP6Blue);il.AddEntry(&bias[1],"R bias, refitted calibration","lp")->SetTextColor(kP6Blue+1);il.Draw();save(c,dir,zs+"_rigidity_bias");

  }
  std::string merge="pdfunite";for(auto&p:pages)merge+=" \""+p+"\"";merge+=" \""+pdf+"\"";int code=gSystem->Exec(merge.c_str());if(code!=0)throw std::runtime_error("PDF merge failed");histout.Close();std::cout<<"DONE: "<<pages.size()<<" single-chart pages"<<std::endl;
}
