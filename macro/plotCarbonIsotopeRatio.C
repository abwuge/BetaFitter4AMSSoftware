#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

void plotCarbonIsotopeRatio(const char *input, const char *output) {
    gROOT->SetStyle("Pub"); gROOT->ForceStyle();
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    std::ifstream in(input); int n; in >> n;
    std::vector<double> x(n), y[3]; for (auto &v:y) v.resize(n);
    for(int i=0;i<n;++i) in>>x[i]>>y[0][i]>>y[1][i]>>y[2][i];
    int nw; in>>nw; std::vector<double> wx(nw),wl(nw),wh(nw),wy(nw);
    std::vector<int> counts(nw),bound(nw);
    for(int i=0;i<nw;++i) in>>wx[i]>>wl[i]>>wh[i]>>wy[i]>>counts[i]>>bound[i];
    int nr;in>>nr;TGraphAsymmErrors refs[4];
    for(int i=0;i<nr;++i){int group;double e,lo,hi,v,el,eh;in>>group>>e>>lo>>hi>>v>>el>>eh;
        auto &g=refs[group];int j=g.GetN();g.SetPoint(j,e,v);g.SetPointError(j,e-lo,hi-e,el,eh);}
    if(!in) throw std::runtime_error("Incomplete ratio input");
    const int colors[]={kP8Blue,kP8Blue+1,kP8Blue+2,kP8Blue+3,kP8Blue+4,kP8Blue+5,kP8Blue+6};
    TCanvas c("carbon_ratio","",1280,720);
    c.SetLeftMargin(.115);c.SetRightMargin(.035);c.SetBottomMargin(.16);c.SetTopMargin(.25);
    c.SetLogx(false);c.SetLogy(false);
    TH1D frame("energy_frame","",100,.04,2.2);frame.SetDirectory(nullptr);
    frame.SetMinimum(0);frame.SetMaximum(.4);
    frame.GetXaxis()->SetTitle("Kinetic energy per nucleon [GeV/n]");
    frame.GetYaxis()->SetTitle("^{13}C / ^{12}C");
    for(auto a:{frame.GetXaxis(),frame.GetYaxis()}){a->SetTitleSize(.046);a->SetLabelSize(.038);}
    frame.GetXaxis()->SetTitleOffset(1.25);frame.GetYaxis()->SetTitleOffset(1.05);
    frame.Draw("AXIS");
    TLegend legend(.12,.775,.97,.94);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetNColumns(2);legend.SetTextFont(42);legend.SetTextSize(.029);
    const char *names[]={"Selected NewCalib + ADC + noise","Native timing + Gaussian R","Native timing + Student-t R"};
    TGraph g[3];
    for(int j=0;j<3;++j){g[j]=TGraph(n,x.data(),y[j].data());g[j].SetLineColor(colors[j]);g[j].SetLineWidth(j==1?4:3);g[j].SetLineStyle(1);g[j].Draw("L SAME");legend.AddEntry(&g[j],names[j],"l")->SetTextColor(colors[j]);}
    const char *referenceNames[]={"Near-Earth balloon/space (1973-2000)","Voyager 1/2, heliospheric (1994/1996)","Voyager 1, VLISM (Cummings 2025)","Balloon 1977 (Buffington 1978)"};
    for(int j=0;j<4;++j){refs[j].SetMarkerStyle(21+j);refs[j].SetMarkerSize(1.3);refs[j].SetMarkerColor(colors[j+3]);refs[j].SetLineColor(colors[j+3]);refs[j].SetLineWidth(2);refs[j].Draw("P SAME");legend.AddEntry(&refs[j],referenceNames[j],"pe")->SetTextColor(colors[j+3]);}
    legend.Draw();
    TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.027);
    text.DrawLatex(.12,.968,"Frozen flight fits: selected-sample population, equal energy; measured R = 2.5-6 GV");
    c.RedrawAxis();std::string base=std::string(output)+"/carbon_ratio_energy";
    c.Print((base+".pdf").c_str());c.Print((base+".png").c_str());
    c.Clear();c.SetTopMargin(.24);
    TH1D wide("wide_frame","",2000,0,105);wide.SetDirectory(nullptr);wide.SetMinimum(0);wide.SetMaximum(.3);
    wide.GetXaxis()->SetTitle("Kinetic energy [GeV/n] (flight: ^{12}C equivalent)");wide.GetYaxis()->SetTitle("^{13}C / ^{12}C");
    for(auto a:{wide.GetXaxis(),wide.GetYaxis()}){a->SetTitleSize(.044);a->SetLabelSize(.038);}
    wide.GetXaxis()->SetTitleOffset(1.25);wide.GetYaxis()->SetTitleOffset(1.08);wide.Draw("AXIS");
    TGraph trend(nw,wx.data(),wy.data());
    trend.SetLineColor(kP6Blue);trend.SetLineWidth(3);trend.Draw("L SAME");
    TGraphAsymmErrors fit,limit;
    for(int i=0;i<nw;++i){auto &g=bound[i]?limit:fit;int j=g.GetN();g.SetPoint(j,wx[i],wy[i]);g.SetPointError(j,wx[i]-wl[i],wh[i]-wx[i],0,0);}
    fit.SetLineColor(kP6Blue);fit.SetMarkerColor(kP6Blue);fit.SetMarkerStyle(20);fit.SetMarkerSize(1.3);fit.SetLineWidth(3);fit.Draw("P SAME");
    limit.SetLineColor(kP6Blue+1);limit.SetMarkerColor(kP6Blue+1);limit.SetMarkerStyle(26);limit.SetMarkerSize(1.8);limit.SetLineWidth(4);limit.Draw("P SAME");
    TLegend leg(.12,.79,.97,.935);leg.SetNColumns(2);leg.SetTextFont(42);leg.SetTextSize(.03);leg.SetBorderSize(0);leg.SetFillStyle(0);
    leg.AddEntry(&fit,"Flight: fitted ratio","pe")->SetTextColor(kP6Blue);
    leg.AddEntry(&limit,"At imposed f_{13} = 0.001 lower bound","pe")->SetTextColor(kP6Blue+1);
    for(int j=0;j<4;++j){refs[j].SetLineColor(kP6Blue+2+j);refs[j].SetMarkerColor(kP6Blue+2+j);refs[j].Draw("P SAME");leg.AddEntry(&refs[j],referenceNames[j],"pe")->SetTextColor(kP6Blue+2+j);}leg.Draw();
    text.DrawLatex(.12,.955,"Flight: same-R bin proportions; reference: equal-energy measurements from CRDB");
    c.RedrawAxis();base=std::string(output)+"/wide_training_ratio_energy";
    c.Print((base+".pdf").c_str());c.Print((base+".png").c_str());
    wide.GetXaxis()->SetRangeUser(0,3.5);c.Modified();c.Update();
    base=std::string(output)+"/wide_training_ratio_energy_low";
    c.Print((base+".pdf").c_str());c.Print((base+".png").c_str());
    // Show the fitted range and only the first bin of the trailing boundary run.
    int lastFree = -1;
    for(int i=0;i<nw;++i) if(!bound[i]) lastFree=i;
    int shown = lastFree+1<nw ? lastFree+2 : nw;
    c.Clear();
    TH1D focus("focus_frame","",100,wl.front(),wh[shown-1]);focus.SetDirectory(nullptr);
    focus.SetMinimum(0);focus.SetMaximum(.3);
    focus.GetXaxis()->SetTitle("Kinetic energy [GeV/n] (flight: ^{12}C equivalent)");
    focus.GetYaxis()->SetTitle("^{13}C / ^{12}C");
    for(auto a:{focus.GetXaxis(),focus.GetYaxis()}){a->SetTitleSize(.044);a->SetLabelSize(.038);}
    focus.GetXaxis()->SetTitleOffset(1.25);focus.GetYaxis()->SetTitleOffset(1.08);focus.Draw("AXIS");
    TGraph shortTrend(shown,wx.data(),wy.data());
    shortTrend.SetLineColor(kP6Blue);shortTrend.SetLineWidth(3);shortTrend.Draw("L SAME");
    TGraphAsymmErrors shortFit,shortLimit;
    for(int i=0;i<shown;++i){auto &g=bound[i]?shortLimit:shortFit;int j=g.GetN();g.SetPoint(j,wx[i],wy[i]);g.SetPointError(j,wx[i]-wl[i],wh[i]-wx[i],0,0);}
    shortFit.SetLineColor(kP6Blue);shortFit.SetMarkerColor(kP6Blue);shortFit.SetMarkerStyle(20);shortFit.SetMarkerSize(1.3);shortFit.SetLineWidth(3);shortFit.Draw("P SAME");
    shortLimit.SetLineColor(kP6Blue+1);shortLimit.SetMarkerColor(kP6Blue+1);shortLimit.SetMarkerStyle(26);shortLimit.SetMarkerSize(1.8);shortLimit.SetLineWidth(3);shortLimit.Draw("P SAME");
    for(int j=0;j<4;++j) refs[j].Draw("P SAME");
    leg.Draw();
    text.DrawLatex(.12,.955,"Flight: same-R proportions; trailing lower-bound bins condensed");
    text.DrawLatex(.58,.62,"R #geq 20 GV: all at f_{13} = 0.001 bound");
    c.RedrawAxis();base=std::string(output)+"/wide_training_ratio_energy_focus";
    c.Print((base+".pdf").c_str());c.Print((base+".png").c_str());

}
