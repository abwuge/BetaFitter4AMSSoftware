#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

void plotExternalCarbonRatio(const char *input, const char *output) {
    gROOT->SetStyle("Pub"); gROOT->ForceStyle();
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    std::ifstream stream(input); int n=0; stream>>n;
    if(n<1) throw std::runtime_error("Empty external isotope ratio");
    std::vector<double> x(n),y(n),error(n),zero(n,0.);
    double ymin=0.,ymax=0.;
    for(int i=0;i<n;++i){stream>>x[i]>>y[i]>>error[i];ymax=std::max(ymax,y[i]+error[i]);ymin=std::min(ymin,y[i]-error[i]);}
    if(!stream) throw std::runtime_error("Incomplete external isotope ratio");
    TCanvas canvas("external_carbon_ratio","",1280,720);
    canvas.SetLeftMargin(.12);canvas.SetRightMargin(.035);
    canvas.SetBottomMargin(.16);canvas.SetTopMargin(.19);
    TGraphErrors graph(n,x.data(),y.data(),zero.data(),error.data());
    graph.SetLineColor(kP6Blue);graph.SetMarkerColor(kP6Blue);
    graph.SetLineWidth(2);graph.SetMarkerStyle(20);graph.SetMarkerSize(.8);
    for(int linear=0;linear<2;++linear){
        canvas.Clear();canvas.SetLogx(!linear);
        TH1D frame(Form("external_ratio_frame_%d",linear),"",100,linear?0.:.1,linear?x.back()*1.05:40.);
        frame.SetDirectory(nullptr);frame.SetMinimum(ymin<0.?ymin-.03*ymax:0.);frame.SetMaximum(ymax*1.15);
        frame.GetXaxis()->SetTitle("Kinetic energy per nucleon [GeV/n]");
        frame.GetYaxis()->SetTitle("^{13}C / ^{12}C");
        for(auto axis:{frame.GetXaxis(),frame.GetYaxis()}){axis->SetTitleSize(.044);axis->SetLabelSize(.038);}
        frame.GetXaxis()->SetTitleOffset(1.25);frame.GetYaxis()->SetTitleOffset(1.15);
        frame.Draw("AXIS");graph.Draw("LP SAME");
        TLegend legend(.13,.845,.94,.91);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.032);
        legend.AddEntry(&graph,"External combined ratio: TOF + RICH NaF + RICH AGL","lp")->SetTextColor(kP6Blue);
        legend.Draw();
        TLatex label;label.SetNDC();label.SetTextFont(42);label.SetTextSize(.029);
        label.DrawLatex(.12,.953,Form("Supplied ratio used for training | %d points | original vertical errors",n));
        canvas.RedrawAxis();
        std::string base=std::string(output)+"/external_ratio_energy"+(linear?"_linear":"");
        canvas.Print((base+".pdf").c_str());canvas.Print((base+".png").c_str());
        canvas.Clear();
    }
}
