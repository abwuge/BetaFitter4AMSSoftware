#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TColor.h>
#include <TLegend.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

void plotRatioMaterialClosure(const char *input,const char *output){
    gROOT->SetStyle("Pub");gROOT->ForceStyle();gStyle->SetOptStat(0);
    std::ifstream stream(input);int groups;stream>>groups;
    TCanvas c("ratio_closure","",1280,720);
    std::string pdf=std::string(output)+"/material_closure.pdf";c.Print((pdf+"[").c_str());
    for(int k=0;k<groups;++k){
        std::string rad,metric;double xmin,xmax,full_mean;int full_n;stream>>rad>>metric>>xmin>>xmax>>full_n>>full_mean;
        std::vector<double> x[2],y[2],err[2],xe[2];int total[2]={0,0};
        double low=0,high=0;
        for(int m=0;m<2;++m){int n;stream>>n;x[m].resize(n);y[m].resize(n);err[m].resize(n);xe[m].resize(n,0);for(int j=0;j<n;++j){int count;stream>>x[m][j]>>y[m][j]>>err[m][j]>>count;total[m]+=count;low=std::min(low,y[m][j]-err[m][j]);high=std::max(high,y[m][j]+err[m][j]);}}
        c.Clear();c.SetLeftMargin(.13);c.SetRightMargin(.04);c.SetBottomMargin(.15);c.SetTopMargin(.25);
        double pad=std::max(.05,.15*(high-low));TH1D frame(Form("mc_frame_%d",k),"",100,xmin,xmax);frame.SetDirectory(nullptr);
        frame.SetMinimum(metric=="mean"?low-pad:0);frame.SetMaximum(high+pad);
        frame.GetXaxis()->SetTitle("#beta_{MC} at L5");frame.GetYaxis()->SetTitle(metric=="mean"?"Mean #Delta#beta_{ref} #times 10^{3}":"W_{68}(#Delta#beta_{ref}) #times 10^{3}");
        for(auto a:{frame.GetXaxis(),frame.GetYaxis()}){a->SetTitleSize(.044);a->SetLabelSize(.038);}frame.GetXaxis()->SetTitleOffset(1.2);frame.GetYaxis()->SetTitleOffset(1.2);frame.Draw("AXIS");
        TGraphErrors *graphs[2];TLegend legend(.13,.87,.96,.92);legend.SetNColumns(2);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.030);
        for(int m=0;m<2;++m){int color=m==0?kP6Blue:kP6Yellow;graphs[m]=new TGraphErrors(x[m].size(),x[m].data(),y[m].data(),xe[m].data(),err[m].data());graphs[m]->SetLineColor(color);graphs[m]->SetMarkerColor(color);graphs[m]->SetMarkerStyle(20+m);graphs[m]->SetLineWidth(3);graphs[m]->Draw("LP SAME");auto entry=legend.AddEntry(graphs[m],Form("^{%d}C | N = %d in shown bins",12+m,total[m]),"lp");entry->SetTextColor(color);}legend.Draw();
        TLine zero(xmin,0,xmax,0);zero.SetLineColor(kP6Red);zero.SetLineStyle(2);if(metric=="mean")zero.Draw();
        TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.030);
        text.DrawLatex(.13,.948,Form("MC / %s | #Delta#beta_{ref} = #beta_{RICH#rightarrow L5} - #beta_{MC}",rad.c_str()));
        text.DrawLatex(.13,.832,"Fixed ^{12}C transport | 2.5 #leq R < 6 GV | source-bootstrap SE | N #geq 15 per bin");
        text.SetTextSize(.027);text.DrawLatex(.13,.788,Form("All selected events: N = %d | mean #Delta#beta_{ref} #times 10^{3} = %.3f | N outside shown bins = %d",full_n,full_mean,full_n-total[0]-total[1]));
        c.RedrawAxis();c.Print(pdf.c_str());c.Print(Form("%s/%s_%s.png",output,rad.c_str(),metric.c_str()));
        c.Clear();for(auto graph:graphs)delete graph;
    }
    c.Print((pdf+"]").c_str());if(!stream)throw std::runtime_error("Incomplete MC closure points");
}
