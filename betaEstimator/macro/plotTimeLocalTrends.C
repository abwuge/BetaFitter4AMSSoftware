#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TColor.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

void plotTimeLocalTrends(const char *input, const char *output, bool center_reference=false, const char *radiator="") {
    gROOT->SetStyle("Pub"); gROOT->ForceStyle();
    gStyle->SetOptStat(0); gStyle->SetOptFit(0);
    std::ifstream in(input); int groups; in>>groups;
    TCanvas c("trends","",1280,720);
    const int colors[]={kP6Blue,kP6Yellow};
    for(int group=0;group<groups;++group) {
        std::string axis; int n;double xmin,xmax;in>>axis>>n>>xmin>>xmax;
        std::vector<double> x(n),width[2],bias[2],gain(n),lower(n),upper(n);
        for(int m=0;m<2;++m){width[m].resize(n);bias[m].resize(n);}
        for(int i=0;i<n;++i)in>>x[i]>>width[0][i]>>width[1][i]>>bias[0][i]>>bias[1][i]>>gain[i]>>lower[i]>>upper[i];
        const bool time=axis=="time"||axis=="yearly";
        std::string base=std::string(output)+"/"+axis;
        const std::string combined=base+"_comparison.pdf";
        if(!center_reference)c.Print((combined+"[").c_str());
        for(int metric=0;metric<3;++metric){
            if(time&&metric<2)continue;
            c.Clear();c.SetLogx(axis=="rigidity");
            c.SetLeftMargin(.12);c.SetRightMargin(.035);c.SetBottomMargin(.16);c.SetTopMargin(.17);
            double ymin=metric==0?1e30:0,ymax=-1e30;
            for(int i=0;i<n;++i){
                if(x[i]<xmin||x[i]>xmax)continue;
                if(metric<2)for(int m=0;m<2;++m){double v=metric==0?width[m][i]:bias[m][i];ymin=std::min(ymin,v);ymax=std::max(ymax,v);}
                else{ymin=std::min(ymin,std::min(gain[i],lower[i]));ymax=std::max(ymax,std::max(gain[i],upper[i]));}
            }
            if(metric>0){ymin=std::min(0.,ymin);ymax=std::max(0.,ymax);}
            double pad=std::max((ymax-ymin)*.12,.2);ymin-=pad;ymax+=pad;if(metric==0)ymin=std::max(0.,ymin);
            TH1D frame(Form("frame_%d_%d",group,metric),"",100,xmin,xmax);frame.SetDirectory(nullptr);
            frame.SetMinimum(ymin);frame.SetMaximum(ymax);
            const char *xtitle=time?"Target period [UTC]":axis=="rigidity"?"Rigidity [GV]":"#beta_{RICH}";
            frame.GetXaxis()->SetTitle(xtitle);
            frame.GetYaxis()->SetTitle(metric==0?"W_{68}(#Delta#beta) #times 10^{3}":metric==1?"Mean #Delta#beta #times 10^{3}":"W_{68} improvement [%]");
            for(auto a:{frame.GetXaxis(),frame.GetYaxis()}){a->SetLabelSize(.038);a->SetTitleSize(.044);}
            frame.GetXaxis()->SetTitleOffset(1.25);frame.GetYaxis()->SetTitleOffset(1.15);
            if(time){frame.GetXaxis()->SetTimeDisplay(1);frame.GetXaxis()->SetTimeOffset(0,"gmt");frame.GetXaxis()->SetTimeFormat("%Y");frame.GetXaxis()->SetNdivisions(508);}
            frame.Draw("AXIS");
            TLegend legend(.12,.865,.965,.92);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextFont(42);legend.SetTextSize(.029);legend.SetNColumns(2);
            std::vector<TGraph*> graphs;
            if(metric<2){
                for(int m=0;m<2;++m){auto &y=metric==0?width[m]:bias[m];auto g=new TGraph(n,x.data(),y.data());g->SetLineColor(colors[m]);g->SetMarkerColor(colors[m]);g->SetLineWidth(3);g->SetMarkerStyle(m==0?20:21);g->SetMarkerSize(.8);g->Draw("LP SAME");auto e=legend.AddEntry(g,m==0?"Linear TOF":center_reference?"External-ratio neighbour":"Neighbour-trained model","lp");e->SetTextColor(colors[m]);graphs.push_back(g);}
            }else{
                if(axis!="time"){
                    auto band=new TGraph(2*n);
                    for(int i=0;i<n;++i){band->SetPoint(i,x[i],lower[i]);band->SetPoint(2*n-1-i,x[i],upper[i]);}
                    band->SetFillColorAlpha(colors[1],.2);band->SetLineWidth(0);band->Draw("F SAME");graphs.push_back(band);
                }
                auto g=new TGraph(n,x.data(),gain.data());g->SetLineColor(colors[1]);g->SetMarkerColor(colors[1]);g->SetLineWidth(time?2:3);g->SetMarkerStyle(21);g->SetMarkerSize(axis=="time"?.35:.8);g->Draw("LP SAME");auto e=legend.AddEntry(g,axis=="time"?"Neighbour-trained / 15 days":axis=="yearly"?"Neighbour-trained / yearly pooled":center_reference?"External-ratio neighbour":"Neighbour-trained model","lp");e->SetTextColor(colors[1]);graphs.push_back(g);
            }
            TLine zero(xmin,0,xmax,0);zero.SetLineColor(colors[0]);zero.SetLineStyle(2);zero.SetLineWidth(1);
            if(metric>0){zero.Draw();if(metric==2){auto e=legend.AddEntry(&zero,"Linear TOF (0%)","l");e->SetTextColor(colors[0]);}}
            legend.Draw();
            TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.029);
            if(metric==2&&axis!="time")text.DrawLatex(.12,.956,center_reference?Form("%s | Band: conditional 95%% run bootstrap",radiator):"Band: conditional 95% run bootstrap");
            else if(metric<2)text.DrawLatex(.12,.956,center_reference?Form("%s | #Delta#beta = #hat{#beta}_{L5} - #beta_{RICH#rightarrow L5} | bins with N #geq 20",radiator):"#Delta#beta = #hat{#beta}_{RICH} - #beta_{RICH}    |    bins with N #geq 20");
            c.RedrawAxis();if(!center_reference)c.Print(combined.c_str());
            const char *suffix=metric==0?"width":metric==1?"bias":"gain";
            std::string individual=base+"_"+suffix;
            c.Print((individual+".pdf").c_str());c.Print((individual+".png").c_str());
            if(metric==0)c.Print((base+"_comparison.png").c_str());
            c.Clear();for(auto g:graphs)delete g;
        }
        if(!center_reference)c.Print((combined+"]").c_str());
    }
    if(!in)throw std::runtime_error("Incomplete trend input");
}
