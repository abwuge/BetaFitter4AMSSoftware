#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TColor.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>

void plotMaterialMCPairs(const char *input, const char *output, double maximum_R=6.) {
    gROOT->SetStyle("Pub");gROOT->ForceStyle();gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    std::ifstream stream(input);int groups=0;stream>>groups;
    TCanvas canvas("material_mc_pairs","",1280,720);
    canvas.SetLeftMargin(.13);canvas.SetRightMargin(.035);canvas.SetBottomMargin(.16);canvas.SetTopMargin(.21);
    for(int group=0;group<groups;++group){
        std::string rad;int n;stream>>rad>>n;
        std::vector<double> x(n),y(n);std::vector<int> isotope(n);
        for(int i=0;i<n;++i)stream>>x[i]>>y[i]>>isotope[i];
        for(int zoom=0;zoom<(rad=="AGL"?2:1);++zoom){
            const double observed_min=std::min(*std::min_element(x.begin(),x.end()),*std::min_element(y.begin(),y.end()));
            const double observed_max=std::max(*std::max_element(x.begin(),x.end()),*std::max_element(y.begin(),y.end()));
            const double xmin=zoom?.955:std::min(rad=="AGL"?.76:.78,observed_min-.005);
            const double xmax=zoom?(maximum_R>6.?1.005:.995):std::max(1.,observed_max+.005);
            std::vector<double> xx[2],yy[2],rr[2];
            double rmin=0.,rmax=0.;
            for(int i=0;i<n;++i){
                if(x[i]<xmin||x[i]>xmax||y[i]<xmin||y[i]>xmax)continue;
                const int m=isotope[i];if(m<0||m>1)throw std::runtime_error("Unknown carbon isotope");
                xx[m].push_back(x[i]);yy[m].push_back(y[i]);rr[m].push_back((y[i]-x[i])*1000.);
                rmin=std::min(rmin,rr[m].back());rmax=std::max(rmax,rr[m].back());
            }
            const int shown=xx[0].size()+xx[1].size();
            for(int residual=0;residual<2;++residual){
                canvas.Clear();
                TH1D frame(Form("pair_frame_%d_%d_%d",group,zoom,residual),"",100,xmin,xmax);frame.SetDirectory(nullptr);
                const double padding=std::max(.15,.08*(rmax-rmin));
                frame.SetMinimum(residual?rmin-padding:xmin);frame.SetMaximum(residual?rmax+padding:xmax);
                frame.GetXaxis()->SetTitle("#beta_{MC} at L5");
                frame.GetYaxis()->SetTitle(residual?"(#beta_{RICH#rightarrow L5} - #beta_{MC}) #times 10^{3}":"#beta_{RICH#rightarrow L5}");
                for(auto axis:{frame.GetXaxis(),frame.GetYaxis()}){axis->SetTitleSize(.044);axis->SetLabelSize(.038);}
                frame.GetXaxis()->SetTitleOffset(1.25);frame.GetYaxis()->SetTitleOffset(1.2);frame.Draw("AXIS");
                TLine ideal(xmin,residual?0.:xmin,xmax,residual?0.:xmax);
                ideal.SetLineColor(kP6Red);ideal.SetLineStyle(2);ideal.SetLineWidth(2);ideal.Draw();
                TGraph graph[2];TLegend legend(.13,.865,.96,.92);legend.SetNColumns(3);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.029);
                for(int m=0;m<2;++m){
                    auto &values=residual?rr[m]:yy[m];const int color=m==0?kP6Blue:kP6Yellow;
                    graph[m]=TGraph(xx[m].size(),xx[m].data(),values.data());
                    graph[m].SetMarkerStyle(maximum_R>6.?1:20+m);graph[m].SetMarkerSize(rad=="NaF"?.55:.8);
                    graph[m].SetMarkerColorAlpha(color,.55);graph[m].Draw("P SAME");
                    legend.AddEntry(&graph[m],Form("^{%d}C | N = %d",12+m,int(xx[m].size())),"p")->SetTextColor(color);
                }
                legend.AddEntry(&ideal,residual?"Zero residual":"y = x","l")->SetTextColor(kP6Red);legend.Draw();
                TLatex label;label.SetNDC();label.SetTextFont(42);label.SetTextSize(.029);
                label.DrawLatex(.13,.954,Form("MC / %s | #beta_{RICH#rightarrow L5} = (1/#beta_{RICH} - #Delta u_{material})^{-1}",rad.c_str()));
                label.DrawLatex(.13,.825,Form("Fixed ^{12}C transport | 2.5 #leq R < %.0f GV | %s | shown %d / %d; outside %d",maximum_R,zoom?"zoom":"full range",shown,n,n-shown));
                canvas.RedrawAxis();
                std::string base=std::string(output)+"/"+rad+(residual?"_mc_residual_pairs":"_mc_beta_pairs")+(zoom?"_zoom":"");
                canvas.Print((base+".png").c_str());canvas.Print((base+".pdf").c_str());canvas.Clear();
            }
        }
    }
    if(!stream)throw std::runtime_error("Incomplete MC pairs");
}
