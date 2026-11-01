#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TF1.h>
#include <TColor.h>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>

namespace material_plot {
double quant(std::vector<double> v,double p){
    std::sort(v.begin(),v.end());double t=p*(v.size()-1);int i=int(t);
    return v[i]+(t-i)*(v[std::min(i+1,int(v.size()-1))]-v[i]);
}
void axes(TH1 *h,const char *x,const char *y){
    h->SetTitle("");h->GetXaxis()->SetTitle(x);h->GetYaxis()->SetTitle(y);
    for(auto a:{h->GetXaxis(),h->GetYaxis()}){a->SetTitleSize(.044);a->SetLabelSize(.037);}
    h->GetXaxis()->SetTitleOffset(1.18);h->GetYaxis()->SetTitleOffset(1.25);
}
}

void plotMaterialComparison(const char *output,bool external_ratio=false){
    using namespace material_plot;
    gROOT->SetStyle("Pub");gROOT->ForceStyle();gStyle->SetOptStat(0);gStyle->SetOptFit(0);
    auto blue=gROOT->GetColor(kP6Blue);
    double stop[3]={0,.45,1},r[3]={1,(1+blue->GetRed())/2,blue->GetRed()},g[3]={1,(1+blue->GetGreen())/2,blue->GetGreen()},b[3]={1,(1+blue->GetBlue())/2,blue->GetBlue()};
    TColor::CreateGradientColorTable(3,stop,r,g,b,100);gStyle->SetNumberContours(100);
    TFile f(Form("%s/comparisons.root",output));
    TCanvas c("material_comparison","",1280,720);
    const std::string pdf=std::string(output)+"/material_comparison.pdf";
    c.Print((pdf+"[").c_str());
    auto save=[&](std::string name){c.RedrawAxis();c.Print(pdf.c_str());c.Print((std::string(output)+"/"+name+".png").c_str());};
    auto setup=[&](bool density=false){c.Clear();c.SetLogy(false);c.SetLogz(density);c.SetRightMargin(density?.14:.04);c.SetLeftMargin(.13);c.SetBottomMargin(.15);c.SetTopMargin(.16);};
    TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.032);
    const int color[2]={kP6Blue,kP6Yellow};
    std::ofstream fits(Form("%s/gaussian_fits.tsv",output));fits<<"sample\tmethod\tn\tmean\tmean_error\tsigma\tsigma_error\tfit_status\n";
    for(std::string radiator:{"NaF","AGL"}){
        TTree *t=nullptr;f.GetObject(("flight_"+radiator).c_str(),t);if(!t)throw std::runtime_error("Missing flight tree");
        double ref,u,raw,rr,nr;t->SetBranchAddress("u_reference",&ref);t->SetBranchAddress("u_new",&u);t->SetBranchAddress("u_raw",&raw);
        t->SetBranchAddress("raw_residual",&rr);t->SetBranchAddress("new_residual",&nr);
        std::vector<double> refs,us,samples[2];
        for(Long64_t i=0;i<t->GetEntries();++i){t->GetEntry(i);refs.push_back(ref);us.push_back(u);samples[0].push_back(1000*rr);samples[1].push_back(1000*nr);}
        double lo=quant(refs,.002),hi=quant(refs,.998),yl=quant(us,.002),yh=quant(us,.998);
        setup(true);TH2D scatter("flight_scatter","",100,lo,hi,100,yl,yh);scatter.SetDirectory(nullptr);
        for(size_t i=0;i<refs.size();++i)scatter.Fill(refs[i],us[i]);
        axes(&scatter,"u_{RICH}-u_{material}","u_{new}");scatter.Draw("COLZ");
        TLine identity(std::max(lo,yl),std::max(lo,yl),std::min(hi,yh),std::min(hi,yh));identity.SetLineColor(color[1]);identity.SetLineStyle(2);identity.Draw();
        text.DrawLatex(.13,.945,Form("Flight / %s   |   N = %lld",radiator.c_str(),t->GetEntries()));
        text.DrawLatex(.13,.89,"u = 1/#beta;  u_{material} = #Delta u_{C12};  central 99.6% axis ranges");save("flight_"+radiator+"_u_comparison");
        setup();double low=std::min(quant(samples[0],.005),quant(samples[1],.005)),high=std::max(quant(samples[0],.995),quant(samples[1],.995));
        TH1D *hist[2];TF1 *fit[2];double peak=0;
        for(int j=0;j<2;++j){
            hist[j]=new TH1D(Form("residual_%s_%d",radiator.c_str(),j),"",120,low,high);hist[j]->SetDirectory(nullptr);
            for(double v:samples[j])hist[j]->Fill(v);peak=std::max(peak,hist[j]->GetMaximum());
            double med=quant(samples[j],.5),sig=.741301109252801*(quant(samples[j],.75)-quant(samples[j],.25));
            TH1D core(Form("core_%d",j),"",120,med-3*sig,med+3*sig);core.SetDirectory(nullptr);for(double v:samples[j])core.Fill(v);
            fit[j]=new TF1(Form("gaus_%s_%d",radiator.c_str(),j),"gaus",med-3*sig,med+3*sig);fit[j]->SetParameters(core.GetMaximum(),med,sig);
            int status=core.Fit(fit[j],"Q0RS");
            fit[j]->SetParameter(0,fit[j]->GetParameter(0)*hist[j]->GetBinWidth(1)/core.GetBinWidth(1));
            fit[j]->SetLineColor(color[j]);fit[j]->SetLineStyle(2);fit[j]->SetLineWidth(3);
            fits<<radiator<<'\t'<<j<<'\t'<<samples[j].size()<<'\t'<<fit[j]->GetParameter(1)<<'\t'<<fit[j]->GetParError(1)<<'\t'<<fit[j]->GetParameter(2)<<'\t'<<fit[j]->GetParError(2)<<'\t'<<status<<'\n';
            hist[j]->SetLineColor(color[j]);hist[j]->SetLineWidth(3);axes(hist[j],"[u_{TOF}-(u_{RICH}-u_{material})] #times 10^{3}","Counts");
        }
        hist[0]->SetMaximum(peak*1.48);hist[0]->Draw("HIST");hist[1]->Draw("HIST SAME");
        for(int j=0;j<2;++j){fit[j]->Draw("SAME");double x=.16+.43*j;text.SetTextColor(color[j]);text.DrawLatex(x,.79,j==0?"Linear TOF":external_ratio?"External-ratio neighbour":"Neighbour-trained");text.DrawLatex(x,.74,Form("#mu = %.3f #pm %.3f",fit[j]->GetParameter(1),fit[j]->GetParError(1)));text.DrawLatex(x,.69,Form("#sigma = %.3f #pm %.3f",fit[j]->GetParameter(2),fit[j]->GetParError(2)));}
        text.SetTextColor(kBlack);text.DrawLatex(.13,.945,Form("Flight / %s   |   N = %lld",radiator.c_str(),t->GetEntries()));text.DrawLatex(.13,.89,"Gaussian core fit; displayed range: pooled central 99%");save("flight_"+radiator+"_residual");
        c.Clear();for(int j=0;j<2;++j){delete hist[j];delete fit[j];}

        TTree *mc=nullptr;f.GetObject(("mc_"+radiator).c_str(),mc);if(!mc)throw std::runtime_error("Missing MC tree");
        double bm,res;mc->SetBranchAddress("beta_mc",&bm);mc->SetBranchAddress("residual",&res);
        std::vector<double> bx,by;for(Long64_t i=0;i<mc->GetEntries();++i){mc->GetEntry(i);bx.push_back(bm);by.push_back(1000*res);}
        double xmin=std::max(.5,quant(bx,0)-.005),xmax=std::min(1.005,quant(bx,1)+.005);
        setup(true);double ymin=quant(by,.005),ymax=quant(by,.995);TH2D closure("closure","",60,xmin,xmax,100,ymin,ymax);closure.SetDirectory(nullptr);
        for(size_t i=0;i<by.size();++i)closure.Fill(bx[i],by[i]);axes(&closure,"#beta_{MC} at L5","(u_{RICH}-u_{material}-u_{MC}) #times 10^{3}");closure.Draw("COLZ");
        TLine zero(xmin,0,xmax,0);zero.SetLineColor(color[1]);zero.SetLineStyle(2);zero.SetLineWidth(2);zero.Draw();
        text.DrawLatex(.13,.945,Form("MC / %s   |   N = %lld",radiator.c_str(),mc->GetEntries()));text.DrawLatex(.13,.89,"Fixed C12 mean material correction; central 99% shown");save("mc_"+radiator+"_material_closure");
        std::vector<double> xx,mean,err,width,xe;
        double trend_min=radiator=="NaF"?xmin:.955,trend_max=radiator=="NaF"?xmax:1.;
        const int bins=radiator=="NaF"?12:9;
        for(int j=0;j<bins;++j){double l=trend_min+(trend_max-trend_min)*j/bins,h=trend_min+(trend_max-trend_min)*(j+1)/bins;std::vector<double> v;double sx=0,sum=0,sum2=0;for(size_t i=0;i<by.size();++i)if(bx[i]>=l&&bx[i]<h){v.push_back(by[i]);sx+=bx[i];sum+=by[i];sum2+=by[i]*by[i];}if(v.size()<15)continue;double n=v.size(),m=sum/n;xx.push_back(sx/n);xe.push_back(0);mean.push_back(m);err.push_back(std::sqrt(std::max(0.,sum2/n-m*m)/n));width.push_back((quant(v,.84)-quant(v,.16))/2);}
        for(int metric=0;metric<2;++metric){
            setup();auto &values=metric==0?mean:width;double a=std::min(0.,*std::min_element(values.begin(),values.end())),z=std::max(0.,*std::max_element(values.begin(),values.end()));if(metric==0)for(size_t i=0;i<values.size();++i){a=std::min(a,values[i]-err[i]);z=std::max(z,values[i]+err[i]);}double pad=std::max(.1,(z-a)*.2);
            TH1D frame("closure_frame","",100,trend_min,trend_max);frame.SetDirectory(nullptr);frame.SetMinimum(metric==0?a-pad:0);frame.SetMaximum(z+pad);axes(&frame,"#beta_{MC} at L5",metric==0?"Mean (u_{RICH}-u_{material}-u_{MC}) #times 10^{3}":"W_{68}(u_{RICH}-u_{material}-u_{MC}) #times 10^{3}");frame.Draw("AXIS");
            TGraphErrors graph(xx.size(),xx.data(),values.data(),xe.data(),metric==0?err.data():xe.data());graph.SetLineColor(color[0]);graph.SetMarkerColor(color[0]);graph.SetMarkerStyle(20);graph.SetLineWidth(3);graph.Draw("LP SAME");if(metric==0){zero.SetX1(trend_min);zero.SetX2(trend_max);zero.Draw();}
            text.DrawLatex(.13,.945,Form("MC / %s   |   N = %lld",radiator.c_str(),mc->GetEntries()));text.DrawLatex(.13,.89,metric==0?"Error bars: standard error of mean; bins with N #geq 15":"W_{68} = (q_{84}-q_{16})/2; bins with N #geq 15");save("mc_"+radiator+(metric==0?"_material_bias":"_material_width"));
        }
    }
    c.Print((pdf+"]").c_str());
}
