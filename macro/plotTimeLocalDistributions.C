#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TLatex.h>
#include <TColor.h>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

void plotTimeLocalDistributions(const char *input, const char *output) {
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    std::ifstream in(input);
    std::ofstream fits(std::string(output)+"/fit_parameters.csv");
    fits << "axis,view,bin,method,n,fit_low,fit_high,mu,mu_error,sigma,sigma_error,status,chi2,ndf\n";
    fits.precision(17);
    int pages; in >> pages;
    TCanvas canvas("distributions", "", 1280, 800);
    const int colors[] = {kP6Blue, kP6Yellow};
    const char *labels[] = {"Linear TOF", "Neighbour-trained model"};
    std::string active;
    for (int page=0; page<pages; ++page) {
        std::string axis, view; int bin, n, nb;
        double low, high, xmin, xmax;
        in >> axis >> view >> bin >> low >> high >> n >> nb >> xmin >> xmax;
        const std::string pdf=std::string(output)+"/"+axis+"_distributions_"+view+".pdf";
        if (pdf!=active) {
            if (!active.empty()) canvas.Print((active+"]").c_str());
            active=pdf; canvas.Print((active+"[").c_str());
        }
        canvas.Clear(); canvas.SetLogy(view=="full");
        canvas.SetLeftMargin(.10); canvas.SetRightMargin(.035);
        canvas.SetBottomMargin(.14); canvas.SetTopMargin(.30);
        std::vector<TH1D*> hist; std::vector<TF1*> curves;
        double maximum=1;
        for (int m=0;m<2;++m) {
            auto h=new TH1D(Form("display_%d_%d",page,m),"",nb,xmin*1000,xmax*1000);
            h->SetDirectory(nullptr); h->SetLineColor(colors[m]); h->SetLineWidth(3);
            for (int j=1;j<=nb;++j) {double count;in>>count; h->SetBinContent(j,count);}
            maximum=std::max(maximum,h->GetMaximum()); hist.push_back(h);
        }
        for (int m=0;m<2;++m) {
            auto h=hist[m];h->SetMinimum(view=="full"?.5:0);
            h->SetMaximum(maximum*(view=="full"?2:1.15));
            h->GetXaxis()->SetTitle("(#hat{#beta}_{RICH} - #beta_{RICH}) #times 10^{3}");
            h->GetYaxis()->SetTitle("Counts");
            for (auto a:{h->GetXaxis(),h->GetYaxis()}) {a->SetLabelSize(.038);a->SetTitleSize(.044);}
            h->GetXaxis()->SetTitleOffset(1.25);h->GetYaxis()->SetTitleOffset(1.0);
            h->Draw(m==0?"HIST":"HIST SAME");
        }
        TLatex text; text.SetNDC(); text.SetTextFont(42); text.SetTextSize(.030);
        text.DrawLatex(.10,.956,axis=="rigidity" ? Form("%.5g #leq R < %.5g GV",low,high) : Form("%.5g #leq #beta_{RICH} < %.5g",low,high));
        text.DrawLatex(.10,.91,Form("N = %d   |   %s range   |   bin width = %.4g #times 10^{-3}",n,view.c_str(),(xmax-xmin)*1000/nb));
        for (int m=0;m<2;++m) {
            int nf; double flo,fhi,median,sigma; in>>nf>>flo>>fhi>>median>>sigma;
            TH1D fitHist(Form("core_%d_%d",page,m),"",nf,flo*1000,fhi*1000);
            fitHist.SetDirectory(nullptr);
            for (int j=1;j<=nf;++j) {double count;in>>count;fitHist.SetBinContent(j,count);fitHist.SetBinError(j,std::sqrt(count));}
            auto f=new TF1(Form("gauss_%d_%d",page,m),"gaus",flo*1000,fhi*1000);
            f->SetParameters(fitHist.GetMaximum(),median*1000,sigma*1000);
            f->SetParLimits(2,sigma*1000*.1,sigma*1000*5);
            auto result=fitHist.Fit(f,"QRSNLI");
            int status=result;
            const bool valid=n>=20 && status==0 && result.Get() && result->CovMatrixStatus()>=2;
            text.SetTextColor(colors[m]); double x=m==0?.10:.55;
            text.DrawLatex(x,.86,labels[m]);
            if (valid) {
                text.DrawLatex(x,.818,Form("#mu = %.4f #pm %.4f  [10^{-3}]",f->GetParameter(1),f->GetParError(1)));
                text.DrawLatex(x,.776,Form("#sigma = %.4f #pm %.4f  [10^{-3}]",f->GetParameter(2),f->GetParError(2)));
            } else text.DrawLatex(x,.818,n<20 ? "Gaussian fit: insufficient entries (N < 20)" : Form("Gaussian fit unavailable (status %d)",status));
            text.SetTextSize(.027);
            text.DrawLatex(x,.733,Form("Gaussian core fit: [%.3f, %.3f] #times 10^{-3}",flo*1000,fhi*1000));
            text.SetTextSize(.030);
            fits<<axis<<','<<view<<','<<bin<<','<<(m==0?"raw":"local")<<','<<n<<','<<flo<<','<<fhi<<',';
            if (valid) fits<<f->GetParameter(1)/1000<<','<<f->GetParError(1)/1000<<','<<f->GetParameter(2)/1000<<','<<f->GetParError(2)/1000;
            else fits<<"nan,nan,nan,nan";
            fits<<','<<(n<20?-1:status)<<','<<f->GetChisquare()<<','<<f->GetNDF()<<'\n';
            f->SetParameter(0,f->GetParameter(0)*(xmax-xmin)/nb/((fhi-flo)/nf));
            f->SetLineColor(colors[m]);f->SetLineStyle(2);f->SetLineWidth(3);
            if (valid) curves.push_back(f); else delete f;
        }
        for (auto f:curves) f->Draw("SAME");
        canvas.RedrawAxis();canvas.Print(pdf.c_str());
        canvas.Clear();for(auto h:hist)delete h;for(auto f:curves)delete f;
    }
    if (!active.empty())canvas.Print((active+"]").c_str());
    if (!in) throw std::runtime_error("Incomplete distribution input");
}
