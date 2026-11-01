#include <TCanvas.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct CoreFit {
    int method, isotope, n, npositive, ncore;
    double beta, mean, meanError, sigma, sigmaError, low, high;
};
void exportPage(TCanvas &canvas, const char *pdf, const char *stem,
                std::vector<TCanvas*> &pages, std::vector<std::string> &stems) {
    canvas.RedrawAxis();
    canvas.Print(pdf);
    canvas.Print((std::string("figures/") + stem + ".png").c_str());
    pages.push_back(static_cast<TCanvas*>(canvas.Clone(stem)));
    stems.push_back(stem);
}
}

void plotBetaMassComparison() {
    std::string line;
    std::vector<CoreFit> betaFits, massFits;
    std::ifstream betaInput("beta_mass_beta_fits.tsv");
    std::getline(betaInput, line);
    while (std::getline(betaInput, line)) {
        std::istringstream stream(line);
        CoreFit r;
        if (!(stream >> r.method >> r.beta >> r.n >> r.mean >> r.meanError
              >> r.sigma >> r.sigmaError >> r.ncore >> r.low >> r.high)) {
            std::cerr << "Invalid beta fit row" << std::endl; gSystem->Exit(1);
        }
        betaFits.push_back(r);
    }
    std::ifstream massInput("beta_mass_peak_fits.tsv");
    std::getline(massInput, line);
    while (std::getline(massInput, line)) {
        std::istringstream stream(line);
        CoreFit r;
        if (!(stream >> r.isotope >> r.method >> r.n >> r.npositive >> r.mean
              >> r.meanError >> r.sigma >> r.sigmaError >> r.ncore >> r.low >> r.high)) {
            std::cerr << "Invalid mass fit row" << std::endl; gSystem->Exit(1);
        }
        massFits.push_back(r);
    }
    if (betaFits.size() != 26 || massFits.size() != 4) {
        std::cerr << "Incomplete fit tables" << std::endl; gSystem->Exit(1);
    }
    gROOT->SetBatch(kTRUE);
    gROOT->SetStyle("Pub"); gROOT->ForceStyle();
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetTitleSize(0.045, "XY"); gStyle->SetLabelSize(0.038, "XY");
    gSystem->mkdir("figures", kTRUE);
    TCanvas canvas("beta_mass_comparison", "Frozen beta and rigidity comparison", 1280, 720);
    canvas.SetLeftMargin(0.13); canvas.SetRightMargin(0.035);
    canvas.SetBottomMargin(0.14); canvas.SetTopMargin(0.10);
    const char *pdf = "figures/beta_mass_comparison.pdf";
    std::vector<TCanvas*> pages;
    std::vector<std::string> stems;
    canvas.Print((std::string(pdf) + "[").c_str());
    const char *betaLabels[] = {"Old #beta (linear TOF)", "New #beta (frozen beta-only)"};
    for (int metric = 0; metric < 2; ++metric) {
        canvas.Clear();
        double ymin = 0, ymax = 0;
        for (const auto &r : betaFits) {
            double value = metric ? r.sigma : r.mean;
            double error = metric ? r.sigmaError : r.meanError;
            ymin = std::min(ymin, value-error); ymax = std::max(ymax, value+error);
        }
        double span = ymax-ymin;
        auto frame = canvas.DrawFrame(.68, metric ? 0 : ymin-.10*span, .98, ymax+.27*span);
        frame->GetXaxis()->SetTitle("#beta_{true} (Tracker L5)");
        frame->GetYaxis()->SetTitle(metric ? "#sigma(#beta_{reco}-#beta_{true}) [10^{-3}]"
                                           : "#mu(#beta_{reco}-#beta_{true}) [10^{-3}]");
        frame->GetXaxis()->SetTitleOffset(1.18); frame->GetYaxis()->SetTitleOffset(1.2);
        TLegend legend(.16, .78, .56, .89);
        legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(.032);
        std::vector<TGraphErrors*> graphs;
        for (int method=0; method<2; ++method) {
            auto graph = new TGraphErrors();
            for (const auto &r : betaFits) if (r.method == method) {
                int i=graph->GetN(); graph->SetPoint(i, r.beta, metric ? r.sigma : r.mean);
                graph->SetPointError(i, 0, metric ? r.sigmaError : r.meanError);
            }
            int color=kP6Blue+method;
            graph->SetLineColor(color); graph->SetMarkerColor(color); graph->SetLineWidth(3);
            graph->SetMarkerStyle(20+method); graph->SetMarkerSize(1.1);
            graph->Draw("LPZ SAME");
            legend.AddEntry(graph,betaLabels[method],"lp")->SetTextColor(color);
            graphs.push_back(graph);
        }
        legend.Draw();
        TLatex text; text.SetNDC(); text.SetTextSize(.038);
        text.DrawLatex(.13,.95,metric ? "C12 + C13: beta residual sigma" : "C12 + C13: beta residual mean");
        text.SetTextSize(.026); text.DrawLatex(.62,.95,"Independent test; Gaussian core, 1#sigma fit errors");
        exportPage(canvas,pdf,metric ? "beta_only_residual_sigma" : "beta_only_residual_mean",pages,stems);
        canvas.Clear(); for (auto graph : graphs) delete graph;
    }
    TH1D *hist[4];
    for (int s=0; s<4; ++s) {
        hist[s]=new TH1D(Form("beta_mass_peak_%d",s),"",140,5,19);
        hist[s]->SetDirectory(nullptr); hist[s]->SetLineColor(kP6Blue+s);
        hist[s]->SetLineWidth(3); hist[s]->SetLineStyle(s<2 ? 1 : 2);
    }
    std::ifstream events("beta_mass_events.tsv"); std::getline(events,line);
    int n=0;
    while (std::getline(events,line)) {
        std::istringstream stream(line); int isotope; double raw, corrected;
        if (!(stream >> isotope >> raw >> corrected) || isotope<0 || isotope>1) {
            std::cerr << "Invalid event row" << std::endl; gSystem->Exit(1);
        }
        hist[2*isotope]->Fill(raw); hist[2*isotope+1]->Fill(corrected); ++n;
    }
    if (n!=20775 || hist[0]->GetEntries()!=9920 || hist[2]->GetEntries()!=10855) {
        std::cerr << "Wrong test coverage" << std::endl; gSystem->Exit(1);
    }
    double ymax=0; for (auto h : hist) ymax=std::max(ymax,h->GetMaximum());
    hist[0]->SetMinimum(0); hist[0]->SetMaximum(1.48*ymax);
    hist[0]->GetXaxis()->SetTitle("m [GeV/c^{2}]"); hist[0]->GetYaxis()->SetTitle("Counts");
    hist[0]->GetXaxis()->SetTitleOffset(1.18); hist[0]->GetYaxis()->SetTitleOffset(1.2);
    hist[0]->Draw("HIST"); for(int s=1;s<4;++s) hist[s]->Draw("HIST SAME");
    const char *labels[]={"C12 raw", "C12 corrected", "C13 raw", "C13 corrected"};
    TLatex text; text.SetNDC(); text.SetTextSize(.034);
    text.DrawLatex(.13,.95,"Mass peaks: frozen beta-only + tracker correction");
    text.SetTextSize(.025); text.DrawLatex(.74,.95,"C12: 9,920; C13: 10,855");
    for (const auto &r : massFits) {
        int s=2*r.isotope+r.method;
        double x=.16+.205*s;
        text.SetTextColor(kP6Blue+s); text.SetTextSize(.030);
        text.DrawLatex(x,.855,labels[s]); text.SetTextSize(.027);
        text.DrawLatex(x,.815,Form("#mu = %.3f #pm %.3f",r.mean,r.meanError));
        text.DrawLatex(x,.775,Form("#sigma = %.3f #pm %.3f",r.sigma,r.sigmaError));
    }
    text.SetTextColor(kBlack); text.SetTextSize(.025);
    text.DrawLatex(.16,.725,"Gaussian core [GeV/c^{2}]; bin width 0.1 GeV/c^{2}");
    exportPage(canvas,pdf,"beta_rigidity_mass_peaks",pages,stems);
    canvas.Print((std::string(pdf)+"]").c_str());
    canvas.Clear(); for(auto h : hist) delete h;
    for (size_t i=0; i<pages.size(); ++i) {
        pages[i]->Print((std::string("figures/")+stems[i]+".pdf").c_str());
        delete pages[i];
    }
    std::cout << "Saved beta mean, sigma, and mass peaks from " << n << " test events" << std::endl;
}
