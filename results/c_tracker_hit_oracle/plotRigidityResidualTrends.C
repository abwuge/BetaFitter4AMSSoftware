#include <TCanvas.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH1.h>
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
struct ResidualBin {
    int isotope, method, n, ncore;
    double low, high, center, mean, meanError, sigma, sigmaError, fitLow, fitHigh;
};
}

void plotRigidityResidualTrends() {
    // Reuse the frozen test sample and its stored conditional Gaussian fits.
    std::ifstream input("residual_trends.tsv");
    std::string line, observable;
    std::getline(input, line);
    std::vector<ResidualBin> rows;
    while (std::getline(input, line)) {
        std::istringstream stream(line);
        ResidualBin row;
        if (!(stream >> observable >> row.isotope >> row.method >> row.low >> row.high
              >> row.center >> row.n >> row.mean >> row.meanError >> row.sigma
              >> row.sigmaError >> row.ncore >> row.fitLow >> row.fitHigh)) {
            std::cerr << "Invalid residual fit row" << std::endl;
            gSystem->Exit(1);
        }
        if (observable == "rigidity") rows.push_back(row);
    }
    if (rows.empty()) {
        std::cerr << "No rigidity fit results" << std::endl;
        gSystem->Exit(1);
    }
    gROOT->SetBatch(kTRUE);
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTitleSize(0.045, "XY");
    gStyle->SetLabelSize(0.038, "XY");
    gSystem->mkdir("figures", kTRUE);
    TCanvas canvas("rigidity_trends", "Rigidity residual trends", 1280, 720);
    canvas.SetLeftMargin(0.12);
    canvas.SetRightMargin(0.035);
    canvas.SetBottomMargin(0.14);
    canvas.SetTopMargin(0.09);
    const char *labels[] = {"C12 raw", "C12 corrected", "C13 raw", "C13 corrected"};
    const char *pdf = "figures/rigidity_residual_trends.pdf";
    std::vector<TCanvas*> pages;
    canvas.Print((std::string(pdf) + "[").c_str());
    for (int metric = 0; metric < 2; ++metric) {
        canvas.Clear();
        double ymin = 0, ymax = 0;
        for (const auto &r : rows) {
            const double value = metric ? r.sigma : r.mean;
            const double error = metric ? r.sigmaError : r.meanError;
            ymin = std::min(ymin, value-error);
            ymax = std::max(ymax, value+error);
        }
        const double span = ymax-ymin;
        auto frame = canvas.DrawFrame(0.68, metric ? 0 : ymin-0.08*span,
                                      0.98, ymax+0.30*span);
        frame->GetXaxis()->SetTitle("#beta_{true} (Tracker L5)");
        frame->GetYaxis()->SetTitle(metric ? "#sigma[100 ln(R/R_{true})] [%]"
                                           : "#mu[100 ln(R/R_{true})] [%]");
        frame->GetXaxis()->SetTitleOffset(1.18);
        frame->GetYaxis()->SetTitleOffset(1.18);
        TLegend legend(0.15, 0.765, 0.58, 0.895);
        legend.SetNColumns(2);
        legend.SetBorderSize(0);
        legend.SetFillStyle(0);
        legend.SetTextSize(0.033);
        std::vector<TGraphErrors*> graphs;
        for (int series=0; series<4; ++series) {
            auto graph = new TGraphErrors();
            for (const auto &r : rows) {
                if (r.isotope != series/2 || r.method != series%2) continue;
                int i = graph->GetN();
                graph->SetPoint(i, r.center, metric ? r.sigma : r.mean);
                graph->SetPointError(i, 0, metric ? r.sigmaError : r.meanError);
            }
            if (!graph->GetN()) {
                std::cerr << "Missing series " << labels[series] << std::endl;
                gSystem->Exit(1);
            }
            const int color = kP6Blue+series;
            graph->SetLineColor(color);
            graph->SetMarkerColor(color);
            graph->SetLineWidth(3);
            graph->SetLineStyle(series<2 ? 1 : 2);
            graph->SetMarkerStyle(series<2 ? 20 : 21);
            graph->SetMarkerSize(1.05);
            graph->Draw("LPZ SAME");
            legend.AddEntry(graph, labels[series], "lp")->SetTextColor(color);
            graphs.push_back(graph);
        }
        legend.Draw();
        TLatex text;
        text.SetNDC();
        text.SetTextSize(0.039);
        text.DrawLatex(0.12, 0.95, metric ? "Rigidity residual: core resolution"
                                        : "Rigidity residual: mean bias");
        text.SetTextSize(0.027);
        text.DrawLatex(0.62, 0.95, "N = 20,775; Gaussian core, fit errors");
        canvas.RedrawAxis();
        canvas.Print(pdf);
        std::string stem = metric ? "figures/rigidity_residual_sigma" : "figures/rigidity_residual_mean";
        canvas.Print((stem+".png").c_str());
        pages.push_back(static_cast<TCanvas*>(canvas.Clone(metric ? "sigma_page" : "mean_page")));
        canvas.Clear();
        for (auto graph : graphs) delete graph;
    }
    canvas.Print((std::string(pdf) + "]").c_str());
    for (int metric=0; metric<2; ++metric) {
        pages[metric]->Print(metric ? "figures/rigidity_residual_sigma.pdf"
                                    : "figures/rigidity_residual_mean.pdf");
        delete pages[metric];
    }
    std::cout << "Saved two rigidity trend pages from " << rows.size() << " stored fits" << std::endl;
}
