#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
const char *stages[] = {"tkadc_base", "tkadc_adcshape", "tkadc_adcstrip", "tkadc_morph", "tkadc_shape_only"};
const char *labels[] = {"Tracker baseline", "+ ADC shape", "+ ADC shape + strip address", "+ ADC shape + strip + Q ratios", "New features only"};
struct Gain { std::string comparison, stage; int bin; double direct, median, low, high; };
void legend_entry(TLegend &legend, TObject *object, const char *label, int color, const char *option) {
    auto entry = legend.AddEntry(object, label, option);
    entry->SetTextColor(color);
}
}

void plot_tkadc_validation() {
    gROOT->SetBatch(kTRUE);
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.045, "XY");
    gStyle->SetLabelSize(0.037, "XY");
    gSystem->mkdir("figures", kTRUE);
    std::ifstream input("tkadc_gain_plot.tsv");
    std::string line;
    std::getline(input, line);
    std::vector<Gain> rows;
    while (std::getline(input, line)) {
        std::istringstream stream(line);
        Gain row;
        if (stream >> row.comparison >> row.stage >> row.bin >> row.direct >> row.median >> row.low >> row.high) rows.push_back(row);
    }
    if (rows.size() != 36) { std::cerr << "Incomplete gain table" << std::endl; gSystem->Exit(1); }
    TFile file("tkadc_test_predictions.root");
    auto tree = dynamic_cast<TTree*>(file.Get("events"));
    if (!tree || tree->GetEntries() != 9920) { std::cerr << "Wrong test sample" << std::endl; gSystem->Exit(1); }
    TCanvas canvas("tkadc_canvas", "C12 strip/cluster validation", 1280, 720);
    canvas.SetLeftMargin(0.12); canvas.SetRightMargin(0.035);
    canvas.SetTopMargin(0.10); canvas.SetBottomMargin(0.13);
    const char *pdf = "figures/tkadc_validation.pdf";
    canvas.Print((std::string(pdf) + "[").c_str());
    const double centers[] = {2.75, 3.5, 4.5, 5.5};
    for (std::string comparison : {"raw", "base"}) {
        canvas.Clear(); canvas.SetLogy(false);
        double lo = 0, hi = 0;
        for (const auto &row : rows) if (row.comparison == comparison && !(comparison == "base" && row.stage == "tkadc_shape_only")) {
            lo = std::min(lo, row.low); hi = std::max(hi, row.high);
        }
        double span = std::max(hi-lo, 2.0);
        auto frame = canvas.DrawFrame(2.45, lo-0.15*span, 6.05, hi+0.50*span);
        frame->SetTitle(comparison == "raw" ? "C12: rigidity width vs raw reconstruction" : "C12: added information vs tracker baseline");
        frame->GetXaxis()->SetTitle("R_{reco} [GV]");
        frame->GetYaxis()->SetTitle(comparison == "raw" ? "W68 narrowing vs raw [%]" : "W68 narrowing vs tracker baseline [%]");
        TLegend legend(0.15, 0.73, 0.90, 0.88);
        legend.SetNColumns(2); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.031);
        std::vector<TGraphAsymmErrors*> graphs;
        for (int s=0; s<5; ++s) {
            if (comparison == "base" && (s == 0 || s == 4)) continue;
            auto graph = new TGraphAsymmErrors();
            for (const auto &row : rows) if (row.comparison == comparison && row.stage == stages[s]) {
                int point = graph->GetN();
                graph->SetPoint(point, centers[row.bin], row.median);
                graph->SetPointError(point, 0, 0, row.median-row.low, row.high-row.median);
            }
            int color = kP6Blue+s+1;
            graph->SetLineColor(color); graph->SetMarkerColor(color);
            graph->SetLineWidth(3); graph->SetMarkerStyle(20+s); graph->SetMarkerSize(1.2);
            graph->Draw("LP SAME"); legend_entry(legend, graph, labels[s], color, "lp"); graphs.push_back(graph);
        }
        legend.Draw();
        TLatex text; text.SetNDC(); text.SetTextSize(0.028);
        text.DrawLatex(0.62, 0.93, "N = 9920; source bootstrap 95% CI");
        text.SetTextSize(0.033);
        text.DrawLatex(0.12, 0.93, comparison == "raw" ? "C12 rigidity: width" : "C12 rigidity: added information");
        canvas.Print(pdf);
        std::string stem = "figures/tkadc_gain_"+comparison;
        canvas.Print((stem+".png").c_str());
        for (auto graph : graphs) delete graph;
    }
    float R, R_true, u;
    float prediction[5];
    tree->SetBranchAddress("R", &R); tree->SetBranchAddress("R_true", &R_true); tree->SetBranchAddress("u", &u);
    for (int s=0; s<5; ++s) tree->SetBranchAddress(stages[s], &prediction[s]);
    for (int observable=0; observable<2; ++observable) {
        canvas.Clear(); canvas.SetLogy(true);
        std::vector<TH1D*> histograms;
        for (int m=0; m<5; ++m) {
            auto hist = new TH1D(Form("tkadc_hist_%d_%d", observable, m), "", observable == 0 ? 160 : 110,
                                 observable == 0 ? -0.4 : 40, observable == 0 ? 0.4 : 260);
            hist->SetDirectory(nullptr);
            hist->SetLineColor(kP6Blue+m); hist->SetLineWidth(3); histograms.push_back(hist);
        }
        for (Long64_t entry=0; entry<tree->GetEntries(); ++entry) {
            tree->GetEntry(entry);
            for (int m=0; m<5; ++m) {
                double corrected = R * std::exp(m == 0 ? 0.0 : -double(prediction[m-1]));
                double value = observable == 0 ? std::log(corrected/R_true) : std::pow(6.0*corrected, 2)*(double(u)*u-1);
                histograms[m]->Fill(value);
            }
        }
        double maxcount=1;
        for (auto hist : histograms) maxcount = std::max(maxcount, hist->GetMaximum());
        auto first=histograms[0];
        first->SetTitle(observable == 0 ? "C12: tracker rigidity residuals" : "C12: mass squared, measured TOF beta");
        first->GetXaxis()->SetTitle(observable == 0 ? "ln(R/R_{true})" : "m^{2} [GeV^{2}/c^{4}]");
        first->GetYaxis()->SetTitle("Counts"); first->SetMinimum(0.8); first->SetMaximum(maxcount*8);
        first->Draw("HIST");
        for (int m=1; m<5; ++m) histograms[m]->Draw("HIST SAME");
        TLegend legend(0.62, 0.70, 0.96, 0.89);
        legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.028);
        const char *methods[] = {"Raw", "Tracker baseline", "+ ADC shape", "+ ADC shape + strip", "+ ADC shape + strip + Q ratios"};
        for (int m=0; m<5; ++m) legend_entry(legend, histograms[m], methods[m], kP6Blue+m, "l");
        legend.Draw();
        TLatex text; text.SetNDC(); text.SetTextSize(0.030); text.DrawLatex(0.14, 0.84, "N = 9920");
        text.SetTextSize(0.033); text.DrawLatex(0.12, 0.93, observable == 0 ? "C12 rigidity residuals" : "C12 mass peaks, measured TOF beta");
        canvas.Print(pdf);
        std::string stem = observable == 0 ? "figures/tkadc_rigidity_residual" : "figures/tkadc_mass2_peaks";
        canvas.Print((stem+".png").c_str());
        for (auto hist : histograms) delete hist;
    }
    canvas.Print((std::string(pdf) + "]").c_str());
    std::cout << "VALIDATION_FIGURES_DONE" << std::endl;
}
