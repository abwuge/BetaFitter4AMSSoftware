#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {

struct Species {
    const char *directory;
    const char *stem;
    const char *label;
    int charge;
};

struct Event {
    double generatedRigidity;
    double reconstructedRigidity;
    double invRigidityResidual;
    double linearInvBetaResidual;
    double nonlinearInvBetaResidual;
    double linearMass;
    double nonlinearMass;
    double trueMass;
};

struct RobustValue {
    bool valid = false;
    int count = 0;
    double median = 0;
    double sigma68 = 0;
};

struct BinSummary {
    double low = 0;
    double high = 0;
    RobustValue rigidity;
    RobustValue invRigidity;
    RobustValue linearInvBeta;
    RobustValue nonlinearInvBeta;
    RobustValue linearMass;
    RobustValue nonlinearMass;
    RobustValue trueMass;
};

double quantile(const std::vector<double> &sorted, double probability)
{
    if (sorted.empty())
        return std::numeric_limits<double>::quiet_NaN();
    const double position = probability * (sorted.size() - 1);
    const size_t lower = static_cast<size_t>(std::floor(position));
    const size_t upper = static_cast<size_t>(std::ceil(position));
    if (lower == upper)
        return sorted[lower];
    return sorted[lower] + (position - lower) *
        (sorted[upper] - sorted[lower]);
}

RobustValue summarize(std::vector<double> values)
{
    RobustValue result;
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](double value) { return !std::isfinite(value); }),
                 values.end());
    result.count = values.size();
    if (values.size() < 4)
        return result;
    std::sort(values.begin(), values.end());
    result.median = quantile(values, 0.5);
    result.sigma68 = 0.5 * (quantile(values, 0.84) - quantile(values, 0.16));
    result.valid = std::isfinite(result.median) &&
                   std::isfinite(result.sigma68);
    return result;
}

bool supportedCharge(int charge)
{
    return charge == 2 || charge == 6 || charge == 8;
}

std::vector<Event> loadMatchedEvents(const Species &species,
                                     const std::string &rdstDirectory,
                                     const std::string &betaDirectory,
                                     double rigidityMin,
                                     double rigidityMax)
{
    const std::string rdstName = rdstDirectory + "/" + species.directory +
                                 "/" + species.stem + ".root";
    const std::string betaName = betaDirectory + "/" + species.directory +
                                 "/" + species.stem + ".root";
    std::unique_ptr<TFile> rdstFile(TFile::Open(rdstName.c_str(), "READ"));
    std::unique_ptr<TFile> betaFile(TFile::Open(betaName.c_str(), "READ"));
    if (!rdstFile || rdstFile->IsZombie() || !betaFile || betaFile->IsZombie()) {
        std::cerr << "Cannot open matched inputs for " << species.label << std::endl;
        return {};
    }
    TTree *rdst = dynamic_cast<TTree *>(rdstFile->Get("amstreea"));
    TTree *beta = dynamic_cast<TTree *>(betaFile->Get("betaTree"));
    if (!rdst || !beta) {
        std::cerr << "Missing amstreea or betaTree for " << species.label << std::endl;
        return {};
    }

    const char *requiredRdst[] = {"mch", "mmass", "mevmom1", "tk_q", "tk_qin",
                                  "tk_edep", "tk_pos", "tk_rigidity"};
    for (const char *branch : requiredRdst) {
        if (!rdst->GetBranch(branch)) {
            std::cerr << "Missing required getter-based MC branch " << branch
                      << " in " << rdstName << std::endl;
            return {};
        }
    }

    float mch = 0;
    float mmass = 0;
    float mevmom1[21]{};
    float tkQ[2]{};
    float tkQin[2][3]{};
    float tkEdep[9]{};
    float tkPos[9][3]{};
    float tkRigidity = 0;
    rdst->SetBranchAddress("mch", &mch);
    rdst->SetBranchAddress("mmass", &mmass);
    rdst->SetBranchAddress("mevmom1", mevmom1);
    rdst->SetBranchAddress("tk_q", tkQ);
    rdst->SetBranchAddress("tk_qin", tkQin);
    rdst->SetBranchAddress("tk_edep", tkEdep);
    rdst->SetBranchAddress("tk_pos", tkPos);
    rdst->SetBranchAddress("tk_rigidity", &tkRigidity);

    double betaMc = 0;
    double betaRigidity = 0;
    double betaLinear = 0;
    double betaNonlinear = 0;
    double betaInnerRigidity = 0;
    int betaCharge = 0;
    beta->SetBranchAddress("mcBeta", &betaMc);
    beta->SetBranchAddress("innerRigidity", &betaInnerRigidity);
    beta->SetBranchAddress("betaRigidity", &betaRigidity);
    beta->SetBranchAddress("linearBeta", &betaLinear);
    beta->SetBranchAddress("nonlinearBeta", &betaNonlinear);
    beta->SetBranchAddress("Z", &betaCharge);

    std::vector<Event> events;
    Long64_t betaEntry = 0;
    Long64_t alignmentMismatches = 0;
    for (Long64_t entry = 0; entry < rdst->GetEntries(); ++entry) {
        rdst->GetEntry(entry);
        if (mevmom1[10] == -1000)
            continue;
        const int reconstructedCharge = static_cast<int>(
            (tkQin[0][2] < 2.5 ? tkQ[1] : tkQin[0][2]) + 0.5);
        if (!supportedCharge(reconstructedCharge))
            continue;
        bool trackerComplete = true;
        for (int layer = 1; layer <= 7; ++layer)
            trackerComplete = trackerComplete && tkEdep[layer] > 0 &&
                              std::isfinite(tkEdep[layer]) &&
                              std::isfinite(tkPos[layer][0]) &&
                              std::isfinite(tkPos[layer][1]) &&
                              std::isfinite(tkPos[layer][2]);
        if (!trackerComplete || !(tkPos[4][2] > 0) || !(tkPos[5][2] < 0))
            continue;
        if (betaEntry >= beta->GetEntries()) {
            std::cerr << "betaTree ended before amstreea alignment for "
                      << species.label << std::endl;
            return {};
        }
        beta->GetEntry(betaEntry++);

        const double trueBeta = mevmom1[10] /
            std::sqrt(mevmom1[10] * mevmom1[10] + mmass * mmass);
        if (betaCharge != reconstructedCharge ||
            std::abs(betaInnerRigidity - tkRigidity) >
                1e-6 * std::max(1.0, std::abs(static_cast<double>(tkRigidity))) ||
            std::abs(betaMc - trueBeta) > 1e-6)
            ++alignmentMismatches;

        const int trueCharge = static_cast<int>(std::lround(std::abs(mch)));
        if (trueCharge != species.charge || reconstructedCharge != species.charge ||
            !std::isfinite(tkRigidity) || tkRigidity == 0 ||
            !std::isfinite(betaLinear) || !std::isfinite(betaNonlinear) ||
            !(betaLinear > 0 && betaLinear < 1) ||
            !(betaNonlinear > 0 && betaNonlinear < 1) || mch == 0)
            continue;
        const double generatedRigidity = std::abs(mevmom1[10] / mch);
        if (!(generatedRigidity >= rigidityMin && generatedRigidity < rigidityMax))
            continue;

        const double momentum = std::abs(tkRigidity) * trueCharge;
        Event event;
        event.generatedRigidity = generatedRigidity;
        event.reconstructedRigidity = tkRigidity;
        event.invRigidityResidual = 1.0 / tkRigidity - 1.0 / generatedRigidity;
        event.linearInvBetaResidual = 1.0 / betaLinear - 1.0 / trueBeta;
        event.nonlinearInvBetaResidual = 1.0 / betaNonlinear - 1.0 / trueBeta;
        event.linearMass = momentum *
            std::sqrt(1.0 / (betaLinear * betaLinear) - 1.0);
        event.nonlinearMass = momentum *
            std::sqrt(1.0 / (betaNonlinear * betaNonlinear) - 1.0);
        event.trueMass = mmass;
        events.push_back(event);
    }
    if (betaEntry != beta->GetEntries()) {
        std::cerr << "Unconsumed betaTree entries for " << species.label << ": "
                  << (beta->GetEntries() - betaEntry) << std::endl;
        return {};
    }
    if (alignmentMismatches != 0) {
        std::cerr << "Event alignment mismatches for " << species.label << ": "
                  << alignmentMismatches << std::endl;
        return {};
    }
    std::cout << species.label << ": matched " << betaEntry
              << " beta entries, selected " << events.size()
              << " events in generated-rigidity range" << std::endl;
    return events;
}

std::vector<BinSummary> binEvents(const std::vector<Event> &events, int bins,
                                  double rigidityMin, double rigidityMax)
{
    std::vector<std::vector<double>> rigidity(bins), invRigidity(bins),
        linearInvBeta(bins), nonlinearInvBeta(bins), linearMass(bins),
        nonlinearMass(bins), trueMass(bins);
    const double width = (rigidityMax - rigidityMin) / bins;
    for (const Event &event : events) {
        int bin = static_cast<int>((event.generatedRigidity - rigidityMin) / width);
        if (bin < 0 || bin >= bins)
            continue;
        rigidity[bin].push_back(event.reconstructedRigidity);
        invRigidity[bin].push_back(event.invRigidityResidual);
        linearInvBeta[bin].push_back(event.linearInvBetaResidual);
        nonlinearInvBeta[bin].push_back(event.nonlinearInvBetaResidual);
        linearMass[bin].push_back(event.linearMass);
        nonlinearMass[bin].push_back(event.nonlinearMass);
        trueMass[bin].push_back(event.trueMass);
    }
    std::vector<BinSummary> output(bins);
    for (int bin = 0; bin < bins; ++bin) {
        output[bin].low = rigidityMin + bin * width;
        output[bin].high = output[bin].low + width;
        output[bin].rigidity = summarize(rigidity[bin]);
        output[bin].invRigidity = summarize(invRigidity[bin]);
        output[bin].linearInvBeta = summarize(linearInvBeta[bin]);
        output[bin].nonlinearInvBeta = summarize(nonlinearInvBeta[bin]);
        output[bin].linearMass = summarize(linearMass[bin]);
        output[bin].nonlinearMass = summarize(nonlinearMass[bin]);
        output[bin].trueMass = summarize(trueMass[bin]);
    }
    return output;
}

TGraphErrors *makeGraph(const std::vector<BinSummary> &bins,
                        const RobustValue BinSummary::*member,
                        Color_t color, Style_t marker)
{
    auto *graph = new TGraphErrors();
    for (const BinSummary &bin : bins) {
        const RobustValue &value = bin.*member;
        if (!value.valid)
            continue;
        const int point = graph->GetN();
        graph->SetPoint(point, 0.5 * (bin.low + bin.high), value.median);
        graph->SetPointError(point, 0, value.sigma68);
    }
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetLineWidth(2);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(1.35);
    return graph;
}

void setAxes(TMultiGraph &multi, const char *title, const char *yTitle,
             double rigidityMin, double rigidityMax)
{
    multi.SetTitle(title);
    multi.GetXaxis()->SetTitle("Generated rigidity at AMS center [GV]");
    multi.GetYaxis()->SetTitle(yTitle);
    multi.GetXaxis()->SetLimits(rigidityMin, rigidityMax);
    multi.GetXaxis()->SetTitleSize(0.046);
    multi.GetYaxis()->SetTitleSize(0.046);
    multi.GetXaxis()->SetLabelSize(0.038);
    multi.GetYaxis()->SetLabelSize(0.038);
    multi.GetYaxis()->SetTitleOffset(1.38);
}

void preparePad()
{
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.10);
    gPad->SetGridy();
}

void savePage(TCanvas &canvas, const TString &base, const TString &pngName,
              const TString &pdfName, int page, int pages)
{
    canvas.SaveAs(base + "/" + pngName);
    TString target = base + "/" + pdfName;
    if (page == 0)
        target += "(";
    else if (page == pages - 1)
        target += ")";
    canvas.Print(target);
}

} // namespace

void plotMcRigidityTruthComparison(
    const char *rdstDirectory =
        "/lustre02/user/yachen/public/hxwu/amsd69n_getrigidity_mc_small",
    const char *betaDirectory = "results/mc_getrigidity_small",
    const char *outputDirectory = "macro/results/mc_truth",
    double rigidityMin = 0.8,
    double rigidityMax = 5.0)
{
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetEndErrorSize(7);
    gSystem->mkdir(outputDirectory, true);

    const Species species[] = {
        {"He4", "82259493", "^{4}He", 2},
        {"C12", "82334649", "^{12}C", 6},
        {"O16", "82357004", "^{16}O", 8},
    };
    constexpr int speciesCount = sizeof(species) / sizeof(species[0]);
    std::vector<Event> events[speciesCount];
    size_t minimumEntries = std::numeric_limits<size_t>::max();
    for (int index = 0; index < speciesCount; ++index) {
        events[index] = loadMatchedEvents(species[index], rdstDirectory,
                                          betaDirectory, rigidityMin, rigidityMax);
        if (events[index].empty())
            return;
        minimumEntries = std::min(minimumEntries, events[index].size());
    }
    // Keep roughly six events per bin in the smallest species sample.  The
    // preliminary MC requested here is intentionally small, so a high fixed
    // minimum would suppress all robust summaries for C and O.
    const int bins = std::clamp(static_cast<int>(minimumEntries / 6), 3, 14);
    const double binWidth = (rigidityMax - rigidityMin) / bins;
    std::cout << "Using " << bins << " bins (width " << binWidth
              << " GV), based on minimum selected sample " << minimumEntries
              << std::endl;

    std::vector<BinSummary> summaries[speciesCount];
    for (int index = 0; index < speciesCount; ++index)
        summaries[index] = binEvents(events[index], bins, rigidityMin, rigidityMax);

    const Color_t rigidityColor = TColor::GetColor("#45364B");
    const Color_t truthColor = TColor::GetColor("#788AA3");
    const Color_t linearColor = TColor::GetColor("#447604");
    const Color_t nonlinearColor = TColor::GetColor("#E55812");
    const TString outputBase(outputDirectory);
    const TString pdfName("mc_truth_comparison.pdf");
    TCanvas canvas("mcTruthCanvas", "", 2100, 720);
    canvas.Divide(speciesCount, 1, 0.001, 0.001);

    // Page 1: reconstructed rigidity versus generated center truth.
    for (int index = 0; index < speciesCount; ++index) {
        canvas.cd(index + 1);
        preparePad();
        auto *graph = makeGraph(summaries[index], &BinSummary::rigidity,
                                rigidityColor, 20);
        auto *multi = new TMultiGraph();
        multi->Add(graph, "P");
        multi->Draw("A");
        setAxes(*multi, Form("%s: rigidity response", species[index].label),
                "Reconstructed rigidity from GetRigidity [GV]",
                rigidityMin, rigidityMax);
        multi->SetMinimum(rigidityMin);
        multi->SetMaximum(rigidityMax);
        auto *truth = new TLine(rigidityMin, rigidityMin, rigidityMax, rigidityMax);
        truth->SetLineColor(truthColor);
        truth->SetLineStyle(7);
        truth->SetLineWidth(2);
        truth->Draw("SAME");
        auto *legend = new TLegend(0.18, 0.74, 0.58, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.032);
        legend->AddEntry(graph, "median #pm #sigma_{68}", "lp");
        legend->AddEntry(truth, "R_{rec}=R_{gen}", "l");
        legend->Draw();
    }
    savePage(canvas, outputBase, "rigidity_response.png", pdfName, 0, 4);

    // Page 2: inverse-rigidity residual around generated center truth.
    canvas.Clear();
    canvas.Divide(speciesCount, 1, 0.001, 0.001);
    for (int index = 0; index < speciesCount; ++index) {
        canvas.cd(index + 1);
        preparePad();
        auto *graph = makeGraph(summaries[index], &BinSummary::invRigidity,
                                rigidityColor, 20);
        auto *multi = new TMultiGraph();
        multi->Add(graph, "P");
        multi->Draw("A");
        setAxes(*multi, Form("%s: inverse-rigidity residual", species[index].label),
                "1/R_{rec} - 1/R_{gen} [GV^{-1}]", rigidityMin, rigidityMax);
        auto *truth = new TLine(rigidityMin, 0, rigidityMax, 0);
        truth->SetLineColor(truthColor);
        truth->SetLineStyle(7);
        truth->SetLineWidth(2);
        truth->Draw("SAME");
    }
    savePage(canvas, outputBase, "inverse_rigidity_residual.png", pdfName, 1, 4);

    // Page 3: linear and nonlinear inverse-beta residuals.
    canvas.Clear();
    canvas.Divide(speciesCount, 1, 0.001, 0.001);
    for (int index = 0; index < speciesCount; ++index) {
        canvas.cd(index + 1);
        preparePad();
        auto *linear = makeGraph(summaries[index], &BinSummary::linearInvBeta,
                                 linearColor, 20);
        auto *nonlinear = makeGraph(summaries[index],
                                    &BinSummary::nonlinearInvBeta,
                                    nonlinearColor, 21);
        auto *multi = new TMultiGraph();
        multi->Add(linear, "P");
        multi->Add(nonlinear, "P");
        multi->Draw("A");
        setAxes(*multi, Form("%s: inverse-beta residual", species[index].label),
                "1/#beta_{rec} - 1/#beta_{MC,center}", rigidityMin, rigidityMax);
        auto *truth = new TLine(rigidityMin, 0, rigidityMax, 0);
        truth->SetLineColor(truthColor);
        truth->SetLineStyle(7);
        truth->SetLineWidth(2);
        truth->Draw("SAME");
        auto *legend = new TLegend(0.18, 0.73, 0.54, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.032);
        legend->AddEntry(linear, "Linear #beta", "lp");
        legend->AddEntry(nonlinear, "Nonlinear #beta", "lp");
        legend->Draw();
    }
    savePage(canvas, outputBase, "inverse_beta_residual.png", pdfName, 2, 4);

    // Page 4: rigidity-derived mass using linear and nonlinear beta.
    canvas.Clear();
    canvas.Divide(speciesCount, 1, 0.001, 0.001);
    for (int index = 0; index < speciesCount; ++index) {
        canvas.cd(index + 1);
        preparePad();
        auto *linear = makeGraph(summaries[index], &BinSummary::linearMass,
                                 linearColor, 20);
        auto *nonlinear = makeGraph(summaries[index], &BinSummary::nonlinearMass,
                                    nonlinearColor, 21);
        auto *truthGraph = makeGraph(summaries[index], &BinSummary::trueMass,
                                     truthColor, 24);
        auto *multi = new TMultiGraph();
        multi->Add(linear, "P");
        multi->Add(nonlinear, "P");
        multi->Add(truthGraph, "L");
        multi->Draw("A");
        setAxes(*multi, Form("%s: reconstructed mass", species[index].label),
                "Rigidity-derived mass [GeV/c^{2}]", rigidityMin, rigidityMax);
        auto *legend = new TLegend(0.18, 0.69, 0.56, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.030);
        legend->AddEntry(linear, "Linear #beta", "lp");
        legend->AddEntry(nonlinear, "Nonlinear #beta", "lp");
        legend->AddEntry(truthGraph, "Generated mass", "l");
        legend->Draw();
    }
    savePage(canvas, outputBase, "mass_response.png", pdfName, 3, 4);

    std::ofstream table((outputBase + "/mc_truth_summary.tsv").Data());
    table << "species\tbin_low_GV\tbin_high_GV\tcount"
             "\tRrec_median_GV\tRrec_sigma68_GV"
             "\tdelta_invR_median_GVinv\tdelta_invR_sigma68_GVinv"
             "\tdelta_invBeta_linear_median\tdelta_invBeta_linear_sigma68"
             "\tdelta_invBeta_nonlinear_median\tdelta_invBeta_nonlinear_sigma68"
             "\tmass_linear_median_GeV\tmass_linear_sigma68_GeV"
             "\tmass_nonlinear_median_GeV\tmass_nonlinear_sigma68_GeV"
             "\ttrue_mass_GeV\n";
    table << std::setprecision(9);
    for (int index = 0; index < speciesCount; ++index) {
        for (const BinSummary &bin : summaries[index]) {
            table << species[index].label << '\t' << bin.low << '\t' << bin.high
                  << '\t' << bin.rigidity.count << '\t'
                  << bin.rigidity.median << '\t' << bin.rigidity.sigma68 << '\t'
                  << bin.invRigidity.median << '\t' << bin.invRigidity.sigma68 << '\t'
                  << bin.linearInvBeta.median << '\t' << bin.linearInvBeta.sigma68 << '\t'
                  << bin.nonlinearInvBeta.median << '\t' << bin.nonlinearInvBeta.sigma68 << '\t'
                  << bin.linearMass.median << '\t' << bin.linearMass.sigma68 << '\t'
                  << bin.nonlinearMass.median << '\t' << bin.nonlinearMass.sigma68 << '\t'
                  << bin.trueMass.median << '\n';
        }
    }

    TFile output((outputBase + "/mc_truth_comparison.root").Data(), "RECREATE");
    for (int index = 0; index < speciesCount; ++index) {
        auto *directory = output.mkdir(species[index].directory);
        directory->cd();
        makeGraph(summaries[index], &BinSummary::rigidity, rigidityColor, 20)
            ->Write("rigidity_response");
        makeGraph(summaries[index], &BinSummary::invRigidity, rigidityColor, 20)
            ->Write("inverse_rigidity_residual");
        makeGraph(summaries[index], &BinSummary::linearInvBeta, linearColor, 20)
            ->Write("inverse_beta_linear");
        makeGraph(summaries[index], &BinSummary::nonlinearInvBeta, nonlinearColor, 21)
            ->Write("inverse_beta_nonlinear");
        makeGraph(summaries[index], &BinSummary::linearMass, linearColor, 20)
            ->Write("mass_linear");
        makeGraph(summaries[index], &BinSummary::nonlinearMass, nonlinearColor, 21)
            ->Write("mass_nonlinear");
        makeGraph(summaries[index], &BinSummary::trueMass, truthColor, 24)
            ->Write("mass_truth");
    }
    output.Close();
    std::cout << "Wrote " << outputBase
              << "/{mc_truth_comparison.pdf,*.png,mc_truth_summary.tsv,mc_truth_comparison.root}"
              << std::endl;
}
