#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
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
#include <memory>
#include <string>

namespace {

struct Species {
    int charge;
    const char *label;
    double mass;
};

struct PeakFit {
    bool valid = false;
    int status = -1;
    double mean = 0;
    double meanError = 0;
    double sigma = 0;
    double sigmaError = 0;
};

PeakFit fitPeak(TH1D *hist, const char *name, double referenceMass)
{
    PeakFit result;
    const int first = hist->FindBin(0.72 * referenceMass);
    const int last = hist->FindBin(1.28 * referenceMass);
    int peakBin = first;
    for (int bin = first + 1; bin <= last; ++bin)
        if (hist->GetBinContent(bin) > hist->GetBinContent(peakBin))
            peakBin = bin;

    const double mode = hist->GetBinCenter(peakBin);
    const double firstHalfWidth = 0.12 * referenceMass;
    TF1 firstFit(Form("%s_first", name), "gaus", mode - firstHalfWidth,
                 mode + firstHalfWidth);
    firstFit.SetParameters(hist->GetBinContent(peakBin), mode,
                           0.06 * referenceMass);
    int status = hist->Fit(&firstFit, "QNR");
    if (status != 0 || !std::isfinite(firstFit.GetParameter(1)) ||
        !std::isfinite(firstFit.GetParameter(2)) || firstFit.GetParameter(2) <= 0)
        return result;

    const double center = firstFit.GetParameter(1);
    const double sigma = std::abs(firstFit.GetParameter(2));
    const double halfWidth = std::clamp(2.2 * sigma,
                                        0.025 * referenceMass,
                                        0.14 * referenceMass);
    TF1 finalFit(name, "gaus", center - halfWidth, center + halfWidth);
    finalFit.SetParameters(firstFit.GetParameter(0), center, sigma);
    status = hist->Fit(&finalFit, "QNR");

    result.status = status;
    result.mean = finalFit.GetParameter(1);
    result.meanError = finalFit.GetParError(1);
    result.sigma = std::abs(finalFit.GetParameter(2));
    result.sigmaError = finalFit.GetParError(2);
    result.valid = status == 0 && std::isfinite(result.mean) &&
                   std::isfinite(result.sigma) && result.sigma > 0 &&
                   result.mean > 0.65 * referenceMass &&
                   result.mean < 1.35 * referenceMass;
    return result;
}

void styleHistogram(TH1D *hist, Color_t color, int lineStyle)
{
    hist->SetStats(false);
    hist->SetLineColor(color);
    hist->SetLineWidth(3);
    hist->SetLineStyle(lineStyle);
    hist->GetXaxis()->SetTitle("Rigidity-derived mass  [GeV/c^{2}]");
    hist->GetYaxis()->SetTitle("Normalized events / (GeV/c^{2})");
    hist->GetXaxis()->SetTitleSize(0.050);
    hist->GetYaxis()->SetTitleSize(0.050);
    hist->GetXaxis()->SetLabelSize(0.042);
    hist->GetYaxis()->SetLabelSize(0.042);
    hist->GetYaxis()->SetTitleOffset(1.45);
}

} // namespace

void plotRigidityMassPeaks(
    const char *inputName = "results/14.root",
    const char *outputDirectory = "macro/results/data_p000-009",
    double rigidityMin = 1.0,
    double rigidityMax = 3.0)
{
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gSystem->mkdir(outputDirectory, true);

    std::unique_ptr<TFile> input(TFile::Open(inputName, "READ"));
    if (!input || input->IsZombie()) {
        std::cerr << "Cannot open " << inputName << std::endl;
        return;
    }
    TTree *tree = dynamic_cast<TTree *>(input->Get("betaTree"));
    if (!tree || tree->GetEntries() == 0) {
        std::cerr << "Missing or empty betaTree in " << inputName << std::endl;
        return;
    }

    // Dominant A/Z ~= 2 species match the mass hypothesis used by BetaFitter.
    const Species species[] = {
        {2, "^{4}He", 3.7273794},
        {6, "^{12}C", 11.174864},
        {8, "^{16}O", 14.895109},
    };
    constexpr int nSpecies = sizeof(species) / sizeof(species[0]);
    const Color_t linearColor = TColor::GetColor("#45364B");
    const Color_t nonlinearColor = TColor::GetColor("#5F6077");
    const Color_t referenceColor = TColor::GetColor("#788AA3");

    double rigidity = 0;
    double linearBeta = 0;
    double nonlinearBeta = 0;
    int charge = 0;
    tree->SetBranchAddress("innerRigidity", &rigidity);
    tree->SetBranchAddress("linearBeta", &linearBeta);
    tree->SetBranchAddress("nonlinearBeta", &nonlinearBeta);
    tree->SetBranchAddress("Z", &charge);

    std::unique_ptr<TH1D> linear[nSpecies];
    std::unique_ptr<TH1D> nonlinear[nSpecies];
    Long64_t selected[nSpecies]{};
    for (int index = 0; index < nSpecies; ++index) {
        const double low = 0.55 * species[index].mass;
        const double high = 1.45 * species[index].mass;
        linear[index] = std::make_unique<TH1D>(
            Form("hMassLinearZ%d", species[index].charge), "", 180, low, high);
        nonlinear[index] = std::make_unique<TH1D>(
            Form("hMassNonlinearZ%d", species[index].charge), "", 180, low, high);
        linear[index]->Sumw2();
        nonlinear[index]->Sumw2();
    }

    const Long64_t entries = tree->GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);
        if (!std::isfinite(rigidity) || rigidity <= rigidityMin ||
            rigidity >= rigidityMax || !std::isfinite(linearBeta) ||
            !std::isfinite(nonlinearBeta) || linearBeta <= 0 ||
            linearBeta >= 1 || nonlinearBeta <= 0 || nonlinearBeta >= 1)
            continue;

        int index = -1;
        for (int candidate = 0; candidate < nSpecies; ++candidate)
            if (charge == species[candidate].charge)
                index = candidate;
        if (index < 0)
            continue;

        const double momentum = std::abs(rigidity) * charge;
        const double linearMass = momentum *
            std::sqrt(1.0 / (linearBeta * linearBeta) - 1.0);
        const double nonlinearMass = momentum *
            std::sqrt(1.0 / (nonlinearBeta * nonlinearBeta) - 1.0);
        if (!std::isfinite(linearMass) || !std::isfinite(nonlinearMass))
            continue;
        linear[index]->Fill(linearMass);
        nonlinear[index]->Fill(nonlinearMass);
        ++selected[index];
    }

    PeakFit linearFit[nSpecies];
    PeakFit nonlinearFit[nSpecies];
    for (int index = 0; index < nSpecies; ++index) {
        if (linear[index]->Integral() > 0)
            linear[index]->Scale(1.0 / linear[index]->Integral(), "width");
        if (nonlinear[index]->Integral() > 0)
            nonlinear[index]->Scale(1.0 / nonlinear[index]->Integral(), "width");
        styleHistogram(linear[index].get(), linearColor, 1);
        styleHistogram(nonlinear[index].get(), nonlinearColor, 2);
        linearFit[index] = fitPeak(linear[index].get(),
                                   Form("fitLinearZ%d", species[index].charge),
                                   species[index].mass);
        nonlinearFit[index] = fitPeak(
            nonlinear[index].get(),
            Form("fitNonlinearZ%d", species[index].charge), species[index].mass);
    }

    TCanvas canvas("massCanvas", "", 2100, 720);
    canvas.Divide(nSpecies, 1, 0.001, 0.001);
    for (int index = 0; index < nSpecies; ++index) {
        canvas.cd(index + 1);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.04);
        gPad->SetBottomMargin(0.14);
        gPad->SetTopMargin(0.09);
        gPad->SetGridy();

        const double maximum = 1.18 * std::max(linear[index]->GetMaximum(),
                                               nonlinear[index]->GetMaximum());
        linear[index]->SetMaximum(maximum);
        linear[index]->SetMinimum(0);
        linear[index]->SetTitle(Form("%s (Z=%d), %.1f < R < %.1f GV",
                                     species[index].label, species[index].charge,
                                     rigidityMin, rigidityMax));
        linear[index]->Draw("HIST");
        nonlinear[index]->Draw("HIST SAME");

        TLine *reference = new TLine(species[index].mass, 0,
                                     species[index].mass, maximum);
        reference->SetLineColor(referenceColor);
        reference->SetLineWidth(3);
        reference->SetLineStyle(7);
        reference->Draw("SAME");

        TLegend *legend = new TLegend(0.18, 0.72, 0.49, 0.89);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.027);
        legend->AddEntry(linear[index].get(), "Linear #beta", "l");
        legend->AddEntry(nonlinear[index].get(), "Nonlinear #beta", "l");
        legend->AddEntry(reference, "Known mass", "l");
        legend->Draw();

        TPaveText *metrics = new TPaveText(0.53, 0.69, 0.96, 0.90, "NDC");
        metrics->SetBorderSize(0);
        metrics->SetFillStyle(0);
        metrics->SetTextAlign(12);
        metrics->SetTextSize(0.027);
        metrics->AddText(Form("%s, %.1f<R<%.1f GV, N=%lld",
                              species[index].label, rigidityMin, rigidityMax,
                              selected[index]));
        if (linearFit[index].valid) {
            const double bias = 100.0 *
                (linearFit[index].mean / species[index].mass - 1.0);
            metrics->AddText(Form("L: #mu=%.4g, bias=%+.2f%%, #sigma/#mu=%.2f%%",
                                  linearFit[index].mean, bias,
                                  100.0 * linearFit[index].sigma /
                                      linearFit[index].mean));
        } else {
            metrics->AddText("Linear: peak fit failed");
        }
        if (nonlinearFit[index].valid) {
            const double bias = 100.0 *
                (nonlinearFit[index].mean / species[index].mass - 1.0);
            metrics->AddText(Form("NL: #mu=%.4g, bias=%+.2f%%, #sigma/#mu=%.2f%%",
                                  nonlinearFit[index].mean, bias,
                                  100.0 * nonlinearFit[index].sigma /
                                      nonlinearFit[index].mean));
        } else {
            metrics->AddText("Nonlinear: peak fit failed");
        }
        metrics->Draw();
        gPad->RedrawAxis();
    }

    const TString outputBase = TString(outputDirectory) + "/rigidity_mass_peaks";
    canvas.SaveAs(outputBase + ".png");
    const TString convertCommand = Form(
        "convert '%s.png' '%s.pdf' >/dev/null 2>&1",
        outputBase.Data(), outputBase.Data());
    if (gSystem->Exec(convertCommand) != 0)
        canvas.SaveAs(outputBase + ".pdf");

    TFile histogramFile(outputBase + ".root", "RECREATE");
    for (int index = 0; index < nSpecies; ++index) {
        linear[index]->Write();
        nonlinear[index]->Write();
    }
    histogramFile.Close();

    std::ofstream summary((outputBase + ".tsv").Data());
    summary << "species\tZ\tknown_mass_GeV\tmethod\tselected_events"
               "\tfit_valid\tmu_GeV\tmu_error_GeV\tsigma_GeV"
               "\tsigma_error_GeV\tbias_percent\tresolution_percent\n";
    summary << std::setprecision(9);
    for (int index = 0; index < nSpecies; ++index) {
        const PeakFit fits[] = {linearFit[index], nonlinearFit[index]};
        const char *methods[] = {"linear", "nonlinear"};
        for (int method = 0; method < 2; ++method) {
            const PeakFit &fit = fits[method];
            summary << species[index].label << '\t' << species[index].charge
                    << '\t' << species[index].mass << '\t' << methods[method]
                    << '\t' << selected[index] << '\t' << fit.valid << '\t'
                    << fit.mean << '\t' << fit.meanError << '\t' << fit.sigma
                    << '\t' << fit.sigmaError << '\t'
                    << (fit.valid ? 100.0 * (fit.mean / species[index].mass - 1.0)
                                  : 0.0)
                    << '\t'
                    << (fit.valid ? 100.0 * fit.sigma / fit.mean : 0.0) << '\n';
        }
    }

    std::cout << "Wrote " << outputBase << ".{png,pdf,root,tsv}" << std::endl;
}
