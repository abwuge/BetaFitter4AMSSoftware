#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
constexpr int kNLayers = 4;
constexpr int kNCharges = 4;
constexpr int kCharges[kNCharges] = {1, 2, 6, 8};
constexpr double kC = 29.9792458; // cm/ns
constexpr double kIqrToSigma = 0.741301109252801;
constexpr int kColors[kNCharges] = {kP6Blue, kP6Yellow, kP6Red, kP6Violet};
constexpr int kMarkers[kNCharges] = {20, 21, 22, 23};

struct LayerObjects
{
    TH1D *residual = nullptr;
    TProfile *vsBeta = nullptr;
    TProfile *vsEdep = nullptr;
    TProfile *vsPosition = nullptr;
    Long64_t entries = 0;
    double mean = 0;
    double meanError = 0;
    double sigma = 0;
    double sigmaError = 0;
};

bool addValidationFiles(TChain &chain, const char *sourcePath, int charge)
{
    std::ifstream input(sourcePath);
    if (!input)
    {
        std::cerr << "Cannot open source list: " << sourcePath << std::endl;
        return false;
    }

    std::string line;
    int added = 0;
    while (std::getline(input, line))
    {
        const std::string::size_type first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#')
            continue;
        std::istringstream values(line);
        std::string role, path;
        int z = 0;
        if (!(values >> role >> z >> path))
            continue;
        if (role == "validation" && z == charge)
        {
            chain.Add(path.c_str());
            ++added;
        }
    }
    std::cout << "Z=" << charge << ": added " << added << " validation files" << std::endl;
    return added > 0;
}

bool fitCore(TH1D *hist, double &mean, double &meanError,
             double &sigma, double &sigmaError)
{
    if (!hist || hist->GetEntries() < 100)
        return false;

    const double probabilities[3] = {0.25, 0.50, 0.75};
    double q[3] = {};
    hist->GetQuantiles(3, q, probabilities);
    const double robustSigma = kIqrToSigma * (q[2] - q[0]);
    if (!(robustSigma > 0) || !std::isfinite(robustSigma))
        return false;

    const double lo = std::max(hist->GetXaxis()->GetXmin(), q[1] - 3 * robustSigma);
    const double hi = std::min(hist->GetXaxis()->GetXmax(), q[1] + 3 * robustSigma);
    TF1 gaussian(Form("fit_%s", hist->GetName()), "gaus", lo, hi);
    gaussian.SetParameters(hist->GetMaximum(), q[1], robustSigma);
    const TFitResultPtr result = hist->Fit(&gaussian, "SQNR");
    if (!result.Get() || result->Status() != 0)
        return false;

    mean = gaussian.GetParameter(1);
    meanError = gaussian.GetParError(1);
    sigma = std::abs(gaussian.GetParameter(2));
    sigmaError = gaussian.GetParError(2);
    return std::isfinite(mean) && std::isfinite(sigma) && sigma > 0;
}

void styleProfile(TProfile *profile, int index)
{
    profile->SetLineColor(kColors[index]);
    profile->SetMarkerColor(kColors[index]);
    profile->SetMarkerStyle(kMarkers[index]);
    profile->SetMarkerSize(0.85);
    profile->SetLineWidth(2);
}

void drawZeroLine(double xmin, double xmax)
{
    TLine *line = new TLine(xmin, 0, xmax, 0);
    line->SetLineColor(kGray + 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();
}

void drawPageHeader(int layer, const char *subtitle)
{
    TLatex text;
    text.SetNDC();
    text.SetTextFont(42);
    text.SetTextSize(0.040);
    text.DrawLatex(0.12, 0.94, Form("TOF S%d timing residual", layer + 1));
    text.SetTextSize(0.031);
    text.DrawLatex(0.12, 0.895, subtitle);
}

void printMultipage(TCanvas &canvas, const TString &pdf, int page, int pages)
{
    if (page == 0)
        canvas.Print(pdf + "(");
    else if (page + 1 == pages)
        canvas.Print(pdf + ")");
    else
        canvas.Print(pdf);
}

} // namespace

void plotLayerResiduals(const char *sourcePath = "config/sources.conf",
                        const char *outputDir = "results/layer_residuals",
                        double betaMin = 0.6)
{
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gSystem->mkdir(outputDir, true);

    std::array<std::array<LayerObjects, kNLayers>, kNCharges> objects;

    for (int iz = 0; iz < kNCharges; ++iz)
    {
        const int targetCharge = kCharges[iz];
        TChain chain("amstreea");
        if (!addValidationFiles(chain, sourcePath, targetCharge))
            continue;

        float mch = 0;
        float mmass = 0;
        float mevmom1[21] = {};
        float tofTl[4] = {};
        float tofLength[4] = {};
        float tofPos[4][3] = {};
        float tofEdep[4] = {};
        int tofBpatt = 0;
        float tkQ[2] = {};
        float tkQin[2][3] = {};
        float tkEdep[9] = {};
        float tkPos[9][3] = {};

        chain.SetBranchStatus("*", 0);
        const char *branches[] = {"mch", "mmass", "mevmom1", "tof_tl", "tof_leng",
                                  "tof_pos", "tof_edep", "tof_bpatt", "tk_q", "tk_qin",
                                  "tk_edep", "tk_pos"};
        for (const char *branch : branches)
            chain.SetBranchStatus(branch, 1);
        chain.SetBranchAddress("mch", &mch);
        chain.SetBranchAddress("mmass", &mmass);
        chain.SetBranchAddress("mevmom1", mevmom1);
        chain.SetBranchAddress("tof_tl", tofTl);
        chain.SetBranchAddress("tof_leng", tofLength);
        chain.SetBranchAddress("tof_pos", tofPos);
        chain.SetBranchAddress("tof_edep", tofEdep);
        chain.SetBranchAddress("tof_bpatt", &tofBpatt);
        chain.SetBranchAddress("tk_q", tkQ);
        chain.SetBranchAddress("tk_qin", tkQin);
        chain.SetBranchAddress("tk_edep", tkEdep);
        chain.SetBranchAddress("tk_pos", tkPos);

        for (int layer = 0; layer < kNLayers; ++layer)
        {
            LayerObjects &obj = objects[iz][layer];
            obj.residual = new TH1D(Form("residual_Z%d_S%d", targetCharge, layer + 1), "",
                                    280, -0.35, 0.35);
            obj.vsBeta = new TProfile(Form("beta_Z%d_S%d", targetCharge, layer + 1), "",
                                      24, betaMin, 1.0);
            obj.vsEdep = new TProfile(Form("edep_Z%d_S%d", targetCharge, layer + 1), "",
                                      30, -0.4, 1.4);
            obj.vsPosition = new TProfile(Form("position_Z%d_S%d", targetCharge, layer + 1), "",
                                          32, -130, 130);
            obj.residual->SetDirectory(nullptr);
            obj.vsBeta->SetDirectory(nullptr);
            obj.vsEdep->SetDirectory(nullptr);
            obj.vsPosition->SetDirectory(nullptr);
            styleProfile(obj.vsBeta, iz);
            styleProfile(obj.vsEdep, iz);
            styleProfile(obj.vsPosition, iz);
        }

        const Long64_t entries = chain.GetEntries();
        for (Long64_t entry = 0; entry < entries; ++entry)
        {
            chain.GetEntry(entry);
            if (tofBpatt != 4444)
                continue;

            const int reconstructedCharge = static_cast<int>((tkQin[0][2] < 2.5 ? tkQ[1] : tkQin[0][2]) + 0.5);
            const int mcCharge = static_cast<int>(mch + 0.5);
            if (reconstructedCharge != targetCharge || mcCharge != targetCharge)
                continue;

            // Match the raw-event quality gate used by the independent estimator study.
            bool trackerComplete = true;
            for (int trackerLayer = 1; trackerLayer <= 7; ++trackerLayer)
            {
                trackerComplete = trackerComplete && tkEdep[trackerLayer] > 0 &&
                                  std::isfinite(tkEdep[trackerLayer]) &&
                                  std::isfinite(tkPos[trackerLayer][0]) &&
                                  std::isfinite(tkPos[trackerLayer][1]) &&
                                  std::isfinite(tkPos[trackerLayer][2]);
            }
            if (!trackerComplete || !(tkPos[4][2] > 0) || !(tkPos[5][2] < 0))
                continue;
            if (!(mmass > 0) || mevmom1[10] <= 0 || !std::isfinite(mevmom1[10]))
                continue;
            const double momentum = mevmom1[10];
            const double mcBeta = momentum / std::sqrt(momentum * momentum + mmass * mmass);
            if (!(mcBeta > betaMin && mcBeta <= 1.0) || !std::isfinite(mcBeta))
                continue;

            double x[kNLayers] = {};
            double tMean = 0;
            double xMean = 0;
            bool valid = true;
            for (int layer = 0; layer < kNLayers; ++layer)
            {
                x[layer] = tofLength[layer] / kC;
                valid = valid && std::isfinite(x[layer]) && std::isfinite(tofTl[layer]) &&
                        std::isfinite(tofPos[layer][0]) && std::isfinite(tofPos[layer][1]) &&
                        std::isfinite(tofEdep[layer]) && tofEdep[layer] > 0;
                xMean += x[layer];
                tMean += tofTl[layer];
            }
            if (!valid)
                continue;
            xMean /= kNLayers;
            tMean /= kNLayers;

            double numerator = 0;
            double denominator = 0;
            for (int layer = 0; layer < kNLayers; ++layer)
            {
                const double dx = x[layer] - xMean;
                numerator += dx * (tofTl[layer] - tMean);
                denominator += dx * dx;
            }
            if (!(denominator > 0) || !std::isfinite(denominator))
                continue;
            const double inverseBeta = numerator / denominator;
            if (!(inverseBeta > 0) || !std::isfinite(inverseBeta))
                continue;

            for (int layer = 0; layer < kNLayers; ++layer)
            {
                const double residual = tofTl[layer] - tMean - inverseBeta * (x[layer] - xMean);
                const double normalizedDeposit = tofEdep[layer] / (targetCharge * targetCharge);
                const double logDeposit = std::log10(normalizedDeposit);
                // From the observed bar orientation: S1/S4 bars run along x; S2/S3 along y.
                const double alongBar = (layer == 0 || layer == 3) ? tofPos[layer][0] : tofPos[layer][1];
                if (!std::isfinite(residual) || !std::isfinite(logDeposit) || !std::isfinite(alongBar))
                    continue;

                LayerObjects &obj = objects[iz][layer];
                obj.residual->Fill(residual);
                obj.vsBeta->Fill(mcBeta, residual);
                obj.vsEdep->Fill(logDeposit, residual);
                obj.vsPosition->Fill(alongBar, residual);
                ++obj.entries;
            }
        }

        for (int layer = 0; layer < kNLayers; ++layer)
        {
            LayerObjects &obj = objects[iz][layer];
            fitCore(obj.residual, obj.mean, obj.meanError, obj.sigma, obj.sigmaError);
        }
    }

    const TString output(outputDir);
    TFile rootOutput(output + "/layer_residuals.root", "RECREATE");
    for (int iz = 0; iz < kNCharges; ++iz)
        for (int layer = 0; layer < kNLayers; ++layer)
        {
            objects[iz][layer].residual->Write();
            objects[iz][layer].vsBeta->Write();
            objects[iz][layer].vsEdep->Write();
            objects[iz][layer].vsPosition->Write();
        }
    rootOutput.Close();

    // 1. Residual distributions.
    {
        TCanvas canvas("distribution_canvas", "", 1280, 720);
        const TString pdf = output + "/residual_distribution.pdf";
        for (int layer = 0; layer < kNLayers; ++layer)
        {
            canvas.Clear();
            canvas.SetLeftMargin(0.11);
            canvas.SetRightMargin(0.04);
            canvas.SetBottomMargin(0.12);
            canvas.SetTopMargin(0.12);

            double maximum = 0;
            for (int iz = 0; iz < kNCharges; ++iz)
            {
                TH1D *hist = objects[iz][layer].residual;
                const double integral = hist->Integral();
                if (integral > 0)
                    hist->Scale(1.0 / integral);
                hist->SetLineColor(kColors[iz]);
                hist->SetLineWidth(3);
                hist->SetFillStyle(0);
                maximum = std::max(maximum, hist->GetMaximum());
            }

            TH1D *frame = static_cast<TH1D *>(objects[0][layer].residual->Clone(Form("frame_dist_S%d", layer + 1)));
            frame->Reset();
            frame->GetXaxis()->SetTitle("Orthogonal timing residual [ns]");
            frame->GetYaxis()->SetTitle("Normalized entries");
            frame->GetXaxis()->SetRangeUser(-0.24, 0.24);
            frame->SetMinimum(0);
            frame->SetMaximum(maximum * 1.28);
            frame->Draw("AXIS");

            TLegend legend(0.62, 0.62, 0.94, 0.86);
            legend.SetBorderSize(0);
            legend.SetFillStyle(0);
            legend.SetTextSize(0.030);
            Long64_t sharedN = 0;
            for (int iz = 0; iz < kNCharges; ++iz)
            {
                LayerObjects &obj = objects[iz][layer];
                obj.residual->Draw("HIST SAME");
                legend.AddEntry(obj.residual,
                                Form("Z=%d: #mu=%.1f ps, #sigma=%.1f ps",
                                     kCharges[iz], obj.mean * 1000, obj.sigma * 1000),
                                "l");
                sharedN += obj.entries;
            }
            legend.Draw();
            drawPageHeader(layer, Form("0.6 < #beta_{MC} #leq 1; unit-area shapes; total N = %lld", sharedN));
            canvas.RedrawAxis();
            printMultipage(canvas, pdf, layer, kNLayers);
            canvas.SaveAs(output + Form("/residual_distribution_S%d.png", layer + 1));
            delete frame;
        }
    }

    auto drawProfiles = [&](const char *stem, const char *xTitle,
                            TProfile *(LayerObjects::*member),
                            double xMin, double xMax, double yAbs,
                            const char *subtitle)
    {
        TCanvas canvas(Form("%s_canvas", stem), "", 1280, 720);
        const TString pdf = output + Form("/%s.pdf", stem);
        for (int layer = 0; layer < kNLayers; ++layer)
        {
            canvas.Clear();
            canvas.SetLeftMargin(0.11);
            canvas.SetRightMargin(0.04);
            canvas.SetBottomMargin(0.12);
            canvas.SetTopMargin(0.12);

            TH1D frame(Form("frame_%s_S%d", stem, layer + 1), "", 100, xMin, xMax);
            frame.GetXaxis()->SetTitle(xTitle);
            frame.GetYaxis()->SetTitle("Mean timing residual [ns]");
            frame.SetMinimum(-yAbs);
            frame.SetMaximum(yAbs);
            frame.Draw("AXIS");
            drawZeroLine(xMin, xMax);

            TLegend legend(0.79, 0.70, 0.94, 0.87);
            legend.SetBorderSize(0);
            legend.SetFillStyle(0);
            legend.SetTextSize(0.031);
            Long64_t totalN = 0;
            for (int iz = 0; iz < kNCharges; ++iz)
            {
                TProfile *profile = objects[iz][layer].*member;
                profile->Draw("E1 SAME");
                legend.AddEntry(profile, Form("Z=%d", kCharges[iz]), "lp");
                totalN += objects[iz][layer].entries;
            }
            legend.Draw();
            drawPageHeader(layer, Form("%s; total N = %lld", subtitle, totalN));
            canvas.RedrawAxis();
            printMultipage(canvas, pdf, layer, kNLayers);
            canvas.SaveAs(output + Form("/%s_S%d.png", stem, layer + 1));
        }
    };

    drawProfiles("residual_vs_beta", "#beta_{MC}", &LayerObjects::vsBeta,
                 betaMin, 1.0, 0.030, "0.6 < #beta_{MC} #leq 1; profile mean #pm statistical error");
    drawProfiles("residual_vs_edep", "log_{10}[(TOF deposit / Z^{2}) / MeV]", &LayerObjects::vsEdep,
                 -0.4, 1.4, 0.040, "Measured layer deposit; profile mean #pm statistical error");
    drawProfiles("residual_vs_position", "Hit position along bar [cm]", &LayerObjects::vsPosition,
                 -130, 130, 0.040, "S1/S4: x coordinate; S2/S3: y coordinate");

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "\n# Layer residual Gaussian core [ns]\n";
    std::cout << "# Z layer entries mean meanErr sigma sigmaErr\n";
    for (int iz = 0; iz < kNCharges; ++iz)
        for (int layer = 0; layer < kNLayers; ++layer)
        {
            const LayerObjects &obj = objects[iz][layer];
            std::cout << kCharges[iz] << " " << (layer + 1) << " " << obj.entries << " "
                      << obj.mean << " " << obj.meanError << " "
                      << obj.sigma << " " << obj.sigmaError << "\n";
        }
}
