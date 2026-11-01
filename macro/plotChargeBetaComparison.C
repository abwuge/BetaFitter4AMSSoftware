#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "beta_fitter_macro_utils.h"

namespace
{
constexpr int kNBetaBins = 28;
constexpr int kNMethods = 3;
constexpr double kResidualMin = -0.5;
constexpr double kResidualMax = 1.5;
constexpr double kIqrToSigma = 0.741301109252801;

struct FitSeries
{
    std::vector<double> beta;
    std::vector<double> mean;
    std::vector<double> meanError;
    std::vector<double> sigma;
    std::vector<double> sigmaError;
};

struct ChargeResult
{
    int charge = 0;
    double betaMin = 0;
    double betaMax = 1;
    std::array<FitSeries, kNMethods> methods;
    std::array<TH1D *, kNMethods> lowBetaResiduals = {};
    Long64_t commonEntries = 0;
    Long64_t lowBetaEntries = 0;
};

bool fitCore(TH1D *hist, double &mean, double &meanError,
             double &sigma, double &sigmaError, int minimumEntries,
             double *amplitude = nullptr)
{
    if (!hist || hist->GetEntries() < minimumEntries)
        return false;

    const double probabilities[3] = {0.25, 0.50, 0.75};
    double quantiles[3] = {};
    hist->GetQuantiles(3, quantiles, probabilities);

    const double robustSigma = kIqrToSigma * (quantiles[2] - quantiles[0]);
    if (!(robustSigma > 0) || !std::isfinite(robustSigma))
        return false;

    const double fitMin = std::max(kResidualMin, quantiles[1] - 3 * robustSigma);
    const double fitMax = std::min(kResidualMax, quantiles[1] + 3 * robustSigma);
    TF1 gaussian(Form("gaussian_%s", hist->GetName()), "gaus", fitMin, fitMax);
    gaussian.SetParameters(hist->GetMaximum(), quantiles[1], robustSigma);
    gaussian.SetLineWidth(3);
    gaussian.SetLineStyle(2);

    const TFitResultPtr result = hist->Fit(&gaussian, "SQNR");
    if (!result.Get() || result->Status() != 0)
        return false;

    mean = result->Parameter(1);
    meanError = result->ParError(1);
    sigma = std::abs(result->Parameter(2));
    sigmaError = result->ParError(2);
    if(amplitude) *amplitude=result->Parameter(0);
    return std::isfinite(mean) && std::isfinite(meanError) &&
           std::isfinite(sigma) && std::isfinite(sigmaError) &&
           sigma < 0.08 && meanError < 0.02 &&
           sigmaError < 0.5 * sigma;
}

ChargeResult analyzeCharge(const char *fileName, int charge, int minimumEntries)
{
    ChargeResult output;
    output.charge = charge;

    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: cannot open " << fileName << std::endl;
        return output;
    }

    TTree *tree = static_cast<TTree *>(file->Get("betaTree"));
    if (!tree || tree->GetEntries() == 0)
    {
        std::cerr << "Error: missing or empty betaTree in " << fileName << std::endl;
        return output;
    }

    output.betaMin = std::max(0.25, tree->GetMinimum("mcBeta"));
    output.betaMax = std::min(1.0, tree->GetMaximum("mcBeta"));
    const double betaWidth = (output.betaMax - output.betaMin) / kNBetaBins;

    std::array<std::array<TH1D *, kNBetaBins>, kNMethods> histograms;
    for (int method = 0; method < kNMethods; ++method)
        for (int bin = 0; bin < kNBetaBins; ++bin)
            histograms[method][bin] = new TH1D(
                Form("hResidual_Z%d_M%d_B%d", charge, method, bin), "",
                800, kResidualMin, kResidualMax);
    for (int method=0;method<kNMethods;++method)
    {
        output.lowBetaResiduals[method]=new TH1D(Form("lowBeta_Z%d_M%d",charge,method),"",320,-0.16,0.16);
        output.lowBetaResiduals[method]->SetDirectory(nullptr);
    }

    double mcBeta = 0;
    double linearBeta = 0;
    double nonlinearBeta = 0;
    double hybridBeta = 0;
    int measuredCharge = charge, mcCharge=charge;
    const bool hasHybrid=tree->GetBranch("hybridBeta") &&
        tree->GetEntries("hybridBeta>0 && TMath::Finite(hybridBeta)")>0;
    tree->SetBranchAddress("mcBeta", &mcBeta);
    tree->SetBranchAddress("linearBeta", &linearBeta);
    tree->SetBranchAddress("nonlinearBeta", &nonlinearBeta);
    if(hasHybrid) tree->SetBranchAddress("hybridBeta", &hybridBeta);
    if(tree->GetBranch("Z")) tree->SetBranchAddress("Z",&measuredCharge);
    if(tree->GetBranch("mcCharge")) tree->SetBranchAddress("mcCharge",&mcCharge);

    const Long64_t entries = tree->GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry)
    {
        tree->GetEntry(entry);
        if(measuredCharge!=charge || mcCharge!=charge) continue;
        if (!(mcBeta > 0) || !std::isfinite(mcBeta) ||
            mcBeta < output.betaMin || mcBeta > output.betaMax)
            continue;

        int bin = static_cast<int>((mcBeta - output.betaMin) / betaWidth);
        bin = std::max(0, std::min(kNBetaBins - 1, bin));
        const double reconstructed[kNMethods] = {nonlinearBeta, linearBeta, hybridBeta};
        const int activeMethods=hasHybrid ? 3:2;
        bool common=true;
        for(int method=0;method<activeMethods;++method)
            common=common && reconstructed[method]>0 && std::isfinite(reconstructed[method]);
        if(!common) continue;
        ++output.commonEntries;
        if(mcBeta<0.9) ++output.lowBetaEntries;
        for (int method = 0; method < activeMethods; ++method)
        {
            if (!(reconstructed[method] > 0) || !std::isfinite(reconstructed[method]))
                continue;
            const double residual = 1.0 / reconstructed[method] - 1.0 / mcBeta;
            if (std::isfinite(residual))
                histograms[method][bin]->Fill(residual);
            if(mcBeta<0.9) output.lowBetaResiduals[method]->Fill(residual);
        }
    }

    for (int method = 0; method < kNMethods; ++method)
    {
        FitSeries &series = output.methods[method];
        for (int bin = 0; bin < kNBetaBins; ++bin)
        {
            double mean = 0;
            double meanError = 0;
            double sigma = 0;
            double sigmaError = 0;
            if (!fitCore(histograms[method][bin], mean, meanError, sigma, sigmaError,
                         minimumEntries))
                continue;

            series.beta.push_back(output.betaMin + (bin + 0.5) * betaWidth);
            series.mean.push_back(mean);
            series.meanError.push_back(meanError);
            series.sigma.push_back(sigma);
            series.sigmaError.push_back(sigmaError);
        }
    }

    std::cout << "Z=" << charge << ": " << entries << " entries, "
              << output.methods[0].beta.size() << " nonlinear bins, "
              << output.methods[1].beta.size() << " linear bins, "
              << output.methods[2].beta.size() << " hybrid bins, common N="
              << output.commonEntries << std::endl;
    file->Close();
    return output;
}

TGraphErrors *makeGraph(const FitSeries &series, bool drawSigma,
                        int color, int markerStyle)
{
    const std::vector<double> &values = drawSigma ? series.sigma : series.mean;
    const std::vector<double> &errors = drawSigma ? series.sigmaError : series.meanError;
    TGraphErrors *graph = new TGraphErrors(
        series.beta.size(), series.beta.data(), values.data(), nullptr, errors.data());
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(0.85);
    graph->SetLineWidth(2);
    return graph;
}

std::pair<double, double> graphRange(const ChargeResult &result, bool sigma)
{
    double minimum = sigma ? 0.0 : std::numeric_limits<double>::max();
    double maximum = sigma ? 0.0 : -std::numeric_limits<double>::max();
    for (int method = 0; method < kNMethods; ++method)
    {
        const FitSeries &series = result.methods[method];
        const std::vector<double> &values = sigma ? series.sigma : series.mean;
        const std::vector<double> &errors = sigma ? series.sigmaError : series.meanError;
        for (size_t i = 0; i < values.size(); ++i)
        {
            minimum = std::min(minimum, values[i] - errors[i]);
            maximum = std::max(maximum, values[i] + errors[i]);
        }
    }

    if (sigma)
        return std::make_pair(0.0, std::max(0.01, maximum * 1.18));

    minimum = std::min(minimum, 0.0);
    maximum = std::max(maximum, 0.0);
    const double padding = std::max(0.005, 0.12 * (maximum - minimum));
    return std::make_pair(minimum - padding, maximum + padding);
}
} // namespace

void plotChargeBetaComparison(
    const char *z2File = "results/11.root",
    const char *z6File = "results/12.root",
    const char *z8File = "results/13.root",
    const char *outputName = nullptr,
    int minimumEntries = 300,
    const char *z1File = nullptr,
    bool reportPages = false)
{
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);
    std::vector<const char *> files = {z2File, z6File, z8File};
    std::vector<int> charges = {2, 6, 8};
    if(z1File && z1File[0])
    { files.insert(files.begin(),z1File); charges.insert(charges.begin(),1); }
    const int nCharges=charges.size();
    std::vector<ChargeResult> results(nCharges);
    for (int index = 0; index < nCharges; ++index)
        results[index] = analyzeCharge(files[index], charges[index], minimumEntries);

    TString resolvedOutput = BetaFitterMacro::output_path(
        outputName, "beta_comparison", "charge_beta_center_comparison.pdf");
    TString pngOutput = resolvedOutput;
    if (pngOutput.EndsWith(".pdf"))
        pngOutput.ReplaceAll(".pdf", ".png");
    else
        pngOutput += ".png";

    const int nonlinearColor = kP6Blue;
    const int linearColor = TColor::GetColor("#222222");
    const int hybridColor = kP6Red;
    TCanvas canvas("chargeBetaComparison", "", 500*nCharges, 860);
    canvas.Divide(nCharges, 2, 0.002, 0.002);

    for (int chargeIndex = 0; chargeIndex < nCharges; ++chargeIndex)
    {
        const ChargeResult &result = results[chargeIndex];
        for (int row = 0; row < 2; ++row)
        {
            canvas.cd(row * nCharges + chargeIndex + 1);
            gPad->SetLeftMargin(0.18);
            gPad->SetRightMargin(0.035);
            gPad->SetBottomMargin(row == 1 ? 0.15 : 0.11);
            gPad->SetTopMargin(row == 0 ? 0.12 : 0.06);
            gPad->SetGridx();
            gPad->SetGridy();

            const bool sigma = row == 1;
            const std::pair<double, double> range = graphRange(result, sigma);
            TH1F *frame = gPad->DrawFrame(result.betaMin, range.first,
                                          result.betaMax, range.second);
            frame->SetTitle(row == 0 ? Form("Z = %d", result.charge) : "");
            frame->GetXaxis()->SetTitle("#beta_{MC} at L5 (center ref.)");
            frame->GetYaxis()->SetTitle(sigma ? "Core #sigma(1/#beta residual)"
                                               : "Core #mu(1/#beta residual)");
            frame->GetXaxis()->SetTitleSize(0.052);
            frame->GetYaxis()->SetTitleSize(0.052);
            frame->GetXaxis()->SetLabelSize(0.043);
            frame->GetYaxis()->SetLabelSize(0.043);
            frame->GetYaxis()->SetTitleOffset(1.55);

            TLine zero(result.betaMin, 0, result.betaMax, 0);
            if (!sigma)
            {
                zero.SetLineColor(TColor::GetColor("#788AA3"));
                zero.SetLineStyle(2);
                zero.SetLineWidth(2);
                zero.DrawClone("SAME");
            }

            TGraphErrors *nonlinear = makeGraph(result.methods[0], sigma,
                                                 nonlinearColor, 20);
            TGraphErrors *linear = makeGraph(result.methods[1], sigma,
                                              linearColor, 21);
            nonlinear->Draw("LP SAME");
            linear->Draw("LP SAME");
            TGraphErrors *hybrid=makeGraph(result.methods[2],sigma,hybridColor,22);
            if(hybrid->GetN()) hybrid->Draw("LP SAME");

            if (row == 0)
            {
                TLatex chargeLabel;
                chargeLabel.SetNDC();
                chargeLabel.SetTextFont(42);
                chargeLabel.SetTextSize(0.060);
                chargeLabel.SetTextAlign(22);
                chargeLabel.DrawLatex(0.52, 0.94, Form("Z = %d, N = %lld", result.charge,result.commonEntries));
            }

            if (row == 0 && chargeIndex == 0)
            {
                TLegend *legend = new TLegend(0.58, 0.67, 0.94, 0.84);
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                legend->SetTextSize(0.039);
                legend->AddEntry(nonlinear, "Nonlinear", "lp");
                legend->AddEntry(linear, "Linear", "lp");
                if(hybrid->GetN()) legend->AddEntry(hybrid,"Hybrid (1 cm)","lp");
                legend->Draw();
            }
        }
    }

    canvas.SaveAs(resolvedOutput);
    canvas.SaveAs(pngOutput);
    // Re-layout the exact fitted panels for presentation; keep the event
    // selection, fit results and original overview untouched.
    if(reportPages && nCharges==4)
        for(int row=0;row<2;++row)
        {
            TCanvas page(Form("reportCore%d",row),"",1600,1000);
            for(int index=0;index<nCharges;++index)
            {
                page.cd();
                auto pad=static_cast<TPad *>(canvas.GetPad(row*nCharges+index+1)->Clone(
                    Form("reportCorePad%d_%d",row,index)));
                const double x=(index%2)*0.5, y=(1-index/2)*0.5;
                pad->SetPad(x,y,x+0.5,y+0.5);
                pad->SetTopMargin(0.12);
                pad->SetBottomMargin(0.17);
                pad->Draw();pad->cd();
                if(row==1)
                {
                    TLatex label;label.SetNDC();label.SetTextFont(42);
                    label.SetTextSize(0.055);label.SetTextAlign(22);
                    label.DrawLatex(0.52,0.94,Form("Z = %d, N = %lld",results[index].charge,results[index].commonEntries));
                    if(index==0)
                    {
                        auto legend=new TLegend(0.57,0.68,0.93,0.85);
                        legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(0.037);
                        const int legendColors[3]={nonlinearColor,linearColor,hybridColor};
                        const char *legendLabels[3]={"Nonlinear","Linear","Hybrid (1 cm)"};
                        for(int method=0;method<3;++method)
                            legend->AddEntry(makeGraph(results[index].methods[method],true,legendColors[method],20+method),legendLabels[method],"lp");
                        legend->Draw();
                    }
                }
            }
            TString pdf=resolvedOutput;pdf.ReplaceAll(".pdf",row==0 ? "_mean.pdf":"_sigma.pdf");
            TString png=pdf;png.ReplaceAll(".pdf",".png");
            page.SaveAs(pdf);page.SaveAs(png);
        }
    TString rootOutput=resolvedOutput; rootOutput.ReplaceAll(".pdf",".root");
    TFile numerical(rootOutput,"RECREATE");
    for(const auto &result:results)
        for(int method=0;method<kNMethods;++method)
        {
            result.lowBetaResiduals[method]->Write();
            auto mean=makeGraph(result.methods[method],false,kBlack,20);
            auto sigma=makeGraph(result.methods[method],true,kBlack,20);
            mean->Write(Form("mean_Z%d_M%d",result.charge,method));
            sigma->Write(Form("sigma_Z%d_M%d",result.charge,method));
        }
    numerical.Close();
    TCanvas residualCanvas("hybridResiduals","",500*nCharges,580);
    residualCanvas.Divide(nCharges,1,0.005,0.005);
    const int colors[3]={nonlinearColor,linearColor,hybridColor};
    const char *labels[3]={"Nonlinear","Linear","Hybrid"};
    for(int index=0;index<nCharges;++index)
    {
        residualCanvas.cd(index+1);
        gPad->SetLeftMargin(0.18);gPad->SetRightMargin(0.035);
        gPad->SetBottomMargin(0.14);gPad->SetTopMargin(0.08);
        auto &result=results[index];
        double maximum=0;
        for(auto hist:result.lowBetaResiduals) maximum=std::max(maximum,hist->GetMaximum());
        auto frame=gPad->DrawFrame(-0.08,0,0.08,maximum*1.8);
        frame->GetXaxis()->SetTitle("1/#beta_{reco} - 1/#beta_{MC}");
        frame->GetYaxis()->SetTitle("Events / 0.001");
        frame->GetXaxis()->SetTitleSize(0.047);frame->GetYaxis()->SetTitleSize(0.047);
        frame->GetXaxis()->SetLabelSize(0.040);frame->GetYaxis()->SetLabelSize(0.040);
        frame->GetYaxis()->SetTitleOffset(1.65);
        TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(0.037);
        text.DrawLatex(0.20,0.96,Form("Z = %d, #beta_{MC}<0.9, N = %lld",result.charge,result.lowBetaEntries));
        for(int method=0;method<kNMethods;++method)
        {
            auto hist=result.lowBetaResiduals[method];
            if(!hist->GetEntries()) continue;
            hist->SetLineColor(colors[method]);hist->SetLineWidth(2);hist->Draw("HIST SAME");
            double mean,em,sigma,es,amplitude;
            if(fitCore(hist,mean,em,sigma,es,50,&amplitude))
            {
                TF1 fit(Form("lowfit_Z%d_M%d",result.charge,method),"gaus",mean-3*sigma,mean+3*sigma);
                fit.SetParameters(amplitude,mean,sigma);
                fit.SetLineColor(colors[method]);fit.SetLineStyle(2);fit.SetLineWidth(3);fit.DrawCopy("SAME");
                text.SetTextColor(colors[method]);
                text.DrawLatex(0.23,0.87-method*0.105,labels[method]);
                text.DrawLatex(0.23,0.825-method*0.105,Form("#mu=%.4f#pm%.4f, #sigma=%.4f#pm%.4f",mean,em,sigma,es));
                std::cout<<"LOW_BETA Z="<<result.charge<<" "<<labels[method]<<" N="<<result.lowBetaEntries<<" mean="<<mean<<" meanError="<<em<<" sigma="<<sigma<<" sigmaError="<<es<<std::endl;
            }
        }
    }
    TString residualPdf=resolvedOutput;residualPdf.ReplaceAll(".pdf","_residuals.pdf");
    TString residualPng=residualPdf;residualPng.ReplaceAll(".pdf",".png");
    residualCanvas.SaveAs(residualPdf);residualCanvas.SaveAs(residualPng);
    if(reportPages && nCharges==4)
    {
        TCanvas page("reportLowBeta","",1600,1000);
        for(int index=0;index<nCharges;++index)
        {
            page.cd();
            auto pad=static_cast<TPad *>(residualCanvas.GetPad(index+1)->Clone(Form("reportResidualPad%d",index)));
            const double x=(index%2)*0.5, y=(1-index/2)*0.5;
            pad->SetPad(x,y,x+0.5,y+0.5);pad->Draw();
        }
        TString pdf=resolvedOutput;pdf.ReplaceAll(".pdf","_low_beta.pdf");
        TString png=pdf;png.ReplaceAll(".pdf",".png");
        page.SaveAs(pdf);page.SaveAs(png);
    }
    std::cout << "Saved comparison to " << resolvedOutput << " and " << pngOutput << std::endl;
}

void plotTrackerEnergyLossComparison(
    const char *oldZ6File = "results/12.root",
    const char *oldZ8File = "results/13.root",
    const char *trackerZ6File = "results/tk_edep_20260731/beta_C12_tracker.root",
    const char *trackerZ8File = "results/tk_edep_20260731/beta_O16_tracker.root",
    const char *outputName = nullptr,
    int trackerMinimumEntries = 20)
{
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    const std::array<const char *, 2> oldFiles = {oldZ6File, oldZ8File};
    const std::array<const char *, 2> trackerFiles = {trackerZ6File, trackerZ8File};
    const std::array<int, 2> charges = {6, 8};
    std::array<ChargeResult, 2> results;
    for (int index = 0; index < 2; ++index)
    {
        const ChargeResult oldResult = analyzeCharge(oldFiles[index], charges[index], 300);
        const ChargeResult trackerResult =
            analyzeCharge(trackerFiles[index], charges[index], trackerMinimumEntries);
        results[index] = trackerResult;
        results[index].betaMin = std::min(oldResult.betaMin, trackerResult.betaMin);
        results[index].betaMax = std::max(oldResult.betaMax, trackerResult.betaMax);
        results[index].methods[0] = trackerResult.methods[0];
        results[index].methods[1] = oldResult.methods[0];
    }

    TString resolvedOutput = BetaFitterMacro::output_path(
        outputName, "beta_comparison", "tracker_energy_loss_beta_comparison.pdf");
    TString pngOutput = resolvedOutput;
    if (pngOutput.EndsWith(".pdf"))
        pngOutput.ReplaceAll(".pdf", ".png");
    else
        pngOutput += ".png";

    const int trackerColor = TColor::GetColor("#45364B");
    const int oldColor = TColor::GetColor("#C45B46");
    TCanvas canvas("trackerEnergyLossComparison", "", 1100, 860);
    canvas.Divide(2, 2, 0.002, 0.002);
    for (int chargeIndex = 0; chargeIndex < 2; ++chargeIndex)
    {
        const ChargeResult &result = results[chargeIndex];
        for (int row = 0; row < 2; ++row)
        {
            canvas.cd(row * 2 + chargeIndex + 1);
            gPad->SetLeftMargin(chargeIndex == 0 ? 0.16 : 0.12);
            gPad->SetRightMargin(0.035);
            gPad->SetBottomMargin(row == 1 ? 0.15 : 0.11);
            gPad->SetTopMargin(row == 0 ? 0.12 : 0.06);
            gPad->SetGridx();
            gPad->SetGridy();

            const bool sigma = row == 1;
            const std::pair<double, double> range = graphRange(result, sigma);
            TH1F *frame = gPad->DrawFrame(result.betaMin, range.first,
                                          result.betaMax, range.second);
            frame->GetXaxis()->SetTitle("#beta_{MC} at AMS center");
            frame->GetYaxis()->SetTitle(sigma ? "Core #sigma(1/#beta residual)"
                                               : "Core #mu(1/#beta residual)");
            frame->GetXaxis()->SetTitleSize(0.052);
            frame->GetYaxis()->SetTitleSize(0.052);
            frame->GetXaxis()->SetLabelSize(0.043);
            frame->GetYaxis()->SetLabelSize(0.043);
            frame->GetYaxis()->SetTitleOffset(chargeIndex == 0 ? 1.45 : 1.20);

            if (!sigma)
            {
                TLine zero(result.betaMin, 0, result.betaMax, 0);
                zero.SetLineColor(TColor::GetColor("#788AA3"));
                zero.SetLineStyle(2);
                zero.SetLineWidth(2);
                zero.Draw("SAME");
            }

            TGraphErrors *tracker = makeGraph(result.methods[0], sigma, trackerColor, 20);
            TGraphErrors *old = makeGraph(result.methods[1], sigma, oldColor, 21);
            tracker->Draw("LP SAME");
            old->Draw("LP SAME");

            if (row == 0)
            {
                TLatex chargeLabel;
                chargeLabel.SetNDC();
                chargeLabel.SetTextFont(42);
                chargeLabel.SetTextSize(0.060);
                chargeLabel.SetTextAlign(22);
                chargeLabel.DrawLatex(0.52, 0.94, Form("Z = %d", result.charge));
            }
            if (row == 0 && chargeIndex == 0)
            {
                TLegend *legend = new TLegend(0.52, 0.67, 0.92, 0.84);
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                legend->SetTextSize(0.040);
                legend->AddEntry(tracker, "Tracker-layer", "lp");
                legend->AddEntry(old, "Old constant-#beta", "lp");
                legend->Draw();
            }
        }
    }
    canvas.SaveAs(resolvedOutput);
    canvas.SaveAs(pngOutput);
    std::cout << "Saved tracker comparison to " << resolvedOutput << " and " << pngOutput << std::endl;
}
