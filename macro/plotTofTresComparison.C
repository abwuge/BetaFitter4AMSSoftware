#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TFitResult.h>
#include <array>

/**
 * Draw plots of tof_tres[4] time residual values for different TOF paddles
 * and segments within each paddle
 *
 * @param fileName Input ROOT file path
 * @param outputName Output PDF file name
 */
void plotTofTresComparison(
    std::string fileName = "test.root",
    const char *outputName = "test_tof_tres_comparison.pdf")
{
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetEndErrorSize(20);

    if (fileName.size() < 5 || fileName.substr(fileName.size() - 5) != ".root")
        fileName += ".root";

    TFile *file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        fileName = "results/" + fileName;
        std::cout << "Trying to open file: " << fileName << std::endl;
        file = TFile::Open(fileName.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "Error: Unable to find/open file " << fileName << std::endl;
            return;
        }
    }

    TTree *tree = (TTree *)file->Get("betaDiff");
    if (!tree)
    {
        std::cerr << "Error: Could not find betaDiff tree in " << fileName << std::endl;
        file->Close();
        delete file;
        return;
    }

    if (tree->GetEntries() == 0)
    {
        std::cerr << "Error: Empty tree, no data to plot" << std::endl;
        file->Close();
        delete file;
        return;
    }

    float tof_tres[4]{};    
    float mevcoo1[21][3]{}; 

    tree->SetBranchAddress("tof_tres", tof_tres);
    tree->SetBranchAddress("mevcoo1", mevcoo1);

    // Define TOF stations configuration
    const int nPaddles[4] = {8, 8, 10, 8};
    const int mevcooIndex[4] = {4, 6, 15, 16};
    const bool isXAxisParallel[4] = {true, false, false, true};
    const float paddleWidth = 12.0;
    const char* stationNames[4] = {"S1", "S2", "S3", "S4"};

    // Define segments within each paddle
    const int nSegments = 10;  // 10 segments of 12cm each
    const float paddleLength = 120.0;
    const float segmentLength = paddleLength / nSegments;

    const int nBinsX = 100;
    double minTres = -0.2;
    double maxTres = 0.2;

    // Create histograms for each paddle
    std::vector<std::vector<TH1F*>> hTofTres;
    
    // Create histograms for segments within each paddle
    std::vector<std::vector<std::vector<TH1F*>>> hTofTresSegments;
    
    // Initialize histograms
    for (int station = 0; station < 4; ++station) {
        std::vector<TH1F*> stationHists;
        std::vector<std::vector<TH1F*>> stationSegmentHists;
        
        for (int paddle = 0; paddle < nPaddles[station]; ++paddle) {
            // Main paddle histogram
            TH1F* hist = new TH1F(Form("hTofTres_%s_P%d", stationNames[station], paddle),
                                 Form("TOF Time Residual - %s Paddle %d;Time Residual [ns];Entries", 
                                      stationNames[station], paddle+1),
                                 nBinsX, minTres, maxTres);
            
            int color = paddle + 1;
            if (color >= 10) color += 10;
            hist->SetLineColor(color);
            hist->SetLineWidth(2);
            
            stationHists.push_back(hist);
            
            // Segment histograms
            std::vector<TH1F*> paddleSegmentHists;
            for (int segment = 0; segment < nSegments; ++segment) {
                TH1F* segHist = new TH1F(Form("hTofTres_%s_P%d_S%d", stationNames[station], paddle, segment),
                                       Form("TOF Time Residual - %s Paddle %d Segment %d;Time Residual [ns];Entries", 
                                            stationNames[station], paddle+1, segment+1),
                                       nBinsX, minTres, maxTres);
                paddleSegmentHists.push_back(segHist);
            }
            stationSegmentHists.push_back(paddleSegmentHists);
        }
        
        hTofTres.push_back(stationHists);
        hTofTresSegments.push_back(stationSegmentHists);
    }

    std::cout << "Processing entries: " << tree->GetEntries() << std::endl;

    for (int i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);

        for (int station = 0; station < 4; ++station)
        {
            if (mevcoo1[mevcooIndex[station]][0] == -1000 || 
                mevcoo1[mevcooIndex[station]][1] == -1000 || 
                mevcoo1[mevcooIndex[station]][2] == -1000)
                continue;
            
            // Determine paddle number
            int paddle = -1;
            
            // Get position along the paddle
            float alongPos = 0.0;
            
            if (isXAxisParallel[station]) {
                float yPos = mevcoo1[mevcooIndex[station]][1];
                
                for (int p = 0; p < nPaddles[station] - 1; ++p) {
                    float boundary = (p - nPaddles[station]/2 + 1) * paddleWidth;
                    if (yPos < boundary) {
                        paddle = p;
                        break;
                    }
                }
                
                if (paddle < 0) {
                    paddle = nPaddles[station] - 1;
                }
                
                // X position determines segment along the paddle length
                alongPos = mevcoo1[mevcooIndex[station]][0];
            } else {
                float xPos = mevcoo1[mevcooIndex[station]][0];
                
                for (int p = 0; p < nPaddles[station] - 1; ++p) {
                    float boundary = (p - nPaddles[station]/2 + 1) * paddleWidth;
                    if (xPos < boundary) {
                        paddle = p;
                        break;
                    }
                }
                
                if (paddle < 0) {
                    paddle = nPaddles[station] - 1;
                }
                
                // Y position determines segment along the paddle length
                alongPos = mevcoo1[mevcooIndex[station]][1];
            }
            
            if (paddle >= 0 && paddle < nPaddles[station]) {
                // Fill paddle histogram
                hTofTres[station][paddle]->Fill(tof_tres[station]);
                
                // Determine segment (0-9) based on position along paddle
                // Map alongPos from [-60cm, 60cm] to [0, nSegments-1]
                int segment = static_cast<int>((alongPos + paddleLength/2) / segmentLength);
                if (segment < 0) segment = 0;
                if (segment >= nSegments) segment = nSegments - 1;
                
                // Fill segment histogram
                hTofTresSegments[station][paddle][segment]->Fill(tof_tres[station]);
            }
        }
    }

    TCanvas *canvas = new TCanvas("canvas", "", 3508, 2480);
    canvas->SetLeftMargin(0.16);
    canvas->SetRightMargin(0.12);
    canvas->SetGridx();
    canvas->SetGridy();

    canvas->Print(Form("%s[", outputName));

    std::array<int, 10> paddleColors = {kRed, kBlue, kGreen+2, kMagenta, kCyan+2, 
                                      kOrange+7, kViolet-3, kSpring-5, kTeal-1, kYellow+2};
    
    // Draw paddle histograms (same as before)
    for (int station = 0; station < 4; ++station)
    {
        TLegend *legend = new TLegend(0.65, 0.50, 0.85, 0.85);
        legend->SetBorderSize(0);
        
        double maxHeight = 0;
        
        for (int paddle = 0; paddle < nPaddles[station]; ++paddle) {
            if (hTofTres[station][paddle]->GetMaximum() > maxHeight) {
                maxHeight = hTofTres[station][paddle]->GetMaximum();
            }
        }
        
        for (int paddle = 0; paddle < nPaddles[station]; ++paddle) {
            hTofTres[station][paddle]->SetLineColor(paddleColors[paddle % paddleColors.size()]);
            hTofTres[station][paddle]->SetLineWidth(2);
            
            if (paddle == 0) {
                hTofTres[station][paddle]->SetMaximum(maxHeight * 1.1);
                hTofTres[station][paddle]->Draw("HIST");
            } else {
                hTofTres[station][paddle]->Draw("HIST SAME");
            }
            
            legend->AddEntry(hTofTres[station][paddle], Form("Paddle %d", paddle+1), "l");
        }
        
        legend->Draw();
        
        TPaveText *stationText = new TPaveText(0.2, 0.92, 0.8, 0.98, "NDC");
        stationText->SetFillColor(0);
        stationText->SetBorderSize(0);
        stationText->AddText(Form("TOF Time Residual - %s", stationNames[station]));
        stationText->Draw();
        
        canvas->Print(outputName);
        
        delete legend;
        delete stationText;
        
        // Individual paddle plots with Gaussian fits
        for (int paddle = 0; paddle < nPaddles[station]; ++paddle)
        {
            canvas->Clear();
            
            if (hTofTres[station][paddle]->GetEntries() < 10) {
                continue;
            }
            
            hTofTres[station][paddle]->Draw("HIST");
            
            TF1 *fGaus = nullptr;
            TPaveText *statText = new TPaveText(0.65, 0.65, 0.85, 0.85, "NDC");
            statText->SetFillColor(0);
            statText->SetBorderSize(0);
            statText->AddText(Form("Entries: %.0f", hTofTres[station][paddle]->GetEntries()));
            statText->AddText(Form("Mean: %.4f", hTofTres[station][paddle]->GetMean()));
            statText->AddText(Form("RMS: %.4f", hTofTres[station][paddle]->GetRMS()));
            
            if (hTofTres[station][paddle]->GetEntries() >= 50) {
                fGaus = new TF1(Form("fGaus_%s_P%d", stationNames[station], paddle), 
                               "gaus", 
                               hTofTres[station][paddle]->GetXaxis()->GetXmin(), 
                               hTofTres[station][paddle]->GetXaxis()->GetXmax());
                fGaus->SetLineColor(kRed);
                hTofTres[station][paddle]->Fit(fGaus, "Q");
                fGaus->Draw("SAME");
                
                statText->AddText(Form("Fit Mean: %.4f", fGaus->GetParameter(1)));
                statText->AddText(Form("Fit Sigma: %.4f", fGaus->GetParameter(2)));
            }
            
            statText->Draw();
            
            TPaveText *paddleTitle = new TPaveText(0.2, 0.92, 0.8, 0.98, "NDC");
            paddleTitle->SetFillColor(0);
            paddleTitle->SetBorderSize(0);
            paddleTitle->AddText(Form("TOF %s, Paddle %d", stationNames[station], paddle+1));
            paddleTitle->Draw();
            
            canvas->Print(outputName);
            
            delete paddleTitle;
            delete statText;
            delete fGaus;
        }
    }
    
    // Calculate means and errors for each segment in each paddle
    std::vector<std::vector<std::vector<double>>> segmentMeans(4);
    std::vector<std::vector<std::vector<double>>> segmentErrors(4);
    
    for (int station = 0; station < 4; ++station) {
        segmentMeans[station].resize(nPaddles[station]);
        segmentErrors[station].resize(nPaddles[station]);
        
        for (int paddle = 0; paddle < nPaddles[station]; ++paddle) {
            segmentMeans[station][paddle].resize(nSegments, 0);
            segmentErrors[station][paddle].resize(nSegments, 0);
            
            for (int segment = 0; segment < nSegments; ++segment) {
                TH1F* hist = hTofTresSegments[station][paddle][segment];
                
                if (hist->GetEntries() >= 10) {
                    // Use Gaussian fit if enough entries
                    if (hist->GetEntries() >= 50) {
                        TF1* fSegGaus = new TF1(Form("fSegGaus_%s_P%d_S%d", stationNames[station], paddle, segment),
                                              "gaus", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
                        hist->Fit(fSegGaus, "Q");
                        segmentMeans[station][paddle][segment] = fSegGaus->GetParameter(1);
                        segmentErrors[station][paddle][segment] = fSegGaus->GetParError(1);
                        delete fSegGaus;
                    } else {
                        // Use histogram stats if not enough for fitting
                        segmentMeans[station][paddle][segment] = hist->GetMean();
                        segmentErrors[station][paddle][segment] = hist->GetRMS() / sqrt(hist->GetEntries());
                    }
                }
            }
        }
    }
    
    // Draw segment graphs for each paddle
    for (int station = 0; station < 4; ++station) {
        for (int paddle = 0; paddle < nPaddles[station]; ++paddle) {
            // Skip if no data
            bool hasData = false;
            for (int segment = 0; segment < nSegments; ++segment) {
                if (segmentMeans[station][paddle][segment] != 0) {
                    hasData = true;
                    break;
                }
            }
            
            if (!hasData) continue;
            
            // Prepare data for graph
            double segPos[nSegments];      // Segment position (center of segment)
            double segPosErr[nSegments];   // Position error (zero)
            double validMeans[nSegments];  // Mean values
            double validErrors[nSegments]; // Error values
            int validCount = 0;
            
            for (int segment = 0; segment < nSegments; ++segment) {
                if (segmentMeans[station][paddle][segment] != 0) {
                    // Position from -60cm to +60cm
                    segPos[validCount] = (segment + 0.5) * segmentLength - paddleLength/2;
                    segPosErr[validCount] = 0;
                    validMeans[validCount] = segmentMeans[station][paddle][segment];
                    validErrors[validCount] = segmentErrors[station][paddle][segment];
                    validCount++;
                }
            }
            
            if (validCount < 2) continue; // Need at least 2 points for a line
            
            // Create graph
            TGraphErrors* grSegments = new TGraphErrors(validCount, segPos, validMeans, segPosErr, validErrors);
            grSegments->SetTitle(Form("TOF Time Residual vs Position - %s Paddle %d", 
                                     stationNames[station], paddle+1));
            
            // Set axis labels based on paddle orientation
            if (isXAxisParallel[station]) {
                grSegments->GetXaxis()->SetTitle("Position along paddle (X) [cm]");
            } else {
                grSegments->GetXaxis()->SetTitle("Position along paddle (Y) [cm]");
            }
            
            grSegments->GetYaxis()->SetTitle("Time Residual Mean [ns]");
            grSegments->SetMarkerStyle(20);
            grSegments->SetMarkerSize(1.5);
            grSegments->SetMarkerColor(kBlue);
            grSegments->SetLineColor(kBlue);
            grSegments->SetLineWidth(2);
            
            // Find appropriate y-axis range
            double yMin = *std::min_element(validMeans, validMeans + validCount);
            double yMax = *std::max_element(validMeans, validMeans + validCount);
            double yPadding = (yMax - yMin) * 0.2;
            if (yPadding < 0.01) yPadding = 0.01; // Minimum padding
            
            canvas->Clear();
            grSegments->Draw("APL");
            grSegments->GetXaxis()->SetRangeUser(-paddleLength/2 - 10, paddleLength/2 + 10);
            grSegments->GetYaxis()->SetRangeUser(yMin - yPadding, yMax + yPadding);
            
            // Add title
            TPaveText *segTitle = new TPaveText(0.2, 0.92, 0.8, 0.98, "NDC");
            segTitle->SetFillColor(0);
            segTitle->SetBorderSize(0);
            segTitle->AddText(Form("TOF Time Residual vs Position - %s Paddle %d", 
                                  stationNames[station], paddle+1));
            segTitle->Draw();
            
            // Add info about orientation
            TPaveText *infoText = new TPaveText(0.65, 0.75, 0.85, 0.85, "NDC");
            infoText->SetFillColor(0);
            infoText->SetBorderSize(0);
            infoText->AddText(Form("Paddle orientation: %s", 
                                  isXAxisParallel[station] ? "X-Parallel" : "Y-Parallel"));
            infoText->AddText(Form("Valid segments: %d/%d", validCount, nSegments));
            infoText->Draw();
            
            canvas->Print(outputName);
            
            delete grSegments;
            delete segTitle;
            delete infoText;
        }
    }
    
    // Generate summary plot showing sigma vs paddle for each station
    std::vector<std::vector<double>> means(4);
    std::vector<std::vector<double>> sigmas(4);
    std::vector<std::vector<double>> errors(4);
    
    for (int station = 0; station < 4; ++station) {
        means[station].resize(nPaddles[station], 0);
        sigmas[station].resize(nPaddles[station], 0);
        errors[station].resize(nPaddles[station], 0);
        
        for (int paddle = 0; paddle < nPaddles[station]; ++paddle) {
            if (hTofTres[station][paddle]->GetEntries() < 10) {
                continue;
            }
            
            TF1 *fGaus = new TF1(Form("fGaus_summary_%s_P%d", stationNames[station], paddle), 
                                "gaus", 
                                hTofTres[station][paddle]->GetXaxis()->GetXmin(), 
                                hTofTres[station][paddle]->GetXaxis()->GetXmax());
            
            if (hTofTres[station][paddle]->GetEntries() >= 50) {
                hTofTres[station][paddle]->Fit(fGaus, "Q");
                means[station][paddle] = fGaus->GetParameter(1);
                sigmas[station][paddle] = fGaus->GetParameter(2);
                errors[station][paddle] = fGaus->GetParError(1);
            } else {
                means[station][paddle] = hTofTres[station][paddle]->GetMean();
                sigmas[station][paddle] = hTofTres[station][paddle]->GetRMS();
                errors[station][paddle] = sigmas[station][paddle] / sqrt(hTofTres[station][paddle]->GetEntries());
            }
            
            delete fGaus;
        }
    }
    
    for (int station = 0; station < 4; ++station) {
        bool hasData = false;
        for (int paddle = 0; paddle < nPaddles[station]; ++paddle) {
            if (means[station][paddle] != 0) {
                hasData = true;
                break;
            }
        }
        
        if (!hasData) continue;
        
        double paddleNums[nPaddles[station]];
        double validMeans[nPaddles[station]];
        double validSigmas[nPaddles[station]];
        double validErrors[nPaddles[station]];
        int validCount = 0;
        
        for (int paddle = 0; paddle < nPaddles[station]; ++paddle) {
            if (means[station][paddle] != 0) {
                paddleNums[validCount] = paddle + 1;
                validMeans[validCount] = means[station][paddle];
                validSigmas[validCount] = sigmas[station][paddle];
                validErrors[validCount] = errors[station][paddle];
                validCount++;
            }
        }
        
        if (validCount == 0) continue;
        
        // Create summary plots
        TPaveText *titleText = nullptr;
        
        // Mean plot
        TGraphErrors *grMeans = new TGraphErrors(validCount, paddleNums, validMeans, nullptr, validErrors);
        grMeans->SetTitle(Form("Mean TOF Time Residual - %s", stationNames[station]));
        grMeans->GetXaxis()->SetTitle("Paddle Number");
        grMeans->GetYaxis()->SetTitle("Time Residual Mean [ns]");
        grMeans->SetMarkerStyle(20);
        grMeans->SetMarkerSize(2.0);
        grMeans->SetMarkerColor(kBlue);
        grMeans->SetLineColor(kBlue);
        grMeans->SetLineWidth(2);
        
        double meanMin = *std::min_element(validMeans, validMeans + validCount);
        double meanMax = *std::max_element(validMeans, validMeans + validCount);
        double meanPadding = (meanMax - meanMin) * 0.2;
        grMeans->GetYaxis()->SetRangeUser(meanMin - meanPadding, meanMax + meanPadding);
        grMeans->GetXaxis()->SetLimits(0.5, nPaddles[station] + 0.5);
        
        canvas->Clear();
        grMeans->Draw("APL");
        
        titleText = new TPaveText(0.2, 0.92, 0.8, 0.98, "NDC");
        titleText->SetFillColor(0);
        titleText->SetBorderSize(0);
        titleText->AddText(Form("Mean TOF Time Residual - %s", stationNames[station]));
        titleText->Draw();
        
        canvas->Print(outputName);
        
        // Sigma plot
        TGraph *grSigmas = new TGraph(validCount, paddleNums, validSigmas);
        grSigmas->SetTitle(Form("Sigma of TOF Time Residual - %s", stationNames[station]));
        grSigmas->GetXaxis()->SetTitle("Paddle Number");
        grSigmas->GetYaxis()->SetTitle("Time Residual Sigma [ns]");
        grSigmas->SetMarkerStyle(24);
        grSigmas->SetMarkerSize(2.0);
        grSigmas->SetMarkerColor(kRed);
        grSigmas->SetLineColor(kRed);
        grSigmas->SetLineWidth(2);
        
        double sigmaMax = *std::max_element(validSigmas, validSigmas + validCount);
        grSigmas->GetYaxis()->SetRangeUser(0, sigmaMax * 1.2);
        grSigmas->GetXaxis()->SetLimits(0.5, nPaddles[station] + 0.5);
        
        canvas->Clear();
        grSigmas->Draw("APL");
        
        titleText->Clear();
        titleText->AddText(Form("Sigma of TOF Time Residual - %s", stationNames[station]));
        titleText->Draw();
        
        canvas->Print(outputName);
        
        delete grMeans;
        delete grSigmas;
        delete titleText;
    }
    
    canvas->Print(Form("%s]", outputName));

    std::cout << "TOF time residual plots saved to: " << outputName << std::endl;

    // Cleanup
    for (int station = 0; station < 4; ++station) {
        for (int paddle = 0; paddle < nPaddles[station]; ++paddle) {
            delete hTofTres[station][paddle];
            
            for (int segment = 0; segment < nSegments; ++segment) {
                delete hTofTresSegments[station][paddle][segment];
            }
        }
    }

    file->Close();
    delete file;
}