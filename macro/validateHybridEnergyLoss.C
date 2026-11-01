#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <TSystem.h>

namespace
{
void printResidualSummary(const char *name, const std::vector<double> &values)
{
    if (values.empty())
    {
        std::cout << name << " 0 nan nan nan nan" << std::endl;
        return;
    }

    double sum = 0;
    double sum2 = 0;
    for (std::vector<double>::const_iterator it = values.begin();
         it != values.end(); ++it)
    {
        sum += *it;
        sum2 += *it * *it;
    }
    const double mean = sum / values.size();
    const double rms = std::sqrt(std::max(0.0, sum2 / values.size() - mean * mean));

    std::vector<double> sorted(values);
    std::sort(sorted.begin(), sorted.end());
    const std::size_t middle = sorted.size() / 2;
    const std::size_t low = static_cast<std::size_t>(0.16 * (sorted.size() - 1));
    const std::size_t high = static_cast<std::size_t>(0.84 * (sorted.size() - 1));
    const double median = sorted[middle];
    const double central68 = 0.5 * (sorted[high] - sorted[low]);

    std::cout << name << " " << values.size() << " " << mean << " " << rms
              << " " << median << " " << central68 << std::endl;
}
}

// Residual convention: 1/beta_reco - 1/beta_MC.  Output columns are
// method, finite entries, mean, RMS, median, central-68% half width.
void validateHybridEnergyLoss(const char *inputPath,
                              double betaLow = 0,
                              double betaHigh = 1.000001)
{
    TFile input(inputPath, "READ");
    TTree *tree = static_cast<TTree *>(input.Get("betaTree"));
    if (!tree)
    {
        std::cerr << "Error: betaTree is missing from " << inputPath << std::endl;
        return;
    }

    double mcBeta = 0;
    double linearBeta = 0;
    double nonlinearBeta = 0;
    double hybridBeta = 0;
    tree->SetBranchAddress("mcBeta", &mcBeta);
    tree->SetBranchAddress("linearBeta", &linearBeta);
    tree->SetBranchAddress("nonlinearBeta", &nonlinearBeta);
    tree->SetBranchAddress("hybridBeta", &hybridBeta);

    std::vector<double> linear;
    std::vector<double> nonlinear;
    std::vector<double> hybrid;
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
    {
        tree->GetEntry(entry);
        if (!(mcBeta >= betaLow && mcBeta < betaHigh) ||
            !(mcBeta > 0) || !std::isfinite(mcBeta))
            continue;
        const double truth = 1.0 / mcBeta;
        if (linearBeta > 0 && std::isfinite(linearBeta))
            linear.push_back(1.0 / linearBeta - truth);
        if (nonlinearBeta > 0 && std::isfinite(nonlinearBeta))
            nonlinear.push_back(1.0 / nonlinearBeta - truth);
        if (hybridBeta > 0 && std::isfinite(hybridBeta))
            hybrid.push_back(1.0 / hybridBeta - truth);
    }

    std::cout << "# " << inputPath << " beta=[" << betaLow << ","
              << betaHigh << ")" << std::endl;
    std::cout << "# method entries mean rms median central68" << std::endl;
    std::cout << std::setprecision(10);
    printResidualSummary("linear", linear);
    printResidualSummary("nonlinear", nonlinear);
    printResidualSummary("hybrid", hybrid);
}

// Inputs must be reconstructions of the same source file in the same order.
void checkHybridStepConvergence(const char *fullStepPath,
                               const char *halfStepPath)
{
    TFile full(fullStepPath),half(halfStepPath);
    auto a=static_cast<TTree *>(full.Get("betaTree"));
    auto b=static_cast<TTree *>(half.Get("betaTree"));
    if(!a || !b || a->GetEntries()!=b->GetEntries())
    { std::cerr<<"Convergence inputs are not paired"<<std::endl;gSystem->Exit(1);return; }
    double truth[2],linear[2],nonlinear[2],hybrid[2];
    int status[2],steps[2];
    TTree *trees[2]={a,b};
    for(int i=0;i<2;++i)
    {
        trees[i]->SetBranchAddress("mcBeta",&truth[i]);
        trees[i]->SetBranchAddress("linearBeta",&linear[i]);
        trees[i]->SetBranchAddress("nonlinearBeta",&nonlinear[i]);
        trees[i]->SetBranchAddress("hybridBeta",&hybrid[i]);
        trees[i]->SetBranchAddress("hybridFitStatus",&status[i]);
        trees[i]->SetBranchAddress("hybridSteps",&steps[i]);
    }
    const auto same=[](double x,double y){return x==y || (std::isnan(x)&&std::isnan(y));};
    std::vector<double> difference,absolute;
    Long64_t changedReference=0;
    for(Long64_t entry=0;entry<a->GetEntries();++entry)
    {
        a->GetEntry(entry);b->GetEntry(entry);
        if(!same(truth[0],truth[1]) || !same(linear[0],linear[1]) ||
           !same(nonlinear[0],nonlinear[1])) ++changedReference;
        if(status[0] || status[1] || !(hybrid[0]>0 && hybrid[1]>0) ||
           !std::isfinite(hybrid[0]) || !std::isfinite(hybrid[1])) continue;
        double delta=1.0/hybrid[1]-1.0/hybrid[0];
        difference.push_back(delta);absolute.push_back(std::abs(delta));
    }
    std::sort(absolute.begin(),absolute.end());
    std::cout<<"STEP_CONVERGENCE "<<fullStepPath<<" entries="<<a->GetEntries()
             <<" referenceChanges="<<changedReference<<" steps="<<steps[0]<<":"<<steps[1]<<std::endl;
    printResidualSummary("half_minus_full",difference);
    if(!absolute.empty())
        std::cout<<"absolute_delta median="<<absolute[absolute.size()/2]
                 <<" p95="<<absolute[static_cast<size_t>(0.95*(absolute.size()-1))]
                 <<" max="<<absolute.back()<<std::endl;
    if(changedReference || absolute.empty()) gSystem->Exit(1);
}
