#include <TFile.h>
#include <TTree.h>

#include <iostream>

void makeValidationSample(const char *inputPath,
                          const char *outputPath = "results/validation_zeta_scope/mc_sample.root",
                          Long64_t entries = 10000)
{
    TFile input(inputPath, "READ");
    TTree *tree = static_cast<TTree *>(input.Get("amstreea"));
    if (!tree)
    {
        std::cerr << "Error: amstreea is missing from " << inputPath << std::endl;
        return;
    }

    TFile output(outputPath, "RECREATE");
    TTree *sample = tree->CopyTree("", "", entries, 0);
    sample->Write();
    std::cout << "Sample entries: " << sample->GetEntries() << std::endl;
}
