#include <string>
#include <fstream>
#include <iostream>

#include <TFile.h>
#include <TTree.h>

void selectHadd(std::string inputFileList = "input_Z2.list", const char *outputFile = "test_sel_Z2.root")
{
    if (inputFileList.size() < 13 || inputFileList.substr(inputFileList.size() - 5) != ".list")
        inputFileList = "input_Z" + inputFileList + ".list";

    std::ifstream inputFile(inputFileList);
    if (!inputFile)
    {
        std::cerr << "Error: Unable to open input file " << inputFileList << std::endl;
        return;
    }

    TFile *output = new TFile(outputFile, "RECREATE");
    if (!output || output->IsZombie())
    {
        std::cerr << "Error: Unable to create output file " << outputFile << std::endl;
        inputFile.close();
        return;
    }

    TTree *outputTree = new TTree("amstrees", "Selected amstreea");
    float mevmom1[21]{};
    float mevcoo1[21][3]{};
    float tof_betah{};

    outputTree->Branch("mevmom1", mevmom1, "mevmom1[21]/F");
    outputTree->Branch("mevcoo1", mevcoo1, "mevcoo1[21][3]/F");
    outputTree->Branch("tof_betah", &tof_betah, "tof_betah/F");

    std::string line;
    while (std::getline(inputFile, line))
    {
        if (line.empty() || line[0] == '#')
            continue;

        std::cout << "Processing file: " << line << std::endl;

        TFile *input = TFile::Open(line.c_str(), "READ");
        if (!input || input->IsZombie())
        {
            std::cerr << "Error: Unable to open file " << line << std::endl;
            continue;
        }

        TTree *inputTree = nullptr;
        input->GetObject("amstreea", inputTree);
        if (!inputTree)
        {
            std::cerr << "Error: amstreea not found in " << line << std::endl;
            input->Close();
            continue;
        }

        inputTree->SetBranchStatus("*", false);
        inputTree->SetBranchStatus("mevmom1", true);
        inputTree->SetBranchStatus("mevcoo1", true);
        inputTree->SetBranchStatus("tof_betah", true);

        inputTree->SetBranchAddress("mevmom1", mevmom1);
        inputTree->SetBranchAddress("mevcoo1", mevcoo1);
        inputTree->SetBranchAddress("tof_betah", &tof_betah);

        outputTree->CopyEntries(inputTree);

        std::cout << "Processed " << inputTree->GetEntries() << " entries\n";
    }

    output->Write();
    output->Close();
    inputFile.close();
}