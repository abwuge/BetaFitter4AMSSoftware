#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Tools.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TMath.h>
#include <TCut.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void train_tracker_bdt_list(const char *input_path="tracker_features.root",
                            const char *output_dir="tmva_ablation",
                            const char *stage="subset",
                            const char *feature_list="subset.txt",
                            int rbin=0,
                            int ntrees=300)
{
    TMVA::Tools::Instance();
    gROOT->SetBatch(kTRUE);
    gSystem->mkdir(output_dir,kTRUE);

    std::vector<std::string> vars;
    std::ifstream fin(feature_list);
    std::string line;
    while (std::getline(fin,line)) {
        if (!line.empty()) vars.push_back(line);
    }
    if (vars.empty()) { std::cerr<<"empty feature list "<<feature_list<<std::endl; return; }

    TFile *input=TFile::Open(input_path,"READ");
    if(!input||input->IsZombie()){std::cerr<<"cannot open "<<input_path<<std::endl;return;}
    TTree *all=dynamic_cast<TTree*>(input->Get("events"));
    if(!all){std::cerr<<"missing events"<<std::endl;return;}
    TString train_cut=Form("role<2 && bin==%d",rbin);
    TString test_cut=Form("role==2 && bin==%d",rbin);
    gROOT->cd();
    TTree *train=all->CopyTree(train_cut);
    TTree *test=all->CopyTree(test_cut);

    TString tag=Form("%s_bin%d",stage,rbin);
    TString output_name=TString(output_dir)+"/"+tag+".root";
    TFile *output=TFile::Open(output_name,"RECREATE");
    TString factory_name="TrackerOracle_"+tag;
    TMVA::Factory *factory=new TMVA::Factory(factory_name,output,"!V:Silent:!Color:!DrawProgressBar:AnalysisType=Regression");
    TString dataset_name="ds_"+tag;
    TMVA::DataLoader *loader=new TMVA::DataLoader(dataset_name);
    for(const auto &v:vars) loader->AddVariable(v.c_str(),'F');
    loader->AddSpectator("event_index","event_index","",'I');
    loader->AddSpectator("label","label","",'I');
    loader->AddSpectator("source","source","",'I');
    loader->AddSpectator("R","R","GV",'F');
    loader->AddSpectator("R_true","R_true","GV",'F');
    loader->AddSpectator("u","u","",'F');
    loader->AddSpectator("u_true","u_true","",'F');
    loader->AddTarget("target");
    loader->AddRegressionTree(train,1.0,TMVA::Types::kTraining);
    loader->AddRegressionTree(test,1.0,TMVA::Types::kTesting);
    TCut empty="";
    loader->PrepareTrainingAndTestTree(empty,"nTrain_Regression=0:nTest_Regression=0:SplitMode=Block:NormMode=None:!V");

    int nvars=vars.size();
    int use_nvars=TMath::Max(16,TMath::Min(64,int(2.0*TMath::Sqrt(double(nvars))+0.5)));
    TString opts=Form("!H:!V:NTrees=%d:MinNodeSize=1.0%%:BoostType=Grad:Shrinkage=0.05:"
                      "UseBaggedBoost=True:BaggedSampleFraction=0.65:nCuts=32:MaxDepth=3:"
                      "RegressionLossFunctionBDTG=Huber:HuberQuantile=0.70:"
                      "UseRandomisedTrees=True:UseNvars=%d:UsePoissonNvars=False",
                      ntrees,use_nvars);
    std::cout<<"TRAIN "<<tag<<" nvars="<<nvars<<" random_nvars="<<use_nvars
             <<" ntrain="<<train->GetEntries()<<" ntest="<<test->GetEntries()<<std::endl;
    factory->BookMethod(loader,TMVA::Types::kBDT,"BDTG",opts);
    factory->TrainAllMethods(); factory->TestAllMethods(); factory->EvaluateAllMethods();
    output->Close();
    delete factory; delete loader; delete train; delete test;
    input->Close(); delete input; delete output;
    std::cout<<"DONE "<<output_name<<std::endl;
}
