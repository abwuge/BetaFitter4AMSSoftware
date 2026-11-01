#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Tools.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCut.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void train_reverse_bdt(const char *input_path="reverse_features.root",
                       const char *output_dir="tmva_reverse",
                       const char *stage="reverse",
                       const char *feature_list="features_reverse_eloss.txt",
                       int rbin=0, int ntrees=250)
{
  TMVA::Tools::Instance(); gROOT->SetBatch(kTRUE); gSystem->mkdir(output_dir,kTRUE);
  std::vector<std::string> vars; std::ifstream fin(feature_list); std::string line;
  while(std::getline(fin,line)) if(!line.empty()) vars.push_back(line);
  TFile *input=TFile::Open(input_path,"READ"); TTree *all=(TTree*)input->Get("events");
  TString train_cut=Form("role==0 && bin==%d",rbin);
  TString test_cut =Form("role>0 && bin==%d",rbin);
  gROOT->cd(); TTree *train=all->CopyTree(train_cut); TTree *test=all->CopyTree(test_cut);
  TString tag=Form("%s_bin%d",stage,rbin); TString outname=TString(output_dir)+"/"+tag+".root";
  TFile *output=TFile::Open(outname,"RECREATE");
  TMVA::Factory *factory=new TMVA::Factory("Reverse_"+tag,output,"!V:Silent:!Color:!DrawProgressBar:AnalysisType=Regression");
  TMVA::DataLoader *loader=new TMVA::DataLoader("ds_"+tag);
  for(const auto &v:vars) loader->AddVariable(v.c_str(),'F');
  loader->AddSpectator("event_index","event_index","",'I');
  loader->AddSpectator("role","role","",'I');
  loader->AddSpectator("source","source","",'I');
  loader->AddSpectator("R","R","GV",'F');
  loader->AddSpectator("R_true","R_true","GV",'F');
  loader->AddTarget("target");
  loader->AddRegressionTree(train,1.0,TMVA::Types::kTraining);
  loader->AddRegressionTree(test,1.0,TMVA::Types::kTesting);
  TCut empty="";
  loader->PrepareTrainingAndTestTree(empty,"nTrain_Regression=0:nTest_Regression=0:SplitMode=Block:NormMode=None:!V");
  int use_nvars=vars.size()<64?vars.size():64;
  TString opts=Form("!H:!V:NTrees=%d:MinNodeSize=1.0%%:BoostType=Grad:Shrinkage=0.05:"
                    "UseBaggedBoost=True:BaggedSampleFraction=0.65:nCuts=32:MaxDepth=3:"
                    "RegressionLossFunctionBDTG=Huber:HuberQuantile=0.70:"
                    "UseRandomisedTrees=True:UseNvars=%d:UsePoissonNvars=False",ntrees,use_nvars);
  std::cout<<"TRAIN "<<tag<<" vars="<<vars.size()<<" train="<<train->GetEntries()<<" valid+test="<<test->GetEntries()<<std::endl;
  factory->BookMethod(loader,TMVA::Types::kBDT,"BDTG",opts);
  factory->TrainAllMethods(); factory->TestAllMethods(); factory->EvaluateAllMethods();
  output->Close(); delete factory; delete loader; delete train; delete test; input->Close(); delete input; delete output;
  std::cout<<"DONE "<<outname<<std::endl;
}
