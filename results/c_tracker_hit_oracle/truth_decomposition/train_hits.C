#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Tools.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TCut.h>
#include <fstream>
#include <string>
void train_hits(const char*stage="adc"){
 TMVA::Tools::Instance();gROOT->SetBatch(1);
 TFile input("hit_features.root");TTree*all=(TTree*)input.Get("events");gROOT->cd();
 TTree*tr=all->CopyTree("role==0");TTree*te=all->CopyTree("role!=0");
 TFile out(Form("hit_%s.root",stage),"RECREATE");
 TMVA::Factory fac(Form("Hit_%s",stage),&out,"!V:Silent:!Color:!DrawProgressBar:AnalysisType=Regression");
 TMVA::DataLoader loader(Form("ds_hit_%s",stage));std::ifstream f(Form("hit_%s.txt",stage));std::string s;
 while(std::getline(f,s))if(!s.empty())loader.AddVariable(s.c_str(),'F');
 loader.AddSpectator("hit_index",'I');loader.AddSpectator("role",'I');loader.AddSpectator("source",'I');loader.AddSpectator("true_error_um",'F');loader.AddTarget("target");
 loader.AddRegressionTree(tr,1.,TMVA::Types::kTraining);loader.AddRegressionTree(te,1.,TMVA::Types::kTesting);
 loader.PrepareTrainingAndTestTree(TCut(""),"nTrain_Regression=0:nTest_Regression=0:SplitMode=Block:NormMode=None:!V");
 fac.BookMethod(&loader,TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=200:MinNodeSize=1.0%:BoostType=Grad:Shrinkage=0.05:UseBaggedBoost=True:BaggedSampleFraction=0.65:nCuts=32:MaxDepth=3:RegressionLossFunctionBDTG=Huber:HuberQuantile=0.70");
 fac.TrainAllMethods();fac.TestAllMethods();fac.EvaluateAllMethods();out.Close();
}
