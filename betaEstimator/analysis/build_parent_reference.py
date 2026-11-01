#!/usr/bin/env python3
"""Generate a parent-bound raw reader while retaining the native output schema."""
from pathlib import Path
import hashlib
import json
import subprocess
from research_conditional_mc import ROOT,OUT,stamp
from research_flight_conditional import write,digest

PREFIX=r'''
#include <map>
#include <memory>
class ParentReferenceReader {
 std::unique_ptr<TFile> file;
 TTree *tree;
 AMSEventR *event;
 TBranch *header;
 std::map<std::pair<unsigned int,unsigned int>,Long64_t> index;
public:
 ParentReferenceReader(const char *path):file(TFile::Open(path,"READ")),tree(0),event(0),header(0) {
  if(!file || file->IsZombie())return;
  tree=dynamic_cast<TTree*>(file->Get("AMSRoot"));if(!tree)return;
  event=new AMSEventR;
  tree->SetBranchStatus("*",0);
  const char *names[]={"ev.fHeader","ev.fMCEventg","ev.fTofMCCluster","ev.fTrMCCluster","ev.fBetaH","ev.fTrTrack"};
  for(int j=0;j<6;++j){UInt_t found=0;tree->SetBranchStatus(names[j],1,&found);if(found!=1){tree=0;return;}}
  if(tree->SetBranchAddress("ev.",&event)<0){tree=0;return;}
  header=tree->GetBranch("ev.fHeader");if(!header){tree=0;return;}
  for(Long64_t entry=0;entry<tree->GetEntries();++entry){
   if(header->GetEntry(entry)<=0){std::cerr<<"HEADER_FAILURE "<<entry<<std::endl;tree=0;return;}
   std::pair<unsigned int,unsigned int> key(event->Run(),event->Event());
   if(index.count(key)){std::cerr<<"DUPLICATE_IDENTITY\n";tree=0;return;}
   index[key]=entry;
  }
 }
 bool valid()const{return tree!=0;}
 AMSEventR *GetEventFast(unsigned int run,unsigned int number){
  auto it=index.find(std::make_pair(run,number));if(it==index.end())return 0;
  if(tree->GetEntry(it->second)<=0)return 0;
  return event;
 }
};
'''

def main():
 out=OUT/'reference';source=ROOT/'analysis/conditional_mc_reference_parent.C'
 if source.exists():raise FileExistsError(str(source))
 original=(ROOT/'analysis/conditional_mc_reference.C').read_text()
 begin=original.index('    AMSChain chain;');end=original.index('    TFile output',begin)
 replacement='    ParentReferenceReader chain(argv[1]);if(!chain.valid())return 3;\n'
 changed=original[:begin]+replacement+original[end:]
 changed=changed.replace('int main(int argc,char **argv)',PREFIX+'\nint main(int argc,char **argv)',1)
 changed=changed.replace('AMSPoint p=v->GetXgl();for(int k=0;k<3;++k)tk_position.push_back(p[k]);',
                         'for(int k=0;k<3;++k)tk_position.push_back(nan); // No Tracker geometry initialization in this reader.')
 source.write_text(changed)
 build=(ROOT/'analysis/build_conditional_mc_reference_branches.sh').read_text().replace('conditional_mc_reference_branches.C','conditional_mc_reference_parent.C').replace('reference/extract_branches','reference/extract_parent')
 script=ROOT/'analysis/build_conditional_mc_reference_parent.sh';script.write_text(build);script.chmod(0o755)
 result=subprocess.run(['bash',str(script)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
 (out/'parent_build.log').write_text(result.stdout)
 if result.returncode:raise RuntimeError(result.stdout)
 write(out/'parent_decoder.json',dict(created_utc=stamp(),source_sha256=digest(source),executable_sha256=digest(out/'extract_parent'),
   base_extractor_sha256=digest(ROOT/'analysis/conditional_mc_reference.C'),build_script_sha256=digest(script),
   read_branches=['Header','MCEventg','TofMCCluster','TrMCCluster','BetaH','TrTrack'],
   deliberately_unread_output='tk_position is NaN because the Tracker geometry database is not initialized; L5 primary acceptance coordinates remain native MCEventg values',
   reference_selection_changed=False,model_parameters_changed=False,validated=False))
 print('PARENT_DECODER_BUILT',flush=True)

if __name__=='__main__':main()
