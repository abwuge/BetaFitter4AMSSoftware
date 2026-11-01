#!/usr/bin/env python3
"""Build a direct ROOT-branch reference reader alongside the AMSChain audit."""
import hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
source=root/'analysis/conditional_mc_reference.C'
if hashlib.sha256(source.read_bytes()).hexdigest()!='a845a52e88c69a09ad05d4b5402c5796f9a2929ad8dc929dddddff3949f0629f':
    raise ValueError('Unexpected original reference reader')
text=source.read_text()
reader=r'''
#include <map>
#include <memory>
class DirectReferenceReader {
    std::unique_ptr<TFile> file;
    TTree *tree;
    HeaderR *header;
    std::vector<MCEventgR> *mc;
    std::vector<TofMCClusterR> *tof;
    std::vector<TrMCClusterR> *tk;
    TBranch *bh,*bm,*bt,*bk;
    std::map<std::pair<unsigned int,unsigned int>,Long64_t> index;
public:
    DirectReferenceReader(const char *path):file(TFile::Open(path,"READ")),tree(0),header(0),mc(0),tof(0),tk(0),bh(0),bm(0),bt(0),bk(0) {
        if(!file || file->IsZombie())return;
        tree=dynamic_cast<TTree*>(file->Get("AMSRoot"));if(!tree)return;
        tree->SetBranchStatus("*",0);
        const char *names[]={"ev.fHeader","ev.fMCEventg","ev.fTofMCCluster","ev.fTrMCCluster"};
        for(int i=0;i<4;++i)tree->SetBranchStatus(names[i],1);
        if(tree->SetBranchAddress(names[0],&header)<0 || tree->SetBranchAddress(names[1],&mc)<0 ||
           tree->SetBranchAddress(names[2],&tof)<0 || tree->SetBranchAddress(names[3],&tk)<0){tree=0;return;}
        bh=tree->GetBranch(names[0]);bm=tree->GetBranch(names[1]);bt=tree->GetBranch(names[2]);bk=tree->GetBranch(names[3]);
        if(!bh || !bm || !bt || !bk){tree=0;return;}
        for(Long64_t entry=0;entry<tree->GetEntries();++entry){
            if(bh->GetEntry(entry)<=0 || !header){tree=0;return;}
            std::pair<unsigned int,unsigned int> key(header->Run,header->Event);
            if(index.count(key)){std::cerr<<"Duplicate raw identity\n";tree=0;return;}
            index[key]=entry;
        }
    }
    bool valid() const{return tree!=0;}
    DirectReferenceReader *GetEventFast(unsigned int run,unsigned int event){
        auto it=index.find(std::make_pair(run,event));if(it==index.end())return 0;
        Long64_t e=it->second;
        if(bh->GetEntry(e)<=0 || bm->GetEntry(e)<=0 || bt->GetEntry(e)<=0 || bk->GetEntry(e)<=0)return 0;
        return this;
    }
    unsigned int Run()const{return header->Run;}
    unsigned int Event()const{return header->Event;}
    int nMCEventg()const{return mc?mc->size():0;}
    int nTofMCCluster()const{return tof?tof->size():0;}
    int nTrMCCluster()const{return tk?tk->size():0;}
    MCEventgR *pMCEventg(int i){return mc && i>=0 && i<(int)mc->size()?&(*mc)[i]:0;}
    TofMCClusterR *pTofMCCluster(int i){return tof && i>=0 && i<(int)tof->size()?&(*tof)[i]:0;}
    TrMCClusterR *pTrMCCluster(int i){return tk && i>=0 && i<(int)tk->size()?&(*tk)[i]:0;}
};
'''
text=text.replace('#undef main','#undef main\n'+reader)
begin=text.index('    AMSChain chain;')
end=text.index('    TFile output',begin)
text=text[:begin]+'    DirectReferenceReader chain(argv[1]);if(!chain.valid())return 3;\n'+text[end:]
text=text.replace('AMSEventR *ev=chain.GetEventFast','DirectReferenceReader *ev=chain.GetEventFast')
old='''            track_exists=ev->pTrTrack(requested_track)!=0;
            BetaHR *betah=ev->pBetaH(requested_beta);beta_exists=betah!=0;if(betah)beta_track=betah->iTrTrack();'''
if text.count(old)!=1:raise ValueError('Original association block differs')
text=text.replace(old,'            track_exists=beta_exists=beta_track=-1; // Measurement branches are not reopened.')
text=text.replace('raw_file_or_list targets output.root','raw.root targets output.root')
text=text.replace('// Extract primary reference records after measured-only MC models are frozen.','// Read only raw reference branches; bypass AMSEventR detector initialization.')
output=root/'analysis/conditional_mc_reference_branches.C'
with output.open('x') as stream:stream.write(text)
print(output,hashlib.sha256(output.read_bytes()).hexdigest())
