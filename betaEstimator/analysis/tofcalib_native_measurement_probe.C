// Reuse the native producer's measurement methods. Emit only missing branches
// and association checks; existing NEWTOF measurements remain authoritative.
#define main native_producer_main
#include "../../amsd69n/analysis_amsd69n.C"
#undef main

int main(int argc, char **argv) {
    if (argc != 4) { std::cerr << "usage: supplement raw.root targets.txt output.root\n"; return 2; }
    AMSChain chain;
    std::string input = argv[1];
    if (input.size() >= 5 && input.substr(input.size()-5) == ".root") {
        if (chain.Add(argv[1]) <= 0) return 3;
    } else {
        std::ifstream files(argv[1]);
        std::string path;
        int count=0;
        while (std::getline(files, path)) {
            if (path.empty()) continue;
            if (chain.Add(path.c_str()) <= 0) return 3;
            ++count;
        }
        if (!count) return 3;
    }
    AMSAnalysis analysis;
    analysis.LoadHeader();
    TFile output(argv[3], "RECREATE");
    if (output.IsZombie()) return 4;
    TTree tree("beta_extra", "Native measured supplements keyed by actual run/event");
    unsigned int run, event;
    int found, selected;
    tree.Branch("run", &run, "run/i");
    tree.Branch("event", &event, "event/i");
    tree.Branch("found", &found, "found/I");
    tree.Branch("selected", &selected, "selected/I");
    tree.Branch("itrtrack", &analysis.itrtrack, "itrtrack/I");
    tree.Branch("ibetah", &analysis.ibetah, "ibetah/I");
    tree.Branch("irich", &analysis.irich, "irich/I");
    tree.Branch("tk_rigidityk", analysis.tk_rigidityk, "tk_rigidityk[4][7]/F");
    tree.Branch("tk_erigidity1", analysis.tk_erigidity1, "tk_erigidity1[3][7]/F");
    tree.Branch("tof_bpatt", &analysis.tof_bpatt, "tof_bpatt/I");
    tree.Branch("rich_itrtrack", &analysis.rich_itrtrack, "rich_itrtrack/I");
    tree.Branch("rich_ebeta1", &analysis.rich_ebeta1, "rich_ebeta1/F");
    tree.Branch("tof_tl", analysis.tof_tl, "tof_tl[4]/F");
    tree.Branch("tof_leng", analysis.tof_leng, "tof_leng[4]/F");
    tree.Branch("tof_edep", analysis.tof_edep, "tof_edep[4]/F");
    tree.Branch("tof_pos", analysis.tof_pos, "tof_pos[4][3]/F");
    tree.Branch("tof_barid", analysis.tof_barid, "tof_barid[4]/I");
    tree.Branch("tk_ql", analysis.tk_ql, "tk_ql[9]/F");
    std::ifstream targets(argv[2]);
    if (!targets) return 5;
    int total=0, missing=0;
    while (targets >> run >> event) {
        AMSEventR *ev = chain.GetEventFast(run, event);
        found = ev && ev->Run()==run && ev->Event()==event;
        selected=0;
        if (found) {
            ev->RecordRTIRun();
            analysis.unbiasetag0=1;
            analysis.unbiasetag=0;
            analysis.istkaddhit=0;
            analysis.InitEvent(ev);
            selected=analysis.Select();
        } else { ++missing; }
        tree.Fill();
        if (++total % 100 == 0) std::cout << "Supplement " << total << " missing " << missing << std::endl;
    }
    output.cd();
    tree.Write();
    output.Close();
    std::cout << "DONE requested=" << total << " missing=" << missing << std::endl;
    return missing ? 6 : 0;
}
