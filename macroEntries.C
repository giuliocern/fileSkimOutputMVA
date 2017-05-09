#include <iostream>
#include <vector>
#include <string>
#include <fstream>


void macroEntries() {

    std::vector<std::string> fileNames;

    fileNames.push_back("DYJetstoLL_amc_0J");
    fileNames.push_back("DYJetstoLL_amc_1J");
    fileNames.push_back("DYJetstoLL_amc_2J");
    fileNames.push_back("DYJetsToLL_M-50");
    fileNames.push_back("DY1JetsToLL_M-50");
    fileNames.push_back("DY2JetsToLL_M-50");
    fileNames.push_back("DY3JetsToLL_M-50");
    fileNames.push_back("DY4JetsToLL_M-50");
    ofstream out_entries;
    out_entries.open("entries.txt"); 





    for (int i = 0; i < fileNames.size(); ++i) {
        std::string thisFileName = fileNames[i] + "_v25_reskim.root";
        TFile * f = new TFile (thisFileName.c_str());
        TH1F * histo_entries = (TH1F*) f->Get("CountWeighted");
        float count = histo_entries->GetBinContent(1);

        TTree * my_tree;
        f->GetObject("tree", my_tree);
        int TEnt = my_tree->GetEntries();

        out_entries << " "<< fileNames[i] << "\t\ttree entries:  " <<  TEnt << "\t\tCount events:  " << count << "\t\tratio:  " << TEnt/count << endl;

    }


    out_entries.close();

}
