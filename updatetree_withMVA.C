#define TMVAapply_all_cxx
//#include "TMVAapply_all.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
//#include "/afs/cern.ch/user/g/gimandor/private/CMSSW_8_0_25/src/nadyaCode/TMVA-v4.2.0/Factory.h"
//#include "/afs/cern.ch/user/g/gimandor/private/CMSSW_8_0_25/src/nadyaCode/TMVA-v4.2.0/Tools.h"
//#include "/afs/cern.ch/user/g/gimandor/private/CMSSW_8_0_25/src/nadyaCode/TMVA-v4.2.0/Reader.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
//#include "/afs/cern.ch/work/n/nchernya/VBFZll/plotter/EWcorr.C"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.h"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/rochcor2016.cc"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.cc"
//#include "/mnt/t3nfs01/data01/shome/nchernya/VBFZll/plotter/RoccoR.h"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>


const int njets = 30;




typedef struct {
   Float_t eta[njets];
   Float_t pt[njets];
   Float_t JEC_corr[njets];
   Float_t JEC_corr_up[njets];
   Float_t JEC_corr_down[njets];
   Float_t JER_corr[njets];
   Float_t JER_corr_up[njets];
   Float_t JER_corr_down[njets];
   Float_t phi[njets];
	Float_t mass[njets];
	Float_t btagCSV[njets];
        Float_t btagCMVA[njets];
	Int_t nsoft;
	Float_t soft_pt[njets];
	Float_t soft_eta[njets];
	Float_t soft_mass[njets];
	Float_t qgl[njets];
	Int_t nsoft2;
	Int_t nsoft5;
	Int_t nsoft10;
        Int_t nsoftActivityEWKJets;
        Int_t EWKnsoft2;
	Int_t EWKnsoft5;
	Int_t EWKnsoft10;
	Int_t id[njets];
	Int_t puId[njets];
	Float_t HTsoft;
	Int_t partonFlavour[njets];
	Float_t EWKHTsoft;
	Float_t EWKsoft_pt[njets];
	Float_t EWKsoft_eta[njets];
	Float_t EWKsoft_phi[njets];
	Float_t EWKsoft_mass[njets];
	Float_t pt_regVBF[njets];	
	Float_t ptd[njets];
	Float_t axis2[njets];
	Int_t mult[njets];
	Float_t leadTrackPt[njets];
	Float_t blike_VBF[njets];
} Jets;


typedef struct {

        Float_t ll_mass;
//         Float_t Zll_pt;
   	Float_t Mqq; 
	Float_t DeltaEtaQQ;
  	Float_t q1_eta;
        Float_t met_pt;
        Float_t EWKHTsoft;
        Float_t btagCMVA;
        Float_t btagCSV;
	Float_t softJet3_pt;
        Float_t cosThetaStar;
	Float_t cosThetaStarAbs;

	Float_t qq_pt;
	Float_t Jet3_pt;
	Float_t axis2_jet1;
	Float_t axis2_jet2;
	Float_t Jet2q_pt; 
	Float_t Jet2q_leadTrackPt;
	Float_t Jet1q_pt;
	Float_t Jet1q_leadTrackPt;
	Float_t ll_zstar;
	Float_t RptHard;
	Float_t ll_pt;
	Float_t qgl_1q;
	Float_t qgl_2q;
        
        Int_t softActivityEWK_njets2;
        Int_t softActivityEWK_njets5;
        Int_t softActivityEWK_njets10;

        float weightMVA;	

        
}TMVAstruct;

void updatetree_SingleFile(std::string fileName) {

    TFile *f = new TFile(fileName.c_str(),"update");
    TTree *T = (TTree*)f->Get("tree");


    Jets Jet;
    int nJets=0;
    float met_pt=0.;
    int EWKnsoft5=0;
    int nselLeptons=0;
	float selLeptons_pt[30], selLeptons_eta[30], selLeptons_phi[30], selLeptons_mass[30], selLeptons_SF_IdCutLoose[30], selLeptons_SF_IdCutTight[30], selLeptons_SF_IsoLoose[30], selLeptons_SF_IsoTight[30],selLeptons_SF_trk_eta[30], selLeptons_SF_HLT_RunD4p2[30],selLeptons_SF_HLT_RunD4p3[30], selLeptons_relIso04[30], selLeptons_relIso03[30], selLeptons_eleSieie[30], selLeptons_eleHoE[30], selLeptons_eleDEta[30],selLeptons_eleDPhi[30], selLeptons_eleEcalClusterIso[30], selLeptons_eleHcalClusterIso[30],selLeptons_dr03TkSumPt[30] ;
	int selLeptons_charge[30], selLeptons_pdgId[30], selLeptons_looseIdPOG[30], selLeptons_trackerLayers[30],  selLeptons_eleMVAIdSppring16GenPurp[30]; 



    T->SetBranchAddress("nJet",&nJets);
    T->SetBranchAddress("Jet_pt",Jet.pt);
    T->SetBranchAddress("Jet_eta",Jet.eta);
    T->SetBranchAddress("Jet_phi",Jet.phi);
	T->SetBranchAddress("Jet_mass",Jet.mass);
	T->SetBranchAddress("Jet_btagCSV",Jet.btagCSV);
	T->SetBranchAddress("Jet_btagCMVA",Jet.btagCMVA);
	T->SetBranchAddress("Jet_id",Jet.id);
	T->SetBranchAddress("Jet_puId",Jet.puId);
    T->SetBranchAddress("Jet_qgl",Jet.qgl);

	T->SetBranchAddress("met_pt",&met_pt);
	T->SetBranchAddress("softActivityEWK_njets5",&Jet.EWKnsoft5);
	T->SetBranchAddress("softActivityEWK_HT",&Jet.EWKHTsoft);


	T->SetBranchAddress("nselLeptons",&nselLeptons);
	T->SetBranchAddress("selLeptons_pt",selLeptons_pt);
	T->SetBranchAddress("selLeptons_eta",selLeptons_eta);
	T->SetBranchAddress("selLeptons_phi",selLeptons_phi);
	T->SetBranchAddress("selLeptons_mass",selLeptons_mass);
	T->SetBranchAddress("selLeptons_charge",selLeptons_charge);
	T->SetBranchAddress("selLeptons_pdgId",selLeptons_pdgId);
	T->SetBranchAddress("selLeptons_looseIdPOG",selLeptons_looseIdPOG);
	T->SetBranchAddress("selLeptons_relIso04",selLeptons_relIso04);
	T->SetBranchAddress("selLeptons_relIso03",selLeptons_relIso03);
    T->SetBranchAddress("selLeptons_eleMVAIdSppring16GenPurp",selLeptons_eleMVAIdSppring16GenPurp); 
	T->SetBranchAddress("selLeptons_trackerLayers", selLeptons_trackerLayers);
	T->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	T->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	T->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	T->SetBranchAddress("selLeptons_eleSieie",selLeptons_eleSieie);
	T->SetBranchAddress("selLeptons_eleHoE",selLeptons_eleHoE);
	T->SetBranchAddress("selLeptons_eleDEta",selLeptons_eleDEta);
	T->SetBranchAddress("selLeptons_eleDPhi",selLeptons_eleDPhi);
	T->SetBranchAddress("selLeptons_eleEcalClusterIso",selLeptons_eleEcalClusterIso);
	T->SetBranchAddress("selLeptons_eleHcalClusterIso",selLeptons_eleHcalClusterIso);
	T->SetBranchAddress("selLeptons_dr03TkSumPt",selLeptons_dr03TkSumPt);



    TMVA::Reader *reader = new TMVA::Reader("Silent");
    float BDT_VBF = 2.;
    TBranch *branchBDT_VBF = T->Branch("BDT_VBF",&BDT_VBF,"BDT_VBF/F");





    float temp_softActivityEWK_njets5 = 0.;
    TMVAstruct TMVA;
    reader->AddVariable("Mqq",&TMVA.Mqq);
    reader->AddVariable("DeltaEtaQQ",&TMVA.DeltaEtaQQ);
    reader->AddVariable("q1_eta",&TMVA.q1_eta);
    reader->AddVariable("ll_pt",&TMVA.ll_pt);
    reader->AddVariable("ll_mass",&TMVA.ll_mass);
    reader->AddVariable("met_pt",&TMVA.met_pt);
    reader->AddVariable("EWKHTsoft",&TMVA.EWKHTsoft);
    reader->AddVariable("qq_pt",&TMVA.qq_pt);
    reader->AddVariable("RptHard",&TMVA.RptHard);
    reader->AddVariable("softActivityEWK_njets5",&temp_softActivityEWK_njets5);
//    reader->AddVariable("softActivityEWK_njets5",&TMVA.softActivityEWK_njets5);
    reader->AddVariable("btagCMVA",&TMVA.btagCSV);
    reader->AddVariable("btagCSV",&TMVA.btagCSV);
    reader->AddVariable("qgl_1q",&TMVA.qgl_1q);
    reader->AddVariable("qgl_2q",&TMVA.qgl_2q);
    reader->AddVariable("cosThetaStar",&TMVA.cosThetaStar);
    reader->AddVariable("cosThetaStarAbs",&TMVA.cosThetaStarAbs);
    reader->BookMVA("BDTG", "Classification_BDTG.weights.xml");







    Long64_t nentries = T->GetEntries();
    for (Long64_t i=0;i<nentries;i++) {
        if (i%1000 == 0) std::cout << "Writing " << i << "th output" << std::endl;
//        if (i > 4000) break;
        T->GetEntry(i);

        float maxBTagCSV = 0.;
        float maxBTagCMVA = 0.;
        BDT_VBF = -2.;



//////////////////Jets////////////////
		TLorentzVector Qjet1;
		TLorentzVector Qjet2;
		TLorentzVector qq;
		int good_jets = 0;
		vector<TLorentzVector> jets_pv;
		vector<int> jets_indices;
		float jet3_pt = 0;
		float jet3_eta;
        for (int j=0;j<nJets;j++){
            TLorentzVector jet0;
            if (!((Jet.id[j]>2)&&(Jet.puId[j]>0))) continue;
            
            jet0.SetPtEtaPhiM(Jet.pt[j],Jet.eta[j],Jet.phi[j],Jet.mass[j]);

            if (maxBTagCSV < Jet.btagCSV[j]) maxBTagCSV = Jet.btagCSV[j];
            if (maxBTagCMVA < Jet.btagCMVA[j]) maxBTagCMVA = Jet.btagCMVA[j];        
            
            bool condition = false;
            if (good_jets < 2) condition = true;
            else if (jet0.Eta() < max(jets_pv[0].Eta(), jets_pv[1].Eta()) && jet0.Eta() > min(jets_pv[0].Eta(), jets_pv[1].Eta())) condition = true;  
            if (condition) {
                jets_pv.push_back(jet0);
                jets_indices.push_back(j);
                good_jets++;
            }
        }
        if (good_jets<2) continue;
        
        
        Int_t EWKnsoft2_toSubtract = 0;
        Int_t EWKnsoft5_toSubtract = 0;
        Int_t EWKnsoft10_toSubtract = 0;


        Qjet1 = jets_pv[0];
        Qjet2 = jets_pv[1];
		if (good_jets>=3) {
			jet3_pt=jets_pv[2].Pt();
			jet3_eta=jets_pv[2].Eta();
		}
		qq=Qjet1+Qjet2;
		TMVA.Mqq = qq.M();
		Float_t qq_pt = qq.Pt();
		Float_t qqDeltaEta = TMath::Abs(Qjet1.Eta()-Qjet2.Eta());
		Float_t qqDeltaPhi = TMath::Abs(Qjet1.DeltaPhi(Qjet2));



////////////////////leptons////////////////
		TLorentzVector lepton1;
		TLorentzVector lepton2;
		TLorentzVector Zll;
		int idx_1stLepton = 0;
		int idx_2ndLepton = 0;
		int count_l=0;

		    for (int i_lep=0; i_lep<nselLeptons;i_lep++ ){
			    if (!((selLeptons_looseIdPOG[i_lep]>0) && (selLeptons_relIso04[i_lep]<0.25) && (TMath::Abs(selLeptons_pdgId[i_lep])==13 )) ) continue;
                if ((count_l==1) && (selLeptons_charge[idx_1stLepton]*selLeptons_charge[i_lep] > 0)) continue;
			    if (count_l==1) {
				    idx_2ndLepton=i_lep;
				    lepton2.SetPtEtaPhiM(selLeptons_pt[idx_2ndLepton], selLeptons_eta[idx_2ndLepton], selLeptons_phi[idx_2ndLepton], selLeptons_mass[idx_2ndLepton]);
				    count_l++;
				    break;
			    }
			    if (count_l==0) {
				    idx_1stLepton=i_lep;
				    lepton1.SetPtEtaPhiM(selLeptons_pt[idx_1stLepton], selLeptons_eta[idx_1stLepton], selLeptons_phi[idx_1stLepton], selLeptons_mass[idx_1stLepton]);
				    count_l++;
			    }
		    }

        if (count_l < 2) continue;

        Zll = lepton1 + lepton2;
		float Zll_ystar = Zll.Rapidity() - (Qjet1.Rapidity() + Qjet2.Rapidity()) ;
		Float_t DeltaRelQQ = (Qjet1+Qjet2).Pt()/( Qjet1.Pt()+Qjet2.Pt()); 

		TMVA.ll_zstar = TMath::Abs( Zll_ystar/ (Qjet1.Rapidity()-Qjet2.Rapidity() )) ;
		TMVA.DeltaEtaQQ = TMath::Abs(Qjet1.Eta()-Qjet2.Eta());
		TMVA.RptHard = (Qjet1+Qjet2+ Zll).Pt()/( Qjet1.Pt()+Qjet2.Pt() + Zll.Pt()); 



        TVector3 BoostVector_toMuonRestFrame = Zll.BoostVector();
        TLorentzVector positiveMuon_newSys;
        TLorentzVector positiveMuon_OldSys;
        if (selLeptons_charge[idx_1stLepton] > 0) positiveMuon_newSys = lepton1;
        else positiveMuon_newSys = lepton2;
        positiveMuon_OldSys = positiveMuon_newSys;
        positiveMuon_newSys.Boost(-BoostVector_toMuonRestFrame);
        TVector3 H_direction  = Zll.Vect();
        TVector3 mu_direction = positiveMuon_newSys.Vect();
        TMVA.cosThetaStar = H_direction.Dot(mu_direction)/H_direction.Mag()/mu_direction.Mag();
        TMVA.cosThetaStarAbs = TMath::Abs(TMVA.cosThetaStar);


        TMVA.btagCMVA = maxBTagCMVA;
        TMVA.btagCSV = maxBTagCSV;
        TMVA.softActivityEWK_njets5 = Jet.EWKnsoft5;// - EWKnsoft5_toSubtract;
        TMVA.qgl_1q = Jet.qgl[jets_indices[0]];
        TMVA.qgl_2q = Jet.qgl[jets_indices[1]];
        TMVA.q1_eta = jets_pv[0].Eta();
        TMVA.qq_pt = qq_pt;
        TMVA.ll_pt = Zll.Pt();
        TMVA.ll_mass = Zll.M();
        TMVA.met_pt = met_pt;
        TMVA.Jet3_pt = jet3_pt;
        TMVA.EWKHTsoft = Jet.EWKHTsoft;

        temp_softActivityEWK_njets5 = (float) TMVA.softActivityEWK_njets5;




        BDT_VBF = reader->EvaluateMVA("BDTG");



        if (BDT_VBF > 0.99) {
            std::cout << "Writing " << i << "th output" << std::endl;
            std::cout << "ll_mass: " << TMVA.ll_mass  << " \tll_pt: " << TMVA.ll_pt << " \tmet_pt: " << TMVA.met_pt << "  \tMqq: " << TMVA.Mqq << " \tDeltaEtaQQ: " << TMVA.DeltaEtaQQ << " \tRptHard: " << TMVA.RptHard  <<  std::endl;
            std::cout  << "EWKnsoft5: " << TMVA.softActivityEWK_njets5 << " \tjet3_pt: " << TMVA.Jet3_pt << " \tEWKHTsoft: " << TMVA.EWKHTsoft << " \tqq_pt: " << TMVA.qq_pt << " \tqgl_1q: " << TMVA.qgl_1q << " \tqgl_2q: " << TMVA.qgl_2q <<  std::endl;
            std::cout << "maxBTagCMVA: " << TMVA.btagCMVA << " \tmaxBTagCSV: " << TMVA.btagCSV  << " \tcosThetaStar: " << TMVA.cosThetaStar << " \tcosThetaStarAbs: " << TMVA.cosThetaStarAbs << " \tBDT_VBF: " << BDT_VBF <<  std::endl;
        }

//        BDT_VBF = 1.;
//        if (BDT_VBF > -1.5)
        branchBDT_VBF->Fill();


   }


   f->cd();

   T->Write();
   delete f;
}



void updatetree_withMVA() {

std::vector<std::string> file_names;
//file_names.push_back("../updateTreeTest/DY0JetsToLL_M-50");
//file_names.push_back("../updateTreeTest/DY1JetsToLL_M-50");

//file_names.push_back("DY0JetsToLL_M-50");
//file_names.push_back("DY1JetsToLL_M-50");
//file_names.push_back("DY2JetsToLL_M-50");
//file_names.push_back("DY3JetsToLL_M-50");
//file_names.push_back("DY4JetsToLL_M-50");
file_names.push_back("DYJetsToLL_M-50");

//file_names.push_back("DYJetstoLL_amc_0J");
//file_names.push_back("DYJetstoLL_amc_1J");
//file_names.push_back("DYJetstoLL_amc_2J");
//file_names.push_back("DYJetstoLL_madgraph");

//file_names.push_back("ST_s");
//file_names.push_back("ST_tW_antitop");
//file_names.push_back("ST_tW_top");
//file_names.push_back("ST_t_antitop");
//file_names.push_back("ST_t_top");

//file_names.push_back("TT");
//file_names.push_back("WW");
//file_names.push_back("WZ");
//file_names.push_back("ZZ");
//file_names.push_back("WJetsToLnu_madgraph");

//file_names.push_back("GluGlu_HToMuMu");
//file_names.push_back("VBF_HToMuMu");





for(int n=0; n<file_names.size(); ++n)
    updatetree_SingleFile(file_names[n]+"_v25_reskim.root");


//updatetree_SingleFile("SingleMuon_reminiaod_v25.root");


}







