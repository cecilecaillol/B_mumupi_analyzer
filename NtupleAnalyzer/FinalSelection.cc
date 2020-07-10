#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include "TMultiGraph.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <stdio.h>
#include <TF1.h>
#include <TDirectoryFile.h>
#include <TRandom3.h>
#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TKey.h"
#include "TTree.h"
#include "THashList.h"
#include "THStack.h"
#include "TPaveLabel.h"
#include "TFile.h"
#include "tr_Tree.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFunctor.h"

using namespace std;

#include <fstream>

int main(int argc, char** argv) {

    std::string input = *(argv + 1);
    std::string output = *(argv + 2);
    std::string sample = *(argv + 3);
    std::string name = *(argv + 4);

    TFile *f_Double = new TFile(input.c_str());
    TH1F* nbevt = (TH1F*) f_Double->Get("h_counter");
    float ngen = nbevt->GetBinContent(1);

    float xs=1.0; float weight=1.0; float luminosity=1.0;
    weight=luminosity*xs/ngen;


    //std::cout.setf(ios::fixed, ios::floatfield);
    //std::cout.precision(4);

    TTree *arbre = (TTree*) f_Double->Get("tr_tree");
    arbre->SetBranchAddress("run", &run);
    arbre->SetBranchAddress("lumi", &lumi);
    arbre->SetBranchAddress("evt", &evt);
    arbre->SetBranchAddress("genpt_1", &genpt_1);
    arbre->SetBranchAddress("geneta_1", &geneta_1);
    arbre->SetBranchAddress("genphi_1", &genphi_1);
    arbre->SetBranchAddress("genq_1", &genq_1);
    arbre->SetBranchAddress("genpt_2", &genpt_2);
    arbre->SetBranchAddress("geneta_2", &geneta_2);
    arbre->SetBranchAddress("genphi_2", &genphi_2);
    arbre->SetBranchAddress("genq_2", &genq_2);
    arbre->SetBranchAddress("genpt_3", &genpt_3);
    arbre->SetBranchAddress("geneta_3", &geneta_3);
    arbre->SetBranchAddress("genphi_3", &genphi_3);
    arbre->SetBranchAddress("genq_3", &genq_3);
    arbre->SetBranchAddress("pt_1", &pt_1);
    arbre->SetBranchAddress("pterr_1", &pterr_1);
    arbre->SetBranchAddress("eta_1", &eta_1);
    arbre->SetBranchAddress("phi_1", &phi_1);
    arbre->SetBranchAddress("q_1", &q_1);
    arbre->SetBranchAddress("looseID_1", &looseID_1);
    arbre->SetBranchAddress("mediumID_1", &mediumID_1);
    arbre->SetBranchAddress("pt_2", &pt_2);
    arbre->SetBranchAddress("pterr_2", &pterr_2);
    arbre->SetBranchAddress("eta_2", &eta_2);
    arbre->SetBranchAddress("phi_2", &phi_2);
    arbre->SetBranchAddress("q_2", &q_2);
    arbre->SetBranchAddress("looseID_2", &looseID_2);
    arbre->SetBranchAddress("mediumID_2", &mediumID_2);
    arbre->SetBranchAddress("pt_3", &pt_3);
    arbre->SetBranchAddress("eta_3", &eta_3);
    arbre->SetBranchAddress("phi_3", &phi_3);
    arbre->SetBranchAddress("q_3", &q_3);
    arbre->SetBranchAddress("vtxprob_12", &vtxprob_12);
    arbre->SetBranchAddress("vtxprob_23", &vtxprob_23);
    arbre->SetBranchAddress("vtxprob_13", &vtxprob_13);
    arbre->SetBranchAddress("vtxprob_123", &vtxprob_123);
    arbre->SetBranchAddress("HLT_Mu7_IP4", &HLT_Mu7_IP4);
    arbre->SetBranchAddress("HLT_Mu8_IP3", &HLT_Mu8_IP3);
    arbre->SetBranchAddress("HLT_Mu8_IP5", &HLT_Mu8_IP5);
    arbre->SetBranchAddress("HLT_Mu9_IP5", &HLT_Mu9_IP5);
    arbre->SetBranchAddress("HLT_Mu9_IP6", &HLT_Mu9_IP6);
    arbre->SetBranchAddress("HLT_Mu12_IP6", &HLT_Mu12_IP6);
    arbre->SetBranchAddress("m123", &m123);
    arbre->SetBranchAddress("validFraction_1", &validFraction_1);
    arbre->SetBranchAddress("normalizedChi2_1", &normalizedChi2_1);
    arbre->SetBranchAddress("chi2LocalPosition_1", &chi2LocalPosition_1);
    arbre->SetBranchAddress("trkKink_1", &trkKink_1);
    arbre->SetBranchAddress("segmentCompatibility_1", &segmentCompatibility_1);
    arbre->SetBranchAddress("softMvaValue_1", &softMvaValue_1);
    arbre->SetBranchAddress("validFraction_2", &validFraction_2);
    arbre->SetBranchAddress("normalizedChi2_2", &normalizedChi2_2);
    arbre->SetBranchAddress("chi2LocalPosition_2", &chi2LocalPosition_2);
    arbre->SetBranchAddress("trkKink_2", &trkKink_2);
    arbre->SetBranchAddress("segmentCompatibility_2", &segmentCompatibility_2);
    arbre->SetBranchAddress("softMvaValue_2", &softMvaValue_2);

   TH1F* cutflow=new TH1F("cutflow","cutflow",10,0,10);
   TH1F* h_m123_SS = new TH1F("h_m123_SS","h_m123_SS",100,0,100); h_m123_SS->Sumw2();
   TH1F* h_m123_OS = new TH1F("h_m123_OS","h_m123_OS",100,0,100); h_m123_OS->Sumw2();

   ofstream myfile_OS;
   myfile_OS.open (name+"_OS.txt");
   ofstream myfile_SS;
   myfile_SS.open (name+"_SS.txt");

   bool do_cutflow=true;

   Int_t nentries_wtn = (Int_t) arbre->GetEntries();
   for (Int_t i = 0; i < nentries_wtn; i++) {
        arbre->GetEntry(i);
        if (i % 100000 == 0) fprintf(stdout, "\r  Processed events: %8d of %8d ", i, nentries_wtn);
        fflush(stdout);

        TLorentzVector mymu1;
        mymu1.SetPtEtaPhiM(pt_1, eta_1, phi_1, 0.105);
        TLorentzVector mymu2;
        mymu2.SetPtEtaPhiM(pt_2, eta_2, phi_2, 0.105);
        TLorentzVector mypion;
        mypion.SetPtEtaPhiM(pt_3, eta_3, phi_3, 0.14);

	// Separation between objects and minimum pt/eta
        if (mymu1.DeltaR(mymu2)<0.02 or mymu1.DeltaR(mypion)<0.02 or mymu2.DeltaR(mypion)<0.02) continue;
	if (fabs(mymu1.Eta())<2.4 or fabs(mymu2.Eta())<2.4 or fabs(mypion.Eta())<2.4) continue;
	if (mymu1.Pt()<2 or mymu2.Pt()<2 or mypion.Pt()<2) continue;
	if (do_cutflow) cutflow->Fill(0.5);

	// Trigger
	if (!HLT_Mu7_IP4 && !HLT_Mu8_IP3 && !HLT_Mu8_IP5 && !HLT_Mu9_IP5 && !HLT_Mu9_IP6 && !HLT_Mu12_IP6) continue;
	if (do_cutflow) cutflow->Fill(1.5);

	// Muon ID
	if (!mediumID_1 or !mediumID_2) continue;
	if (do_cutflow) cutflow->Fill(2.5);

	// Not all the same charge
	if (fabs(q_1+q_2+q_3)==3) continue;
        if (do_cutflow) cutflow->Fill(3.5);
	bool is_OS=(q_1*q_2<0);
        bool is_SS=(q_1*q_2>0);

	// Vertex probability
	if (vtxprob_123<0.01) continue;
	if (do_cutflow) cutflow->Fill(4.5);
	if (do_cutflow && is_SS) cutflow->Fill(5.5);

	// Fill 
	float mass_123=(mymu1+mymu2+mypion).M();
	if (is_SS){
	   h_m123_SS->Fill(mass_123,weight);
           myfile_SS<<mass_123<<endl;
	}
        if (is_OS){
           h_m123_OS->Fill(mass_123,weight);
           myfile_OS<<mass_123<<endl;
       	}
				
    } // end of loop over events

    TFile *fout = TFile::Open(output.c_str(), "RECREATE");
    fout->cd();
    h_m123_OS->Write();
    cutflow->Write();
    h_m123_SS->Write();

    fout->Close();

    myfile_OS.close();
    myfile_SS.close();

} 


