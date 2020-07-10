// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/UtilAlgos/interface/MCMatchSelector.h"
#include "CommonTools/UtilAlgos/interface/MatchByDRDPt.h"
#include "CommonTools/UtilAlgos/interface/MatchLessByDPt.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace reco;
using namespace std;
using namespace edm;

class miniaod_ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit miniaod_ntuplizer(const edm::ParameterSet&);
      ~miniaod_ntuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      EDGetTokenT<GenParticleCollection> genToken_;
      EDGetTokenT<pat::MuonCollection> muonToken_;
      //EDGetTokenT<VertexCollection> vtxToken_;
      EDGetTokenT<TriggerResults> triggerBits_;
      //EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      //EDGetTokenT<TrackCollection> trackToken_;
      EDGetTokenT<pat::PackedCandidateCollection> pfToken_;

      // Declare the histograms
      TH1D *h_count;
      
      TTree *tr;
      float run, lumi, evt=-1;
      float genpt_1, geneta_1, genphi_1, genq_1=-1; // generated leading muon
      float genpt_2, geneta_2, genphi_2, genq_2=-1; // generated subleading muon
      float genpt_3, geneta_3, genphi_3, genq_3=-1; // generated pion 
      float pt_1, pterr_1, eta_1, phi_1, q_1, softID_1, looseID_1, mediumID_1, iso_1=-1; // leading muon
      float pt_2, pterr_2, eta_2, phi_2, q_2, softID_2, looseID_2, mediumID_2, iso_2=-1; // subleading muon
      float pt_3, eta_3, phi_3, q_3, iso_3=-1; // pion 
      float vtxprob_12, vtxprob_23, vtxprob_13, vtxprob_123=-1; // vertex fit probability
      int HLT_Mu7_IP4, HLT_Mu8_IP3, HLT_Mu8_IP5, HLT_Mu9_IP5, HLT_Mu9_IP6, HLT_Mu12_IP6=0; // B parking triggers
      int match_Mu7_IP4_1, match_Mu8_IP3_1, match_Mu8_IP5_1, match_Mu9_IP5_1, match_Mu9_IP6_1, match_Mu12_IP6_1=0;
      int match_Mu7_IP4_2, match_Mu8_IP3_2, match_Mu8_IP5_2, match_Mu9_IP5_2, match_Mu9_IP6_2, match_Mu12_IP6_2=0;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
miniaod_ntuplizer::miniaod_ntuplizer(const edm::ParameterSet& iConfig):
    genToken_(consumes<GenParticleCollection>(iConfig.getParameter<InputTag>("genParticles"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<InputTag>("muons"))),
    //vtxToken_(consumes<VertexCollection>(iConfig.getParameter<InputTag>("vertices"))),
    triggerBits_(consumes<TriggerResults>(iConfig.getParameter<InputTag>("bits"))),
    //triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<InputTag>("objects"))),
    //trackToken_(consumes<TrackCollection>(iConfig.getParameter<InputTag>("trks"))),
    pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
   //now do what ever initialization is needed
   //give the binning of the histograms
   usesResource("TFileService");
   edm::Service<TFileService> fs;

   h_count = fs->make<TH1D>("h_count" , "" , 10 , 0 , 10);

   tr = fs->make<TTree>("tr_tree", "tr_tree");

   tr->Branch("run", &run, "run/F");
   tr->Branch("lumi", &lumi, "lumi/F");
   tr->Branch("evt", &evt, "evt/F");
   tr->Branch("genpt_1", &genpt_1, "genpt_1/F");
   tr->Branch("geneta_1", &geneta_1, "geneta_1/F");
   tr->Branch("genphi_1", &genphi_1, "genphi_1/F");
   tr->Branch("genq_1", &genq_1, "genq_1/F");
   tr->Branch("genpt_2", &genpt_2, "genpt_2/F");
   tr->Branch("geneta_2", &geneta_2, "geneta_2/F");
   tr->Branch("genphi_2", &genphi_2, "genphi_2/F");
   tr->Branch("genq_2", &genq_2, "genq_2/F");
   tr->Branch("genpt_3", &genpt_3, "genpt_3/F");
   tr->Branch("geneta_3", &geneta_3, "geneta_3/F");
   tr->Branch("genphi_3", &genphi_3, "genphi_3/F");
   tr->Branch("genq_3", &genq_3, "genq_3/F");
   tr->Branch("pt_1", &pt_1, "pt_1/F");
   tr->Branch("pterr_1", &pterr_1, "pterr_1/F");
   tr->Branch("eta_1", &eta_1, "eta_1/F");
   tr->Branch("phi_1", &phi_1, "phi_1/F");
   tr->Branch("q_1", &q_1, "q_1/F");
   tr->Branch("iso_1", &iso_1, "iso_1/F");
   tr->Branch("softID_1", &softID_1, "softID_1/F");
   tr->Branch("looseID_1", &looseID_1, "looseID_1/F");
   tr->Branch("mediumID_1", &mediumID_1, "mediumID_1/F");
   tr->Branch("pt_2", &pt_2, "pt_2/F");
   tr->Branch("pterr_2", &pterr_2, "pterr_2/F");
   tr->Branch("eta_2", &eta_2, "eta_2/F");
   tr->Branch("phi_2", &phi_2, "phi_2/F");
   tr->Branch("q_2", &q_2, "q_2/F");
   tr->Branch("iso_2", &iso_2, "iso_2/F");
   tr->Branch("softID_2", &softID_2, "softID_2/F");
   tr->Branch("looseID_2", &looseID_2, "looseID_2/F");
   tr->Branch("mediumID_2", &mediumID_2, "mediumID_2/F");
   tr->Branch("pt_3", &pt_3, "pt_3/F");
   tr->Branch("eta_3", &eta_3, "eta_3/F");
   tr->Branch("phi_3", &phi_3, "phi_3/F");
   tr->Branch("q_3", &q_3, "q_3/F");
   tr->Branch("iso_3", &iso_3, "iso_3/F");
   tr->Branch("vtxprob_12", &vtxprob_12, "vtxprob_12/F");
   tr->Branch("vtxprob_23", &vtxprob_23, "vtxprob_23/F");
   tr->Branch("vtxprob_13", &vtxprob_13, "vtxprob_13/F");
   tr->Branch("vtxprob_123", &vtxprob_123, "vtxprob_123/F");
   tr->Branch("HLT_Mu7_IP4", &HLT_Mu7_IP4, "HLT_Mu7_IP4/I");
   tr->Branch("HLT_Mu8_IP3", &HLT_Mu8_IP3, "HLT_Mu8_IP3/I");
   tr->Branch("HLT_Mu8_IP5", &HLT_Mu8_IP5, "HLT_Mu8_IP5/I");
   tr->Branch("HLT_Mu9_IP5", &HLT_Mu9_IP5, "HLT_Mu9_IP5/I");
   tr->Branch("HLT_Mu9_IP6", &HLT_Mu9_IP6, "HLT_Mu9_IP6/I");
   tr->Branch("HLT_Mu12_IP6", &HLT_Mu12_IP6, "HLT_Mu12_IP6/I");
   tr->Branch("match_Mu7_IP4_1", &match_Mu7_IP4_1, "match_Mu7_IP4_1/I");
   tr->Branch("match_Mu8_IP3_1", &match_Mu8_IP3_1, "match_Mu8_IP3_1/I");
   tr->Branch("match_Mu8_IP5_1", &match_Mu8_IP5_1, "match_Mu8_IP5_1/I");
   tr->Branch("match_Mu9_IP5_1", &match_Mu9_IP5_1, "match_Mu9_IP5_1/I");
   tr->Branch("match_Mu9_IP6_1", &match_Mu9_IP6_1, "match_Mu9_IP6_1/I");
   tr->Branch("match_Mu12_IP6_1", &match_Mu12_IP6_1, "match_Mu12_IP6_1/I");
   tr->Branch("match_Mu7_IP4_2", &match_Mu7_IP4_2, "match_Mu7_IP4_2/I");
   tr->Branch("match_Mu8_IP3_2", &match_Mu8_IP3_2, "match_Mu8_IP3_2/I");
   tr->Branch("match_Mu8_IP5_2", &match_Mu8_IP5_2, "match_Mu8_IP5_2/I");
   tr->Branch("match_Mu9_IP5_2", &match_Mu9_IP5_2, "match_Mu9_IP5_2/I");
   tr->Branch("match_Mu9_IP6_2", &match_Mu9_IP6_2, "match_Mu9_IP6_2/I");
   tr->Branch("match_Mu12_IP6_2", &match_Mu12_IP6_2, "match_Mu12_IP6_2/I");
}


miniaod_ntuplizer::~miniaod_ntuplizer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
miniaod_ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   // Event counter
   //h_count->Fill(0.5,1.0);

   // Get all collections
   Handle<GenParticleCollection> genParticles;
   iEvent.getByToken(genToken_, genParticles);

   /*Handle<VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if (vertices->empty()) return; // skip the event if no PV found
   const reco::Vertex &PV = vertices->front();*/

   Handle<TriggerResults> triggerBits;
   //Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerBits_, triggerBits);
   //iEvent.getByToken(triggerObjects_, triggerObjects);

   /*Handle<TrackCollection> trks;
   iEvent.getByToken(trackToken_, trks);*/

   // #####################
   // Analyze gen particles
   // #####################
   TLorentzVector mygenmuon_1;
   mygenmuon_1.SetPtEtaPhiM(0,0,0,0);
   TLorentzVector mygenmuon_2;
   mygenmuon_2.SetPtEtaPhiM(0,0,0,0);
   TLorentzVector mygenpion;
   mygenpion.SetPtEtaPhiM(0,0,0,0);

   for(size_t i = 2; i < genParticles->size()-2; ++i) { // Loop over generated particles
     const GenParticle & p = (*genParticles)[i];
     const GenParticle & p2 = (*genParticles)[i+1];
     const GenParticle & p3 = (*genParticles)[i+2];
     const Candidate * mom = p.mother();
     if (fabs(p.pdgId())==13 and fabs(mom->pdgId())==521 and p.pdgId()==p2.pdgId()){ // find the generated muons
	  mygenmuon_1.SetPtEtaPhiM(p.pt(),p.eta(),p.phi(),p.mass());
	  mygenmuon_2.SetPtEtaPhiM(p2.pt(),p2.eta(),p2.phi(),p2.mass());
          mygenpion.SetPtEtaPhiM(p3.pt(),p3.eta(),p3.phi(),p3.mass());
          genq_1=p.charge();
          genq_2=p2.charge();
	  genq_3=p3.charge();
     }
     //std::cout<<p.pdgId()<<" ("<<mom->pdgId()<<"), ";
   }
   //std::cout<<std::endl;

   // Assign the generated variables
   genpt_1=mygenmuon_1.Pt(); geneta_1=mygenmuon_1.Eta(); genphi_1=mygenmuon_1.Phi();
   genpt_2=mygenmuon_2.Pt(); geneta_2=mygenmuon_2.Eta(); genphi_2=mygenmuon_2.Phi();
   genpt_3=mygenpion.Pt(); geneta_3=mygenpion.Eta(); genphi_3=mygenpion.Phi();

   // ######################
   // Get list of reco muons
   // ######################
   Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);

   vector<int> muon_indices;
   TLorentzVector myrecomu1;
   TLorentzVector myrecomu2;
   myrecomu1.SetPtEtaPhiM(0,0,0,0);
   myrecomu2.SetPtEtaPhiM(0,0,0,0);
   TLorentzVector mytmpmu;
   reco::TrackRef trk1;
   reco::TrackRef trk2;
   int k=0;
   for (const pat::Muon &mu : *muons) {
      mytmpmu.SetPtEtaPhiM(mu.pt(),mu.eta(),mu.phi(),0.105);
      if (mytmpmu.Pt()>2 and fabs(mytmpmu.Eta())<2.4 and mu.isGlobalMuon()) muon_indices.push_back(k);
      k++; 
   }

   // ######################
   // Get list of reco pions
   // ###################### 

   edm::Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);
   vector<int> pion_indices;
   reco::Track trk3;

   for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
        const pat::PackedCandidate &pf = (*pfs)[i];
	if (fabs(pf.pdgId()==211) and pf.pt()>2 and fabs(pf.eta())<2.4) pion_indices.push_back(i);
   }

   // ############################################
   // Fill events with at least 2 muons and 1 pion
   // ############################################

   if (muon_indices.size()>=2 and pion_indices.size()>=1){

      // ###########################################
      // Vertex fitting and choice of best 3 objects
      // ###########################################
      //int index_bestmu1, index_bestmu2, index_bestpion=-1;
      int index_bestmu1=muon_indices[0];
      int index_bestmu2=muon_indices[1];
      int index_bestpion=pion_indices[0];

      double fv_tC = -1;
      double fv_dOF = -1;
      double fv_Prob = -1;
      double old_prob=-1;

      for (unsigned int i_mu1=0; i_mu1<muon_indices.size()-1; ++i_mu1){
        for (unsigned int i_mu2=i_mu1+1; i_mu2<muon_indices.size(); ++i_mu2){
	  for (unsigned int i_pi=0; i_pi<pion_indices.size(); ++i_pi){
          //const pat::PackedCandidate &pf = (*pfs)[i];
            const pat::Muon &muon1 = (*muons)[muon_indices[i_mu1]];
            const pat::Muon &muon2 = (*muons)[muon_indices[i_mu2]];
	    const pat::PackedCandidate &pion = (*pfs)[pion_indices[i_pi]];
	    if (pion.hasTrackDetails() && !muon1.innerTrack().isNull() && !muon2.innerTrack().isNull()){ // all tracks exist and the objects are in the same area
              trk1 = muon1.innerTrack();
              trk2 = muon2.innerTrack();
              trk3 = pion.pseudoTrack();
              std::vector<reco::TransientTrack> t_trks;
              edm::ESHandle<TransientTrackBuilder> theB;
              iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
              t_trks.push_back(theB->build(trk1));
              t_trks.push_back(theB->build(trk2));
              t_trks.push_back(theB->build(trk3));
              KalmanVertexFitter kvf;
              TransientVertex fv = kvf.vertex(t_trks);
              if(fv.isValid()){
                fv_tC = fv.totalChiSquared();
                fv_dOF = fv.degreesOfFreedom();
                fv_Prob = TMath::Prob(fv_tC,(int)fv_dOF);
              }
	    }
	    else{ // bad candidate 
	      fv_Prob = 0;
	    }
            if (fv_Prob>old_prob) {
	       index_bestmu1=muon_indices[i_mu1]; 
	       index_bestmu2=muon_indices[i_mu2]; 
	       index_bestpion=pion_indices[i_pi];}
	       old_prob=fv_Prob;
          }
        }
      }

      // ##########################
      // Fill best object variables
      // ##########################

      run=iEvent.eventAuxiliary().run();
      lumi=iEvent.eventAuxiliary().luminosityBlock();
      evt=iEvent.eventAuxiliary().event();

      const pat::Muon &mymu1 = (*muons)[index_bestmu1];
      const pat::Muon &mymu2 = (*muons)[index_bestmu2];
      const pat::PackedCandidate &mypion = (*pfs)[index_bestpion];
      pt_1=mymu1.pt();
      pt_2=mymu2.pt();
      pterr_1=mymu1.muonBestTrack()->ptError();
      pterr_2=mymu2.muonBestTrack()->ptError();
      eta_1=mymu1.eta();
      eta_2=mymu2.eta();
      phi_1=mymu1.phi();
      phi_2=mymu2.phi();
      q_1=mymu1.charge();
      q_2=mymu2.charge();
      //softID_1=mymu1.isSoftMuon();
      //softID_2=mymu2.isSoftMuon();
      looseID_1=mymu1.isLooseMuon();
      looseID_2=mymu2.isLooseMuon();
      mediumID_1=mymu1.isMediumMuon();
      mediumID_2=mymu2.isMediumMuon();

      pt_3=mypion.pt();
      eta_3=mypion.eta();
      phi_3=mypion.phi();
      q_3=mypion.charge();

      // Compute vertex probabilities for best combination
      if (mypion.hasTrackDetails() && !mymu1.innerTrack().isNull() && !mymu2.innerTrack().isNull()){ // all tracks exist and the objects are in the same area
          trk1 = mymu1.innerTrack();
          trk2 = mymu2.innerTrack();
          trk3 = mypion.pseudoTrack();
          std::vector<reco::TransientTrack> t_trks;
          edm::ESHandle<TransientTrackBuilder> theB;
          iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
          t_trks.push_back(theB->build(trk1));
          t_trks.push_back(theB->build(trk2));
          t_trks.push_back(theB->build(trk3));
          KalmanVertexFitter kvf;
          TransientVertex fv = kvf.vertex(t_trks);
          if(fv.isValid()){
            fv_tC = fv.totalChiSquared();
            fv_dOF = fv.degreesOfFreedom();
            vtxprob_123 = TMath::Prob(fv_tC,(int)fv_dOF);
          }
      }
      if (mypion.hasTrackDetails() && !mymu1.innerTrack().isNull()){
          trk1 = mymu1.innerTrack();
          trk3 = mypion.pseudoTrack();
          std::vector<reco::TransientTrack> t_trks;
          edm::ESHandle<TransientTrackBuilder> theB;
          iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
          t_trks.push_back(theB->build(trk1));
          t_trks.push_back(theB->build(trk3));
          KalmanVertexFitter kvf;
          TransientVertex fv = kvf.vertex(t_trks);
          if(fv.isValid()){
            fv_tC = fv.totalChiSquared();
            fv_dOF = fv.degreesOfFreedom();
            vtxprob_13 = TMath::Prob(fv_tC,(int)fv_dOF);
          }
      }
      if (!mymu1.innerTrack().isNull() && !mymu2.innerTrack().isNull()){
          trk1 = mymu1.innerTrack();
          trk2 = mymu2.innerTrack();
          std::vector<reco::TransientTrack> t_trks;
          edm::ESHandle<TransientTrackBuilder> theB;
          iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
          t_trks.push_back(theB->build(trk1));
          t_trks.push_back(theB->build(trk2));
          KalmanVertexFitter kvf;
          TransientVertex fv = kvf.vertex(t_trks);
          if(fv.isValid()){
            fv_tC = fv.totalChiSquared();
            fv_dOF = fv.degreesOfFreedom();
            vtxprob_12 = TMath::Prob(fv_tC,(int)fv_dOF);
          }
      }
      if (mypion.hasTrackDetails() && !mymu2.innerTrack().isNull()){
          trk2 = mymu2.innerTrack();
          trk3 = mypion.pseudoTrack();
          std::vector<reco::TransientTrack> t_trks;
          edm::ESHandle<TransientTrackBuilder> theB;
          iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
          t_trks.push_back(theB->build(trk2));
          t_trks.push_back(theB->build(trk3));
          KalmanVertexFitter kvf;
          TransientVertex fv = kvf.vertex(t_trks);
          if(fv.isValid()){
            fv_tC = fv.totalChiSquared();
            fv_dOF = fv.degreesOfFreedom();
            vtxprob_23 = TMath::Prob(fv_tC,(int)fv_dOF);
          }
      }

      // ########
      // Triggers
      // ########

      const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
      for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
          if (names.triggerName(i).find("HLT_Mu7_IP4_part")<140) HLT_Mu7_IP4=triggerBits->accept(i);
          if (names.triggerName(i).find("HLT_Mu8_IP3_part")<140) HLT_Mu8_IP3=triggerBits->accept(i);
          if (names.triggerName(i).find("HLT_Mu8_IP5_part")<140) HLT_Mu8_IP5=triggerBits->accept(i);
          if (names.triggerName(i).find("HLT_Mu9_IP5_part")<140) HLT_Mu9_IP5=triggerBits->accept(i);
          if (names.triggerName(i).find("HLT_Mu9_IP6_part")<140) HLT_Mu9_IP6=triggerBits->accept(i);
          if (names.triggerName(i).find("HLT_Mu12_IP6_part")<140) HLT_Mu12_IP6=triggerBits->accept(i);
      }

      /*TLorentzVector mytrgobj;
      mytrgobj.SetPtEtaPhiM(0,0,0,0);
      TLorentzVector myrecomuon1;
      myrecomuon1.SetPtEtaPhiM(mymu1.pt(),mymu1.eta(),mymu1.phi(),mymu1.mass());
      TLorentzVector myrecomuon2;
      myrecomuon2.SetPtEtaPhiM(mymu2.pt(),mymu2.eta(),mymu2.phi(),mymu2.mass());

      for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
         obj.unpackPathNames(names);
         mytrgobj.SetPtEtaPhiM(obj.pt(),obj.eta(),obj.phi(),0);
         if (HLT_Mu7_IP4 && mytrgobj.DeltaR(myrecomuon1)<0.03 && mytrgobj.Pt()>0.9*myrecomuon1.Pt() && mytrgobj.Pt()<1.1*myrecomuon1.Pt() && obj.hasPathName( "HLT_Mu7_IP4_part", false, true )) match_Mu7_IP4_1=1; // FIXME need the exact name
         if (HLT_Mu7_IP4 && mytrgobj.DeltaR(myrecomuon2)<0.03 && mytrgobj.Pt()>0.9*myrecomuon2.Pt() && mytrgobj.Pt()<1.1*myrecomuon2.Pt() && obj.hasPathName( "HLT_Mu7_IP4_part", false, true )) match_Mu7_IP4_2=1; // FIXME need the exact name
         if (HLT_Mu8_IP3 && mytrgobj.DeltaR(myrecomuon1)<0.03 && mytrgobj.Pt()>0.9*myrecomuon1.Pt() && mytrgobj.Pt()<1.1*myrecomuon1.Pt() && obj.hasPathName( "HLT_Mu8_IP3_part", false, true )) match_Mu8_IP3_1=1; // FIXME need the exact name
         if (HLT_Mu8_IP3 && mytrgobj.DeltaR(myrecomuon2)<0.03 && mytrgobj.Pt()>0.9*myrecomuon2.Pt() && mytrgobj.Pt()<1.1*myrecomuon2.Pt() && obj.hasPathName( "HLT_Mu8_IP3_part", false, true )) match_Mu8_IP3_2=1; // FIXME need the exact name
         if (HLT_Mu8_IP5 && mytrgobj.DeltaR(myrecomuon1)<0.03 && mytrgobj.Pt()>0.9*myrecomuon1.Pt() && mytrgobj.Pt()<1.1*myrecomuon1.Pt() && obj.hasPathName( "HLT_Mu8_IP5_part", false, true )) match_Mu8_IP5_1=1; // FIXME need the exact name
         if (HLT_Mu8_IP5 && mytrgobj.DeltaR(myrecomuon2)<0.03 && mytrgobj.Pt()>0.9*myrecomuon2.Pt() && mytrgobj.Pt()<1.1*myrecomuon2.Pt() && obj.hasPathName( "HLT_Mu8_IP5_part", false, true )) match_Mu8_IP5_2=1; // FIXME need the exact name
         if (HLT_Mu9_IP5 && mytrgobj.DeltaR(myrecomuon1)<0.03 && mytrgobj.Pt()>0.9*myrecomuon1.Pt() && mytrgobj.Pt()<1.1*myrecomuon1.Pt() && obj.hasPathName( "HLT_Mu9_IP5_part", false, true )) match_Mu9_IP5_1=1; // FIXME need the exact name
         if (HLT_Mu9_IP5 && mytrgobj.DeltaR(myrecomuon2)<0.03 && mytrgobj.Pt()>0.9*myrecomuon2.Pt() && mytrgobj.Pt()<1.1*myrecomuon2.Pt() && obj.hasPathName( "HLT_Mu9_IP5_part", false, true )) match_Mu9_IP5_2=1; // FIXME need the exact name
         if (HLT_Mu9_IP6 && mytrgobj.DeltaR(myrecomuon1)<0.03 && mytrgobj.Pt()>0.9*myrecomuon1.Pt() && mytrgobj.Pt()<1.1*myrecomuon1.Pt() && obj.hasPathName( "HLT_Mu9_IP6_part", false, true )) match_Mu9_IP6_1=1; // FIXME need the exact name
         if (HLT_Mu9_IP6 && mytrgobj.DeltaR(myrecomuon2)<0.03 && mytrgobj.Pt()>0.9*myrecomuon2.Pt() && mytrgobj.Pt()<1.1*myrecomuon2.Pt() && obj.hasPathName( "HLT_Mu9_IP6_part", false, true )) match_Mu9_IP6_2=1; // FIXME need the exact name
         if (HLT_Mu12_IP6 && mytrgobj.DeltaR(myrecomuon1)<0.03 && mytrgobj.Pt()>0.9*myrecomuon1.Pt() && mytrgobj.Pt()<1.1*myrecomuon1.Pt() && obj.hasPathName( "HLT_Mu12_IP6_part", false, true )) match_Mu12_IP6_1=1; // FIXME need the exact name
         if (HLT_Mu12_IP6 && mytrgobj.DeltaR(myrecomuon2)<0.03 && mytrgobj.Pt()>0.9*myrecomuon2.Pt() && mytrgobj.Pt()<1.1*myrecomuon2.Pt() && obj.hasPathName( "HLT_Mu12_IP6_part", false, true )) match_Mu12_IP6_2=1; // FIXME need the exact name
      }*/

      tr->Fill();
   }
}

// ------------ method called once each job just before starting event loop  ------------
void 
miniaod_ntuplizer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
miniaod_ntuplizer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
miniaod_ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(miniaod_ntuplizer);
