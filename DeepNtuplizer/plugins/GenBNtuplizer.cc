// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include <boost/core/demangle.hpp>

//ROOT includes
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include "TBranch.h"
#include <string>
#include <vector>


//CMSSW includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

// for ivf
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"

#include "DataFormats/BTauReco/interface/PixelClusterTagInfo.h"

#include <algorithm>
#include <iterator>
#include <map>

//trash?

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#if defined( __GXX_EXPERIMENTAL_CXX0X__)
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#endif

//struct MagneticField;


class GenBNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit GenBNtuplizer(const edm::ParameterSet&);
  ~GenBNtuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  edm::EDGetTokenT<edm::View<pat::Jet> >      jetToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> >      qjetToken_;
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  bool WriteJetPart=false;
  bool WritePFcands=false;


  TTree *tree_;
  std::vector<float> m_qjet_pt;
  std::vector<float> m_qjet_eta;
  std::vector<float> m_qjet_phi;
  std::vector<float> m_qjet_mass;
  std::vector<float> m_qjet_energy;
  std::vector<float> m_qjet_et;
  std::vector<float> m_qjet_ntrk;
  std::vector<float> m_qjet_vx;
  std::vector<float> m_qjet_vy;
  std::vector<float> m_qjet_vz;

  std::vector<float> m_jet_pt;
  std::vector<float> m_jet_eta;
  std::vector<float> m_jet_phi;
  std::vector<float> m_jet_mass;
  std::vector<float> m_jet_energy;
  std::vector<float> m_jet_et;
  std::vector<float> m_jet_ntrk;
  std::vector<float> m_jet_npf;
  std::vector<float> m_jet_vx;
  std::vector<float> m_jet_vy;
  std::vector<float> m_jet_vz;
  std::vector<int> m_jetpart_jetIdx;
  std::vector<float> m_jetpart_pt;
  std::vector<float> m_jetpart_eta;
  std::vector<float> m_jetpart_phi;
  std::vector<float> m_jetpart_charge;
  std::vector<int> m_jetpart_pdgId;


  std::vector<int> m_pf_pdgId;
  std::vector<float> m_pf_pt;
  std::vector<float> m_pf_eta;
  std::vector<float> m_pf_phi;
  std::vector<float> m_pf_charge;

  std::vector<float> m_sv_pt;
  std::vector<float> m_sv_eta;
  std::vector<float> m_sv_phi;
  std::vector<float> m_sv_mass;
  std::vector<float> m_sv_energy;
  std::vector<float> m_sv_et;
  std::vector<float> m_sv_ntrk;
  std::vector<float> m_sv_vx;
  std::vector<float> m_sv_vy;
  std::vector<float> m_sv_vz;


  std::vector<float> m_genquark_pdgId;
  std::vector<float> m_genquark_pt;
  std::vector<float> m_genquark_eta;
  std::vector<float> m_genquark_phi;
  std::vector<float> m_genquark_status;
  std::vector<float> m_genquark_mother;
  std::vector<float> m_genquark_vx;
  std::vector<float> m_genquark_vy;
  std::vector<float> m_genquark_vz;
  std::vector<float> m_genquark_charge;

  std::vector<float> m_genmeson_pdgId;
  std::vector<float> m_genmeson_pt;
  std::vector<float> m_genmeson_eta;
  std::vector<float> m_genmeson_phi;
  std::vector<float> m_genmeson_status;
  std::vector<float> m_genmeson_vx;
  std::vector<float> m_genmeson_vy;
  std::vector<float> m_genmeson_vz;
  std::vector<float> m_genmeson_nfinals;
  std::vector<float> m_genmeson_ncharge;
  std::vector<float> m_genmeson_frcharge;
  std::vector<float> m_genmeson_charge;
  std::vector<float> m_genmeson_jet1;
  std::vector<float> m_genmeson_jet1_dr;
  std::vector<float> m_genmeson_jet2;
  std::vector<float> m_genmeson_jet2_dr;
  std::vector<float> m_genmeson_jet3;
  std::vector<float> m_genmeson_jet3_dr;
  std::vector<float> m_genmeson_qjet1;
  std::vector<float> m_genmeson_qjet1_dr;
  std::vector<float> m_genmeson_qjet2;
  std::vector<float> m_genmeson_qjet2_dr;
  std::vector<float> m_genmeson_qjet3;
  std::vector<float> m_genmeson_qjet3_dr;
  std::vector<float> m_genmeson_sv1;
  std::vector<float> m_genmeson_sv1_dr;
  std::vector<float> m_genmeson_sv2;
  std::vector<float> m_genmeson_sv2_dr;
  std::vector<float> m_genmeson_sv3;
  std::vector<float> m_genmeson_sv3_dr;

  std::vector<float> m_gendaughters_pdgId;
  std::vector<float> m_gendaughters_mesonIdx;

  std::vector<float> m_genpart_pdgId;
  std::vector<float> m_genpart_pt;
  std::vector<float> m_genpart_eta;
  std::vector<float> m_genpart_phi;
  std::vector<float> m_genpart_status;
  std::vector<float> m_genpart_vx;
  std::vector<float> m_genpart_vy;
  std::vector<float> m_genpart_vz;
  std::vector<float> m_genpart_mother;
  std::vector<float> m_genpart_Bmeson;
  std::vector<float> m_genpart_Bmeson_pt;
  std::vector<float> m_genpart_Bmeson_eta;
  std::vector<float> m_genpart_Bmeson_phi;
  std::vector<float> m_genpart_BmesonIdx;


  int m_qjets=0;
  int m_genquarks=0;
  int m_jets=0;  
  int m_svs=0;
};


GenBNtuplizer::GenBNtuplizer(const edm::ParameterSet& iConfig):
  jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
  qjetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("pfChargedJets"))),
  svToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secVertices"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparts"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfcands"))),
  WriteJetPart((iConfig.getParameter<bool>("writeJetPart"))),
  WritePFcands((iConfig.getParameter<bool>("writePFcands")))
{


}


GenBNtuplizer::~GenBNtuplizer()
{
}


// ------------ method called for each event  ------------
void
GenBNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  m_qjet_pt.clear();
  m_qjet_eta.clear();
  m_qjet_phi.clear();
  m_qjet_mass.clear();
  m_qjet_energy.clear();
  m_qjet_et.clear();
  m_qjet_vx.clear();
  m_qjet_vy.clear();
  m_qjet_vz.clear();
  m_qjet_ntrk.clear();

  m_jet_pt.clear();
  m_jet_eta.clear();
  m_jet_phi.clear();
  m_jet_mass.clear();
  m_jet_energy.clear();
  m_jet_et.clear();
  m_jet_vx.clear();
  m_jet_vy.clear();
  m_jet_vz.clear();
  m_jet_ntrk.clear();
  m_jet_npf.clear();
  m_jetpart_jetIdx.clear();
  m_jetpart_pt.clear();
  m_jetpart_eta.clear();
  m_jetpart_phi.clear();
  m_jetpart_charge.clear();
  m_jetpart_pdgId.clear();


  m_pf_pdgId.clear();
  m_pf_pt.clear();
  m_pf_eta.clear();
  m_pf_phi.clear();
  m_pf_charge.clear();
  

  m_sv_pt.clear();
  m_sv_eta.clear();
  m_sv_phi.clear();
  m_sv_mass.clear();
  m_sv_energy.clear();
  m_sv_et.clear();
  m_sv_vx.clear();
  m_sv_vy.clear();
  m_sv_vz.clear();
  m_sv_ntrk.clear();

  m_genquark_pdgId.clear();;
  m_genquark_pt.clear();;
  m_genquark_eta.clear();
  m_genquark_phi.clear();
  m_genquark_status.clear();
  m_genquark_mother.clear();
  m_genquark_vx.clear();
  m_genquark_vy.clear();
  m_genquark_vz.clear();
  m_genquark_charge.clear();

  m_genmeson_pdgId.clear();
  m_genmeson_pt.clear();
  m_genmeson_eta.clear();
  m_genmeson_phi.clear();
  m_genmeson_status.clear();
  m_genmeson_vx.clear();
  m_genmeson_vy.clear();
  m_genmeson_vz.clear();
  m_genmeson_nfinals.clear();
  m_genmeson_ncharge.clear();
  m_genmeson_frcharge.clear();
  m_genmeson_charge.clear();
  m_genmeson_qjet1.clear();
  m_genmeson_qjet1_dr.clear();
  m_genmeson_qjet2.clear();
  m_genmeson_qjet2_dr.clear();
  m_genmeson_qjet3.clear();
  m_genmeson_qjet3_dr.clear();
  m_genmeson_jet1.clear();
  m_genmeson_jet1_dr.clear();
  m_genmeson_jet2.clear();
  m_genmeson_jet2_dr.clear();
  m_genmeson_jet3.clear();
  m_genmeson_jet3_dr.clear();
  m_genmeson_sv1.clear();
  m_genmeson_sv1_dr.clear();
  m_genmeson_sv2.clear();
  m_genmeson_sv2_dr.clear();
  m_genmeson_sv3.clear();
  m_genmeson_sv3_dr.clear();

  m_gendaughters_pdgId.clear();
  m_gendaughters_mesonIdx.clear();

  m_genpart_pdgId.clear();
  m_genpart_pt.clear();
  m_genpart_eta.clear();
  m_genpart_phi.clear();
  m_genpart_mother.clear();
  m_genpart_Bmeson.clear();
  m_genpart_Bmeson_pt.clear();
  m_genpart_Bmeson_eta.clear();
  m_genpart_Bmeson_phi.clear();
  m_genpart_BmesonIdx.clear();


  m_qjets=0;
  m_genquarks=0;
  m_jets=0;
  m_svs=0;


  //global info

  edm::Handle<std::vector<reco::VertexCompositePtrCandidate> > secvertices;
  iEvent.getByToken(svToken_, secvertices);

  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetToken_, jets);

  edm::Handle<edm::View<pat::Jet> > qjets;
  iEvent.getByToken(qjetToken_, qjets);

  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfToken_, pfcands);
  //std::cout<<std::endl<<"here"<<std::endl;

   int igen = -1;
   for (auto Bmeson_iter = genParticlesHandle->begin(); Bmeson_iter != genParticlesHandle->end(); ++Bmeson_iter) {
     if ( fabs(Bmeson_iter->pdgId())<500 || fabs(Bmeson_iter->pdgId())>600) continue;
     const reco::GenParticle* mother =  static_cast< const reco::GenParticle*> (Bmeson_iter->mother());
     
     if (mother==NULL) continue; 
     if ( Bmeson_iter->numberOfDaughters()==0) continue;

//     const reco::GenParticle*  daughter =static_cast< const reco::GenParticle*> (Bmeson_iter->daughter(0));

     //std::cout<<" gen "<<Bmeson_iter->pdgId()<<" pt "<<Bmeson_iter->pt()<<" eta "<<Bmeson_iter->eta()<<" phi "<<Bmeson_iter->phi()<<" mother "<<mother->pdgId()<<std::endl;

     bool Skip=false;
     for (unsigned int ib=0; ib<Bmeson_iter->numberOfDaughters(); ib++){
      // std::cout<<" --daughter "<<ib<<" pdg  "<<Bmeson_iter->daughter(ib)->pdgId()<<" pt "<<Bmeson_iter->daughter(ib)->pt()<<" eta "<<Bmeson_iter->daughter(ib)->eta()<<" phi "<<Bmeson_iter->daughter(ib)->phi()<<std::endl;
       if ( fabs(Bmeson_iter->daughter(ib)->pdgId())>500 && fabs(Bmeson_iter->daughter(ib)->pdgId())<600)
         Skip=true;       
     }
     if (Skip) {
    //   std::cout<<"skip1"<<std::endl;
       continue;
     }

     const reco::GenParticle* quark_iter =  static_cast< const reco::GenParticle*> (Bmeson_iter->mother());

     bool Read=true; 
     if (fabs(quark_iter->pdgId())==5 && fabs( (quark_iter->mother())->pdgId() )>9000000)  Read=false;
    // if (Read) continue;
     Skip=false;
     while (Read && !Skip){
        quark_iter = static_cast<const reco::GenParticle*> (quark_iter->mother());
        
        if (quark_iter==NULL){ Read=false; Skip=true; break;}
        if (quark_iter->mother()==NULL ){ Read=false; Skip=true; break;}
        if (fabs(quark_iter->pdgId())==5 && fabs( (quark_iter->mother())->pdgId() )>9000000)  Read=false;
     }
     if (Skip){
  //      std::cout<<"skip2"<<std::endl;
        continue;
     }
     igen++;
    // std::cout<<" quark "<<quark_iter->pdgId()<<"  "<< (quark_iter->mother())->pdgId()<<std::endl;
     m_genmeson_pdgId.push_back(Bmeson_iter->pdgId());
     m_genmeson_pt.push_back(Bmeson_iter->pt());
     m_genmeson_eta.push_back(Bmeson_iter->eta());
     m_genmeson_phi.push_back(Bmeson_iter->phi());
     m_genmeson_status.push_back(Bmeson_iter->status());
     m_genmeson_vx.push_back(Bmeson_iter->vx());
     m_genmeson_vy.push_back(Bmeson_iter->vy());
     m_genmeson_vz.push_back(Bmeson_iter->vz());
     m_genmeson_nfinals.push_back(Bmeson_iter->numberOfDaughters());
     m_genmeson_charge.push_back(Bmeson_iter->charge());

     float ncharge=0;
     float frcharge=0;
     for (unsigned int ib=0; ib<Bmeson_iter->numberOfDaughters(); ib++){
         ncharge += fabs(Bmeson_iter->daughter(ib)->charge());
         frcharge += fabs(Bmeson_iter->daughter(ib)->charge())*Bmeson_iter->daughter(ib)->et();
         m_gendaughters_pdgId.push_back(Bmeson_iter->daughter(ib)->pdgId());
         m_gendaughters_mesonIdx.push_back(igen);
     }
     m_genmeson_ncharge.push_back(ncharge);
     m_genmeson_frcharge.push_back(frcharge/Bmeson_iter->et());
    
     m_genquark_pdgId.push_back(quark_iter->pdgId());
     m_genquark_pt.push_back(quark_iter->pt());
     m_genquark_eta.push_back(quark_iter->eta());
     m_genquark_phi.push_back(quark_iter->phi());
     m_genquark_status.push_back(quark_iter->status());
     m_genquark_vx.push_back(quark_iter->vx());
     m_genquark_vy.push_back(quark_iter->vy());
     m_genquark_vz.push_back(quark_iter->vz());
     m_genquark_mother.push_back(quark_iter->mother()->pdgId());
     m_genquark_charge.push_back(quark_iter->charge());
      
     // match stuff
     // trackjet
     std::vector<std::pair<float,int>> qjets_dr_idx;
     int nqjet=-1;
     for(const pat::Jet& qjet: *qjets){
        nqjet+=1;
        if (1.0< TMath::Sqrt(reco::deltaR2(qjet.eta(),qjet.phi(),Bmeson_iter->eta(),Bmeson_iter->phi()) ) )
           continue;
        qjets_dr_idx.push_back( std::make_pair(TMath::Sqrt(reco::deltaR2(qjet.eta(),qjet.phi(),Bmeson_iter->eta(),Bmeson_iter->phi())),nqjet)  );
     }

     std::sort(qjets_dr_idx.begin(), qjets_dr_idx.end(), [](auto &left, auto &right) {   return left.first < right.first;  });

     if (qjets_dr_idx.size()>0){
        m_genmeson_qjet1.push_back(qjets_dr_idx[0].second);
        m_genmeson_qjet1_dr.push_back(qjets_dr_idx[0].first);
     } else {
        m_genmeson_qjet1.push_back(-1);
        m_genmeson_qjet1_dr.push_back(100);
     }
     if (qjets_dr_idx.size()>1){
        m_genmeson_qjet2.push_back(qjets_dr_idx[1].second);
        m_genmeson_qjet2_dr.push_back(qjets_dr_idx[1].first);
     } else{
        m_genmeson_qjet2.push_back(-1);
        m_genmeson_qjet2_dr.push_back(100);
     }
     if (qjets_dr_idx.size()>2){
        m_genmeson_qjet3.push_back(qjets_dr_idx[2].second);
        m_genmeson_qjet3_dr.push_back(qjets_dr_idx[2].first);
     } else {
        m_genmeson_qjet3.push_back(-1);
        m_genmeson_qjet3_dr.push_back(100);
     }

     // std jet
     std::vector<std::pair<float,int>> jets_dr_idx;
     int njet=-1;
     for(const pat::Jet& jet: *jets){
        njet+=1;
        if (1.0< TMath::Sqrt(reco::deltaR2(jet.eta(),jet.phi(),Bmeson_iter->eta(),Bmeson_iter->phi()) ) )
           continue;
        jets_dr_idx.push_back( std::make_pair(TMath::Sqrt(reco::deltaR2(jet.eta(),jet.phi(),Bmeson_iter->eta(),Bmeson_iter->phi())),njet)  );
     }

     std::sort(jets_dr_idx.begin(), jets_dr_idx.end(), [](auto &left, auto &right) {   return left.first < right.first;  });

     if (jets_dr_idx.size()>0){
        m_genmeson_jet1.push_back(jets_dr_idx[0].second);
        m_genmeson_jet1_dr.push_back(jets_dr_idx[0].first);
     } else {
        m_genmeson_jet1.push_back(-1);
        m_genmeson_jet1_dr.push_back(100);
     }
     if (jets_dr_idx.size()>1){
        m_genmeson_jet2.push_back(jets_dr_idx[1].second);
        m_genmeson_jet2_dr.push_back(jets_dr_idx[1].first);
     } else{
        m_genmeson_jet2.push_back(-1);
        m_genmeson_jet2_dr.push_back(100);
     }
     if (jets_dr_idx.size()>2){
        m_genmeson_jet3.push_back(jets_dr_idx[2].second);
        m_genmeson_jet3_dr.push_back(jets_dr_idx[2].first);
     } else {
        m_genmeson_jet3.push_back(-1);
        m_genmeson_jet3_dr.push_back(100);
     }

     std::vector<std::pair<float,int>> sv_dr_idx;
     int nsv=-1;
     for (const reco::VertexCompositePtrCandidate &sv : *secvertices){
        nsv+=1;
        if (1.0< TMath::Sqrt(reco::deltaR2(sv.eta(),sv.phi(),Bmeson_iter->eta(),Bmeson_iter->phi()) ) )
           continue;
        sv_dr_idx.push_back( std::make_pair(TMath::Sqrt(reco::deltaR2(sv.eta(),sv.phi(),Bmeson_iter->eta(),Bmeson_iter->phi())),nsv)  );
     }

     std::sort(sv_dr_idx.begin(), sv_dr_idx.end(), [](auto &left, auto &right) {   return left.first < right.first;  });

     if (sv_dr_idx.size()>0){
        m_genmeson_sv1.push_back(sv_dr_idx[0].second);
        m_genmeson_sv1_dr.push_back(sv_dr_idx[0].first);
     } else {
        m_genmeson_sv1.push_back(-1);
        m_genmeson_sv1_dr.push_back(100);
     }
     if (sv_dr_idx.size()>1){
        m_genmeson_sv2.push_back(sv_dr_idx[1].second);
        m_genmeson_sv2_dr.push_back(sv_dr_idx[1].first);
     } else{
        m_genmeson_sv2.push_back(-1);
        m_genmeson_sv2_dr.push_back(100);
     }
     if (sv_dr_idx.size()>2){
        m_genmeson_sv3.push_back(sv_dr_idx[2].second);
        m_genmeson_sv3_dr.push_back(sv_dr_idx[2].first);
     } else {
        m_genmeson_sv3.push_back(-1);
        m_genmeson_sv3_dr.push_back(100);
     }
   
   }	   


  // genpart   
  for (auto gens_iter = genParticlesHandle->begin(); gens_iter != genParticlesHandle->end(); ++gens_iter) {
     if ( (11!=abs(gens_iter->pdgId()) and abs(gens_iter->pdgId())!=13 
           and abs(gens_iter->pdgId())!=15) and
          (abs(gens_iter->pdgId())!=211 and abs(gens_iter->pdgId())!=321)   
        ) 
           continue;
    // std::cout<<" reverse pdg "<<gens_iter->pdgId()<<std::endl;
     const reco::GenParticle* mother =  static_cast<const reco::GenParticle*> (gens_iter->mother());
     if (mother==NULL) 
        continue;
     float tmp_mother_pdg = mother->pdgId();
    // std::cout<<" mother1 "<<mother->pdgId()<<std::endl;
     bool FromB=false; 
     while (!FromB) {
        mother =  static_cast< const reco::GenParticle*> (mother->mother());
        if ( mother==NULL) break;
      //std::cout<<" motherX "<<mother->pdgId()<<std::endl;
        if (abs(mother->pdgId())>500 and abs(mother->pdgId())<600){
           FromB=true;
          // std::cout<<"  "<<std::distance(genParticlesHandle->begin(), it); <<std::endl;
        }
    
     }
     if (!FromB) continue;
     m_genpart_pdgId.push_back(gens_iter->pdgId());
     m_genpart_pt.push_back(gens_iter->pt());
     m_genpart_eta.push_back(gens_iter->eta());
     m_genpart_phi.push_back(gens_iter->phi());
     m_genpart_mother.push_back(tmp_mother_pdg);
     m_genpart_Bmeson.push_back(mother->pdgId());
     m_genpart_Bmeson_pt.push_back(mother->pt());
     m_genpart_Bmeson_eta.push_back(mother->eta());
     m_genpart_Bmeson_phi.push_back(mother->phi());
     //m_genpart_BmesonIdx.clear();

     //std::cout<<"saved"<<std::endl;
  }


  // loop over the jets
  for(const pat::Jet& jet: *qjets){
    m_qjets+=1;
    m_qjet_pt.push_back(jet.pt());
    m_qjet_eta.push_back(jet.eta());
    m_qjet_phi.push_back(jet.phi());
    m_qjet_mass.push_back(jet.mass());
    m_qjet_energy.push_back(jet.energy());
    m_qjet_et.push_back(jet.et());
    m_qjet_vx.push_back(jet.vx());
    m_qjet_vy.push_back(jet.vy());
    m_qjet_vz.push_back(jet.vz());
    m_qjet_ntrk.push_back(jet.chargedMultiplicity());
  } // end of looping over the jets

  // loop over the jets
  for(const pat::Jet& jet: *jets){
    m_jet_pt.push_back(jet.pt());
    m_jet_eta.push_back(jet.eta());
    m_jet_phi.push_back(jet.phi());
    m_jet_mass.push_back(jet.mass());
    m_jet_energy.push_back(jet.energy());
    m_jet_et.push_back(jet.et());
    m_jet_vx.push_back(jet.vx());
    m_jet_vy.push_back(jet.vy());
    m_jet_vz.push_back(jet.vz());
    m_jet_ntrk.push_back(jet.chargedMultiplicity() );
    m_jet_npf.push_back(jet.numberOfDaughters() );
    if (WriteJetPart){
      for (unsigned int idx=0; idx< jet.numberOfDaughters(); idx++){
          m_jetpart_pt.push_back(jet.daughter(idx)->pt());
          m_jetpart_eta.push_back(jet.daughter(idx)->eta());
          m_jetpart_phi.push_back(jet.daughter(idx)->phi());
          m_jetpart_charge.push_back(jet.daughter(idx)->charge());
          m_jetpart_jetIdx.push_back(m_jets);
          m_jetpart_pdgId.push_back(jet.daughter(idx)->pdgId());
      }
    }
    m_jets+=1;
//    std::vector< reco::PFCandidatePtr > pfs = jet.getPFConstituents () ;
  } // end of looping over the jets


  // loop over the sv
  for (const reco::VertexCompositePtrCandidate &sv : *secvertices){
    m_svs+=1;
    m_sv_pt.push_back(sv.pt());
    m_sv_eta.push_back(sv.eta());
    m_sv_phi.push_back(sv.phi());
    m_sv_mass.push_back(sv.mass());
    m_sv_energy.push_back(sv.energy());
    m_sv_et.push_back(sv.et());
    m_sv_vx.push_back(sv.vx());
    m_sv_vy.push_back(sv.vy());
    m_sv_vz.push_back(sv.vz());
    m_sv_ntrk.push_back(sv.numberOfDaughters());
  } // end of looping over the jets


  if (WritePFcands){
    for (auto & pfcand: *pfcands){
       m_pf_pt.push_back(pfcand.pt());
       m_pf_eta.push_back(pfcand.eta());
       m_pf_phi.push_back(pfcand.phi());
       m_pf_charge.push_back(pfcand.charge());
       m_pf_pdgId.push_back(pfcand.pdgId());
    }
 }

  tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void
GenBNtuplizer::beginJob()
{
   edm::Service<TFileService> fs;

   tree_=(fs->make<TTree>("tree" ,"tree" ));
   tree_->Branch("nqjet", &m_qjets);
   tree_->Branch("qjet_pt", &m_qjet_pt);
   tree_->Branch("qjet_eta", &m_qjet_eta);
   tree_->Branch("qjet_phi", &m_qjet_phi);
   tree_->Branch("qjet_mass", &m_qjet_mass);
   tree_->Branch("qjet_energy", &m_qjet_energy);
   tree_->Branch("qjet_et", &m_qjet_et);
   tree_->Branch("qjet_vx", &m_qjet_vx);
   tree_->Branch("qjet_vy", &m_qjet_vy);
   tree_->Branch("qjet_vz", &m_qjet_vz);
   tree_->Branch("qjet_ntrk", &m_qjet_ntrk);

   tree_->Branch("njet", &m_jets);
   tree_->Branch("jet_pt", &m_jet_pt);
   tree_->Branch("jet_eta", &m_jet_eta);
   tree_->Branch("jet_phi", &m_jet_phi);
   tree_->Branch("jet_mass", &m_jet_mass);
   tree_->Branch("jet_energy", &m_jet_energy);
   tree_->Branch("jet_et", &m_jet_et);
   tree_->Branch("jet_vx", &m_jet_vx);
   tree_->Branch("jet_vy", &m_jet_vy);
   tree_->Branch("jet_vz", &m_jet_vz);
   tree_->Branch("jet_ntrk", &m_jet_ntrk);
   tree_->Branch("jet_npf", &m_jet_npf);
   tree_->Branch("jetpart_jetIdx", &m_jetpart_jetIdx);

   tree_->Branch("jetpart_pt", &m_jetpart_pt);
   tree_->Branch("jetpart_eta", &m_jetpart_eta);
   tree_->Branch("jetpart_phi", &m_jetpart_phi);
   tree_->Branch("jetpart_charge", &m_jetpart_charge);
   tree_->Branch("jetpart_pdgId", &m_jetpart_pdgId);

   tree_->Branch("nsv", &m_svs);
   tree_->Branch("sv_pt", &m_sv_pt);
   tree_->Branch("sv_eta", &m_sv_eta);
   tree_->Branch("sv_phi", &m_sv_phi);
   tree_->Branch("sv_mass", &m_sv_mass);
   tree_->Branch("sv_energy", &m_sv_energy);
   tree_->Branch("sv_et", &m_sv_et);
   tree_->Branch("sv_vx", &m_sv_vx);
   tree_->Branch("sv_vy", &m_sv_vy);
   tree_->Branch("sv_vz", &m_sv_vz);
   tree_->Branch("sv_ntrk", &m_sv_ntrk);

   tree_->Branch("pf_pdgId", &m_pf_pdgId);
   tree_->Branch("pf_pt", &m_pf_pt);
   tree_->Branch("pf_eta", &m_pf_eta);
   tree_->Branch("pf_phi", &m_pf_phi);
   tree_->Branch("pf_charge", &m_pf_charge);

   tree_->Branch("genquark_pdgId", &m_genquark_pdgId);
   tree_->Branch("genquark_pt", &m_genquark_pt);
   tree_->Branch("genquark_eta", &m_genquark_eta);
   tree_->Branch("genquark_phi", &m_genquark_phi);
   tree_->Branch("genquark_status", &m_genquark_status);
   tree_->Branch("genquark_mother", &m_genquark_mother);
   tree_->Branch("genquark_vx", &m_genquark_vx);
   tree_->Branch("genquark_vy", &m_genquark_vy);
   tree_->Branch("genquark_vz", &m_genquark_vz);   
   tree_->Branch("genquark_charge", &m_genquark_charge);

   tree_->Branch("genmeson_pdgId", &m_genmeson_pdgId);
   tree_->Branch("genmeson_pt", &m_genmeson_pt);
   tree_->Branch("genmeson_eta", &m_genmeson_eta);
   tree_->Branch("genmeson_phi", &m_genmeson_phi);
   tree_->Branch("genmeson_status", &m_genmeson_status);
   tree_->Branch("genmeson_vx", &m_genmeson_vx);
   tree_->Branch("genmeson_vy", &m_genmeson_vy);
   tree_->Branch("genmeson_vz", &m_genmeson_vz);   
   tree_->Branch("genmeson_nfinals", &m_genmeson_nfinals);
   tree_->Branch("genmeson_ncharge", &m_genmeson_ncharge);
   tree_->Branch("genmeson_frcharge", &m_genmeson_frcharge);
   tree_->Branch("genmeson_charge", &m_genmeson_charge);
   tree_->Branch("genmeson_qjet1", &m_genmeson_qjet1);
   tree_->Branch("genmeson_qjet2", &m_genmeson_qjet2);
   tree_->Branch("genmeson_qjet3", &m_genmeson_qjet3);
   tree_->Branch("genmeson_qjet1_dr", &m_genmeson_qjet1_dr);
   tree_->Branch("genmeson_qjet2_dr", &m_genmeson_qjet2_dr);
   tree_->Branch("genmeson_qjet3_dr", &m_genmeson_qjet3_dr);
   tree_->Branch("genmeson_jet1", &m_genmeson_jet1);
   tree_->Branch("genmeson_jet2", &m_genmeson_jet2);
   tree_->Branch("genmeson_jet3", &m_genmeson_jet3);
   tree_->Branch("genmeson_jet1_dr", &m_genmeson_jet1_dr);
   tree_->Branch("genmeson_jet2_dr", &m_genmeson_jet2_dr);
   tree_->Branch("genmeson_jet3_dr", &m_genmeson_jet3_dr);
   tree_->Branch("genmeson_sv1", &m_genmeson_sv1);
   tree_->Branch("genmeson_sv2", &m_genmeson_sv2);
   tree_->Branch("genmeson_sv3", &m_genmeson_sv3);
   tree_->Branch("genmeson_sv1_dr", &m_genmeson_sv1_dr);
   tree_->Branch("genmeson_sv2_dr", &m_genmeson_sv2_dr);
   tree_->Branch("genmeson_sv3_dr", &m_genmeson_sv3_dr);


   tree_->Branch("gendaughters_pdgId", &m_gendaughters_pdgId);  
   tree_->Branch("gendaughters_mesonIdx", &m_gendaughters_mesonIdx);

   tree_->Branch("genpart_pdgId", &m_genpart_pdgId);
   tree_->Branch("genpart_pt", &m_genpart_pt);
   tree_->Branch("genpart_eta", &m_genpart_eta);
   tree_->Branch("genpart_phi", &m_genpart_phi);
   tree_->Branch("genpart_mother", &m_genpart_mother);
   tree_->Branch("genpart_Bmeson", &m_genpart_Bmeson);
   tree_->Branch("genpart_Bmeson_pt", &m_genpart_Bmeson_pt);
   tree_->Branch("genpart_Bmeson_eta", &m_genpart_Bmeson_eta);
   tree_->Branch("genpart_Bmeson_phi", &m_genpart_Bmeson_phi);
   tree_->Branch("genpart_BmesonIdx", &m_genpart_BmesonIdx);

}

// ------------ method called once each job just after ending the event loop  ------------
void
GenBNtuplizer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenBNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenBNtuplizer);
