/*
 * ntuple_V0Ks.cc
 *
 *  Created on: 3rd August 2023
 *      Author: Alexandre De Moor
 */


#include "../interface/ntuple_V0Ks.h"
// for ivf
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "TVector3.h"

const reco::Vertex * ntuple_V0Ks::spvp_;

ntuple_V0Ks::ntuple_V0Ks(std::string prefix, double jetR):ntuple_content(jetR),v0ks_num_(0){
    prefix_ = prefix;
}
ntuple_V0Ks::~ntuple_V0Ks(){}


void ntuple_V0Ks::getInput(const edm::ParameterSet& iConfig){

}

void ntuple_V0Ks::initBranches(TTree* tree){
  // SV candidates
  addBranch(tree,(prefix_+"n_v0ks").c_str()           ,&v0ks_num_         ,(prefix_+"v0ks_num_/I").c_str()     );
  addBranch(tree,(prefix_+"nv0ks").c_str()            ,&nv0ks_            ,(prefix_+"nv0ks_/F").c_str()         );
  addBranch(tree,(prefix_+"v0ks_pt").c_str()          ,&v0ks_pt_          ,(prefix_+"v0ks_pt_["+prefix_+"v0ks_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0ks_px").c_str()          ,&v0ks_px_          ,(prefix_+"v0ks_px_["+prefix_+"v0ks_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0ks_py").c_str()          ,&v0ks_py_          ,(prefix_+"v0ks_py_["+prefix_+"v0ks_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0ks_pz").c_str()          ,&v0ks_pz_          ,(prefix_+"v0ks_pz_["+prefix_+"v0ks_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0ks_eta").c_str()         ,&v0ks_eta_         ,(prefix_+"v0ks_eta_["+prefix_+"v0ks_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0ks_phi").c_str()         ,&v0ks_phi_         ,(prefix_+"v0ks_phi_["+prefix_+"v0ks_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0ks_e").c_str()           ,&v0ks_e_           ,(prefix_+"v0ks_e_["+prefix_+"v0ks_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0ks_etarel").c_str()      ,&v0ks_etarel_      ,(prefix_+"v0ks_etarel_["+prefix_+"v0ks_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0ks_phirel").c_str()      ,&v0ks_phirel_      ,(prefix_+"v0ks_phirel_["+prefix_+"v0ks_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0ks_deltaR").c_str()      ,&v0ks_deltaR_      ,(prefix_+"v0ks_deltaR_["+prefix_+"v0ks_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0ks_mass").c_str()        ,&v0ks_mass_        ,(prefix_+"v0ks_mass_["+prefix_+"v0ks_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0ks_ntracks").c_str()     ,&v0ks_ntracks_     ,(prefix_+"v0ks_ntracks_["+prefix_+"v0ks_num_]/F").c_str()     );
  addBranch(tree,(prefix_+"v0ks_chi2").c_str()        ,&v0ks_chi2_        ,(prefix_+"v0ks_chi2_["+prefix_+"v0ks_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0ks_ndf").c_str()         ,&v0ks_ndf_         ,(prefix_+"v0ks_ndf_["+prefix_+"v0ks_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0ks_normchi2").c_str()    ,&v0ks_normchi2_    ,(prefix_+"v0ks_normchi2_["+prefix_+"v0ks_num_]/F").c_str()     );
  addBranch(tree,(prefix_+"v0ks_dxy").c_str()         ,&v0ks_dxy_         ,(prefix_+"v0ks_dxy_["+prefix_+"v0ks_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0ks_dxyerr").c_str()      ,&v0ks_dxyerr_      ,(prefix_+"v0ks_dxyerr_["+prefix_+"v0ks_num_]/F").c_str()      );
  addBranch(tree,(prefix_+"v0ks_dxysig").c_str()      ,&v0ks_dxysig_      ,(prefix_+"v0ks_dxysig_["+prefix_+"v0ks_num_]/F").c_str()      );
  addBranch(tree,(prefix_+"v0ks_d3d").c_str()         ,&v0ks_d3d_         ,(prefix_+"v0ks_d3d_["+prefix_+"v0ks_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0ks_d3derr").c_str()      ,&v0ks_d3derr_      ,(prefix_+"v0ks_d3err_["+prefix_+"v0ks_num_]/F").c_str()       );
  addBranch(tree,(prefix_+"v0ks_d3dsig").c_str()      ,&v0ks_d3dsig_      ,(prefix_+"v0ks_d3dsig_["+prefix_+"v0ks_num_]/F").c_str()      );
  addBranch(tree,(prefix_+"v0ks_costhetasvpv").c_str(),&v0ks_costhetasvpv_,(prefix_+"v0ks_costhetasvpv_["+prefix_+"v0ks_num_]/F").c_str());
  addBranch(tree,(prefix_+"v0ks_enratio").c_str()     ,&v0ks_enratio_     ,(prefix_+"v0ks_enratio_["+prefix_+"v0ks_num_]/F").c_str());

}

void ntuple_V0Ks::readSetup(const edm::EventSetup& iSetup){

}

void ntuple_V0Ks::readEvent(const edm::Event& iEvent){


}


bool ntuple_V0Ks::compareDxyDxyErr(const reco::VertexCompositePtrCandidate &sva,const reco::VertexCompositePtrCandidate &svb){
    reco::Vertex pv = *spvp_;
    float adxy = ntuple_V0Ks::vertexDxy(sva,pv).value();
    float bdxy = ntuple_V0Ks::vertexDxy(svb,pv).value();
    float aerr = ntuple_V0Ks::vertexDxy(sva,pv).error();
    float berr = ntuple_V0Ks::vertexDxy(svb,pv).error();

    float asig = ntuple_V0Ks::catchInfs(adxy/aerr,0.);
    float bsig = ntuple_V0Ks::catchInfs(bdxy/berr,0.);
    return bsig<asig;
}

bool ntuple_V0Ks::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){

    const float jet_uncorr_e = -1;// jet.correctedJet("Uncorrected").energy();
    const reco::Vertex & pv = vertices()->at(0);
    GlobalVector jetRefTrackDir(jet.px(),jet.py(),jet.pz());

    v0ks_num_ = 0;

    reco::VertexCompositePtrCandidateCollection cpvtx=*V0ks();
    spvp_ =   & vertices()->at(0);
    std::sort(cpvtx.begin(),cpvtx.end(),ntuple_V0Ks::compareDxyDxyErr);

    float etasign=1;
    etasign++; //avoid unused warning
    if(jet.eta()<0)etasign=-1;

    double jet_radius = jetR();
    if (jet_radius<0){
      // subjets: use maxDR(subjet, pfcand)
      for (unsigned idau=0; idau<jet.numberOfDaughters(); ++idau){
        double dR = reco::deltaR(*jet.daughter(idau), jet);
        if (dR>jet_radius)
          jet_radius = dR;
      }
    }

    for (const reco::VertexCompositePtrCandidate &sv : cpvtx) {

        if (reco::deltaR(sv,jet)>jet_radius) { continue; }
        if((int)max_sv>v0ks_num_){

            v0ks_pt_[v0ks_num_]           = sv.pt();
            v0ks_px_[v0ks_num_]           = sv.px();
            v0ks_py_[v0ks_num_]           = sv.py();
            v0ks_pz_[v0ks_num_]           = sv.pz();
            v0ks_eta_[v0ks_num_]          = sv.eta();
            v0ks_phi_[v0ks_num_]          = sv.phi();
            v0ks_etarel_[v0ks_num_]       = catchInfsAndBound(fabs(sv.eta()-jet.eta())-0.5,0,-2,0);
            v0ks_phirel_[v0ks_num_]       = catchInfsAndBound(fabs(reco::deltaPhi(sv.phi(),jet.phi()))-0.5,0,-2,0);
            v0ks_deltaR_[v0ks_num_]       = catchInfsAndBound(fabs(reco::deltaR(sv,jet))-0.5,0,-2,0);
            v0ks_mass_[v0ks_num_]         = sv.mass();
            v0ks_ntracks_[v0ks_num_]      = sv.numberOfDaughters();
            v0ks_chi2_[v0ks_num_]         = sv.vertexChi2();
            v0ks_ndf_[v0ks_num_]          = sv.vertexNdof();
            v0ks_normchi2_[v0ks_num_]     = catchInfsAndBound(v0ks_chi2_[v0ks_num_]/v0ks_ndf_[v0ks_num_],1000,-1000,1000);
            v0ks_dxy_[v0ks_num_]          = vertexDxy(sv,pv).value();
            v0ks_dxyerr_[v0ks_num_]       = catchInfsAndBound(vertexDxy(sv,pv).error()-2,0,-2,0);
            v0ks_dxysig_[v0ks_num_]       = catchInfsAndBound(v0ks_dxy_[v0ks_num_]/vertexDxy(sv,pv).error() ,0,-1,800);
            v0ks_d3d_[v0ks_num_]          = vertexD3d(sv,pv).value();
            v0ks_d3derr_[v0ks_num_]       = catchInfsAndBound(vertexD3d(sv,pv).error()-2,0,-2,0);
            v0ks_d3dsig_[v0ks_num_]       = catchInfsAndBound(vertexD3d(sv,pv).value()/vertexD3d(sv,pv).error(),0,-1,800);
            v0ks_costhetasvpv_[v0ks_num_] = vertexDdotP(sv,pv); // the pointing angle (i.e. the angle between the sum of the momentum
            // of the tracks in the SV and the flight direction betwen PV and SV)

            v0ks_enratio_[v0ks_num_]      = sv.energy()/jet_uncorr_e;
            v0ks_e_[v0ks_num_]            = sv.energy();

            v0ks_num_++;
        }
    } // end of looping over the secondary vertices
    nv0ks_=v0ks_num_;

    return true;
}

///helpers seldomly touched
Measurement1D ntuple_V0Ks::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

Measurement1D ntuple_V0Ks::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

float ntuple_V0Ks::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv)  {
    reco::Candidate::Vector p = sv.momentum();
    reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
    return p.Unit().Dot(d.Unit());
}
