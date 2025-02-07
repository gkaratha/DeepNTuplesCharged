//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Feb  1 17:30:01 2025 by ROOT version 6.30/03
// from TTree tree/tree
// found on file: output1p5cPFCvsPuppiRecluster_0.root
//////////////////////////////////////////////////////////

#ifndef jet_tree_h
#define jet_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class jet_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nqjet;
   vector<float>   *qjet_pt;
   vector<float>   *qjet_eta;
   vector<float>   *qjet_phi;
   vector<float>   *qjet_mass;
   vector<float>   *qjet_energy;
   vector<float>   *qjet_et;
   vector<float>   *qjet_vx;
   vector<float>   *qjet_vy;
   vector<float>   *qjet_vz;
   vector<float>   *qjet_ntrk;
   Int_t           njet;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_mass;
   vector<float>   *jet_energy;
   vector<float>   *jet_et;
   vector<float>   *jet_vx;
   vector<float>   *jet_vy;
   vector<float>   *jet_vz;
   vector<float>   *jet_ntrk;
   vector<float>   *jet_npf;
   vector<int>     *jetpart_jetIdx;
   vector<float>   *jetpart_pt;
   vector<float>   *jetpart_eta;
   vector<float>   *jetpart_phi;
   vector<float>   *jetpart_charge;
   vector<int>     *jetpart_pdgId;
   Int_t           nsv;
   vector<float>   *sv_pt;
   vector<float>   *sv_eta;
   vector<float>   *sv_phi;
   vector<float>   *sv_mass;
   vector<float>   *sv_energy;
   vector<float>   *sv_et;
   vector<float>   *sv_vx;
   vector<float>   *sv_vy;
   vector<float>   *sv_vz;
   vector<float>   *sv_ntrk;
   vector<int>     *pf_pdgId;
   vector<float>   *pf_pt;
   vector<float>   *pf_eta;
   vector<float>   *pf_phi;
   vector<float>   *pf_charge;
   vector<float>   *genquark_pdgId;
   vector<float>   *genquark_pt;
   vector<float>   *genquark_eta;
   vector<float>   *genquark_phi;
   vector<float>   *genquark_status;
   vector<float>   *genquark_mother;
   vector<float>   *genquark_vx;
   vector<float>   *genquark_vy;
   vector<float>   *genquark_vz;
   vector<float>   *genquark_charge;
   vector<float>   *genmeson_pdgId;
   vector<float>   *genmeson_pt;
   vector<float>   *genmeson_eta;
   vector<float>   *genmeson_phi;
   vector<float>   *genmeson_status;
   vector<float>   *genmeson_vx;
   vector<float>   *genmeson_vy;
   vector<float>   *genmeson_vz;
   vector<float>   *genmeson_nfinals;
   vector<float>   *genmeson_ncharge;
   vector<float>   *genmeson_frcharge;
   vector<float>   *genmeson_charge;
   vector<float>   *genmeson_qjet1;
   vector<float>   *genmeson_qjet2;
   vector<float>   *genmeson_qjet3;
   vector<float>   *genmeson_qjet1_dr;
   vector<float>   *genmeson_qjet2_dr;
   vector<float>   *genmeson_qjet3_dr;
   vector<float>   *genmeson_jet1;
   vector<float>   *genmeson_jet2;
   vector<float>   *genmeson_jet3;
   vector<float>   *genmeson_jet1_dr;
   vector<float>   *genmeson_jet2_dr;
   vector<float>   *genmeson_jet3_dr;
   vector<float>   *genmeson_sv1;
   vector<float>   *genmeson_sv2;
   vector<float>   *genmeson_sv3;
   vector<float>   *genmeson_sv1_dr;
   vector<float>   *genmeson_sv2_dr;
   vector<float>   *genmeson_sv3_dr;
   vector<float>   *gendaughters_pdgId;
   vector<float>   *gendaughters_mesonIdx;
   vector<float>   *genpart_pdgId;
   vector<float>   *genpart_pt;
   vector<float>   *genpart_eta;
   vector<float>   *genpart_phi;
   vector<float>   *genpart_mother;
   vector<float>   *genpart_Bmeson;
   vector<float>   *genpart_Bmeson_pt;
   vector<float>   *genpart_Bmeson_eta;
   vector<float>   *genpart_Bmeson_phi;
   vector<float>   *genpart_BmesonIdx;

   // List of branches
   TBranch        *b_nqjet;   //!
   TBranch        *b_qjet_pt;   //!
   TBranch        *b_qjet_eta;   //!
   TBranch        *b_qjet_phi;   //!
   TBranch        *b_qjet_mass;   //!
   TBranch        *b_qjet_energy;   //!
   TBranch        *b_qjet_et;   //!
   TBranch        *b_qjet_vx;   //!
   TBranch        *b_qjet_vy;   //!
   TBranch        *b_qjet_vz;   //!
   TBranch        *b_qjet_ntrk;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_mass;   //!
   TBranch        *b_jet_energy;   //!
   TBranch        *b_jet_et;   //!
   TBranch        *b_jet_vx;   //!
   TBranch        *b_jet_vy;   //!
   TBranch        *b_jet_vz;   //!
   TBranch        *b_jet_ntrk;   //!
   TBranch        *b_jet_npf;   //!
   TBranch        *b_jetpart_jetIdx;   //!
   TBranch        *b_jetpart_pt;   //!
   TBranch        *b_jetpart_eta;   //!
   TBranch        *b_jetpart_phi;   //!
   TBranch        *b_jetpart_charge;   //!
   TBranch        *b_jetpart_pdgId;   //!
   TBranch        *b_nsv;   //!
   TBranch        *b_sv_pt;   //!
   TBranch        *b_sv_eta;   //!
   TBranch        *b_sv_phi;   //!
   TBranch        *b_sv_mass;   //!
   TBranch        *b_sv_energy;   //!
   TBranch        *b_sv_et;   //!
   TBranch        *b_sv_vx;   //!
   TBranch        *b_sv_vy;   //!
   TBranch        *b_sv_vz;   //!
   TBranch        *b_sv_ntrk;   //!
   TBranch        *b_pf_pdgId;   //!
   TBranch        *b_pf_pt;   //!
   TBranch        *b_pf_eta;   //!
   TBranch        *b_pf_phi;   //!
   TBranch        *b_pf_charge;   //!
   TBranch        *b_genquark_pdgId;   //!
   TBranch        *b_genquark_pt;   //!
   TBranch        *b_genquark_eta;   //!
   TBranch        *b_genquark_phi;   //!
   TBranch        *b_genquark_status;   //!
   TBranch        *b_genquark_mother;   //!
   TBranch        *b_genquark_vx;   //!
   TBranch        *b_genquark_vy;   //!
   TBranch        *b_genquark_vz;   //!
   TBranch        *b_genquark_charge;   //!
   TBranch        *b_genmeson_pdgId;   //!
   TBranch        *b_genmeson_pt;   //!
   TBranch        *b_genmeson_eta;   //!
   TBranch        *b_genmeson_phi;   //!
   TBranch        *b_genmeson_status;   //!
   TBranch        *b_genmeson_vx;   //!
   TBranch        *b_genmeson_vy;   //!
   TBranch        *b_genmeson_vz;   //!
   TBranch        *b_genmeson_nfinals;   //!
   TBranch        *b_genmeson_ncharge;   //!
   TBranch        *b_genmeson_frcharge;   //!
   TBranch        *b_genmeson_charge;   //!
   TBranch        *b_genmeson_qjet1;   //!
   TBranch        *b_genmeson_qjet2;   //!
   TBranch        *b_genmeson_qjet3;   //!
   TBranch        *b_genmeson_qjet1_dr;   //!
   TBranch        *b_genmeson_qjet2_dr;   //!
   TBranch        *b_genmeson_qjet3_dr;   //!
   TBranch        *b_genmeson_jet1;   //!
   TBranch        *b_genmeson_jet2;   //!
   TBranch        *b_genmeson_jet3;   //!
   TBranch        *b_genmeson_jet1_dr;   //!
   TBranch        *b_genmeson_jet2_dr;   //!
   TBranch        *b_genmeson_jet3_dr;   //!
   TBranch        *b_genmeson_sv1;   //!
   TBranch        *b_genmeson_sv2;   //!
   TBranch        *b_genmeson_sv3;   //!
   TBranch        *b_genmeson_sv1_dr;   //!
   TBranch        *b_genmeson_sv2_dr;   //!
   TBranch        *b_genmeson_sv3_dr;   //!
   TBranch        *b_gendaughters_pdgId;   //!
   TBranch        *b_gendaughters_mesonIdx;   //!
   TBranch        *b_genpart_pdgId;   //!
   TBranch        *b_genpart_pt;   //!
   TBranch        *b_genpart_eta;   //!
   TBranch        *b_genpart_phi;   //!
   TBranch        *b_genpart_mother;   //!
   TBranch        *b_genpart_Bmeson;   //!
   TBranch        *b_genpart_Bmeson_pt;   //!
   TBranch        *b_genpart_Bmeson_eta;   //!
   TBranch        *b_genpart_Bmeson_phi;   //!
   TBranch        *b_genpart_BmesonIdx;   //!

   jet_tree(TTree *tree=0);
   virtual ~jet_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

jet_tree::jet_tree(TTree *tree) : fChain(0) 
{
   Init(tree);
}

jet_tree::~jet_tree()
{
   if (!fChain) return;
}

Int_t jet_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t jet_tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void jet_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   qjet_pt = 0;
   qjet_eta = 0;
   qjet_phi = 0;
   qjet_mass = 0;
   qjet_energy = 0;
   qjet_et = 0;
   qjet_vx = 0;
   qjet_vy = 0;
   qjet_vz = 0;
   qjet_ntrk = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_mass = 0;
   jet_energy = 0;
   jet_et = 0;
   jet_vx = 0;
   jet_vy = 0;
   jet_vz = 0;
   jet_ntrk = 0;
   jet_npf = 0;
   jetpart_jetIdx = 0;
   jetpart_pt = 0;
   jetpart_eta = 0;
   jetpart_phi = 0;
   jetpart_charge = 0;
   jetpart_pdgId = 0;
   sv_pt = 0;
   sv_eta = 0;
   sv_phi = 0;
   sv_mass = 0;
   sv_energy = 0;
   sv_et = 0;
   sv_vx = 0;
   sv_vy = 0;
   sv_vz = 0;
   sv_ntrk = 0;
   pf_pdgId = 0;
   pf_pt = 0;
   pf_eta = 0;
   pf_phi = 0;
   pf_charge = 0;
   genquark_pdgId = 0;
   genquark_pt = 0;
   genquark_eta = 0;
   genquark_phi = 0;
   genquark_status = 0;
   genquark_mother = 0;
   genquark_vx = 0;
   genquark_vy = 0;
   genquark_vz = 0;
   genquark_charge = 0;
   genmeson_pdgId = 0;
   genmeson_pt = 0;
   genmeson_eta = 0;
   genmeson_phi = 0;
   genmeson_status = 0;
   genmeson_vx = 0;
   genmeson_vy = 0;
   genmeson_vz = 0;
   genmeson_nfinals = 0;
   genmeson_ncharge = 0;
   genmeson_frcharge = 0;
   genmeson_charge = 0;
   genmeson_qjet1 = 0;
   genmeson_qjet2 = 0;
   genmeson_qjet3 = 0;
   genmeson_qjet1_dr = 0;
   genmeson_qjet2_dr = 0;
   genmeson_qjet3_dr = 0;
   genmeson_jet1 = 0;
   genmeson_jet2 = 0;
   genmeson_jet3 = 0;
   genmeson_jet1_dr = 0;
   genmeson_jet2_dr = 0;
   genmeson_jet3_dr = 0;
   genmeson_sv1 = 0;
   genmeson_sv2 = 0;
   genmeson_sv3 = 0;
   genmeson_sv1_dr = 0;
   genmeson_sv2_dr = 0;
   genmeson_sv3_dr = 0;
   gendaughters_pdgId = 0;
   gendaughters_mesonIdx = 0;
   genpart_pdgId = 0;
   genpart_pt = 0;
   genpart_eta = 0;
   genpart_phi = 0;
   genpart_mother = 0;
   genpart_Bmeson = 0;
   genpart_Bmeson_pt = 0;
   genpart_Bmeson_eta = 0;
   genpart_Bmeson_phi = 0;
   genpart_BmesonIdx = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nqjet", &nqjet, &b_nqjet);
   fChain->SetBranchAddress("qjet_pt", &qjet_pt, &b_qjet_pt);
   fChain->SetBranchAddress("qjet_eta", &qjet_eta, &b_qjet_eta);
   fChain->SetBranchAddress("qjet_phi", &qjet_phi, &b_qjet_phi);
   fChain->SetBranchAddress("qjet_mass", &qjet_mass, &b_qjet_mass);
   fChain->SetBranchAddress("qjet_energy", &qjet_energy, &b_qjet_energy);
   fChain->SetBranchAddress("qjet_et", &qjet_et, &b_qjet_et);
   fChain->SetBranchAddress("qjet_vx", &qjet_vx, &b_qjet_vx);
   fChain->SetBranchAddress("qjet_vy", &qjet_vy, &b_qjet_vy);
   fChain->SetBranchAddress("qjet_vz", &qjet_vz, &b_qjet_vz);
   fChain->SetBranchAddress("qjet_ntrk", &qjet_ntrk, &b_qjet_ntrk);
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_mass", &jet_mass, &b_jet_mass);
   fChain->SetBranchAddress("jet_energy", &jet_energy, &b_jet_energy);
   fChain->SetBranchAddress("jet_et", &jet_et, &b_jet_et);
   fChain->SetBranchAddress("jet_vx", &jet_vx, &b_jet_vx);
   fChain->SetBranchAddress("jet_vy", &jet_vy, &b_jet_vy);
   fChain->SetBranchAddress("jet_vz", &jet_vz, &b_jet_vz);
   fChain->SetBranchAddress("jet_ntrk", &jet_ntrk, &b_jet_ntrk);
   fChain->SetBranchAddress("jet_npf", &jet_npf, &b_jet_npf);
   fChain->SetBranchAddress("jetpart_jetIdx", &jetpart_jetIdx, &b_jetpart_jetIdx);
   fChain->SetBranchAddress("jetpart_pt", &jetpart_pt, &b_jetpart_pt);
   fChain->SetBranchAddress("jetpart_eta", &jetpart_eta, &b_jetpart_eta);
   fChain->SetBranchAddress("jetpart_phi", &jetpart_phi, &b_jetpart_phi);
   fChain->SetBranchAddress("jetpart_charge", &jetpart_charge, &b_jetpart_charge);
   fChain->SetBranchAddress("jetpart_pdgId", &jetpart_pdgId, &b_jetpart_pdgId);
   fChain->SetBranchAddress("nsv", &nsv, &b_nsv);
   fChain->SetBranchAddress("sv_pt", &sv_pt, &b_sv_pt);
   fChain->SetBranchAddress("sv_eta", &sv_eta, &b_sv_eta);
   fChain->SetBranchAddress("sv_phi", &sv_phi, &b_sv_phi);
   fChain->SetBranchAddress("sv_mass", &sv_mass, &b_sv_mass);
   fChain->SetBranchAddress("sv_energy", &sv_energy, &b_sv_energy);
   fChain->SetBranchAddress("sv_et", &sv_et, &b_sv_et);
   fChain->SetBranchAddress("sv_vx", &sv_vx, &b_sv_vx);
   fChain->SetBranchAddress("sv_vy", &sv_vy, &b_sv_vy);
   fChain->SetBranchAddress("sv_vz", &sv_vz, &b_sv_vz);
   fChain->SetBranchAddress("sv_ntrk", &sv_ntrk, &b_sv_ntrk);
   fChain->SetBranchAddress("pf_pdgId", &pf_pdgId, &b_pf_pdgId);
   fChain->SetBranchAddress("pf_pt", &pf_pt, &b_pf_pt);
   fChain->SetBranchAddress("pf_eta", &pf_eta, &b_pf_eta);
   fChain->SetBranchAddress("pf_phi", &pf_phi, &b_pf_phi);
   fChain->SetBranchAddress("pf_charge", &pf_charge, &b_pf_charge);
   fChain->SetBranchAddress("genquark_pdgId", &genquark_pdgId, &b_genquark_pdgId);
   fChain->SetBranchAddress("genquark_pt", &genquark_pt, &b_genquark_pt);
   fChain->SetBranchAddress("genquark_eta", &genquark_eta, &b_genquark_eta);
   fChain->SetBranchAddress("genquark_phi", &genquark_phi, &b_genquark_phi);
   fChain->SetBranchAddress("genquark_status", &genquark_status, &b_genquark_status);
   fChain->SetBranchAddress("genquark_mother", &genquark_mother, &b_genquark_mother);
   fChain->SetBranchAddress("genquark_vx", &genquark_vx, &b_genquark_vx);
   fChain->SetBranchAddress("genquark_vy", &genquark_vy, &b_genquark_vy);
   fChain->SetBranchAddress("genquark_vz", &genquark_vz, &b_genquark_vz);
   fChain->SetBranchAddress("genquark_charge", &genquark_charge, &b_genquark_charge);
   fChain->SetBranchAddress("genmeson_pdgId", &genmeson_pdgId, &b_genmeson_pdgId);
   fChain->SetBranchAddress("genmeson_pt", &genmeson_pt, &b_genmeson_pt);
   fChain->SetBranchAddress("genmeson_eta", &genmeson_eta, &b_genmeson_eta);
   fChain->SetBranchAddress("genmeson_phi", &genmeson_phi, &b_genmeson_phi);
   fChain->SetBranchAddress("genmeson_status", &genmeson_status, &b_genmeson_status);
   fChain->SetBranchAddress("genmeson_vx", &genmeson_vx, &b_genmeson_vx);
   fChain->SetBranchAddress("genmeson_vy", &genmeson_vy, &b_genmeson_vy);
   fChain->SetBranchAddress("genmeson_vz", &genmeson_vz, &b_genmeson_vz);
   fChain->SetBranchAddress("genmeson_nfinals", &genmeson_nfinals, &b_genmeson_nfinals);
   fChain->SetBranchAddress("genmeson_ncharge", &genmeson_ncharge, &b_genmeson_ncharge);
   fChain->SetBranchAddress("genmeson_frcharge", &genmeson_frcharge, &b_genmeson_frcharge);
   fChain->SetBranchAddress("genmeson_charge", &genmeson_charge, &b_genmeson_charge);
   fChain->SetBranchAddress("genmeson_qjet1", &genmeson_qjet1, &b_genmeson_qjet1);
   fChain->SetBranchAddress("genmeson_qjet2", &genmeson_qjet2, &b_genmeson_qjet2);
   fChain->SetBranchAddress("genmeson_qjet3", &genmeson_qjet3, &b_genmeson_qjet3);
   fChain->SetBranchAddress("genmeson_qjet1_dr", &genmeson_qjet1_dr, &b_genmeson_qjet1_dr);
   fChain->SetBranchAddress("genmeson_qjet2_dr", &genmeson_qjet2_dr, &b_genmeson_qjet2_dr);
   fChain->SetBranchAddress("genmeson_qjet3_dr", &genmeson_qjet3_dr, &b_genmeson_qjet3_dr);
   fChain->SetBranchAddress("genmeson_jet1", &genmeson_jet1, &b_genmeson_jet1);
   fChain->SetBranchAddress("genmeson_jet2", &genmeson_jet2, &b_genmeson_jet2);
   fChain->SetBranchAddress("genmeson_jet3", &genmeson_jet3, &b_genmeson_jet3);
   fChain->SetBranchAddress("genmeson_jet1_dr", &genmeson_jet1_dr, &b_genmeson_jet1_dr);
   fChain->SetBranchAddress("genmeson_jet2_dr", &genmeson_jet2_dr, &b_genmeson_jet2_dr);
   fChain->SetBranchAddress("genmeson_jet3_dr", &genmeson_jet3_dr, &b_genmeson_jet3_dr);
   fChain->SetBranchAddress("genmeson_sv1", &genmeson_sv1, &b_genmeson_sv1);
   fChain->SetBranchAddress("genmeson_sv2", &genmeson_sv2, &b_genmeson_sv2);
   fChain->SetBranchAddress("genmeson_sv3", &genmeson_sv3, &b_genmeson_sv3);
   fChain->SetBranchAddress("genmeson_sv1_dr", &genmeson_sv1_dr, &b_genmeson_sv1_dr);
   fChain->SetBranchAddress("genmeson_sv2_dr", &genmeson_sv2_dr, &b_genmeson_sv2_dr);
   fChain->SetBranchAddress("genmeson_sv3_dr", &genmeson_sv3_dr, &b_genmeson_sv3_dr);
   fChain->SetBranchAddress("gendaughters_pdgId", &gendaughters_pdgId, &b_gendaughters_pdgId);
   fChain->SetBranchAddress("gendaughters_mesonIdx", &gendaughters_mesonIdx, &b_gendaughters_mesonIdx);
   fChain->SetBranchAddress("genpart_pdgId", &genpart_pdgId, &b_genpart_pdgId);
   fChain->SetBranchAddress("genpart_pt", &genpart_pt, &b_genpart_pt);
   fChain->SetBranchAddress("genpart_eta", &genpart_eta, &b_genpart_eta);
   fChain->SetBranchAddress("genpart_phi", &genpart_phi, &b_genpart_phi);
   fChain->SetBranchAddress("genpart_mother", &genpart_mother, &b_genpart_mother);
   fChain->SetBranchAddress("genpart_Bmeson", &genpart_Bmeson, &b_genpart_Bmeson);
   fChain->SetBranchAddress("genpart_Bmeson_pt", &genpart_Bmeson_pt, &b_genpart_Bmeson_pt);
   fChain->SetBranchAddress("genpart_Bmeson_eta", &genpart_Bmeson_eta, &b_genpart_Bmeson_eta);
   fChain->SetBranchAddress("genpart_Bmeson_phi", &genpart_Bmeson_phi, &b_genpart_Bmeson_phi);
   fChain->SetBranchAddress("genpart_BmesonIdx", &genpart_BmesonIdx, &b_genpart_BmesonIdx);
   Notify();
}

Bool_t jet_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void jet_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t jet_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef jet_tree_cxx
