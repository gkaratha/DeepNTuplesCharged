

int EtaPtCoverage(){
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  std::vector<TString> names;
  TChain * cc =new TChain("deepntuplizer/tree");
  cc->Add("/eos/cms/store/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v3/CRAB_UserFiles/PFC_Signal_cascade_m100_31_13_04_25/250413_155917/0000/*1.root");
  cc->Add("/eos/cms/store/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v3/CRAB_UserFiles/PFC_Signal_cascade_m220_67_13_04_25/250415_133413/0000/*1.root");
  cc->Add("/eos/cms/store/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v3/CRAB_UserFiles/PFC_Signal_chain_m70_dm20_13_04_25/250413_155615/0000/*1.root");
  cc->Add("/eos/cms/store/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v3/CRAB_UserFiles/PFC_Signal_chain_m70_dm8_13_04_25/250413_155726/0000/*1.root");
//  cc->Add("/eos/cms/store/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v2/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/PFC_QCD_pt15to7k_ext1/250326_153555/0000/*1.root");


  std::vector<TString> pt_branches={"jet_pt"};
  std::vector<TString> eta_branches={"jet_eta"};
  std::vector<TString> name={"jet_newlabel"};
  TString Blabel= "isMatchedB2d || isMatchedB3d || isMatchedBMore3d  ";
  TString Clabel= "isMatchedC2d || isMatchedC3d || isMatchedCMore3d";
  TString BKGlabel= "isNotMatched";


for (int idx=0; idx<name.size(); idx++){
  TH2F * hpt_etaB = new TH2F("hpt_etaB"+name[idx],"",10,-2.5,2.5,8,0,40);
  cc->Draw(pt_branches[idx]+":"+eta_branches[idx]+">>hpt_etaB"+name[idx],Blabel);
  cout<<"hpt_etaB "<<hpt_etaB->Integral()<<endl;
  TCanvas * c1 = new TCanvas("c1","",800,600);
  c1->SetRightMargin(0.2);
  hpt_etaB->Draw("COLZ");
  hpt_etaB->GetXaxis()->SetTitle("eta"); 
  hpt_etaB->GetYaxis()->SetTitle("pt");   
  c1->SaveAs("hpteta_B"+name[idx]+".png");


  TH2F * hpt_etaUDS = new TH2F("hpt_etaUDS"+name[idx],"",10,-2.5,2.5,8,0,40);
  cc->Draw(pt_branches[idx]+":"+eta_branches[idx]+">>hpt_etaUDS"+name[idx],BKGlabel);
  cout<<"hpt_etaUDSG "<<hpt_etaUDS->Integral()<<endl;
  TCanvas * c2 = new TCanvas("c2","",800,600);
  c2->SetRightMargin(0.2);
  hpt_etaUDS->Draw("COLZ");
  hpt_etaUDS->GetXaxis()->SetTitle("eta");   
  hpt_etaUDS->GetYaxis()->SetTitle("pt"); 
  c2->SaveAs("hpteta_UDS"+name[idx]+".png");


  TH2F * hpt_etaC = new TH2F("hpt_etaC"+name[idx],"",10,-2.5,2.5,8,0,40);
  cc->Draw(pt_branches[idx]+":"+eta_branches[idx]+">>hpt_etaC"+name[idx],Clabel);
  cout<<"hpt_etaC "<<hpt_etaC->Integral()<<endl;
  TCanvas * c4 = new TCanvas("c4","",800,600);
  c4->SetRightMargin(0.2);
  hpt_etaC->Draw("COLZ");
  hpt_etaC->GetXaxis()->SetTitle("eta");   
  hpt_etaC->GetYaxis()->SetTitle("pt"); 
  c4->SaveAs("hpteta_C"+name[idx]+".png");

}


  return 0;
}
