

int GetNjets(){

  std::vector<TString> names;
  TChain * cc1 =new TChain("deepntuplizer/tree");
  cc1->Add("/eos/cms/store/cmst3/user/gkaratha/SoftMultiJet/DeepNtuples_v2/CRAB_UserFiles/PFC_Signal_cascade_m100_12_03_25/250312_143450/0000/*.root");
  names.push_back("cascade m100");
//  cc1->Add("/eos/cms/store/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v2/QCD_Bin-PT-15to7000_Par-PT-flat2022_TuneCP5_13p6TeV_pythia8/PFC_QCD_pt15to7k/250326_153724/0000/output_pfc_0_1.root");
 // names.push_back("QCD");


  TChain * cc2 =new TChain("deepntuplizer/tree");
  cc2->Add("/eos/cms/store/cmst3/user/gkaratha/SoftMultiJet/DeepNtuples_v2/CRAB_UserFiles/PFC_Signal_cascade_m220_12_03_25/250312_143736/0000/*.root");
  names.push_back("cascade m220");

  TChain * cc3 =new TChain("deepntuplizer/tree");
  cc3->Add("/eos/cms/store/cmst3/user/gkaratha/SoftMultiJet/DeepNtuples_v2/CRAB_UserFiles/PFC_Signal_chain_m70dm20_12_03_25/250312_141701/0000/*.root");
  names.push_back("chain dm20");

  TChain * cc4 =new TChain("deepntuplizer/tree");
  cc4->Add("/eos/cms/store/cmst3/user/gkaratha/SoftMultiJet/DeepNtuples_v2/CRAB_UserFiles/PFC_Signal_chain_m70dm8_12_03_25/250312_141934/0000/*.root");
  names.push_back("chain dm8");

  //TChain * cc5 =new TChain("deepntuplizer/tree");
  //cc5->Add("/eos/cms/store/cmst3/user/gkaratha/SoftMultiJet/DeepNtuples_v2/DYto2E-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/PFC_DY_2E4Jets/250318_144011/0000/*.root");
  //names.push_back("DY");


  int ismp=0;
  for (auto cc: {cc1,cc2,cc3,cc4}){
    TH1F * hmb = new TH1F("hmb_"+names[ismp],"",1,0,2);
    TH1F * hmc = new TH1F("hmc_"+names[ismp],"",1,0,2);
    TH1F * hmt = new TH1F("hmt_"+names[ismp],"",1,0,2);
    TH1F * hml = new TH1F("hml_"+names[ismp],"",1,0,2);
    
    cc->Draw("1>>hmb_"+names[ismp],"(isB || isBB || isLeptonicB || isGBB || isLeptonicB_C)");
    cc->Draw("1>>hmc_"+names[ismp],"(isC || isCC || isGCC)");
    cc->Draw("1>>hml_"+names[ismp]," (isU || isD || isS || isG)");
//    cc->Draw("1>>hmt_"+names[ismp]," (isTaup1h0p_ || isTaup1h1p_ || isTaup1h2p_ || isTaup3h0p_ || isTaup3h1p_ || isTaum1h0p_ || isTaum1h1p_ || isTaum1h2p_ || isTaum3h0p_ || isTaum3h1p_)");

    cout<<"sample "<<names[ismp]<<" B "<<hmb->Integral()<<" C "<<hmc->Integral()<<" tau "<<hmt->Integral()<<" light "<<hml->Integral()<<endl;


    TH1F * hmb_b = new TH1F("hmb_b_"+names[ismp],"",1,0,2);
    TH1F * hmb_bb = new TH1F("hmb_bb_"+names[ismp],"",1,0,2);
    TH1F * hmb_lep = new TH1F("hmb_lep_"+names[ismp],"",1,0,2);
    TH1F * hmbc_lep = new TH1F("hmbc_lep_"+names[ismp],"",1,0,2);

    cc->Draw("1>>hmb_b_"+names[ismp],"(isB )");
    cc->Draw("1>>hmb_bb_"+names[ismp],"(isBB)");
    cc->Draw("1>>hmb_lep_"+names[ismp],"(isLeptonicB)");
    cc->Draw("1>>hmbc_lep_"+names[ismp],"(isLeptonicB_C)");


    cout<<" single B "<<hmb_b->Integral()<<" double B "<<hmb_bb->Integral()<<" lep B "<<hmb_lep->Integral()<<" lep B->C "<<hmbc_lep->Integral()<<endl;


    ismp+=1;

  }



  return 0;
}
