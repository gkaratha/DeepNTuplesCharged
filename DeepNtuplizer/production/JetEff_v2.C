#include "JetEff_helper.h"
#include "jet_tree.h"




int JetEff_v2(){

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  
  TChain * cc_puppi = new TChain("genbanalizer/tree");
  cc_puppi->Add("outputSignalBkgPupCvsPuppiRecluster_0.root");
  
  TChain * cc_pfq0 = new TChain("genbanalizer/tree");
  cc_pfq0->Add("outputSignalBkgPVQual4_pt0p75_dz1PFCvsPFRecluster_0.root");
  
  TChain * cc_pfq4 = new TChain("genbanalizer/tree");
  cc_pfq4->Add("outputSignalBkgPVQual4PFCvsSlimmed_0.root");
  
  jet_tree tree_puppi;
  tree_puppi.Init(cc_puppi);
  
  jet_tree tree_pfq0;
  tree_pfq0.Init(cc_pfq0);
  
  jet_tree tree_pfq4;
  tree_pfq4.Init(cc_pfq4);
  
  
  TH1F* hdr_bb = new TH1F("hdr_bb","",100,0,3); 
  std::vector<TH1F*> hdr_jets = {new TH1F("hdr_jet_puppi","",100,0,1), new TH1F("hdr_jet_pf","",100,0,1) };
  std::vector<TH1F*> hdr_qjets = {new TH1F("hdr_jet_pupc","",100,0,1), new TH1F("hdr_jet_pfc","",100,0,1), new TH1F("hdr_jet_pfcq4","",100,0,1) };
  
  TH1F* hdr_sv = new TH1F("hdr_sv","",100,0,1);
  TH1F* heffpt_den = new TH1F("heffpt_den","",100,0,30);
  TH1F* heffeta_den = new TH1F("heffeta_den","",100,-2.3,2.3);
  
  vector<TH1F*> heff_jetpt_nums = { new TH1F("heffpt_puppi","",100,0,30), new TH1F("heffpt_pf","",100,0,30)};
  vector<TH1F*> heff_jeteta_nums = { new TH1F("heffeta_puppi","",100,-2.3,2.3), new TH1F("heffeta_pf","",100,-2.3,2.3)};
  
  vector<TH1F*> heff_qjetpt_nums = { new TH1F("heffpt_pupc","",100,0,30), new TH1F("heffpt_pfc","",100,0,30),new TH1F("heffpt_pfcq4","",100,0,30)};
  vector<TH1F*> heff_qjeteta_nums = { new TH1F("heffeta_pupc","",100,-2.3,2.3), new TH1F("heffeta_pfc","",100,-2.3,2.3), new TH1F("heffeta_pfcq4","",100,-2.3,2.3)};
  
  
  vector<TH1F*> hrate_jetpt = { new TH1F("hratept_puppi","",100,0,30), new TH1F("hratept_pf","",100,0,30)};
  
  vector<TH1F*> hrate_qjetpt = { new TH1F("hratept_pupc","",100,0,30), new TH1F("hratept_pfc","",100,0,30),new TH1F("hratept_pfcq4","",100,0,30)};
  
  
  
  bool CheckMatch=false;
  bool CheckGenContentInTrees=false;
  int evt_to_check=100;
  bool CheckGenFilter=true;
  
  
  
  if (CheckGenContentInTrees ){
     int ijk=1;
     for (jet_tree tree: {tree_pfq0,tree_pfq4}){
         if (tree.GetEntries()!= tree_puppi.GetEntries())
            std::cout<<"Difference in evts between 0 and "<<ijk<<" entries "<<tree_puppi.GetEntries()<<"  "<<tree.GetEntries()<<std::endl;
         ijk+=1;
     }
  
     std::vector<float> genBpt1,genBeta1,genBphi1;
     for (int ievt=0; ievt<evt_to_check; ievt++){
        tree_puppi.GetEntry(ievt);
        for (int ib=0; ib<tree_puppi.genquark_pt->size(); ib++){
            genBpt1.push_back(tree_puppi.genmeson_pt->at(ib));
            genBeta1.push_back(tree_puppi.genmeson_eta->at(ib));
            genBphi1.push_back(tree_puppi.genmeson_phi->at(ib));
        }
     }
     int itr=1;
     for (auto tree:{tree_pfq0,tree_pfq4} ){
        std::vector<float> genBpt,genBeta,genBphi;
        for (int ievt=0; ievt<evt_to_check; ievt++){
           tree.GetEntry(ievt);
           for (int ib=0; ib<tree.genquark_pt->size(); ib++){
              genBpt.push_back(tree.genmeson_pt->at(ib));
              genBeta.push_back(tree.genmeson_eta->at(ib));
              genBphi.push_back(tree.genmeson_phi->at(ib));
           }
        }
        if (genBpt != genBpt1 || genBeta!= genBeta1 || genBphi != genBphi1)
           std::cout<<"tree "<<itr<<" not the same as 0"<<endl;
        itr+=1;
     }
  
  }
  
  
  //////////////////////////////////// Efficiency /////////////////////////////
  // get gen objects and denominator
  std::vector< std::vector< int > > idx_genB_allevts;
  std::vector< std::vector< float > > pt_genB_allevts;

  for (int ievt=0; ievt<tree_puppi.GetEntries(); ievt++){
     tree_puppi.GetEntry(ievt); 
     std::vector<int> idx_genB_thisevt;
     std::vector<float> pt_genB_thisevt;
     // main Loop for B (B1)
     for (int ib=0; ib<tree_puppi.genmeson_pt->size(); ib++){
        if (tree_puppi.genquark_mother->at(ib)<1.0e+06)
           continue;
        if (fabs(tree_puppi.genmeson_eta->at(ib))>2.2) 
           continue;
        if (tree_puppi.genmeson_pt->at(ib)<1)
           continue;
        float minDR_bb=1000;
        // ib2 search for close by B to B1
        for (int ib2=0; ib2<tree_puppi.genmeson_pt->size(); ib2++){
            if (ib==ib2) continue;
            if (tree_puppi.genquark_mother->at(ib2)<1.0e+06)
                continue;   
            if (minDR_bb<deltaR(tree_puppi.genmeson_eta->at(ib2),tree_puppi.genmeson_phi->at(ib2),tree_puppi.genmeson_eta->at(ib),tree_puppi.genmeson_phi->at(ib) ) )
               continue;
            minDR_bb=deltaR(tree_puppi.genmeson_eta->at(ib2),tree_puppi.genmeson_phi->at(ib2),tree_puppi.genmeson_eta->at(ib),tree_puppi.genmeson_phi->at(ib) ); 
       } // end ib2
       hdr_bb->Fill(minDR_bb);
       // remove close by B pairs
       if (minDR_bb<0.8) 
          continue;
       heffpt_den->Fill(tree_puppi.genmeson_pt->at(ib));
       heffeta_den->Fill(tree_puppi.genmeson_eta->at(ib));
       pt_genB_thisevt.push_back(tree_puppi.genmeson_pt->at(ib));
       idx_genB_thisevt.push_back(ib);
    }
    idx_genB_allevts.push_back(idx_genB_thisevt);
    pt_genB_allevts.push_back(pt_genB_thisevt);
  }
  PlotHisto(hdr_bb, "jeff_drbb", "minimum DR between B", false);
  cout<<" gen done "<<endl;
  
  // standard objects
  int itree=0;
  for (auto tree: {&tree_puppi, &tree_pfq0}){
     for (int ievt=0; ievt<tree->GetEntries(); ievt++){
       tree->GetEntry(ievt);
       std::vector<int> igenBs = idx_genB_allevts[ievt];
       std::vector<float> pt_genBs = pt_genB_allevts[ievt];
       int idx_xcheck=0;
       for (auto igenB: igenBs){
          if (CheckGenFilter && pt_genBs[idx_xcheck]!=tree->genmeson_pt->at(igenB) )
             std::cout<<"error in gen filter saved "<<pt_genBs[idx_xcheck]<<" onfly "<<tree_puppi.genmeson_pt->at(igenB)<<endl;
          idx_xcheck+=1; 
          hdr_jets[itree]->Fill(tree->genmeson_jet1_dr->at(igenB));
          if (tree->genmeson_jet1_dr->at(igenB)>0.4)
             continue;
          heff_jetpt_nums[itree]->Fill(tree->genmeson_pt->at(igenB)); 
          heff_jeteta_nums[itree]->Fill(tree->genmeson_eta->at(igenB));
       }
     }
     itree+=1;
  }
  
  PlotHistos(hdr_jets, "jeff_dr_stdjet", "minDR(gen,reco)", {"puppi","PF","Slimmed"}, false, 0, 700);


  double sgn_pf[]={heff_jetpt_nums[1]->Integral()};
  double sgn_puppi[]={heff_jetpt_nums[0]->Integral()};

  
  heff_jetpt_nums[0]->Divide(heffpt_den);
  heff_jetpt_nums[1]->Divide(heffpt_den);
  //heff_jetpt_nums[2]->Divide(heffpt_den);
  PlotHistos(heff_jetpt_nums, "jeff_effpt_stdjet", "pT(gen)", {"Puppi","PF"}, false,0,1.05);
  
  heff_jeteta_nums[0]->Divide(heffeta_den);
  heff_jeteta_nums[1]->Divide(heffeta_den);
  PlotHistos(heff_jeteta_nums, "jeff_effeta_stdjet", "eta(gen)", {"Puppi","PF"}, false,0,1.05);
  cout<<" stdjet efficiency done "<<endl;
  
 
  // track objects
  itree=0;
  for (auto tree: {&tree_puppi, &tree_pfq0, &tree_pfq4}){
     for (int ievt=0; ievt<tree->GetEntries(); ievt++){
       tree->GetEntry(ievt);
       std::vector<int> igenBs = idx_genB_allevts[ievt];
       std::vector<float> pt_genBs = pt_genB_allevts[ievt];
       for (auto igenB: igenBs){
          hdr_qjets[itree]->Fill(tree->genmeson_qjet1_dr->at(igenB));
          if (tree->genmeson_qjet1_dr->at(igenB)>0.4)
             continue;
          heff_qjetpt_nums[itree]->Fill(tree->genmeson_pt->at(igenB));
          heff_qjeteta_nums[itree]->Fill(tree->genmeson_eta->at(igenB));
       }
     }
     itree+=1;
  }
  
  PlotHistos(hdr_qjets, "jeff_dr_qjet", "minDR(gen,reco)", {"PupC","PFC","PFC-q4"}, false, 0, 700);
  
  double sgn_pfc[]={heff_qjetpt_nums[1]->Integral()};
  double sgn_pfcq4[]={heff_qjetpt_nums[2]->Integral()};
  double sgn_pupc[]={heff_qjetpt_nums[0]->Integral()};


  heff_qjetpt_nums[0]->Divide(heffpt_den);
  heff_qjetpt_nums[1]->Divide(heffpt_den);
  heff_qjetpt_nums[2]->Divide(heffpt_den);
  PlotHistos(heff_qjetpt_nums, "jeff_effpt_qjet", "pT(gen)", {"PupC","PFC","PFC-q4"}, false,0,1.05);
  
  heff_qjeteta_nums[0]->Divide(heffeta_den);
  heff_qjeteta_nums[1]->Divide(heffeta_den);
  heff_qjeteta_nums[2]->Divide(heffeta_den);
  PlotHistos(heff_qjeteta_nums, "jeff_effeta_qjet", "eta(gen)", {"PupC","PFC","PFC-q4"}, false,0,1.05);
  
  cout<<" trkjet efficiency done "<<endl;
  
  
  //Combine jet performance
  PlotHistos({heff_jetpt_nums[1],heff_qjetpt_nums[1],heff_qjetpt_nums[2]}, "jeff_eff_pt", "pt(gen)", {"PF","PFC","PFC-q4"}, false,0.6,1.05);
  PlotHistos({heff_jeteta_nums[1],heff_qjeteta_nums[1],heff_qjeteta_nums[2]}, "jeff_eff_eta", "eta(gen)", {"PF","PFC","PFC-q4"}, false,0.6,1.05);
  
  
  /////////////////////////////////// RATE //////////////////////////////////////
  //standard
  itree=0;
  for (auto tree: {&tree_puppi, &tree_pfq0}){
    for (int ievt=0; ievt<tree->GetEntries(); ievt++){
       tree->GetEntry(ievt);
       for (int ij=0; ij<tree->njet; ij++){
          if (fabs(tree->jet_eta->at(ij))>2.5)
             continue;
          bool Skip=false;
          for ( int igen=0; igen<tree->genmeson_pt->size(); igen++ ){
             if (deltaR(tree->genmeson_eta->at(igen),tree->genmeson_phi->at(igen),tree->jet_eta->at(ij),tree->jet_phi->at(ij))>0.4) 
                continue;
             Skip=true;
             break;
          }
          if (Skip) continue;
          hrate_jetpt[itree]->Fill(tree->jet_pt->at(ij));
       }   
    }
    itree+=1;
  }
  
  double rate_pf[]={hrate_jetpt[1]->Integral()};
  double rate_puppi[]={hrate_jetpt[0]->Integral()};


  TH1F* hrate_puppi = PlotIntegrate(hrate_jetpt[0]);
  TH1F* hrate_pf = PlotIntegrate(hrate_jetpt[1]);
  PlotHistos({hrate_puppi, hrate_pf}, "jeff_rate_stdjet", "pT(jet)", {"Puppi","PF"}, true,0,10000000);
  
  cout<<" stdjet rate done "<<endl;
  
  
  itree=0;
  for (auto tree: {&tree_puppi, &tree_pfq0, &tree_pfq4}){
    for (int ievt=0; ievt<tree->GetEntries(); ievt++){
       tree->GetEntry(ievt);
       for (int ij=0; ij<tree->nqjet; ij++){
          if (fabs(tree->qjet_eta->at(ij))>2.5)
             continue;
          bool Skip=false;
          for ( int igen=0; igen<tree->genmeson_pt->size(); igen++ ){
             if (deltaR(tree->genmeson_eta->at(igen),tree->genmeson_phi->at(igen),tree->qjet_eta->at(ij),tree->qjet_phi->at(ij))>0.4) 
                continue;
             Skip=true;
             break;
          }
          if (Skip) continue;
          hrate_qjetpt[itree]->Fill(tree->qjet_pt->at(ij));
       }   
    }
    itree+=1;
  }
  
  double rate_pfc[]={hrate_qjetpt[1]->Integral()};
  double rate_pupc[]={hrate_qjetpt[0]->Integral()};
  double rate_pfcq4[]={hrate_qjetpt[2]->Integral()};


  TH1F* hrate_pupc = PlotIntegrate(hrate_qjetpt[0]);
  TH1F* hrate_pfc = PlotIntegrate(hrate_qjetpt[1]);
  TH1F* hrate_pfcq4 = PlotIntegrate(hrate_qjetpt[2]);
  
  PlotHistos({hrate_pupc, hrate_pfc, hrate_pfcq4}, "jeff_rate_qjet", "pT(jet)", {"PupC","PFC","PFC-q4"}, true,0,10000000);
  cout<<" trkjet rate done "<<endl;
  
  PlotHistos({hrate_pf, hrate_pfc, hrate_pfcq4}, "jeff_rate_jet", "pT(jet)", {"PF","PFC","PFC-q4"}, true,0,10000000);
  
  rate_pfc[0]/=rate_pf[0];
  rate_pupc[0]/=rate_pf[0];
  rate_pfcq4[0]/=rate_pf[0];
  rate_puppi[0]/=rate_pf[0];

  sgn_pfc[0]/=sgn_pf[0];
  sgn_pupc[0]/=sgn_pf[0];
  sgn_pfcq4[0]/=sgn_pf[0];
  sgn_puppi[0]/=sgn_pf[0];

  cout<<" rate pfc "<<rate_pfc[0]<<" pupc "<<rate_pupc[0]<<" pfc-q4 "<<rate_pfcq4[0]<<" puppi "<<rate_puppi[0]<<endl;
  cout<<" eff pfc "<<sgn_pfc[0]<<" pupc "<<sgn_pupc[0]<<" pfc-q4 "<<sgn_pfcq4[0]<<" puppi "<<sgn_puppi[0]<<endl;

  float axis_graph[]={0,1};

  TGraph * gr_pfc = new TGraph(1,rate_pfc,sgn_pfc);
  TGraph * gr_pupc = new TGraph(1,rate_pupc,sgn_pupc);
  TGraph * gr_puppi = new TGraph(1,rate_puppi,sgn_puppi);
  TGraph * gr_pfcq4 = new TGraph(1,rate_pfcq4,sgn_pfcq4);

  TGraph * gr_axis = new TGraph(2,axis_graph,axis_graph);
  
  TCanvas * c1 = new TCanvas("c1","",800,600);

  gr_pfc->SetMarkerColor(1);
  gr_pupc->SetMarkerColor(2);
  gr_puppi->SetMarkerColor(3);
  gr_pfcq4->SetMarkerColor(4);

  gr_pfc->SetMarkerStyle(8);
  gr_pupc->SetMarkerStyle(8);
  gr_puppi->SetMarkerStyle(8);
  gr_pfcq4->SetMarkerStyle(8);

  gr_pfc->SetLineColor(1);
  gr_pupc->SetLineColor(2);
  gr_puppi->SetLineColor(3);
  gr_pfcq4->SetLineColor(4);
  gr_axis->SetLineColor(0);

  gr_axis->Draw("AP");
  gr_pfc->Draw("P sames");
  gr_pupc->Draw("P sames");
  gr_puppi->Draw("P sames");
  gr_pfcq4->Draw("P sames");

  gr_axis->SetTitle(";Rate wrt PF jets ;Efficiency wrt PF jets");
  gr_axis->GetXaxis()->SetRangeUser(0,1);
  gr_axis->GetYaxis()->SetRangeUser(0,1);
    
//  gr_pfc->GetXaxis()->SetMaximum(1);
  c1->SaveAs("jeteff_graph.png");


 

  return 0;
}
