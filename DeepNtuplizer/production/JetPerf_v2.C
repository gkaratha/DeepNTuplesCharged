#include "JetEff_helper.h"
#include "jet_tree.h"




int JetPerf_v2(){

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  
  
  TChain * cc = new TChain("genbanalizer/tree");
//  cc->Add("outputSignalBkgPVQual4PFCvsPFRecluster_0.root");
  cc->Add("outputSignalBkgPVQual4_pt0p75PFCvsPFRecluster_0.root");
  
 
  jet_tree tree;
  tree.Init(cc);
  
  
  TH1F* hdr = new TH1F("hdr","",100,0,3); 
  TH1F* hgen_num = new TH1F("hgen_num","",20,0,20);
  TH1F* hmatched_num = new TH1F("hmatched_num","",20,0,20);
  TH2F* hmatched_2d = new TH2F("hmatched_2d","",10,0,10,10,0,10);  
  TH2F* hmatched_2d_bellow5 = new TH2F("hmatched_2d_below5","",10,0,10,10,0,10);
  TH2F* hmatched_2d_bellow10 = new TH2F("hmatched_2d_below10","",10,0,10,10,0,10);
  TH2F* hmatched_2d_above10 = new TH2F("hmatched_2d_above10","",10,0,10,10,0,10);
 
  TH1F* heffpt_den = new TH1F("heffpt_den","",50,0,100);
  TH1F* heffpt_num = new TH1F("heffpt_num","",50,0,100);

  TH1F* heffnd_den = new TH1F("heffnd_den","",10,0,10);
  TH1F* heffnd_num = new TH1F("heffnd_num","",10,0,10);
 
  TH1F* hpt_b = new TH1F("hpt_b","",100,0,20);
  TH1F* hpt_nob = new TH1F("hpt_nob","",100,0,20);
  
  TH2F* hfrac_nd = new TH2F("hfrac_nd","",10,0,1,100,0,100);
  TH2F* hfrac_nd_bellow4 = new TH2F("hfrac_nd_bellow4","",10,0,1,100,0,100);
  TH2F* hfrac_nd_above3 = new TH2F("hfrac_nd_above3","",10,0,1,100,0,100);


  TH2F* heff_nd_pt_den = new TH2F("heff_nd_pt_den","",10,0,10,50,0,50);
  TH2F* heff_nd_pt_num = new TH2F("heff_nd_pt_num","",10,0,10,50,0,50);


  
  for (int ievt=0; ievt<cc->GetEntries(); ievt++){
     tree.GetEntry(ievt);
    // if (ievt>0) break; 
     if (ievt%1000==0) 
        std::cout<<"-- ievt "<<ievt<<std::endl;
     
     std::vector<float> part_to_B= connectPartToB(tree.genpart_Bmeson,tree.genpart_Bmeson_pt,tree.genpart_Bmeson_eta,tree.genpart_Bmeson_phi);
     
     std::vector<std::vector<std::vector<float>> >  Bdecay_daughters = CreateBdecayChain(part_to_B, tree.genpart_Bmeson_pt->size(), tree.genpart_pt,tree.genpart_eta, tree.genpart_phi, tree.genpart_pdgId );
      std::vector< std::vector<float> > B_daughters_pt = Bdecay_daughters[0];
      std::vector< std::vector<float> > B_daughters_eta = Bdecay_daughters[1];
      std::vector< std::vector<float> > B_daughters_phi = Bdecay_daughters[2]; 
      std::vector< std::vector<float> > B_daughters_pdgId = Bdecay_daughters[3];

      std::vector<std::vector< std::vector<int> >> B_daughters_match =  MatchBdaughtersToTracks(B_daughters_eta, B_daughters_phi, tree.qjetpart_eta, tree.qjetpart_phi );
      std::vector<std::vector<int>>  B_daughter_matchedTrkIdx = B_daughters_match[0];
      std::vector<std::vector<int>>  B_daughter_matchedDaughterIdx = B_daughters_match[1];
   
    
      for (int ib=0; ib<B_daughters_pt.size(); ib++){
         if (B_daughters_pt[ib].size()==0) 
            continue;
   //      if (tree.genpart_Bmeson_pt->at(ib)<5) 
     //       continue;
         heffpt_den->Fill(tree.genpart_Bmeson_pt->at(ib));
         int nfound = B_daughter_matchedTrkIdx[ib].size();
         if (tree.genpart_Bmeson_pt->at(ib)>5){
            hgen_num->Fill(B_daughters_pt[ib].size());
            hmatched_num->Fill(nfound);
            hmatched_2d->Fill(B_daughters_pt[ib].size(),nfound);
            heffnd_den->Fill(B_daughters_pt[ib].size());
            if (nfound>0) heffnd_num->Fill(B_daughters_pt[ib].size());
         }
         heff_nd_pt_den->Fill(B_daughters_pt[ib].size(),tree.genpart_Bmeson_pt->at(ib));
         if (tree.genpart_Bmeson_pt->at(ib)<10 && tree.genpart_Bmeson_pt->at(ib)>5)
            hmatched_2d_bellow5->Fill(B_daughters_pt[ib].size(),nfound);
         else if (tree.genpart_Bmeson_pt->at(ib)<20 && tree.genpart_Bmeson_pt->at(ib)>10)
            hmatched_2d_bellow10->Fill(B_daughters_pt[ib].size(),nfound);
         else if (tree.genpart_Bmeson_pt->at(ib)>20)
            hmatched_2d_above10->Fill(B_daughters_pt[ib].size(),nfound);
         if (nfound>0){
            heffpt_num->Fill(tree.genpart_Bmeson_pt->at(ib));
            heff_nd_pt_num->Fill(B_daughters_pt[ib].size(),tree.genpart_Bmeson_pt->at(ib));
         }
      }
      
      std::vector<std::vector<int>> Bconst=OrderConstituents(tree.qjetpart_qjetIdx);
      std::vector< std::vector< std::vector<int> >> Bconst_matched = MatchJetConstToBdaughters( Bconst,  B_daughter_matchedTrkIdx, B_daughter_matchedDaughterIdx);
      std::vector<std::vector<int>> Bconst_matched_genB = Bconst_matched[0];
      std::vector<std::vector<int>> Bconst_matched_trackIdx = Bconst_matched[1];
      std::vector<std::vector<int>> Bconst_matched_daughterIdx= Bconst_matched[2];

      for (int ib=0; ib<Bconst.size(); ib++){
        //if (Bconst_matched_trackIdx[ib].size()<2)
        //   continue;
        hfrac_nd->Fill((Bconst_matched_trackIdx[ib].size()*1.0)/Bconst[ib].size(),Bconst[ib].size());
        if (Bconst_matched_trackIdx[ib].size()<4)
           hfrac_nd_bellow4->Fill((Bconst_matched_trackIdx[ib].size()*1.0)/Bconst[ib].size(),Bconst[ib].size());
        if (Bconst_matched_trackIdx[ib].size()>3)
           hfrac_nd_above3->Fill((Bconst_matched_trackIdx[ib].size()*1.0)/Bconst[ib].size(),Bconst[ib].size());
        for(int itrk=0; itrk<Bconst[ib].size(); itrk++){ 
           if ( std::find(Bconst_matched_trackIdx[ib].begin(), Bconst_matched_trackIdx[ib].end(), Bconst[ib][itrk]) == Bconst_matched_trackIdx[ib].end() )
               hpt_nob->Fill(tree.qjetpart_pt->at(Bconst[ib][itrk]));
           else
               hpt_b->Fill(tree.qjetpart_pt->at(Bconst[ib][itrk]));
        }
      }
  }

  heffpt_num->Divide(heffpt_den);
  heffnd_num->Divide(heffnd_den);
  heff_nd_pt_num->Divide(heff_nd_pt_den);
  

  TCanvas * cndau = new TCanvas("cndau","c1",800,600);
  hgen_num->Draw();
  hmatched_num->SetLineColor(2);
  hmatched_num->Draw("sames");
  hgen_num->GetXaxis()->SetTitle("N_{daughter}");
  hgen_num->GetYaxis()->SetTitle("Events");
  cndau->SetLogy();
  cndau->SaveAs("jetperf_Ndaughters.png");
  
  TCanvas * c2 = new TCanvas("c2","c1",800,600);
  hdr->Draw();
  c2->SetLogy();
  c2->SaveAs("jetperf_mindr.png");

  TCanvas * cndau2d = new TCanvas("cndau2d","c1",800,600);
  hmatched_2d->Draw("COLZ");
  hmatched_2d->GetXaxis()->SetTitle("N_{daughter}(gen)");
  hmatched_2d->GetYaxis()->SetTitle("N_{daughter}(matched)");
  cndau2d->SetLogz();
  cndau2d->SaveAs("jetperf_Ndaughters2d.png");
  
  TCanvas * cndau2d_bellow5 = new TCanvas("cndau2d_bellow5","c1",800,600);
  cndau2d_bellow5->SetRightMargin(0.15);
  hmatched_2d_bellow5->Draw("COLZ");
  hmatched_2d_bellow5->GetXaxis()->SetTitle("N_{daughter}(gen)");
  hmatched_2d_bellow5->GetYaxis()->SetTitle("N_{daughter}(matched)");
  cndau2d_bellow5->SetLogz();
  cndau2d_bellow5->SaveAs("jetperf_Ndaughters2d_bellow10.png");

  TCanvas * cndau2d_bellow10 = new TCanvas("cndau2d_bellow10","c1",800,600);
  cndau2d_bellow10->SetRightMargin(0.15);
  hmatched_2d_bellow10->Draw("COLZ");
  hmatched_2d_bellow10->GetXaxis()->SetTitle("N_{daughter}(gen)");
  hmatched_2d_bellow10->GetYaxis()->SetTitle("N_{daughter}(matched)");
  cndau2d_bellow10->SetLogz();
  cndau2d_bellow10->SaveAs("jetperf_Ndaughters2d_bellow20.png");

  TCanvas * cndau2d_above10 = new TCanvas("cndau2d_above10","c1",800,600);
  cndau2d_above10->SetRightMargin(0.15);
  hmatched_2d_above10->Draw("COLZ");
  hmatched_2d_above10->GetXaxis()->SetTitle("N_{daughter}(gen)");
  hmatched_2d_above10->GetYaxis()->SetTitle("N_{daughter}(matched)");
  cndau2d_above10->SetLogz();
  cndau2d_above10->SaveAs("jetperf_Ndaughters2d_above20.png");

  TCanvas * ceffpt = new TCanvas("ceffpt","c1",800,600);
  heffpt_num->Draw();
  heffpt_num->GetXaxis()->SetTitle("p_{T}(B)");
  heffpt_num->GetYaxis()->SetTitle("Efficiency");
  ceffpt->SaveAs("jetperf_effpt.png");

  TCanvas * ceffnd = new TCanvas("ceffnd","c1",800,600);
  heffnd_num->Draw();
  heffnd_num->GetXaxis()->SetTitle("N_{daughter}(gen)");
  heffnd_num->GetYaxis()->SetTitle("Efficiency");
  ceffnd->SaveAs("jetperf_effnd.png");
  
  TCanvas * cpt = new TCanvas("cpt","c1",800,600);
  hpt_nob->Scale(1./hpt_nob->Integral());
  hpt_b->Scale(1./hpt_b->Integral());
  hpt_b->Draw("HIST");
  hpt_nob->SetLineColor(2);
  hpt_nob->Draw("HIST sames");
  hpt_b->GetYaxis()->SetTitle("Events");
  hpt_b->GetXaxis()->SetTitle("p_{T} (trk)");
  cpt->SetLogy();
  cpt->SaveAs("jetperf_pt_daughters.png");


  TCanvas * cfrac2d = new TCanvas("cfrac2d","c1",800,600);
  cfrac2d->SetRightMargin(0.15);
  hfrac_nd->Draw("COLZ");
  hfrac_nd->GetYaxis()->SetTitle("p_{T}(B)");
  hfrac_nd->GetXaxis()->SetTitle("N_{daughter}(matched)/N_{daughter}(reco)");
  cfrac2d->SetLogz();
  cfrac2d->SaveAs("jetperf_frac2d.png");

  TCanvas * ceff2d = new TCanvas("ceff2d","c1",800,600);
  ceff2d->SetRightMargin(0.15);
  heff_nd_pt_num->Draw("COLZ");
  heff_nd_pt_num->GetXaxis()->SetTitle("N_{daughters}(Gen)");
  heff_nd_pt_num->GetYaxis()->SetTitle("p_{T}(B)");
  ceff2d->SaveAs("jetperf_eff2d.png");  


  TCanvas * cfrac2d_bellow4 = new TCanvas("cfrac2d_bellow4","c1",800,600);
  cfrac2d_bellow4->SetRightMargin(0.15);
  hfrac_nd_bellow4->Draw("COLZ");
  hfrac_nd_bellow4->GetYaxis()->SetTitle("p_{T}(B)");
  hfrac_nd_bellow4->GetXaxis()->SetTitle("N_{daughter}(matched)/N_{daughter}(reco)");
  cfrac2d_bellow4->SetLogz();
  cfrac2d_bellow4->SaveAs("jetperf_frac2d_bellow4.png");

  TCanvas * cfrac2d_above3 = new TCanvas("cfrac2d_above3","c1",800,600);
  cfrac2d_above3->SetRightMargin(0.15);
  hfrac_nd_above3->Draw("COLZ");
  hfrac_nd_above3->GetYaxis()->SetTitle("p_{T}(B)");
  hfrac_nd_above3->GetXaxis()->SetTitle("N_{daughter}(matched)/N_{daughter}(reco)");
  cfrac2d_above3->SetLogz();
  cfrac2d_above3->SaveAs("jetperf_frac2d_above3.png");



  return 0;
}
