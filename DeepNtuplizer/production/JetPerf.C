#include "JetEff_helper.h"
#include "jet_tree.h"




int JetPerf(){

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  
  
  TChain * cc = new TChain("genbanalizer/tree");
  cc->Add("outputSignalBkgPVQual4PFCvsPFRecluster_0.root");
  
 
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
 
  TH1F* hpt_b = new TH1F("hpt_b","",50,0,10);
  TH1F* hpt_nob = new TH1F("hpt_nob","",50,0,10);
  
  TH2F* hfrac_nd = new TH2F("hfrac_nd","",10,0,1,100,0,100);





  
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

      std::vector<int> matched_tracks; 
      std::vector<int> matched_trackB;
      std::vector<int> matched_tracks_gen;


    
      for (int ib=0; ib<B_daughters_pt.size(); ib++){
         if (B_daughters_pt[ib].size()>0) {
            hgen_num->Fill(B_daughters_pt[ib].size());
            heffpt_den->Fill(tree.genpart_Bmeson_pt->at(ib));
            heffnd_den->Fill(B_daughters_pt[ib].size());
         }
         int nfound=0;
         for (int id=0; id<B_daughters_pt[ib].size(); id++){
            float minDR=1000;
            int index=-1;
            for (int itrk=0; itrk<tree.qjetpart_qjetIdx->size(); itrk++){
               if (minDR < deltaR(B_daughters_eta[ib][id],B_daughters_phi[ib][id],tree.qjetpart_eta->at(itrk),tree.qjetpart_phi->at(itrk)) )
                  continue;
               minDR = deltaR(B_daughters_eta[ib][id],B_daughters_phi[ib][id],tree.qjetpart_eta->at(itrk),tree.qjetpart_phi->at(itrk));
               index = itrk;
            }
	    hdr->Fill(minDR);
            if (minDR>0.03) 
               continue;
             nfound+=1;
             matched_tracks.push_back(index);
             matched_trackB.push_back(ib);
             matched_tracks_gen.push_back(id);
         }
         hmatched_num->Fill(nfound);
         hmatched_2d->Fill(B_daughters_pt[ib].size(),nfound);
         if (tree.genpart_Bmeson_pt->at(ib)<5)
            hmatched_2d_bellow5->Fill(B_daughters_pt[ib].size(),nfound);
         else if (tree.genpart_Bmeson_pt->at(ib)<10)
            hmatched_2d_bellow10->Fill(B_daughters_pt[ib].size(),nfound);
         else
            hmatched_2d_above10->Fill(B_daughters_pt[ib].size(),nfound);
         if (nfound>0){
            heffpt_num->Fill(tree.genpart_Bmeson_pt->at(ib));
            heffnd_num->Fill(B_daughters_pt[ib].size());
         }
      }
      
      std::vector<std::vector<int>> Bconst;
      std::vector<int> temp_Bconst;
      std::vector<std::vector<int>> Bconst_matched;
      std::vector<int> temp_Bconst_matched;
      int last_jet=-1;
      for (int itrk=0; itrk<tree.qjetpart_qjetIdx->size(); itrk++){
          if (temp_Bconst.size()>0 && last_jet != tree.qjetpart_qjetIdx->at(itrk) ){
             Bconst.push_back(temp_Bconst);
             temp_Bconst.clear();  
             Bconst_matched.push_back(temp_Bconst_matched);
             temp_Bconst_matched.clear();
          }
          temp_Bconst.push_back(itrk);
          last_jet=tree.qjetpart_qjetIdx->at(itrk);
          if ( std::find(matched_tracks.begin(),matched_tracks.end(),itrk) != matched_tracks.end() ){
              hpt_b->Fill(tree.qjetpart_pt->at(itrk));
              temp_Bconst_matched.push_back(itrk);
          }
          else
             hpt_nob->Fill(tree.qjetpart_pt->at(itrk));  
      } 

      for (int ib=0; ib<Bconst.size(); ib++){
          hfrac_nd->Fill((Bconst_matched[ib].size()*1.0)/Bconst[ib].size(),Bconst[ib].size());
      }
  }

  heffpt_num->Divide(heffpt_den);
  heffnd_num->Divide(heffnd_den);

  TCanvas * cndau = new TCanvas("cndau","c1",800,600);
  hgen_num->Draw();
  hmatched_num->SetLineColor(2);
  hmatched_num->Draw("sames");
  cndau->SetLogy();
  cndau->SaveAs("jetperf_Ndaughters.png");
  
  TCanvas * c2 = new TCanvas("c2","c1",800,600);
  hdr->Draw();
  c2->SetLogy();
  c2->SaveAs("jetperf_mindr.png");

  TCanvas * cndau2d = new TCanvas("cndau2d","c1",800,600);
  hmatched_2d->Draw("COLZ");
  cndau2d->SetLogz();
  cndau2d->SaveAs("jetperf_Ndaughters2d.png");
  
  TCanvas * cndau2d_bellow5 = new TCanvas("cndau2d_bellow5","c1",800,600);
  hmatched_2d_bellow5->Draw("COLZ");
  cndau2d_bellow5->SetLogz();
  cndau2d_bellow5->SaveAs("jetperf_Ndaughters2d_bellow5.png");

  TCanvas * cndau2d_bellow10 = new TCanvas("cndau2d_bellow10","c1",800,600);
  hmatched_2d_bellow10->Draw("COLZ");
  cndau2d_bellow10->SetLogz();
  cndau2d_bellow10->SaveAs("jetperf_Ndaughters2d_bellow10.png");

  TCanvas * cndau2d_above10 = new TCanvas("cndau2d_above10","c1",800,600);
  hmatched_2d_above10->Draw("COLZ");
  cndau2d_above10->SetLogz();
  cndau2d_above10->SaveAs("jetperf_Ndaughters2d_above10.png");

  TCanvas * ceffpt = new TCanvas("ceffpt","c1",800,600);
  heffpt_num->Draw();
  ceffpt->SaveAs("jetperf_effpt.png");

  TCanvas * ceffnd = new TCanvas("ceffnd","c1",800,600);
  heffnd_num->Draw();
  ceffnd->SaveAs("jetperf_effnd.png");
  
  TCanvas * cpt = new TCanvas("cpt","c1",800,600);
  hpt_nob->Scale(1./hpt_nob->Integral());
  hpt_b->Scale(1./hpt_b->Integral());
  hpt_b->Draw("HIST");
  hpt_nob->SetLineColor(2);
  hpt_nob->Draw("HIST sames");
  cpt->SetLogy();
  cpt->SaveAs("jetperf_pt_daughters.png");


  TCanvas * cfrac2d = new TCanvas("cfrac2d","c1",800,600);
  hfrac_nd->Draw("COLZ");
  cfrac2d->SetLogz();
  cfrac2d->SaveAs("jetperf_frac2d.png");


  return 0;
}
