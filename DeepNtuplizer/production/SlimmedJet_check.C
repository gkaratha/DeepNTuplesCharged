#include "jet_tree.h"

int SlimmedJet_check(){


TChain * cc_slimmed = new TChain("genbanalizer/tree");
TChain * cc_puppi = new TChain("genbanalizer/tree");


cc_slimmed->Add("outputSelectedPatPFCvsSlimmed_0.root");
cc_puppi->Add("outputSelectedPatPFCvsPuppiRecluster_0.root");


jet_tree tree_slimmed;
tree_slimmed.Init(cc_slimmed);

jet_tree tree_puppi;
tree_puppi.Init(cc_puppi);

TH1F* hpt_pup = new TH1F("hpt_pup","",100,0,400);
TH1F* hpt_slm = new TH1F("hpt_slm","",100,0,400);
TH1F* heta_pup = new TH1F("heta_pup","",100,-4,4);
TH1F* heta_slm = new TH1F("heta_slm","",100,-4,4);
TH1F* hphi_pup = new TH1F("hphi_pup","",100,-3.5,3.5);
TH1F* hphi_slm = new TH1F("hphi_slm","",100,-3.5,3.5);

TH1F* hnchg_pup = new TH1F("hnchg_pup","",20,0,20);
TH1F* hnchg_slm = new TH1F("hnchg_slm","",20,0,20);

TH1F* hnpf_pup = new TH1F("hnpf_pup","",30,0,30);
TH1F* hnpf_slm = new TH1F("hnpf_slm","",30,0,30);

TH1F* hdr = new TH1F("hdr","",100,0,1);

TH1F* hdpf = new TH1F("hdpf","",26,-5,20);
TH1F* hdchg = new TH1F("hdchg","",12,-2,10);

TH1F* hjppt_pup = new TH1F("hjppt_pup","",100,0,25);
TH1F* hjppt_slm = new TH1F("hjppt_slm","",100,0,25);
TH1F* hjpeta_pup = new TH1F("hjpeta_pup","",100,-4,4);
TH1F* hjpeta_slm = new TH1F("hjpeta_slm","",100,-4,4);
TH1F* hjpphi_pup = new TH1F("hjpphi_pup","",100,-3.5,3.5);
TH1F* hjpphi_slm = new TH1F("hjpphi_slm","",100,-3.5,3.5);
TH1F* hjppdg_pup = new TH1F("hjppdg_pup","",5,0,5);
TH1F* hjppdg_slm = new TH1F("hjppdg_slm","",5,0,5);

TH1F* hjpcharge_pup = new TH1F("hjpchg_pup","",3,-1.5,1.5);
TH1F* hjpcharge_slm = new TH1F("hjpchg_slm","",3,-1.5,1.5);
TH1F* hjppt_charged_pup = new TH1F("hjppt_chg_pup","",100,0,25);
TH1F* hjppt_charged_slm = new TH1F("hjppt_chg_slm","",100,0,25);
TH1F* hjppt_neutral_pup = new TH1F("hjppt_ntr_pup","",100,0,25);
TH1F* hjppt_neutral_slm = new TH1F("hjppt_ntr_slm","",100,0,25);

TH1F* hjppt_mu_pup = new TH1F("hjppt_mu_pup","",100,0,25);
TH1F* hjppt_mu_slm = new TH1F("hjppt_mu_slm","",100,0,25);
TH1F* hjppt_el_pup = new TH1F("hjppt_el_pup","",100,0,25);
TH1F* hjppt_el_slm = new TH1F("hjppt_el_slm","",100,0,25);
TH1F* hjppt_pi_pup = new TH1F("hjppt_pi_pup","",100,0,25);
TH1F* hjppt_pi_slm = new TH1F("hjppt_pi_slm","",100,0,25);
TH1F* hjppt_pi0_pup = new TH1F("hjppt_pi0_pup","",100,0,25);
TH1F* hjppt_pi0_slm = new TH1F("hjppt_pi0_slm","",100,0,25);
TH1F* hjppt_g_pup = new TH1F("hjppt_g_pup","",100,0,25);
TH1F* hjppt_g_slm = new TH1F("hjppt_g_slm","",100,0,25);


TH1F* hpfpt_mu = new TH1F("hpfpt_mu","",100,0,25);
TH1F* hpfpt_el = new TH1F("hpfpt_el","",100,0,25);
TH1F* hpfpt_pi = new TH1F("hpfpt_pi","",100,0,25);
TH1F* hpfpt_pi0 = new TH1F("hpfpt_pi0","",100,0,25);
TH1F* hpfpt_g = new TH1F("hpfpt_g","",100,0,25);

TH1F* hdr2 = new TH1F("hdr2","",100,0,1);


for (int ievt=0; ievt<cc_puppi->GetEntries(); ievt++){
   tree_puppi.GetEntry(ievt);
   tree_slimmed.GetEntry(ievt);

   for (int irjet=0; irjet<tree_puppi.jet_pt->size(); irjet++){
       hpt_pup->Fill(tree_puppi.jet_pt->at(irjet)); 
       heta_pup->Fill(tree_puppi.jet_eta->at(irjet));
       hphi_pup->Fill(tree_puppi.jet_phi->at(irjet));
       hnchg_pup->Fill(tree_puppi.jet_ntrk->at(irjet));
       hnpf_pup->Fill(tree_puppi.jet_npf->at(irjet));
   }

   for (int irjet=0; irjet<tree_slimmed.jet_pt->size(); irjet++){
       hpt_slm->Fill(tree_slimmed.jet_pt->at(irjet));
       heta_slm->Fill(tree_slimmed.jet_eta->at(irjet));
       hphi_slm->Fill(tree_slimmed.jet_phi->at(irjet));
       hnchg_slm->Fill(tree_slimmed.jet_ntrk->at(irjet));
       hnpf_slm->Fill(tree_slimmed.jet_npf->at(irjet));
   }

   std::vector<int> losts_idx; 
   for (int islm=0; islm<tree_slimmed.jet_pt->size(); islm++){
      TLorentzVector vslim;
      vslim.SetPtEtaPhiM(tree_slimmed.jet_pt->at(islm),tree_slimmed.jet_eta->at(islm),tree_slimmed.jet_phi->at(islm),0);
      float minDR=1000;
      int idx=-1;
      for (int ipup=0; ipup<tree_puppi.jet_pt->size(); ipup++){
         TLorentzVector vpup;
         vpup.SetPtEtaPhiM(tree_puppi.jet_pt->at(ipup),tree_puppi.jet_eta->at(ipup),tree_puppi.jet_phi->at(ipup),0);     
         if (minDR<vpup.DeltaR(vslim)) 
            continue;
         minDR=vpup.DeltaR(vslim);
         idx = ipup;
     }
     hdr->Fill(minDR);   
    
     if (minDR>0.1)
        continue;
     hdpf->Fill(tree_slimmed.jet_npf->at(islm)-tree_puppi.jet_npf->at(idx));
     hdchg->Fill(tree_slimmed.jet_ntrk->at(islm)-tree_puppi.jet_ntrk->at(idx));
     //std::cout<<"slimmed pt "<<tree_slimmed.jet_pt->at(islm)<<" eta "<<tree_slimmed.jet_eta->at(islm)<<" phi "<<tree_slimmed.jet_phi->at(islm)<<" npf "<<tree_slimmed.jet_npf->at(islm)<<endl;
     for (int jpart=0; jpart<tree_slimmed.jetpart_pt->size(); jpart++){
         if ( tree_slimmed.jetpart_jetIdx->at(jpart) <islm ) continue;
         if ( tree_slimmed.jetpart_jetIdx->at(jpart) >islm ) break;
         hjppt_slm->Fill(tree_slimmed.jetpart_pt->at(jpart));
         hjpeta_slm->Fill(tree_slimmed.jetpart_eta->at(jpart));
         hjpphi_slm->Fill(tree_slimmed.jetpart_phi->at(jpart));
         hjpcharge_slm->Fill(tree_slimmed.jetpart_charge->at(jpart));
         if (tree_slimmed.jetpart_charge->at(jpart)==0)
             hjppt_neutral_slm->Fill(tree_slimmed.jetpart_pt->at(jpart));
         else
             hjppt_charged_slm->Fill(tree_slimmed.jetpart_pt->at(jpart));
         if (fabs(tree_slimmed.jetpart_pdgId->at(jpart))==11)
            hjppt_el_slm->Fill(tree_slimmed.jetpart_pt->at(jpart));
         if (fabs(tree_slimmed.jetpart_pdgId->at(jpart))==13)
            hjppt_mu_slm->Fill(tree_slimmed.jetpart_pt->at(jpart));
         if (fabs(tree_slimmed.jetpart_pdgId->at(jpart))==22)
            hjppt_g_slm->Fill(tree_slimmed.jetpart_pt->at(jpart));
         if (fabs(tree_slimmed.jetpart_pdgId->at(jpart))==211)
            hjppt_pi_slm->Fill(tree_slimmed.jetpart_pt->at(jpart));
         if (fabs(tree_slimmed.jetpart_pdgId->at(jpart))==130)
            hjppt_pi0_slm->Fill(tree_slimmed.jetpart_pt->at(jpart));
       //  std::cout<<" --daughter "<<jpart<<" pt "<<tree_slimmed.jetpart_pt->at(jpart)<<" eta "<<tree_slimmed.jetpart_eta->at(jpart)<<" phi  "<<tree_slimmed.jetpart_phi->at(jpart)<<" pdg "<<tree_slimmed.jetpart_pdgId->at(jpart)<<endl;
     }
     //std::cout<<"puppi pt "<<tree_puppi.jet_pt->at(idx)<<" eta "<<tree_puppi.jet_eta->at(idx)<<" phi "<<tree_puppi.jet_phi->at(idx)<<" npf "<<tree_puppi.jet_npf->at(idx)<<endl;
     for (int jpart=0; jpart<tree_puppi.jetpart_pt->size(); jpart++){
         if ( tree_puppi.jetpart_jetIdx->at(jpart) <idx ) continue;
         if ( tree_puppi.jetpart_jetIdx->at(jpart) >idx ) break;
         hjppt_pup->Fill(tree_puppi.jetpart_pt->at(jpart));
         hjpeta_pup->Fill(tree_puppi.jetpart_eta->at(jpart));
         hjpphi_pup->Fill(tree_puppi.jetpart_phi->at(jpart));
         hjpcharge_pup->Fill(tree_puppi.jetpart_charge->at(jpart));
         if (tree_puppi.jetpart_charge->at(jpart)==0) 
             hjppt_neutral_pup->Fill(tree_puppi.jetpart_pt->at(jpart));
         else 
             hjppt_charged_pup->Fill(tree_puppi.jetpart_pt->at(jpart));
         if (fabs(tree_puppi.jetpart_pdgId->at(jpart))==11)
            hjppt_el_pup->Fill(tree_puppi.jetpart_pt->at(jpart));
         if (fabs(tree_puppi.jetpart_pdgId->at(jpart))==13)
            hjppt_mu_pup->Fill(tree_puppi.jetpart_pt->at(jpart));
         if (fabs(tree_puppi.jetpart_pdgId->at(jpart))==22)
            hjppt_g_pup->Fill(tree_puppi.jetpart_pt->at(jpart));
         if (fabs(tree_puppi.jetpart_pdgId->at(jpart))==211)
            hjppt_pi_pup->Fill(tree_puppi.jetpart_pt->at(jpart));
         if (fabs(tree_puppi.jetpart_pdgId->at(jpart))==130)
            hjppt_pi0_pup->Fill(tree_puppi.jetpart_pt->at(jpart));
       //  std::cout<<" --daughter "<<jpart<<" pt "<<tree_puppi.jetpart_pt->at(jpart)<<" eta "<<tree_puppi.jetpart_eta->at(jpart)<<" phi "<<tree_puppi.jetpart_phi->at(jpart)<<" pdg "<<tree_puppi.jetpart_pdgId->at(jpart)<<endl;
     }
     for (int jpart=0; jpart<tree_slimmed.jetpart_pt->size(); jpart++){
         if ( tree_slimmed.jetpart_jetIdx->at(jpart) <islm ) continue;
         if ( tree_slimmed.jetpart_jetIdx->at(jpart) >islm ) break;    
         bool Found=false;
         for (int jpart2=0; jpart2<tree_puppi.jetpart_pt->size(); jpart2++){
            if ( tree_puppi.jetpart_jetIdx->at(jpart2) <idx ) continue;
            if ( tree_puppi.jetpart_jetIdx->at(jpart2) >idx ) break;
            if ( tree_slimmed.jetpart_pt->at(jpart) != tree_puppi.jetpart_pt->at(jpart2) || 
                 tree_slimmed.jetpart_eta->at(jpart) != tree_puppi.jetpart_eta->at(jpart2) || 
                tree_slimmed.jetpart_phi->at(jpart) != tree_puppi.jetpart_phi->at(jpart2)  ) 
               continue;
               Found= true;
         }
         if (!Found) 
            losts_idx.push_back(jpart);
     }
 }

//  for( int ik=0; ik<losts_vec.size(); ik++){
  //   cout<<"lost pt "<<losts_vec[ik].Pt()<<" eta "<<losts_vec[ik].Eta()<<" phi "<<losts_vec[ik].Phi()<<endl;
 // }
 //for (int ipf=0; ipf<tree_puppi.pf_pt->size(); ipf++){
   //   cout<<"pf pt "<<tree_puppi.pf_pt->at(ipf)<<" eta "<<tree_puppi.pf_eta->at(ipf)<<" phi "<<tree_puppi.pf_phi->at(ipf)<<" pdg "<<tree_puppi.pf_pdgId->at(ipf)<<endl;
  //}

  for( int ik : losts_idx){
    TLorentzVector lost_vec;
    lost_vec.SetPtEtaPhiM(tree_slimmed.jetpart_pt->at(ik),tree_slimmed.jetpart_eta->at(ik),tree_slimmed.jetpart_phi->at(ik),0);
   // cout<<"lost pt "<<lost_vec.Pt()<<" eta "<<lost_vec.Eta()<<" phi "<<lost_vec.Phi()<<endl;
    float minDR=100;
    for (int ipf=0; ipf<tree_puppi.pf_pt->size(); ipf++){
        TLorentzVector vpf;
        vpf.SetPtEtaPhiM(tree_puppi.pf_pt->at(ipf),tree_puppi.pf_eta->at(ipf),tree_puppi.pf_phi->at(ipf),0); 
        if (minDR<lost_vec.DeltaR(vpf))
            continue;
        minDR=lost_vec.DeltaR(vpf);
     //   cout<<"  --- matched to lost pt "<<tree_puppi.pf_pt->at(ipf)<<" eta "<<tree_puppi.pf_eta->at(ipf)<<" phi "<<tree_puppi.pf_phi->at(ipf)<<endl;
    } 
    hdr2->Fill(minDR);
   /* if (fabs(tree_slimmed.pf_pdgId->at(ipf))==11)
        hpfpt_el->Fill(tree_slimmed.pf_pt->at(ipf));
    if (fabs(tree_slimmed.pf_pdgId->at(ipf))==13)
        hpfpt_mu->Fill(tree_slimmed.pf_pt->at(ipf));
    if (fabs(tree_slimmed.pf_pdgId->at(ipf))==22)
        hpfpt_g->Fill(tree_slimmed.pf_pt->at(ipf));
    if (fabs(tree_slimmed.pf_pdgId->at(ipf))==211)
        hpfpt_pi->Fill(tree_slimmed.pf_pt->at(ipf));
    if (fabs(tree_slimmed.pf_pdgId->at(ipf))==130)
        hpfpt_pi0->Fill(tree_slimmed.pf_pt->at(ipf));*/
 }

}




 TLegend * leg = new TLegend(0.7,0.7,1,1);
 leg->AddEntry(hpt_slm,"Slimmed");
 leg->AddEntry(hpt_pup,"Puppi");


 TCanvas * cpt_reco = new TCanvas("cpt_reco","",800,600);
 hpt_slm->Draw();
 hpt_slm->SetLineColor(2);
 hpt_pup->SetLineColor(1);
 hpt_pup->Draw("sames");
 leg->Draw("sames");
 hpt_slm->GetXaxis()->SetTitle("pT(jet)");
 cpt_reco->SetLogy();
 cpt_reco->SaveAs("xcheck_pt_slim.png");

 TCanvas * ceta_reco = new TCanvas("ceta_reco","",800,600);
 heta_slm->Draw();
 heta_slm->SetLineColor(2);
 heta_pup->SetLineColor(1);
 heta_pup->Draw("sames");
 leg->Draw("sames");
 heta_slm->GetXaxis()->SetTitle("eta(jet)");
 ceta_reco->SetLogy();
 ceta_reco->SaveAs("xcheck_eta_slim.png");

 TCanvas * cphi_reco = new TCanvas("cphi_reco","",800,600);
 hphi_slm->Draw();
 hphi_slm->SetLineColor(2);
 hphi_pup->SetLineColor(1);
 hphi_pup->Draw("sames");
 leg->Draw("sames");
 hphi_slm->GetXaxis()->SetTitle("phi(jet)");
 cphi_reco->SetLogy();
 cphi_reco->SaveAs("xcheck_phi_slim.png");

 TCanvas * cnchg_reco = new TCanvas("cnchg_reco","",800,600);
 hnchg_slm->Draw();
 hnchg_slm->SetLineColor(2);
 hnchg_pup->SetLineColor(1);
 hnchg_pup->Draw("sames");
 leg->Draw("sames");
 hnchg_slm->GetXaxis()->SetTitle("N(CHG)");
 cnchg_reco->SetLogy();
 cnchg_reco->SaveAs("xcheck_nchg_slim.png");


 TCanvas * cnpf_reco = new TCanvas("cnpf_reco","",800,600);
 hnpf_slm->Draw();
 hnpf_slm->SetLineColor(2);
 hnpf_pup->SetLineColor(1);
 hnpf_pup->Draw("sames");
 leg->Draw("sames");
 hnpf_slm->GetXaxis()->SetTitle("N(PF)");
 cnpf_reco->SetLogy();
 cnpf_reco->SaveAs("xcheck_npf_slim.png");

 TCanvas * cdr = new TCanvas("cdr","",800,600);
 hdr->Draw();
 hdr->GetXaxis()->SetTitle("DR");
 cdr->SetLogy();
 cdr->SaveAs("xcheck_dr.png");

 TCanvas * cdtrk = new TCanvas("cdtrk","",800,600);
 hdchg->Draw();
 hdchg->GetXaxis()->SetTitle("Ntrk(Slim)-Ntrk(Puppi)");
 cdtrk->SetLogy();
 cdtrk->SaveAs("xcheck_dtrk.png");

 TCanvas * cdpf = new TCanvas("cdpf","",800,600);
 hdpf->Draw();
 hdpf->GetXaxis()->SetTitle("Npf(Slim)-Npf(Puppi)");
 cdpf->SetLogy();
 cdpf->SaveAs("xcheck_dpf.png");

 TCanvas * cjppt_reco = new TCanvas("cjppt_reco","",800,600);
 hjppt_slm->Draw();
 hjppt_slm->SetLineColor(2);
 hjppt_pup->SetLineColor(1);
 hjppt_pup->Draw("sames");
 leg->Draw("sames");
 hjppt_slm->GetXaxis()->SetTitle("Jet const. pT");
 cjppt_reco->SetLogy();
 cjppt_reco->SaveAs("xcheck_jppt_slim.png");

 TCanvas * cjpeta_reco = new TCanvas("cjpeta_reco","",800,600);
 hjpeta_slm->Draw();
 hjpeta_slm->SetLineColor(2);
 hjpeta_pup->SetLineColor(1);
 hjpeta_pup->Draw("sames");
 leg->Draw("sames");
 hjpeta_slm->GetXaxis()->SetTitle("Jet const. eta");
 cjpeta_reco->SetLogy();
 cjpeta_reco->SaveAs("xcheck_jpeta_slim.png");

 TCanvas * cjpphi_reco = new TCanvas("cjpphi_reco","",800,600);
 hjpphi_slm->Draw();
 hjpphi_slm->SetLineColor(2);
 hjpphi_pup->SetLineColor(1);
 hjpphi_pup->Draw("sames");
 leg->Draw("sames");
 hjpphi_slm->GetXaxis()->SetTitle("Jet const. phi");
 cjpphi_reco->SetLogy();
 cjpphi_reco->SaveAs("xcheck_jpphi_slim.png");


 TCanvas * cjp_charge = new TCanvas("cjp_charge","",800,600);
 hjpcharge_slm->Draw();
 hjpcharge_slm->SetLineColor(2);
 hjpcharge_pup->SetLineColor(1);
 hjpcharge_pup->Draw("sames");
 leg->Draw("sames");
 hjpcharge_slm->GetXaxis()->SetTitle("Jet const. charge");
 //cjp_charge->SetLogy();
 hjpcharge_slm->SetMinimum(0);
 cjp_charge->SaveAs("xcheck_jpcharge_slim.png");


 TCanvas * cjppt_charged = new TCanvas("cjppt_charged","",800,600);
 hjppt_charged_slm->Draw();
 hjppt_charged_slm->SetLineColor(2);
 hjppt_charged_pup->SetLineColor(1);
 hjppt_charged_pup->Draw("sames");
 leg->Draw("sames");
 hjppt_charged_slm->GetXaxis()->SetTitle("Jet charged const. pT");
 cjppt_charged->SetLogy();
 cjppt_charged->SaveAs("xcheck_jppt_charged_slim.png");

 TCanvas * cjppt_neutral = new TCanvas("cjppt_neutral","",800,600);
 hjppt_neutral_slm->Draw();
 hjppt_neutral_slm->SetLineColor(2);
 hjppt_neutral_pup->SetLineColor(1);
 hjppt_neutral_pup->Draw("sames");
 leg->Draw("sames");
 hjppt_neutral_slm->GetXaxis()->SetTitle("Jet neutral const. pT");
 cjppt_neutral->SetLogy();
 cjppt_neutral->SaveAs("xcheck_jppt_neutral_slim.png");


 TCanvas * cjppt_el = new TCanvas("cjppt_el","",800,600);
 hjppt_el_slm->Draw();
 hjppt_el_slm->SetLineColor(2);
 hjppt_el_pup->SetLineColor(1);
 hjppt_el_pup->Draw("sames");
 leg->Draw("sames");
 hjppt_el_slm->GetXaxis()->SetTitle("Jet e const. pT");
 cjppt_el->SetLogy();
 cjppt_el->SaveAs("xcheck_jppt_el_slim.png");

 TCanvas * cjppt_mu = new TCanvas("cjppt_mu","",800,600);
 hjppt_mu_slm->Draw();
 hjppt_mu_slm->SetLineColor(2);
 hjppt_mu_pup->SetLineColor(1);
 hjppt_mu_pup->Draw("sames");
 leg->Draw("sames");
 hjppt_mu_slm->GetXaxis()->SetTitle("Jet mu const. pT");
 cjppt_mu->SetLogy();
 cjppt_mu->SaveAs("xcheck_jppt_mu_slim.png");


 TCanvas * cjppt_g = new TCanvas("cjppt_g","",800,600);
 hjppt_g_slm->Draw();
 hjppt_g_slm->SetLineColor(2);
 hjppt_g_pup->SetLineColor(1);
 hjppt_g_pup->Draw("sames");
 leg->Draw("sames");
 hjppt_g_slm->GetXaxis()->SetTitle("Jet g const. pT");
 cjppt_g->SetLogy();
 cjppt_g->SaveAs("xcheck_jppt_g_slim.png");


 TCanvas * cjppt_pi = new TCanvas("cjppt_pi","",800,600);
 hjppt_pi_slm->Draw();
 hjppt_pi_slm->SetLineColor(2);
 hjppt_pi_pup->SetLineColor(1);
 hjppt_pi_pup->Draw("sames");
 leg->Draw("sames");
 hjppt_pi_slm->GetXaxis()->SetTitle("Jet pi const. pT");
 cjppt_pi->SetLogy();
 cjppt_pi->SaveAs("xcheck_jppt_pi_slim.png");

 TCanvas * cjppt_pi0 = new TCanvas("cjppt_pi0","",800,600);
 hjppt_pi0_slm->Draw();
 hjppt_pi0_slm->SetLineColor(2);
 hjppt_pi0_pup->SetLineColor(1);
 hjppt_pi0_pup->Draw("sames");
 leg->Draw("sames");
 hjppt_pi0_slm->GetXaxis()->SetTitle("Jet pi0 const. pT");
 cjppt_pi0->SetLogy();
 cjppt_pi0->SaveAs("xcheck_jppt_pi0_slim.png");

 TCanvas * cpfpt_mu = new TCanvas("cpfpt_mu","",800,600);
 hpfpt_mu->Draw();
 cpfpt_mu->SetLogy();
 cpfpt_mu->SaveAs("xcheck_pfpt_mu_slim.png");

 TCanvas * cpfpt_el = new TCanvas("cpfpt_el","",800,600);
 hpfpt_el->Draw();
 cpfpt_el->SetLogy();
 cpfpt_el->SaveAs("xcheck_pfpt_el_slim.png");

 TCanvas * cpfpt_g = new TCanvas("cpfpt_g","",800,600);
 hpfpt_g->Draw();
 cpfpt_g->SetLogy();
 cpfpt_g->SaveAs("xcheck_pfpt_g_slim.png");

 TCanvas * cpfpt_pi = new TCanvas("cpfpt_pi","",800,600);
 hpfpt_pi->Draw();
 cpfpt_pi->SetLogy();
 cpfpt_pi->SaveAs("xcheck_pfpt_pi_slim.png");

 TCanvas * cpfpt_pi0 = new TCanvas("cpfpt_pi0","",800,600);
 hpfpt_pi0->Draw();
 cpfpt_pi0->SetLogy();
 cpfpt_pi0->SaveAs("xcheck_pfpt_pi0_slim.png");


TCanvas * cdr2 = new TCanvas("cdr2","",800,600);
 hdr2->Draw();
 hdr2->GetXaxis()->SetTitle("DR");
 cdr2->SetLogy();
 cdr2->SaveAs("xcheck_dr2.png");


 return 0;
}
