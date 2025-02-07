#include "jet_tree.h"

int JetEff(){


TChain * cc = new TChain("genbanalizer/tree");

cc->Add("output_0.root");

jet_tree tree;
tree.Init(cc);

TH1F* hdr_bb = new TH1F("hdr_bb","",100,0,3); 
TH1F* hdr_jet = new TH1F("hdr_jet","",100,0,1);
TH1F* hdr_qjet = new TH1F("hdr_qjet","",100,0,1);
TH1F* hdr_sv = new TH1F("hdr_sv","",100,0,1);
TH1F* heff_den = new TH1F("heff_den","",100,0,30);
TH1F* heff_jet = new TH1F("heff_jet","",100,0,30);
TH1F* heff_qjet = new TH1F("heff_qjet","",100,0,30);
TH1F* heff_sv = new TH1F("heff_sv","",100,0,30);

TH1F* hrate_jet = new TH1F("hrate_jet","",100,0,30);
TH1F* hrate_qjet = new TH1F("hrate_qjet","",100,0,30);
TH1F* hrate_sv = new TH1F("hrate_sv","",100,0,30);
TH1F* hpt_jet = new TH1F("hpt_jet","",100,0,30);
TH1F* hpt_qjet = new TH1F("hpt_qjet","",100,0,30);
TH1F* hpt_sv = new TH1F("hpt_sv","",100,0,30);

float avgpt_jet[30]={0},avgpt_qjet[30]={0},avgpt_sv[30]={0};
float avgpt_njet[30]={0},avgpt_nqjet[30]={0},avgpt_nsv[30]={0};

bool CheckMatch=true;


for (int ievt=0; ievt<cc->GetEntries(); ievt++){
   tree.GetEntry(ievt);
   if (ievt%1000==0) std::cout<<"ievt "<<ievt<<std::endl;
   std::vector<int> used;
   for (int ib=0; ib<tree.genquark_pt->size(); ib++){
      if (tree.genquark_mother->at(ib)<1.0e+06)
         continue;
      TLorentzVector vb1;
      vb1.SetPtEtaPhiM(tree.genmeson_pt->at(ib), tree.genmeson_eta->at(ib), tree.genmeson_phi->at(ib), 5.28);
      if (fabs(vb1.Eta())>2.2) 
         continue;

      float minDR_bb=1000;
      for (int ib2=0; ib2<tree.genquark_pt->size(); ib2++){
          if (ib==ib2) continue;
          if (tree.genquark_mother->at(ib2)<1.0e+06)
            continue;	 
         TLorentzVector vb2;
         vb2.SetPtEtaPhiM(tree.genmeson_pt->at(ib2), tree.genmeson_eta->at(ib2), tree.genmeson_phi->at(ib2), 5.28);
//	 if (vb1==vb2) continue;
         
	 if (minDR_bb>vb2.DeltaR(vb1))
            minDR_bb=vb2.DeltaR(vb1);
	  if (vb1==vb2){
	     cout<<"ib1 "<<ib<<"  "<<vb1.Pt()<<"  "<<vb1.Eta()<<"  "<<vb1.Phi()<<endl;
	     cout<<"ib2 "<<ib2<<"  "<<vb2.Pt()<<"  "<<vb2.Eta()<<"  "<<vb2.Phi()<<endl;
	  }
      }
      hdr_bb->Fill(minDR_bb);
      if (minDR_bb<0.8) 
         continue;
   
      if (CheckMatch){
        float minDR_sv=1000;
        int index_sv=-1;
        for (int ir=0; ir<tree.sv_pt->size(); ir++){
           TLorentzVector rvec;
           rvec.SetPtEtaPhiM(tree.sv_pt->at(ir),tree.sv_eta->at(ir),tree.sv_phi->at(ir),0);   
           if (minDR_sv>rvec.DeltaR(vb1) ){
               minDR_sv=rvec.DeltaR(vb1);
               index_sv=ir;
            }
        }
        if (tree.genmeson_sv1->at(ib)!=index_sv && minDR_sv<1.0)
          std::cout<<" SV: ntuple "<<tree.genmeson_sv1_dr->at(ib)<<" this "<<minDR_sv<<endl;

        float minDR_jet=1000;
        int index_jet=-1;
        for (int ir=0; ir<tree.jet_pt->size(); ir++){
           TLorentzVector rvec;
           rvec.SetPtEtaPhiM(tree.jet_pt->at(ir),tree.jet_eta->at(ir),tree.jet_phi->at(ir),0);
           if (minDR_jet>rvec.DeltaR(vb1) ){
               minDR_jet=rvec.DeltaR(vb1);
               index_jet=ir;
           }
        }
        if (tree.genmeson_jet1->at(ib)!=index_jet && minDR_jet<1.0)
             std::cout<<" Jet: ntuple "<<tree.genmeson_jet1_dr->at(ib)<<" this "<<minDR_jet<<endl;

        float minDR_qjet=1000;
        int index_qjet=-1;
        for (int ir=0; ir<tree.qjet_pt->size(); ir++){
           TLorentzVector rvec;
           rvec.SetPtEtaPhiM(tree.qjet_pt->at(ir),tree.qjet_eta->at(ir),tree.qjet_phi->at(ir),0);
           if (minDR_qjet>rvec.DeltaR(vb1) ){
               minDR_qjet=rvec.DeltaR(vb1);
               index_qjet=ir;
           }
        }
        if (tree.genmeson_qjet1->at(ib)!=index_qjet && minDR_qjet<1.0)
            std::cout<<" qjet: ntuple "<<tree.genmeson_qjet1_dr->at(ib)<<" this "<<minDR_qjet<<endl;
     }

      hdr_qjet -> Fill(tree.genmeson_qjet1_dr->at(ib));
      hdr_jet -> Fill(tree.genmeson_jet1_dr->at(ib));
      hdr_sv -> Fill(tree.genmeson_sv1_dr->at(ib));
      heff_den -> Fill(tree.genmeson_pt->at(ib));
      if (tree.genmeson_sv1_dr->at(ib)<0.4){
 	  heff_sv -> Fill(tree.genmeson_pt->at(ib));   
          hpt_sv-> Fill(tree.sv_pt->at(tree.genmeson_sv1->at(ib))); 
          if (std::find(used.begin(),used.end(),tree.genmeson_sv1->at(ib)) != used.end())
             cout<<"skt "<<tree.genmeson_sv1->at(ib)<<endl;
          used.push_back(tree.genmeson_sv1->at(ib));
      }
      if (tree.genmeson_jet1_dr->at(ib)<0.4){
          heff_jet -> Fill(tree.genmeson_pt->at(ib));
          hpt_jet-> Fill(tree.jet_pt->at(tree.genmeson_jet1->at(ib)));
      }
      if (tree.genmeson_qjet1_dr->at(ib)<0.4){
          heff_qjet -> Fill(tree.genmeson_pt->at(ib));     
          hpt_qjet-> Fill(tree.qjet_pt->at(tree.genmeson_qjet1->at(ib)));
       }

       for (int ik=0; ik<30; ik++){
          if ( tree.genmeson_pt->at(ib)<ik || (ik+1)<tree.genmeson_pt->at(ib))
             continue;
          if (tree.genmeson_sv1_dr->at(ib)<0.4) {
             avgpt_nsv[ik]+=1;
             avgpt_sv[ik]+=tree.sv_pt->at(tree.genmeson_sv1->at(ib));
          }
          if (tree.genmeson_jet1_dr->at(ib)<0.4) {
             avgpt_njet[ik]+=1;
             avgpt_jet[ik]+=tree.jet_pt->at(tree.genmeson_jet1->at(ib));
          }
          if (tree.genmeson_qjet1_dr->at(ib)<0.4) {
             avgpt_nqjet[ik]+=1;
             avgpt_qjet[ik]+=tree.qjet_pt->at(tree.genmeson_qjet1->at(ib));
          }
       }
   }


   /// rate
   for (int ir=0; ir<tree.qjet_pt->size(); ir++){
     TLorentzVector rvec;
     rvec.SetPtEtaPhiM(tree.qjet_pt->at(ir),tree.qjet_eta->at(ir),tree.qjet_phi->at(ir),0);
     bool Skip=false;
     for (int ib=0; ib<tree.genmeson_pt->size(); ib++){
        TLorentzVector gvec;
        gvec.SetPtEtaPhiM(tree.genmeson_pt->at(ib), tree.genmeson_eta->at(ib), tree.genmeson_phi->at(ib), 5.28);
	if (gvec.DeltaR(rvec)>0.4)
           continue;		
        Skip=true;		
	break;
     }   
     if (Skip) continue;
     hrate_qjet->Fill(rvec.Pt());
   }


   for (int ir=0; ir<tree.jet_pt->size(); ir++){
     TLorentzVector rvec;
     rvec.SetPtEtaPhiM(tree.jet_pt->at(ir),tree.jet_eta->at(ir),tree.jet_phi->at(ir),0);
     bool Skip=false;
     for (int ib=0; ib<tree.genmeson_pt->size(); ib++){
        TLorentzVector gvec;
        gvec.SetPtEtaPhiM(tree.genmeson_pt->at(ib), tree.genmeson_eta->at(ib), tree.genmeson_phi->at(ib), 5.28);
	if (gvec.DeltaR(rvec)>0.4)
           continue;		
        Skip=true;		
	break;
     }   
     if (Skip) continue;
     hrate_jet->Fill(rvec.Pt());
   }

   for (int ir=0; ir<tree.sv_pt->size(); ir++){
     TLorentzVector rvec;
     rvec.SetPtEtaPhiM(tree.sv_pt->at(ir),tree.sv_eta->at(ir),tree.sv_phi->at(ir),0);
     bool Skip=false;
     for (int ib=0; ib<tree.genmeson_pt->size(); ib++){
        TLorentzVector gvec;
        gvec.SetPtEtaPhiM(tree.genmeson_pt->at(ib), tree.genmeson_eta->at(ib), tree.genmeson_phi->at(ib), 5.28);
	if (gvec.DeltaR(rvec)>0.4)
           continue;		
        Skip=true;		
	break;
     }   
     if (Skip) continue;
     hrate_sv->Fill(rvec.Pt());
   }

}


 cout<<" pf jet: eff "<<hpt_jet->Integral(0,101)/heff_den->Integral(0,101)<<" rate "<<hrate_jet->Integral(0,101)<<endl;
 cout<<" sv: eff "<<hpt_sv->Integral(0,101)/heff_den->Integral(0,101)<<" rate "<<hrate_sv->Integral(0,101)<<endl;
 cout<<" pfc jet: eff "<<hpt_qjet->Integral(0,101)/heff_den->Integral(0,101)<<" rate "<<hrate_qjet->Integral(0,101)<<endl;

 TH1F* hcrate_sv = (TH1F*) hrate_sv->Clone();
 TH1F* hcrate_jet = (TH1F*) hrate_jet->Clone();
 TH1F* hcrate_qjet = (TH1F*) hrate_qjet->Clone();
 
 float eff_jet[100],rate_jet[100],eff_qjet[100],rate_qjet[100],eff_sv[100],rate_sv[100];

 for (int ibin=0; ibin<hcrate_sv->GetNbinsX(); ibin++){
   hcrate_sv->SetBinContent(ibin,hrate_sv->Integral(ibin,101));
   hcrate_jet->SetBinContent(ibin,hrate_jet->Integral(ibin,101));
   hcrate_qjet->SetBinContent(ibin,hrate_qjet->Integral(ibin,101));
   rate_jet[ibin]=hrate_jet->Integral(ibin,101)/hrate_jet->Integral(0,101);
   rate_qjet[ibin]=hrate_qjet->Integral(ibin,101)/hrate_jet->Integral(0,101);
   rate_sv[ibin]=hrate_sv->Integral(ibin,101)/hrate_jet->Integral(0,101);
   eff_jet[ibin]=hpt_jet->Integral(ibin,101)/heff_den->Integral(0,101);
   eff_qjet[ibin]=hpt_qjet->Integral(ibin,101)/heff_den->Integral(0,101);
   eff_sv[ibin]=hpt_sv->Integral(ibin,101)/heff_den->Integral(0,101);
 } 

 heff_sv->Divide(heff_den);
 heff_jet->Divide(heff_den);
 heff_qjet->Divide(heff_den);

 TGraph * gr_jet = new TGraph(100,eff_jet,rate_jet);
 TGraph * gr_qjet = new TGraph(100,eff_qjet,rate_qjet);
 TGraph * gr_sv = new TGraph(100,eff_sv,rate_sv);

 float avgpt_x[30]={0};

 for (int ik=0; ik<30; ik++){
     avgpt_x[ik]=(ik+0.5);
     if (avgpt_nqjet[ik]>0)
        avgpt_qjet[ik]/=avgpt_nqjet[ik];
     if (avgpt_njet[ik]>0)
        avgpt_jet[ik]/=avgpt_njet[ik];
     if (avgpt_nsv[ik]>0)
        avgpt_sv[ik]/=avgpt_nsv[ik];
  }

 TGraph * gr_avgpt_jet = new TGraph(30,avgpt_x,avgpt_jet);
 TGraph * gr_avgpt_qjet = new TGraph(30,avgpt_x,avgpt_qjet);
 TGraph * gr_avgpt_sv = new TGraph(30,avgpt_x,avgpt_sv);


 gStyle->SetOptStat(0);

 TCanvas * cdr_bb = new TCanvas("cdr_bb","",800,600);
 hdr_bb->Draw();
 hdr_bb->GetXaxis()->SetTitle("DeltaR(b1,b2)");
 cdr_bb->SaveAs("eff_hdr_bb.png");

 TLegend * leg = new TLegend(0.7,0.3,1,0.6);
 leg->AddEntry(hdr_jet,"PF jet");
 leg->AddEntry(hdr_qjet,"PFC jet");
 leg->AddEntry(hdr_sv,"SV");

 TCanvas * cdr_reco = new TCanvas("cdr_reco","",800,600);
 hdr_jet->Draw();
 hdr_jet->SetLineColor(1);
 hdr_qjet->SetLineColor(2);
 hdr_sv->SetLineColor(3);
 hdr_qjet->Draw("sames"); 
 hdr_sv->Draw("sames");
 leg->Draw("sames");
 hdr_jet->GetXaxis()->SetTitle("DR(B,reco)");
 cdr_reco->SaveAs("eff_hdr_reco.png");

 TCanvas * ceff = new TCanvas("ceff","",800,600);
 heff_sv->Draw();
 heff_sv->SetLineColor(3);
 heff_jet->SetLineColor(1);
 heff_qjet->SetLineColor(2);
 heff_jet->Draw("sames");
 heff_qjet->Draw("sames");
 leg->Draw("sames");
 heff_sv->GetXaxis()->SetTitle("B meson pT");
 ceff->SaveAs("eff_vs.png");

 TCanvas * crate = new TCanvas("crate","",800,600);
 hcrate_jet->Draw();
 hcrate_jet->SetLineColor(1);
 hcrate_qjet->SetLineColor(2);
 hcrate_sv->SetLineColor(3);
 hcrate_qjet->Draw("sames");
 hcrate_sv->Draw("sames");
 leg->Draw("sames");
 crate ->SetLogy();
 hcrate_jet->GetYaxis()->SetTitle("Rate (# reco obj)");
 hcrate_jet->GetXaxis()->SetTitle("pT(Reco obj)");
 crate->SaveAs("eff_rate.png");

 TCanvas * cgr = new TCanvas("cgr","",800,600);
 gr_jet->SetLineColor(1);
 gr_jet->SetMarkerColor(1);
 gr_qjet->SetLineColor(2);
 gr_qjet->SetMarkerColor(2);
 gr_sv->SetLineColor(3);
 gr_sv->SetMarkerColor(3);
 cgr->SetLogy();
 gr_jet->Draw("AP*");
 gr_jet->GetYaxis()->SetRangeUser(0.001,1.2);
 gr_jet->GetYaxis()->SetTitle("Rate (# reco obj)");
 gr_jet->GetXaxis()->SetTitle("Efficiency");
 gr_qjet->Draw("P* sames");
 gr_sv->Draw("P* sames");
 leg->Draw("sames");
 cgr->SaveAs("eff_roc.png");


 TCanvas * cavg = new TCanvas("cavg","",800,600);
 gr_avgpt_jet->SetLineColor(1);
 gr_avgpt_jet->SetMarkerColor(1);
 gr_avgpt_qjet->SetLineColor(2);
 gr_avgpt_qjet->SetMarkerColor(2);
 gr_avgpt_sv->SetLineColor(3);
 gr_avgpt_sv->SetMarkerColor(3);
 gr_avgpt_jet->Draw("AP*");
 gr_avgpt_jet->GetYaxis()->SetTitle("<avg reco pT>");
 gr_avgpt_jet->GetXaxis()->SetTitle("gen B pT");
 gr_avgpt_qjet->Draw("P* sames");
 gr_avgpt_sv->Draw("P* sames");
 gr_avgpt_jet->GetYaxis()->SetRangeUser(0,80);
// leg->Draw("sames");
 cavg->SaveAs("eff_avg.png");


 return 0;
}
