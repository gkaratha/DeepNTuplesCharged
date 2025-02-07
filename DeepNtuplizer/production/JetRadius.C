#include "jet_tree.h"

int JetRadius(){

TChain * cc = new TChain("genbanalizer/tree");

cc->Add("output_0.root");

jet_tree tree;
tree.Init(cc);

TH1F * hdaughter_pdg = new TH1F("hdaughter_pdg","",5,1,6);
TH1F * hdaughters = new TH1F("hdaughters","",8,0,8);
TH1F * hdaughters_acc = new TH1F("hdaughters_acc","",8,0,8);
TH1F * hdr = new TH1F("hdr","",100,0,3);
TH1F * hdaughter_dir_lep = new TH1F("hdaughter_dir_lep","",2,0,2);
TH1F * hdaughter_dir_hdr = new TH1F("hdaughter_dir_hdr","",2,0,2);
TH2F * hpt_dr = new TH2F("hpt_dr","",40,0,40,20,0,3);
TH2F * hvispt_dr = new TH2F("hvispt_dr","",40,0,40,20,0,3);

float ntrk1p2[4]={0},ntrk1[4]={0},ntrk0p8[4]={0},ntrk0p6[4]={0},ntrk0p4[4]={0};
float avg[4]={0};



for (int ievt=0; ievt<cc->GetEntries(); ievt++){
   tree.GetEntry(ievt);
   if (ievt%1000==0) std::cout<<"ievt "<<ievt<<std::endl;
   std::vector<float> Bmeson_id, Bmeson_pt, Bmeson_eta, Bmeson_phi;
   std::vector<std::vector<float>> Bdaughters_id, Bdaughters_pt, Bdaughters_eta, Bdaughters_phi, Bdaughters_mom;
   std::vector<float> part_to_B;

   for (int ipart=0; ipart<tree.genpart_pdgId->size(); ipart++){
     if (Bmeson_id.size()==0){
        std::vector<float> part_pt{tree.genpart_pt->at(ipart)},
                               part_eta{tree.genpart_eta->at(ipart)},
                               part_phi{tree.genpart_phi->at(ipart)},
                               part_pdg{tree.genpart_pdgId->at(ipart)},
                               part_mom{tree.genpart_mother->at(ipart)};
         Bmeson_id.push_back(tree.genpart_Bmeson->at(ipart));
         Bmeson_pt.push_back(tree.genpart_Bmeson_pt->at(ipart));
         Bmeson_eta.push_back(tree.genpart_Bmeson_eta->at(ipart));
         Bmeson_phi.push_back(tree.genpart_Bmeson_phi->at(ipart));
         part_to_B.push_back(Bmeson_id.size()-1);
         Bdaughters_id.push_back(part_pdg);
         Bdaughters_pt.push_back(part_pt);
         Bdaughters_eta.push_back(part_eta);
         Bdaughters_phi.push_back(part_phi);
         Bdaughters_mom.push_back(part_mom);
         continue;
     }	     
     int foundBidx=-1;
     for (int ifound=0; ifound<Bmeson_id.size(); ifound++){
        if (Bmeson_id[ifound] != tree.genpart_Bmeson->at(ipart) ||
            Bmeson_pt[ifound] != tree.genpart_Bmeson_pt->at(ipart) ||
            Bmeson_eta[ifound] != tree.genpart_Bmeson_eta->at(ipart) ||
	    Bmeson_phi[ifound] != tree.genpart_Bmeson_phi->at(ipart) ) 
	      continue;
        foundBidx=ifound;
        break;
     }
     if (foundBidx==-1){
        std::vector<float> part_pt{tree.genpart_pt->at(ipart)},
		               part_eta{tree.genpart_eta->at(ipart)},
			       part_phi{tree.genpart_phi->at(ipart)},
			       part_pdg{tree.genpart_pdgId->at(ipart)},
			       part_mom{tree.genpart_mother->at(ipart)};	
        Bmeson_id.push_back(tree.genpart_Bmeson->at(ipart));
        Bmeson_pt.push_back(tree.genpart_Bmeson_pt->at(ipart));
        Bmeson_eta.push_back(tree.genpart_Bmeson_eta->at(ipart));
        Bmeson_phi.push_back(tree.genpart_Bmeson_phi->at(ipart));
        part_to_B.push_back(Bmeson_id.size()-1);
        Bdaughters_id.push_back(part_pdg);
        Bdaughters_pt.push_back(part_pt);
        Bdaughters_eta.push_back(part_eta);
        Bdaughters_phi.push_back(part_phi);
        Bdaughters_mom.push_back(part_mom);
     } else{
        part_to_B.push_back(foundBidx);
        Bdaughters_id[foundBidx].push_back(tree.genpart_pdgId->at(ipart));
        Bdaughters_pt[foundBidx].push_back(tree.genpart_pt->at(ipart));
        Bdaughters_eta[foundBidx].push_back(tree.genpart_eta->at(ipart));
        Bdaughters_phi[foundBidx].push_back(tree.genpart_phi->at(ipart));
        Bdaughters_mom[foundBidx].push_back(tree.genpart_mother->at(ipart));
     }
      

     // std::cout<<"ipart "<<ipart<<" pt "<<tree.genpart_Bmeson_pt->at(ipart)<<" eta "<<tree.genpart_Bmeson_eta->at(ipart) <<" phi "<<tree.genpart_Bmeson_phi->at(ipart)<<" pdg "<<tree.genpart_Bmeson->at(ipart)<<std::endl;	   
//      if (Bmeson_id.size()>0 && std::find(Bmeson_pt.begin(),Bmeson_pt.end(),tree.genpart_Bmeson_pt[ipart]) != Bmeson_pt.end() && std::find(Bmeson_eta.begin(),Bmeson_eta.end(),tree.genpart_Bmeson_eta[ipart]) != Bmeson_eta.end() && std::find(Bmeson_pt.begin(),Bmeson_phi.end(),tree.genpart_Bmeson_phi[ipart]) != Bmeson_phi.end()
   }

   //std::cout<<Bdaughters_id.size()<<endl;
   for (int ib=0; ib<Bdaughters_id.size(); ib++){
       hdaughters->Fill(Bdaughters_id[ib].size());
     //  std::cout<<Bdaughters_id[ib].size()<<std::endl;
       TLorentzVector vb;
       vb.SetPtEtaPhiM(Bmeson_pt[ib],Bmeson_eta[ib],Bmeson_phi[ib],5.28);
       TLorentzVector sum_vec;
       float ntrk_dr1p2=0,ntrk_dr1=0,ntrk_dr0p8=0,ntrk_dr0p6=0,ntrk_dr0p4=0;
       int nacc=0;
       float maxDR=0;
       for (int id=0; id<Bdaughters_id[ib].size(); id++){
          if (Bdaughters_pt[ib][id]< 0.5 || fabs(Bdaughters_eta[ib][id])>2.4)
	     continue;
          nacc++;	  
	  TLorentzVector vd;
          vd.SetPtEtaPhiM(Bdaughters_pt[ib][id],Bdaughters_eta[ib][id],Bdaughters_phi[ib][id],0.139);
	  sum_vec+=vd;
	  if (maxDR< vb.DeltaR(vd))
	      maxDR= vb.DeltaR(vd);
	  if (vb.DeltaR(vd)<1.2) ntrk_dr1p2+=1;
	  if (vb.DeltaR(vd)<1.0) ntrk_dr1+=1;
          if (vb.DeltaR(vd)<0.8) ntrk_dr0p8+=1;
	  if (vb.DeltaR(vd)<0.6) ntrk_dr0p6+=1;
	  if (vb.DeltaR(vd)<0.4) ntrk_dr0p4+=1;
          if(fabs(Bdaughters_id[ib][id])==11) {		  
	     hdaughter_pdg->Fill(1);
	     if ( fabs(Bdaughters_mom[ib][id])>500 )
	        hdaughter_dir_lep->Fill(0);
	     else
		hdaughter_dir_lep->Fill(1);     
	  } else if(fabs(Bdaughters_id[ib][id])==13){
             hdaughter_pdg->Fill(2);
	     if ( fabs(Bdaughters_mom[ib][id])>500 )
                hdaughter_dir_lep->Fill(0);
             else
                hdaughter_dir_lep->Fill(1);
	  } else{
             hdaughter_pdg->Fill(3);
	     if ( fabs(Bdaughters_mom[ib][id])>500 )
                hdaughter_dir_hdr->Fill(0);
             else
                hdaughter_dir_hdr->Fill(1);
	  }
           	 
       }
       if (nacc>1) {
          hdr->Fill(maxDR);
	  hpt_dr->Fill(Bmeson_pt[ib],maxDR);
	  hvispt_dr->Fill(sum_vec.Pt(),maxDR);
	  if (sum_vec.Pt()<5) {
             ntrk1p2[0]+=ntrk_dr1p2/nacc,ntrk1[0]+=ntrk_dr1/nacc,ntrk0p8[0]+=ntrk_dr0p8/nacc,ntrk0p6[0]+=ntrk_dr0p6/nacc,ntrk0p4[0]+=ntrk_dr0p4/nacc;
	     avg[0]+=1;
	  }else if (sum_vec.Pt()<10){
             ntrk1p2[1]+=ntrk_dr1p2/nacc,ntrk1[1]+=ntrk_dr1/nacc,ntrk0p8[1]+=ntrk_dr0p8/nacc,ntrk0p6[1]+=ntrk_dr0p6/nacc,ntrk0p4[1]+=ntrk_dr0p4/nacc;
	     avg[1]+=1;
	  }else if (sum_vec.Pt()<20){
             ntrk1p2[2]+=ntrk_dr1p2/nacc,ntrk1[2]+=ntrk_dr1/nacc,ntrk0p8[2]+=ntrk_dr0p8/nacc,ntrk0p6[2]+=ntrk_dr0p6/nacc,ntrk0p4[2]+=ntrk_dr0p4/nacc;
	     avg[2]+=1;
	  }else{
             ntrk1p2[3]+=ntrk_dr1p2/nacc,ntrk1[3]+=ntrk_dr1/nacc,ntrk0p8[3]+=ntrk_dr0p8/nacc,ntrk0p6[3]+=ntrk_dr0p6/nacc,ntrk0p4[3]+=ntrk_dr0p4/nacc;
	     avg[3]+=1;
	  }
       }  
       hdaughters_acc->Fill(nacc);
   }
 }
 float xaxis[]={2.5,7.5,15,30};
 for(int i=0; i<4; i++){
    ntrk1p2[i]/=avg[i],ntrk1[i]/=avg[i],ntrk0p8[i]/=avg[i],ntrk0p6[i]/=avg[i],ntrk0p4[i]/=avg[i];
 }
 TGraph * gr_ntrk1p2 = new TGraph(4,xaxis,ntrk1p2);
 TGraph * gr_ntrk1 = new TGraph(4,xaxis,ntrk1);
 TGraph * gr_ntrk0p8 = new TGraph(4,xaxis,ntrk0p8);
 TGraph * gr_ntrk0p6 = new TGraph(4,xaxis,ntrk0p6);
 TGraph * gr_ntrk0p4 = new TGraph(4,xaxis,ntrk0p4);

 gStyle->SetOptStat(0);
 TLegend * leg3 = new TLegend(0.7,0.2,1,0.4);
 leg3->AddEntry(gr_ntrk1p2,"DR<1.2");
 leg3->AddEntry(gr_ntrk1,"DR<1.0");
 leg3->AddEntry(gr_ntrk0p8,"DR<0.8");
 leg3->AddEntry(gr_ntrk0p6,"DR<0.6");
 leg3->AddEntry(gr_ntrk0p4,"DR<0.4");
 TCanvas * cntrk_dr = new TCanvas("cntrk_dr","",800,600);
 gr_ntrk1p2->SetLineColor(1);
 gr_ntrk1->SetLineColor(2);
 gr_ntrk0p8->SetLineColor(3);
 gr_ntrk0p6->SetLineColor(4);
 gr_ntrk0p4->SetLineColor(5);
 gr_ntrk1p2->Draw("AL");
 gr_ntrk1->Draw("L sames");
 gr_ntrk0p8->Draw("L sames");
 gr_ntrk0p6->Draw("L sames");
 gr_ntrk0p4->Draw("L sames");
 gr_ntrk1p2->GetXaxis()->SetTitle("Sum pT(charged)");
 gr_ntrk1p2->GetYaxis()->SetTitle("Fraction inside cone");
 gr_ntrk1p2->GetYaxis()->SetRangeUser(0,1.1);
 leg3->Draw("sames");
 cntrk_dr->SaveAs("cntrk_dr.png");
 

 TCanvas * cdpdg = new TCanvas("dpdg","",800,600);
 hdaughter_pdg->Draw("HIST");
 hdaughter_pdg->GetXaxis()->SetTitle("PDG");
 cdpdg->SaveAs("cdaughter_pdg.png");

 TLegend * leg = new TLegend(0.7,0.7,1,1);
 leg->AddEntry(hdaughters,"all");
 leg->AddEntry(hdaughters_acc,"pT>0.5 & #eta<2.4");
 TCanvas * cdaughter = new TCanvas("cdaughter","",800,600);
 hdaughters->Draw("HIST");
 hdaughters_acc->SetLineColor(2);
 hdaughters_acc->Draw("HIST sames");
 hdaughters->GetXaxis()->SetTitle("# charged daughter");
 leg->Draw("sames");
 cdaughter->SaveAs("cdaughter.png");

 TCanvas * cdaughter_dr = new TCanvas("cdaughter_dr","",800,600);
 hdr->Draw("HIST");
 hdr->GetXaxis()->SetTitle("max DR");
 cdaughter_dr->SaveAs("cdaughter_maxdr.png");

 TLegend * leg2 = new TLegend(0.7,0.7,1,1);
 leg2->AddEntry(hdaughter_dir_lep,"leptons");
 leg2->AddEntry(hdaughter_dir_hdr,"hadrons");
 TCanvas * cdaughter_dir = new TCanvas("cdaughter_dir","",800,600);
 hdaughter_dir_lep->Scale(1./hdaughter_dir_lep->Integral());
 hdaughter_dir_hdr->Scale(1./hdaughter_dir_hdr->Integral());
 hdaughter_dir_hdr->SetLineColor(2);
 hdaughter_dir_lep->GetYaxis()->SetRangeUser(0,1);
 hdaughter_dir_lep->Draw("HIST");
 hdaughter_dir_hdr->Draw("HIST sames");
 hdaughter_dir_lep->GetXaxis()->SetTitle("Indirect decay");
 leg2->Draw("sames");
 cdaughter_dir->SaveAs("cdaughter_dir.png");

 TCanvas * cpt_dr = new TCanvas("cpt_dr","",800,600);
 hpt_dr->Draw("COLZ");
 hpt_dr->GetYaxis()->SetTitle("DR");
 hpt_dr->GetXaxis()->SetTitle("pT(B)");
 cpt_dr->SaveAs("cpt_dr.png");

 TCanvas * cvispt_dr = new TCanvas("cvispt_dr","",800,600);
 hvispt_dr->Draw("COLZ");
 hvispt_dr->GetYaxis()->SetTitle("DR");
 hvispt_dr->GetXaxis()->SetTitle("Sum pT(charged)");
 cvispt_dr->SaveAs("cvispt_dr.png");
 
 return 0;
}
