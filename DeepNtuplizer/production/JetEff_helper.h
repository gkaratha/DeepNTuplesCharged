#include "jet_tree.h"

double deltaR(double eta1, double phi1, double eta2, double phi2) {
    double dEta = eta1 - eta2;
    double dPhi = TMath::Abs(phi1 - phi2);
    if (dPhi > TMath::Pi()) dPhi = 2 * TMath::Pi() - dPhi;
    return TMath::Sqrt(dEta * dEta + dPhi * dPhi);
}




std::vector< std::vector< std::vector<int> >> MatchJetConstToBdaughters( std::vector<std::vector<int>> Bconst,  std::vector<std::vector<int>>  B_daughter_matchedTrkIdx, std::vector<std::vector<int>>  B_daughter_matchedDaughterIdx){
  std::vector<std::vector<int>> Bconst_matched_genB;
  std::vector<std::vector<int>> Bconst_matched_trackIdx;
  std::vector<std::vector<int>> Bconst_matched_daughterIdx;
  for(auto atBconst: Bconst){
     std::vector<int> tmp_Bconst_matched_genB;
     std::vector<int> tmp_Bconst_matched_trackIdx;
     std::vector<int> tmp_Bconst_matched_daughterIdx;
     for(int itrk=0; itrk<atBconst.size(); itrk++){
        for(int igenB=0; igenB<B_daughter_matchedTrkIdx.size(); igenB++){
           if (std::find(B_daughter_matchedTrkIdx[igenB].begin(),B_daughter_matchedTrkIdx[igenB].end(),atBconst[itrk]) == B_daughter_matchedTrkIdx[igenB].end()  )
              continue;
           tmp_Bconst_matched_genB.push_back(igenB);
           tmp_Bconst_matched_trackIdx.push_back(itrk);
           tmp_Bconst_matched_daughterIdx.push_back( B_daughter_matchedDaughterIdx[igenB][std::find(B_daughter_matchedTrkIdx[igenB].begin(),B_daughter_matchedTrkIdx[igenB].end(),atBconst[itrk]) - B_daughter_matchedTrkIdx[igenB].begin()] );
           break;
        }
     }
     Bconst_matched_genB.push_back(tmp_Bconst_matched_genB);
     Bconst_matched_trackIdx.push_back(tmp_Bconst_matched_trackIdx);
     Bconst_matched_daughterIdx.push_back(tmp_Bconst_matched_daughterIdx);
  }

  return {Bconst_matched_genB,Bconst_matched_trackIdx,Bconst_matched_daughterIdx};
}


std::vector<std::vector< std::vector<int> >> MatchBdaughtersToTracks(std::vector< std::vector<float>> B_daughters_eta, std::vector< std::vector<float>> B_daughters_phi, std::vector<float>* qjetpart_eta, std::vector<float> * qjetpart_phi ){

  std::vector<std::vector<int>> B_daughter_matchedTrkIdx;
  std::vector<std::vector<int>> B_daughter_matchedDaughterIdx; 
  for (int ib=0; ib<B_daughters_eta.size(); ib++){
      
      std::vector<int> tmp_B_daughter_matchedTrkIdx;
      std::vector<int> tmp_B_daughter_matchedDaughterIdx;
      for (int id=0; id<B_daughters_eta[ib].size(); id++){
          float minDR=1000;
          int index=-1;
          for (int itrk=0; itrk<qjetpart_eta->size(); itrk++){
              if (minDR < deltaR(B_daughters_eta[ib][id],B_daughters_phi[ib][id],qjetpart_eta->at(itrk),qjetpart_phi->at(itrk)) )
                  continue;
              minDR = deltaR(B_daughters_eta[ib][id],B_daughters_phi[ib][id],qjetpart_eta->at(itrk),qjetpart_phi->at(itrk));
              index = itrk;
          }
          if (minDR>0.03)
             continue;
          tmp_B_daughter_matchedTrkIdx.push_back(index);
          tmp_B_daughter_matchedDaughterIdx.push_back(id);    
          
      }
      B_daughter_matchedTrkIdx.push_back(tmp_B_daughter_matchedTrkIdx);
      B_daughter_matchedDaughterIdx.push_back(tmp_B_daughter_matchedDaughterIdx);
  }
  
  return {B_daughter_matchedTrkIdx, B_daughter_matchedDaughterIdx};

}

 


std::vector<std::vector<int>> OrderConstituents(std::vector<int> * qjetpart_qjetIdx){
    std::vector<std::vector<int>> Bconst;
    std::vector<int> temp_Bconst;
    int last_jet=-1;
    for (int itrk=0; itrk<qjetpart_qjetIdx->size(); itrk++){
      if (temp_Bconst.size()>0 && last_jet != qjetpart_qjetIdx->at(itrk) ){
           Bconst.push_back(temp_Bconst);
           temp_Bconst.clear();
      }
      temp_Bconst.push_back(itrk);
      last_jet=qjetpart_qjetIdx->at(itrk);
    }
    return Bconst;
}



// put all daughters in order of B mesons
std::vector<std::vector<std::vector<float>> > CreateBdecayChain(std::vector<float> part_to_B, int nBmeson, std::vector<float>* genpart_pt,std::vector<float>* genpart_eta,std::vector<float>* genpart_phi, std::vector<float>* genpart_pdgId ) {


      std::vector< std::vector<float> > B_daughters_pt, B_daughters_eta,
                                        B_daughters_phi, B_daughters_pdgId;

      std::vector<float> already_stored;
      std::vector<float> B_pt;

      for ( int ig=0; ig<nBmeson; ig++){
          if ( already_stored.size()>0 && std::find(already_stored.begin(), already_stored.end(), part_to_B[ig])!= already_stored.end() )
              continue;
          already_stored.push_back( part_to_B[ig] );
          std::vector<float> tmp_daughters_pt,tmp_daughters_eta,tmp_daughters_phi, tmp_daughters_pdgId;
          tmp_daughters_pt.push_back(genpart_pt->at(ig));
          tmp_daughters_eta.push_back(genpart_eta->at(ig));
          tmp_daughters_phi.push_back(genpart_phi->at(ig));
          tmp_daughters_pdgId.push_back(genpart_pdgId->at(ig));

          for (int ig2=ig+1; ig2<nBmeson; ig2++){
              if (part_to_B[ig] != part_to_B[ig2])
                 continue;
              tmp_daughters_pt.push_back(genpart_pt->at(ig2));
              tmp_daughters_eta.push_back(genpart_eta->at(ig2));
              tmp_daughters_phi.push_back(genpart_phi->at(ig2));
              tmp_daughters_pdgId.push_back(genpart_pdgId->at(ig2));
          }

          B_daughters_pt.push_back(tmp_daughters_pt);
          B_daughters_eta.push_back(tmp_daughters_eta);
          B_daughters_phi.push_back(tmp_daughters_phi);
          B_daughters_pdgId.push_back(tmp_daughters_pdgId);
      }
  return {B_daughters_pt, B_daughters_eta, B_daughters_phi, B_daughters_pdgId };
}

//maps genparts to mother B
std::vector<float> connectPartToB(std::vector<float>* genpart_Bmeson, std::vector<float>* genpart_Bmeson_pt,std::vector<float>* genpart_Bmeson_eta,std::vector<float>* genpart_Bmeson_phi ) {

   std::vector<float> Bmeson_id, Bmeson_pt, Bmeson_eta, Bmeson_phi;
   std::vector<float> part_to_B;

   for (int ipart=0; ipart<genpart_Bmeson->size(); ipart++){
     if (Bmeson_id.size()==0){
         Bmeson_id.push_back(genpart_Bmeson->at(ipart));
         Bmeson_pt.push_back(genpart_Bmeson_pt->at(ipart));
         Bmeson_eta.push_back(genpart_Bmeson_eta->at(ipart));
         Bmeson_phi.push_back(genpart_Bmeson_phi->at(ipart));
         part_to_B.push_back(Bmeson_id.size()-1);
         continue;
     }
     int foundBidx=-1;
     for (int ifound=0; ifound<Bmeson_id.size(); ifound++){
        if (Bmeson_id[ifound] != genpart_Bmeson->at(ipart) ||
            Bmeson_pt[ifound] != genpart_Bmeson_pt->at(ipart) ||
            Bmeson_eta[ifound] != genpart_Bmeson_eta->at(ipart) ||
            Bmeson_phi[ifound] != genpart_Bmeson_phi->at(ipart) )
              continue;
        foundBidx=ifound;
        break;
     }
     if (foundBidx==-1){
        Bmeson_id.push_back(genpart_Bmeson->at(ipart));
        Bmeson_pt.push_back(genpart_Bmeson_pt->at(ipart));
        Bmeson_eta.push_back(genpart_Bmeson_eta->at(ipart));
        Bmeson_phi.push_back(genpart_Bmeson_phi->at(ipart));
        part_to_B.push_back(Bmeson_id.size()-1);
     } else{
        part_to_B.push_back(foundBidx);
     }
   }  
   return part_to_B;
}





std::pair<int,float> GenRecoMatchIdxDr(int nreco,std::vector<float>* reco_eta, std::vector<float>* reco_phi, float gen_eta,float gen_phi){
   float minDR=1000;
   int index=-1;
   TLorentzVector gvec;
   gvec.SetPtEtaPhiM(0,gen_eta,gen_phi,0);
   for (int ir=0; ir<nreco; ir++){
       TLorentzVector rvec;
       rvec.SetPtEtaPhiM(0,reco_eta->at(ir),reco_phi->at(ir),0);
       if (minDR<rvec.DeltaR(gvec) )
          continue;
       minDR=rvec.DeltaR(gvec);
       index=ir;
   }
   return std::make_pair(index,minDR);
}


void PlotHisto(TH1F* histo, TString name, TString xaxis, bool LogY){
  TCanvas * c1 = new TCanvas("c"+name,"",800,600);
  histo->Draw();
  histo->GetYaxis()->SetTitle("Events");
  histo->GetXaxis()->SetTitle(xaxis);
  histo->SetLineWidth(3);
  if (LogY) 
     c1->SetLogy();
  c1->SaveAs(name+".png");
}

void PlotHistos(vector<TH1F*> histos, TString name, TString xaxis, vector<TString> legs, bool LogY, float minY=-1, float maxY=-1){
  TCanvas * c1 = new TCanvas("c"+name,"",800,600);
  TLegend * leg = new TLegend(0.7,0.7,1,1);
  for(int idx=0; idx<histos.size(); idx++){
      histos[idx]->SetLineWidth(2);
      histos[idx]->SetLineColor(idx+1);
      leg ->AddEntry(histos[idx],legs[idx]);
      if (idx==0){
         histos[idx]->Draw();
         histos[idx]->GetXaxis()->SetTitle(xaxis);
      } else{
        histos[idx]->Draw("sames");
      }
        
  }
  leg->Draw("sames");
  if (LogY)
     c1->SetLogy();
  if (minY>0) histos[0]->SetMinimum(minY);
  if (maxY>0) histos[0]->SetMaximum(maxY);

  c1->SaveAs(name+".png");
}


TH1F*  PlotIntegrate(TH1F* histo){
  TH1F* hcopy = (TH1F*) histo->Clone();
  for (int ibin=1; ibin< histo->GetNbinsX()+1; ibin++){
    hcopy->SetBinContent(ibin, histo->Integral(ibin,histo->GetNbinsX()+1) );
  }
  return hcopy;
}


