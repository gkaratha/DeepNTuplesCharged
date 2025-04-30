import ROOT as rt

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(0)

cc = rt.TChain("deepntuplizer/tree")
cc.Add("/eos/cms/store/cmst3/group/softJets/gkaratha/SoftMultiJet/DeepNtuples_v3/CRAB_UserFiles/PFC_Signal_chain_m70_dm20_13_04_25/250413_155615/0000/*1.root")


#Blabel = "(isB || isBB || isLeptonicB || isGBB || isLeptonicB_C)"
#Clabel = "(isC || isCC || isGCC)"
#BKGlabel = "(isU || isD || isG || isS)"
#name="oldlabel"

Blabel = "(isMatchedB2d || isMatchedB3d || isMatchedBMore3d)"
Clabel = "(isMatchedC2d || isMatchedC3d || isMatchedCMore3d)"
BKGlabel = "(isNotMatched)"
name="newlabel"


cut="1"

obj_plots = {
            #### SV variables
            "sv_etarel":{"min":-0.7,"max":0,"name":"svdeta_jet"},\
            "sv_phirel":{"min":-0.7,"max":0,"name":"svdphi_jet"},\
            "sv_eta":{"min":-2.5,"max":2.5,"name":"sv_eta"},\
            "sv_phi":{"min":-3.2,"max":3.2,"name":"sv_phi"},\
            "sv_pt":{"min":0,"max":10,"name":"sv_pt"},\
            "sv_e":{"min":0,"max":10,"name":"sv_e"},\
            "sv_enratio":{"min":-100,"max":0,"name":"sv_enratio"},\
            "sv_costhetasvpv":{"min":-100,"max":0,"name":"sv_cos"},\
            "sv_d3d":{"min":0,"max":3.0,"name":"sv_d3d"},\
            "sv_d3dsig":{"min":0,"max":20.0,"name":"sv_d3dsig"},\
            "sv_dxy":{"min":0,"max":0.7,"name":"sv_dxy"},\
            "sv_dxysig":{"min":0,"max":20.0,"name":"sv_dxysig"},\
            "sv_normchi2":{"min":0,"max":5,"name":"sv_normchi2"},\
            "sv_chi2":{"min":0,"max":15,"name":"sv_chi2"},\
            "sv_mass":{"min":0,"max":10,"name":"sv_mass"},\
            "sv_ntracks":{"min":0,"max":10,"name":"sv_ntrk"},\
            "sv_deltaR":{"min":-0.6,"max":0,"name":"sv_dr"},\
            ### global vars
            "jet_pt":{"min":0,"max":50,"name":"glb_jetpt"},\
            "jet_eta":{"min":-2.5,"max":2.5,"name":"glb_jeteta"},\
            "npv":{"min":20,"max":100,"name":"glb_jetnpv"},\
            "n_Cpfcand":{"min":0,"max":50,"name":"glb_jet_ncpf"},\
            "nsv":{"min":0,"max":10,"name":"glb_jetnsv"},\
            "TagVarCSV_trackSumJetEtRatio":{"min":0,"max":2,"name":"glb_csvjet_etratio"},\
            "TagVarCSV_trackSumJetDeltaR":{"min":0,"max":0.3,"name":"glb_csvjet_dr"},\
            "TagVarCSV_vertexCategory":{"min":0,"max":5,"name":"glb_csvjet_vtx_cat"},\
            "TagVarCSV_trackSip2dValAboveCharm":{"min":-0.05,"max":0.05,"name":"glb_csvjet_ip2d"},\
             "TagVarCSV_trackSip2dSigAboveCharm":{"min":-10,"max":10,"name":"glb_csvjet_sip2d"},\
             "TagVarCSV_trackSip3dValAboveCharm":{"min":-0.05,"max":0.05,"name":"glb_csvjet_ip3d"},\
            "TagVarCSV_trackSip3dSigAboveCharm":{"min":-10,"max":10,"name":"glb_csvjet_sip3d"},\
            "TagVarCSV_jetNSelectedTracks":{"min":0,"max":10,"name":"glb_csvtrk_nseltracks"},\
            "TagVarCSV_jetNTracksEtaRel":{"min":0,"max":10,"name":"glb_csvjet_etarel"},\
            #### PF vars
            "Cpfcan_etarel":{"min":-0.5,"max":0.1,"name":"cpfdeta"},\
            "Cpfcan_phirel":{"min":-0.5,"max":0.1,"name":"cpfdphi"},\
            "Cpfcan_pt":{"min":0,"max":30,"name":"cpf_pt"},\
            "Cpfcan_e":{"min":0,"max":30,"name":"cpf_e"},\
            "Cpfcan_eta":{"min":-2.5,"max":2.5,"name":"cpf_eta"},\
            "Cpfcan_phi":{"min":-3.2,"max":3.2,"name":"cpf_phi"},\
            "Cpfcan_quality":{"min":0,"max":10,"name":"cpf_qual"},\
            "Cpfcan_VTX_ass":{"min":3,"max":8,"name":"cpf_vtxas"},\
            "Cpfcan_chi2":{"min":0,"max":100,"name":"cpf_chi2"},\
            "Cpfcan_puppiw":{"min":0,"max":1.5,"name":"cpf_puppi"},\
            "Cpfcan_drminsv":{"min":-0.4,"max":0,"name":"cpf_drminsv"},\
            "Cpfcan_ptrel":{"min":-2,"max":0,"name":"cpf_ptrel"},\
            "Cpfcan_BtagPf_trackEtaRel":{"min":0,"max":6,"name":"cpf_trketarel"},\
            "Cpfcan_BtagPf_trackPtRel":{"min":0,"max":4,"name":"cpf_trkptrel"},\
            "Cpfcan_BtagPf_trackPPar":{"min":0,"max":20,"name":"cpf_trkppar"},\
            "Cpfcan_BtagPf_trackDeltaR":{"min":0,"max":1,"name":"cpf_trkdr"},\
            "Cpfcan_BtagPf_trackPParRatio":{"min":0,"max":1,"name":"cpf_trkpparatio"},\
            "Cpfcan_BtagPf_trackSip2dVal":{"min":-0.5,"max":0.5,"name":"cpf_trksip2d"},\
            "Cpfcan_BtagPf_trackSip2dSig":{"min":0,"max":100,"name":"cpf_trksip2dsig"},\
            "Cpfcan_BtagPf_trackSip3dVal":{"min":-1,"max":5,"name":"cpf_trksip3d"},\
            "Cpfcan_BtagPf_trackSip3dSig":{"min":0,"max":100,"name":"cpf_trksip3dsig"},\
            "Cpfcan_BtagPf_trackJetDistVal":{"min":-10,"max":5,"name":"cpf_trkjetdist"},\

            ##### pair variables
            "pair_pca_distance":{"min":0,"max":20,"name":"pair_pca"},\
            "pair_pca_significance":{"min":0,"max":100,"name":"pair_pcasig"},\
            "pair_pcaSeed_x2":{"min":-100,"max":100,"name":"pair_pca_x"},\
            "pair_pcaSeed_y2":{"min":-100,"max":100,"name":"pair_pca_y"},\
            "pair_pcaSeed_z2":{"min":-400,"max":400,"name":"pair_pca_z"},\
            "pair_pcaSeed_xerr2":{"min":0,"max":0.15,"name":"pair_pca_xerr"},\
            "pair_pcaSeed_yerr2":{"min":0,"max":0.15,"name":"pair_pca_yerr"},\
            "pair_pcaSeed_zerr2":{"min":0,"max":0.09,"name":"pair_pca_zerr"},\
            "pair_dotprod1":{"min":-1,"max":1.5,"name":"pair_dotprod_trkksi"},\
            "pair_pca_dist2":{"min":0,"max":100,"name":"pair_pca_pv"},\
            "pair_dotprod12_2D":{"min":-1,"max":1.5,"name":"pair_dotprod_trk1trk2"},\
            "pair_dotprod12_2DV":{"min":-1,"max":1.5,"name":"pair_dotprod_ksi1ksi2"},\
            "pair_dotprod12_3D":{"min":-1,"max":1.5,"name":"pair_dotprod_trk1trk2_3d"},\
            "pair_dotprod12_3DV":{"min":-1,"max":1.5,"name":"pair_dotprod_ksi1ksi2_3d"},\
            "pair_pca_jetAxis_dist":{"min":0,"max":10,"name":"pair_dist_k_jet"},\
            "pair_pca_jetAxis_dotprod":{"min":0,"max":1.2,"name":"pair_dorprod_zi_jet"},\
            "pair_pca_jetAxis_dEta":{"min":0,"max":10,"name":"pair_deta_jetaxis"},\
            "pair_pca_jetAxis_dPhi":{"min":0,"max":5,"name":"pair_dphi_jetaxis"}
           }


for plot in obj_plots.keys():
    ctmp = rt.TCanvas("c"+plot,"",800,600)
    htmp_b = rt.TH1F("b"+obj_plots[plot]["name"],"",50,obj_plots[plot]["min"],obj_plots[plot]["max"])
    htmp_c = rt.TH1F("c"+obj_plots[plot]["name"],"",50,obj_plots[plot]["min"],obj_plots[plot]["max"])
    htmp_uds = rt.TH1F("uds"+obj_plots[plot]["name"],"",50,obj_plots[plot]["min"],obj_plots[plot]["max"])
    cc.Draw(plot+">>b"+obj_plots[plot]["name"],cut+" && "+Blabel)
    cc.Draw(plot+">>c"+obj_plots[plot]["name"],cut+" && "+Clabel)
    cc.Draw(plot+">>uds"+obj_plots[plot]["name"],cut+" && "+BKGlabel)
    leg = rt.TLegend(0.8,0.8,1,1)
    leg.AddEntry(htmp_b,"from B")
    leg.AddEntry(htmp_c,"from C")
    leg.AddEntry(htmp_uds,"from light")

    ihtmp=0
    for htmp in [htmp_b,htmp_c,htmp_uds]: 
      htmp.SetLineColor(1+ihtmp)
      htmp.SetLineWidth(3)
      if htmp.Integral()>0:
         htmp.Scale(1./htmp.Integral());
      else:
         print("empty histo "+plot)
      if ihtmp==0:
         htmp.GetXaxis().SetTitle(obj_plots[plot]["name"])
         htmp.Draw("HIST")
      else:
         htmp.Draw("HIST sames")
      ihtmp+=1
    leg.Draw("sames")
    ctmp.SaveAs("input_cmpr_"+name+"_"+obj_plots[plot]["name"]+".png")
    ctmp.SetLogy()
    ctmp.SaveAs("input_cmpr_"+name+"_"+obj_plots[plot]["name"]+"_log.png")

'''

    ########################### only B ####################################
    ctmpb = rt.TCanvas("cb"+plot,"",800,600)
    htmpb_b = rt.TH1F("b_b"+obj_plots[plot]["name"],"",50,obj_plots[plot]["min"],obj_plots[plot]["max"])
    htmpb_bb = rt.TH1F("b_bb"+obj_plots[plot]["name"],"",50,obj_plots[plot]["min"],obj_plots[plot]["max"])
    htmpb_gbb = rt.TH1F("b_gbb"+obj_plots[plot]["name"],"",50,obj_plots[plot]["min"],obj_plots[plot]["max"])
    htmpb_lepb = rt.TH1F("b_lepb"+obj_plots[plot]["name"],"",50,obj_plots[plot]["min"],obj_plots[plot]["max"])
    cc.Draw(plot+">>b_b"+obj_plots[plot]["name"],cut+" && isB")
    cc.Draw(plot+">>b_bb"+obj_plots[plot]["name"],cut+" && isBB")
    cc.Draw(plot+">>b_gbb"+obj_plots[plot]["name"],cut+" && isGBB")
    cc.Draw(plot+">>b_lepb"+obj_plots[plot]["name"],cut+" && isLeptonicB")

    legb = rt.TLegend(0.8,0.8,1,1)
    legb.AddEntry(htmpb_b,"cat B: B")
    legb.AddEntry(htmpb_bb,"cat B: BB")
    legb.AddEntry(htmpb_gbb,"cat B: gBB")
    legb.AddEntry(htmpb_lepb,"cat B: lep B")


    ihtmp=0
    for htmp in [htmpb_b,htmpb_bb,htmpb_gbb,htmpb_lepb]:
      htmp.SetLineColor(1+ihtmp)
      htmp.SetLineWidth(3)
      if htmp.Integral()>0:
         htmp.Scale(1./htmp.Integral());
      else:
         print("B plot: empty histo "+plot+" #"+str(ihtmp))
      if ihtmp==0:
         htmp.GetXaxis().SetTitle(obj_plots[plot]["name"])
         htmp.Draw("HIST")
      else:
         htmp.Draw("HIST sames")
      ihtmp+=1
    legb.Draw("sames")
    ctmpb.SaveAs("Binput_cmpr_"+obj_plots[plot]["name"]+".png")
    ctmpb.SetLogy()
    ctmpb.SaveAs("Binput_cmpr_"+obj_plots[plot]["name"]+"_log.png")

'''
