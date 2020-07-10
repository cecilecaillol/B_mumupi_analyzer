      float run, lumi, evt=-1;
      float genpt_1, geneta_1, genphi_1, genq_1=-1; // generated leading muon
      float genpt_2, geneta_2, genphi_2, genq_2=-1; // generated subleading muon
      float genpt_3, geneta_3, genphi_3, genq_3=-1; // generated pion 
      float pt_1, pterr_1, eta_1, phi_1, q_1, looseID_1, mediumID_1=-1; // leading muon
      float pt_2, pterr_2, eta_2, phi_2, q_2, looseID_2, mediumID_2=-1; // subleading muon
      float pt_3, eta_3, phi_3, q_3=-1; // pion 
      float vtxprob_12, vtxprob_23, vtxprob_13, vtxprob_123=-1; // vertex fit probability
      int HLT_Mu7_IP4, HLT_Mu8_IP3, HLT_Mu8_IP5, HLT_Mu9_IP5, HLT_Mu9_IP6, HLT_Mu12_IP6=0; // B parking triggers
      int match_Mu7_IP4_1, match_Mu8_IP3_1, match_Mu8_IP5_1, match_Mu9_IP5_1, match_Mu9_IP6_1, match_Mu12_IP6_1=0;
      int match_Mu7_IP4_2, match_Mu8_IP3_2, match_Mu8_IP5_2, match_Mu9_IP5_2, match_Mu9_IP6_2, match_Mu12_IP6_2=0;
      float validFraction_1, normalizedChi2_1, chi2LocalPosition_1, trkKink_1, segmentCompatibility_1, softMvaValue_1=-1;
      float validFraction_2, normalizedChi2_2, chi2LocalPosition_2, trkKink_2, segmentCompatibility_2, softMvaValue_2=-1;
      float m123=-1;

