{
//...
#include "HFlowEPres.h"
   //...

   Bool_t useEPcorr = kTRUE; // kFALSE;
   if (AnalysisParameter.Contains("EPcorr")) useEPcorr = kTRUE;

   HFlowEPres     *EPresCorr = NULL;
   TString         EPresCorrFile;
   vector<Float_t> EPres;
   Int_t           nHarmomics = 6;

   EPres.resize(nHarmomics);

   if (AnalysisParameter.Contains("apr12")) {
      // paraFileCentrality  = "/lustre/hades/user/bkardan/param/centrality_epcorr_apr12_gen8_2019_02_pass30.root";
      // occupancyWeightFile = "/lustre/hades/user/bkardan/weights/weights_PID0_day108_cent5.root";
      // occupancyNewFile    =
      // "/lustre/hades/user/bkardan/weights/flow2024_08_apr12_gen8_pass7_day107-109_PID0_Occupancy.root";
      // sectorFilelist      = "/lustre/hades/dst/apr12/gen8/sector_selection/FileListHadron.list";
      
      EPresCorrFile = "EPHistQA_apr12_gen8_pass37_EPResolution_TOFRPC_FWCharge_ParameterFileCal.root";
      //EPresCorrFile = "EPHistQA_mar19_ag158ag_3200A_gen6_pass11_EPResolution_TOFRPC_FWCharge_ParameterFileCal.root";
      gHades->setBeamTimeID(kApr12);
      sBeamtime = "apr12";
      // eBeam = 1230;  // in MeV/A
      // dEdxCorr.setDefaultPar("apr12");
   }

   if (useEPcorr) {
      if (EPresCorrFile != "") {
         EPresCorr     = new HFlowEPres("HFlowEPres", "HFlowEPres");
         fileEPresCorr = TFile::Open(EPresCorrFile.Data(), "READ");
         if (fileEPresCorr->IsOpen()) {
            std::cout << "weight file open: " << EPresCorrFile << std::endl;
            EPresCorr->LoadHistograms(fileEPresCorr);
         } else {
            cout << " WARNING: the EPres - list was not found." << endl;
            return 0;
         }
      } else
         return 0;
   }

   //---- Event Schleife

   // EP
   Float_t ep  = evtChara->getEventPlane(eEPcorr);
   Float_t epA = evtChara->getEventPlane(eEPcorr, 1);
   Float_t epB = evtChara->getEventPlane(eEPcorr, 2);
   if (ep == -1 || epA == -1 || epB == -1) return kFALSE; // if no EP and SubEP reconstructed skip event

   // Centrality
   Int_t centralityClass = evtChara->getCentralityClass(eCentEst, eCentClass);
   if (centralityClass < 1 || centralityClass > nCentrality) return kFALSE;
   Float_t cent = evtChara->getCentralityPercentile(eCentEst);
   // FIXME
   if (cent < 0. || cent > 70.) return kFALSE;

   Int_t TOFRPC       = evtChara->getCentralityEstimator(HParticleEvtChara::kTOFRPC);
   Int_t FWSumChargeZ = evtChara->getCentralityEstimator(HParticleEvtChara::kFWSumChargeZ);

   //-----------------------------------------------------------------
   for (Int_t iHarm = 0; iHarm < EPres.size(); iHarm++) EPres[iHarm] = 1.;
   if (EPresCorr) {
      for (Int_t iHarm = 0; iHarm < EPres.size(); iHarm++) {
         Double_t temp = EPresCorr->GetEPresCorr(iHarm + 1, TOFRPC, FWSumChargeZ);
         if (temp > 0. && temp < 0.99) EPres[iHarm] = temp;
         // cout << "EP : " << iHarm+1 << "  res:" << EPres[iHarm] << "  TOFRPC:" << TOFRPC << "  FW:"<< FWSumChargeZ <<
         // endl;
      }
   }

   //-----------------------------------------------------------------
   // Teilchen Schleife:
   //.....

   for (int iHarm = 0; iHarm < (int)nHarmomics; ++iHarm) {
      Float_t Vn = TMath::Cos((iHarm + 1.) * deltaPhi);
      Float_t Sn = TMath::Cos((iHarm + 1.) * deltaPhi);
      if (EPres[iHarm] > 0 && EPres[iHarm] < 1.) {
         // if(Vn > EPres[iHarm]) cout << "Vn : " << Vn <<" Harm: "<< iHarm+1 << "  res:" << EPres[iHarm] << "  TOFRPC:"
         // << TOFRPC << "  FW:"<< FWSumChargeZ <<  endl;
         Vn /= EPres[iHarm];
         Sn /= EPres[iHarm];
      }
      // Float_t weightEP = weight;
      // if(EPres[iHarm] > 0 && EPres[iHarm] < 1.) weightEP /= EPres[iHarm];
      // cout << "weightEP : " << weightEP <<" Harm: "<< iHarm+1 << "  res:" << EPres[iHarm] << "  TOFRPC:" << TOFRPC <<
      // "  FW:"<< FWSumChargeZ <<  endl;
      if (fProfile_PtYCent_FlowCos[iHarm][zent]) fProfile_PtYCent_FlowCos[iHarm]->Fill(Y, pt, V2, Vn, weight);
      if (fProfile_PtYCent_FlowSin[iHarm][zent]) fProfile_PtYCent_FlowSin[iHarm]->Fill(Y, pt, V2, Sn, weight);
      // if(fProfile_PtPzCent_FlowCos[iHarm]) fProfile_PtPzCent_FlowCos[iHarm]->Fill(pz, pt,  cent, Vn, weight);
   }

   return 0;
}
