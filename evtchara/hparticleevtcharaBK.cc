//_HADES_CLASS_DESCRIPTION 
////////////////////////////////////////////////////////////////////////////
//*-- AUTHOR : B. Kardan  04.06.2019
//*-- VERSION : 0.72
//
// HParticleEvtCharaBK
//
// Purpose: EventCharacterization
// - Centrality from Hit(TOF,RPC,FW) and Track Estimators in Data/Sim
// - QVector and Phi (Reaction Plane Estimate) from FW
// - Event-weight for downscaled Events PT2
//
//--------------------------------------------------------------------------
// Usage:
//
//  - input files can be found at : example /cvmfs/hades.gsi.de/param/eventchara/
//
//  - to define the ParameterFile where Classes and Estimators are stored use
//    setParameterFile("/cvmfs/hades.gsi.de/param/eventchara/ParameterFile.root")
//    
//  - to print the definition of estimator & class use
//    printCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k10);
//  
//  - to get the CentralityClass of an event (with estimator and class definition) use 
//    getCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k10);
//
//  - to get the CentralityPercentile of an event (only estimator is needed) use 
//    getCentralityPercentile(HParticleEvtCharaBK::kTOFRPC);
//
//  - to get the EventWeight to re-weight downscaled events use
//    getEventWeight();
//    
//  - to get the EventPlane with ReCentering use
//    getEventPlane(HParticleEvtCharaBK::kReCentering);
//  
//--------------------------------------------------------------------------  
//   Estimators:
//   
//     TOFRPC                   - (default) TOF and RPC hit multiplicity in time-window
//     TOF                      - TOF hit multiplicity in time-window (0 < tof < 35ns)
//     RPC                      - RPC hit multiplicity in time-window (0 < tof < 25ns)
//     TOFRPCtot                - total TOF and RPC hit multiplicity in event-window
//     SelectedParticleCand     - selected Particle multiplicity
//     PrimaryParticleCand      - primary Particle multiplicity
//     SelectedParticleCandCorr - selected Particle multiplicity corrected by the 
//                                running mean and scaled with the reference mean 
//                                (selTrack * referenceMean/<selTrack>)
//     SelectedParticleCandSecCorr
//                               - selected Particle multiplicity corrected by the 
//                                running mean and scaled with the reference mean 
//                                (selTrack * referenceMean/<selTrack>)
//     SelectedParticleCandNorm - selected Particle multiplicity normalized by 
//                                the running mean (selTrack/<selTrack>)
//     FWSumChargeSpec          - sum of charge (dE/dx in a.u.) of Spectator in 
//                                FW acceptance with beta~0.9
//     FWSumChargeZ             - sum of charge (int Z till 6 with fixed cuts in dE/dx)
//                                of Spectator in FW acceptance with beta~0.9
//
//   Classes:
//     2                   - 2% classes
//     5                   - 5% classes
//     10                  - (default) 10% classes
//     13                  - 13% classes
//     15                  - 15% classes
//     20                  - 20% classes
//     FOPI                - FOPI centrality classes
//   
//   
//   EventPlaneCorrection:
//     kNoCorrection            - EP only selection of spectator candidates in FW,
//                                no Correction
//     kSmear                   - smearing of x and y of each FW hit inside cell size
//     kShiftFW                 - global shift x and y with the centre of gravity of
//                                all FW hits in all events (per class/day)
//     kScalingFW               - global scale x and yith the sigma of the distribution
//                                all FW hits in all events (per class/day)
//     kReCentering             - re-centering of QVector with <Qx> and <Qy> (calc. evt-by-evt)
//     kScaling                 - scaling of QVector with the sigma of Qx and Qy (calc. evt-by-evt)
//     kRotation                - rotation of EP via residual Fourier harmonics 
//                                after re-centering and scaling
//   
//--------------------------------------------------------------------------
// quick how-To:
//  - for an example see demo-macro:
//       /scripts/batch/GE/loopDST/loopDST.C
//
//  - include header-file
//      #include "HParticleEvtCharaBK.h"
//  
//  - before eventLoop:
//       HParticleEvtCharaBK evtChara;
//           ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_apr12_gen8_2018_06.root";
//  
//       if(isSimulation)
//           ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_sim_au1230au_gen8a_UrQMD_minbias_2018_06.root";
//    
//       if(!evtChara.setParameterFile(ParameterfileCVMFS)){
//           cout << "Parameterfile not found !!! " << endl;return kFALSE;
//       }
//  
//       if(!evtChara.init()) {cout << "HParticleEvtCharaBK not init!!! " << endl;return kFALSE;}
//  
//       evtChara.printCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k5);
//       evtChara.printCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k10);
//  
//  - inside event-loop:
//  
//    Int_t   CentralityClassTOFRPC = evtChara.getCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k10);
//    10% Centrality-Classes:       1(0-10%) - 5(40-50%) ... 0:Overflow max:Underflow
//  
//    Float_t CentralityTOFRPC      = evtChara.getCentralityPercentile(HParticleEvtCharaBK::kTOFRPC);
//    CentralityPercentile:         0 - 100% percentile of the total cross section
//                                      101% Over-,Underflow or Outlier Events
//  
//    Float_t EventPlane            = evtChara.getEventPlane(HParticleEvtCharaBK::kReCentering|HParticleEvtCharaBK::kScaling);
//    EventPlane:                   0 - 2pi (in rad)  re-centered & scaled EP
//                                 -1   no EP could be reconstructed
//
////////////////////////////////////////////////////////////////////////////

#include "hparticleevtcharaBK.h"
ClassImp(HParticleEvtCharaBK)

// const char* HParticleEvtCharaBK::StringEventPlaneParameter[HParticleEvtCharaBK::kNumEventPlaneParameter] = {
//     "FWx"  
//    ,"FWy"  
//    ,"Qx"  
//    ,"Qy"  
//    ,"FWxWCharge"  
//    ,"FWyWCharge"  
//    ,"QxWCharge"  
//    ,"QyWCharge"  
//    ,"FourierC1"  
//    ,"FourierC2"  
//    ,"FourierC3"  
//    ,"FourierC4"  
//    ,"FourierS1"  
//    ,"FourierS2"  
//    ,"FourierS3"  
//    ,"FourierS4"  
//    ,"Resolution"  
//    ,"ResolutionWCharge"  
//    ,"Chi" 
// };


HParticleEvtCharaBK::HParticleEvtCharaBK(const Text_t* name,const Text_t* title)
   : HReconstructor(name,title) {
        fParticleEvtInfoCat   = 0;
        fParameterFile        = "centrality.root";
        fPostFix              = "";
        isSimulation          = kFALSE;
        fReferenceMeanSelTrack= 36.44;
        fDayOfYear            = 0;
        fEventPlaneCorrectionFlag = kNoCorrection;
        fQVectorCalcDone = kFALSE;
        //maxFWCells = 302;
        
        useFWCut.resize(kNumFWCutValues);
        useFWCut[0]=kFALSE;  //kBetaCuts  
        useFWCut[1]=kTRUE;   //kTimeCuts  
        useFWCut[2]=kTRUE;   //kChargeCuts
        excludeNoisyFWcells = kTRUE;
        
        fRandom       = new TRandom2();
        
        vQPhi.resize(3);

        fCentralityPercentileHist.resize(kNumCentralityEstimator);
        fCentralityHist.resize(kNumCentralityEstimator);
        fEstimatorHist.resize(kNumCentralityEstimator);
        for (int centEst = 0; centEst < (int)kNumCentralityEstimator; ++centEst) fCentralityHist[centEst].resize(kNumCentralityClass);
        fEventPlaneCorrectionHist.resize(kNumEventPlaneParameter);
        fFWCutValuesHist.resize(kNumFWCutValues);
        for (int cutValue = 0; cutValue < (int)kNumFWCutValues; ++cutValue) fFWCutValuesHist[cutValue].resize(MAXFWCELLS);
        
        
        fFWminBeta.resize(3);
        fFWmaxBeta.resize(3);
        fFWminCharge.resize(3);
        fFWminBeta[0]=0.84;   fFWmaxBeta[0]=1.2;    fFWminCharge[0]=80; // small cells  
        fFWminBeta[1]=0.85;   fFWmaxBeta[1]=1.2;    fFWminCharge[1]=84; // medium cells
        fFWminBeta[2]=0.8;    fFWmaxBeta[2]=1.2;    fFWminCharge[2]=88; // large cells
        
        fFWChargeCuts.resize(6);
        fFWChargeCuts[5]=386;
        fFWChargeCuts[4]=339;
        fFWChargeCuts[3]=298;
        fFWChargeCuts[2]=241;
        fFWChargeCuts[1]=175;
        fFWChargeCuts[0]=84 ;
        
    }

HParticleEvtCharaBK::~HParticleEvtCharaBK() {
}

Bool_t HParticleEvtCharaBK::init()
{
    if(gHades){
        HRecEvent* evt =(HRecEvent*) gHades->getCurrentEvent();
        if(evt){

            HCategory* catKin=HCategoryManager::getCategory(catGeantKine,kFALSE,"catGeantKine");
            if(catKin) {isSimulation=kTRUE; Info("init()","GeantKine found - is Simulation");}
            else       {isSimulation=kFALSE;Info("init()","GeantKine not found - is not Simulation");}


            fParticleEvtInfoCat =  HCategoryManager::getCategory(catParticleEvtInfo,kTRUE,"catParticleEvtInfo, from HParticleEvtCharaBK::init()");
            if(!fParticleEvtInfoCat) { Info("init()","No catParticleEvtInfo in input!"); return kFALSE;}

            fParticleCandCat =  HCategoryManager::getCategory(catParticleCand,kTRUE,"catParticleCand, from HParticleEvtCharaBK::init()");
            if(!fParticleCandCat) { Info("init()","No catParticleCand in input!");}

            fCatWallHit = HCategoryManager::getCategory(catWallHit,kTRUE,"catWallHit, from HParticleEvtCharaBK::init()");
            if(!fCatWallHit) { Info("init()","No catWallHit in input!");}

        } else {
            Error("init()","NO EventStructure found!");
            return kFALSE;
    }

    } else {
        Error("init()","NO HADES found!");
        return kFALSE;
    }
    // read parameter file
    if(fParameterFile){
        if(!loadParameterFile()) return kFALSE;
    } else {
        Error("init()","NO Parameterfile found!");
        return kFALSE;
    }
    loadDayOfYear();
    
    return kTRUE;
}

Bool_t HParticleEvtCharaBK::reinit()
{
    loadDayOfYear();
    return kTRUE;
}

Int_t HParticleEvtCharaBK::loadDayOfYear()
{
    TString tempFileName = "";
    if(gLoop){
        gLoop->isNewFile(tempFileName);
    }
    else if(gHades){
        if(gHades->getRuntimeDb()){
            tempFileName = gHades->getRuntimeDb()->getCurrentFileName();
        }
        else if(gHades->getDataSource()){
            tempFileName = gHades->getDataSource()->getCurrentFileName();
        }
    }
    if(tempFileName==""){
        //Error("loadDayOfYear()","File/Day not found! using default day108");
        fDayOfYear = 108; // default FIXME  in apr12
        if(isSimulation) fDayOfYear = 0;
        return fDayOfYear;
    }
    fDayOfYear = HTime::getDayFileName(HTime::stripFileName(tempFileName,kTRUE,kTRUE));
    //Info("loadDayOfYear()",">>> Day of Year: %d",fDayOfYear);
    
    if( fDayOfYear>365 && isSimulation){
        fDayOfYear = 0;
    }
    
    return fDayOfYear;
}

Bool_t HParticleEvtCharaBK::setParameterFile(TString ParameterFile)
{
    fParameterFile = ParameterFile;
    TString path = gSystem->ExpandPathName(fParameterFile.Data());
    if (gSystem->AccessPathName(path)) {
        Error("init()","File  %s does not exist!",path.Data());
        return kFALSE;
    }
    return kTRUE;
}

Bool_t HParticleEvtCharaBK::loadParameterFile()
{
    // read parameter file
    TString path = gSystem->ExpandPathName(fParameterFile.Data());
    if (gSystem->AccessPathName(path)) {
        Error("loadParameterFile()","File  %s does not exist!",path.Data());
        return kFALSE;
    } else {
        fFile  = new TFile(path,"OPEN");
        TObject* ParameterFileVersion = (TObject*) fFile->Get("HParticleEvtCharaVersion");
        if(ParameterFileVersion){
            TString Version(ParameterFileVersion->GetTitle());
            Float_t fVersion = Version.Atof();
            if(fVersion < getVersion()){
                Error("loadParameterFile()","File  %s is out-dated! Needed Version: %02.1f -- in file: %02.1f",path.Data(), getVersion(), fVersion);
                return kFALSE;
            }
            cout<<"\n--------------------------------------------------------------------------------------"  << endl;
            Info("loadParameterFile()",">>> Parameter input file (ver. %02.1f) : %s",fVersion,path.Data());
            
            cout<<"\n--------------------------------------------------------------------------------------"  << endl;
        }
        else{
            Error("loadParameterFile()","In File  %s no Version information found!",path.Data());
            return kFALSE;
        }
    }
    loadCentralityEstimatorHist();
    loadEventPlaneCorrectionHist();
    loadFWCutValuesHist();
    return kTRUE;
}

Bool_t HParticleEvtCharaBK::saveParameterFile()
{
    // read parameter file
    TString path = gSystem->ExpandPathName(fParameterFile.Data());
    //if (gSystem->AccessPathName(path)) {
    //    Error("loadParameterFile()","File  %s does not exist!",path.Data());
    //    return kFALSE;
    //} else {
    fFile  = new TFile(path,"CREATE");
    if(!fFile) return kFALSE;
    cout << "Version of HParticleEvtCharaBK:  " << getVersion() << endl;
    TString version=Form("%02.1f",getVersion());
    TNamed OutputVersion("HParticleEvtCharaVersion",version.Data()); 
    OutputVersion.Write();
    saveCentralityEstimatorHist();
    saveEventPlaneCorrectionHist();
    saveFWCutValuesHist();
    return kTRUE;
}

Bool_t HParticleEvtCharaBK::loadCentralityEstimatorHist()
{

    cout<<"\n--------------------------------------------------------------------------------------"  << endl;
    Info("loadCentralityEstimatorHist()","Calibration for Centrality Estimator and Classes loading:");
    Int_t n=0,m=0;
    for (int centEst = 0; centEst < (int) kNumCentralityEstimator; ++centEst){
       printf("\n>>>> %31s %s:  ",getStringCentralityEstimator(centEst).Data() ,fPostFix.Data());
       TString temp;
       if(fPostFix.CompareTo("")==0) temp = Form("%s_percentile",getStringCentralityEstimator(centEst).Data() );
       else                       temp = Form("%s_%s_percentile",getStringCentralityEstimator(centEst).Data(), fPostFix.Data() );
       fCentralityPercentileHist[centEst] = (TH1F*) fFile->FindObjectAny( temp.Data() );
       if(fCentralityPercentileHist[centEst]){
            printf("percentile");
            m++;
        }
       
       for (int centC = 0; centC < (int) kNumCentralityClass; ++centC){
           TString temp2;
           if(fPostFix.CompareTo("")==0) temp2 = Form("%s_%s_fixedCuts",getStringCentralityEstimator(centEst).Data(), getStringCentralityClass(centC).Data() );
           else                       temp2 = Form("%s_%s_%s_fixedCuts",getStringCentralityEstimator(centEst).Data(), fPostFix.Data(), getStringCentralityClass(centC).Data() );
           fCentralityHist[centEst][centC] = (TH1F*) fFile->FindObjectAny( temp2.Data() );
           if(fCentralityHist[centEst][centC]){
                printf(" %5s",getStringCentralityClass(centC).Data() );
                n++;
           }
       }
   }
   cout<< "\n\nCentrality Percentile Curves      #"<< m << endl;
   cout<< "Centrality Estimator and Classes  #"<< n << endl;
   cout<< "--------------------------------------------------------------------------------------"  << endl;

   return 0;
}

Bool_t HParticleEvtCharaBK::saveCentralityEstimatorHist()
{
    if(!fFile) return kFALSE;
    fFile->mkdir("EstimatorHist");
    fFile->cd("/EstimatorHist/");

    cout<<"\n--------------------------------------------------------------------------------------"  << endl;
    Info("saveCentralityEstimatorHist()","Estimator Hist saving:");
    Int_t m=0;
    for (int centEst = 0; centEst < (int) kNumCentralityEstimator; ++centEst){
       if(fEstimatorHist[centEst]){
            fEstimatorHist[centEst]->Write();
            m++;
        }
   }
   cout<< "\n\nEstimator Hist      #"<< m << "saved" << endl;
   cout<< "--------------------------------------------------------------------------------------"  << endl;

    fFile->mkdir("Centrality");
    fFile->cd("/Centrality/");
    
    cout<<"\n--------------------------------------------------------------------------------------"  << endl;
    Info("saveCentralityEstimatorHist()","Calibration for Centrality Estimator and Classes saving:");
    Int_t n=0;
    m=0;
    for (int centEst = 0; centEst < (int) kNumCentralityEstimator; ++centEst){
       printf("\n>>>> %31s %s:  ",getStringCentralityEstimator(centEst).Data() ,fPostFix.Data());
       if(fCentralityPercentileHist[centEst]){
            fCentralityPercentileHist[centEst]->Write();
            printf("percentile");
            m++;
        }
       
       for (int centC = 0; centC < (int) kNumCentralityClass; ++centC){
           if(fCentralityHist[centEst][centC]){
                fCentralityHist[centEst][centC]->Write();
                printf(" %5s",getStringCentralityClass(centC).Data() );
                n++;
           }
       }
   }
   cout<< "\n\nCentrality Percentile Curves      #"<< m << "saved" << endl;
   cout<< "Centrality Estimator and Classes  #"<< n << "saved" << endl;
   cout<< "--------------------------------------------------------------------------------------"  << endl;
   return 0;
}


Bool_t HParticleEvtCharaBK::loadEventPlaneCorrectionHist()
{

    cout<<"\n--------------------------------------------------------------------------------------"  << endl;
    Info("loadEventPlaneCorrectionHist()","Calibration for EventPlane Correction Histogramms loading:");
    Int_t n=0;
    for (int epParam = 0; epParam < (int) kNumEventPlaneParameter; ++epParam){
       printf("\n>>>> %27s :  ",getStringEventPlaneParameter(epParam).Data() );
       TString temp;
       temp = Form("EPcorr_%s_Day_Centrality",getStringEventPlaneParameter(epParam).Data() );
       fEventPlaneCorrectionHist[epParam] = (TProfile2D*) fFile->FindObjectAny( temp.Data() );   //FIXME
       if(fEventPlaneCorrectionHist[epParam]){
            //for (int dim = 0; dim < fEventPlaneCorrectionHist[qAxis]->GetDimension(); ++dim){
                printf("%4s(%2d bins) ",fEventPlaneCorrectionHist[epParam]->GetXaxis()->GetTitle(),
                                  fEventPlaneCorrectionHist[epParam]->GetNbinsX());
                printf("%4s(%2d bins) ",fEventPlaneCorrectionHist[epParam]->GetYaxis()->GetTitle(),
                                  fEventPlaneCorrectionHist[epParam]->GetNbinsY());
                printf("%4s(%2d bins) ",fEventPlaneCorrectionHist[epParam]->GetZaxis()->GetTitle(),
                                  fEventPlaneCorrectionHist[epParam]->GetNbinsZ());
            n++;
        }
   }

   cout<< "\n\nEventPlane Correction Histogramms  #"<< n << endl;
   cout<< "--------------------------------------------------------------------------------------"  << endl;

   return 0;
}

Bool_t HParticleEvtCharaBK::addEventPlaneCorrectionHist(TProfile2D *hist, UInt_t epParam)
{
    if(!hist)return kFALSE;
    if(epParam>=kNumEventPlaneParameter) return kFALSE;
    if(epParam==kChi)  return kTRUE;
    
    TString temp;
    temp = Form("EPcorr_%s_Day_Centrality",getStringEventPlaneParameter(epParam).Data() );
    hist->SetNameTitle(temp.Data(),temp.Data());

    if(epParam==kResolution){
         fEventPlaneCorrectionHist[kResolution] = (TH2D*) makeEPresolution(hist);

         temp = Form("EPcorr_%s_Day_Centrality",getStringEventPlaneParameter(kChi).Data() );
         hist->SetNameTitle(temp.Data(),temp.Data());
         fEventPlaneCorrectionHist[kChi] = (TH2D*) makeEPresolution(hist, kTRUE);
    }
    else{
         fEventPlaneCorrectionHist[epParam] = (TProfile2D*) hist;
    }
    return kTRUE;
}

Bool_t HParticleEvtCharaBK::saveEventPlaneCorrectionHist()
{
    if(!fFile) return kFALSE;
    fFile->mkdir("EPcorr");
    fFile->cd("/EPcorr/");
    
    cout<<"\n--------------------------------------------------------------------------------------"  << endl;
    Info("saveEventPlaneCorrectionHist()","Calibration for EventPlane Correction Histogramms saving:");
    Int_t n=0;
    for (int epParam = 0; epParam < (int) kNumEventPlaneParameter; ++epParam){
       printf("\n>>>> %27s :  ",getStringEventPlaneParameter(epParam).Data() );
       if(fEventPlaneCorrectionHist[epParam]){
           fEventPlaneCorrectionHist[epParam]->Write();
           printf("%4s(%2d bins) ",fEventPlaneCorrectionHist[epParam]->GetXaxis()->GetTitle(),
                             fEventPlaneCorrectionHist[epParam]->GetNbinsX());
           printf("%4s(%2d bins) ",fEventPlaneCorrectionHist[epParam]->GetYaxis()->GetTitle(),
                             fEventPlaneCorrectionHist[epParam]->GetNbinsY());
           printf("%4s(%2d bins) ",fEventPlaneCorrectionHist[epParam]->GetZaxis()->GetTitle(),
                             fEventPlaneCorrectionHist[epParam]->GetNbinsZ());
            n++;
        }
   }

   cout<< "\n\nEventPlane Correction Histogramms  #"<< n << "saved" << endl;
   cout<< "--------------------------------------------------------------------------------------"  << endl;

   return 0;
}

Bool_t HParticleEvtCharaBK::loadFWCutValuesHist()
{
    cout<<"\n--------------------------------------------------------------------------------------"  << endl;
    Info("loadFWCutValuesHist()","Histogramms with CutValues for FW Hits loading:");
    Int_t n=0;
    for (int ep = 0; ep < (int) kNumFWCutValues; ++ep){
       Int_t m = 0;
       printf("\n>>>> %27s :  ",getStringFWCutValues(ep).Data() );
       for (int iCell = 0; iCell < MAXFWCELLS ; ++iCell){
           TString temp = Form("FWCuts_%s_cell%d",getStringFWCutValues(ep).Data(), iCell );
           fFWCutValuesHist[ep][iCell] = (TH1F*) fFile->FindObjectAny( temp.Data() );   //FIXME
           if(fFWCutValuesHist[ep][iCell]){
               if(m==0){
                   printf("%4s(%2d bins) ",fFWCutValuesHist[ep][iCell]->GetXaxis()->GetTitle(),
                                     fFWCutValuesHist[ep][iCell]->GetNbinsX());
                   printf("%4s(%2d bins) ",fFWCutValuesHist[ep][iCell]->GetYaxis()->GetTitle(),
                                     fFWCutValuesHist[ep][iCell]->GetNbinsY());
               }
               n++; m++;
            }
       }
       printf("active cells:%2d ", m);
       if(m==0){
            useFWCut[ep] = kFALSE;
            printf("  not used for FWCut!");
       }
   }

   cout<< "\n\nHistogramms with CutValues for FW Hits  #"<< n << endl;
   cout<< "--------------------------------------------------------------------------------------"  << endl;

   return 0;
}

Bool_t HParticleEvtCharaBK::addFWCutValuesHist(TH1 *hist, Int_t cell, UInt_t eFWCut)
{
    if(!hist)return kFALSE;
    if(eFWCut>=kNumFWCutValues) return kFALSE;
    if(cell<0 || cell >= MAXFWCELLS) return kFALSE;
    
    TString temp = Form("FWCuts_%s_cell%d",getStringFWCutValues(eFWCut).Data(), cell );
    hist->SetNameTitle(temp.Data(),temp.Data());
    fFWCutValuesHist[eFWCut][cell] = (TH1F*) hist;
    return kTRUE;
}

Bool_t HParticleEvtCharaBK::saveFWCutValuesHist()
{
    if(!fFile) return kFALSE;
    fFile->mkdir("FWCutValue");
    fFile->cd("/FWCutValue/");
    
    cout<<"\n--------------------------------------------------------------------------------------"  << endl;
    Info("saveFWCutValuesHist()","Histogramms with CutValues for FW Hits saving:");
    Int_t n=0;
    for (int ep = 0; ep < (int) kNumFWCutValues; ++ep){
       Int_t m = 0;
       printf("\n>>>> %27s :  ",getStringFWCutValues(ep).Data() );
       for (int iCell = 0; iCell < (int) MAXFWCELLS; ++iCell){
           if(fFWCutValuesHist[ep][iCell]){
               TString temp = Form("FWCuts_%s_cell%d",getStringFWCutValues(ep).Data(), iCell );
               fFWCutValuesHist[ep][iCell]->SetNameTitle(temp.Data(),temp.Data());
               fFWCutValuesHist[ep][iCell]->Write();
               if(m==0){
                   printf("%4s(%2d bins) ",fFWCutValuesHist[ep][iCell]->GetXaxis()->GetTitle(),
                                     fFWCutValuesHist[ep][iCell]->GetNbinsX());
                   printf("%4s(%2d bins) ",fFWCutValuesHist[ep][iCell]->GetYaxis()->GetTitle(),
                                     fFWCutValuesHist[ep][iCell]->GetNbinsY());
               }
               n++; m++;
            }
       }
       printf("active cells:%2d ", m);

   }

   cout<< "\n\ntotal Histogramms with CutValues for FW Hits  #"<< n << "saved" << endl;
   cout<< "--------------------------------------------------------------------------------------"  << endl;

   return 0;
}


Int_t HParticleEvtCharaBK::execute(){
    return 0;
}

Int_t HParticleEvtCharaBK::print(){
    return 0;
}

void HParticleEvtCharaBK::reset()
{
    // Reset.
}
//----------------------------------------------------------------------
Bool_t HParticleEvtCharaBK::isNewEvent()
{
    if(gHades->getCurrentEvent()->getHeader()->getEventSeqNumber() == currentEventSeqNumber){
        return kFALSE;
    } 
    else{
        currentEventSeqNumber = gHades->getCurrentEvent()->getHeader()->getEventSeqNumber();
        loadDayOfYear();
        return kTRUE;
    }
}

Float_t HParticleEvtCharaBK::getEventWeight()
{
    if(isSimulation) return 1;
// ---- PT3 Event without PT2 Events --------------
    if(gHades->getCurrentEvent()->getHeader()->getTBit() == 8192) return 1;
// ---- PT2 Event with PT3 Events --------------
    if(gHades->getCurrentEvent()->getHeader()->getTBit() == 12288)return 1;
// ---- PT2 Event without PT3 Events ----------
    if(gHades->getCurrentEvent()->getHeader()->getTBit() == 4096){
        if(fDayOfYear<=106) return 4; //see logbook entry 15.4.2012 ~23.52h  from MultSectra: 3.7;
        else                return 8; //   7.8;
    }
    return 0;
}


//----------------------------------------------------------------------

Int_t HParticleEvtCharaBK::GetFWmoduleSize(HWallHitSim* wall)
{
    Int_t i = wall->getCell();
    
    if ( i < 0 || i > 301 ) return -1;
    else if ( i==65 || i==66 || i==77 || i==78 ) return -1;      //beam hole
    else if ( i<144 ) return 1;    // small modules
    else if ( i<208 ) return 2;    // middle modules
    else if ( i<302 ) return 3;    // large modules
    return -1;
}


// Bool_t HParticleEvtCharaBK::PassesCutsFW(HWallHitSim* wall)
// {
//     if(!wall) return kFALSE;
//     Float_t betaFW  =  wall->getDistance() / wall->getTime()/ 299.792458;
//     
//      if(!isSimulation){
//          if( (GetFWmoduleSize(wall) == 1
//               &&  wall->getCharge()>fFWminCharge[0]&&betaFW>fFWminBeta[0] && betaFW <fFWmaxBeta[0] ) //small
//     	  || (GetFWmoduleSize(wall) == 2
//               &&  wall->getCharge()>fFWminCharge[1]&&betaFW>fFWminBeta[1] && betaFW <fFWmaxBeta[1] ) //medium
//     	  || (GetFWmoduleSize(wall) == 3
//               &&  wall->getCharge()>fFWminCharge[2]&&betaFW>fFWminBeta[2] && betaFW <fFWmaxBeta[2] )  //large
//               ) return kTRUE;	  
//      }
//      else
//      {
//          if( (GetFWmoduleSize(wall) == 1  
//               &&  betaFW>fFWminBeta[0] && betaFW <fFWmaxBeta[0] ) //small
//           || (GetFWmoduleSize(wall) == 2  
//               &&  betaFW>fFWminBeta[1] && betaFW <fFWmaxBeta[1] ) //medium
//           || (GetFWmoduleSize(wall) == 3  
//               &&  betaFW>fFWminBeta[2] && betaFW <fFWmaxBeta[2] )  //large
//               ) return kTRUE;
//      }
//      return kFALSE;
// }

Bool_t HParticleEvtCharaBK::PassesCutsFW(HWallHitSim* wall)
{
    for (int ep = 0; ep < (int) kNumFWCutValues; ++ep){
        if(useFWCut[ep] && !PassesCutsFW(wall, ep) ) return kFALSE;
    }
    return kTRUE;
}

Bool_t HParticleEvtCharaBK::PassesCutsFW(HWallHitSim* wall, UInt_t eFWCut)
{
    if(!wall) return kFALSE;
    Int_t cell      =  wall->getCell();
    if(cell<0 || cell>= MAXFWCELLS) return kFALSE;
    
    if(excludeNoisyFWcells){
        if(cell ==20 || cell ==53 || cell ==54 
            || cell ==64 || cell ==67 || cell ==79) return kFALSE;      //noisy cells
    }
    
    Float_t value =0;
    if(eFWCut==kBetaCuts  ) value  =  wall->getDistance() / wall->getTime()/ 299.792458;
    if(eFWCut==kTimeCuts  ) value  =  wall->getTime();
    if(eFWCut==kChargeCuts) value  =  wall->getCharge();
    

    if(fFWCutValuesHist[eFWCut][cell]){
        if(fFWCutValuesHist[eFWCut][cell]->GetBinContent(fFWCutValuesHist[eFWCut][cell]->FindBin(value)) > 0) return kTRUE;
    }
    return kFALSE;
}

Int_t HParticleEvtCharaBK::getFWCharge(HWallHitSim* wall)
{
    Float_t Charge = 0;
    Int_t Z = 0;

    if(!wall) return 0;
    if(!isSimulation) Charge = wall->getCharge();
    else              Charge = 93.*pow(wall->getCharge(),0.46-0.006*sqrt(wall->getCharge()));  // parametrization from R.Holzmann
    Int_t cell = wall->getCell();
    if(cell<0 || cell>= MAXFWCELLS) return 0;

    if(fFWCutValuesHist[kChargeCuts][cell]){
        Z = fFWCutValuesHist[kChargeCuts][cell]->GetBinContent(fFWCutValuesHist[kChargeCuts][cell]->FindBin(Charge) );
    }
    if(Z>1) Z=Z-1;
    
    return Z;
}
/*
Int_t HParticleEvtCharaBK::getFWCharge(HWallHitSim* wall)
{
    Float_t Charge = 0;
    if(!wall) return kFALSE;
    if(!isSimulation) Charge = wall->getCharge();
    else              Charge = 93.*pow(wall->getCharge(),0.46-0.006*sqrt(wall->getCharge()));  // parametrization from R.Holzmann

    if(Charge>=fFWChargeCuts[5])       return 6;
    else if(Charge>=fFWChargeCuts[4])  return 5;
    else if(Charge>=fFWChargeCuts[3])  return 4;
    else if(Charge>=fFWChargeCuts[2])  return 3;
    else if(Charge>=fFWChargeCuts[1])  return 2;
    else if(Charge>=fFWChargeCuts[0])  return 1;
    else return 0;
    
}
*/

TH1F* HParticleEvtCharaBK::getEventPlaneCorrectionHist(UInt_t flag) const
{
    if(flag<kNumEventPlaneParameter) return (TH1F*)fEventPlaneCorrectionHist[flag];
    else return 0;
}


Float_t HParticleEvtCharaBK::getCorrection(UInt_t flag)
{
    Float_t fCorrection = 0.;
    Double_t cent = getCentralityPercentile();
    if(getEventPlaneCorrectionHist(flag)) 
        fCorrection = getEventPlaneCorrectionHist(flag)->GetBinContent(getEventPlaneCorrectionHist(flag)->FindBin(cent, fDayOfYear) );
    return fCorrection;
}

Float_t HParticleEvtCharaBK::getCorrectionError(UInt_t flag)
{
    Float_t fCorrectionError = 1.;
    Double_t cent = getCentralityPercentile();
    if(getEventPlaneCorrectionHist(flag)) 
        fCorrectionError = getEventPlaneCorrectionHist(flag)->GetBinError(getEventPlaneCorrectionHist(flag)->FindBin(cent, fDayOfYear) );
    return fCorrectionError;
}

Float_t HParticleEvtCharaBK::getSmearValue(Int_t size){
    if(size==1) return ( fRandom->Rndm(1) - 0.5 )*40.;
    else if(size==2) return ( fRandom->Rndm(1) - 0.5 )*80.;
    else if(size==3) return ( fRandom->Rndm(1) - 0.5 )*160.;
    return 0;
}

Float_t HParticleEvtCharaBK::getThetaWeight(HWallHitSim* wall, Float_t min, Float_t max){
     Float_t theta = wall->getTheta();
     if(theta>max) return 0.;    // max Theta 
     if(theta<min) return 1.;    // min Theta
     Float_t wtheta = 1. -( (theta-min)/(max-min) );   
     if(wtheta>0 && wtheta < 1.) return wtheta;
     else return 0.;
}

void HParticleEvtCharaBK::fillHitArray()
{
    //  -----------------------------------------------------------------
    //fill Hits into HitArray for QVector
    fQVectorCalcDone = kFALSE; //reset QVectors
    fEventPlaneCorrectionFlag = kNoCorrection;
    Int_t nHits = arrayOfHits.size();
    for (Int_t i=0; i<nHits; i++){
        if (!arrayOfHits[i]) continue;
        delete arrayOfHits[i];
        arrayOfHits[i] = NULL;
    }
    iFWHitvector.clear();
    arrayOfHits.clear();
    
    if(!fCatWallHit) return;
    
    HWallHitSim *wall = 0;  // includes HWallHit
    Float_t wallXOrg=0,wallYOrg=0,wallZOrg=0;
    Float_t wallX=0,wallY=0,wallZ=0;
    Float_t weight1 = 1.;
    Float_t weight2 = 1.;
    for(Int_t i = 0; i < fCatWallHit->getEntries(); i ++ ){
        wall  = (HWallHitSim*) HCategoryManager::getObject(wall,fCatWallHit,i);
        if(!PassesCutsFW(wall)) continue;
         iFWHitvector.push_back(i);
        
         wall->getXYZLab(wallXOrg,wallYOrg,wallZOrg);
        
        //if(isFlagSet(kSmear, statusflag)){
        //      Int_t size = GetFWmoduleSize(wall);
        //      wallX += getSmearValue(size);
        //      wallY += getSmearValue(size);
        //}
        
        wallX = (wallXOrg - getCorrection(kFWx)) / getCorrectionError(kFWx);
        wallY = (wallYOrg - getCorrection(kFWy)) / getCorrectionError(kFWy);
        weight1 = getFWCharge(wall);
        weight2 = getThetaWeight(wall);

        //FIXME   catch x=0 y=0 
        if(wallX==0 || wallY==0) continue;
        //if(weight1<=0 || weight1>1) continue;
        //if(weight2<=0 || weight2>1) continue;

        SimpleQVector* HitVector = new SimpleQVector();
        HitVector->Set(wallX,wallY,weight1, weight2); 
        HitVector->SetOrg(wallXOrg,wallYOrg); 


        //HitVector.Weight = weight;
        arrayOfHits.push_back(HitVector);
        
    } 
    nHits = arrayOfHits.size();
    if(nHits<1) return;
    //shuffle
    std::random_shuffle ( arrayOfHits.begin(), arrayOfHits.end() );
    //Flag two subevents in given shuffled sub-events
    for (Int_t i=0; i<nHits; i++){
        //SimpleQVector* hit = static_cast<SimpleQVector*>(arrayOfHits->At(i));
        if (!arrayOfHits[i]) continue;
        if((i+1)*2 > nHits)  arrayOfHits[i]->SetSubEvent(1);
        else                 arrayOfHits[i]->SetSubEvent(2);
        
       // printf("i:%d  Phi: %.3f, Weight1: %.3f, Weight2: %.3f, SubEvent: %d",
       //          i,  arrayOfHits[i]->Phi(),arrayOfHits[i]->Weight1(),arrayOfHits[i]->Weight2(),arrayOfHits[i]->SubEvent() );
       // printf("\n");
    }
    return;
}

void HParticleEvtCharaBK::printHitArray()
{
    //  -----------------------------------------------------------------
    if(isNewEvent()){
        fillHitArray();
    }
    //print Hits of HitArray
    Int_t nHits = arrayOfHits.size();
    printf( "printHitArray #Hits = %d\n", nHits );
    if(nHits<1) return;
    for (Int_t i=0; i<nHits; i++){
        if (arrayOfHits[i]) arrayOfHits[i]->print();
    }
    return;
}

vector<Int_t> HParticleEvtCharaBK::getFWhits(){
    if(isNewEvent()){
        fillHitArray();
    }
    return iFWHitvector;
}

Bool_t HParticleEvtCharaBK::fillQVectors(UInt_t statusflag, UInt_t nHarmonic)
{
    //  -----------------------------------------------------------------
    //HWallHit  into QVector
 
    fEventPlaneCorrectionFlag = kNoCorrection;
    fQVectorCalcDone = kFALSE; 

    vQPhi[0] = -1;
    vQPhi[1] = -1;
    vQPhi[2] = -1;

    Double_t dQX  = 0.,dQY  = 0.;
    Double_t dQXA = 0.,dQYA = 0.;
    Double_t dQXB = 0.,dQYB = 0.;

    Double_t dQXShift  = 0.,dQYShift  = 0.;
    Double_t dQXScale  = 1.,dQYScale  = 1.;

    Double_t phi = -1;
    Double_t w = 1.;

    Double_t sumOfWeights=0;
    //Double_t nHits=0;

    Int_t nHits = arrayOfHits.size();
    if(nHits<4) return kFALSE;  // at least two hits per sub-event for  
    for (Int_t i=0; i<nHits; i++){
        w = 1.;
        if(isFlagSet(kShiftFW, statusflag)) phi = arrayOfHits[i]->Phi();
        else                                phi = arrayOfHits[i]->PhiOrg();

        if(isFlagSet(kWeightCharge, statusflag)) w = arrayOfHits[i]->Weight1();
        if(isFlagSet(kWeightTheta, statusflag))  w *= arrayOfHits[i]->Weight2();
        
        dQX += w * TMath::Cos(nHarmonic*phi);
        dQY += w * TMath::Sin(nHarmonic*phi);
        
        if(arrayOfHits[i]->SubEvent() == 1) {
            dQXA += w * TMath::Cos(nHarmonic*phi);
            dQYA += w * TMath::Sin(nHarmonic*phi);
        }
        else if(arrayOfHits[i]->SubEvent() == 2) {
            dQXB += w * TMath::Cos(nHarmonic*phi);
            dQYB += w * TMath::Sin(nHarmonic*phi);
        }
            sumOfWeights +=w;
    }
//   printf("nHits: %d dQX: %.3f, dQY: %.3f \t  A:  dQX: %.3f, dQY: %.3f \t B:  dQX: %.3f, dQY: %.3f \t b weight: %.3f \n",
//                 nHits, dQX, dQY , dQXA, dQYA ,dQXB, dQYB, sumOfWeights );


        if(isFlagSet(kReCentering, statusflag)){
            if(isFlagSet(kShiftFW, statusflag)){
                if(isFlagSet(kWeightCharge, statusflag)){ dQXShift  = getCorrection(kQx2WCharge); dQYShift  = getCorrection(kQy2WCharge);}
                else{                                     dQXShift  = getCorrection(kQx2);        dQYShift  = getCorrection(kQy2);}
            }
            else{
                if(isFlagSet(kWeightCharge, statusflag)){ dQXShift  = getCorrection(kQxWCharge); dQYShift  = getCorrection(kQyWCharge);}
                else{                                     dQXShift  = getCorrection(kQx);        dQYShift  = getCorrection(kQy);}
            }
            
            dQX  -= dQXShift;
            dQY  -= dQYShift;
            dQXA -= dQXShift;
            dQYA -= dQYShift;
            dQXB -= dQXShift;
            dQYB -= dQYShift;
            if(isFlagSet(kScaling, statusflag)){
                if(isFlagSet(kShiftFW, statusflag)){
                    if(isFlagSet(kWeightCharge, statusflag)){dQXScale  = getCorrectionError(kQx2WCharge); dQYScale  = getCorrectionError(kQy2WCharge);}
                    else{                                    dQXScale  = getCorrectionError(kQx2); dQYScale  = getCorrectionError(kQy2);}
                }
                else{
                    if(isFlagSet(kWeightCharge, statusflag)){dQXScale  = getCorrectionError(kQxWCharge); dQYScale  = getCorrectionError(kQyWCharge);}
                    else{                                    dQXScale  = getCorrectionError(kQx); dQYScale  = getCorrectionError(kQy);}
                }
                
                dQX  /= dQXScale;
                dQY  /= dQYScale;
                dQXA /= dQXScale;
                dQYA /= dQYScale;
                dQXB /= dQXScale;
                dQYB /= dQYScale;
            }
        }
    //    printf("n nHits: %d dQXShift: %.3f, dQYShift: %.3f \t  dQXScale: %.3f, dQYScale: %.3f \t Corr: %s \n",
    //               nHits, dQXShift, dQYShift , dQXScale, dQYScale, getStringEventPlaneCorrection(statusflag).Data() );
    //    printf("n nHits: %d dQX: %.3f, dQY: %.3f \t A:  dQX: %.3f, dQY: %.3f \t B:  dQX: %.3f, dQY: %.3f \t Corr: %s \n",
    //               nHits, dQX, dQY , dQXA, dQYA ,dQXB, dQYB, getStringEventPlaneCorrection(statusflag).Data() );
    //    printf("\n");

        if(dQX==0 && dQX==0)   return kFALSE;
        if(dQXA==0 && dQXA==0) return kFALSE;
        if(dQXB==0 && dQXB==0) return kFALSE;



   //vQ ->Set(dQX ,dQY );
   //vQA->Set(dQXA,dQYA);
   //vQB->Set(dQXB,dQYB);
    
   if(isFlagSet(kRotation, statusflag)){
       Float_t corrPsi = getCorrectionPhi(getPhi(dQX ,dQY ));
       vQPhi[0] = TVector2::Phi_0_2pi( getPhi(dQX ,dQY ) + corrPsi );
       vQPhi[1] = TVector2::Phi_0_2pi( getPhi(dQXA,dQYA) + corrPsi );
       vQPhi[2] = TVector2::Phi_0_2pi( getPhi(dQXB,dQYB) + corrPsi );
   }
   else{
       vQPhi[0] = getPhi(dQX ,dQY );
       vQPhi[1] = getPhi(dQXA,dQYA);
       vQPhi[2] = getPhi(dQXB,dQYB);
   }
    
    fEventPlaneCorrectionFlag = statusflag;
    fQVectorCalcDone = kTRUE; 
    
    return kTRUE;
}

void HParticleEvtCharaBK::printQVectors()
{
    //  -----------------------------------------------------------------
    if(fQVectorCalcDone){ printf("QVector Calculation Done\n");}
    else{
        printf("QVector Calculation NOT Done !!!!!!! retry now\n");
        if(!fillQVectors(fEventPlaneCorrectionFlag)) {printf("QVector Calculation failed !!!!!!!\n"); return ;}
    }

    printf("##### vQPhi: %.3f \t vQPhiA: %.3f \t vQPhiB: %.3f \t\t Corr: %s #####\n",
               vQPhi[0], vQPhi[1], vQPhi[2],  getStringEventPlaneCorrection(fEventPlaneCorrectionFlag).Data() );
    return;
}


Float_t HParticleEvtCharaBK::getCorrectionPhi(Float_t phi){
    
  // for(UInt_t k = 0; k < charm; k++){
  //   cpphi += An[k]*cos((Double_t)(k+1) * val);
  //   cpphi += Bn[k]*sin((Double_t)(k+1) * val);
  // } 

    //used Formular in Footnoot1 of Phys.Rev. C58 (1998) 1671-1678
    //also appendix of Phys.Rev. C56 (1997) 3254-3264
    // kFourierCn  is equivalent to Bn
    // kFourierSn  is equivalent to An

    Float_t corrPsi = phi;
    corrPsi =    2./1.*(getCorrection(kFourierC1)*TMath::Sin(1.*phi)) 
               + 2./2.*(getCorrection(kFourierC2)*TMath::Sin(2.*phi)) 
               + 2./3.*(getCorrection(kFourierC3)*TMath::Sin(3.*phi)) 
               + 2./4.*(getCorrection(kFourierC4)*TMath::Sin(4.*phi))
               + 2./5.*(getCorrection(kFourierC5)*TMath::Sin(5.*phi))
               + 2./6.*(getCorrection(kFourierC6)*TMath::Sin(6.*phi))
               + 2./7.*(getCorrection(kFourierC7)*TMath::Sin(7.*phi))
               + 2./8.*(getCorrection(kFourierC8)*TMath::Sin(8.*phi))

               - 2./1.*(getCorrection(kFourierS1)*TMath::Cos(1.*phi))
               - 2./2.*(getCorrection(kFourierS2)*TMath::Cos(2.*phi))
               - 2./3.*(getCorrection(kFourierS3)*TMath::Cos(3.*phi))
               - 2./4.*(getCorrection(kFourierS4)*TMath::Cos(4.*phi))
               - 2./5.*(getCorrection(kFourierS5)*TMath::Cos(5.*phi))
               - 2./6.*(getCorrection(kFourierS6)*TMath::Cos(6.*phi))
               - 2./7.*(getCorrection(kFourierS7)*TMath::Cos(7.*phi))
               - 2./8.*(getCorrection(kFourierS8)*TMath::Cos(8.*phi));
    return corrPsi;
}

Float_t HParticleEvtCharaBK::getEventPlane(UInt_t statusflag, UInt_t SubEvent, UInt_t nHarmonic)
{
    if(isNewEvent()){
        fillHitArray();
    }
    if(!fQVectorCalcDone){
        if(!fillQVectors(statusflag)) return -1;
    }
    
    if(fQVectorCalcDone && fEventPlaneCorrectionFlag == statusflag){
        if(SubEvent<3)  return vQPhi[SubEvent];
        else return -1;
    }
    else{
        if(!fillQVectors(statusflag, nHarmonic)) return -1;
        if(fQVectorCalcDone && fEventPlaneCorrectionFlag == statusflag){
            if(SubEvent<3)  return vQPhi[SubEvent];
            else return -1;
        }
    }
    return -1;

}

Float_t HParticleEvtCharaBK::getEventPlaneWeight(UInt_t statusflag, UInt_t SubEvent, UInt_t nHarmonic)
{
    return 1.;
}

Float_t HParticleEvtCharaBK::getEventPlaneParameter(UInt_t e, Bool_t corr)
{
         
    if(isNewEvent()){
        fillHitArray();
    }
    if(e==kFWx|| e==kFWy)                      return -999;  

    if(corr){
         if(e==kQx  || e==kQy)                                          fillQVectors(kReCentering|kScaling);
         if(e==kQxWCharge  || e==kQyWCharge)                            fillQVectors(kReCentering|kScaling|kWeightCharge);
         if(e==kQx2 || e==kQy2 || e==kResolution)                       fillQVectors(kShiftFW|kReCentering|kScaling);
         if(e==kQx2WCharge || e==kQy2WCharge || e==kResolutionWCharge)  fillQVectors(kShiftFW|kReCentering|kScaling|kWeightCharge);
         if(e>=kFourierC1 && e<=kFourierS8)                             fillQVectors(kShiftFW|kWeightCharge|kRotation);
    }
    else{
         if(e==kQx  || e==kQy)                                          fillQVectors(kNoCorrection);
         if(e==kQxWCharge  || e==kQyWCharge)                            fillQVectors(kWeightCharge);
         if(e==kQx2 || e==kQy2|| e==kResolution)                        fillQVectors(kShiftFW);
         if(e==kQx2WCharge || e==kQy2WCharge || e==kResolutionWCharge)  fillQVectors(kShiftFW|kWeightCharge);
         if(e>=kFourierC1 && e<=kFourierS8)                             fillQVectors(kShiftFW|kWeightCharge);
    }
    if(vQPhi[0]==-1) return -999;
    
    if(e==kQx || e==kQxWCharge || e==kQx2 || e==kQx2WCharge)      return TMath::Cos(    vQPhi[0] );
    if(e==kQy || e==kQyWCharge || e==kQy2 || e==kQy2WCharge)      return TMath::Sin(    vQPhi[0] );
    if(e==kFourierC1)                return TMath::Cos(    vQPhi[0] );
    if(e==kFourierC2)                return TMath::Cos(2.* vQPhi[0] );
    if(e==kFourierC3)                return TMath::Cos(3.* vQPhi[0] );
    if(e==kFourierC4)                return TMath::Cos(4.* vQPhi[0] );
    if(e==kFourierC5)                return TMath::Cos(5.* vQPhi[0] );
    if(e==kFourierC6)                return TMath::Cos(6.* vQPhi[0] );
    if(e==kFourierC7)                return TMath::Cos(7.* vQPhi[0] );
    if(e==kFourierC8)                return TMath::Cos(8.* vQPhi[0] );
    if(e==kFourierS1)                return TMath::Sin(    vQPhi[0] );
    if(e==kFourierS2)                return TMath::Sin(2.* vQPhi[0] );
    if(e==kFourierS3)                return TMath::Sin(3.* vQPhi[0] );
    if(e==kFourierS4)                return TMath::Sin(4.* vQPhi[0] );
    if(e==kFourierS5)                return TMath::Sin(5.* vQPhi[0] );
    if(e==kFourierS6)                return TMath::Sin(6.* vQPhi[0] );
    if(e==kFourierS7)                return TMath::Sin(7.* vQPhi[0] );
    if(e==kFourierS8)                return TMath::Sin(8.* vQPhi[0] );
    if(e==kResolution || e==kResolutionWCharge){  //FIXME  only 1 higher harmonics not included!!!
       return TMath::Cos(vQPhi[1])*TMath::Cos(vQPhi[2]) + TMath::Sin(vQPhi[1])*TMath::Sin(vQPhi[2]);
    }
    return -999;
}


//----------------------------------------------------------------------
Int_t HParticleEvtCharaBK::getCentralityClass(TString estimator)
{
    // legacy code
    // Return centrality class, default in 5% of total cross section with estimator
    // or with preset classes like FOPI {6.3%, 21.0%, 30.9%}

    if (estimator.CompareTo("TOFRPCtimecut")==0 || estimator.CompareTo("TOFRPC5")==0) 
        return getCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k5);
    else if (estimator.CompareTo("TOFRPCtimecutFOPI")==0 || estimator.CompareTo("TOFRPCFOPI")==0) 
        return getCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::kFOPI);
    else if (estimator.CompareTo("TOFRPCtimecut10")==0 || estimator.CompareTo("TOFRPC10")==0) 
        return getCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k10);
    else {
        Error("getCentralityClass()","No CentralityEstimator defined!");
        return 0;
    }
}

Float_t HParticleEvtCharaBK::getCentralityEstimator(UInt_t centE)
{
    HParticleEvtInfo *event_info = (HParticleEvtInfo*)fParticleEvtInfoCat->getObject(0);
    if(!event_info) {
        Error("getCentralityClass()","No HParticleEvtInfo");
    return 0;
    }
    
    if(centE     ==kTOFRPC)                      return event_info->getSumTofMultCut() + event_info->getSumRpcMultHitCut();
    else if(centE==kTOF)                         return event_info->getSumTofMultCut();
    else if(centE==kRPC)                         return event_info->getSumRpcMultHitCut();
    else if(centE==kTOFRPCtot)                   return event_info->getSumTofMult() + event_info->getSumRpcMultHit();
    else if(centE==kTOFtot)                      return event_info->getSumTofMult();
    else if(centE==kRPCtot)                      return event_info->getSumRpcMultHit();
    else if(centE==kSelectedParticleCand)        return event_info->getSumSelectedParticleCandMult();
    else if(centE==kSelectedParticleCandCorr)    return event_info->getSumSelectedParticleCandMult()*(fReferenceMeanSelTrack/event_info->getMeanMult());    
    else if(centE==kSelectedParticleCandNorm)    return event_info->getSumSelectedParticleCandMult()/event_info->getMeanMult();   
    else if(centE==kSelectedParticleCandSecNorm) return getSelectedParticleCandSecNorm();    
    else if(centE==kSelectedParticleCandCorrPerWire) return getSelectedParticleCandCorrPerWire();   
    else if(centE==kPrimaryParticleCand)         return event_info->getSumPrimaryParticleCandMult();
    else if(centE==kMdcWires)                    return event_info->getMdcWires();
    else if(centE==kMdcWiresOuterMod)            return getMdcWiresOuterMod();
    else if(centE==kFWSumChargeSpec)             return getFWSumChargeSpec();
    else if(centE==kFWSumChargeZ)                return getFWSumZ();
    else if(centE==kDirectivity)                 return getDirectivity();
    else if(centE==kRatioEtEz)                   return getRatioEtEz();
    else if(centE==kEt)                          return getEt();
    else return 0;
}

TH1F* HParticleEvtCharaBK::getCentralityClassHist(UInt_t centE, UInt_t centC) const
{
    if(centE<kNumCentralityEstimator && centC<kNumCentralityClass) return (TH1F*)fCentralityHist[centE][centC];
    else return 0;
}

TH1F* HParticleEvtCharaBK::getCentralityPercentileHist(UInt_t centE) const
{
    if(centE<kNumCentralityEstimator) return (TH1F*)fCentralityPercentileHist[centE];
    else return 0;
}

Int_t HParticleEvtCharaBK::getCentralityClass(UInt_t centE, UInt_t centC)
{
    // Return centrality class, default in 5% of total cross section with estimator
    // or with preset classes like FOPI {6.3%, 21.0%, 30.9%}
    if(centE>=kNumCentralityEstimator){
        return 101.;
    }
    else if(!getCentralityClassHist(centE, centC)){
        //Error("getCentralityClass()","No CentralityEstimator defined!");
        return 0;
    }
    else {
        Int_t buffer = getCentralityEstimator(centE);
        Int_t bin = getCentralityClassHist(centE, centC)->FindBin(buffer);
        return getCentralityClassHist(centE, centC)->GetBinContent(bin);
    }
    return 0;
    
}

Bool_t HParticleEvtCharaBK::printCentralityClass(TString estimator)
{
    // legacy code
    // print all CentralityClasses in the Estimator

    if (estimator.CompareTo("TOFRPCtimecut")==0 || estimator.CompareTo("TOFRPC5")==0) 
        return printCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k5);
    else if (estimator.CompareTo("TOFRPCtimecutFOPI")==0 || estimator.CompareTo("TOFRPCFOPI")==0) 
        return printCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::kFOPI);
    else if (estimator.CompareTo("TOFRPCtimecut10")==0 || estimator.CompareTo("TOFRPC10")==0) 
        return printCentralityClass(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k10);
    else {
         Error("printCentralityClass()","Sorry. printCentralityClass() for %s not implemented yet!",estimator.Data());
         return kFALSE;
     }

    return kFALSE;
}

Bool_t HParticleEvtCharaBK::printCentralityClass(UInt_t centE, UInt_t centC)
{
    // print all CentralityClasses in the Estimator
    cout << endl;
    cout<<"---------------------------------------------------------------------------------------------"  << endl;
    cout<<"Centrality Classes for "<< getStringCentralityEstimator(centE).Data() << "  with " << getStringCentralityClass(centC).Data()  <<" bins" << endl;
    return printCentralityClass(getCentralityClassHist(centE, centC)); 
}

Bool_t HParticleEvtCharaBK::printCentralityClass(TH1* htemp)
{
    // print all CentralityClasses in the Estimator
    if(!htemp){
        Error("printCentralityClass()","Sorry. printCentralityClass() for ... not implemented yet!");
        return kFALSE;
    }
    else{
        cout<<"---------------------------------------------------------------------------------------------"  << endl;
        cout<< "# of Classes: "<< htemp->GetNbinsX()-2 << endl;
        printf(" Class   lowEdge -   upEdge  Centrality[%%]  BinWidth[%%]  real CentralityBin[%%]  BinCenter[%%]\n");
        Float_t pcent = 0;
        for(Int_t i = htemp->GetNbinsX(); i>0 ; --i) {
            printf("  %2.f :  %8.2f - %8.2f %13s %13.3f    %8.2f - %8.2f %13.2f\n",
                      htemp->GetBinContent(i),
                      (htemp->GetXaxis())->GetBinLowEdge(i),
                      (htemp->GetXaxis())->GetBinUpEdge(i),
                      (htemp->GetXaxis())->GetBinLabel(i),
                      htemp->GetBinError(i),
                      pcent,
                      pcent+htemp->GetBinError(i),
                      pcent+0.5*htemp->GetBinError(i) );
            pcent += htemp->GetBinError(i);
        }
        cout<<"---------------------------------------------------------------------------------------------"  << endl;
        return kTRUE;
    }
    return kFALSE;
}

TH2D* HParticleEvtCharaBK::makeEPresolution(TProfile2D *hist, Bool_t calcChi)
{
    if(!hist) return 0;
    TString sName=hist->GetName();   
    sName.Append("_res");
    
    Int_t nxBins = hist->GetNbinsX();
    Int_t nyBins = hist->GetNbinsY();
    
    TH2D* hresolution  = (TH2D*) new TH2D(sName.Data(), sName.Data(), nxBins, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(),
                                                                      nyBins, hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax());
    (hresolution->GetXaxis())->SetTitle(hist->GetXaxis()->GetTitle());
    (hresolution->GetYaxis())->SetTitle(hist->GetYaxis()->GetTitle());
    (hresolution->GetZaxis())->SetTitle(hist->GetZaxis()->GetTitle());

    Int_t bin,binx,biny;
    for (binx =0;binx<=nxBins+1;binx++) {
       for (biny =0;biny<=nyBins+1;biny++) {
           bin   = biny*(nxBins+2) + binx;
           
           Double_t dQaQb        = hist->GetBinContent(bin);
           Double_t dEntriesQaQb = hist->GetBinEntries(bin);
           if( dQaQb <= 0 || dEntriesQaQb < 1 ){
               hresolution->SetBinContent(bin, 0. );
               hresolution->SetBinError(  bin, 0. );
                continue;
           }
           Double_t dSpreadQaQb  = hist->GetBinError(bin);
           Double_t dV = TMath::Sqrt(dQaQb);
           printf("\nbin=%d binX=%d binY=%d QaQb = %f  +- %f    V = %f\n",bin, binx, biny, dQaQb, dSpreadQaQb , dV);
           Double_t dChi = FindXi(dV,1e-6);
           dV = ComputeResolution( TMath::Sqrt2()*dChi );
           printf("An estimate of the event plane resolution is: %f\n", dV );
           Double_t dVerr = 0.;
           if(dQaQb > 0.) dVerr = (1./(2.*pow(dQaQb,0.5)))*dSpreadQaQb;
           Double_t dChiErr = FindXi(dVerr,1e-6);
           printf("An estimate chi of the event plane is: %f +- %f\n", dChi, dChiErr );
           printf("R:(subevents) = %f +- %f\n",dV,dVerr);
           if(calcChi){
               hresolution->SetBinContent(binx, biny, dChi );
               hresolution->SetBinError(  binx, biny, dChiErr);
           }
           else{
               hresolution->SetBinContent(binx, biny, dV );
               hresolution->SetBinError(  binx, biny, dVerr );
           }
       }
    }
    TString name=hist->GetName();   
    hresolution->SetNameTitle(name.Data(),name.Data());
    hresolution->SetEntries(hist->GetEntries());
    return hresolution;

}
TH1* HParticleEvtCharaBK::makeEPresolution(TH3 *hist, Bool_t calcChi)
{
    if(!hist) return 0;
    TString sName=hist->GetName();   
    sName.Append("_res");
    
    TAxis* axis = hist->GetZaxis();
    axis->SetRange(1, axis->GetNbins());
    TH1* temp  = hist->Project3D("xy_1");

    axis->SetRange(axis->FindBin(TMath::PiOver2()), axis->GetNbins());
    TH1* temp2  = hist->Project3D("xy_2");    
    temp2->Divide(temp);

    Int_t nxBins = temp2->GetNbinsX();
    Int_t nyBins = temp2->GetNbinsY();
    
    TH2D* hresolution  = (TH2D*) new TH2D(sName.Data(), sName.Data(), nxBins, temp2->GetXaxis()->GetXmin(), temp2->GetXaxis()->GetXmax(),
                                                                      nyBins, temp2->GetYaxis()->GetXmin(), temp2->GetYaxis()->GetXmax());
    (hresolution->GetXaxis())->SetTitle(hist->GetXaxis()->GetTitle());
    (hresolution->GetYaxis())->SetTitle(hist->GetYaxis()->GetTitle());
    (hresolution->GetZaxis())->SetTitle(hist->GetZaxis()->GetTitle());
    
    Int_t bin,binx,biny;
    
    for (binx =0;binx<=nxBins+1;binx++) {
       for (biny =0;biny<=nyBins+1;biny++) {
           bin   = biny*(nxBins+2) + binx;
           
           Double_t ratio        = temp2->GetBinContent(binx,biny);
           printf("\nbin=%d binX=%d binY=%d  Ratio=%f\n",bin, binx, biny, ratio);
           if( ratio <= 0 ){
               continue;
           }
           Double_t chisq = -2.*TMath::Log(2.*ratio);
           Double_t dChi = sqrt(chisq);
           printf("An estimate chi of the event plane is: %f\n", dChi );
           Double_t dV = ComputeResolution( TMath::Sqrt2()*dChi );
           printf("An estimate of the event plane resolution is: %f\n", dV );
           Double_t dVerr = 0.;
           printf("R:(subevents) = %f +- %f\n",dV,dVerr);

           if(calcChi){
               hresolution->SetBinContent(binx, biny, dChi );
               //hresolution->SetBinError(  binx, biny, dChiErr);
           }
           else{

               hresolution->SetBinContent(binx, biny, dV );
               hresolution->SetBinError(  binx, biny, dVerr );
           }
       }
    }
    
    TString name=hist->GetName();   
    hresolution->SetNameTitle(name.Data(),name.Data());
    hresolution->SetEntries(hist->GetEntries());
    return hresolution;
}

//--------------------------------------------------------------------
Double_t HParticleEvtCharaBK::ModifiedBesselI(Int_t n, Double_t x) const
{
  // compute half-integer modified Bessel functions
  // order: n>0, for n>5/2, interpolation is used (improve this by using recurrence!!!)
  const Double_t FACTOR = 0.797884561;  //FIXME
  if (n<0) return 0;
  if (x<1e-7) return 0;

  if (n==0)       return FACTOR*sqrt(x)*TMath::SinH(x)/x;   // 1/2
  else if (n==1)  return FACTOR*sqrt(x)*( -TMath::SinH(x)/(x*x) + TMath::CosH(x)/x );   // 3/2
  else if (n==2)  return FACTOR*sqrt(x)*( (3./(x*x)+1.)*TMath::SinH(x)/x - 3.*TMath::CosH(x)/(x*x) );   // 5/2
  return 0.5*(TMath::BesselI(n,x)+TMath::BesselI(n+1,x));   // use average of integer-order Bessel
}

//--------------------------------------------------------------------
Double_t HParticleEvtCharaBK::ComputeResolution( Double_t x, Int_t n ) const
{
  // Computes resolution for Event Plane method
  if(x > 51.3) {
    printf("Warning: Estimation of total resolution might be WRONG. Please check!");
    return 0.99981;
  }
  if (n<1) return 0;
  Double_t n1 = (n-1)/2;      
  Double_t n2 = (n+1)/2;
  
  Double_t a = x*x/2;  //   in formula (8) of Ollitrault arXiv:nucl-ex/9711003v2
                       //   it is x*x/2 for  full EP resolution here sub-event resolution
  Double_t b = TMath::Exp(-a);
  if (n==1)        b *= TMath::BesselI0(a)+TMath::BesselI1(a);
  else if(n%2==1)  b *= TMath::BesselI(n1,a)+TMath::BesselI(n2,a);
  else             b *= ModifiedBesselI(n1, a) + ModifiedBesselI(n2, a);
  return TMath::Sqrt(TMath::Pi())/2*x*b; 
}


//--------------------------------------------------------------------
Double_t HParticleEvtCharaBK::FindXi( Double_t res, Double_t prec, Int_t n ) const
{
  // Computes x(res) for Event Plane method
  if(res > 0.99981) {
    printf("Warning: Resolution for subEvent is high. You reached the precision limit.");
    return 51.3;
  }
  int nSteps =0;
  Double_t xtmp=0, xmin=0, xmax=51.3, rtmp=0, delta=2*prec;
  while( delta > prec ) {
    xtmp = 0.5*(xmin+xmax);
    rtmp = ComputeResolution(xtmp, n);
    delta = TMath::Abs( res-rtmp );
    if(rtmp>res) xmax = xtmp;
    if(rtmp<res) xmin = xtmp;
    nSteps++;
  }
  return xtmp;
}
//--------------------------------------------------------------------
Bool_t HParticleEvtCharaBK::addEstimatorHist(TH1F *hist, Float_t fractionXsection, UInt_t centE, Int_t direction)
{
    if(!hist)return kFALSE;
    if(fractionXsection<0.1) return kFALSE;
    if(centE>=kNumCentralityEstimator) return kFALSE;
    
    TString temp;
    if(fPostFix.CompareTo("")==0) temp = Form("%s",getStringCentralityEstimator(centE).Data() );
    else                       temp = Form("%s_%s",getStringCentralityEstimator(centE).Data(), fPostFix.Data() );
    
    hist->SetName(temp.Data());
    hist->SetTitle(temp.Data());
    fEstimatorHist[centE]= (TH1F*) hist;
    fCentralityPercentileHist[centE] = (TH1F*) makePercentiles(hist, fractionXsection, direction);
    for (int centC = 0; centC < (int) kNumCentralityClass; ++centC){
        fCentralityHist[centE][centC] = (TH1F*) makeClasses(hist, fractionXsection, centC, direction);
        if(fCentralityHist[centE][centC]) printCentralityClass(centE, centC);
    }
    return kTRUE;
}

TH1F* HParticleEvtCharaBK::makePercentiles(TH1F *htemp, Float_t fractionXsection, Int_t direction)
{
  if(!htemp) return 0;  

  TH1F *hpercent  = (TH1F*) htemp->Clone("hpercent");
  TString name=htemp->GetName();   
  name.Append("_percentile");
  hpercent->SetNameTitle(name.Data(),name.Data());
  hpercent->Reset();
  Float_t totIntegral = fractionXsection / htemp->Integral(1,htemp->GetNbinsX()); 
  if(direction<0){
      for (int ibin=1; ibin<=htemp->GetNbinsX(); ibin++){
          hpercent->SetBinContent(ibin, totIntegral * htemp->Integral(ibin,htemp->GetNbinsX()) );
      } 
  }
  else if(direction>0){
      for (int ibin=1; ibin<=htemp->GetNbinsX(); ibin++){
          hpercent->SetBinContent(ibin, totIntegral * htemp->Integral(1,ibin) );
      } 
  }
  return hpercent;
}

TH1F* HParticleEvtCharaBK::makeClasses(TH1F *htemp, Float_t fractionXsection, UInt_t centC, Int_t direction)
{

    if(!htemp) return 0;
    if(fractionXsection<0.1) return 0;
    Int_t nClasses = getCentralityClassNbins(centC);
    Double_t integral = htemp->Integral();
    Double_t norm = integral / fractionXsection;
    Float_t* PercentileArray = getCentralityClassArray(centC);
    std::vector<Double_t> binEdge;
    std::vector<Double_t> xSection;
    std::vector<TString>  fLabels; 


    Int_t start = 1;
    Int_t stop = htemp->GetNbinsX();
    if(direction < 0){
        direction = -1;
        start = htemp->GetNbinsX()-1;
        stop = 0;
        binEdge.push_back(htemp->GetBinLowEdge(htemp->GetNbinsX()));
    }
    else{
        direction = 1;
        binEdge.push_back(htemp->GetBinLowEdge(1));
    }
    Int_t bin = start;
    Float_t lxs = 0;
    Float_t txs = 0;


    // Find edge and starting point at begining
    while(1){
        lxs += htemp->GetBinContent(bin);
        Double_t pxs = lxs/norm;
        if( pxs>0.001 ){
          if (direction>0)   binEdge.push_back(htemp->GetBinLowEdge(bin+1));
          else               binEdge.push_back(htemp->GetBinLowEdge(bin));
          xSection.push_back(pxs);
          bin += direction;
          break;
        }
        bin += direction;
    }

    Float_t totxs = 0;
    for (Int_t i = 0; i < nClasses; ++i) {
        lxs = 0 ;
        totxs += PercentileArray[i];
        while(1){
            lxs += htemp->GetBinContent(bin);
            txs += htemp->GetBinContent(bin);
            Double_t pxs = txs/norm;
            Double_t tdiff = (txs+htemp->GetBinContent(bin+direction))/norm;
            if( ( pxs>totxs )
                || (TMath::Abs(pxs-totxs)<=TMath::Abs(tdiff-totxs)) )
            {
              if (direction>0)  binEdge.push_back(htemp->GetBinLowEdge(bin+1));
              else              binEdge.push_back(htemp->GetBinLowEdge(bin));
              xSection.push_back(lxs/norm);
              bin += direction;
              break;
            }
            bin += direction;
            if ( (direction>0 && bin>=stop) || (direction<0 && bin<=stop ) || (txs>=integral) ) break;
        }
        if(totxs>fractionXsection || (direction>0 && bin>=stop) || (direction<0 && bin<=stop ) || (txs>=integral)) break;
    }

      fLabels.push_back("overflow");
      Double_t totXsection =0;
      for(std::vector<double>::size_type index = 1; index < xSection.size()-1; ++index){
          fLabels.push_back(Form("%02.0f-%02.0f",round(totXsection),round(totXsection+xSection[index])) );
          totXsection += xSection[index];
      }
      fLabels.push_back("underflow");
 
      std::reverse(fLabels.begin(),fLabels.end());
      std::reverse(binEdge.begin(), binEdge.end());
      std::reverse(xSection.begin(), xSection.end());

      Double_t xlowbins[binEdge.size()];
      std::copy(binEdge.begin(), binEdge.end(), xlowbins);

      TString name = htemp->GetTitle();
      name = Form("%s_%s_fixedCuts", name.Data(), getStringCentralityClass(centC).Data());
      TH1F *hfixedCuts = new TH1F(name.Data(), name.Data(), binEdge.size()-1, xlowbins); 

      for(std::vector<double>::size_type bin = 0; bin < fLabels.size(); ++bin){
          (hfixedCuts->GetXaxis())->SetBinLabel(bin+1,fLabels[bin]);
          hfixedCuts->SetBinContent(bin+1, fLabels.size()-bin-1);
          hfixedCuts->SetBinError(bin+1, xSection[bin]);
      }
  return hfixedCuts;
}



Int_t HParticleEvtCharaBK::getNbins(TString estimator)
{
    // legacy code
    // Number of Bins (CentralityClasses and Over- and Underflow) in the Estimator

    if (estimator.CompareTo("TOFRPCtimecut")==0 || estimator.CompareTo("TOFRPC5")==0) 
        return getNbins(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k5);
    else if (estimator.CompareTo("TOFRPCtimecutFOPI")==0 || estimator.CompareTo("TOFRPCFOPI")==0) 
        return getNbins(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::kFOPI);
    else if (estimator.CompareTo("TOFRPCtimecut10")==0 || estimator.CompareTo("TOFRPC10")==0) 
        return getNbins(HParticleEvtCharaBK::kTOFRPC, HParticleEvtCharaBK::k10);
    else { Error("getNbins()","Sorry. getNbins() for %s not implemented yet!",estimator.Data()); return 0;}
}
Int_t HParticleEvtCharaBK::getNbins(UInt_t centE, UInt_t centC)
{
    // Number of Bins (CentralityClasses and Over- and Underflow) in the Estimator
    TH1F *htemp = getCentralityClassHist(centE, centC); 
    if(htemp) return htemp->GetNbinsX();
    else return 0;
}

Int_t HParticleEvtCharaBK::getCentralityClassNbins(UInt_t centC)
{
    if(centC==k1040)      return 6;
    else if(centC==kFOPI) return 6;
    Float_t binSize = getCentralityClassBinSize(centC);
    if(binSize>0) return round((100./binSize));
    //    cout << "binSize: "<< binSize  <<"  nBins: " <<nBins << " " <<endl;
    return 0;
}

Float_t* HParticleEvtCharaBK::getCentralityClassArray(UInt_t centC)
{
    // 
    Float_t binSize = getCentralityClassBinSize(centC);
    Int_t nBins = getCentralityClassNbins(centC);
    if (nBins<1) return 0;
    
    Float_t* arr = new Float_t[nBins];
    if(binSize>0){
        Float_t xs = 0.;
        for(Int_t i =0;  i < nBins ; i++){
            xs+=binSize;
            if(xs < 100.){
                arr[i] = (Float_t) binSize;
            }
            else{
                arr[i] = 100. - (binSize*i);
            }
        }
    }
    else if(centC==k1040) {
        Double_t fxs[6] = {10.,30.,30, 10.,10.,10.};
        for(Int_t i =0;  i < nBins ; i++) arr[i] = fxs[i];
    }
    else if(centC==kFOPI) {
        Double_t fxs[6] = {6.3,14.7,9.9,10.,10.,10.};
        for(Int_t i =0;  i < nBins ; i++) arr[i] = fxs[i];
    }
    
    return arr;
}

Float_t* HParticleEvtCharaBK::getUpEdgeArray(UInt_t centE, UInt_t centC)
{
    // 
    TH1F *htemp = getCentralityClassHist(centE, centC);
    
    if(htemp){
        Int_t nBins = htemp->GetNbinsX();
        Float_t* arr = new Float_t[nBins];
        for(Int_t i =0;  i < nBins ; i++) {
            //printf(" %8.2f\n",(htemp->GetXaxis())->GetBinLowEdge(nBins-i));
            arr[i] = (Float_t) (htemp->GetXaxis())->GetBinLowEdge(nBins-i);
        }
         return arr;
    }
    else return 0;
}

Float_t* HParticleEvtCharaBK::getBinCenterArray(UInt_t centE, UInt_t centC)
{
    // 
    TH1F *htemp = getCentralityClassHist(centE, centC);
    
    if(htemp){
        Int_t nBins = htemp->GetNbinsX();
        Float_t* arr = new Float_t[nBins];
        Float_t pcent = 0;
        for(Int_t i =0;  i < nBins ; i++) {
                   pcent += 0.5*htemp->GetBinError(i);
            //       printf(" %8.2f\n",pcent);
            arr[nBins-i] = (Float_t) pcent;
        }
         return arr;
    }
    else return 0;
}

vector<TString> HParticleEvtCharaBK::getLabelArray(UInt_t centE, UInt_t centC)
{
    // 
    TH1F *htemp = getCentralityClassHist(centE, centC); 
    vector<TString> ar;
    if(htemp){
        Int_t nBins = htemp->GetNbinsX();
        for(Int_t i =1;  i < nBins-1 ; i++) {
            //printf(" %s\n",(htemp->GetXaxis())->GetBinLabel(nBins-i));
            TString prv = TString( (htemp->GetXaxis())->GetBinLabel(nBins-i) );
            ar.push_back(prv);
        }
    }
    return ar;
}

Int_t HParticleEvtCharaBK::getCentralityClass5(TString estimator)
{
    // Return centrality class 5% of total cross section
    return 1+(Int_t)getCentralityPercentile(estimator)/5.;
}
Int_t HParticleEvtCharaBK::getCentralityClass10(TString estimator)
{
    // Return centrality class 5% of total cross section
    return 1+(Int_t)getCentralityPercentile(estimator)/10.;
}

Float_t HParticleEvtCharaBK::getCentralityPercentile(TString estimator)
{
    // legacy code
    return getCentralityPercentile(getCentralityEstimatorFromString(estimator));
}


Float_t HParticleEvtCharaBK::getCentralityPercentile(UInt_t centE)
{

    if(centE>=kNumCentralityEstimator){
        return 101.;
    }
    else if(!getCentralityPercentileHist(centE)){
      //  Error("getCentralityPercentile()","No CentralityEstimator defined!");
        return 101.;
    }
    else {
        Int_t buffer = getCentralityEstimator(centE);
        Int_t bin = getCentralityPercentileHist(centE)->FindBin(buffer);
        return getCentralityPercentileHist(centE)->GetBinContent(bin);
    }

    return 101.;
}

Int_t HParticleEvtCharaBK::getMdcWiresOuterMod()
{
    HParticleEvtInfo *event_info = (HParticleEvtInfo*)fParticleEvtInfoCat->getObject(0);
    if(!event_info) {
        //Error("getCentralityClass()","No HParticleEvtInfo");
        return 0;
    }
    Int_t sum = 0;
    for(Int_t s=0;s<6;s++) sum += event_info->getMdcWiresMod(s,2)+event_info->getMdcWiresMod(s,3);  //outer mdc
    return sum;
}


Float_t HParticleEvtCharaBK::getSelectedParticleCandSecNorm(){
    
    if(!gLoop) return 0;
    
    HParticleEvtInfo *event_info = (HParticleEvtInfo*)fParticleEvtInfoCat->getObject(0);
    if(!event_info) {
       // Error("getCentralityClass()","No HParticleEvtInfo");
        return 0;
    }
    
    Float_t SelectedParticleCandSecCorr=0;
    Int_t nGoodSectors=0;
    for (Int_t s = 0; s < 6; ++s ){
        if(gLoop->goodSector(s)){
            SelectedParticleCandSecCorr += event_info->getSelectedParticleCandMult(s)/event_info->getMeanMult(s);
            nGoodSectors++;
        }
    }
    if(nGoodSectors>0) SelectedParticleCandSecCorr /= nGoodSectors;
    return SelectedParticleCandSecCorr;
}

Float_t HParticleEvtCharaBK::getSelectedParticleCandCorrPerWire(){
    
    HParticleEvtInfo *event_info = (HParticleEvtInfo*)fParticleEvtInfoCat->getObject(0);
    if(!event_info) {
        //Error("getCentralityClass()","No HParticleEvtInfo");
        return 0;
    }
    // Correction function with parameters provided by O.&V.Pechenova
    Float_t w0 = -0.660682;
    Float_t w1 = -0.0652827;
    Float_t w2 = -0.660682;
    Float_t w3 = -0.0652827;
    Float_t w4 = -0.655548;
    Float_t w5 = -0.00547515;
    
    Int_t   nRawWires = 0;
    Int_t   nSelMult = 0;
    Float_t WirePerTracks = 0;
    Float_t SelectedParticleCandCorrPerWire = 0;
    Float_t EffWirePerTrack = 0;
    for(Int_t s=0;s<6;s++) {
        nSelMult  = event_info->getSelectedParticleCandMult(s);
        nRawWires = event_info->getMdcWiresMod(s,2)+event_info->getMdcWiresMod(s,3);  //outer mdc
        if(nSelMult>0){
            WirePerTracks = Float_t(nRawWires)/nSelMult;
            EffWirePerTrack = 1/(exp(w0+w1*WirePerTracks) + exp(w2+w3*WirePerTracks) + exp(w4+w5*WirePerTracks));
            SelectedParticleCandCorrPerWire += nSelMult*EffWirePerTrack;
        }
    }
    return SelectedParticleCandCorrPerWire;
}

Float_t HParticleEvtCharaBK::getFWSumChargeSpec(){

    Float_t MultFWSumChargeSpec=0;

    //HWallHit
    if(!fCatWallHit) return -1;
    HWallHitSim *wall = 0;
    for(Int_t i = 0; i < fCatWallHit->getEntries(); i ++ ){
        wall  = (HWallHitSim*) HCategoryManager::getObject(wall,fCatWallHit,i);
        if(!wall) continue;
        if(!PassesCutsFW(wall)) continue;
        if(!isSimulation) MultFWSumChargeSpec += wall->getCharge();
        else              MultFWSumChargeSpec += 93.*pow(wall->getCharge(),0.46-0.006*sqrt(wall->getCharge()));  // parametrization from R.Holzmann
     }

    return MultFWSumChargeSpec;
}

Float_t HParticleEvtCharaBK::getFWSumZ(Int_t minZ, Int_t maxZ, UInt_t SubEvent){

    if(isNewEvent()){
        fillHitArray();
    }
    Float_t Z = 0.;
    Float_t w = 1.;

    Int_t nHits = arrayOfHits.size();
    if(nHits<1) return 0;  // at least one hit
    for (Int_t i=0; i<nHits; i++){
        w  = arrayOfHits[i]->Weight1();
        //w2 *= arrayOfHits[i]->Weight2();
        if(w>maxZ) w=maxZ;
        if(w<minZ) continue;

        if(SubEvent == 0) {
            Z += w;
        }
        else if(SubEvent == 1 && arrayOfHits[i]->SubEvent() == 1) {
            Z += w;
        }
        else if(SubEvent == 2 && arrayOfHits[i]->SubEvent() == 2) {
            Z += w;
        }
    }
    return Z;
}

Float_t HParticleEvtCharaBK::getEt(){

    Float_t Et=0;
    if(!fParticleCandCat) return -1;
    //if(!fParticleCandCat) {Error("calcFWSumChargeSpec()","No calcEt");return kFALSE;}
    HParticleCand* particle_cand = 0;

    for(Int_t i = 0; i < fParticleCandCat->getEntries(); i ++ ){
        particle_cand  = (HParticleCand*) HCategoryManager::getObject(particle_cand,fParticleCandCat,i);
        if(!particle_cand->isFlagBit(Particle::kIsUsed)) continue;
        if(particle_cand->isFakeRejected() )            continue;
        if(particle_cand->getSystemUsed() == -1)        continue;
        if(!particle_cand->isFlagAND(4,                               
                Particle::kIsAcceptedHitInnerMDC,
                Particle::kIsAcceptedHitOuterMDC,
                Particle::kIsAcceptedHitMETA,
                Particle::kIsAcceptedRK)  )       continue;
//        if(!particle_cand->getPID()==14)       continue;
          if(particle_cand->getPID()==-1)       continue;
//        if(!particle_cand->getPID()==8)       continue;
//        if(!particle_cand->getPID()==9)       continue;
//        if( particle_cand->getMomentumOrg() == particle_cand->getMomentum() ) continue;
        if( TMath::Abs(particle_cand->Rapidity()-0.74)>0.5 && particle_cand->Pt()>300 ){  //
            Et              += particle_cand->Et();
        }
    //-----------------------------------------------------------------
    } // end particle loop

    return Et;

}

Float_t HParticleEvtCharaBK::getRatioEtEz(){
    
    Float_t RatioEtEz=0;
    Float_t Et=0;
    Float_t Ez=0;
    if(!fParticleCandCat) return -1;
    //if(!fParticleCandCat) {Error("calcRatioEtEz()","No fParticleCandCat");return kFALSE;}
    HParticleCand* particle_cand = 0;

    for(Int_t i = 0; i < fParticleCandCat->getEntries(); i ++ ){
        particle_cand  = (HParticleCand*) HCategoryManager::getObject(particle_cand,fParticleCandCat,i);
        if(!particle_cand->isFlagBit(Particle::kIsUsed)) continue;
        if(particle_cand->isFakeRejected() )            continue;
        if(particle_cand->getSystemUsed() == -1)        continue;
        if(!particle_cand->isFlagAND(4,                               
                Particle::kIsAcceptedHitInnerMDC,
                Particle::kIsAcceptedHitOuterMDC,
                Particle::kIsAcceptedHitMETA,
                Particle::kIsAcceptedRK)  )       continue;
        if(particle_cand->getPID()==-1)       continue;
        //if(!particle_cand->getPID()==14)       continue;
        //if(!particle_cand->getPID()==8)       continue;
        //if(!particle_cand->getPID()==9)       continue;
        //if( particle_cand->getMomentumOrg() == particle_cand->getMomentum() ) continue;
//      if( TMath::Abs(particle_cand->Rapidity()-0.74)>0.5 && particle_cand->Pt()>300 ){  //
            Ez              += particle_cand->E()*particle_cand->CosTheta();   // pt2 == 0 ? 0 : E()*E() * pt2/(pt2+Z()*Z());
            Et              += particle_cand->Et();
            //cout<<"  Sum Et:  "<<Et<<"  Et:  "<<tlv1.Et()<<endl;
//      }
    //-----------------------------------------------------------------
    } // end particle loop
    if(Ez){
        RatioEtEz = Et/Ez;
    }
    else {
        RatioEtEz=0;
    }
    return RatioEtEz;


}


Float_t HParticleEvtCharaBK::getDirectivity(){

    Float_t Directivity=0;
    TLorentzVector QVector = TLorentzVector();
    Int_t QVectorT = 0;
        if(!fParticleCandCat) return -1;
    //if(!fParticleCandCat) {Error("calcDirectivity()","No fParticleCandCat");return kFALSE;}
    HParticleCand* particle_cand = 0;

    for(Int_t i = 0; i < fParticleCandCat->getEntries(); i ++ ){
        particle_cand  = (HParticleCand*) HCategoryManager::getObject(particle_cand,fParticleCandCat,i);
        if(!particle_cand->isFlagBit(Particle::kIsUsed)) continue;
        if(particle_cand->isFakeRejected() )            continue;
        if(particle_cand->getSystemUsed() == -1)        continue;
        if(!particle_cand->isFlagAND(4,                               
                Particle::kIsAcceptedHitInnerMDC,
                Particle::kIsAcceptedHitOuterMDC,
                Particle::kIsAcceptedHitMETA,
                Particle::kIsAcceptedRK)  )       continue;
        if(particle_cand->getPID()==-1)           continue;
        if( particle_cand->Rapidity() > 0.74 ){  //
            QVector  += (TLorentzVector)* particle_cand;
            QVectorT +=  particle_cand->Pt();
        }
        else {  //
            QVector  -= (TLorentzVector)* particle_cand;
            QVectorT += particle_cand->Pt();
        }



    //-----------------------------------------------------------------
    } // end particle loop

    if(QVectorT>0)Directivity = QVector.Pt() / QVectorT;      //  = |Sum(vector_p_{t,i})| / Sum(|p_{t,i}|)   for 0>Y_cm
    else Directivity = 0;

    return Directivity;

}

