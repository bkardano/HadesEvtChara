#ifndef  HPARTICLEEVTCHARABK_H
#define  HPARTICLEEVTCHARABK_H

#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TNamed.h"
#include "TH1F.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "vector"
#include "TNtuple.h"
#include "TVector2.h"
#include "TRandom2.h"

#include "hades.h"
#include "htool.h"
#include "hphysicsconstants.h"
#include "hrootsource.h"
#include "hiterator.h"
#include "hloop.h"
#include "htime.h"
#include "hruntimedb.h"
#include "htaskset.h"

#include "haddef.h"
#include "heventheader.h"
#include "hreconstructor.h"
#include "hparticledef.h"
#include "hparticlestructs.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"

#include "hcategory.h"
#include "hcategorymanager.h"
#include "hparticleevtinfo.h"

#include "walldef.h"
#include "hwallhitsim.h"
//#include "hwallraw.h"
#include "henergylosscorrpar.h"

using namespace std;
class HParticleCand;
class HCategory;

class HParticleEvtCharaBK : public HReconstructor {

public:


    // -------------------------------------------------------------------------
    //  flags for HParticleEvtCharaBK, CentralityEstimator
    enum eCentralityEstimator { 
        kTOFRPC                          =  0,
        kTOF                                 ,
        kRPC                                 ,
        kTOFRPCtot                           ,
        kTOFtot                              ,
        kRPCtot                              ,
        kSelectedParticleCand                ,
        kSelectedParticleCandCorr            ,
        kSelectedParticleCandNorm            ,
        kSelectedParticleCandSecNorm         ,
        kSelectedParticleCandCorrPerWire     ,
        kPrimaryParticleCand                 ,
        kMdcWires                            ,
        kMdcWiresOuterMod                    ,
        kFWSumChargeSpec                     ,
        kFWSumChargeZ                        ,
        kDirectivity                         ,
        kRatioEtEz                           ,
        kEt                                  ,
//        kTOFRPCvsFWSumChargeZ                ,
        kNumCentralityEstimator              
    };

    Int_t getCentralityEstimatorBinSize(UInt_t e)
    {
        if(e     ==kTOFRPC)                          return    500;
        else if(e==kTOF)                             return    500;
        else if(e==kRPC)                             return    500;
        else if(e==kTOFRPCtot)                       return    500;
        else if(e==kTOFtot)                          return    500;
        else if(e==kRPCtot)                          return    500;
        else if(e==kSelectedParticleCand)            return    500;
        else if(e==kSelectedParticleCandCorr)        return    500;
        else if(e==kSelectedParticleCandSecNorm)     return    500;
        else if(e==kSelectedParticleCandNorm)        return    500;
        else if(e==kSelectedParticleCandCorrPerWire) return    500;
        else if(e==kPrimaryParticleCand)             return    500;
        else if(e==kMdcWires)                        return   1000;
        else if(e==kMdcWiresOuterMod)                return   1000;
        else if(e==kFWSumChargeSpec)                 return   1000;
        else if(e==kFWSumChargeZ)                    return    500;
        else if(e==kDirectivity)                     return   1000;
        else if(e==kRatioEtEz)                       return   1000;
        else if(e==kEt)                              return   1000;
        else return 0;
    }

    Int_t getCentralityEstimatorBinMax(UInt_t e)
    {
        if(e     ==kTOFRPC)                          return     500; 
        else if(e==kTOF)                             return     500;
        else if(e==kRPC)                             return     500;
        else if(e==kTOFRPCtot)                       return     500;
        else if(e==kTOFtot)                          return     500;
        else if(e==kRPCtot)                          return     500;
        else if(e==kSelectedParticleCand)            return     500;
        else if(e==kSelectedParticleCandCorr)        return     500;
        else if(e==kSelectedParticleCandSecNorm)     return      50;
        else if(e==kSelectedParticleCandNorm)        return      50;
        else if(e==kSelectedParticleCandCorrPerWire) return     500;
        else if(e==kPrimaryParticleCand)             return     500;
        else if(e==kMdcWires)                        return   10000; 
        else if(e==kMdcWiresOuterMod)                return   10000; 
        else if(e==kFWSumChargeSpec)                 return   10000; 
        else if(e==kFWSumChargeZ)                    return     500; 
        else if(e==kDirectivity)                     return       3; 
        else if(e==kRatioEtEz)                       return       3; 
        else if(e==kEt)                              return   30000; 
        else return 0;
    }

    TString getStringCentralityEstimator(UInt_t e)
    {
        if(e     ==kTOFRPC)                         return "TOFRPC";
        else if(e==kTOF)                            return "TOF";
        else if(e==kRPC)                            return "RPC";
        else if(e==kTOFRPCtot)                      return "TOFRPCtot";
        else if(e==kTOFtot)                         return "TOFtot";
        else if(e==kRPCtot)                         return "RPCtot";
        else if(e==kSelectedParticleCand)           return "SelectedParticleCand";
        else if(e==kSelectedParticleCandCorr)       return "SelectedParticleCandCorr";
        else if(e==kSelectedParticleCandSecNorm)    return "SelectedParticleCandSecNorm";
        else if(e==kSelectedParticleCandNorm)       return "SelectedParticleCandNorm";
        else if(e==kSelectedParticleCandCorrPerWire)return "SelectedParticleCandCorrPerWire";
        else if(e==kPrimaryParticleCand)            return "PrimaryParticleCand";
        else if(e==kMdcWires)                       return "MdcWires";
        else if(e==kMdcWiresOuterMod)               return "MdcWiresOuterMod";
        else if(e==kFWSumChargeSpec)                return "FWSumChargeSpec";
        else if(e==kFWSumChargeZ)                   return "FWSumChargeZ";
        else if(e==kDirectivity)                    return "Directivity";
        else if(e==kRatioEtEz)                      return "RatioEtEz";
        else if(e==kEt)                             return "Et";
        else return "";
    }

    UInt_t getCentralityEstimatorFromString(TString s)
    {
             if(s.CompareTo("TOFRPC")==0 || s.CompareTo("TOFRPCtimecut")==0)    return kTOFRPC;
        else if(s.CompareTo("TOF")==0 || s.CompareTo("TOFtimecut")==0)          return kTOF;
        else if(s.CompareTo("RPC")==0 || s.CompareTo("RPCtimecut")==0)          return kRPC;
        else if(s.CompareTo("TOFRPCtot")==0)                                    return kTOFRPCtot;
        else if(s.CompareTo("TOFtot")==0)                                       return kTOFtot;
        else if(s.CompareTo("RPCtot")==0)                                       return kRPCtot;
        else if(s.CompareTo("SelectedParticleCand")==0)                         return kSelectedParticleCand;
        else if(s.CompareTo("SelectedParticleCandCorr")==0)                     return kSelectedParticleCandCorr;
        else if(s.CompareTo("SelectedParticleCandSecNorm")==0)                  return kSelectedParticleCandSecNorm;
        else if(s.CompareTo("SelectedParticleCandNorm")==0)                     return kSelectedParticleCandNorm;
        else if(s.CompareTo("SelectedParticleCandCorrPerWire")==0)              return kSelectedParticleCandCorrPerWire;
        else if(s.CompareTo("PrimaryParticleCand")==0)                          return kPrimaryParticleCand;
        else if(s.CompareTo("MdcWires")==0)                                     return kMdcWires;
        else if(s.CompareTo("MdcWiresOuterMod")==0)                             return kMdcWiresOuterMod;
        else if(s.CompareTo("FWSumChargeSpec")==0)                              return kFWSumChargeSpec;
        else if(s.CompareTo("FWSumChargeZ")==0)                                 return kFWSumChargeZ;
        else if(s.CompareTo("Directivity")==0)                                  return kDirectivity;
        else if(s.CompareTo("RatioEtEz")==0)                                    return kRatioEtEz;
        else if(s.CompareTo("Et")==0)                                           return kEt;
//        else if(s.CompareTo("TOFRPCvsFWSumChargeZ")==0)                                    return kTOFRPCvsFWSumChargeZ;
        else return kNumCentralityEstimator;
    }


    //  CentralityClass
    enum eCentralityClass { 
        k10                 = 0,        // (default) 10% classes
        k5                     ,        // 5% classes
        k2                     ,        // 2% classes
        k1                     ,        // 1% classes
        k13                    ,        // 13% classes 
        k15                    ,        // 15% classes
        k20                    ,        // 20% classes
        k1040                  ,        // 10-40% class
        kFOPI                  ,        // FOPI centrality classes
        kNumCentralityClass    
    };

    TString getStringCentralityClass(UInt_t e)
    {
             if(e==k10)                  return "10pc";
        else if(e==k5)                   return "5pc";
        else if(e==k2)                   return "2pc";
        else if(e==k1)                   return "1pc";
        else if(e==k13)                  return "13pc";
        else if(e==k15)                  return "15pc";
        else if(e==k20)                  return "20pc";
        else if(e==k1040)                return "1040class";
        else if(e==kFOPI)                return "FOPI";
        else return "";
    }

    Float_t getCentralityClassBinSize(UInt_t e)
    {
             if(e==k10)                  return 10. ;
        else if(e==k5)                   return  5. ;
        else if(e==k2)                   return  2. ;
        else if(e==k1)                   return  1. ;
        else if(e==k13)                  return 13. ;
        else if(e==k15)                  return 15. ;
        else if(e==k20)                  return 20. ;
        else if(e==k1040)                return  0. ;
        else if(e==kFOPI)                return  0. ;
        else return 0;
    }

    enum eEventPlaneCorrection {         // additive
        kNoCorrection           = 0x00,  //     
        kSmear                  = 0x01,  //    1
        kShiftFW                = 0x02,  //    2 
        kScaleFW                = 0x04,  //    2 
        kReCenter               = 0x08,  //    4
        kWidthEqualization      = 0x10,  //    8
        kFlatten                = 0x20,  //   16
        kWeightCharge           = 0x40,  //   32
        kWeightTheta            = 0x80, 
         

       kScalingFW = kScaleFW,
       kReCentering = kReCenter, 
       kScaling = kWidthEqualization,
       kFlattening = kFlatten,
       kRotation = kFlatten,
        
        kDefault                = kShiftFW|kWeightCharge|kRotation
    };

    TString getStringEventPlaneCorrection(UInt_t e)
    {
        TString s="";
        if(e==kNoCorrection)    return "NoCorrection";

        if(isFlagSet(kSmear,e)          )  s +="Smear";
        if(isFlagSet(kShiftFW,e)        )  s +="ShifFW";
        if(isFlagSet(kScaleFW,e)        )  s +="ScalFW";
        if(isFlagSet(kReCenter,e)       )  s +="Recent";
        if(isFlagSet(kScaling,e)        )  s +="Scale";
        if(isFlagSet(kFlattening,e)     )  s +="Flatt";
        if(isFlagSet(kWeightCharge,e)   )  s +="WCharge";
        if(isFlagSet(kWeightTheta,e)    )  s +="WTheta";
        return s;
    }

    //  EventPlaneCorrection Parameter
    enum eEventPlaneParameter {     
        kFWx                    = 0,
        kFWy                    ,
        kQx                     ,
        kQy                     ,
        kQxWCharge              ,
        kQyWCharge              ,
        kQx2                    ,
        kQy2                    ,
        kQx2WCharge             ,
        kQy2WCharge             ,
        kFourierC1              ,
        kFourierC2              ,
        kFourierC3              ,
        kFourierC4              ,
        kFourierC5              ,
        kFourierC6              ,
        kFourierC7              ,
        kFourierC8              , 
        kFourierS1              ,
        kFourierS2              ,
        kFourierS3              ,
        kFourierS4              ,
        kFourierS5              ,
        kFourierS6              ,
        kFourierS7              ,
        kFourierS8              ,
        kResolution             ,
        kResolutionWCharge      ,
        kChi                    ,
        kNumEventPlaneParameter 
    };
    //TString StringEventPlaneParameter[] = {
    //    "FWx"  
    //   ,"FWy"  
    //   ,"Qx"  
    //   ,"Qy"  
    //   ,"FWxWCharge"  
    //   ,"FWyWCharge"  
    //   ,"QxWCharge"  
    //   ,"QyWCharge"  
    //   ,"FourierC1"  
    //   ,"FourierC2"  
    //   ,"FourierC3"  
    //   ,"FourierC4"  
    //   ,"FourierS1"  
    //   ,"FourierS2"  
    //   ,"FourierS3"  
    //   ,"FourierS4"  
    //   ,"Resolution"  
    //   ,"ResolutionWCharge"  
    //   ,"Chi" 
    //};
     
    TString getStringEventPlaneParameter(UInt_t e)
    {
           if(e==kFWx)                  return "FWx"  ;
      else if(e==kFWy              )    return "FWy"  ;
      else if(e==kQx               )    return "Qx"  ;
      else if(e==kQy               )    return "Qy"  ;
      else if(e==kQxWCharge        )    return "QxWCharge"  ;
      else if(e==kQyWCharge        )    return "QyWCharge"  ;
      else if(e==kQx2               )    return "Qx2"  ;
      else if(e==kQy2               )    return "Qy2"  ;
      else if(e==kQx2WCharge        )    return "Qx2WCharge"  ;
      else if(e==kQy2WCharge        )    return "Qy2WCharge"  ;

      else if(e==kFourierC1        )    return "FourierC1"  ;
      else if(e==kFourierC2        )    return "FourierC2"  ;
      else if(e==kFourierC3        )    return "FourierC3"  ;
      else if(e==kFourierC4        )    return "FourierC4"  ;
      else if(e==kFourierC5        )    return "FourierC5"  ;
      else if(e==kFourierC6        )    return "FourierC6"  ;
      else if(e==kFourierC7        )    return "FourierC7"  ;
      else if(e==kFourierC8        )    return "FourierC8"  ;
      else if(e==kFourierS1        )    return "FourierS1"  ;
      else if(e==kFourierS2        )    return "FourierS2"  ;
      else if(e==kFourierS3        )    return "FourierS3"  ;
      else if(e==kFourierS4        )    return "FourierS4"  ;
      else if(e==kFourierS5        )    return "FourierS5"  ;
      else if(e==kFourierS6        )    return "FourierS6"  ;
      else if(e==kFourierS7        )    return "FourierS7"  ;
      else if(e==kFourierS8        )    return "FourierS8"  ;
      else if(e==kResolution       )    return "Resolution"  ;
      else if(e==kResolutionWCharge)    return "ResolutionWCharge";
      else if(e==kChi              )    return "Chi" ;
   //     if(e<kNumEventPlaneParameter) return StringEventPlaneParameter[e];
        else return "";
    }

    UInt_t getEventPlaneParameterFromString(TString s)
    {
         if(s.CompareTo("FWx")==0)          return kFWx;
    else if(s.CompareTo("FWy")==0)          return kFWy;
    else if(s.CompareTo("Qx")==0)           return kQx;
    else if(s.CompareTo("Qy")==0)           return kQy;
    else if(s.CompareTo("QxWCharge")==0)    return kQxWCharge;
    else if(s.CompareTo("QyWCharge")==0)    return kQyWCharge;
    else if(s.CompareTo("Qx2")==0)          return kQx2;
    else if(s.CompareTo("Qy2")==0)          return kQy2;
    else if(s.CompareTo("Qx2WCharge")==0)   return kQx2WCharge;
    else if(s.CompareTo("Qy2WCharge")==0)   return kQy2WCharge;
    else if(s.CompareTo("FourierC1")==0)    return kFourierC1;
    else if(s.CompareTo("FourierC2")==0)    return kFourierC2;
    else if(s.CompareTo("FourierC3")==0)    return kFourierC3;
    else if(s.CompareTo("FourierC4")==0)    return kFourierC4;
    else if(s.CompareTo("FourierC5")==0)    return kFourierC5;
    else if(s.CompareTo("FourierC6")==0)    return kFourierC6;
    else if(s.CompareTo("FourierC7")==0)    return kFourierC7;
    else if(s.CompareTo("FourierC8")==0)    return kFourierC8;
    else if(s.CompareTo("FourierS1")==0)    return kFourierS1;
    else if(s.CompareTo("FourierS2")==0)    return kFourierS2;
    else if(s.CompareTo("FourierS3")==0)    return kFourierS3;
    else if(s.CompareTo("FourierS4")==0)    return kFourierS4;
    else if(s.CompareTo("FourierS5")==0)    return kFourierS5;
    else if(s.CompareTo("FourierS6")==0)    return kFourierS6;
    else if(s.CompareTo("FourierS7")==0)    return kFourierS7;
    else if(s.CompareTo("FourierS8")==0)    return kFourierS8;
    else if(s.CompareTo("Resolution")==0)   return kResolution;
    else if(s.CompareTo("ResolutionWCharge")==0)  return kResolutionWCharge;
    else if(s.CompareTo("Chi")==0)          return kChi;
    else return kNumEventPlaneParameter;
    }



  // static const char* getStringEventPlaneParameter(UInt_t e){
  //     if(e<kNumEventPlaneParameter) return StringEventPlaneParameter[e];
  //     else return "";
  // }
  // //{ return (e<kNumEventPlaneParameter)?StringEventPlaneParameter[e]:""; }
  // 
  // UInt_t getEventPlaneParameterFromString(const char* valname)
  // {
  //   TString name(valname);
  //   for(UInt_t i=0; i<kNumEventPlaneParameter; i++) {
  //     if(!name.CompareTo(StringEventPlaneParameter[i])) return i;
  //   }
  //   return kNumEventPlaneParameter; // alias none
  // }

    // -------------------------------------------------------------------------
    // helper class
    class SimpleQVector {
        public:
            Int_t   fSubEvent;
            Float_t fPhi;
            Float_t fPhiOrg;
            Float_t fWeight1;
            Float_t fWeight2;
            UInt_t  fFlag;
            
            SimpleQVector(Int_t se=0, Float_t phi=-1., Float_t w1=1., Float_t w2=1., UInt_t f=kNoCorrection){
                fSubEvent = se;
                fPhi      = phi;
                fWeight1  = w1;
                fWeight2  = w2;
                fFlag     = f;
            }
            
            void SetSubEvent(Int_t se)                { fSubEvent = se;}
            void SetFlag(UInt_t flag)                  { fFlag = flag;}
            void Set(Float_t x, Float_t y, Float_t w1=1., Float_t w2=1.)
                                              { fPhi = TVector2::Phi_0_2pi( TMath::Pi()+TMath::ATan2(-y,-x) );fWeight1 = w1;fWeight2 = w2;}
            void SetOrg(Float_t x, Float_t y) { fPhiOrg = TVector2::Phi_0_2pi( TMath::Pi()+TMath::ATan2(-y,-x) );}

//            void Set(Float_t x, Float_t y)            { fPhi = TVector2::Phi_0_2pi( TMath::Pi()+TMath::ATan2(-y,-x) );}
//            void Set(Float_t x, Float_t y, Float_t w) { fPhi = TVector2::Phi_0_2pi( TMath::Pi()+TMath::ATan2(-y,-x) );fWeight1 = w;}

            void Rotate(Float_t psi) { fPhi = TVector2::Phi_0_2pi( fPhi+psi );}

            
            UInt_t  Flag()     {return fFlag;}   
            Int_t   SubEvent() {return fSubEvent;}   
            Float_t Weight1()  {return fWeight1;}   
            Float_t Weight2()  {return fWeight2;}   

            Float_t Phi()     {return fPhi;}   
            Float_t PhiOrg()  {return fPhiOrg;}   

            Float_t X()   {return TMath::Cos(fPhi);}   
            Float_t Y()   {return TMath::Sin(fPhi);}
            Float_t XOrg()   {return TMath::Cos(fPhiOrg);}   
            Float_t YOrg()   {return TMath::Sin(fPhiOrg);}
            
            void print(){
                //print stuff
                printf("Phi: %.5f, PhiOrg: %.5f, Weight1: %.3f, Weight2: %.3f, SubEvent: %d \n",
                           fPhi,fPhiOrg,fWeight1,fWeight2,fSubEvent);
            }
        
            void Clear()
            {
                fSubEvent  = 0;
                fPhi       = -1.;
                fPhiOrg    = -1.;
                fWeight1   = 1.;
                fWeight2   = 1.;
                fFlag      = kNoCorrection;
            }
    };

#define MAXFWCELLS 302

    //  FWCutValues
    enum eFWCut {     
        kBetaCuts            = 0,
        kTimeCuts            = 1,
        kChargeCuts          = 2,
        kNumFWCutValues        = 3
    };
    TString getStringFWCutValues(UInt_t e)
    {
             if(e==kBetaCuts    )   return "BetaCuts"  ;
        else if(e==kTimeCuts    )   return "TimeCuts"  ;
        else if(e==kChargeCuts  )   return "ChargeCuts"  ;
        else return "";
    }
    

    HParticleEvtCharaBK(const Text_t* name="HParticleEvtCharaBK",const Text_t* title="HParticleEvtCharaBK"); /// constructor
    ~HParticleEvtCharaBK();  /// destructor
    Float_t getVersion() {return 0.7;}  // get parameterfile version, to avoid loading old format 
    Int_t   getNumFWCells() {return MAXFWCELLS;}
    Int_t   execute();
    Bool_t  init();
    Bool_t  reinit();
    Bool_t  finalize(){ return kTRUE;}

    //setter:
    void    setPostFix(TString postFix)              {fPostFix = postFix;}
    Bool_t  setParameterFile(TString ParameterFile);
    void    setUseFWCut(UInt_t eFWCut, Bool_t b=kTRUE)  {useFWCut[eFWCut] = b;}
    void    setExcludeNoisyFWcells(Bool_t b=kTRUE)      {excludeNoisyFWcells = b;}

    Int_t   loadDayOfYear();
    Bool_t  loadParameterFile();
    Bool_t  saveParameterFile();
    Bool_t  loadCentralityEstimatorHist();
    Bool_t  addEstimatorHist(TH1F* hist, Float_t fractionXsection=100., UInt_t centE=kTOFRPC, Int_t direction=-1);
    Bool_t  saveCentralityEstimatorHist();

    Bool_t  loadEventPlaneCorrectionHist();
    Bool_t  addEventPlaneCorrectionHist(TProfile2D *hist, UInt_t epParam=kQx);
    Bool_t  saveEventPlaneCorrectionHist();

    Bool_t  loadFWCutValuesHist();
    Bool_t  addFWCutValuesHist(TH1 *hist, Int_t cell, UInt_t eFWCut);
    Bool_t  saveFWCutValuesHist();

    void    setReferenceMean(Float_t referenceMean)       {fReferenceMeanSelTrack = referenceMean;}

    //getter:
    Bool_t   isNewEvent();
    Float_t  getEventWeight();
    
    Int_t    GetFWmoduleSize(HWallHitSim* wall);
    Bool_t   PassesCutsFW(HWallHitSim* wall);
    Bool_t   PassesCutsFW(HWallHitSim* wall, UInt_t eFWCut);

    TH1F*    getEventPlaneCorrectionHist(UInt_t flag) const;
    Float_t  getCorrection(UInt_t flag=kQx);
    Float_t  getCorrectionError(UInt_t flag=kQx);
    Int_t    getFWCharge(HWallHitSim* wall);
    Float_t  getSmearValue(Int_t size=1);
    Float_t  getThetaWeight(HWallHitSim* wall, Float_t min=3.5, Float_t max=8.);

    vector<Int_t> getFWhits();
    void     fillHitArray();
    void     printHitArray();
    void     printQVectors();

    Bool_t    fillQVectors(UInt_t EPcorr=kNoCorrection, UInt_t nHarmonic=1);  
//    TVector2* getQVector(UInt_t EPcorr=kNoCorrection, UInt_t SubEvent=0, UInt_t nHarmonic=1);  

    Float_t  getCorrectionPhi(Float_t phi);
    Float_t  getEventPlane(UInt_t EPcorr=kNoCorrection, UInt_t SubEvent=0, UInt_t nHarmonic=1);  
//    Float_t  getQVectorMagnitude(UInt_t EPcorr=kNoCorrection, UInt_t SubEvent=0, UInt_t nHarmonic=1);
    Float_t  getEventPlaneWeight(UInt_t EPcorr=kNoCorrection, UInt_t SubEvent=0, UInt_t nHarmonic=1);
    Float_t  getEventPlaneParameter(UInt_t flag=kQx, Bool_t corr=kFALSE);
    
    
    Int_t   getCentralityClass(TString estimator);                       // Lagacy function should not be used
    Int_t   getCentralityClass(UInt_t centE=kTOFRPC, UInt_t centC=k10);  // Return centrality classes with adjusted fixed Cuts
    Float_t getCentralityEstimator(UInt_t centE=kTOFRPC);
    TH1F*   getCentralityClassHist(UInt_t centE=kTOFRPC, UInt_t centC=k10) const;

    Bool_t  printCentralityClass(TString estimator);
    Bool_t  printCentralityClass(UInt_t centE=kTOFRPC, UInt_t centC=k10);
    Bool_t  printCentralityClass(TH1* htemp);


    Float_t getCentralityPercentile(TString estimator);          // Return centrality percentile - Lagacy function should not be used
    Float_t getCentralityPercentile(UInt_t centE=kTOFRPC);       // Return centrality percentile [0.0-100.0%] of total cross section - Error gives 101% back
    TH1F*   getCentralityPercentileHist(UInt_t centE=kTOFRPC) const;

    Int_t   getCentralityClass5(TString estimator);              // Return centrality 5%  classes binned from centrality percentile
    Int_t   getCentralityClass10(TString estimator);             // Return centrality 10% classes binned from centrality percentile

    //maker
    TH2D*    makeEPresolution(TProfile2D* hist, Bool_t calcChi=kFALSE);
    TH1*     makeEPresolution(TH3* hist, Bool_t calcChi=kFALSE);
    Double_t ModifiedBesselI( Int_t n , Double_t x) const;
    Double_t ComputeResolution( Double_t x, Int_t n=1 ) const;
    Double_t FindXi( Double_t res, Double_t prec=1e-6, Int_t n=1) const;

    Int_t    findEdgeBin(TH1F *htemp);
    TH1F*    makePercentiles(TH1F* hist, Float_t fractionXsection=100., Int_t direction=-1);
    //TH1F*    makePercentiles(TH2F* hist, Float_t fractionXsection=100., Int_t direction=-1); //FIXME
    Int_t    getCentralityClassNbins(UInt_t centC=k10);
    Float_t* getCentralityClassArray(UInt_t centC=k10);
    TH1F*    makeClasses(TH1F *h, Float_t fractionXsection=100., UInt_t centC=k10, Int_t direction=-1);
    
    //helper
    Int_t    getNbins(TString estimator);
    Int_t    getNbins(UInt_t centE=kTOFRPC, UInt_t centC=k10);
    Float_t* getUpEdgeArray(UInt_t centE=kTOFRPC, UInt_t centC=k10);
    Float_t* getBinCenterArray(UInt_t centE=kTOFRPC, UInt_t centC=k10);
    vector<TString> getLabelArray(UInt_t centE=kTOFRPC, UInt_t centC=k10);

    // estimator getter:
    Int_t   getMdcWiresOuterMod();
    Float_t getSelectedParticleCandSecNorm();
    Float_t getSelectedParticleCandCorrPerWire();
    Float_t getFWSumChargeSpec();
    Float_t getFWSumZ(Int_t minZ=1, Int_t maxZ=99, UInt_t SubEvent=0);
    Float_t getEt();
    Float_t getRatioEtEz();
    Float_t getDirectivity();

    void    reset();
    Int_t   print();
    
private:

   // const char* StringEventPlaneParameter[kNumEventPlaneParameter];   // variable names



    HCategory* fParticleCandCat;     // ParticleCand category
    HCategory* fParticleEvtInfoCat;  // ParticleEvtInfo category
    HCategory* fCatWallHit;          // ForwardWall category
    
    TString    fParameterFile;          // location and name of ParameterFile
    TFile*     fFile;
    
    TRandom*  fRandom;
    
    // -------------------------------------------------------------------------
    UInt_t  fDayOfYear;                 // day of the year extracted from filename (sim=0)
    UInt_t  currentEventSeqNumber;      // running seq. num. of event
    
    // flags
    TString fPostFix;
    Bool_t  isSimulation;               // for simulation if catGeantKine is available

    Float_t fReferenceMeanSelTrack;     // track running mean reference point
    
    //FIXME fixed cuts not needed
    vector<Float_t>  fFWminBeta;
    vector<Float_t>  fFWmaxBeta;
    vector<Float_t>  fFWminCharge;
    vector<Float_t>  fFWChargeCuts;    

    vector<vector<TH1*> > fCentralityHist;
    vector<TH1*>          fCentralityPercentileHist;
    vector<TH1*>          fEstimatorHist;
    vector<TH2*>          fEventPlaneCorrectionHist;
    vector<vector<TH1*> > fFWCutValuesHist;
    vector<Bool_t>        useFWCut;
    Bool_t                excludeNoisyFWcells;
    
    vector<SimpleQVector*> arrayOfHits;
    vector<Int_t>         iFWHitvector;   // list of selected and shuffled FW hit-num. 
    vector<Float_t>       vQPhi;    
    


    UInt_t fEventPlaneCorrectionFlag;
    Bool_t fQVectorCalcDone;
    
    inline Bool_t  isFlagSet(UInt_t flag, UInt_t status){ return (flag==(status&flag));}
    inline Float_t getPhi(Float_t x, Float_t y){ return (x==0 && y==0) ? -1 : TVector2::Phi_0_2pi( TMath::Pi()+TMath::ATan2(-y,-x) );}

    ClassDef(HParticleEvtCharaBK,3)

};

#endif
