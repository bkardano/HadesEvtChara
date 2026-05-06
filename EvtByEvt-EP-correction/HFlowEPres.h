#ifndef HFLOWEPRES_H
#define HFLOWEPRES_H

#include "hades.h"
#include "htool.h"
#include "htaskset.h"
#include "hruntimedb.h"
#include "hrecevent.h"
#include "hreconstructor.h"
#include "hcategorymanager.h"
#include "hcategory.h"

#include "hparticlestructs.h"

//--------category definitions---------
#include "hparticledef.h"
#include "hstartdef.h"
#include "hgeantdef.h"
//-------------------------------------

//-------objects-----------------------
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticleevtinfo.h"
#include "hstart2hit.h"
#include "hstart2cal.h"
#include "hgeantkine.h"
#include "hparticlegeantevent.h"
#include "htboxchan.h"
#include "htofhit.h"
#include "hrpchit.h"

//-------------------------------------
#include "tofdef.h"
#include "rpcdef.h"
#include "walldef.h"
#include "hwallhitsim.h"
#include "henergylosscorrpar.h"

//-------------------------------------
// root stuff
#include "TString.h"
#include "TFile.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TCanvas.h"


//-------------------------------------
// standard stuff
#include <iostream>
#include <cstdlib>
using namespace std;
//-------------------------------------


// HFlowEPCorr:
// Description: Class to organise event-by-event EP resolution correction
// author: Behruz Kardan 24.09.2024


class HFlowEPres : public HReconstructor
{
    // put your global vars here


private:

    //-------------------------------------
    // config 
    TString   Beamtime;  //FIXME  with  gHades->getBeamTimeID()==kApr12; kMar19)
    Bool_t    isSimulation;  // will be evaluated by HGeantKine Cat later
    Bool_t    debug;
    Bool_t    makeQA;

    Int_t     fieldConfig;
    Int_t     nHarm;
    //--------------------------------------------------

    TString AnalysisParaCode;   
    Bool_t  useDSTPID;
    Bool_t  useMassPID;
    Bool_t  usePolarityPID;

    //--------------------------------------------------        
    vector<TH2*>   fHistEPresCorr;
    //--------------------------------------------------

public:
    HFlowEPres();
    HFlowEPres(const Text_t *name = "HFlowEPres",const Text_t *title ="HFlowEPres"){
        AnalysisParaCode = name;
        
        isSimulation = kFALSE;
        debug        = kFALSE;
        makeQA       = kFALSE;

        //Beamtime="mar19";//apr12 //cg
        fieldConfig = 1;   // 1 -nom   2- rev
        nHarm =8;

    }

    virtual ~HFlowEPres(){};
    
    Bool_t      init(void)    {return 0; }
    Bool_t      reinit(void)  {return 0; }
    Bool_t      finalize(void){return 0; }
    Int_t       execute(void) {return 0; }
    
    void        SetBeamtime(TString f)              { Beamtime = f;};
    void        SetDebug(Bool_t deb)                { debug = deb; }
    void        SetMakeQA(Bool_t deb)               { makeQA = deb; }
    void        SetAnalysisParaCode(TString s)      { AnalysisParaCode = s;}
    void        SetFieldConfig(Int_t i)             { fieldConfig = i;};
    //--------------------------------------------------

    //-----------------------------------------------------------------------
    Bool_t LoadHistograms(TFile *outputFile)
    {
       TString sName;
       // get the pointers to all output histograms before calling Finish()
       if (outputFile) {
          // Get the common histograms from the output list

          if (!outputFile->IsOpen()) {
             cout << ">>>>>WARNING: the file does not exist!!! " << endl;
             return kFALSE;
          }
          
          cout << " |||  Loading 2D-Histogramm \n" << endl;
          fHistEPresCorr.resize(nHarm);
          for(Int_t iHarm=0; iHarm<nHarm; iHarm++){
              sName = Form("EPResolution_TOFRPC_FWCharge_ShifFWFlattWCharge_EPresPara_R%d", iHarm+1);
              fHistEPresCorr[iHarm] = (TH2*) outputFile->FindObjectAny(sName.Data());
              if(fHistEPresCorr[iHarm]){
                  fHistEPresCorr[iHarm]->Print();
              }
              else{
                  cout << "No "<< sName << " found! I am out!!!" << endl;
                  return kFALSE;
              }
          }
       }
       else { 
          cout << "outputFile pointer is empty" << endl;
          return kFALSE;
       }
       //fHistList->Print();
       return kTRUE;
    }
    
    //--------------------------------------------------
    Double_t GetEPresCorr(Int_t iHarm, Int_t TOFRPC, Int_t FW){
        // track cuts would go here
        if(iHarm >0 && iHarm < nHarm+1 && fHistEPresCorr[iHarm-1]){
            return fHistEPresCorr[iHarm-1]->GetBinContent(fHistEPresCorr[iHarm-1]->FindBin(TOFRPC, FW));
            //Double_t error  = fHistEPresCorr[iHarm-1]->GetBinError(fHistEPresCorr[iHarm-1]->FindBin(TOFRPC, FW));
        }
        return -1;
    }
    //-----------------------------------------------------------------------
    
    inline Bool_t isFlagSet(UInt_t flag, UInt_t status){ return (flag==(status&flag));}
    //ClassDef(HFlowEPres, 0);
};
// -------------------------------------------------------------------------

#endif