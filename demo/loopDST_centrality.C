#include "hades.h"
#include "hloop.h"
#include "htime.h"
#include "hcategory.h"
#include "hcategorymanager.h"
#include "hparticlecand.h"
#include "hparticletracksorter.h"
#include "hparticlebooker.h"
#include "hparticletool.h"
#include "hparticledef.h"
#include "hparticleevtinfo.h"
#include "henergylosscorrpar.h"
#include "hphysicsconstants.h"

#include "hparticleevtchara.h"

#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TProfile.h"


#include <iostream>
#include <vector>
using namespace std;

Int_t loopDST_centrality(
              TString infileList="/lustre/nyx/hades/dst/apr12/gen8/108/root/be1210819405908.hld_dst_apr12.root",
          TString outfile="test.root",Int_t nEvents=1000)
{
    //  infileList : comma seprated file list "file1.root,file2.root" or "something*.root"
    //  outfile    : optional (not used here) , used to store hists in root file
    //  nEvents    : number of events to processed. if  nEvents < entries or < 0 the chain will be processed

    Bool_t isSimulation = kFALSE;

    //-------------------------------------------------
    // create loop obejct and hades
    HLoop loop(kTRUE);
    //-------------------------------------------------

    //-------------------------------------------------
    // list of all files with working sectors
    if(!isSimulation) loop.readSectorFileList("/lustre/nyx/hades/dst/apr12/gen8/sector_selection/FileListHadron.list",kFALSE,kFALSE);
    //-------------------------------------------------


    //-------------------------------------------------
    // reading input files and decalring containers
     Bool_t ret =kFALSE;
    if(infileList.Contains(",")){
    ret = loop.addMultFiles(infileList);      // file1,file2,file3
    } else{
    ret = loop.addFiles(infileList); // myroot*.root
    }

    if(ret == 0) {
    cout<<"READBACK: ERROR : cannot find inputfiles : "<<infileList.Data()<<endl;
    return 1;
    }

    //if(!loop.setInput("")) {   // all input categories
    //if(!loop.setInput("-*,+HParticleCand,+HParticleEvtInfo")) {
    if(!loop.setInput("-*,+HParticleEvtInfo, +HParticleCand, +HWallHit")) {   // for FWSumChargeSpec HWallHit needed
    cout<<"READBACK: ERROR : cannot read input !"<<endl;
    exit(1);
    } // read all categories
    loop.printCategories();
    loop.printChain();
    //-------------------------------------------------


    //-------------------------------------------------
    //parameters
    HEnergyLossCorrPar dEdxCorr;
    dEdxCorr.setDefaultPar("apr12");
    //-------------------------------------------------

    //-------------------------------------------------
    //Event Chara
    HParticleEvtChara evtChara;

    TString ParameterfileCVMFS;
    //     ParameterfileCVMFS = "/cvmfs/hades.gsi.de/param/eventchara/centrality_hydra2_4.9i_19012016.root";
    //   ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_apr12_gen8_2018_07_pass10.root";

       ParameterfileCVMFS = "/cvmfs/hades.gsi.de/param/eventchara/centrality_epcorr_apr12_gen8_2018_07.root";

    if(isSimulation)
           ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_sim_au1230au_gen8a_UrQMD_minbias_2018_06.root";
    
    if(!evtChara.setParameterFile(ParameterfileCVMFS)){
           cout << "Centrality Parameterfile not found !!! " << endl;
           return kFALSE;  
    }
    if(!evtChara.init()) {
           cout << "HParticleEvtChara not init!!! " << endl;
           return kFALSE;
    }

    evtChara.printCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k5);
    evtChara.printCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);

    Int_t eCentEst = HParticleEvtChara::kTOFRPC;
    Int_t eEPcorrection = HParticleEvtChara::kShiftFW;
    cout << "eEPcorrection: " << eEPcorrection << endl;

    Int_t eEPcorrection2 = HParticleEvtChara::kReCentering;
    cout << "eEPcorrection2: " << eEPcorrection2 << endl;
    
    Int_t eEPcorrectionNoCorr = HParticleEvtChara::kNoCorrection;
    TString sName;
    TList*  fHistList = new TList();
   
    Int_t nBinCent;
    Double_t XaxisCent[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85};
    nBinCent = sizeof(XaxisCent) / sizeof(XaxisCent[0]) - 1;
    
    //-------------------------------------------------

    //#################################################
    //#################################################
    // add your histograms here
    TFile* out = new TFile(outfile.Data(),"RECREATE");
    out->cd();

    // centrality-hist:
    TH1F*  hTOFRPC            = new TH1F("TOFRPC","TOFRPC",260,0,260);
    fHistList->Add(hTOFRPC);
    
    TH1F*  hCentralityTOFRPC            = new TH1F("CentralityTOFRPC","Centrality Bins TOFRPC",22,0,22);
    fHistList->Add(hCentralityTOFRPC);
    TH1F*  hCentralityTOFRPCPercentile  = new TH1F("CentralityTOFRPCPercentile","Centrality Percentile Bins TOFRPC",nBinCent,XaxisCent);
    fHistList->Add(hCentralityTOFRPCPercentile);

    TH1F*  hEventPlaneNoCorr  = new TH1F("EventPlaneNoCorr","EventPlaneNoCorr",180,0, TMath::TwoPi());
    fHistList->Add(hEventPlaneNoCorr);    
    TH1F*  hEventPlane  = new TH1F("EventPlane","EventPlane",180,0, TMath::TwoPi());
    fHistList->Add(hEventPlane);
    
    TH1F*  hEventPlane2  = new TH1F("EventPlane2","EventPlane2",180,0, TMath::TwoPi());
    fHistList->Add(hEventPlane2);
    
    sName = Form("ProtonsVs%s_percentile", evtChara.getStringCentralityEstimator(eCentEst).Data());
    TProfile* hProtonsVsCentralityPercentile = new TProfile(sName.Data(), sName.Data(), nBinCent, XaxisCent, "s" );
    fHistList->Add(hProtonsVsCentralityPercentile);

    sName = Form("PionsVs%s_percentile", evtChara.getStringCentralityEstimator(eCentEst).Data());
    TProfile* hPionsVsCentralityPercentile = new TProfile(sName.Data(), sName.Data(), nBinCent, XaxisCent, "s" );
    fHistList->Add(hPionsVsCentralityPercentile);



   //#################################################
   //#################################################


    //-------------------------------------------------
    // input data
    HCategory* candCat    = (HCategory*)HCategoryManager::getCategory(catParticleCand);
    HCategory* evtInfoCat = (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo);
    //-------------------------------------------------


    //-------------------------------------------------
    // event loop starts here
    Int_t entries = loop.getEntries();
    if(nEvents < entries && nEvents >= 0 ) entries = nEvents;
    TString filename;
    Int_t sectors [6];
    Float_t event_weight = 1;
    Float_t fEventPlaneNoCorr = 0;
    Float_t fEventPlane = 0;
    Int_t nProtons = 0;
    Int_t nPions   = 0;

    for (Int_t i = 0; i < entries; i++) {
    Int_t nbytes =  loop.nextEvent(i);             // get next event. categories will be cleared before
    if(nbytes <= 0) { cout<<nbytes<<endl; break; } // last event reached
        if(i%1000 == 0) cout<<"event "<<i<< "\r" << flush;

    loop.getSectors(sectors); // fill sector array

    if(loop.isNewFile(filename)){
       if(!isSimulation) filename = HTime::stripFileName(filename,kTRUE,kFALSE);
       cout << "Day: " << evtChara.loadDayOfYear() << " filename:  " << filename << endl;
    }

    //-------------------------------------------------
        // summary event info object
    HParticleEvtInfo* evtInfo=0;
    evtInfo = HCategoryManager::getObject(evtInfo,evtInfoCat,0 );

    if(evtInfo&&!evtInfo->isGoodEvent(
                    //Particle::kGoodTRIGGER|     // no trigger selection to also use PT2 events
                      Particle::kGoodVertexClust|
                      Particle::kGoodVertexCand|
                      Particle::kGoodSTART//|
                   //   Particle::kNoPileUpSTART|   // min. bias selection only vertex cuts
                   //   Particle::kNoVETO|
                   //   Particle::kGoodSTARTVETO|
                   //   Particle::kGoodSTARTMETA
                     )) continue;
    
    event_weight = evtChara.getEventWeight();   // event_weight dependent if PT2(down-scaled) or PT3
    //-------------------------------------------------

    hTOFRPC->Fill(evtChara.getCentralityEstimator(HParticleEvtChara::kTOFRPC), event_weight);
    hCentralityTOFRPC->Fill(evtChara.getCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k5), event_weight);
    hCentralityTOFRPCPercentile->Fill(evtChara.getCentralityPercentile(HParticleEvtChara::kTOFRPC), event_weight);
    fEventPlaneNoCorr = evtChara.getEventPlane(HParticleEvtChara::kNoCorrection,0);
    if(fEventPlaneNoCorr!=1) hEventPlaneNoCorr->Fill(fEventPlaneNoCorr, event_weight);
    fEventPlane = evtChara.getEventPlane(HParticleEvtChara::kShiftFW);
    if(fEventPlane!=1) hEventPlane->Fill(fEventPlane, event_weight);
    fEventPlane = evtChara.getEventPlane(HParticleEvtChara::kReCentering);
    if(fEventPlane!=1) hEventPlane2->Fill(fEventPlane, event_weight);


    //cout << "EventCorr:  " << fEventPlane << "  noEPcorrections: " << fEventPlaneNoCorr << "  eventweight:  " << event_weight <<  endl;
    /*
    nProtons = 0;
    nPions   = 0;
    //-------------------------------------------------
        // loop over particle candidates in event
    if(candCat){

        Int_t size = candCat->getEntries();
        HParticleCand* cand=0;
        for(Int_t j = 0; j < size; j++)
        {
        cand = HCategoryManager::getObject(cand,candCat,j);
        if(cand) {

            if(!loop.goodSector(cand->getSector()))  continue;  // skipp inactive sectors
            if(!cand->isFlagBit(kIsUsed)) continue;
            if(cand->getPID()<0)          continue;
            //#################################################
            //#################################################
            // fill your histograms here

            if(cand->getPID()==14) nProtons++;
            if(cand->getPID()==8 || cand->getPID()==9) ++nPions;
            //#################################################
            //#################################################
            //-------------------------------------------------
        }
        } // end cand loop
    } // end cand cat
    
    hProtonsVsCentralityPercentile->Fill(evtChara.getCentralityPercentile(eCentEst), nProtons, event_weight);
    hPionsVsCentralityPercentile->Fill(evtChara.getCentralityPercentile(eCentEst), nPions, event_weight);
   */
    //-------------------------------------------------

    } // end eventloop

    //#################################################
    //#################################################
    // write your histograms here
    out->cd();
    
    fHistList->Write();


    out->Save();
    out->Close();
    //#################################################
    //#################################################

    delete gHades;
    return 0;
}
