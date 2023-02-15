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
//#include "HParticleEvtChara.h"

#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TProfile.h"


#include <iostream>
#include <vector>
using namespace std;

Int_t loopDST_full(
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
    //
	ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_apr12_gen8_2019_02_pass30.root";
    
    //due to the overlap of day126  there is dedidacate param-file for the reverse-field runs
    //if(isRevFieldDATA)   ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_apr12_gen8_revfield_2019_02_pass29.root";
    
    //if you need centrality-values for sim UrQMD gen8a you find them here
    //due missing primary light nuclei in UrQMD... there is no correction on the FW-EventPlane in this version!
    if(isSimulation)     ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_sim_au1230au_gen8a_UrQMD_minbias_2019_04_pass0.root";
    //if(isSimulation)     ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_sim_au1230au_gen9vertex_UrQMD_minbias_2019_04_pass0.root";

    if(!evtChara.setParameterFile(ParameterfileCVMFS)){
           cout << "Centrality Parameterfile not found !!! " << endl;
           return kFALSE;  
    }
    evtChara.setPostFix("");
    if(!evtChara.init()) {
           cout << "HParticleEvtChara not init!!! " << endl;
           return kFALSE;
    }

    evtChara.printCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k5);
    evtChara.printCentralityClass(HParticleEvtChara::kTOFRPC, HParticleEvtChara::k10);

   
   Int_t nBinCent;
   Double_t XaxisCent[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85};
   nBinCent = sizeof(XaxisCent) / sizeof(XaxisCent[0]) - 1;
   
    vector<vector<TH1*> >  fCentralityHist;
    vector<TH1*>       fCentralityHistPercentile;
    vector<TH1*>       fEstimatorHist;

    vector<TProfile*>  fProtonsVsCentralityPercentile;
    vector<TProfile*>  fPionsVsCentralityPercentile;
    TList*  fHistList = new TList();
    TList*  fHistListEstimator = new TList();
    TList*  fHistListPID = new TList();
    //-------------------------------------------------

    //#################################################
    //#################################################
    // add your histograms here
    TFile* out = new TFile(outfile.Data(),"RECREATE");
    out->cd();


    // automated all available centrality-hist:
    fCentralityHist.resize(HParticleEvtChara::kNumCentralityEstimator);
    for (int centEst = 0; centEst < (int)HParticleEvtChara::kNumCentralityEstimator; ++centEst) fCentralityHist[centEst].resize(HParticleEvtChara::kNumCentralityClass);
    fCentralityHistPercentile.resize(HParticleEvtChara::kNumCentralityEstimator);
    fEstimatorHist.resize(HParticleEvtChara::kNumCentralityEstimator);
    fProtonsVsCentralityPercentile.resize(HParticleEvtChara::kNumCentralityEstimator);
    fPionsVsCentralityPercentile.resize(HParticleEvtChara::kNumCentralityEstimator);
    TString sName;
    
    for (int centEst = 0; centEst < (int)HParticleEvtChara::kNumCentralityEstimator; ++centEst){
       sName = Form("%s_percentile", evtChara.getStringCentralityEstimator(centEst).Data());
       fCentralityHistPercentile[centEst] = new TH1F(sName.Data(), sName.Data(), nBinCent, XaxisCent );
       fHistList->Add(fCentralityHistPercentile[centEst]);
       
       sName = Form("ProtonsVs%s_percentile", evtChara.getStringCentralityEstimator(centEst).Data());
       fProtonsVsCentralityPercentile[centEst] = new TProfile(sName.Data(), sName.Data(), nBinCent, XaxisCent, "s" );
       fHistListPID->Add(fProtonsVsCentralityPercentile[centEst]);
       
       sName = Form("PionsVs%s_percentile", evtChara.getStringCentralityEstimator(centEst).Data());
       fPionsVsCentralityPercentile[centEst] = new TProfile(sName.Data(), sName.Data(), nBinCent, XaxisCent, "s" );
       fHistListPID->Add(fPionsVsCentralityPercentile[centEst]);

       sName = Form("%s", evtChara.getStringCentralityEstimator(centEst).Data());
       fEstimatorHist[centEst] =  new TH1F(sName.Data(), sName.Data(), evtChara.getCentralityEstimatorBinSize(centEst), 0, evtChara.getCentralityEstimatorBinMax(centEst));
       fHistListEstimator->Add(fEstimatorHist[centEst]);
       
       for (int centC = 0; centC < (int)HParticleEvtChara::kNumCentralityClass; ++centC){
           sName = Form("%s_%s", evtChara.getStringCentralityEstimator(centEst).Data(),  evtChara.getStringCentralityClass(centC).Data());
           Int_t nBins = evtChara.getNbins(centEst, centC);
           if(nBins>0){
               fCentralityHist[centEst][centC] = new TH1F(sName.Data(), sName.Data(), nBins, 0, nBins );
               fHistList->Add(fCentralityHist[centEst][centC]);
           }
       }
    }



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
    Int_t nProtons = 0;
    Int_t nPions   = 0;

    for (Int_t i = 0; i < entries; i++) {
    Int_t nbytes =  loop.nextEvent(i);             // get next event. categories will be cleared before
    if(nbytes <= 0) { cout<<nbytes<<endl; break; } // last event reached
        if(i%1000 == 0) cout<<"event "<<i<< "\r" << flush;

    loop.getSectors(sectors); // fill sector array

    if(loop.isNewFile(filename)){
       if(!isSimulation) filename = HTime::stripFileName(filename,kTRUE,kFALSE);
    }

    //-------------------------------------------------
        // summary event info object
    HParticleEvtInfo* evtInfo=0;
        evtInfo = HCategoryManager::getObject(evtInfo,evtInfoCat,0 );

    if(evtInfo&&!evtInfo->isGoodEvent(
                    //Particle::kGoodTRIGGER|     // no trigger selection to also use PT2 events
                      Particle::kGoodVertexClust|
                      Particle::kGoodVertexCand|
                      Particle::kGoodSTART|
                      Particle::kNoPileUpSTART|
                      Particle::kNoVETO|
                      Particle::kGoodSTARTVETO|
                      Particle::kGoodSTARTMETA
                     )) continue;
    
    event_weight = evtChara.getEventWeight();   // event_weight dependent if PT2(down-scaled) or PT3
    //-------------------------------------------------

    for (int centEst = 0; centEst < (int)HParticleEvtChara::kNumCentralityEstimator; ++centEst){
       if(fCentralityHistPercentile[centEst])fCentralityHistPercentile[centEst]->Fill(evtChara.getCentralityPercentile(centEst), event_weight);
       if(fEstimatorHist[centEst])fEstimatorHist[centEst]->Fill(evtChara.getCentralityEstimator(centEst), event_weight);
       for (int centC = 0; centC < (int)HParticleEvtChara::kNumCentralityClass; ++centC){
           if(fCentralityHist[centEst][centC]) fCentralityHist[centEst][centC]->Fill(evtChara.getCentralityClass(centEst, centC), event_weight);
       }
    }
    

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
    
    for (int centEst = 0; centEst < (int)HParticleEvtChara::kNumCentralityEstimator; ++centEst){
            if(fProtonsVsCentralityPercentile[centEst])
                   fProtonsVsCentralityPercentile[centEst]->Fill(evtChara.getCentralityPercentile(centEst), nProtons, event_weight);
            if(fPionsVsCentralityPercentile[centEst])
                   fPionsVsCentralityPercentile[centEst]->Fill(evtChara.getCentralityPercentile(centEst), nPions, event_weight);
    }
    //-------------------------------------------------

    } // end eventloop

    //#################################################
    //#################################################
    // write your histograms here
    out->cd();

    out->cd();
    out->mkdir("Estimator");
    out->cd("/Estimator/");
    fHistListEstimator->Write();

    out->mkdir("Centrality");
    out->cd("/Centrality/");
    fHistList->Write();

    out->mkdir("PID");
    out->cd("/PID/");
    fHistListPID->Write();

    out->Save();
    out->Close();
    //#################################################
    //#################################################

    delete gHades;
    return 0;
}
