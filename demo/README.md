 quick how-To:

1.  login  to kronos.hpc.gsi.de
2.  be sure to use the right defall.sh

      . ./defall.sh

3.  that MYHADDIR=... is set to your private include/ & lib/-directory
    check with

     echo ${MYHADDIR}
    or
     export

3.1 if not, edit you ./defall.sh and add a directory after MYHADDIR=
    something like:
    
     export MYHADDIR=/lustre/nyx/hades/user/USERNAME/lib/hydra

4.  be sure that in your defall.sh ROOT_INCLUDE_PATH=... is set to the root-include-
    and your private include-directory   
    check with

     echo ${ROOT_INCLUDE_PATH}
    or
     export
     
     it should be something like this:
     /cvmfs/hades.gsi.de/install/5.34.34/hydra2-4.9xx/include:/lustre/nyx/hades/user/USERNAME/lib/hydra/include

4.1 if not, edit you ./defall.sh and add the directory after ROOT_INCLUDE_PATH=
    something like:
    
     export ROOT_INCLUDE_PATH=${HADDIR}/include:${MYHADDIR}/include


4.2  be sure that MYHADDIR contain lib/libEvtChara.so

      ls ${MYHADDIR}/lib

    and include/hparticleevtcharaBK.h

      ls ${MYHADDIR}/include

5.  if not, copy the folowing directory to any clean & calm place:

      mkdir ./CleanCalmPlace
      cd ./CleanCalmPlace
      svn co https://subversion.gsi.de/hades/analysis/bkardan/HParticleEvtChara/  ./

5.1  go to there and run

      cd ./evtchara
      make
      make install

6.  got to point 4.1 and check if libEvtChara.so exist
 
7.  add to analysis code the following parts:

7.1. include header-file

        #include "hparticleevtcharaBK.h" 

7.2. if you use Hloop check that following catagories are loaded:

     HParticleEvtInfo
     HParticleCand
     HWallHit

         if(!loop.setInput("-*,+HParticleEvtInfo, +HParticleCand, +HWallHit")){
             cout<<"READBACK: ERROR : cannot read input !"<<endl;
         }
  
7.3. before eventLoop:

         HParticleEvtCharaBK evtChara;
         TString ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_apr12_gen8_2019_02_pass30.root";
         
         //due to the overlap of day126  there is dedidacate param-file for the reverse-field runs
         if(isRevFieldDATA)   ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_apr12_gen8_revfield_2019_02_pass29.root";
         
         //if you need centrality-values for sim UrQMD gen8a you find them here
         //due missing primary light nuclei in UrQMD... there is no correction on the FW-EventPlane in this version!
         if(isSimulation)     ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_sim_au1230au_gen8a_UrQMD_minbias_2019_04_pass0.root";
         //if(isSimulation)     ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_sim_au1230au_gen9vertex_UrQMD_minbias_2019_04_pass0.root";

         if(!evtChara.setParameterFile(ParameterfileCVMFS)){
             cout << "Parameterfile not found !!! " << endl;return kFALSE;
         }
  
         if(!evtChara.init()) {
             cout << "HParticleEvtCharaBK not init!!! " << endl;return kFALSE;
         }
         
         Int_t eCentEst   = HParticleEvtCharaBK::kTOFRPC;
         Int_t eCentClass = HParticleEvtCharaBK::k10;
         Int_t eEPcorr    = HParticleEvtCharaBK::kDefault;
         cout << "\t selected EPcorrection method is:  "  << evtChara.getStringEventPlaneCorrection(eEPcorr) << endl;
         
         
         evtChara.printCentralityClass(eCentEst, eCentClass);
  
7.4. inside event-loop:
  
          Float_t  event_weight = evtChara.getEventWeight();   // event_weight dependent if PT2(down-scaled) or PT3
          
          Int_t   CentralityClassTOFRPC = evtChara.getCentralityClass(eCentEst, eCentClass);
          // 10% Centrality-Classes:       1(0-10%) - 5(40-50%) ... 0:Overflow max:Underflow
          
          Float_t CentralityTOFRPC      = evtChara.getCentralityPercentile(eCentClass);
          // CentralityPercentile:         0 - 100% percentile of the total cross section
          //                                   101% Over-,Underflow or Outlier Events
          
          Float_t EventPlane            = evtChara.getEventPlane(eEPcorr);
          //EventPlane:                   0 - 2pi (in rad)  re-centered & scaled EP
          //                             -1   no EP could be reconstructed
          
          //for EP-resolution use the sub-event:
          Float_t EventPlaneA            = evtChara.getEventPlane(eEPcorr,1);
          Float_t EventPlaneB            = evtChara.getEventPlane(eEPcorr,2);

          //check if the EP(and subEvent EP) could be reconstructed, most of the case less than 4 Hits in FW
          if(EventPlane  == -1)  continue;
          if(EventPlaneA  == -1 || EventPlaneB  == -1 ) continue;

          check if you have any further selections on FW-hit multiplicites in you own code,
          this will bias the result!
          
7.5.  to calculate the EP-resolution one option is the following:
      (Ollitrault method.  - approximation valid for high chi-values)
      inside event-loop fill the difference of the two sub-event EP
      into a histogramm for each centrality class:
      
             Float_t deltaPhiAB =  TMath::Abs(TVector2::Phi_mpi_pi(EventPlaneA-EventPlaneB));
             h->Fill(deltaPhiAB);

      after your event-analysis you can finalize your results by calculating the EP-resolution for each centrality class
      by the ratio of event N(pi/2 < deltaEP_AB < pi)/N(0 < deltaEP_AB < pi):
      
          Int_t i90        = h->FindBin(0.5*TMath::Pi());
          Int_t i180       = h->FindBin(TMath::Pi());
          Double_t ratio   = h->Integral(i90,i180) /h->Integral(1,i180);
          Double_t dChi    = sqrt(-2.*TMath::Log(2.*ratio));
          Double_t dEPres1 = HParticleEvtCharaBK::ComputeResolution( dChi , 1);    // for first harmonic
          Double_t dEPres2 = HParticleEvtCharaBK::ComputeResolution( dChi , 2);    // for second ...
          printf("An estimate of the event plane resolution first order is: %f\n", dEPres1 );

7.6.  the second option will be described here soon here ;)
      (Voloshin method - precise iterative resolution search)
      inside event-loop fill the Cosine of the difference of the two sub-event EP
      into a histogramm for each centrality class:
     
          Float_t CosDeltaPhiAB =  TMath::Cos(EventPlaneA-EventPlaneB);
          h2->Fill(CosDeltaPhiAB);

      after your event-analysis you can finalize your results by calculating the EP-resolution for each centrality class
      by the mean-value:

          Double_t MeanCosDeltaPhiAB  = h2->GetMean();        //< cos n delta_phi>;
          Double_t dV = TMath::Sqrt(MeanCosDeltaPhiAB);
          Double_t dChiSub = FindXi(dV,1e-6, 1);              // Sub-Event Resolution-Parameter
          Double_t dChi    = TMath::Sqrt2()*dChiSub;
          Double_t dEPres1 = HParticleEvtCharaBK::ComputeResolution( dChi , 1);    // for first harmonic
          Double_t dEPres2 = HParticleEvtCharaBK::ComputeResolution( dChi , 2);    // for second ...
          printf("An estimate of the event plane resolution first order is: %f\n", dEPres1 );


7.7. to correct your flow-values use 1./dEPres1

8.  check if you added the lib into your Makefile in the analysis-dir
    change there this line:

           HYDRA_LIBS    += -lDst -lEvtChara

9.  compile your executable :  make

10. test your macro with data and if you real have the right parameterfile
    with right values for centrality and EPcorections, check the output of 
    
      evtChara.printCentralityClass(eCentEst, eCentClass);

11. and now happy plotting ;)

