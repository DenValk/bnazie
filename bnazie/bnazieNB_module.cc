/*
 *Low Energy Metrics Analysis program by Jason Stock of South Dakota School of Mines and Technology
 *September 2015.
 *
 *
 * ToDo:
 * Add TTree to show parent particle data (EvtNum, PDG num, Initial Energy, Initial Position).
 * Fix TNtuples to actually store data.
 *
 */

#ifndef bnazieNB_Module
#define bnazieNB_Module

// LArSoft includes
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/RecoBase/Hit.h"
//#include "RecoBase/OpFlash.h"
#include "lardata/RecoBase/OpHit.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <http://root.cern.ch/root/html532/ClassIndex.html>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TProfile.h"

// C++ Includes
#include <map>
#include <iomanip>
#include <sstream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

namespace bnazieNB {

  class bnazieNB : public art::EDAnalyzer
  {
  public:
 
    explicit bnazieNB(fhicl::ParameterSet const& parameterSet);

    virtual void beginJob();
    virtual void reconfigure(fhicl::ParameterSet const& parameterSet) override;
    virtual void analyze (const art::Event& event) override;
    virtual void endJob();
	

  private:

    std::string fPositionsLabel;
    std::string fHitLabel;
    std::string fCheatLabel;
    std::string fPhotonsLabel;

    TTree* fMCData;  
    TTree* fMetricsNtuple;  
    TTree* fCheatMetricsNtuple;  
    TProfile* fPeakAmpHist;
    TProfile* fMaxPeakAmpHist;
    TProfile* fCheatPeakAmpHist;
    TProfile* fMaxCheatPeakAmpHist;
    TProfile* fSummedIntegralHist; 
    TProfile* fCheatSummedIntegralHist;
    TProfile* fNumHitsHist;
    TProfile* fCheatNumHitsHist;
    TProfile* fPEvXHist;
    TH1D*     fPEHist;
    TH1D*     fPENHist;
    TH3D*     partsTotal;
    TH3D*     partsMissed;
//    TH1D*     fTPEHist;
    
    Int_t fEvent, fPdgCode, numDetected, cNumDetected, numGenerated, totHits, numPhotoHit;    
    Double_t xPos, fPeakAmp, fMaxPeakAmp, fSumIntegral, totPeakAmp, totCheatPeakAmp, totIntegral, fCheatPeakAmp, fCheatMaxPeakAmp, fCheatSumIntegral;
    Int_t fNumHits;
    Int_t fCheatNumHits;



  }; // class bnazieNB

  bnazieNB::bnazieNB(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  void bnazieNB::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
 
    fMCData        = tfs->make<TTree>("fMCData", "MCInfo");
    fMCData->Branch("EvtNum", &fEvent, "EvtNum/I");
    fMCData->Branch("XPos", &xPos, "XPos/D");
    fMCData->Branch("PdgNum", &fPdgCode, "PdgCode/I");
 //   fMCData->Branch("Momentum", &fMCmomentum, "MCmomentum/D");
   // fTimeHist     = tfs->make<TH1D>("timehist",";Histogram of Times;",4000, 0, 4000);
   
    fMetricsNtuple        = tfs->make<TTree>("fMetrics", "Metrics");
    fMetricsNtuple->Branch("EvtNum", &fEvent, "EvtNum/I");
    fMetricsNtuple->Branch("XPos", &xPos, "XPos/D");
    fMetricsNtuple->Branch("MaxPeakAmp", &fMaxPeakAmp, "MaxPeakAmp/D");
    fMetricsNtuple->Branch("PeakAmp", &fPeakAmp, "PeakAmp/D");
    fMetricsNtuple->Branch("SumIntegral", &fSumIntegral, "SumIntegral/D");
    fMetricsNtuple->Branch("NumHits", &fNumHits, "NumHits/I");

    fCheatMetricsNtuple        = tfs->make<TTree>("fCheatMetrics", "CheatMetrics");
    fCheatMetricsNtuple->Branch("EvtNum", &fEvent, "EvtNum/I");
    fCheatMetricsNtuple->Branch("XPos", &xPos, "XPos/D");
    fCheatMetricsNtuple->Branch("CheatMaxPeakAmp", &fCheatMaxPeakAmp, "CheatMaxPeakAmp/D");
    fCheatMetricsNtuple->Branch("CheatPeakAmp", &fCheatPeakAmp, "CheatPeakAmp/D");
    fCheatMetricsNtuple->Branch("CheatSumIntegral", &fCheatSumIntegral, "CheatSumIntegral/D");
    fCheatMetricsNtuple->Branch("CheatNumHits", &fCheatNumHits, "CheatNumHits/I");

    fPeakAmpHist              = tfs->make<TProfile>("peakAmpHist",        ";Profile of Peak Amplitudes for each event vs Distance from APA (cm);",                      31, -55, 255.0);
    fCheatPeakAmpHist         = tfs->make<TProfile>("cheatPeakAmpHist",   ";Profile of Cheated Peak Amplitudes for each event vs Distance from APA (cm);",              31, -55, 255.0);
    fMaxPeakAmpHist           = tfs->make<TProfile>("maxPeakAmpHist",     ";Profile of Max Peak Amplitudes for each event vs Distance from APA (cm);",                      31, -55, 255.0);
    fMaxCheatPeakAmpHist      = tfs->make<TProfile>("cheatMaxPeakAmpHist",";Profile of Max Cheated Peak Amplitudes for each event vs Distance from APA (cm);",              31, -55, 255.0);
    fSummedIntegralHist       = tfs->make<TProfile>("summedIntHist",      ";Profile of the Summed Integral Charge (DUC) for each event vs Distance from APA (cm);",         31, -55, 255.0);
    fCheatSummedIntegralHist  = tfs->make<TProfile>("cheatSummedIntHist", ";Profile of the Cheated Summed Integral Charge (DUC) for each event vs Distance from APA (cm);", 31, -55, 255.0);
    fNumHitsHist              = tfs->make<TProfile>("numHitsHist",        ";Number of hits detected per event vs Distance from APA (cm);",                                31, -55, 255.0);
    fCheatNumHitsHist         = tfs->make<TProfile>("cheatNumHitsHist",   ";Profile of the Cheated Number of hits detected per event vs Distance from APA (cm);",                         31, -55, 255.0);
    fPEvXHist                 = tfs->make<TProfile>("pevXHist",           ";Profile of Detected Photons VS Distance from APA (cm);",                         31, -55, 255.0);
    fPEHist                   = tfs->make< TH1D >  ("PEHist",             ";Histogram of the Photon Hits from recob::OpHits_opflash__Reco.obj.fPE() ;",                 1000, 0, 2);
    fPENHist                  = tfs->make< TH1D >  ("PENHist",            ";Histogram of the Photon Hits from recob::OpHits_opflash__Reco.obj.fPE() (Set to 1 for PE>.3) ;",                 31, -55, 255);
    partsTotal                = tfs->make< TH3D >  ("partsAllPlot",       ";TH3 of the origin for all Particles Generated (Only valid for 1 particle per event simulations).;", 301, -50.5,  250.5, 251, -100.5, 150.5, 186, -15.5, 170.5);
    partsMissed               = tfs->make< TH3D >  ("partsMissedPlot",    ";TH3 of the origin for all Undetected Particles (Only valid for 1 particle per event simulations);", 301, -50.5,  250.5, 251, -100.5, 150.5, 186, -15.5, 170.5);
//    fTPEHist                  = tfs->make< TH1D >  ("TPEHist",            ";Histogram of the Photon Hits from recob::OpFlash_opflash__Reco.obj with TotalPE();",        10000, 0, 2000);

    ///////////////////////SET HISTOGRAM LABELS/////////////////
    fPeakAmpHist->SetXTitle("Distance from APA (cm)");
    fPeakAmpHist->SetYTitle("Peak Amplitude");
    fCheatPeakAmpHist->SetXTitle("Distance from APA (cm)");
    fCheatPeakAmpHist->SetYTitle("Cheated Peak Amplitude");
    fMaxPeakAmpHist->SetXTitle("Distance from APA (cm)");
    fMaxPeakAmpHist->SetYTitle("Max Peak Amplitude");
    fMaxCheatPeakAmpHist->SetXTitle("Distance from APA (cm)");
    fMaxCheatPeakAmpHist->SetYTitle("Max Cheated Peak Amplitude");
    fSummedIntegralHist->SetXTitle("Distance from APA (cm)");
    fSummedIntegralHist->SetYTitle("Summed Integral Charge per Event (Deposited Units Charge)");
    fCheatSummedIntegralHist->SetXTitle("Distance from APA (cm)");
    fCheatSummedIntegralHist->SetYTitle("Cheated Summed Integral Charge per Event (Deposited Units Charge)");
    fNumHitsHist->SetXTitle("Distance from APA (cm)");
    fCheatNumHitsHist->SetXTitle("Distance from APA (cm)");
	
    partsTotal->SetXTitle("X Position ( cm )");
    partsMissed->SetXTitle("X Position ( cm )");
    partsTotal->SetYTitle("Y Position ( cm )");
    partsMissed->SetYTitle("Y Position ( cm )");	
    partsTotal->SetZTitle("Z Position ( cm )");
    partsMissed->SetZTitle("Z Position ( cm )");	
	


    fPEHist->SetXTitle("Number of detected Photons");
    fPEHist->SetYTitle("Number of occurences");
    fPEvXHist->SetXTitle("Distance from APA");
    fPEvXHist->SetYTitle("Number of PE's");
//    fTPEHist->SetXTitle("Number of detected Photons");
//    fTPEHist->SetYTitle("Number of occurences");

   

 
  }
   

   //-----------------------------------------------------------------------
  void bnazieNB::endJob()
  {
	std::cout<< "\n================================================================================\n\n";
    std::cout << "Hit35t Parameters|\n    SummedIntegral/ParticleDetected : " << totIntegral/((double)numDetected) 
	          << "\n    NumHits per particle : " << ((double)totHits)/((double)numGenerated) 
			  << "\n    Avg Max Peak Amp per detected particel: " << ((double)totPeakAmp/(double)numDetected)
			  << "\n    Avg Max Peak Amp per generated particel: " << ((double)totPeakAmp/(double)numGenerated);
    std::cout << "\n\nCheated Parameters:\n    Average CheatedPeak Amplitude (of all detected particles). : " 
                << totCheatPeakAmp/((double)cNumDetected) 
              << "    Avg CheatedPeakAmplitude per particle: " << totCheatPeakAmp/numGenerated <<"\n\n";
    std::cout << "OpHit Parameters|  Number of detected photons/ParticlesProduced : " << ((double)numPhotoHit)/((double)numGenerated) 
              << "\n\nNumber of Detected Particles/Number of Particles Produced: " << ((double)numDetected/(double)numGenerated) 
              << "\n\n==============================================\n\n";
  }
  //-----------------------------------------------------------------------
  void bnazieNB::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fPositionsLabel         = parameterSet.get< std::string >("TruthLabel");
    fHitLabel               = parameterSet.get< std::string >("HitLabel");
    fCheatLabel             = parameterSet.get< std::string >("CheatLabel");
    fPhotonsLabel           = parameterSet.get< std::string >("PhotLabel");
    numDetected =0;
    cNumDetected=0;
    numGenerated=0;
    totHits=0;
    numPhotoHit=0;
    xPos=0.0;
    fPeakAmp=0.0;
    fMaxPeakAmp=0.0;
    fSumIntegral=0.0;
    totPeakAmp=0.0;
    totCheatPeakAmp=0.0;
    totIntegral=0.0;
  }

  //-----------------------------------------------------------------------
  void bnazieNB::analyze(const art::Event& evt) 
  {
	Bool_t detected = false;
	TLorentzVector posVec;	  

    art::Handle< std::vector< simb::MCParticle > > truthHandle;
    std::vector< art::Ptr < simb::MCParticle > > truthlist;
    if (evt.getByLabel(fPositionsLabel, truthHandle) )
      art::fill_ptr_vector(truthlist, truthHandle);
    

    art::Handle< std::vector< recob::Hit > > hitHandle;
    std::vector< art::Ptr< recob::Hit > > hitlist;
    if (evt.getByLabel(fHitLabel, hitHandle) )
      art::fill_ptr_vector(hitlist, hitHandle);

    art::Handle< std::vector< recob::Hit > > hitCheatHandle;
    std::vector< art::Ptr< recob::Hit > > hitCheatlist;
    if (evt.getByLabel(fCheatLabel, hitCheatHandle) )
      art::fill_ptr_vector(hitCheatlist, hitCheatHandle);

//    art::Handle< std::vector< recob::OpFlash > > opflashHandle;
//    std::vector< art::Ptr< recob::OpFlash > > opflashlist;
//    if (evt.getByLabel(fPhotonsLabel, opflashHandle) )
//      art::fill_ptr_vector(opflashlist, opflashHandle);

    art::Handle< std::vector< recob::OpHit > > ophitHandle;
    std::vector< art::Ptr< recob::OpHit > > ophitlist;
    if (evt.getByLabel(fPhotonsLabel, ophitHandle) )
      art::fill_ptr_vector(ophitlist, ophitHandle);

    
     //get first x position (starting point, or generation point)
    fEvent = evt.id().event();
    //This section needs to be rewritten to make a map of Hits to their coresponding MCParticles. Then I can get the position where each Hit started instead of the position for the first hit of each event.
    if(!truthlist.empty()){ 
     art::Ptr< simb::MCParticle > MCDat = *truthlist.begin();
     posVec = MCDat->Position(0);
     xPos = posVec(0);
     fPdgCode = MCDat->PdgCode();
     fMCData->Fill();
	numGenerated++;
    }


     /////////Det Reco Section

     std::vector<double> peakAmps;
     double sumIntegral = 0.0;
     fNumHits = 0;
     for ( std::vector< art::Ptr < recob::Hit > >::iterator hitIt = hitlist.begin(); hitIt != hitlist.end(); ++hitIt){
       art::Ptr<recob::Hit> thisHit = *hitIt;
       peakAmps.push_back(thisHit->PeakAmplitude());
       fPeakAmp += thisHit->PeakAmplitude();
       sumIntegral += thisHit->Integral();
       ++fNumHits;
       totHits++;
     }
     if(!peakAmps.empty())
     { 
       fMaxPeakAmp = *std::max_element(peakAmps.begin(), peakAmps.end()); 
	   totPeakAmp += fMaxPeakAmp;
	   detected = true;
	   numDetected++;
//     }
//     else if(peakAmps.empty()) { fMaxPeakAmp = 0.0; }
       fSumIntegral = sumIntegral;
       totIntegral += fSumIntegral;
       fMetricsNtuple->Fill();
       fMaxPeakAmpHist->Fill(xPos, fMaxPeakAmp);
       fPeakAmpHist->Fill(xPos, fPeakAmp);
       fSummedIntegralHist->Fill(xPos, fSumIntegral);
       fNumHitsHist->Fill(xPos, fNumHits);
     } //Modified to write only when peakamps is not empty.

	 
	 
	 //Fill Appropriate TH3Ds.
	partsTotal->Fill(posVec(0), posVec(1), posVec(2)); //always fill Total, including 0,0,0 pileup
	if(detected==true && !truthlist.empty()){
	  partsMissed->Fill(posVec(0), posVec(1), posVec(2)); //Only fill missed when there was a particle, and it didn't make a hit35t signla.
	}
	
     //CHEAT SECTION!!!!!!!!!!!!!!1

     std::vector<double> cheatPeakAmps;
     double cheatSumIntegral = 0.0;
     fCheatNumHits = 0;
     for ( std::vector< art::Ptr < recob::Hit > >::iterator hitIt = hitCheatlist.begin(); hitIt != hitCheatlist.end(); ++hitIt){
       art::Ptr<recob::Hit> thisHit = *hitIt;
       cheatPeakAmps.push_back(thisHit->PeakAmplitude());
       cheatSumIntegral += thisHit->Integral();
       ++fCheatNumHits;
      }
     if(!cheatPeakAmps.empty())
     { 
       fCheatMaxPeakAmp = *std::max_element(cheatPeakAmps.begin(), cheatPeakAmps.end()); 
	   totCheatPeakAmp += fCheatMaxPeakAmp;
	   cNumDetected++;
//     }
//     else if(cheatPeakAmps.empty()) { fCheatMaxPeakAmp = 0.0; }
       fCheatSumIntegral = cheatSumIntegral;
       fCheatMetricsNtuple->Fill();
       fMaxCheatPeakAmpHist->Fill(xPos, fCheatMaxPeakAmp);
       fCheatPeakAmpHist->Fill(xPos, totCheatPeakAmp);
       fCheatSummedIntegralHist->Fill(xPos, fCheatSumIntegral);
       fCheatNumHitsHist->Fill(xPos, fCheatNumHits);
     }  //Only fill when there was a particle. Do not write bogus zeros.

     //Photons Histograms Section

     for ( std::vector< art::Ptr <recob::OpHit > >::iterator peIT = ophitlist.begin(); peIT != ophitlist.end(); ++peIT){
       double peVal = 0.0;
       art::Ptr<recob::OpHit> thisPoint = *peIT;
	   numPhotoHit++;
       peVal = thisPoint->PE();
       fPEHist->Fill(peVal);
       if( peVal>0.3 ){
         fPEvXHist->Fill(xPos, peVal);
         fPENHist->Fill(xPos);
       }
     }

//     for ( std::vector< art::Ptr <recob::OpFlash > >::iterator tpeIT = opflashlist.begin(); tpeIT != opflashlist.end(); ++tpeIT){
//       double tpeVal = 0.0;
//       art::Ptr<recob::OpFlash> thisPoint = *tpeIT;
//       tpeVal = thisPoint->TotalPE();
//       fTPEHist->Fill(tpeVal);
//     }
       
//    art::Ptr< recob::OpFlash > tpeIT = *opflashlist.begin();
//    const double tpeVal = tpeIT->TotalPE();
//    fTPEHist->Fill(tpeVal);
       
   

  } // bnazieNB::analyze()
  
  
  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see bnazieNB.fcl for more information.
  DEFINE_ART_MODULE(bnazieNB)

} // namespace bnazieNB

#endif // bnazieNB_Module

