/*
 *Metrics Analysis program by Jason Stock of South Dakota School of Mines and Technology
 *September 2015.
 *
 *
 * This Version of the bnazie analysis uses the backtracker and is built for events with multiple particles. As such, its structure is different than previous analysis used before Dec 2015.
 *
 */

#ifndef bnazie_Module
#define bnazie_Module

// LArSoft includes
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/PhotonBackTracker.h"

//#include "larsim/Simulation/LArG4Parameters.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
//#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/Exception.h"
//#include "fhiclcpp/ParameterSet.h"
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
#include "TProfile3D.h"
#include "TFile.h"

// C++ Includes
#include <map>
#include <iomanip>
#include <sstream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

namespace bnazie {

  class bnazie : public art::EDAnalyzer
  {
  public:
 
    explicit bnazie(fhicl::ParameterSet const& parameterSet);

    virtual void beginJob();
    virtual void reconfigure(fhicl::ParameterSet const& parameterSet) override;
    virtual void analyze (const art::Event& event) override;
    virtual void endJob();
	

  private:

    std::string fPositionsLabel;
    std::string fHitLabel;
    std::string fCheatLabel;
    std::string fPhotonsLabel;
    std::string fRescaleFile;
    std::string fRescaleObj;
    std::string fPBTRLabel;

    TFile* rescaleSourceFile;

    TTree* fMCData;  
    TTree* fMetricsNtuple;  
    TTree* fCheatMetricsNtuple;  

    TProfile* fTSummedIntegralProf; 
    TProfile* fTCheatSummedIntegralProf;
    TProfile* fCSummedIntegralProf; 
    TProfile* fCCheatSummedIntegralProf;
    TProfile* fBSummedIntegralProf; 
    TProfile* fBCheatSummedIntegralProf;
    TProfile* fPeakAmpProf;
    TProfile* fCheatPeakAmpProf;
    TProfile* fSummedIntegralProf; 
    TProfile* fCheatSummedIntegralProf;
    TProfile* fPEvXProf;
    TProfile3D* rescaleSourceProfile;

    TCanvas*  can6plot;

    TH1D*     fPEvXHist;
    TH1D*     fPENHist;
    TH1D*     fNumHitsHist;
    TH1D*     fCheatNumHitsHist;
    TH1D*     xHits;
    TH1D*     yHits;
    TH1D*     zHits;
    TH3D*     hitsDetected;
    
    Int_t fEvent, fPdgCode, numDetected, cNumDetected, numGenerated, totHits, numPhotoHit, numHits;    
    Int_t fNumHits, faultHits;
    Int_t fCheatNumHits;

    Double_t xPos, fPeakAmp, fMaxPeakAmp, fSumIntegral, totPeakAmp, totCheatPeakAmp, totIntegral, fCheatPeakAmp, fCheatMaxPeakAmp, fCheatSumIntegral, fxCorrection, fyCorrection, fzCorrection, fzOffset, fTau, fVel;

    //detector boundaries

    Double_t xLo, xHi, yLo, yHi, zLo, zHi;

    Bool_t fRescale;

    //art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<cheat::BackTracker> bt;
    art::ServiceHandle<cheat::PhotonBackTracker> pbt;   
    //cheat:: bt = lar::providerFrom<cheat::BackTracker>();
    geo::GeometryCore const* geom =  lar::providerFrom<geo::Geometry>();
    //art::TFileService const* tfs  =  lar::providerFrom<art::TFileService>();

  }; // class bnazie

  bnazie::bnazie(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  void bnazie::beginJob()
  {
 
    fMCData        = tfs->make<TTree>("fMCData", "MCInfo");
    fMCData->Branch("EvtNum", &fEvent, "EvtNum/I");
   
    fMetricsNtuple        = tfs->make<TTree>("fMetrics", "Metrics");
    fMetricsNtuple->Branch("EvtNum", &fEvent, "EvtNum/I");
    fMetricsNtuple->Branch("MaxPeakAmp", &fMaxPeakAmp, "MaxPeakAmp/D");
    fMetricsNtuple->Branch("PeakAmp", &fPeakAmp, "PeakAmp/D");
    fMetricsNtuple->Branch("SumIntegral", &fSumIntegral, "SumIntegral/D");
    fMetricsNtuple->Branch("NumHits", &fNumHits, "NumHits/I");

    fCheatMetricsNtuple        = tfs->make<TTree>("fCheatMetrics", "CheatMetrics");
    fCheatMetricsNtuple->Branch("EvtNum", &fEvent, "EvtNum/I");
    fCheatMetricsNtuple->Branch("CheatMaxPeakAmp", &fCheatMaxPeakAmp, "CheatMaxPeakAmp/D");
    fCheatMetricsNtuple->Branch("CheatPeakAmp", &fCheatPeakAmp, "CheatPeakAmp/D");
    fCheatMetricsNtuple->Branch("CheatSumIntegral", &fCheatSumIntegral, "CheatSumIntegral/D");
    fCheatMetricsNtuple->Branch("CheatNumHits", &fCheatNumHits, "CheatNumHits/I");

    fPeakAmpProf               = tfs->make<TProfile>("peakAmpHist",         ";Profile of Peak Amplitudes for each hit vs BackTracker Distance from APA (cm);",                36, -105.0, 255.0);
    fCheatPeakAmpProf          = tfs->make<TProfile>("cheatPeakAmpHist",    ";Profile of Cheated Peak Amplitudes for each hit vs BackTracker Distance from APA (cm);",        36, -105.0, 255.0);
    fTSummedIntegralProf       = tfs->make<TProfile>("tSummedIntHist",      ";Profile of the Summed Integral Charge (DUC) for each event vs Distance from APA (cm);",         36, -105.0, 255.0);
    fTCheatSummedIntegralProf  = tfs->make<TProfile>("tCheatSummedIntHist", ";Profile of the Cheated Summed Integral Charge (DUC) for each event vs Distance from APA (cm);", 36, -105.0, 255.0);
    fCSummedIntegralProf       = tfs->make<TProfile>("cSummedIntHist",      ";Profile of the Summed Integral Charge (DUC) for each event vs Distance from APA (cm);",         36, -105.0, 255.0);
    fCCheatSummedIntegralProf  = tfs->make<TProfile>("cCheatSummedIntHist", ";Profile of the Cheated Summed Integral Charge (DUC) for each event vs Distance from APA (cm);", 36, -105.0, 255.0);
    fBSummedIntegralProf       = tfs->make<TProfile>("bSummedIntHist",      ";Profile of the Summed Integral Charge (DUC) for each event vs Distance from APA (cm);",         36, -105.0, 255.0);
    fBCheatSummedIntegralProf  = tfs->make<TProfile>("bCheatSummedIntHist", ";Profile of the Cheated Summed Integral Charge (DUC) for each event vs Distance from APA (cm);", 36, -105.0, 255.0);
    fSummedIntegralProf        = tfs->make<TProfile>("summedIntHist",       ";Profile of the Summed Integral Charge (DUC) for each event vs Distance from APA (cm);",         36, -105.0, 255.0);
    fCheatSummedIntegralProf   = tfs->make<TProfile>("cheatSummedIntHist",  ";Profile of the Cheated Summed Integral Charge (DUC) for each event vs Distance from APA (cm);", 36, -105.0, 255.0);
    hitsDetected               = tfs->make< TH3D >  ("hitsDetectedPlot",    ";TH3 of the weighted average position for all Detected Hits (As given in Backtracker) ;",        36, -105.0, 255.0, 31, -260.0, 360.0, 33, -205.0, 455.0);
    fNumHitsHist               = tfs->make< TH1D >  ("numHits",             ";The number of hits detected at a distance x from the APA;",                                     36, -105.0, 255.0);
    fCheatNumHitsHist          = tfs->make< TH1D >  ("cheatNumHits",        ";The number of cheated hits detected at a distance x from the APA;",                             36, -105.0, 255.0);
    xHits                      = tfs->make< TH1D >  ("xHits",               ";The number of hits at different positions X;",                                                  36,-105,405);
    yHits                      = tfs->make< TH1D >  ("yHits",               ";The number of hits at different positions Y;",                                                  141, -705,705);
    zHits                      = tfs->make< TH1D >  ("zHits",               ";The number of hits at different positions Z;",                                                  211, -1055,1055);
    fPEvXHist                  = tfs->make< TH1D >  ("nPEvX",               ";The numebr of PEs at different positions X;", 36, -105.0, 255.0);

    if(fRescale==0){can6plot = tfs->make<TCanvas>("VersionValidationPlots", "vv6plotspread", 2550, 3300);}

    ///////////////////////SET HISTOGRAM LABELS/////////////////
    fPeakAmpProf->SetXTitle("Distance from APA (cm)");
    fPeakAmpProf->SetYTitle("Peak Amplitude");

    fCheatPeakAmpProf->SetXTitle("Distance from APA (cm)");
    fCheatPeakAmpProf->SetYTitle("Cheated Peak Amplitude");

    fSummedIntegralProf->SetXTitle("Distance from APA (cm)");
    fSummedIntegralProf->SetYTitle("Summed Integral Charge per Event (Deposited Units Charge)");

    fCheatSummedIntegralProf->SetXTitle("Distance from APA (cm)");
    fCheatSummedIntegralProf->SetYTitle("Cheated Summed Integral Charge per Event (Deposited Units Charge)");

    fNumHitsHist->SetXTitle("Distance from APA (cm)");
    fNumHitsHist->SetYTitle("Number of hits detected.");

    fCheatNumHitsHist->SetXTitle("Distance from APA (cm)");
    fCheatNumHitsHist->SetYTitle("Number of Cheated Hits Detected.");

    hitsDetected->SetXTitle("X Position ( cm )");
    hitsDetected->SetYTitle("Y Position ( cm )");
    hitsDetected->SetZTitle("Z Position ( cm )");

    fPEvXHist->SetXTitle("X Position (cm)");
    fPEvXHist->SetYTitle("Number if PE detected");
  }
  

   //-----------------------------------------------------------------------
  void bnazie::endJob()
  {
/*	std::cout<< "\n================================================================================\n\n";
    if(cNumDetected!=0 && numGenerated!=0 && numDetected!=0){
        std::cout << "HitFD Parameters|\n    SummedIntegral/ParticleDetected : " << totIntegral/((double)numDetected) 
	              << "\n    NumHits per particle : " << ((double)totHits)/((double)numGenerated) 
			      << "\n    Avg Max Peak Amp per detected particel: " << ((double)totPeakAmp/(double)numDetected)
			      << "\n    Avg Max Peak Amp per generated particel: " << ((double)totPeakAmp/(double)numGenerated);
        std::cout << "\n\nCheated Parameters:\n    Average CheatedPeak Amplitude (of all detected particles). : " 
                  << totCheatPeakAmp/((double)cNumDetected) 
                  << "    Avg CheatedPeakAmplitude per particle: " << totCheatPeakAmp/numGenerated <<"\n\n";
    }else{std::cout<<"numGen or cNumDetected is zero.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n";}
    std::cout << "OpHit Parameters|  Number of detected photons/ParticlesProduced : " << ((double)numPhotoHit)/((double)numGenerated) 
              << "\n\nNumber of Detected Particles/Number of Particles Produced: " << ((double)numDetected/(double)numGenerated) 
              << "\n\n==============================================\n\n"
              << "\n\n\n\nFault Hits = : " << faultHits << "\n"
              << "Total Hits = : " << numHits << "\n\n\n\n\n";*/

    //Make a TCanvas and place the histograms in it.
    if(fRescale==0){
        can6plot->SetFixedAspectRatio(1);
        can6plot ->Divide(2,3);
        can6plot ->cd(1);    
        fPeakAmpProf->Draw();
        can6plot -> cd(2);    
        fCheatPeakAmpProf->Draw();
        can6plot -> cd(3);
        fSummedIntegralProf->Draw();
        can6plot -> cd(4);
        fCheatSummedIntegralProf->Draw();
        can6plot -> cd(5);
        fNumHitsHist->Draw();
        can6plot -> cd(6);
        fCheatNumHitsHist->Draw();
        can6plot -> Update();    
        can6plot ->SaveAs("VersionValidationOutput.pdf");
        std::cout<<"can6plot output written";
    }


  }
  //-----------------------------------------------------------------------
  void bnazie::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fPositionsLabel         = parameterSet.get< std::string >("TruthLabel");
    fHitLabel               = parameterSet.get< std::string >("HitLabel");
    fCheatLabel             = parameterSet.get< std::string >("CheatLabel");
    fPhotonsLabel           = parameterSet.get< std::string >("PhotLabel");
    fRescale                = parameterSet.get< Bool_t     >("rescaleOn");
    fPBTRLabel              = parameterSet.get< std::string >("PhotonBackTrackerLabel");
    if(fRescale==true){
      fxCorrection            = parameterSet.get< Double_t    >("xCorrection");
      fyCorrection            = parameterSet.get< Double_t    >("yCorrection");
      fzCorrection            = parameterSet.get< Double_t    >("zCorrection");
      fzOffset                = parameterSet.get< Double_t    >("zOffset");     //Use Correction to align the origin of the scaling file with the origin of the detector geometry. Use the offset to shift to different regions of the detector.
      fRescaleFile            = parameterSet.get< std::string >("RescaleFile");
      fRescaleObj             = parameterSet.get< std::string >("RescaleTP3D");
      rescaleSourceFile       = new TFile( (const char*)(fRescaleFile.c_str() ) );
      rescaleSourceFile->ls();
      rescaleSourceProfile    = ( (TProfile3D*)rescaleSourceFile->Get( (const char*)(fRescaleObj.c_str()) ) );
      fTau                    = parameterSet.get< Double_t >("eLifetime");
      fVel                    = parameterSet.get< Double_t >("eDriftVelocity");

    }
    numDetected =0;
    cNumDetected=0;
    numGenerated=0;
    totHits=0;
    fPeakAmp=0.0;
    fSumIntegral=0.0;
    totCheatPeakAmp=0.0;
    totIntegral=0.0;
    faultHits=0;
    numHits=0;
    xLo=0.0;
    xHi=0.0;
    yLo=0.0;
    yHi=0.0;
    zLo=0.0;
    zHi=0.0;
    //Set boundaries
    for(geo::CryostatID const& cID: geom->IterateCryostatIDs()){
   // for(it = geom->begin_cryostat_id(); it = geom->end_cryostat_id(); it++){
      Double_t* bounds = new double[6];
      geom->CryostatBoundaries(bounds, cID);
      std::cout<< "Bounds xMax " << bounds[1] << "min " << bounds[0] <<"\n";
      xLo = std::min(xLo, bounds[0]);
      xHi = std::max(xHi, bounds[1]);
      yLo = std::min(yLo, bounds[2]);
      yHi = std::max(yHi, bounds[3]);
      zLo = std::min(zLo, bounds[4]);
      zHi = std::max(zHi, bounds[5]);
    }
    
    
  }

  //-----------------------------------------------------------------------
  void bnazie::analyze(const art::Event& evt) 
  {

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

    art::Handle< std::vector< recob::OpHit > > ophitHandle;
    std::vector< art::Ptr< recob::OpHit > > ophitlist;
    if (evt.getByLabel(fPhotonsLabel, ophitHandle) )
      art::fill_ptr_vector(ophitlist, ophitHandle);

    
    fEvent = evt.id().event();

    if(!truthlist.empty()){ 
     art::Ptr< simb::MCParticle > MCDat = *truthlist.begin();
     posVec = MCDat->Position(0);
     fPdgCode = MCDat->PdgCode();
     fMCData->Fill();
	numGenerated++;
    }

     //OpHit Section
     for(const art::Ptr<recob::OpHit>& ptrOpHit: ophitlist){
       try{
         recob::OpHit opHit = *ptrOpHit;
         std::vector<double> xyzOpPos = pbt->OpHitToXYZ(ptrOpHit);
         double nPE = opHit.PE();
         fPEvXHist->Fill(xyzOpPos.at(0), nPE);
       }
       catch(...){
       }
     }
     
     
     //Hit Section

     for ( std::vector< art::Ptr < recob::Hit > >::iterator hitIt = hitlist.begin(); hitIt != hitlist.end(); ++hitIt){
       art::Ptr<recob::Hit> thisHit = *hitIt;
       double fSumIntegral=0.0;
       fPeakAmp = thisHit->PeakAmplitude();
       fSumIntegral = thisHit->Integral();
       std::vector<double> hpos;
       
       // BACKTRACKER CORRECTION!!!!!!!!!!!!!!
       try{
           hpos = bt->HitToXYZ(thisHit);
           fNumHitsHist->Fill(hpos.at(0));
           hitsDetected->Fill(hpos.at(0), hpos.at(1), hpos.at(2));
           if(fRescale==true){
             Double_t weightingFactor = 0.0;
             Double_t sumStepWeight = 0.0;
             Int_t binXorigin = rescaleSourceProfile->GetXaxis()->FindBin(0.01*( 0.0 - fxCorrection ) );
             Int_t binX = rescaleSourceProfile->GetXaxis()->FindBin(0.01*(hpos.at(0) - fxCorrection)); //Factor of 100 hard coded for now. Bad. Profile in meters, Sim in centimeters.
             Int_t binY = rescaleSourceProfile->GetYaxis()->FindBin(0.01*(hpos.at(1) - fyCorrection));
             Int_t binZ = rescaleSourceProfile->GetZaxis()->FindBin(0.01*(hpos.at(2) - fzCorrection + fzOffset));

             while( binX >= binXorigin ){  //Currently only works for hits in the x>0 region.
               //Double_t stepPos = 15*binX+7.5;
               //Double_t stepPos = ( (Doubel_t*)(rescaleSourceProfile->GetXaxis()->GetBinCenter(binX))+fxCorrection ) ;
               Double_t stepPos = ( (Double_t)(rescaleSourceProfile->GetXaxis()->GetBinCenter(binX)) + fxCorrection );
               Double_t impurityLevel = 
                   rescaleSourceProfile->GetBinContent( binX, binY, binZ ); 
               //weightingFactor += exp(-(stepPos-fxCorrection)/(fVel*(2.0-impurityLevel)*fTau))/exp(02.0*(stepPos-fxCorrection)/(fVel*fTau));          
               //sumStepWeight += 1.0/exp(-(stepPos-fxCorrection)/(fVel*fTau));
               weightingFactor += exp(-(stepPos)/(fVel*(2.0-impurityLevel)*fTau))/exp(02.0*(stepPos)/(fVel*fTau));          
               sumStepWeight += 1.0/exp(-(stepPos)/(fVel*fTau));
               binX -= 1;
             }

             if(sumStepWeight!=0){weightingFactor /= sumStepWeight;}else{weightingFactor=0;}
             if(fPeakAmp!=0.0){fPeakAmpProf->Fill(hpos.at(0), fPeakAmp*weightingFactor );}
             fSummedIntegralProf->Fill(hpos.at(0), fSumIntegral*weightingFactor );
             if( hpos.at(1)<-555 ){fBSummedIntegralProf->Fill(hpos.at(0), fSumIntegral*weightingFactor);}
             else if( hpos.at(2)<30 && hpos.at(2)>-30 ){fCSummedIntegralProf->Fill(hpos.at(0), fSumIntegral*weightingFactor);}
             else if( hpos.at(2)>555 ) { fTSummedIntegralProf->Fill(hpos.at(0), fSumIntegral*weightingFactor);}

           }else{
             if(fPeakAmp!=0.0){fPeakAmpProf->Fill(hpos.at(0), fPeakAmp);}
             fSummedIntegralProf->Fill(hpos.at(0), fSumIntegral);
           }  

           fMetricsNtuple->Fill();
           { // The filling of XYZ Histograms (Separate function for easy removal.
             xHits->Fill(hpos.at(0));
             yHits->Fill(hpos.at(1));
             zHits->Fill(hpos.at(2));
           }
           numHits++;
       }
       catch(cet::exception e){
//           std::cout << "ERROR: Backtracker failed. This (hitfd)hit's peak amplitude was: " << fCheatPeakAmp << "\n";
           faultHits++;
       }
     }
	 
	 
     //CHEAT SECTION!!!!

     for ( std::vector< art::Ptr < recob::Hit > >::iterator hitIt = hitCheatlist.begin(); hitIt != hitCheatlist.end(); ++hitIt){
       art::Ptr<recob::Hit> thisHit = *hitIt;
       double fCheatSumIntegral;
       fCheatPeakAmp = thisHit->PeakAmplitude();
       fCheatSumIntegral = thisHit->Integral();
       std::vector<double> hpos;

       try{
           hpos = bt->HitToXYZ(thisHit);
           fCheatNumHitsHist->Fill(hpos.at(0));
           fCheatMetricsNtuple->Fill();
           fCheatPeakAmpProf->Fill(hpos.at(0), fCheatPeakAmp);
           fCheatSummedIntegralProf->Fill(hpos.at(0), fCheatSumIntegral);
           if(fRescale==true){
             Double_t weightingFactor = 0.0;
             Double_t sumStepWeight = 0.0;
             Int_t binXorigin = rescaleSourceProfile->GetXaxis()->FindBin(0.01*( - fxCorrection ) );
             Int_t binX = rescaleSourceProfile->GetXaxis()->FindBin(0.01*(hpos.at(0) - fxCorrection)); //Factor of 100 hard coded for now. Bad. Profile in meters, Sim in centimeters.
             Int_t binY = rescaleSourceProfile->GetYaxis()->FindBin(0.01*(hpos.at(1) - fyCorrection));
             Int_t binZ = rescaleSourceProfile->GetZaxis()->FindBin(0.01*(hpos.at(2) - fzCorrection + fzOffset));
             while( binX >= binXorigin ){
               //Double_t stepPos = 15*binX+7.5;
               //Double_t stepPos = rescaleSourceProfile->GetXaxis()->GetBinCenter(binX);
               Double_t stepPos = ( (Double_t)(rescaleSourceProfile->GetXaxis()->GetBinCenter(binX)) + fxCorrection );
               Double_t impurityLevel = 
                   rescaleSourceProfile->GetBinContent( binX, binY, binZ ); 
               weightingFactor += exp(-(stepPos)/(fVel*(2.0-impurityLevel)*fTau))/exp(02.0*(stepPos)/(fVel*fTau));          
               sumStepWeight += 1.0/exp(-(stepPos)/(fVel*fTau));
               binX -= 1;
             }
             if(sumStepWeight!=0){weightingFactor /= sumStepWeight;}else{weightingFactor=0;}
             if(fCheatPeakAmp!=0.0){fCheatPeakAmpProf->Fill(hpos.at(0), fCheatPeakAmp*weightingFactor );}
             fCheatSummedIntegralProf->Fill(hpos.at(0), fCheatSumIntegral*weightingFactor );
             if( hpos.at(1)<-555 ){fBCheatSummedIntegralProf->Fill(hpos.at(0), fCheatSumIntegral*weightingFactor);}
             else if( hpos.at(2)<30 && hpos.at(2)>-30 ){fCCheatSummedIntegralProf->Fill(hpos.at(0), fCheatSumIntegral*weightingFactor);}
             else if( hpos.at(2)>555 ) { fTCheatSummedIntegralProf->Fill(hpos.at(0), fCheatSumIntegral*weightingFactor);}

           }else{
             if(fCheatPeakAmp!=0.0){fCheatPeakAmpProf->Fill(hpos.at(0), fCheatPeakAmp);}
             fCheatSummedIntegralProf->Fill(hpos.at(0), fCheatSumIntegral);
           }  
        }

       catch(cet::exception e){
//           std::cout << "ERROR: Backtracker failed. This (dcheat)hit's peak amplitude was: " << fCheatPeakAmp << "\n";

       }

     }  

  } // bnazie::analyze()
  
  
  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see bnazie.fcl for more information.
  DEFINE_ART_MODULE(bnazie)

} // namespace bnazie

#endif // bnazie_Module

