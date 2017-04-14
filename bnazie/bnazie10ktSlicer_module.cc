/*
 *Metrics Analysis program by Jason Stock of South Dakota School of Mines and Technology
 *September 2015.
 *
 *
 * This Version of the bnazie10KTSlicer analysis uses the backtracker and is built for events with multiple particles. As such, its structure is different than previous analysis used before Dec 2015.
 *
 */

#ifndef bnazie10KTSlicer_Module
#define bnazie10KTSlicer_Module

// LArSoft includes
#include "larsim/MCCheater/BackTracker.h"
//#include "Simulation/LArG4Parameters.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/OpHit.h"
#include "lardata/RecoBase/OpFlash.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
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
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Utilities/Exception.h"
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


/////////////////////////////ANONYMOUS NAMESPACE
namespace {

  /*int signF (int a){
    if(a==0){ return 0;}
    else if(a<0){ return -1; }
    else{ return 1;
    }
  }*/
  int signF (double a){
    if(a==0){ return 0; }
    else if(a<0){return -1;}
    else{ return 1;}
  }
  double modulo (double a, int b){
    double out =0.0;
    out = signF(a)*fmod(abs(a),((double)abs(b)));
    return out;
  }
  /*double modulo (double a, double b){
    double out =0.0;
    out = signF(a)*fmod(abs(a),abs(b));
    return out;
  }*/
  /*int modulo (int a, int b){
    int out;
    out = signF(a)*fmod(abs(a),abs(b));
    return out;
  }*/


} // end local namespace





namespace bnazie10KTSlicer {

  class bnazie10KTSlicer : public art::EDAnalyzer
  {
  public:
 
    explicit bnazie10KTSlicer(fhicl::ParameterSet const& parameterSet);

    virtual void beginJob();
    virtual void reconfigure(fhicl::ParameterSet const& parameterSet) override;
    virtual void analyze (const art::Event& event) override;
    virtual void endJob();
	

  private:

/*    double modulo (double a, double b);
    double modulo (double a, int b);
    int modulo (int a, int b);
    int signF (int a);
    int signF (double a);*/

    std::string fPositionsLabel;
    std::string fHitLabel;
    std::string fCheatLabel;
    std::string fPhotonsLabel;
    std::string fRescaleFile;
    std::string fRescaleObj;

    TFile* rescaleSourceFile;

    TTree* fMCData;  
    TTree* fMetricsNtuple;  
    TTree* fCheatMetricsNtuple;  

    TProfile* barrelPeak;
//    TProfile* barrelPhot;

    TProfile* yn720_to_n660Hist;
    TProfile* yn660_to_n600Hist;
    TProfile* yn600_to_n540Hist;
    TProfile* yn540_to_n480Hist;
    TProfile* yn480_to_n420Hist;
    TProfile* yn420_to_n360Hist;
    TProfile* yn360_to_n300Hist;
    TProfile* yn300_to_n240Hist;
    TProfile* yn240_to_n180Hist;
    TProfile* yn180_to_n120Hist;
    TProfile* yn120_to_n060Hist;
    TProfile* yn060_to_n000Hist;
    TProfile* y000_to_060Hist;
    TProfile* y060_to_120Hist;
    TProfile* y120_to_180Hist;
    TProfile* y180_to_240Hist;
    TProfile* y240_to_300Hist;
    TProfile* y300_to_360Hist;
    TProfile* y360_to_420Hist;
    TProfile* y420_to_480Hist;
    TProfile* y480_to_540Hist;
    TProfile* y540_to_600Hist;
    TProfile* y600_to_660Hist;
    TProfile* y660_to_720Hist;


    TProfile* fFullVolumeSummedIntegralHist; 


    TProfile* fSummedIntegralHist;
    TProfile* fCheatSummedIntegralHist;
    TProfile* fPeakAmpHist;
    TProfile* fCheatPeakAmpHist;
    TProfile* fPEvXHist;
    TProfile3D* rescaleSourceProfile;

    TCanvas*  can6plot;

    TH1D*     fPEHist;
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
    Int_t xBs, yBs, zBs;

    Double_t xPos, fPeakAmp, fMaxPeakAmp, fSumIntegral, totPeakAmp, totCheatPeakAmp, totIntegral, fCheatPeakAmp, fCheatMaxPeakAmp;
    Double_t fCheatSumIntegral, fxCorrection, fyCorrection, fzCorrection, fzOffset, fTau, fVel, asymAmp;

    //detector boundaries

    Double_t xLo, xHi, yLo, yHi, zLo, zHi;

    Bool_t fRescale;

    //art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<cheat::BackTracker> bt;
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();

  }; // class bnazie10KTSlicer


  
  bnazie10KTSlicer::bnazie10KTSlicer(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  void bnazie10KTSlicer::beginJob()
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



    fFullVolumeSummedIntegralHist       = tfs->make<TProfile>("fullVolumeSummedIntHist",      ";Profile of the Summed Integral Charge (DUC) for each event vs Distance from APA (cm);",           xBs, xLo, xHi);

    barrelPeak             = tfs->make<TProfile>("barrelPeakHist", ";Selected Pipe Near Field Cage. 20x20cm PeakAmp;", xBs, xLo, xHi );
//    barrelPhot             = tfs->make<TProfile>("barrelPhotHist", ";Selected Pipe Near Field Cage. 20x20cm Photons;", xBs, xLo, xHi );

    yn720_to_n660Hist  = tfs->make<TProfile>("n720y",  ";-660y to -600y  PeakAmp;",  xBs, xLo, xHi );
    yn660_to_n600Hist  = tfs->make<TProfile>("n660y",  ";-600y to -540y  PeakAmp;",  xBs, xLo, xHi );
    yn600_to_n540Hist  = tfs->make<TProfile>("n600y",  ";-600y to -540y  PeakAmp;",  xBs, xLo, xHi );
    yn540_to_n480Hist  = tfs->make<TProfile>("n540y",  ";-540y to -480y  PeakAmp;",  xBs, xLo, xHi );
    yn480_to_n420Hist  = tfs->make<TProfile>("n480y",  ";-480y to -420y  PeakAmp;",  xBs, xLo, xHi );
    yn060_to_n000Hist  = tfs->make<TProfile>("n420y",  ";-420y to -360y  PeakAmp;",  xBs, xLo, xHi );
    yn420_to_n360Hist  = tfs->make<TProfile>("n360y",  ";-360y to -300y  PeakAmp;",  xBs, xLo, xHi );
    yn360_to_n300Hist  = tfs->make<TProfile>("n300y",  ";-300y to -240y  PeakAmp;",  xBs, xLo, xHi );
    yn300_to_n240Hist  = tfs->make<TProfile>("n240y",  ";-240y to -180y  PeakAmp;",  xBs, xLo, xHi );
    yn240_to_n180Hist  = tfs->make<TProfile>("n180y",  ";-180y to -120y  PeakAmp;",  xBs, xLo, xHi );
    yn180_to_n120Hist  = tfs->make<TProfile>("n120y",  ";-120y to -060y  PeakAmp;",  xBs, xLo, xHi );
    yn120_to_n060Hist  = tfs->make<TProfile>("n060y",  ";-060y to -000y  PeakAmp;",  xBs, xLo, xHi );
    y000_to_060Hist    = tfs->make<TProfile>("p000y",  ";000y to 060y  PeakAmp;",    xBs, xLo, xHi );
    y060_to_120Hist    = tfs->make<TProfile>("p060y",  ";060y to 120y  PeakAmp;",    xBs, xLo, xHi );
    y120_to_180Hist    = tfs->make<TProfile>("p120y",  ";120y to 180y  PeakAmp;",    xBs, xLo, xHi );
    y180_to_240Hist    = tfs->make<TProfile>("p180y",  ";180y to 240y  PeakAmp;",    xBs, xLo, xHi );
    y240_to_300Hist    = tfs->make<TProfile>("p240y",  ";240y to 300y  PeakAmp;",    xBs, xLo, xHi );
    y300_to_360Hist    = tfs->make<TProfile>("p300y",  ";300y to 360y  PeakAmp;",    xBs, xLo, xHi );
    y360_to_420Hist    = tfs->make<TProfile>("p360y",  ";360y to 420y  PeakAmp;",    xBs, xLo, xHi );
    y420_to_480Hist    = tfs->make<TProfile>("p420y",  ";420y to 480y  PeakAmp;",    xBs, xLo, xHi );
    y480_to_540Hist    = tfs->make<TProfile>("p480y",  ";480y to 540y  PeakAmp;",    xBs, xLo, xHi );
    y540_to_600Hist    = tfs->make<TProfile>("p540y",  ";540y to 600y  PeakAmp;",    xBs, xLo, xHi );
    y600_to_660Hist    = tfs->make<TProfile>("p600y",  ";600y to 660y  PeakAmp;",    xBs, xLo, xHi );
    y660_to_720Hist    = tfs->make<TProfile>("p660y",  ";660y to 720y  PeakAmp;",    xBs, xLo, xHi );

    fPeakAmpHist              = tfs->make<TProfile>("peakAmpHist",        ";Profile of Peak Amplitudes for each hit vs BackTracker Distance from APA (cm);",                  xBs, xLo, xHi);
    fCheatPeakAmpHist         = tfs->make<TProfile>("cheatPeakAmpHist",   ";Profile of Cheated Peak Amplitudes for each hit vs BackTracker Distance from APA (cm);",          xBs, xLo, xHi);
    fSummedIntegralHist       = tfs->make<TProfile>("summedIntHist",      ";Profile of the Summed Integral Charge (DUC) for each event vs Distance from APA (cm);",           xBs, xLo, xHi);
    fCheatSummedIntegralHist  = tfs->make<TProfile>("cheatSummedIntHist", ";Profile of the Cheated Summed Integral Charge (DUC) for each event vs Distance from APA (cm);",   xBs, xLo, xHi);
    fPEvXHist                    = tfs->make<TProfile>("PEvXHist",              ";Profile of the Number of Photo Electrons detected per event vs Distance from APA;",          xBs, xLo, xHi);
    hitsDetected              = tfs->make< TH3D >  ("hitsDetectedPlot",   ";TH3 of the weighted average position for all Detected Hits (As given in Backtracker) ;",          xBs, xLo, xHi, yBs, yLo, yHi, zBs, zLo, zHi);
    fNumHitsHist              = tfs->make< TH1D >  ("numHits", ";The number of hits detected at a distance x from the APA;", xBs, xLo, xHi);
    fCheatNumHitsHist         = tfs->make< TH1D >  ("cheatNumHits", ";The number of cheated hits detected at a distance x from the APA;", xBs, xLo, xHi);
    xHits                     = tfs->make< TH1D >  ("xHits", ";The number of hits at different positions X;", xBs, xLo, xHi);
    yHits                     = tfs->make< TH1D >  ("yHits", ";The number of hits at different positions Y;", yBs, yLo, yHi);
    zHits                     = tfs->make< TH1D >  ("zHits", ";The number of hits ad different positions Z:", zBs, zLo, zHi);

    if(fRescale==0){can6plot = tfs->make<TCanvas>("VersionValidationPlots", "vv6plotspread", 2550, 3300);}

    ///////////////////////SET HISTOGRAM LABELS/////////////////
    fPeakAmpHist->SetXTitle("Distance from APA (cm)");
    fPeakAmpHist->SetYTitle("Peak Amplitude");

    fCheatPeakAmpHist->SetXTitle("Distance from APA (cm)");
    fCheatPeakAmpHist->SetYTitle("Cheated Peak Amplitude");

    fSummedIntegralHist->SetXTitle("Distance from APA (cm)");
    fSummedIntegralHist->SetYTitle("Summed Integral Charge per Event (Deposited Units Charge)");

    fCheatSummedIntegralHist->SetXTitle("Distance from APA (cm)");
    fCheatSummedIntegralHist->SetYTitle("Cheated Summed Integral Charge per Event (Deposited Units Charge)");

    fNumHitsHist->SetXTitle("Distance from APA (cm)");
    fNumHitsHist->SetYTitle("Number of hits detected.");

    fCheatNumHitsHist->SetXTitle("Distance from APA (cm)");
    fCheatNumHitsHist->SetYTitle("Number of Cheated Hits Detected.");

    hitsDetected->SetXTitle("X Position ( cm )");
    hitsDetected->SetYTitle("Y Position ( cm )");
    hitsDetected->SetZTitle("Z Position ( cm )");
  }
  

   //-----------------------------------------------------------------------
  void bnazie10KTSlicer::endJob()
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
        fPeakAmpHist->Draw();
        can6plot -> cd(2);    
        fCheatPeakAmpHist->Draw();
        can6plot -> cd(3);
        fSummedIntegralHist->Draw();
        can6plot -> cd(4);
        fCheatSummedIntegralHist->Draw();
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
  void bnazie10KTSlicer::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fPositionsLabel         = parameterSet.get< std::string >("TruthLabel");
    fHitLabel               = parameterSet.get< std::string >("HitLabel");
    fCheatLabel             = parameterSet.get< std::string >("CheatLabel");
    fPhotonsLabel           = parameterSet.get< std::string >("PhotLabel");
    fRescale                = parameterSet.get< Bool_t      >("RescaleOn");
    asymAmp                 = parameterSet.get< Double_t    >("AsymetryAmplitude");
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
    xBs=0;
    yBs=0;
    zBs=0;
    xLo=0.0;
    xHi=0.0;
    yLo=0.0;
    yHi=0.0;
    zLo=0.0;
    zHi=0.0;
    //Set boundaries
    {
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
      Double_t minDist = std::min( std::min( xHi-xLo, yHi-yLo ), zHi-zLo );
      Int_t divDist = ((int)(minDist))/100;
      xLo = ( xLo - modulo(xLo, divDist) )  - divDist;
      xHi = ( xHi - modulo(xHi, divDist) )  + divDist;
      xBs = ( xHi - xLo )/divDist;
      yLo = ( yLo - modulo(yLo, divDist) )  - divDist;
      yHi = ( yHi - modulo(yHi, divDist) )  + divDist;
      yBs = ( yHi - yLo )/divDist;
      zLo = ( zLo - modulo(zLo, divDist) )  - divDist;
      zHi = ( zHi - modulo(zHi, divDist) )  + divDist;
      zBs = ( zHi - zLo )/divDist;

    }
    
  }

  //-----------------------------------------------------------------------
  void bnazie10KTSlicer::analyze(const art::Event& evt) 
  {

    TLorentzVector posVec;	  
    Double_t nPE = 0.0;
    Double_t aPos = 0.0;
    Int_t tNumHits = 0;

	  
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

    for ( std::vector< art::Ptr <recob::OpHit> >::iterator opHitIt = ophitlist.begin(); opHitIt != ophitlist.end(); ++opHitIt){
      nPE += (*(*opHitIt)).PE();
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
           aPos += hpos.at(0);
           fNumHitsHist->Fill(hpos.at(0));
           hitsDetected->Fill(hpos.at(0), hpos.at(1), hpos.at(2));
           if( fRescale==true ){
             Double_t weightingFactor = 0.0;
             Double_t sumStepWeight = 0.0;
             Int_t binXorigin = rescaleSourceProfile->GetXaxis()->FindBin(0.01*( 0.0 - fxCorrection ) );
             Int_t binX = rescaleSourceProfile->GetXaxis()->FindBin(0.01*(hpos.at(0) - fxCorrection)); //Factor of 100 hard coded for now. Bad. Profile in meters, Sim in centimeters.
             Int_t binY = rescaleSourceProfile->GetYaxis()->FindBin(0.01*(hpos.at(1) - fyCorrection));
             Int_t binZ = rescaleSourceProfile->GetZaxis()->FindBin(0.01*(hpos.at(2) - fzCorrection + fzOffset));
             
             if (binX > binXorigin){
               while( binX >= binXorigin ){  //Currently only works for hits in the x>0 region.
                 //STEP POS!!  Get bin center (in meters), and add the correction (in centimeters). Final value in centimeters.
                 Double_t stepPos = ( 100 * ( (Double_t)(rescaleSourceProfile->GetXaxis()->GetBinCenter(binX)) + (0.01*fxCorrection) ) ) ;
                 Double_t impurityLevel = rescaleSourceProfile->GetBinContent( binX, binY, binZ ); 
                 impurityLevel = ((impurityLevel - 1.0) * asymAmp) + 1.0;
                 weightingFactor += exp(-(stepPos) / (fVel * (2.0 - impurityLevel) * fTau))/exp(-02.0 * (stepPos)/(fVel * fTau));          
                 sumStepWeight += 1.0/exp(-(stepPos)/(fVel*fTau));
                 binX -= 1;
               }
             }else if(binX < binXorigin){
               while( binX <= binXorigin ){  //Currently only works for hits in the x>0 region.
                 //STEP POS!!  Get bin center (in meters), and add the correction (in centimeters). Final value in centimeters.
                 Double_t stepPos = ( 100 * ( (Double_t)(rescaleSourceProfile->GetXaxis()->GetBinCenter(binX)) + (0.01*fxCorrection) ) ) ;
                 Double_t impurityLevel = rescaleSourceProfile->GetBinContent( binX, binY, binZ ); 
                 impurityLevel = ((impurityLevel - 1.0) * asymAmp) + 1.0;
                 weightingFactor += exp(-(stepPos) / (fVel * (2.0 - impurityLevel) * fTau))/exp(-02.0 * (stepPos)/(fVel * fTau));          
                 sumStepWeight += 1.0/exp(-(stepPos)/(fVel*fTau));
                 binX += 1;
               } 
             }

             if(sumStepWeight!=0){weightingFactor /= sumStepWeight;}else{weightingFactor=1;}
             
             //APPLY WEIGHTING FACTOR
             fPeakAmp = fPeakAmp * weightingFactor;
             fSumIntegral = fSumIntegral*weightingFactor;
             
             if(fPeakAmp!=0.0){
                //std::cout<<"Weighting Factor: " << weightingFactor << "\n";
                fPeakAmpHist->Fill(hpos.at(0), fPeakAmp );
             }
             fSummedIntegralHist->Fill(hpos.at(0), fSumIntegral );
             if( (hpos.at(1) < (yHi ) ) && (hpos.at(1)>yLo ) ){
                 fFullVolumeSummedIntegralHist->Fill(hpos.at(0),fSumIntegral);
             }
             //////////////////////////////////// ManySlicers //////////////////////////////////////////
             
             if      ( (hpos.at(1) > -10) && (hpos.at(1) < 10) && (hpos.at(2) > 0) && (hpos.at(2) < 20) ){barrelPeak->Fill(hpos.at(0), fPeakAmp);}
//             if      (hpos.at(1) > -10 && hpos.at(1) < 10 && hpos.at(2) > 0 && hpos.at(2)<20){barrelPhot ->Fill(hpos.at(0), fPeakAmp);}
             
             if (hpos.at(1) >= -720 && hpos.at(1) < -660 ){yn720_to_n660Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -660 && hpos.at(1) < -600 ){yn660_to_n600Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -600 && hpos.at(1) < -540 ){yn600_to_n540Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -540 && hpos.at(1) < -480 ){yn540_to_n480Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -480 && hpos.at(1) < -420 ){yn480_to_n420Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -420 && hpos.at(1) < -360 ){yn420_to_n360Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -360 && hpos.at(1) < -300 ){yn360_to_n300Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -300 && hpos.at(1) < -240 ){yn300_to_n240Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -240 && hpos.at(1) < -180 ){yn240_to_n180Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -180 && hpos.at(1) < -120 ){yn180_to_n120Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -120 && hpos.at(1) < - 60 ){yn120_to_n060Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= - 60 && hpos.at(1) <    0 ){yn060_to_n000Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=    0 && hpos.at(1) <   60 ){y000_to_060Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=   60 && hpos.at(1) <  120 ){y060_to_120Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  120 && hpos.at(1) <  180 ){y120_to_180Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  180 && hpos.at(1) <  240 ){y180_to_240Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  240 && hpos.at(1) <  300 ){y240_to_300Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  300 && hpos.at(1) <  360 ){y300_to_360Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  360 && hpos.at(1) <  420 ){y360_to_420Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  420 && hpos.at(1) <  480 ){y420_to_480Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  480 && hpos.at(1) <  540 ){y480_to_540Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  480 && hpos.at(1) <  540 ){y480_to_540Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  540 && hpos.at(1) <  600 ){y540_to_600Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  600 && hpos.at(1) <  660 ){y600_to_660Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  660 && hpos.at(1) <  720 ){y660_to_720Hist    ->Fill(hpos.at(0), fPeakAmp);}
             //else if (hpos.at(1) >=  540 && hpos.at(1) <= 600 ){y540_to_600Hist  ->Fill(hpos.at(0), fPeakAmp);}

             //////////////////////////////////////////////////////////////////////////////////////////

           }else{
             if(fPeakAmp!=0.0){fPeakAmpHist->Fill(hpos.at(0), fPeakAmp);}
             fSummedIntegralHist->Fill(hpos.at(0), fSumIntegral);

             if      ( (hpos.at(1) > -10) && (hpos.at(1) < 10) && (hpos.at(2) > 0) && (hpos.at(2) < 20) ){barrelPeak->Fill(hpos.at(0), fPeakAmp);}

             if (hpos.at(1) >= -720 && hpos.at(1) < -660 ){yn720_to_n660Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -660 && hpos.at(1) < -600 ){yn660_to_n600Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -600 && hpos.at(1) < -540 ){yn600_to_n540Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -540 && hpos.at(1) < -480 ){yn540_to_n480Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -480 && hpos.at(1) < -420 ){yn480_to_n420Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -420 && hpos.at(1) < -360 ){yn420_to_n360Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -360 && hpos.at(1) < -300 ){yn360_to_n300Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -300 && hpos.at(1) < -240 ){yn300_to_n240Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -240 && hpos.at(1) < -180 ){yn240_to_n180Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -180 && hpos.at(1) < -120 ){yn180_to_n120Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= -120 && hpos.at(1) < - 60 ){yn120_to_n060Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >= - 60 && hpos.at(1) <    0 ){yn060_to_n000Hist  ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=    0 && hpos.at(1) <   60 ){y000_to_060Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=   60 && hpos.at(1) <  120 ){y060_to_120Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  120 && hpos.at(1) <  180 ){y120_to_180Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  180 && hpos.at(1) <  240 ){y180_to_240Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  240 && hpos.at(1) <  300 ){y240_to_300Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  300 && hpos.at(1) <  360 ){y300_to_360Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  360 && hpos.at(1) <  420 ){y360_to_420Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  420 && hpos.at(1) <  480 ){y420_to_480Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  480 && hpos.at(1) <  540 ){y480_to_540Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  480 && hpos.at(1) <  540 ){y480_to_540Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  540 && hpos.at(1) <  600 ){y540_to_600Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  600 && hpos.at(1) <  660 ){y600_to_660Hist    ->Fill(hpos.at(0), fPeakAmp);}
             else if (hpos.at(1) >=  660 && hpos.at(1) <  720 ){y660_to_720Hist    ->Fill(hpos.at(0), fPeakAmp);}


           }  

           fMetricsNtuple->Fill();
           { // The filling of XYZ Histograms (Separate function for easy removal.
             xHits->Fill(hpos.at(0));
             yHits->Fill(hpos.at(1));
             zHits->Fill(hpos.at(2));
           }
           numHits++;
           tNumHits++;
       }
       catch(cet::exception e){
//           std::cout << "ERROR: Backtracker failed. This (hitfd)hit's peak amplitude was: " << fCheatPeakAmp << "\n";
           faultHits++;
       }
     }
	 
     ///////////Plot number of PEs
     if( fRescale==false ){
         aPos = aPos / ( (double)tNumHits );
         fPEvXHist->Fill(aPos, nPE);
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
           if(fRescale==true){
             Double_t weightingFactor = 0.0;
             Double_t sumStepWeight = 0.0;
             Int_t binXorigin = rescaleSourceProfile->GetXaxis()->FindBin(0.01*( - fxCorrection ) );
             Int_t binX = rescaleSourceProfile->GetXaxis()->FindBin(0.01*(hpos.at(0) - fxCorrection)); //Factor of 100 hard coded for now. Bad. Profile in meters, Sim in centimeters.
             Int_t binY = rescaleSourceProfile->GetYaxis()->FindBin(0.01*(hpos.at(1) - fyCorrection));
             Int_t binZ = rescaleSourceProfile->GetZaxis()->FindBin(0.01*(hpos.at(2) - fzCorrection + fzOffset));

             if (binX > binXorigin){
               while( binX >= binXorigin ){  //Currently only works for hits in the x>0 region.
                 //STEP POS!!  Get bin center (in meters), and add the correction (in centimeters). Final value in centimeters.
                 Double_t stepPos = ( 100 * ( (Double_t)(rescaleSourceProfile->GetXaxis()->GetBinCenter(binX)) + (0.01*fxCorrection) ) ) ;
                 Double_t impurityLevel = rescaleSourceProfile->GetBinContent( binX, binY, binZ ); 
                 impurityLevel = ((impurityLevel - 1.0) * asymAmp) + 1.0;
                 weightingFactor += exp(-(stepPos) / (fVel * (2.0 - impurityLevel) * fTau))/exp(-02.0 * (stepPos)/(fVel * fTau));          
                 sumStepWeight += 1.0/exp(-(stepPos)/(fVel*fTau));
                 binX -= 1;
               }
             }else if(binX < binXorigin){
               while( binX <= binXorigin ){  //Currently only works for hits in the x>0 region.
                 //STEP POS!!  Get bin center (in meters), and add the correction (in centimeters). Final value in centimeters.
                 Double_t stepPos = ( 100 * ( (Double_t)(rescaleSourceProfile->GetXaxis()->GetBinCenter(binX)) + (0.01*fxCorrection) ) ) ;
                 Double_t impurityLevel = rescaleSourceProfile->GetBinContent( binX, binY, binZ ); 
                 impurityLevel = ((impurityLevel - 1.0) * asymAmp) + 1.0;
                 weightingFactor += exp(-(stepPos) / (fVel * (2.0 - impurityLevel) * fTau))/exp(-02.0 * (stepPos)/(fVel * fTau));          
                 sumStepWeight += 1.0/exp(-(stepPos)/(fVel*fTau));
                 binX += 1;
               } 
             }

             if(sumStepWeight!=0){weightingFactor /= sumStepWeight;}else{weightingFactor=0;}

             fCheatPeakAmp = fCheatPeakAmp*weightingFactor;
             fCheatSumIntegral = fCheatSumIntegral*weightingFactor;

             if(fCheatPeakAmp!=0.0){
                 fCheatPeakAmpHist->Fill(hpos.at(0), fCheatPeakAmp );
             }
             fCheatSummedIntegralHist->Fill(hpos.at(0), fCheatSumIntegral );

           }else{
             if(fCheatPeakAmp!=0.0){fCheatPeakAmpHist->Fill(hpos.at(0), fCheatPeakAmp);}
             fCheatSummedIntegralHist->Fill(hpos.at(0), fCheatSumIntegral);
           }  
        }

       catch(cet::exception e){
//           std::cout << "ERROR: Backtracker failed. This (dcheat)hit's peak amplitude was: " << fCheatPeakAmp << "\n";

       }

     }  

  } // bnazie10KTSlicer::analyze()
  
  
  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see bnazie10KTSlicer.fcl for more information.
  DEFINE_ART_MODULE(bnazie10KTSlicer)



} // namespace bnazie10KTSlicer

#endif // bnazie10KTSlicer_Module














