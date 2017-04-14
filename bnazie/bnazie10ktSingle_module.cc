/*
 *Metrics Analysis program by Jason Stock of South Dakota School of Mines and Technology
 *June 2016
 *
 *
 * This Version of the bnazie10KTSingle analysis uses the backtracker and is built for events with multiple particles. As such, its structure is different than previous analysis used before Dec 2015.
 *
 */

#ifndef bnazie10KTSingle_Module
#define bnazie10KTSingle_Module

// LArSoft includes
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
  /*int signF (double a){
    if(a==0){ return 0; }
    else if(a<0){return -1;}
    else{ return 1;}
  }*/
  /*double modulo (double a, int b){
    double out =0.0;
    out = signF(a)*fmod(abs(a),((double)abs(b)));
    return out;
  }*/
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





namespace bnazie10KTSingle {

  class bnazie10KTSingle : public art::EDAnalyzer
  {
  public:
 
    explicit bnazie10KTSingle(fhicl::ParameterSet const& parameterSet);

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

    std::string fHitLabel;
    std::string fPhotonsLabel;
    std::string fTruthLabel;

    TH2D* fIntegral2Hist;
    TH2D* fPhotoElec2Hist;

    TH1D* fIntegral1DHist1;
    TH1D* fIntegral1DHist2;
    TH1D* fIntegral1DHist3;
    TH1D* fIntegral1DHist4;
    TH1D* fIntegral1DHist5;
    TH1D* fIntegral1DHist6;
    TH1D* fIntegral1DHist7;
    TH1D* fIntegral1DHist8;
    TH1D* fIntegral1DHist9;
    TH1D* fIntegral1DHist10;
    TH1D* fIntegral1DHist11;
    TH1D* fIntegral1DHist12;
    TH1D* fIntegral1DHist13;
    TH1D* fIntegral1DHist14;
    TH1D* fIntegral1DHist15;
    TH1D* fIntegral1DHist16;
    TH1D* fIntegral1DHist17;
    TH1D* fIntegral1DHist18;
    TH1D* fIntegral1DHist19;
    TH1D* fIntegral1DHist20;

    Double_t fSummedInt;
    Double_t fNPE;

    //detector boundaries

    //art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();

  }; // class bnazie10KTSingle


  
  bnazie10KTSingle::bnazie10KTSingle(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  void bnazie10KTSingle::beginJob()
  {
 
  }
  

   //-----------------------------------------------------------------------
  void bnazie10KTSingle::endJob()
  {

  }
  //-----------------------------------------------------------------------
  void bnazie10KTSingle::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    fHitLabel               = parameterSet.get< std::string >("HitLabel");
    fPhotonsLabel           = parameterSet.get< std::string >("PhotLabel");
    fTruthLabel             = parameterSet.get< std::string >("TruthLabel");
    fSummedInt  = 0.0;
    fNPE        = 0.0;


    fIntegral2Hist             = tfs->make<TH2D>("integralHist", ";Historgram of IntegralCharge;", 21, 122.5, 227.5, 20, -20.5, -0.5 );
    fPhotoElec2Hist            = tfs->make<TH2D>("photoEHist",   ";Historgram of PhotoElectrons;", 21, 122.5, 227.5, 20, -20.5, -0.5 );

    fIntegral1DHist1           = tfs->make<TH1D>("int1", "Histogram of Integrated charge vs X distance. z=-1;", 21, 122.5, 227.5);
    fIntegral1DHist2           = tfs->make<TH1D>("int2", "Histogram of Integrated charge vs X distance; z=-2", 21, 122.5, 227.5);
    fIntegral1DHist3           = tfs->make<TH1D>("int3", "Histogram of Integrated charge vs X distance; z=-3", 21, 122.5, 227.5);
    fIntegral1DHist4           = tfs->make<TH1D>("int4", "Histogram of Integrated charge vs X distance; z=-4", 21, 122.5, 227.5);
    fIntegral1DHist5           = tfs->make<TH1D>("int5", "Histogram of Integrated charge vs X distance; z=-5", 21, 122.5, 227.5);
    fIntegral1DHist6           = tfs->make<TH1D>("int6", "Histogram of Integrated charge vs X distance; z=-6", 21, 122.5, 227.5);
    fIntegral1DHist7           = tfs->make<TH1D>("int7", "Histogram of Integrated charge vs X distance; z=-7", 21, 122.5, 227.5);
    fIntegral1DHist8           = tfs->make<TH1D>("int8", "Histogram of Integrated charge vs X distance; z=-8", 21, 122.5, 227.5);
    fIntegral1DHist9           = tfs->make<TH1D>("int9", "Histogram of Integrated charge vs X distance; z=-9", 21, 122.5, 227.5);
    fIntegral1DHist10          = tfs->make<TH1D>("int10", "Histogram of Integrated charge vs X distance; z=-10", 21, 122.5, 227.5);
    fIntegral1DHist11          = tfs->make<TH1D>("int11", "Histogram of Integrated charge vs X distance; z=-11", 21, 122.5, 227.5);
    fIntegral1DHist12          = tfs->make<TH1D>("int12", "Histogram of Integrated charge vs X distance; z=-12", 21, 122.5, 227.5);
    fIntegral1DHist13          = tfs->make<TH1D>("int13", "Histogram of Integrated charge vs X distance; z=-13", 21, 122.5, 227.5);
    fIntegral1DHist14          = tfs->make<TH1D>("int14", "Histogram of Integrated charge vs X distance; z=-14", 21, 122.5, 227.5);
    fIntegral1DHist15          = tfs->make<TH1D>("int15", "Histogram of Integrated charge vs X distance; z=-15", 21, 122.5, 227.5);
    fIntegral1DHist16          = tfs->make<TH1D>("int16", "Histogram of Integrated charge vs X distance; z=-16", 21, 122.5, 227.5);
    fIntegral1DHist17          = tfs->make<TH1D>("int17", "Histogram of Integrated charge vs X distance; z=-17", 21, 122.5, 227.5);
    fIntegral1DHist18          = tfs->make<TH1D>("int18", "Histogram of Integrated charge vs X distance; z=-18", 21, 122.5, 227.5);
    fIntegral1DHist19          = tfs->make<TH1D>("int19", "Histogram of Integrated charge vs X distance; z=-19", 21, 122.5, 227.5);
    fIntegral1DHist20          = tfs->make<TH1D>("int20", "Histogram of Integrated charge vs X distance; z=-20", 21, 122.5, 227.5);
    
    
  }

  //-----------------------------------------------------------------------
  void bnazie10KTSingle::analyze(const art::Event& evt) 
  {

    TLorentzVector posVec;	  
    fNPE       = 0.0;
    fSummedInt = 0.0;

    art::Handle< std::vector< simb::MCParticle > > truthHandle;
    std::vector< art::Ptr < simb::MCParticle > > truthlist;
    if (evt.getByLabel(fTruthLabel, truthHandle) )
      art::fill_ptr_vector(truthlist, truthHandle);

    art::Handle< std::vector< recob::Hit > > hitHandle;
    std::vector< art::Ptr< recob::Hit > > hitlist;
    if (evt.getByLabel(fHitLabel, hitHandle) )
      art::fill_ptr_vector(hitlist, hitHandle);

    art::Handle< std::vector< recob::OpHit > > ophitHandle;
    std::vector< art::Ptr< recob::OpHit > > ophitlist;
    if (evt.getByLabel(fPhotonsLabel, ophitHandle) )
      art::fill_ptr_vector(ophitlist, ophitHandle);

    posVec=(*(*truthlist.begin())).Position();

     for ( std::vector< art::Ptr <recob::OpHit> >::iterator opHitIt = ophitlist.begin(); opHitIt != ophitlist.end(); ++opHitIt){
        art::Ptr<recob::OpHit> thisHit = *opHitIt;
        fNPE = thisHit->PE();
        if(fNPE!=0){fPhotoElec2Hist ->Fill(posVec.X(), posVec.Z(), fNPE);}

     }

     for ( std::vector< art::Ptr < recob::Hit > >::iterator hitIt = hitlist.begin(); hitIt != hitlist.end(); ++hitIt){
        art::Ptr<recob::Hit> thisHit = *hitIt;
        fSummedInt    = thisHit->Integral();
        std::cout<<"X"<<posVec.X()<<"Z"<<posVec.Z()<<"INT"<<fSummedInt<<"\n";
        if(fSummedInt!=0){fIntegral2Hist  ->Fill(posVec.X(), posVec.Z(), fSummedInt);}
        
        if( posVec.Z()>0 && posVec.Z()<2 ){
                fIntegral1DHist1->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>1 && posVec.Z()<3 ){
                fIntegral1DHist2->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>2 && posVec.Z()<4 ){
                fIntegral1DHist3->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>3 && posVec.Z()<5 ){
                fIntegral1DHist4->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>4 && posVec.Z()<6 ){
                fIntegral1DHist5->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>5 && posVec.Z()<7 ){
                fIntegral1DHist6->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>6 && posVec.Z()<8 ){
                fIntegral1DHist7->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>7 && posVec.Z()<9 ){
                fIntegral1DHist8->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>8 && posVec.Z()<10 ){
                fIntegral1DHist9->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>9 && posVec.Z()<11 ){
                fIntegral1DHist10->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>10 && posVec.Z()<12 ){
                fIntegral1DHist11->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>11 && posVec.Z()<13 ){
                fIntegral1DHist12->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>12 && posVec.Z()<14 ){
                fIntegral1DHist13->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>13 && posVec.Z()<15 ){
                fIntegral1DHist14->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>14 && posVec.Z()<16 ){
                fIntegral1DHist15->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>15 && posVec.Z()<17 ){
                fIntegral1DHist16->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>16 && posVec.Z()<18 ){
                fIntegral1DHist17->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>17 && posVec.Z()<19 ){
                fIntegral1DHist18->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>18 && posVec.Z()<20 ){
                fIntegral1DHist19->Fill(posVec.X(), fSummedInt);
        }else if( posVec.Z()>19 && posVec.Z()<21 ){
                fIntegral1DHist20->Fill(posVec.X(), fSummedInt);
        }
     } 
    

  } // bnazie10KTSingle::analyze()
  
  
  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see bnazie10KTSingle.fcl for more information.
  DEFINE_ART_MODULE(bnazie10KTSingle)



} // namespace bnazie10KTSingle

#endif // bnazie10KTSingle_Module














