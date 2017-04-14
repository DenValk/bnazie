/*
 *Metrics Analysis program by Jason Stock of South Dakota School of Mines and Technology
 *September 2015.
 *
 *
 * This Version of the PhotonCount10KT analysis uses the backtracker and is built for events with multiple particles. As such, its structure is different than previous analysis used before Dec 2015.
 *
 */

#ifndef PhotonCount10KT_Module
#define PhotonCount10KT_Module

// LArSoft includes
#include "lardata/RecoBase/OpHit.h"
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

} // end local namespace





namespace PhotonCount10KT {

  class PhotonCount10KT : public art::EDAnalyzer
  {
  public:
 
    explicit PhotonCount10KT(fhicl::ParameterSet const& parameterSet);

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

    std::string fOpHitLabel;

    TH1D*     fPEHist;
    
    Double_t fNPE;

    //art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();

  }; // class PhotonCount10KT


  
  PhotonCount10KT::PhotonCount10KT(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
  {
    reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  void PhotonCount10KT::beginJob()
  {
 
    fPEHist                   = tfs->make< TH1D >  ("fPEHist", ";The number of photons detected per OpHit;", 200, 0, 2);

    fPEHist -> SetXTitle("Number of PE detected each OpHit");
    fPEHist -> SetYTitle("");

  }
  

   //-----------------------------------------------------------------------
  void PhotonCount10KT::endJob()
  {

  }
  //-----------------------------------------------------------------------
  void PhotonCount10KT::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    // Read parameters from the .fcl file. The names in the arguments
    fOpHitLabel           = parameterSet.get< std::string >("PhotLabel");
    
  }

  //-----------------------------------------------------------------------
  void PhotonCount10KT::analyze(const art::Event& evt) 
  {

    art::Handle< std::vector< recob::OpHit > > opHitHandle;
    std::vector< art::Ptr < recob::OpHit > > opHitList;
    if (evt.getByLabel(fOpHitLabel, opHitHandle) )
      art::fill_ptr_vector(opHitList, opHitHandle);
    

     for ( std::vector< art::Ptr <recob::OpHit> >::iterator opHitIt = opHitList.begin(); opHitIt != opHitList.end(); ++opHitIt){
         art::Ptr<recob::OpHit> thisHit = *opHitIt;
         double fPE = thisHit->PE();
         fPEHist->Fill(fPE);

     }

  } // PhotonCount10KT::analyze()
  
  
  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see PhotonCount10KT.fcl for more information.
  DEFINE_ART_MODULE(PhotonCount10KT)



} // namespace PhotonCount10KT

#endif // PhotonCount10KT_Module














