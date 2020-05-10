// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      SoftbTableProducer
// 
/**\class SoftbTableProducer SoftbTableProducer.cc PhysicsTools/SoftbTableProducer/plugins/SoftbTableProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrea Rizzi
//         Created:  Mon, 28 Aug 2017 09:26:39 GMT
//
//


// system include files
#include <memory>
#include "GeneratorInterface/Core/interface/GenFilterEfficiencyAnalyzer.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "PhysicsTools/SelectorUtils/interface/RunLumiSelector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DetectorDescription/Core/interface/Singleton.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "PhysicsTools/NanoAOD/interface/MatchingUtils.h"
//
// class declaration
//

class SignalSumGenWeight : public edm::stream::EDProducer<> {
   public:
      explicit SignalSumGenWeight(const edm::ParameterSet&);
      ~SignalSumGenWeight() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginStream(edm::StreamID) override;
      void produce(edm::Event&, const edm::EventSetup&) override;
      void endStream() override;

      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override; 
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      const edm::InputTag genFilterInfoTag_;
      const edm::EDGetTokenT<GenFilterInfo> genFilterInfoToken_;
      const std::string  pvName_;
      edm::Handle<GenFilterInfo> genFilter;

};



//
// constructors and destructor
//

SignalSumGenWeight::SignalSumGenWeight(const edm::ParameterSet& params):
    genFilterInfoTag_(params.getParameter<edm::InputTag>("genFilterInfoTag")),
    genFilterInfoToken_(consumes<GenFilterInfo,edm::InLumi>(params.getParameter<edm::InputTag>("genFilterInfoTag"))),
    pvName_(params.getParameter<std::string>("pvName") )
   
{
   produces<float>("sumPassWeights");
   produces<float>("sumWeights");
   produces<nanoaod::FlatTable>("pv");
   produces<edm::PtrVector<reco::Candidate> >();
}


SignalSumGenWeight::~SignalSumGenWeight()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------


void
SignalSumGenWeight::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    using namespace edm;

    auto GenTable = std::make_unique<nanoaod::FlatTable>(1,pvName_,true,false);
    // For SV we fill from here only stuff that cannot be created with the SimpleFlatTableProducer 


    float sumPassWeights_ = genFilter->sumPassWeights();
    float sumWeights_ = genFilter->sumWeights();
        
     GenTable->addColumnValue<float>("SigWeightGenCutPass", genFilter->sumPassWeights(), "Gen Cut weight on     Signal", nanoaod::FlatTable::FloatColumn,10);
     GenTable->addColumnValue<float>("SigWeightGenCut", genFilter->sumWeights(), "Gen Cut weighton Signa    l", nanoaod::FlatTable::FloatColumn,10);

    iEvent.put(std::move(GenTable),"pv");
}

void SignalSumGenWeight::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup) {
    iLumi.getByLabel(genFilterInfoTag_,genFilter);
}

void SignalSumGenWeight::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup) {
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
SignalSumGenWeight::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
SignalSumGenWeight::endStream() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SignalSumGenWeight::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SignalSumGenWeight);
