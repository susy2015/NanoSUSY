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

class SoftbTableProducer : public edm::stream::EDProducer<> {
   public:
      explicit SoftbTableProducer(const edm::ParameterSet&);
      ~SoftbTableProducer() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginStream(edm::StreamID) override;
      void produce(edm::Event&, const edm::EventSetup&) override;
      void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

      const edm::EDGetTokenT<std::vector<reco::Vertex>> pvs_;
      const edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate> > svs_;
      const edm::EDGetTokenT<edm::View<pat::Jet>> srcJet_;
      const StringCutObjectSelector<reco::Candidate> svCut_;
      const StringCutObjectSelector<reco::Vertex> goodPvCut_;
      const std::string goodPvCutString_;
      const std::string  svName_;
      const std::string svDoc_;
      const double dlenMin_,dlenSigMin_;

};



//
// constructors and destructor
//
SoftbTableProducer::SoftbTableProducer(const edm::ParameterSet& params):
    pvs_(consumes<std::vector<reco::Vertex>>( params.getParameter<edm::InputTag>("pvSrc") )),
    svs_(consumes<edm::View<reco::VertexCompositePtrCandidate> >( params.getParameter<edm::InputTag>("svSrc") )),
    srcJet_(consumes<edm::View<pat::Jet>>(params.getParameter<edm::InputTag>("Jetsrc"))),
    svCut_(params.getParameter<std::string>("svCut") , true),
    goodPvCut_(params.getParameter<std::string>("goodPvCut") , true),
    goodPvCutString_(params.getParameter<std::string>("goodPvCut") ),
    svName_(params.getParameter<std::string>("svName") ),
    svDoc_(params.getParameter<std::string>("svDoc") ),
    dlenMin_(params.getParameter<double>("dlenMin") ),
    dlenSigMin_(params.getParameter<double>("dlenSigMin") )
   
{
   produces<nanoaod::FlatTable>("SB");
   produces<edm::PtrVector<reco::Candidate> >();
}


SoftbTableProducer::~SoftbTableProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------


void
SoftbTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    edm::Handle<std::vector<reco::Vertex>> pvsIn;
    iEvent.getByToken(pvs_, pvsIn);
    edm::Handle<edm::View<pat::Jet>> srcJet;
    iEvent.getByToken(srcJet_, srcJet);

    edm::Handle<edm::View<reco::VertexCompositePtrCandidate> > svsIn;
    iEvent.getByToken(svs_, svsIn);
    auto selCandSv = std::make_unique<PtrVector<reco::Candidate>>();
    std::vector<float> dlen,dlenSig,pAngle, dxy,dxySig ;
    std::vector<int> nTrack, jetIdx;
    VertexDistance3D vdist;
    VertexDistanceXY vdistXY;

    size_t i=0;
    unsigned int nJet = srcJet->size();
    const auto & PV0 = pvsIn->front();
    for (const auto & sv : *svsIn) {
       if (svCut_(sv)) {
           Measurement1D dl= vdist.distance(PV0,VertexState(RecoVertex::convertPos(sv.position()),RecoVertex::convertError(sv.error())));
           if(dl.value() > dlenMin_ and dl.significance() > dlenSigMin_){
               edm::Ptr<reco::Candidate> c =  svsIn->ptrAt(i); 
               selCandSv->push_back(c);

               dlen.push_back(dl.value());
               dlenSig.push_back(dl.significance());

               double dx = (PV0.x() - sv.vx()), dy = (PV0.y() - sv.vy()), dz = (PV0.z() - sv.vz());
               double pdotv = (dx * sv.px() + dy*sv.py() + dz*sv.pz())/sv.p();
               pAngle.push_back(std::acos(pdotv));
               Measurement1D d2d = vdistXY.distance(PV0, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
               dxy.push_back(d2d.value());
               dxySig.push_back(d2d.significance());
               nTrack.push_back(c->numberOfDaughters());


               int matchidx = -1;
               for (unsigned int ij = 0; ij<nJet; ij++){
                 auto jet = srcJet->ptrAt(ij);
                 if(matchByCommonSourceCandidatePtr(*c,*jet))
                 {
                   matchidx = int(ij);
                   break;
                 }
               }
               jetIdx.push_back(matchidx);
           }
       }
       i++;
    }


    auto svsTable = std::make_unique<nanoaod::FlatTable>(selCandSv->size(),svName_,false);
    // For SV we fill from here only stuff that cannot be created with the SimpleFlatTableProducer 
    svsTable->addColumn<float>("dlen",dlen,"decay length in cm",nanoaod::FlatTable::FloatColumn,10);
    svsTable->addColumn<float>("dlenSig",dlenSig,"decay length significance",nanoaod::FlatTable::FloatColumn, 10);
    svsTable->addColumn<float>("pAngle",pAngle,"pointing angle, i.e. acos(p_SV * (SV - PV)) ",nanoaod::FlatTable::FloatColumn,10);
    svsTable->addColumn<float>("dxy",dxy,"decay length in cm",nanoaod::FlatTable::FloatColumn,10);
    svsTable->addColumn<float>("dxySig",dxySig,"decay length significance",nanoaod::FlatTable::FloatColumn, 10);
    svsTable->addColumn<int>("nTracks",nTrack,"number of trakcs",nanoaod::FlatTable::IntColumn);
    svsTable->addColumn<int>("JetIdx",jetIdx,"index of the associated jet (-1 if none)",nanoaod::FlatTable::IntColumn);

    iEvent.put(std::move(svsTable),"SB");
    iEvent.put(std::move(selCandSv));
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
SoftbTableProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
SoftbTableProducer::endStream() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SoftbTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SoftbTableProducer);
