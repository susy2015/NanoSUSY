// -*- C++ -*-
//
// Package:    TopTagger/TopTagger
// Class:      TauMVAProducer
// 
/**\class TauMVAProducer TauMVAProducer.cc TopTagger/TopTagger/plugins/TauMVAProducer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Matthew Kilpatrick
//         Created:  Thu, 20 Dec 2018 21:32:56 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "TLorentzVector.h"

class TauMVAProducer : public edm::stream::EDProducer<> 
{
public:
    explicit TauMVAProducer(const edm::ParameterSet&);
    ~TauMVAProducer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    bool  isA(int particleID, int pdgId, bool checkCharge = false);
    int   getNearPhotonIndex(const pat::PackedCandidate* pfc, edm::Handle<std::vector<pat::PackedCandidate> >& pfcands_);
    float getDRNearestTrack(const pat::PackedCandidate* particle, edm::Handle<std::vector<pat::PackedCandidate> >& pfcands_, const float minTrackPt=1.0);
    float computePFIsolation(const pat::PackedCandidate* particle, edm::Handle<std::vector<pat::PackedCandidate> >& pfcands, const float minDR, const float maxDR, const unsigned int isotype=0, const float minNeutralPt=0.5, const float minPUPt=0.5, const float dBeta=0.5);
    float TrackIso(const pat::PackedCandidate* particle, edm::Handle<std::vector<pat::PackedCandidate> >& pfcands_, const float maxDR=0.3, const float deltaZCut=0.1);
     
private:
     
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<std::vector<pat::Jet> >             jetToken_;
    edm::EDGetTokenT<std::vector<pat::PackedCandidate> > pfCandToken_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate> >   pfCandPtrToken_;
    double        candptMin_;
    double        candetaMax_;
    const std::string tauName_;
};

TauMVAProducer::TauMVAProducer(const edm::ParameterSet& params):
    jetToken_(consumes<std::vector<pat::Jet>>( params.getParameter<edm::InputTag>("jetsSrc") )),
    pfCandToken_(consumes<std::vector<pat::PackedCandidate>>( params.getParameter<edm::InputTag>("pfCandSrc") )),
    pfCandPtrToken_(consumes<edm::View<pat::PackedCandidate>>( params.getParameter<edm::InputTag>("pfCandSrc") )),
    candptMin_( params.getParameter<double>("minCandPt") ),
    candetaMax_( params.getParameter<double>("maxCandEta") ),
    tauName_( params.getParameter<std::string>("tauName") )

{
    produces<nanoaod::FlatTable>("pfcands");
    produces<edm::PtrVector<reco::Candidate> >();
}

TauMVAProducer::~TauMVAProducer()
{
 
}

//
// member functions
//

bool TauMVAProducer::isA(int particleID, int pdgId, bool checkCharge)
{
  int id = pdgId;
  if(checkCharge) return (id == particleID);

  return (std::fabs(id) == particleID);
}

int TauMVAProducer::getNearPhotonIndex(const pat::PackedCandidate* pfc, edm::Handle<std::vector<pat::PackedCandidate> >& pfcands_)
{
  const double minPhotonPt = 0.5;
  const double maxPhotonDR = 0.2;

  int photonInd = -1;
  double maxPhotonPT = 0.0;
  int p_gamma = 22;

  for(unsigned int ic = 0; ic < pfcands_->size(); ic++) {
    const pat::PackedCandidate* c = &pfcands_->at(ic);
    if(!isA(p_gamma, c->pdgId())) continue;
    if(c->pt() < minPhotonPt) continue;
    //float dr = deltaR(c->eta(), c->phi(), pfc->eta(), pfc->phi());
    float dr = deltaR(c->p4(), pfc->p4());
    if(dr > maxPhotonDR) continue;
    if(c->pt() > maxPhotonPT) {
      maxPhotonPT = c->pt();
      photonInd = ic;
    }
  }

  return photonInd;

}

float TauMVAProducer::getDRNearestTrack(const pat::PackedCandidate* particle, edm::Handle<std::vector<pat::PackedCandidate> >& pfcands_, const float minTrackPt)
{
  float minDR = 10.0;

  for(unsigned int ic = 0; ic < pfcands_->size(); ic++) {
    const pat::PackedCandidate* cand = &pfcands_->at(ic);
    if(particle == cand) continue;

    if(cand->charge() == 0) continue;
    if(cand->pt() < minTrackPt) continue;

    const float dR = deltaR(particle->p4(), cand->p4());
    if(dR > minDR) continue;
    minDR = dR;
  }

  return minDR;

}

float TauMVAProducer::TrackIso(const pat::PackedCandidate* particle, edm::Handle<std::vector<pat::PackedCandidate> >& pfcands_, const float maxDR, const float deltaZCut)
{
  float absIso=0;

  for (unsigned int ipf = 0; ipf < pfcands_->size(); ipf++) {
    const pat::PackedCandidate* cand = &pfcands_->at(ipf);
    if(particle == cand) continue;
    if(cand->charge() == 0) continue; // skip neutrals
    const float dR = deltaR(particle->p4(), cand->p4());
    if(dR > maxDR) continue;
    if( cand->pt()>=0.0 && fabs(cand->dz()) < deltaZCut) absIso += cand->pt();
  }

  return absIso;
}

// next round: move isolation computations to IsolationVariables
float TauMVAProducer::computePFIsolation(const pat::PackedCandidate* particle, edm::Handle<std::vector<pat::PackedCandidate> >& pfcands_, const float minDR, const float maxDR, const unsigned int isotype, const float minNeutralPt, const float minPUPt, const float dBeta)
{
  float chargedIso = 0.0;
  float neutralIso = 0.0;
  float puIso = 0.0;

  for(unsigned int ic = 0; ic < pfcands_->size(); ic++) {
    const pat::PackedCandidate* cand = &pfcands_->at(ic);
    if(particle == cand) continue;

    const float dR = deltaR(particle->p4(), cand->p4());
    if(dR < minDR || dR > maxDR) continue;

    if(cand->charge() == 0) {
      if(cand->pt() > minNeutralPt) neutralIso += cand->pt();
    } else if(cand->fromPV() > 1) {
      chargedIso += cand->pt();
    } else {
      if(cand->pt() > minPUPt) puIso += cand->pt();
    }
  }

  float isolation = 0.0;

  switch(isotype) {
    case 0 : {
      isolation = chargedIso;
      break;
    }
    case 1 : {
      isolation = (chargedIso + std::max(float(0.0), (neutralIso - (dBeta*puIso))));
      break;
    }
    default : {
      printf("Isolation type not known!\n");
    }
  }

  return isolation;

}


// ------------ method called to produce the data  ------------
void TauMVAProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    edm::Handle<std::vector<pat::Jet> > jets_;
    iEvent.getByToken(jetToken_, jets_);
    edm::Handle<edm::View<pat::PackedCandidate> > pfcandsPtr_;
    iEvent.getByToken(pfCandPtrToken_, pfcandsPtr_);
    edm::Handle<std::vector<pat::PackedCandidate> > pfcands_;
    iEvent.getByToken(pfCandToken_, pfcands_);

    int p_piplus = 211;

    std::vector<float> dz, fromPV, chiso0p1, chiso0p2, chiso0p3, chiso0p4, totiso0p1, totiso0p2, totiso0p3, totiso0p4, trackiso, nearphopt, nearphoeta, nearphophi, nearestTrkDR;
    std::vector<int> contJetIndex;

    auto selCandPf = std::make_unique<PtrVector<reco::Candidate>>();

    for (unsigned int ic = 0; ic < pfcands_->size(); ic++) {
	const pat::PackedCandidate* pfc = &pfcands_->at(ic);
	if (pfc->pt() < candptMin_ || fabs(pfc->eta()) > candetaMax_) continue;
	if(!isA(p_piplus, pfc->pdgId())) continue; // only save charged hadrons unless otherwise specified

	edm::Ptr<reco::Candidate> c = pfcandsPtr_->ptrAt(ic);
	selCandPf->push_back(c);
	float chiso0p1_  = computePFIsolation(pfc, pfcands_, 0.0, 0.1, 0);
	float chiso0p2_  = computePFIsolation(pfc, pfcands_, 0.0, 0.2, 0);
	float chiso0p3_  = computePFIsolation(pfc, pfcands_, 0.0, 0.3, 0);
	float chiso0p4_  = computePFIsolation(pfc, pfcands_, 0.0, 0.4, 0);
	float totiso0p1_ = computePFIsolation(pfc, pfcands_, 0.0, 0.1, 1);
	float totiso0p2_ = computePFIsolation(pfc, pfcands_, 0.0, 0.2, 1);
	float totiso0p3_ = computePFIsolation(pfc, pfcands_, 0.0, 0.3, 1);
	float totiso0p4_ = computePFIsolation(pfc, pfcands_, 0.0, 0.4, 1);
	float trackiso_ = TrackIso(pfc, pfcands_, 0.3,0.1);
	int photonIndex_ = getNearPhotonIndex(pfc, pfcands_);
	float nearesttrkdr_ = getDRNearestTrack(pfc, pfcands_);

        int index = -1;
        for (unsigned int ijet = 0; ijet < jets_->size(); ijet++) {
          const pat::Jet &j = (*jets_)[ijet];
          for(unsigned int id = 0; id < j.numberOfDaughters(); id++) {
            const pat::PackedCandidate* dau = dynamic_cast<const pat::PackedCandidate*>(j.daughter(id));
            if(pfc == dau) {
              index = ijet;
              break;
            }
          }
        }

	dz.push_back(pfc->dz());
	fromPV.push_back(pfc->fromPV());
	chiso0p1.push_back(chiso0p1_);
	chiso0p2.push_back(chiso0p2_);
	chiso0p3.push_back(chiso0p3_);
	chiso0p4.push_back(chiso0p4_);
	totiso0p1.push_back(totiso0p1_);
	totiso0p2.push_back(totiso0p2_);
	totiso0p3.push_back(totiso0p3_);
	totiso0p4.push_back(totiso0p4_);
	trackiso.push_back(trackiso_);
	nearphopt.push_back(photonIndex_ > -1 ? pfcands_->at(photonIndex_).pt() : -10.0);
	nearphoeta.push_back(photonIndex_ > -1 ? pfcands_->at(photonIndex_).eta() : -10.0);
	nearphophi.push_back(photonIndex_ > -1 ? pfcands_->at(photonIndex_).phi() : -10.0);
	nearestTrkDR.push_back(nearesttrkdr_);
	contJetIndex.push_back(index);
       
    }

    auto out = std::make_unique<nanoaod::FlatTable>(selCandPf->size(), tauName_, false); 
    out->setDoc("save pfcand and tau mva variables");
    
    out->addColumn<float>("dz", 		dz	  , "pfcand info dz"			  , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("fromPV", 		fromPV	  ,"pfcand info from Primary Vertex"  	  , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("chiso0p1", 		chiso0p1, "charged hadron isolation with R = 0.1" , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("chiso0p2", 		chiso0p2, "charged hadron isolation with R = 0.2" , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("chiso0p3", 		chiso0p3, "charged hadron isolation with R = 0.3" , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("chiso0p4", 		chiso0p4, "charged hadron isolation with R = 0.4" , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("totiso0p1", 		totiso0p1, "total hadron isolation with R = 0.1"  , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("totiso0p2", 		totiso0p2, "total hadron isolation with R = 0.2"  , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("totiso0p3", 		totiso0p3, "total hadron isolation with R = 0.3"  , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("totiso0p4", 		totiso0p4, "total hadron isolation with R = 0.4"  , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("trackiso", 		trackiso,  "individual track isolation"		  , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("nearphopt", 		nearphopt, "pT of nearest photon to track" 	  , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("nearphoeta", 	nearphoeta, "Eta of nearest photon to track"      , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("nearphophi", 	nearphophi, "Phi of nearest photon to track"      , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<float>("nearestTrkDR", 	nearestTrkDR, "Nearest track in deltaR"           , nanoaod::FlatTable::FloatColumn,10);
    out->addColumn<int>(  "contJetIndex", 	contJetIndex,     "Index of jet containing track" , nanoaod::FlatTable::IntColumn);
  
    iEvent.put(std::move(out),"pfcands");
    iEvent.put(std::move(selCandPf));
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void TauMVAProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void TauMVAProducer::endStream() 
{
}
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TauMVAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauMVAProducer);
