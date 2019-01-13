// -*- C++ -*-
//
// Package:    TopTagger/TopTagger
// Class:      prodIsoTracksProducer
// 
/**\class prodIsoTracksProducer prodIsoTracksProducer.cc 

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Matthew Kilpatrick
//         Created:  Sun, 13 Jan 2019 21:32:56 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <string>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include <DataFormats/Math/interface/deltaR.h>
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"

// additional headers
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TString.h"
#include "TTree.h"
#include <TFile.h>

#include "TLorentzVector.h"

class prodIsoTracksProducer : public edm::stream::EDProducer<> 
{
public:
    explicit prodIsoTracksProducer(const edm::ParameterSet&);
    ~prodIsoTracksProducer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
     
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT< std::vector<reco::Vertex> >VtxTok_;
    edm::EDGetTokenT<edm::View<reco::MET> >MetTok_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> Loose_IsoTrksHandle_Tok_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate> > Loose_IsoTrksHandlePtr_Tok_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> ForVetoIsoTrks_Tok_;
    edm::EDGetTokenT<std::vector<double> >  Loose_Isotrk_IsoVecHandle_Tok_;
    edm::EDGetTokenT<std::vector<double> >Loose_Isotrk_DzpvVecHandle_Tok_;
    std::vector<int> exclPdgIdVec_;
    double dR_, dzcut_;
    double minPt_, isoCut_;
    int loose_nIsoTrks, nIsoTrksForVeto;
    bool debug_;
    const std::string isoTrackName_;
};

prodIsoTracksProducer::prodIsoTracksProducer(const edm::ParameterSet& params):
    VtxTok_                             (consumes<std::vector<reco::Vertex>>     (params.getParameter<edm::InputTag>("vtxSrc"))),
    MetTok_                             (consumes<edm::View<reco::MET>>          (params.getParameter<edm::InputTag>("metSrc"))),
    Loose_IsoTrksHandle_Tok_            (consumes<pat::PackedCandidateCollection>(params.getParameter<edm::InputTag>("loose_isoTrkSrc"))),
    Loose_IsoTrksHandlePtr_Tok_		(consumes<edm::View<pat::PackedCandidate>>( params.getParameter<edm::InputTag>("loose_isoTrkSrc") )),
    ForVetoIsoTrks_Tok_                 (consumes<pat::PackedCandidateCollection>(params.getParameter<edm::InputTag>("forVetoIsoTrkSrc"))),
    exclPdgIdVec_                       (params.getParameter<std::vector<int>>   ("exclPdgIdVec")),
    dR_                                 (params.getParameter<double>             ("dR_ConeSize")),
    dzcut_                              (params.getParameter<double>             ("dz_CutValue")),
    minPt_                              (params.getParameter<double>             ("minPt_PFCandidate")),
    isoCut_                             (params.getParameter<double>             ("isoCut")),
    debug_                              (params.getParameter<bool>               ("debug")),
    isoTrackName_			(params.getParameter<std::string>        ("isoTrackName"))

{
    produces<nanoaod::FlatTable>(isoTrackName_);
    produces<edm::PtrVector<reco::Candidate> >();
}

prodIsoTracksProducer::~prodIsoTracksProducer()
{
 
}

//
// member functions
//


// ------------ method called to produce the data  ------------
void prodIsoTracksProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    edm::Handle< std::vector<reco::Vertex> > vertices;
    iEvent.getByToken(VtxTok_, vertices);
    edm::Handle<edm::View<reco::MET> > met;
    iEvent.getByToken(MetTok_, met);
    edm::Handle<edm::View<pat::PackedCandidate> > loose_isoTrksHandlePtr_;
    iEvent.getByToken(Loose_IsoTrksHandlePtr_Tok_, loose_isoTrksHandlePtr_);
    edm::Handle<pat::PackedCandidateCollection> loose_isoTrksHandle_;
    iEvent.getByToken(Loose_IsoTrksHandle_Tok_, loose_isoTrksHandle_);
    edm::Handle<pat::PackedCandidateCollection> forVetoIsoTrks_;
    iEvent.getByToken(ForVetoIsoTrks_Tok_, forVetoIsoTrks_);

    std::vector<float>  loose_isoTrks_pt;
    std::vector<float>  loose_isoTrks_eta;
    std::vector<float>  loose_isoTrks_phi;
    std::vector<float>  loose_isoTrks_mass;
    std::vector<double>  loose_isoTrks_charge;
    std::vector<double>  loose_isoTrks_dz;
    std::vector<int>  loose_isoTrks_pdgId;
    std::vector<int>  loose_isoTrks_idx;
    std::vector<double>  loose_isoTrks_iso;
    std::vector<double>  loose_isoTrks_mtw;
    std::vector<double>  loose_isoTrks_pfActivity;
    std::vector<int>  forVetoIsoTrks_idx;

    auto selCandPf = std::make_unique<PtrVector<reco::Candidate>>();

    if( loose_isoTrksHandle_.isValid() ) loose_nIsoTrks = loose_isoTrksHandle_->size(); else loose_nIsoTrks =0;
    if( forVetoIsoTrks_.isValid() ) nIsoTrksForVeto = forVetoIsoTrks_->size(); else nIsoTrksForVeto =0;
    
    if( vertices->size() > 0) {
      if( debug_ ) std::cout<<"\nloose_nIsoTrks : "<<loose_nIsoTrks<<"  nIsoTrksForVeto : "<<nIsoTrksForVeto<<std::endl;
    
      for(int is=0; is<loose_nIsoTrks; is++){
         const pat::PackedCandidate isoTrk = (*loose_isoTrksHandle_)[is];
    
         if( isoTrk.charge() == 0 ) continue;
         if( std::isnan(isoTrk.pt()) || std::isinf(isoTrk.pt()) ) continue;
         if( isoTrk.pt() < minPt_ ) continue;
    
         if(isoTrk.charge() != 0){
           double trkiso = 0.0;
           for(int iP=0; iP<loose_nIsoTrks; iP++){
              const pat::PackedCandidate isoTrk_other = (*loose_isoTrksHandle_)[iP];
    
              // don't count the PFCandidate in its own isolation sum
              if( is == iP ) continue;
              if( isoTrk_other.charge() == 0 ) continue;
    
              // cut on dR between the PFCandidates
              double dR = deltaR(isoTrk.p4(), isoTrk_other.p4());
              if( dR > dR_ ) continue;
    
              // cut on the PFCandidate dz
              double dz_other = isoTrk_other.dz();
              if( fabs(dz_other) > dzcut_ ) continue;
              if(fabs(isoTrk_other.dxy()) > 0.2) continue;
              if( std::find( exclPdgIdVec_.begin(), exclPdgIdVec_.end(), isoTrk_other.pdgId() ) != exclPdgIdVec_.end() ) continue;
    
              trkiso += isoTrk_other.pt();
           }
    
           if(fabs(isoTrk.dxy()) > 0.2) continue;
           if( trkiso/isoTrk.pt() > isoCut_ ) continue;
           if( std::abs(isoTrk.dz()) > dzcut_ ) continue;
           
	   edm::Ptr<reco::Candidate> c = loose_isoTrksHandlePtr_->ptrAt(is);
	   selCandPf->push_back(c);
           double mtw = sqrt( 2*( (*met)[0].pt()*isoTrk.pt() -( (*met)[0].px()*isoTrk.px() + (*met)[0].py()*isoTrk.py() ) ) );
           
           loose_isoTrks_pt.push_back(isoTrk.pt());
           loose_isoTrks_eta.push_back(isoTrk.eta());
           loose_isoTrks_phi.push_back(isoTrk.phi());
           loose_isoTrks_mass.push_back(isoTrk.energy());
           loose_isoTrks_charge.push_back(isoTrk.charge());
           loose_isoTrks_dz.push_back(isoTrk.dz());
           loose_isoTrks_pdgId.push_back(isoTrk.pdgId());
           loose_isoTrks_iso.push_back(trkiso);
           loose_isoTrks_mtw.push_back(mtw);
           
           if( debug_ ){
              std::cout<<"  --> is : "<<is<<"  pt/eta/phi/chg : "<<isoTrk.pt()<<"/"<<isoTrk.eta()<<"/"<<isoTrk.phi()<<"/"<<isoTrk.charge()<<"  mtw : "<<mtw<<"  pdgId : "<<(*loose_isoTrksHandle_)[is].pdgId()<<"  dz : " << isoTrk.dz()<<"  iso/pt : "<<trkiso/isoTrk.pt()<<std::endl;
           }  
         }else{
               //neutral particle, set trkiso and dzpv to 9999
         }     
      }  
      if( debug_ ) std::cout<<std::endl;
      
      auto out = std::make_unique<nanoaod::FlatTable>(selCandPf->size(), isoTrackName_, false); 
      out->setDoc("save Iso Track variables");
      

      //out->addColumn<double>("q", 		loose_isoTrks_charge, "charge of track" , nanoaod::FlatTable::FloatColumn,10);
      //out->addColumn<double>("dz", 		loose_isoTrks_dz,     "dz of track" , nanoaod::FlatTable::FloatColumn,10);
      //out->addColumn<int>(   "pdgid", 	        loose_isoTrks_pdgId,  "pdgid of track" , nanoaod::FlatTable::IntColumn);
      out->addColumn<double>("iso", 		loose_isoTrks_iso,    "iso of track" , nanoaod::FlatTable::FloatColumn,10);
      out->addColumn<double>("mtw", 		loose_isoTrks_mtw,    "mtw of track" , nanoaod::FlatTable::FloatColumn,10);
  
      auto single = std::make_unique<nanoaod::FlatTable>(1, isoTrackName_ + "_single", false);
      single->setDoc("save single values for iso track variables");

      single->addColumn<int>("nIsoTrks", loose_nIsoTrks,      "total number of iso tracks", nanoaod::FlatTable::IntColumn, 1);
      single->addColumn<int>("nIsoTrackVeto", nIsoTrksForVeto,"total number of iso tracks for veto", nanoaod::FlatTable::IntColumn, 1);

      iEvent.put(std::move(out),isoTrackName_);
      iEvent.put(std::move(single),isoTrackName_ + "_single");
      iEvent.put(std::move(selCandPf));

    }
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void prodIsoTracksProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void prodIsoTracksProducer::endStream() 
{
}
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void prodIsoTracksProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(prodIsoTracksProducer);
