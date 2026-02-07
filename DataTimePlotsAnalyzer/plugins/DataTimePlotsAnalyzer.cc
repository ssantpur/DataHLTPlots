// -*- C++ -*-
//
// Package:    DataTimePlots/DataTimePlotsAnalyzer
// Class:      DataTimePlotsAnalyzer
//
/**\class DataTimePlotsAnalyzer DataTimePlotsAnalyzer.cc DataTimePlots/DataTimePlotsAnalyzer/plugins/DataTimePlotsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthew Citron <mcitron@ucsb.edu> 10/19/2017
//         Created:  Fri, 20 Aug 2021 19:40:29 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/MuonReco/interface/MuonRecHitCluster.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "TH2.h"
#include "TTree.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class DataTimePlotsAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DataTimePlotsAnalyzer(const edm::ParameterSet&);
  ~DataTimePlotsAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;

    edm::InputTag jetLabel_;
    edm::InputTag jetTimeLabel_;
    edm::InputTag jetEmLabel_;
    edm::InputTag jetCellsLabel_;
    std::string triggerString_;
    edm::EDGetTokenT<reco::CaloJetCollection> jetInputToken;
    edm::EDGetTokenT<edm::ValueMap<float>> jetTimesInputToken;
    edm::EDGetTokenT<edm::ValueMap<unsigned int>> jetCellsInputToken;
    edm::EDGetTokenT<edm::ValueMap<float>> jetEmInputToken;
    bool triggerFired;
    TH1D * triggerFiredHist;
    TH1D * timeHist;
    TH1D * timeHistHighHad;
    TH1D * maxTimeHist;
    TH1D * maxTimeHistHighHad;
    TH1D * secondTimeHist;
    TH1D * secondTimeHistHighHad;
    TH2D * timeVsJetPtHist;
    TH2D * timeVsJetEmHist;
    TH2D * timeVsJetCellsHist;
    TH2D * timeVsJetEtaHist;
    TTree * timeTree;
    unsigned int event_i = 0;
    unsigned int run_i = 0;
    unsigned int lumi_i = 0;
    unsigned int era_i = 0;
    unsigned int bx_i = 0;
    unsigned int store_i = 0;
    unsigned int orbit_i = 0;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
    std::vector<double> * 	v_caloJetHLTPt= new std::vector<double>();
    std::vector<double> * 	v_caloJetHLTEta= new std::vector<double>();
    std::vector<double> *   v_caloJetHLTPhi= new std::vector<double>();
    std::vector<double> *   v_caloJetHLTE= new std::vector<double>();
    std::vector<double> *   v_caloJetHLTEcalTime= new std::vector<double>();
    std::vector<double> *   v_caloJetHLTEcalE= new std::vector<double>();
    std::vector<int> *   v_caloJetHLTEcalCells= new std::vector<int>();
    std::map<const TString, bool> delayedJetHLT;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DataTimePlotsAnalyzer::DataTimePlotsAnalyzer(const edm::ParameterSet& iConfig)
{ 
    edm::Service<TFileService> fs;
    jetLabel_= iConfig.getParameter<edm::InputTag>("jets");
    jetTimeLabel_= iConfig.getParameter<edm::InputTag>("jetTimes");
    jetCellsLabel_= iConfig.getParameter<edm::InputTag>("jetCells");
    jetEmLabel_= iConfig.getParameter<edm::InputTag>("jetEmEnergy");
    triggerString_ = iConfig.getParameter<std::string>("triggerString");
    jetInputToken = consumes<std::vector<reco::CaloJet>>(jetLabel_);
    jetTimesInputToken = consumes<edm::ValueMap<float>>(jetTimeLabel_);
    jetCellsInputToken = consumes<edm::ValueMap<unsigned int>>(jetCellsLabel_);
    jetEmInputToken = consumes<edm::ValueMap<float>>(jetEmLabel_);
    triggerFired = fs->make<TH1D>("triggerFired","",2,0,2);
    timeHist = fs->make<TH1D>("jetTimeHist","",200,-25,25);
    timeHistHighHad = fs->make<TH1D>("jetTimeHistHighHad","",200,-25,25);
    maxTimeHist = fs->make<TH1D>("maxJetTimeHist","",200,-25,25);
    maxTimeHistHighHad = fs->make<TH1D>("maxJetTimeHistHighHad","",200,-25,25);
    secondTimeHist = fs->make<TH1D>("secondJetTimeHist","",200,-25,25);
    secondTimeHistHighHad = fs->make<TH1D>("secondJetTimeHistHighHad","",200,-25,25);
    timeVsJetPtHist = fs->make<TH2D>("jetTimeVsPtHist","",200,-25,25,200,0,400);
    timeVsJetEtaHist = fs->make<TH2D>("jetTimeVsEtaHist","",200,-25,25,60,-1.5,1.5);
    timeVsJetCellsHist = fs->make<TH2D>("jetTimeVsnCellsHist","",200,-25,25,100,0,100);
    timeVsJetEmHist = fs->make<TH2D>("jetTimeVsEmHist","",200,-25,25,200,0,200);
    triggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","MYHLT"));
    std::vector<std::string> delayedJetStrings = {"HLT_HT425",
						  "HLT_HT430_DelayedJet40_DoubleDelay0p5nsInclusive",
						  "HLT_HT430_DelayedJet40_DoubleDelay0p5nsTrackless",
						  "HLT_HT430_DelayedJet40_DoubleDelay0p75nsTrackless",
						  "HLT_HT430_DelayedJet40_DoubleDelay1nsInclusive",
						  "HLT_HT430_DelayedJet40_DoubleDelay1nsTrackless",
						  "HLT_HT430_DelayedJet40_DoubleDelay1p25nsInclusive",
						  "HLT_HT430_DelayedJet40_DoubleDelay1p5nsInclusive",
						  "HLT_HT430_DelayedJet40_SingleDelay0p5nsInclusive",
						  "HLT_HT430_DelayedJet40_SingleDelay0p5nsTrackless",
						  "HLT_HT430_DelayedJet40_SingleDelay1nsInclusive",
						  "HLT_HT430_DelayedJet40_SingleDelay1nsTrackless",
						  "HLT_HT430_DelayedJet40_SingleDelay1p25nsTrackless",
						  "HLT_HT430_DelayedJet40_SingleDelay1p5nsInclusive",
						  "HLT_HT430_DelayedJet40_SingleDelay1p5nsTrackless",
						  "HLT_HT430_DelayedJet40_SingleDelay2nsInclusive",
						  "HLT_HT430_DelayedJet40_SingleDelay2p25nsInclusive",
						  "HLT_HT430_DelayedJet40_SingleDelay2p5nsInclusive",
						  "HLT_L1Tau_DelayedJet40_DoubleDelay0p5nsTrackless",
						  "HLT_L1Tau_DelayedJet40_DoubleDelay0p75nsInclusive",
						  "HLT_L1Tau_DelayedJet40_DoubleDelay1nsTrackless",
						  "HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsInclusive",
						  "HLT_L1Tau_DelayedJet40_DoubleDelay1p25nsTrackless",
						  "HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsInclusive",
						  "HLT_L1Tau_DelayedJet40_DoubleDelay1p5nsTrackless",
						  "HLT_L1Tau_DelayedJet40_DoubleDelay1p75nsInclusive",
						  "HLT_L1Tau_DelayedJet40_SingleDelay2p5nsTrackless",
						  "HLT_L1Tau_DelayedJet40_SingleDelay2p75nsTrackless",
						  "HLT_L1Tau_DelayedJet40_SingleDelay3nsTrackless",
						  "HLT_L1Tau_DelayedJet40_SingleDelay3p5nsInclusive",
						  "HLT_L1Tau_DelayedJet40_SingleDelay3p75nsInclusive",
						  "HLT_L1Tau_DelayedJet40_SingleDelay4nsInclusive",
						  "HLT_L1Tau_DelayedJet40_SingleDelayPos2p5nsInclusive_SatBunch",
						  "HLT_HT430_DelayedJet40_SingleDelayNeg2p5nsInclusive_SatBunch",
						  "HLT_HT430_DelayedJet40_SingleDelayPos0p75nsInclusive_BeamHalo",
						  "HLT_HT430_DelayedJet40_SingleDelayNeg0p75nsInclusive_BeamHalo",
						  "HLT_L1Tau_DelayedJet40_SingleDelayNeg2p5nsInclusive_SatBunch",
						  "HLT_L1Tau_DelayedJet40_SingleDelayNeg0p75nsInclusive_BeamHalo",
						  "HLT_L1Tau_DelayedJet40_SingleDelayPos0p75nsInclusive_BeamHalo",
						  "HLT_PFMET100_PFMHT100_IDTight_PFHT60","HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60","HLT_CaloMET90_NotCleaned"
    };
    for (auto delayedJetHLTIt = delayedJetStrings.begin(); delayedJetHLTIt != delayedJetStrings.end(); delayedJetHLTIt++){
        delayedJetHLT[*delayedJetHLTIt] = false;
    }   
    timeTree = fs->make<TTree>("timeTree","timeTree");
    timeTree->Branch("triggerFired",&triggerFired,"triggerFired/O");
    timeTree->Branch("caloJetHLT_pt",             &v_caloJetHLTPt);
    timeTree->Branch("caloJetHLT_eta",             &v_caloJetHLTEta);
    timeTree->Branch("caloJetHLT_phi",            &v_caloJetHLTPhi);
    timeTree->Branch("caloJetHLT_e",               &v_caloJetHLTE);
    timeTree->Branch("caloJetHLT_ecalE",               &v_caloJetHLTEcalE);
    timeTree->Branch("caloJetHLT_ecalCells",               &v_caloJetHLTEcalCells);
    timeTree->Branch("caloJetHLT_ecalTime",               &v_caloJetHLTEcalTime);
    timeTree->Branch("event",&event_i,"event/i");
    timeTree->Branch("lumi",&lumi_i,"lumi/i");
    timeTree->Branch("era",&era_i,"era/i");
    timeTree->Branch("bx",&bx_i,"bx/i");
    timeTree->Branch("store",&store_i,"store/i");
    timeTree->Branch("orbit",&orbit_i,"orbit/i");
    timeTree->Branch("run",&run_i,"run/i");
    for (auto iDelayedJetHLT = delayedJetHLT.begin(); iDelayedJetHLT != delayedJetHLT.end(); iDelayedJetHLT++){
        timeTree->Branch(iDelayedJetHLT->first,&iDelayedJetHLT->second,iDelayedJetHLT->first+"/O");
    }
  //now do what ever initialization is needed
}

DataTimePlotsAnalyzer::~DataTimePlotsAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void DataTimePlotsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
    v_caloJetHLTPt->clear();
    v_caloJetHLTEta->clear();
    v_caloJetHLTPhi->clear();
    v_caloJetHLTE->clear();
    v_caloJetHLTEcalTime->clear();
    v_caloJetHLTEcalCells->clear();
    v_caloJetHLTEcalE->clear();
    for (auto iDelayedJetHLT = delayedJetHLT.begin(); iDelayedJetHLT != delayedJetHLT.end(); iDelayedJetHLT++){
        iDelayedJetHLT->second = false;
    }
    lumi_i = 0;
    run_i = 0;
    event_i = 0;
    bx_i = 0;
    store_i = 0;
    orbit_i = 0;

    std::unique_ptr<unsigned int >  run   ( new unsigned int(iEvent.id().run()        ) );
    std::unique_ptr<unsigned int >  event ( new unsigned int(iEvent.id().event()      ) );
    std::unique_ptr<unsigned int >  lumi    ( new unsigned int(iEvent.luminosityBlock() ) );
    std::unique_ptr<unsigned int >  bx    ( new unsigned int(iEvent.bunchCrossing() ) );
    std::unique_ptr<unsigned int >  store    ( new unsigned int(iEvent.eventAuxiliary().storeNumber() ) );
    std::unique_ptr<unsigned int >  orbit    ( new unsigned int(iEvent.orbitNumber() ) );

    lumi_i = *lumi;
    run_i = *run;
    event_i = *event;
    bx_i = *bx;
    store_i = *store;
    orbit_i = *orbit;

    edm::Handle<reco::CaloJetCollection> jets;
    iEvent.getByToken(jetInputToken, jets);
    auto const& jetTimes = iEvent.get(jetTimesInputToken);
    auto const& jetCells = iEvent.get(jetCellsInputToken);
    auto const& jetEm = iEvent.get(jetEmInputToken);
    // edm::Handle<edm::ValueMap<int>> jetCells;
    // iEvent.getByToken(jetCellsInputToken, jetCells);
    // edm::Handle<edm::ValueMap<float>> jetEm;
    // iEvent.getByToken(jetEmInputToken, jetEm);
    Handle<edm::TriggerResults> triggerResults; //our trigger result object
    iEvent.getByToken(triggerResultsToken, triggerResults);
    unsigned int njets = 0;
    double maxTime = -25.;
    double maxTimeHad = -25.;
    double secondTime = -25.;
    double secondTimeHad = -25.;
    if (triggerResults.isValid()) {
	const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
	for (unsigned int i = 0; i < triggerResults->size(); ++i) {
	    const auto &trigname = triggerNames.triggerName(i);
	    // if (trigname.find(triggerString_ + "_v") != std::string::npos) {
		// triggerFired   = triggerResults->accept(i);
	    // }
	    for (auto iDelayedJetHLT = delayedJetHLT.begin(); iDelayedJetHLT != delayedJetHLT.end(); iDelayedJetHLT++){
                if (trigname.find(iDelayedJetHLT->first + "_v") != std::string::npos)   iDelayedJetHLT->second   = triggerResults->accept(i);
	    }
	}
    }
    for (auto iterJet = jets->begin(); iterJet != jets->end(); ++iterJet) {
	edm::Ref<vector<reco::CaloJet>> const caloJetRef(jets, std::distance(jets->begin(), iterJet));
	if (iterJet->pt() > 15){
	    v_caloJetHLTPt->push_back(iterJet->pt());
	    v_caloJetHLTPhi->push_back(iterJet->phi());
	    v_caloJetHLTEta->push_back(iterJet->eta());
	    v_caloJetHLTE->push_back(iterJet->energy());
	    v_caloJetHLTEcalTime->push_back((jetTimes)[caloJetRef]);
	    v_caloJetHLTEcalE->push_back((jetEm)[caloJetRef]);
	    v_caloJetHLTEcalCells->push_back((jetCells)[caloJetRef]);
	}
    }
    // for (uint ijet = 0; ijet < jets->size(); ijet++) {
	// auto const& obj = jets->at(ijet);
	// reco::CaloJetRef calojetref(jets, ijet);
	// v_caloJetHLTPt->push_back(obj.pt());
	// v_caloJetHLTPhi->push_back(obj.phi());
	// v_caloJetHLTEta->push_back(obj.eta());
	// v_caloJetHLTE->push_back(obj.energy());
	// // v_caloJetHLTEcalTime->push_back((*jetTimes)[calojetref]);
	// // v_caloJetHLTEcalE->push_back((*jetEm)[calojetref]);
	// // v_caloJetHLTEcalCells->push_back((*jetCells)[calojetref]);
    // }
    timeTree->Fill();
    // ijet = 0;
    // if (true){ //triggerResults.isValid()) {
	// // const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
	// // for (unsigned int i = 0; i < triggerResults->size(); ++i) {
	// //     const auto &trigname = triggerNames.triggerName(i);
	// //     if (trigname.find(triggerString_ + "_v") != std::string::npos) {
	// // 	trigger   = triggerResults->accept(i);
	// //     }
	// //     if (triggerResults->accept(i)){
	// //     // std::cout << trigname << " " << triggerResults->accept(i) << std::endl;
	// //     }
	// // }
	// if (true){
	// // if (trigger){
	//     triggerFiredHist->Fill(1.5);
	//     for (auto const& c : *jets) {
	// 	reco::CaloJetRef calojetref(jets, ijet);
	// 	timeHist->Fill((*jetTimes)[calojetref]);
	// 	if (maxTime < (*jetTimes)[calojetref]) {
	// 	    if (secondTime > -25) secondTime = maxTime;
	// 	    maxTime = (*jetTimes)[calojetref];
	// 	}
	// 	else if (secondTime < (*jetTimes)[calojetref]){
	// 	    secondTime = (*jetTimes)[calojetref];
	// 	}
	// 	if (c.energy()*c.energyFractionHadronic() > 20) {
	// 	    timeHistHighHad->Fill((*jetTimes)[calojetref]);
	// 	    if (maxTimeHad < (*jetTimes)[calojetref]) {
	// 		if (secondTimeHad > -25) secondTimeHad = maxTimeHad;
	// 		maxTimeHad = (*jetTimes)[calojetref];
	// 	    }
	// 	    else if (secondTimeHad < (*jetTimes)[calojetref]){
	// 		secondTimeHad = (*jetTimes)[calojetref];
	// 	    }
	// 		
	// 	}
	// 	timeVsJetPtHist->Fill((*jetTimes)[calojetref],c.pt());
	// 	timeVsJetEmHist->Fill((*jetTimes)[calojetref],(*jetEm)[calojetref]);
	// 	timeVsJetCellsHist->Fill((*jetTimes)[calojetref],(*jetCells)[calojetref]);
	// 	timeVsJetEtaHist->Fill((*jetTimes)[calojetref],c.eta());
	// 	ijet++;
	//     }
	//     maxTimeHist->Fill(maxTime);
	//     maxTimeHistHighHad->Fill(maxTimeHad);
	//     secondTimeHist->Fill(secondTime);
	//     secondTimeHistHighHad->Fill(secondTimeHad);
	// }
	// else{
	//     triggerFiredHist->Fill(0.5);
	// }
    // }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DataTimePlotsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    //Specify that only 'tracks' is allowed
    //To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DataTimePlotsAnalyzer);
