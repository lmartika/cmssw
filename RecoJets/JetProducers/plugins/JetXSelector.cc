/* \class PFJetSelector
 *
 * Selects jets with a configurable string-based cut,
 * and also writes out the constituents of the jet
 * into a separate collection.
 *
 * \author: Sal Rappoccio
 *
 *
 *  Addendum:  Selects jets that contain an X particle, which I'm using for
 * compositeCandidates (Matt)
 *
 * for more details about the cut syntax, see the documentation
 * page below:
 *
 *   https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
 *
 *
 */

#include "FWCore/Framework/interface/one/EDFilter.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

template <class T, typename C = std::vector<typename T::ConstituentTypeFwdPtr>>
class JetXSelector : public edm::one::EDFilter<> {

public:
  typedef std::vector<T> JetsOutput;
  typedef C ConstituentsOutput;

  JetXSelector(edm::ParameterSet const &params)
      : srcToken_(consumes<typename edm::View<T>>(
            params.getParameter<edm::InputTag>("src"))),
        cut_(params.getParameter<std::string>("cut")), filter_(false),
        selector_(cut_) {
    produces<JetsOutput>();
    produces<ConstituentsOutput>("constituents");
  }

  ~JetXSelector() override{};

  virtual void beginJob(){};
  virtual void endJob(){};
  /*
  // Default initialization is for edm::FwdPtr. Specialization (below) is for
  edm::Ptr. typename ConstituentsOutput::value_type const
  initptr(edm::Ptr<pat::PackedCandidate> const& dau) const { return typename
  ConstituentsOutput::value_type(dau, dau);
  }
  */

  virtual bool filter(edm::Event &iEvent,
                      const edm::EventSetup &iSetup) override {

    auto jets = std::make_unique<std::vector<T>>();
    auto candsOut = std::make_unique<ConstituentsOutput>();

    edm::Handle<typename edm::View<T>> h_jets;
    iEvent.getByToken(srcToken_, h_jets);

    // Now set the Ptrs with the orphan handles.
    for (typename edm::View<T>::const_iterator ibegin = h_jets->begin(),
                                               iend = h_jets->end(),
                                               ijet = ibegin;
         ijet != iend; ++ijet) {

      // Check the selection
      bool isJetPlusX = false;
      //std::cout<<" # of daughters "<<ijet->numberOfDaughters()<<std::endl;
      if (selector_(*ijet)) {
        //for (unsigned int ida = 0; ida < ijet->numberOfDaughters(); ++ida) {
	for (const auto& _d : ijet->daughterPtrVector()) {
	  const auto daughter = dynamic_cast<const reco::Candidate*>(_d.get());


          // const reco::PFCandidate* daughter = dynamic_cast<const
          // reco::PFCandidate*>(ijet->daughterPtr(ida).get());
      /*  // was using this one
      const pat::PackedCandidate *daughter =
              dynamic_cast<const pat::PackedCandidate *>(
                  ijet->daughterPtr(ida).get());
      */
	  //std::cout<<" getting mass "<<std::endl;
          // reco::PFCandidate::ParticleType pfCandidateType =
          // daughter->particleId(); int pfCandidateType = daughter->pdgId();
          //float pfCandidateMass = daughter->mass();
	  float pfCandidateMass = daughter->mass();
          if (pfCandidateMass < -0.01){
	    //std::cout<<" hooray!  "<<daughter->mass()<<std::endl;
            isJetPlusX = true;
            break;
          }
        }

        if (!isJetPlusX)
          continue; // Add the jets that pass to the output collection
        jets->push_back(*ijet);
        for (unsigned int ida = 0; ida < ijet->numberOfDaughters(); ++ida) {
          candsOut->push_back(typename ConstituentsOutput::value_type(
              ijet->daughterPtr(ida), ijet->daughterPtr(ida)));
        }
      }
    }
    // put  in Event
    bool pass = jets->size() > 0;
    iEvent.put(std::move(jets));
    iEvent.put(std::move(candsOut), "constituents");

    if (filter_)
      return pass;
    else
      return true;
  }

protected:
  edm::EDGetTokenT<typename edm::View<T>> srcToken_;
  std::string cut_;
  bool filter_;
  StringCutObjectSelector<T> selector_;
};

// typedef JetXSelector<reco::PFJet, std::vector< edm::FwdPtr<reco::PFCandidate>
// > > PFJetXSelector; typedef JetXSelector<reco::PFJet, std::vector<
// edm::FwdPtr<pat::PackedCandidate> >> PFJetXSelector;
typedef JetXSelector<reco::PFJet,
                     std::vector<edm::FwdPtr<pat::PackedCandidate>>>
    PFJetXSelector;
// typedef JetXSelector<reco::PFJet, std::vector< edm::Ptr<pat::PackedCandidate>
// > >  PFJetXSelectorPtr;
typedef JetXSelector<pat::Jet, std::vector<edm::FwdPtr<pat::PackedCandidate>>>
    PatJetXSelector;

DEFINE_FWK_MODULE(PFJetXSelector);
// DEFINE_FWK_MODULE( PFJetXSelectorPtr );
DEFINE_FWK_MODULE(PatJetXSelector);
