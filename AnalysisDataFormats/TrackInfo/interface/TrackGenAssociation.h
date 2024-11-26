#ifndef TrackGenAssociation_h
#define TrackGenAssociation_h

#include <vector>
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefProd.h"

namespace reco {
    // association map
    typedef edm::AssociationMap<edm::OneToOne<std::vector<pat::PackedCandidate>, std::vector<pat::PackedGenParticle>>> TrackGenAssociationCollection;

    typedef TrackGenAssociationCollection::value_type TrackGenAssociation;

    // reference to an object in a collection of SeedMap objects
    typedef edm::Ref<TrackGenAssociationCollection> TrackGenAssociationRef;
   
    // reference to a collection of SeedMap objects
    typedef edm::RefProd<TrackGenAssociationCollection> TrackGenAssociationRefProd;
       
    // vector of references to objects in the same colletion of SeedMap objects
    typedef edm::RefVector<TrackGenAssociationCollection> TrackGenAssociationRefVector;
}

#endif