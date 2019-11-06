/**
 *  @file   larpandoracontent/LArThreeDReco/LArSecondaryInteractions/SecondaryInteractionsAlgorithm.h
 *
 *  @brief  Header file for the secondary interactions algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SECONDARY_INTERACTIONS_ALGORITHM_H
#define LAR_SECONDARY_INTERACTIONS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <map>

namespace lar_content
{

/**
 *  @brief  SecondaryInteractionsAlgorithm class
 */
class SecondaryInteractionsAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SecondaryInteractionsAlgorithm();
    
    pandora::StatusCode Run();

protected:
    typedef std::map<const pandora::CaloHit *const, std::map<const pandora::CaloHit *const, float> > HitSeparationMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the vertex from the input list name
     *
     *  @return pVertex the output vertex
     */
    const pandora::Vertex *GetVertex() const;

    /**
     *  @brief  Collect the pfos from the input pfo list names
     *
     *  @param  inputPfos the vector of input pfos to populate
     */
    void CollectInputPfos(pandora::PfoList &inputPfos) const;

    /**
     *  @brief  Reorder the input hit by nearest neighbour starting with the vertex position, remove any isolated hits
     *
     *  @param  pVertex the input neutrino interaction vertex
     *  @param  caloHitList the input hit list to order
     *  @param  outputHitList the filtered and ordered output hits
     */
    void FilterAndOrderHitsByNearestNeighbor(const pandora::Vertex *const pVertex, const pandora::CaloHitList &caloHitList, pandora::CaloHitList &outputHitList) const;

    /**
     *  @brief  Get the mapping between pairs of hits and their separation, only stored if separation is within isolated hit distance threshold
     *
     *  @param  caloHitList the input list of calo hits
     *  @param  separationMap the output separation map
     */
    void GetHitSeparationMap(const pandora::CaloHitList &caloHitList, HitSeparationMap &separationMap) const;

    /**
     *  @brief  Get the continuous hit segments by linking hits via nearest neighbor
     *
     *  @param  pVertex the input neutrino vertex
     *  @param  caloHitList the input list of hits
     *  @param  separationMap the separation map
     *  @param  continuousSegments the output vector of continuous segments
     */
    void GetContinuousSegments(const pandora::Vertex *const pVertex, const pandora::CaloHitList &caloHitList, const HitSeparationMap &separationMap, std::vector<pandora::CaloHitList> &continuousSegments) const;

    /**
     *  @brief  Get the hit in the input list that is closest to the vertex, if using 2D hits then use projected distances
     *
     *  @param  caloHitList the input list of hits
     *  @param  pVertex the input vertex
     *
     *  @return the closest hit in the input list to the vertex
     */
    const pandora::CaloHit *GetClosestHitToVertex(const pandora::CaloHitList &caloHitList, const pandora::Vertex *const pVertex) const;

    /**
     *  @brief  Get the hits in the input list that feature in the separation map, i.e those that aren't isolated
     *
     *  @param  caloHitList the input list of hits
     *  @param  separationMap the input separation map
     *  @param  nonIsolatedHits the output list of non-isolated hits
     */
    void GetNonIsolatedHits(const pandora::CaloHitList &caloHitList, const HitSeparationMap &separationMap, pandora::CaloHitList &nonIsolatedHits) const;

    /**
     *  @brief  Collect all of the hits that link with the input seed hit and are in the input list of remaining hits
     *
     *  @param  pSeedHit the seed hit
     *  @param  separationMap the separation map
     *  @param  remainingHits the input list of remaining hits, any hits collected in this function will be removed from this list
     *  @param  segment the output list of hits connected to the seed hit
     */
    void CollectHitsInSegment(const pandora::CaloHit *const pSeedHit, const HitSeparationMap &separationMap, pandora::CaloHitList &remainingHits, pandora::CaloHitList &segment) const;

    /**
     *  @brief  Get the nearest neighbor hit that is remaining
     *
     *  @param  separationMap the separation map
     *  @param  remainingHits the input list of remaining hits
     *  @param  pLastHit the last hit that we considered, we use the direction from last->current to define longitudinal
     *  @param  pCurrentHit the input hit from which we want to find the nearest neighbor
     *  @param  pNextHit the nearest neighbor hit of pCurrentHit
     *
     *  @return boolean: true if a nearest neighbor existed within the isolation threshold, false otherwise
     */
    bool GetNearestNeighbor(const HitSeparationMap &separationMap, const pandora::CaloHitList &remainingHits, const pandora::CaloHit *const pLastHit, const pandora::CaloHit *const pCurrentHit, const pandora::CaloHit *&pNextHit) const;

    // TODO doxygen comments
    void GetKinkIndices(const pandora::CaloHitList &segment, std::vector<unsigned int> &kinkIndices) const;
    void GetKinkIndices(const std::vector<pandora::CaloHitList> &continuousSegments, std::vector< std::vector<unsigned int> > &segmentKinkIndices) const;

    pandora::StringVector m_pfoListNames;        ///< The input vector of pfo list names
    std::string           m_vertexListName;      ///< The input list of vectors - must be of size one: the neutrino vertex
    float                 m_isolatedHitDistance; ///< The separation distance from all other hits for to be considered isolated
    float                 m_transverseBias;      ///< The amount by which we bias the transverse coordinate over the longitudinal when finding the nearest neighbor
    int                   m_nSampleHits;         ///< The number of hits to sample either side of a given hit to find a kink
};

} // namespace lar_content

#endif // #ifndef LAR_SECONDARY_INTERACTIONS_ALGORITHM_H
