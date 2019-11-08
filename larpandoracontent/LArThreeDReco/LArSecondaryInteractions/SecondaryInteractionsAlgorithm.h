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
     *  @brief  Reorder the input hit by nearest neighbour starting with the vertex position, remove any isolated hits then return a vector of continuous segments
     *
     *  @param  pVertex the input neutrino interaction vertex
     *  @param  caloHitList the input hit list to order
     *  @param  outputSegments the filtered and ordered output hits arranged into continous segments
     */
    void GetContinuousSegments(const pandora::Vertex *const pVertex, const pandora::CaloHitList &caloHitList, std::vector<pandora::CaloHitList> &outputSegments) const;

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

    /**
     *  @brief  Stitch together continuous segments if their end-points are within some threshold
     *
     *  @param  initialSegments the input vector of continuous segments
     *  @param  stitchedSegments the output vector of stitched segments
     */
    void StitchSegments(const std::vector<pandora::CaloHitList> &initialSegments, std::vector<pandora::CaloHitList> &stitchedSegments) const;

    // TODO doxygen comments
    void GetKinkHits(const pandora::CaloHitList &segment, pandora::CaloHitList &kinkHits) const;
    void GetKinkHits(const std::vector<pandora::CaloHitList> &continuousSegments, pandora::CaloHitList &kinkHits) const;
    float GetKinkAngle(const pandora::CartesianVector &thisHitPos, const pandora::CaloHitList &preHits, const pandora::CaloHitList &postHits) const;
    unsigned int GetIndexWithMaxKinkAngle(const std::vector<std::pair<unsigned int, float> > &hitIndexCosThetaBunch) const;

    pandora::StringVector m_pfoListNames;        ///< The input vector of pfo list names
    std::string           m_vertexListName;      ///< The input list of vectors - must be of size one: the neutrino vertex
    float                 m_isolatedHitDistance; ///< The separation distance from all other hits for to be considered isolated
    float                 m_stitchingThreshold;  ///< The threshold distance within which two segments of hits should be stitched
    float                 m_transverseBias;      ///< The amount by which we bias the transverse coordinate over the longitudinal when finding the nearest neighbor
    int                   m_nSampleHits;         ///< The number of hits to sample either side of a given hit to find a kink
    float                 m_cosAngleThreshold;   ///< The cosine angle below which we identify a kink
};

} // namespace lar_content

#endif // #ifndef LAR_SECONDARY_INTERACTIONS_ALGORITHM_H
