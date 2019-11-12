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

private: 
    /**
     *  @brief  3D split point class
     */
    class SplitPoint3D
    {
    public:
        /**
         *  @brief  Constructor from three 2D hits (must be from different views)
         *
         *  @param  pandora the pandora instance, required to access the geomerty plugin
         *  @param  pHitA the first hit
         *  @param  pHitB the second hit
         *  @param  pHitC the third hit
         */
        SplitPoint3D(const pandora::Pandora &pandora, const pandora::CaloHit *const pHitA, const pandora::CaloHit *const pHitB, const pandora::CaloHit *const pHitC);
        
        /**
         *  @brief  Constructor from two 2D hits (must be from different views)
         *
         *  @param  pandora the pandora instance, required to access the geomerty plugin
         *  @param  pHitA the first hit
         *  @param  pHitB the second hit
         */
        SplitPoint3D(const pandora::Pandora &pandora, const pandora::CaloHit *const pHitA, const pandora::CaloHit *const pHitB);

        /**
         *  @brief  Get the address of the hit with the given view
         *
         *  @param  view the view to search for
         *
         *  @return the hit (A, B or C) with the given view
         */
        const pandora::CaloHit *GetHitWithView(const pandora::HitType &view) const;

        const pandora::CaloHit   *m_pHitA;       ///< the first hit
        const pandora::CaloHit   *m_pHitB;       ///< the second hit
        const pandora::CaloHit   *m_pHitC;       ///< the third hit 
        bool                      m_hasHitC;     ///< if the third hit is used
        pandora::CartesianVector  m_position3D;  ///< the 3d position corresponding to the hits
        float                     m_chi2;        ///< the measure of the extent to which the hit positions are compatible
        float                     m_maxDeltaX;   ///< the largest separation between any of the hits in X

    private:
        /**
         *  @brief  Check that the hit types are unique and either U, V or W
         */
        void CheckValidHitTypes() const;

        /**
         *  @brief  Set the maximum delta x value between the hits
         */
        void SetMaxDeltaX();
    };

    typedef std::map<const pandora::CaloHit *const, std::map<const pandora::CaloHit *const, float> > HitSeparationMap;
    typedef std::map<const pandora::CaloHit *, std::vector<const pandora::CaloHit *> > CaloHitHierarchyMap;

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
     *  @brief  Get the mapping between pairs of hits and their separation, only stored if separation is within the isolated hit distance threshold
     *
     *  @param  caloHitList the input list of calo hits
     *  @param  separationMap the output separation map
     */
    void GetHitSeparationMap(const pandora::CaloHitList &caloHitList, HitSeparationMap &separationMap) const;

    /**
     *  @brief  Get the mapping between pairs of hits and their separation, only stored if separation is within the user specified isolated hit distance threshold
     *
     *  @param  caloHitList the input list of calo hits
     *  @param  separationMap the output separation map
     *  @param  isolatedHitDistance the isolated hit distanc threshold
     */
    void GetHitSeparationMap(const pandora::CaloHitList &caloHitList, HitSeparationMap &separationMap, const float isolatedHitDistance) const;

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
     *  @brief  Get the hit in the input list that is closest to the input point
     *
     *  @param  caloHitList the input list of hits
     *  @param  point the input point
     *
     *  @return the closest hit to the input point
     */
    const pandora::CaloHit *GetClosestHitToPoint(const pandora::CaloHitList &caloHitList, const pandora::CartesianVector &point) const;

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

    /**
     *  @brief  Build the hierarchy of hits
     *
     *  @param  pVertex the neutrino vertex
     *  @param  continuousSegments the input vector of continuous segments
     *  @param  hitHierarchy the ouput hit hierarchy map
     */
    void BuildHitHierarchy(const pandora::Vertex *const pVertex, const std::vector<pandora::CaloHitList> &continuousSegments, CaloHitHierarchyMap &hitHierarchy) const;


    // TODO doxygen comments
    void BuildInitialHitHierarchy(const pandora::CaloHitList &segment, const pandora::CaloHit *const pSeedHit, CaloHitHierarchyMap &hitHierarchy) const;
    void MakeBestHierarchyLink(const HitSeparationMap &separationMap, std::vector<pandora::CaloHitList> &orphanSegments, CaloHitHierarchyMap &hitHierarchy) const;
    void AddSegmentToHitHierarchy(const unsigned int nearestSegmentIndex, const bool matchFront, const pandora::CaloHit *const pNearestHit, std::vector<pandora::CaloHitList> &orphanSegments, CaloHitHierarchyMap &hitHierarchy) const;



    void PrintHitHierarchyMap(const CaloHitHierarchyMap &hitHierarchy, const pandora::CaloHit *const pSeedHit, const unsigned int depth) const;
    void GetSplitHits(const std::vector<pandora::CaloHitList> &continuousSegments, const HitSeparationMap &separationMap, pandora::CaloHitList &splitHits) const;
    void GetKinkHits(const pandora::CaloHitList &segment, pandora::CaloHitList &kinkHits) const;
    float GetKinkAngle(const pandora::CartesianVector &thisHitPos, const pandora::CaloHitList &preHits, const pandora::CaloHitList &postHits) const;
    unsigned int GetIndexWithMaxKinkAngle(const std::vector<std::pair<unsigned int, float> > &hitIndexCosThetaBunch) const;
    void GetBifurcationHits(const std::vector<pandora::CaloHitList> &continuousSegments, const HitSeparationMap &separationMap, pandora::CaloHitList &bifurcationHits) const;

    void Get3DSplitPoints(const pandora::CaloHitList &allHits, const pandora::CaloHitList &splitHits, std::vector<SplitPoint3D> &splitPoints3D) const;
    void OrganiseHitsByView(const pandora::CaloHitList &allHits, pandora::CaloHitList &hitsU, pandora::CaloHitList &hitsV, pandora::CaloHitList &hitsW) const;

    bool MakeThreeViewMatch(pandora::CaloHitList &hitsU, pandora::CaloHitList &hitsV, pandora::CaloHitList &hitsW, std::vector<SplitPoint3D> &splitPoints3D) const;
    bool MakeTwoViewMatch(const pandora::CaloHitList &allHitsU, const pandora::CaloHitList &allHitsV, const pandora::CaloHitList &allHitsW, pandora::CaloHitList &hitsU, pandora::CaloHitList &hitsV, pandora::CaloHitList &hitsW, std::vector<SplitPoint3D> &splitPoints3D) const;
    void GetViableTwoViewMatches(const pandora::CaloHitList &hitsA, const pandora::CaloHitList &hitsB, const pandora::CaloHitList &hitsC, std::vector<SplitPoint3D> &splitPointsAB) const;

    pandora::StringVector m_pfoListNames;               ///< The input vector of pfo list names
    //pandora::StringVector m_clusterListNames;
    std::string           m_vertexListName;             ///< The input list of vectors - must be of size one: the neutrino vertex
    float                 m_isolatedHitDistance;        ///< The separation distance from all other hits for to be considered isolated
    float                 m_stitchingThreshold;         ///< The threshold distance within which two segments of hits should be stitched
    float                 m_transverseBias;             ///< The amount by which we bias the transverse coordinate over the longitudinal when finding the nearest neighbor
    int                   m_nSampleHits;                ///< The number of hits to sample either side of a given hit to find a kink
    float                 m_cosAngleThreshold;          ///< The cosine angle below which we identify a kink
    float                 m_maxMatchDeltaX;             ///< The maxiumum delta X between kink positions to be considered a match between views
    float                 m_maxMatch3ViewChi2;          ///< The maximum chi2 for kinks from three views to be considered a match
    float                 m_twoViewProjectionThreshold; ///< The threshold distance for a 2-view match when we project the 3D position into the remaining view
};

} // namespace lar_content

#endif // #ifndef LAR_SECONDARY_INTERACTIONS_ALGORITHM_H
