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

    typedef std::map<const pandora::ParticleFlowObject *const, pandora::CartesianVector> PfoToCartesianVectorMap;
    typedef std::map<const pandora::CaloHit *const, std::map<const pandora::CaloHit *const, float> > HitSeparationMap;
    typedef std::map<const pandora::CaloHit *, std::vector<const pandora::CaloHit *> > CaloHitHierarchyMap;
    typedef std::map<pandora::HitType, pandora::CaloHitList> ViewToHitsMap;
    typedef std::map<pandora::HitType, CaloHitHierarchyMap> ViewToHitHierarchyMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the vertex from the input list name
     *
     *  @return the output vertex position
     */
    pandora::CartesianVector GetVertexPosition() const;

    /**
     *  @brief  Collect the pfos from the input pfo list
     *
     *  @param  inputPfos the vector of input pfos to populate
     *  @param  pfoListName the name of the pfo list to collect
     */
    void CollectInputPfos(pandora::PfoList &inputPfos, const std::string &pfoListName) const;

    /**
     *  @brief  Get a vertex position for each pfo from which we should seed the kink finding
     *
     *  @param  allPfos the input list of all pfos
     *  @param  vertexPos the neutrino vertex position
     *  @param  pfoToSeedPointMap the output mapping from pfo to seed position
     */
    void GetSeedVertices(const pandora::PfoList &allPfos, const pandora::CartesianVector &vertexPos, PfoToCartesianVectorMap &pfoToSeedPointMap) const;

    /**
     *  @brief  Split the input PFO if a viable split point can be found
     *          When splitting a PFO create two new PFOs to replace the old one, the PFOs to save and delete are tracked in the output lists
     *
     *  @param  pPfo the input PFO
     *  @param  vertexPos the neutrino vertex
     *  @param  pfosToSave the new PFOs to save
     *  @param  pfosToDelete the old PFOs to delete
     */
    void SplitPfo(const pandora::ParticleFlowObject *const pPfo, const pandora::CartesianVector &vertexPos, pandora::PfoList &pfosToSave, pandora::PfoList &pfosToDelete) const;

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
     *  @param  vertexPos the input neutrino vertex
     *  @param  caloHitList the input list of hits
     *  @param  separationMap the separation map
     *  @param  continuousSegments the output vector of continuous segments
     */
    void GetContinuousSegments(const pandora::CartesianVector &vertexPos, const pandora::CaloHitList &caloHitList, const HitSeparationMap &separationMap, std::vector<pandora::CaloHitList> &continuousSegments) const;

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
     *  @param  vertexPos the input vertex
     *
     *  @return the closest hit in the input list to the vertex
     */
    const pandora::CaloHit *GetClosestHitToVertex(const pandora::CaloHitList &caloHitList, const pandora::CartesianVector &vertexPos) const;

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
     *  @param  vertexPos the neutrino vertex
     *  @param  continuousSegments the input vector of continuous segments
     *  @param  hitHierarchy the ouput hit hierarchy map
     */
    void BuildHitHierarchy(const pandora::CartesianVector &vertexPos, const std::vector<pandora::CaloHitList> &continuousSegments, CaloHitHierarchyMap &hitHierarchy) const;

    /**
     *  @brief  Get the 3D direction of the input PFO when fitted using PCA 
     *
     *  @param  pPfo the input PFO
     *  @param  trackDir the output 3D direction
     *
     *  @return boolean: true if a 3D fit was possible
     */
    bool GetTrackDirection(const pandora::ParticleFlowObject *const pPfo, pandora::CartesianVector &trackDir) const;

    /**
     *  @brief  Get the cosine of the angle between the input direction and the wires in the given view
     *
     *  @param  trackDir the input direction
     *  @param  view the wire plane to use
     *
     *  @return the cosine angle between the track direction and the wires
     */
    float GetCosAngleToWire(const pandora::CartesianVector &trackDir, const pandora::HitType &view) const;

    /**
     *  @brief  Filter the input list of segments and only return those in which we should attempt to find a kink
     *
     *  @param  trackDir the input track direction
     *  @param  view the view in which the segments reside
     *  @param  segments the input list of segments
     *
     *  @return the filtered segments
     */
    std::vector<pandora::CaloHitList> FilterSegments(const pandora::CartesianVector &trackDir, const pandora::HitType &view, const std::vector<pandora::CaloHitList> &segments) const;

    /**
     *  @brief  Build the initial hit hierarchy using the hits in the input segment. Use the seed hit as the start of the hierarchy
     *
     *  @param  segment the input segment
     *  @param  pSeedHit the seed of the hierarchy
     *  @param  hitHierarchy the output initial hit hierarchy
     */
    void BuildInitialHitHierarchy(const pandora::CaloHitList &segment, const pandora::CaloHit *const pSeedHit, CaloHitHierarchyMap &hitHierarchy) const;

    /**
     *  @brief  Find the orphan segment that's closest to an existing hit in the hierarchy and add it to the hierarchy
     *
     *  @param  separationMap the separation map between hits
     *  @param  orphanSegments the segments that aren't yet part of the heirarchy
     *  @param  hitHierarchy the hit herarchy to update
     */
    void MakeBestHierarchyLink(const HitSeparationMap &separationMap, std::vector<pandora::CaloHitList> &orphanSegments, CaloHitHierarchyMap &hitHierarchy) const;

    /**
     *  @brief  Add the input segment to the hierarchy
     *
     *  @param  nearestSegmentIndex the index of the orphan segment to add
     *  @param  matchFront if the match is at the front of the orphan segment (to determine if we need to flip the segment before adding)
     *  @param  pNearestHit the nearest hit to the segment we are adding that's already in the hierarchy
     *  @param  orphanSegments the full list of orphan segments
     *  @param  hitHierarchy the hit hierarchy to update
     */
    void AddSegmentToHitHierarchy(const unsigned int nearestSegmentIndex, const bool matchFront, const pandora::CaloHit *const pNearestHit, std::vector<pandora::CaloHitList> &orphanSegments, CaloHitHierarchyMap &hitHierarchy) const;

    /**
     *  @brief  Print the hit hierarchy map
     *
     *  @param  hitHierarchy the hit hierarchy map to print
     *  @param  pSeedHit the hit to seed the printing
     *  @param  depth the number of indents
     */
    void PrintHitHierarchyMap(const CaloHitHierarchyMap &hitHierarchy, const pandora::CaloHit *const pSeedHit, const unsigned int depth) const;

    /**
     *  @brief  Get the hits downstream of the seed hit in the hierarchy
     *
     *  @param  hitHierarchy the hit heirarchy
     *  @param  pSeedHit the seed hit downstream of which we want to collect
     *  @param  downstreamHits the output downstream hits
     */
    void GetDownstreamHits(const CaloHitHierarchyMap &hitHierarchy, const pandora::CaloHit *const pSeedHit, pandora::CaloHitList &downstreamHits) const;

    /**
     *  @brief  Get the hits in the input segments that represent a point of interest at which we might want to split
     *
     *  @param  continuousSegments the input list of continuous segments
     *  @param  separationMap the input separation map between hits
     *  @param  splitHits the output split hits
     */
    void GetSplitHits(const std::vector<pandora::CaloHitList> &continuousSegments, const HitSeparationMap &separationMap, pandora::CaloHitList &splitHits) const;

    /**
     *  @brief  Get the hits in the segment that represent a kink
     *
     *  @param  segment the input segment
     *  @param  kinkHits the output hits at a kink
     */
    void GetKinkHits(const pandora::CaloHitList &segment, pandora::CaloHitList &kinkHits) const;

    /**
     *  @brief  Get the cosine angle between the hits pre and post of the current sample hit
     *
     *  @param  thisHitPos the current sample hit
     *  @param  preHits the hits immediately before the sample hit
     *  @param  postHits the hits immediately after the sample hit
     *
     *  @return the cosine kink angle
     */
    float GetKinkAngle(const pandora::CartesianVector &thisHitPos, const pandora::CaloHitList &preHits, const pandora::CaloHitList &postHits) const;

    /**
     *  @brief  Get the index of the hit in the input bunch with the maximum kink angle
     *
     *  @param  hitIndexCosThetaBunch the input bunch of hits as hitIndex-kinkAngle pairs
     *
     *  @return the index of the hit in the input bunch with the max kink angle
     */
    unsigned int GetIndexWithMaxKinkAngle(const std::vector<std::pair<unsigned int, float> > &hitIndexCosThetaBunch) const;

    /**
     *  @brief  Get the hits in the input segments that represent a bifurcation
     *
     *  @param  continuousSegments the input segments
     *  @param  separationMap the input hit separation map
     *  @param  bifurcationHits the output bifurcation hits
     */
    void GetBifurcationHits(const std::vector<pandora::CaloHitList> &continuousSegments, const HitSeparationMap &separationMap, pandora::CaloHitList &bifurcationHits) const;

    /**
     *  @brief  Get the viable 3D split points from the input 2D split hits 
     *
     *  @param  viewToAllHitsMap the input mapping from view to all hits in that view for the current PFO
     *  @param  viewToSplitHitsMap the input mapping from view to all 2D split hits in that view for the current PFO
     *  @param  splitPoints3D the output viable 3D split points
     */
    void Get3DSplitPoints(const ViewToHitsMap &viewToAllHitsMap, const ViewToHitsMap &viewToSplitHitsMap, std::vector<SplitPoint3D> &splitPoints3D) const;

    /**
     *  @brief  Organise the hits in the input list by view
     *
     *  @param  allHits the input list of all hits
     *  @param  hitsU the output list of hits in the U view
     *  @param  hitsV the output list of hits in the V view
     *  @param  hitsW the output list of hits in the W view
     */
    void OrganiseHitsByView(const pandora::CaloHitList &allHits, pandora::CaloHitList &hitsU, pandora::CaloHitList &hitsV, pandora::CaloHitList &hitsW) const;

    /**
     *  @brief  If possible, form a 3D split point from consistent 2D split hits from all 3 views
     *
     *  @param  hitsU the input list of remaining split hits in the U view
     *  @param  hitsV the input list of remaining split hits in the V view
     *  @param  hitsW the input list of remaining split hits in the W view
     *  @param  splitPoints3D the output 3D split points
     *
     *  @return boolean: true if a 3-view match could be made
     */
    bool MakeThreeViewMatch(pandora::CaloHitList &hitsU, pandora::CaloHitList &hitsV, pandora::CaloHitList &hitsW, std::vector<SplitPoint3D> &splitPoints3D) const;

    /**
     *  @brief  If possible, form a 2D split point from consistent 2D split hits from 2 views
     *
     *  @param  allHitsU the input list of all hits in the U view
     *  @param  allHitsV the input list of all hits in the V view
     *  @param  allHitsW the input list of all hits in the W view
     *  @param  hitsU the input list of remaining split hits in the U view
     *  @param  hitsV the input list of remaining split hits in the V view
     *  @param  hitsW the input list of remaining split hits in the W view
     *  @param  splitPoints3D the output 3D split points
     *
     *  @return boolean: true if a 2-view match could be made
     */
    bool MakeTwoViewMatch(const pandora::CaloHitList &allHitsU, const pandora::CaloHitList &allHitsV, const pandora::CaloHitList &allHitsW, pandora::CaloHitList &hitsU, pandora::CaloHitList &hitsV, pandora::CaloHitList &hitsW, std::vector<SplitPoint3D> &splitPoints3D) const;

    /**
     *  @brief  Get the viable 2-view matches between the input 2D split hits
     *
     *  @param  hitsA the input list of remaining split hits in the first view
     *  @param  hitsB the input list of remaining split hits in the second view
     *  @param  hitsC the input list all hits in the third view
     *  @param  splitPointsAB the output 3D split points made from a match in the first and second view
     */
    void GetViableTwoViewMatches(const pandora::CaloHitList &hitsA, const pandora::CaloHitList &hitsB, const pandora::CaloHitList &hitsC, std::vector<SplitPoint3D> &splitPointsAB) const;

    /**
     *  @brief  Filter the input 3D split points and only keep the ones at which we might want to make a split
     *  
     *  @param  vertexPos the position of the reconstructed neutrino vertex
     *  @param  viewToAllHitsMap the input map from view to all hits in that view
     *  @param  viewToHitHierarchyMap the input map from view to hit hierarchy
     *  @param  splitPoints3D the input 3D split points
     *  @param  filteredSplitPoints3D the output list of filtered split points
     */
    void Filter3DSplitPoints(const pandora::CartesianVector &vertexPos, const ViewToHitsMap &viewToAllHitsMap, const ViewToHitHierarchyMap &viewToHitHierarchyMap, const std::vector<SplitPoint3D> &splitPoints3D, std::vector<SplitPoint3D> &filteredSplitPoints3D) const;

    /**
     *  @brief  Given a list of 3D split points, choose the best candidate for making a split
     *
     *  @param  splitPoints3D the input split points to choose from
     *  @param  viewToAllHitsMap the input map from view to all hits in that view
     *  @param  viewToHitHierarchyMap the input map from view to hit hierarchy
     *
     *  @return the best split point
     */
    SplitPoint3D Select3DSplitPoint(const std::vector<SplitPoint3D> &splitPoints3D, const ViewToHitsMap &viewToAllHitsMap, const ViewToHitHierarchyMap &viewToHitHierarchyMap) const;

    /**
     *  @brief  Get the score for a given split point, higher scores are better split points  
     *
     *  @param  splitPoint the input 3D split point
     *  @param  viewToAllHitsMap the input map from view to all hits in that view
     *  @param  viewToHitHierarchyMap the input map from view to hit hierarchy
     *
     *  @return the score for the input split point
     */
    float GetSplitPointScore(const SplitPoint3D &splitPoint, const ViewToHitsMap &viewToAllHitsMap, const ViewToHitHierarchyMap &viewToHitHierarchyMap) const;

    /**
     *  @brief  Split the input list of calo hits at the split point into upstream and downstream segments 
     *
     *  @param  view the view which the input calo hits are in
     *  @param  caloHits the input calo hits to split
     *  @param  hitHierarchyMap the input hit hierarchy map
     *  @param  chosenSplitPoint the 3D point at which to split
     *  @param  upstreamHitsMap the output upstream hits
     *  @param  downstreamHitsMap the output downstream hits
     */
    void SplitHitListAtPoint(const pandora::HitType view, const pandora::CaloHitList &caloHits, const CaloHitHierarchyMap &hitHierarchyMap, const SplitPoint3D &chosenSplitPoint, pandora::CaloHitList &upstreamHitsMap, pandora::CaloHitList &downstreamHitsMap) const;

    /**
     *  @brief  Split the hits in the current PFO upstream and downstream of the input 3D split point
     *
     *  @param  viewToCaloHitsMap the input map from view to all hits in that view
     *  @param  viewToHitHierarchyMap the input map from view to the hit hierarchy in that view
     *  @param  chosenSplitPoint the 3D point at which to make a split
     *  @param  viewToUpstreamHitsMap the output map from view to upstream hits
     *  @param  viewToDownstreamHitsMap the output map from view to downstream hits
     */
    void SplitPfoHits(const ViewToHitsMap &viewToCaloHitsMap, const ViewToHitHierarchyMap &viewToHitHierarchyMap, const SplitPoint3D &chosenSplitPoint, ViewToHitsMap &viewToUpstreamHitsMap, ViewToHitsMap &viewToDownstreamHitsMap) const;

    /**
     *  @brief  Break the clusters in the input PFO into new upstream and downstream clusters
     *
     *  @param  pPfo the input PFO
     *  @param  viewToUpstreamHitsMap the input map from view to upstream hits
     *  @param  viewToDownstreamHitsMap the input map from view to downstream hits
     *  @param  upstreamClusters the new output upstream clusters
     *  @param  downstreamClusters the new output downstream clusters
     */
    void SplitPfoClusters(const pandora::ParticleFlowObject *const pPfo, const ViewToHitsMap &viewToUpstreamHitsMap, const ViewToHitsMap &viewToDownstreamHitsMap, pandora::ClusterList &upstreamClusters, pandora::ClusterList &downstreamClusters) const;

    /**
     *  @brief  Split the input PFO into two new upstream and downstream PFOs  
     *
     *  @param  pPfo the input PFO
     *  @param  upstreamClusters the input upstream clusters
     *  @param  downstreamClusters the input downstream clusters
     *  @param  pfosToSave the new PFOs before and after the split that we should save
     *  @param  pfosToDelete the old PFOs before splitting that we need to delete
     */
    void SplitPfo(const pandora::ParticleFlowObject *const pPfo, const pandora::ClusterList &upstreamClusters, const pandora::ClusterList &downstreamClusters, pandora::PfoList &pfosToSave, pandora::PfoList &pfosToDelete) const;

    pandora::StringVector m_pfoListNames;               ///< The input pfo list names to modify
    std::string           m_vertexListName;             ///< The input list of vertices - must be of size one: the neutrino vertex
    std::string           m_clusterListNameU;           ///< The input list of clusters in the U view
    std::string           m_clusterListNameV;           ///< The input list of clusters in the V view
    std::string           m_clusterListNameW;           ///< The input list of clusters in the W view
    std::string           m_trackClusterListName;       ///< The input list of 3D track clusters
    std::string           m_showerClusterListName;      ///< The input list of 3D shower clusters
    float                 m_isolatedHitDistance;        ///< The separation distance from all other hits for to be considered isolated
    float                 m_stitchingThreshold;         ///< The threshold distance within which two segments of hits should be stitched
    float                 m_transverseBias;             ///< The amount by which we bias the transverse coordinate over the longitudinal when finding the nearest neighbor
    unsigned int          m_nSampleHits;                ///< The number of hits to sample either side of a given hit to find a kink
    float                 m_cos3DAngleToWireThreshold;  ///< The angular threshold to between a PFO and a wire for a segment to be used in the kink finding
    unsigned int          m_minHitsThreshold;           ///< The minimum number of hits in a segment to be used in the kink finding
    float                 m_cosAngleThreshold;          ///< The cosine angle below which we identify a kink
    float                 m_maxMatchDeltaX;             ///< The maxiumum delta X between kink positions to be considered a match between views
    float                 m_maxMatch3ViewChi2;          ///< The maximum chi2 for kinks from three views to be considered a match
    float                 m_twoViewProjectionThreshold; ///< The threshold distance for a 2-view match when we project the 3D position into the remaining view
    float                 m_cos3DAngleThreshold;        ///< The cosine angle below which we identify a kink in 3D
    float                 m_minVertexDist;              ///< The minimum distance between a 3D split position and the neutrino vertex
};

} // namespace lar_content

#endif // #ifndef LAR_SECONDARY_INTERACTIONS_ALGORITHM_H
