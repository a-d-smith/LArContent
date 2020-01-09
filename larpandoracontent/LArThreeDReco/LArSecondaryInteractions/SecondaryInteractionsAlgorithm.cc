/**
 *  @file   larpandoracontent/LArThreeDReco/LArSecondaryInteractions/SecondaryInteractionsAlgorithm.cc
 *
 *  @brief  Implementation of the secondary interactions algorithm class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArSecondaryInteractions/SecondaryInteractionsAlgorithm.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"

using namespace pandora;

namespace lar_content
{

SecondaryInteractionsAlgorithm::SecondaryInteractionsAlgorithm() :
    m_isolatedHitDistance(2.f),
    m_stitchingThreshold(1.f),
    m_transverseBias(2.f),
    m_nSampleHits(5),
    m_cosAngleThreshold(0.9),
    m_maxMatchDeltaX(0.6f),
    m_maxMatch3ViewChi2(1.f),
    m_twoViewProjectionThreshold(1.8f)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

StatusCode SecondaryInteractionsAlgorithm::Run()
{
    // Get the input collections
    PfoList allPfos;
    for (const auto &pfoListName : m_pfoListNames)
        this->CollectInputPfos(allPfos, pfoListName);
    
    const auto vertexPos = this->GetVertexPosition();

    PfoToCartesianVectorMap pfoToSeedPointMap;
    this->GetSeedVertices(allPfos, vertexPos, pfoToSeedPointMap);

    for (const auto &pfoListName : m_pfoListNames)
    {
        // Get the input pfos in this list
        PfoList inputPfos;
        this->CollectInputPfos(inputPfos, pfoListName);

        // Make a temporary PFO list to work within
        std::string newPfoListName;
        const PfoList *pNewPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNewPfoList, newPfoListName));

        // Where desired, make new PFOs in the current list that represent existing PFOs before and after a split
        PfoList pfosToSave, pfosToDelete;
        for (const auto &pPfo : inputPfos)
        {
            const auto seedPoint = pfoToSeedPointMap.at(pPfo);
            this->SplitPfo(pPfo, seedPoint, pfosToSave, pfosToDelete);
        }

        // Save any new PFOs to the input list
        if (!pfosToSave.empty())
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, pfoListName, pfosToSave));

        // Delete any old PFOs that we have split
        if (pfosToDelete.empty())
            continue;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, pfoListName));
        for (const auto &pPfo : pfosToDelete)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo));
    }

    return STATUS_CODE_SUCCESS;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::GetSeedVertices(const PfoList &allPfos, const CartesianVector &vertexPos, PfoToCartesianVectorMap &pfoToSeedPointMap) const
{
    for (const auto &pPfo : allPfos)
        pfoToSeedPointMap.emplace(pPfo, vertexPos);

    // TODO make this less vomit in code form
    /*    
    const auto wirePitch = LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W);
    const unsigned int fitWindow = 5;
    const auto dummy = std::numeric_limits<float>::max();

    // Cache the distances between all Pfos and the vertex
    std::map<const ParticleFlowObject *const, std::map<const ParticleFlowObject *const, float> > pfoSeparationMap;
    std::map<const ParticleFlowObject *const, std::map<const ParticleFlowObject *const, const CartesianVector*> > pfoMatchPointMap;
    std::map<const ParticleFlowObject *const, float> pfoVertexDistMap;

    // Loop over the PFOs
    auto iterI = allPfos.begin();
    for (unsigned int i = 0; i < allPfos.size(); ++i)
    {
        const auto pPfoI = *iterI;

        // Initialise the map entries
        pfoSeparationMap[pPfoI];
        pfoMatchPointMap[pPfoI];
        pfoVertexDistMap[pPfoI] = std::numeric_limits<float>::max();
        
        // Get the endpoints
        CartesianVector minPosI(dummy, dummy, dummy);
        CartesianVector maxPosI(dummy, dummy, dummy);

        try
        {
            ClusterList clusters3DI;
            LArPfoHelper::GetThreeDClusterList(pPfoI, clusters3DI);

            const ThreeDSlidingFitResult slidingFitResultI(clusters3DI.front(), fitWindow, wirePitch);
            minPosI = slidingFitResultI.GetGlobalMinLayerPosition();
            maxPosI = slidingFitResultI.GetGlobalMaxLayerPosition();
        }
        catch (const StatusCodeException &)
        {
            // If we can't fit then just skip this pfo
            continue;
        }

        // Get the distance to the vertex
        pfoVertexDistMap[pPfoI] = std::min(vertexPos.GetDistanceSquared(minPosI), vertexPos.GetDistanceSquared(maxPosI));

        // Now loop over the other PFOs without double counting any pairs
        auto iterJ = allPfos.begin();
        for (unsigned int j = 0; j < i; ++j)
        {
            const auto pPfoJ = *iterJ;
        
            // Get the endpoints
            CartesianVector minPosJ(dummy, dummy, dummy);
            CartesianVector maxPosJ(dummy, dummy, dummy);

            try
            {
                ClusterList clusters3DJ;
                LArPfoHelper::GetThreeDClusterList(pPfoJ, clusters3DJ);

                const ThreeDSlidingFitResult slidingFitResultJ(clusters3DJ.front(), fitWindow, wirePitch);
                minPosJ = slidingFitResultJ.GetGlobalMinLayerPosition();
                maxPosJ = slidingFitResultJ.GetGlobalMaxLayerPosition();
            }
            catch (const StatusCodeException &)
            {
                // If we can't fit then just skip this pfo
                continue;
            }

            // Get all the different endpoint distances
            float minDist = std::numeric_limits<float>::max();
            for (const auto &posI : {minPosI, maxPosI})
            {
                for (const auto &posJ : {minPosJ, maxPosJ})
                {
                    const auto dist = posI.GetDistanceSquared(posJ);
                    
                    if (dist > minDist)
                        continue;

                    minDist = dist;
                    pfoSeparationMap[pPfoI][pPfoJ] = minDist;
                    pfoSeparationMap[pPfoJ][pPfoI] = minDist;

                    const auto pMatchPoint = new CartesianVector((posI + posJ) * 0.5);
                    const auto pMatchPointCopy = new CartesianVector(*pMatchPoint);
                    pfoMatchPointMap[pPfoI][pPfoJ] = pMatchPoint;
                    pfoMatchPointMap[pPfoJ][pPfoI] = pMatchPointCopy;
                }
            }

            std::advance(iterJ, 1);
        }

        std::advance(iterI, 1);
    }

    // Arrange the PFOs into a crude hierarchy
    auto remainingPfos = allPfos;
    PfoList usedPfos;

    while (!remainingPfos.empty())
    {
        const ParticleFlowObject *pBestChild = nullptr;
        CartesianVector bestMatchPoint(dummy, dummy, dummy);
        
        float minDist = std::numeric_limits<float>::max();

        for (const auto &pPfo : remainingPfos)
        {
            const auto vertexDist = pfoVertexDistMap.at(pPfo);
            const auto separationMap = pfoSeparationMap.at(pPfo);
            const auto matchPointMap = pfoMatchPointMap.at(pPfo);

            if (vertexDist <= minDist)
            {
                minDist = vertexDist;
                pBestChild = pPfo;
                bestMatchPoint = vertexPos;
            }

            for (const auto &pUsedPfo : usedPfos)
            {
                if (separationMap.find(pUsedPfo) == separationMap.end())
                    continue;

                const auto dist = separationMap.at(pUsedPfo);
                if (dist <= minDist)
                {
                    minDist = dist;
                    pBestChild = pPfo;
                    bestMatchPoint = *matchPointMap.at(pUsedPfo);
                }
            }
        }

        // Make the match
        remainingPfos.remove(pBestChild);
        usedPfos.push_back(pBestChild);
        pfoToSeedPointMap.emplace(pBestChild, bestMatchPoint);
    }
    
    // Clean up the mess you just made
    for (auto &entry : pfoMatchPointMap)
    {
        for (const auto &entry2 : entry.second)
        {
            delete entry2.second;
        }

        entry.second.clear();
    }
    */
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::SplitPfo(const ParticleFlowObject *const pPfo, const CartesianVector &vertexPos, PfoList &pfosToSave, PfoList &pfosToDelete) const
{
    std::cout << "Checking if I should split PFO: " << pPfo << std::endl;

    //// BEGIN TEST
    CaloHitList threeDHits;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, threeDHits);
    const auto has3DHits = !threeDHits.empty();

    CartesianVector trackDir(0.f, 0.f, 0.f);
    if (has3DHits)
    {
        CartesianVector centroid(0.f, 0.f, 0.f);
        LArPcaHelper::EigenValues eigenvalues(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors eigenvectors;
        LArPcaHelper::RunPca(threeDHits, centroid, eigenvalues, eigenvectors);
        trackDir = eigenvectors.front();
    }

    if (!has3DHits)
        return;
    //// END TEST

    // Identify viable hits that lie on viable 2D split positions and get the hierarchy of hits
    ViewToHitHierarchyMap viewToHitHierarchyMap;
    ViewToHitsMap viewToCaloHitsMap, viewToSplitHitsMap;
    for (const auto &view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        auto &splitHits = viewToSplitHitsMap[view];

        auto &caloHits = viewToCaloHitsMap[view];
        LArPfoHelper::GetCaloHits(pPfo, view, caloHits);

        HitSeparationMap separationMap;
        this->GetHitSeparationMap(caloHits, separationMap);

        std::vector<CaloHitList> segments;
        this->GetContinuousSegments(vertexPos, caloHits, separationMap, segments);

        auto &hitHierarchy = viewToHitHierarchyMap[view];
        this->BuildHitHierarchy(vertexPos, segments, hitHierarchy);

        //// TEST
        // Filter the segments before finding kinks
        const auto wireDir = CartesianVector(1.f, 0.f, 0.f).GetCrossProduct(LArGeometryHelper::GetWireAxis(this->GetPandora(), view));
        const auto trackWireCosTheta = std::abs(trackDir.GetDotProduct(wireDir));
        
        if (trackWireCosTheta > 0.95)
            continue;

        std::vector<CaloHitList> selectedSegments;
        for (const auto &segment : segments)
        {
            if (segment.size() < 15)
                continue;
       
            /*
            const auto dimension = this->GetDimension(segment);
            std::cout << "      - Dimension: " << dimension << std::endl;
            */

            selectedSegments.push_back(segment);
        }
        //// END TEST

        // Get the hits that sit at a kink or bifurcation position
        this->GetSplitHits(selectedSegments, separationMap, splitHits);
    }

    // Match the 2D split hits between views to find viable a 3D split position
    std::vector<SplitPoint3D> splitPoints3D;
    this->Get3DSplitPoints(viewToCaloHitsMap, viewToSplitHitsMap, splitPoints3D);

    if (splitPoints3D.empty())
        return;

    const auto chosenSplitPoint = this->Select3DSplitPoint(splitPoints3D, viewToCaloHitsMap, viewToHitHierarchyMap);
    
    //// BEGIN DEBUG
    /*
    for (const auto &view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        const auto hits = viewToCaloHitsMap.at(view);
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, "Hits", RED));
    }

    const auto point3D = chosenSplitPoint.m_position3D;
    const auto pointU = LArGeometryHelper::ProjectPosition(this->GetPandora(), point3D, TPC_VIEW_U);
    const auto pointV = LArGeometryHelper::ProjectPosition(this->GetPandora(), point3D, TPC_VIEW_V);
    const auto pointW = LArGeometryHelper::ProjectPosition(this->GetPandora(), point3D, TPC_VIEW_W);

    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pointU, "SplitU", GREEN, 3));
    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pointV, "SplitV", GREEN, 3));
    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pointW, "SplitW", GREEN, 3));

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    */
    //// END DEBUG


    // Add the 3D hits the map
    auto &hits3D = viewToCaloHitsMap[TPC_3D];
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, hits3D);

    // Split the hits in each view upstream and downstream of the split position using the hierarchy
    ViewToHitsMap viewToUpstreamHitsMap, viewToDownstreamHitsMap;
    this->SplitPfoHits(viewToCaloHitsMap, viewToHitHierarchyMap, chosenSplitPoint, viewToUpstreamHitsMap, viewToDownstreamHitsMap);

    //// BEGIN DEBUG
    /*
    for (const auto &view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        const auto upstreamHits = viewToUpstreamHitsMap.at(view);
        const auto downstreamHits = viewToDownstreamHitsMap.at(view);
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &upstreamHits, "Upstream", RED));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &downstreamHits, "Downstream", BLUE));
    }
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    */
    //// END DEBUG

    // Now make the clusters and PFOs for these upstream and downstream segments
    ClusterList upstreamClusters, downstreamClusters;
    this->SplitPfoClusters(pPfo, viewToUpstreamHitsMap, viewToDownstreamHitsMap, upstreamClusters, downstreamClusters);
    this->SplitPfo(pPfo, upstreamClusters, downstreamClusters, pfosToSave, pfosToDelete); 
    std::cout << "Split made!" << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CartesianVector SecondaryInteractionsAlgorithm::GetVertexPosition() const
{
    const VertexList *pVertexList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));

    if (pVertexList->size() != 1)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    return pVertexList->front()->GetPosition();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::CollectInputPfos(PfoList &inputPfos, const std::string &pfoListName) const
{
    const PfoList *pPfoList = nullptr;
    PfoList collectedPfos;

    if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) == STATUS_CODE_SUCCESS)
        collectedPfos.insert(collectedPfos.end(), pPfoList->begin(), pPfoList->end());

    LArPfoHelper::GetAllDownstreamPfos(collectedPfos, inputPfos);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::GetHitSeparationMap(const CaloHitList &caloHitList, HitSeparationMap &separationMap) const
{
    this->GetHitSeparationMap(caloHitList, separationMap, m_isolatedHitDistance);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::GetHitSeparationMap(const CaloHitList &caloHitList, HitSeparationMap &separationMap, const float isolatedHitDistance) const
{
    const auto isolatedHitDistanceSquared = isolatedHitDistance * isolatedHitDistance;

    unsigned int i = 0;
    for (const auto &pHitI : caloHitList)
    {
        unsigned int j = 0;
        for (const auto &pHitJ : caloHitList)
        {
            if (j >= i)
                continue;

            const auto distSquared = (pHitI->GetPositionVector() - pHitJ->GetPositionVector()).GetMagnitudeSquared();

            if (distSquared < isolatedHitDistanceSquared)
            {
                separationMap[pHitI][pHitJ] = distSquared;
                separationMap[pHitJ][pHitI] = distSquared;
            }

            j++;
        }
        i++;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::GetContinuousSegments(const CartesianVector &vertexPos, const CaloHitList &caloHitList, const HitSeparationMap &separationMap, std::vector<CaloHitList> &continuousSegments) const
{
    // We don't need to order 0 or 1 hit
    if (caloHitList.size() < 2)
    {
        continuousSegments.push_back(caloHitList);
        return;
    }

    CaloHitList remainingHits = caloHitList;
    std::vector<CaloHitList> initialSegments;

    while (!remainingHits.empty())
    {
        const auto pSeedHit = this->GetClosestHitToVertex(remainingHits, vertexPos);

        CaloHitList segment;
        this->CollectHitsInSegment(pSeedHit, separationMap, remainingHits, segment);
        initialSegments.push_back(segment);
    }

    this->StitchSegments(initialSegments, continuousSegments);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *SecondaryInteractionsAlgorithm::GetClosestHitToPoint(const CaloHitList &caloHitList, const CartesianVector &point) const
{
    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float minDistanceSquared = std::numeric_limits<float>::max();
    const CaloHit *pClosestHit = nullptr;
    for (const auto &pHit : caloHitList)
    {
        const auto distSquared = pHit->GetPositionVector().GetDistanceSquared(point);
        if (distSquared > minDistanceSquared)
            continue;

        minDistanceSquared = distSquared;
        pClosestHit = pHit;
    }

    if (!pClosestHit)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return pClosestHit;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *SecondaryInteractionsAlgorithm::GetClosestHitToVertex(const CaloHitList &caloHitList, const CartesianVector &vertexPos) const
{
    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // ATTN here we assume that all hits in the input list are of the same type as the first
    const auto hitType = caloHitList.front()->GetHitType();
    const auto is2D = (hitType == TPC_VIEW_U || hitType == TPC_VIEW_V || hitType == TPC_VIEW_W);
    const auto is3D = (hitType == TPC_3D);
    
    if (!is2D && !is3D)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Project the 3D vertex position into the desired view - if applicable
    const auto vertexPosProjected = (is2D ? LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPos, hitType) : vertexPos);

    return this->GetClosestHitToPoint(caloHitList, vertexPosProjected);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::CollectHitsInSegment(const CaloHit *const pSeedHit, const HitSeparationMap &separationMap, CaloHitList &remainingHits, CaloHitList &segment) const
{
    const CaloHit *pLastHit = pSeedHit;
    const CaloHit *pCurrentHit = pSeedHit;

    while (true)
    {
        segment.push_back(pCurrentHit);
        remainingHits.remove(pCurrentHit);

        const CaloHit *pNextHit = nullptr;
        if (!this->GetNearestNeighbor(separationMap, remainingHits, pLastHit, pCurrentHit, pNextHit))
            break;

        pCurrentHit = pNextHit;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SecondaryInteractionsAlgorithm::GetNearestNeighbor(const HitSeparationMap &separationMap, const CaloHitList &remainingHits, const CaloHit *const pLastHit, const CaloHit *const pCurrentHit, const CaloHit *&pNextHit) const
{
    // Make sure we have an entry for the current hit in the separation map
    const auto iter = separationMap.find(pCurrentHit);
    if (iter == separationMap.end())
        return false;

    // Get the direction unit vector from the last to the current hit
    const auto diff = (pCurrentHit->GetPositionVector() - pLastHit->GetPositionVector());
    CartesianVector dir(0.f, 0.f, 1.f);
    bool hasDir = false;
    if (diff.GetMagnitudeSquared() > std::pow(std::numeric_limits<float>::epsilon(), 2))
    {
        hasDir = true;
        dir = diff.GetUnitVector();
    }
    
    // Find the neighbor hit with the minimum score
    bool foundNeighbor = false;
    float minScore = std::numeric_limits<float>::max();
    for (const auto &entry : iter->second)
    {
        const auto pHit = entry.first;
        const auto distanceSquared = entry.second;

        if (std::find(remainingHits.begin(), remainingHits.end(), pHit) == remainingHits.end())
            continue;

        auto score = distanceSquared;
        if (hasDir)
        {
            const auto d = pHit->GetPositionVector() - pCurrentHit->GetPositionVector();
            const auto longitudinalDist = d.GetDotProduct(dir);
            const auto transverseDist = d.GetMagnitudeSquared() - longitudinalDist * longitudinalDist;

            score = m_transverseBias*transverseDist*transverseDist + longitudinalDist*longitudinalDist;
        }

        if (score > minScore)
            continue;

        foundNeighbor = true;
        minScore = score;
        pNextHit = pHit;
    }

    return foundNeighbor;
}
        
// -----------------------------------------------------------------------------------------------------------------------------------------
    
void SecondaryInteractionsAlgorithm::StitchSegments(const std::vector<CaloHitList> &initialSegments, std::vector<CaloHitList> &stitchedSegments) const
{
    stitchedSegments = initialSegments;
    bool stitchMade = true;
    const auto thresholdSquared = m_stitchingThreshold * m_stitchingThreshold;

    while (stitchMade)
    {
        stitchMade = false;
        unsigned int stitchIndexI = std::numeric_limits<unsigned int>::max();
        unsigned int stitchIndexJ = std::numeric_limits<unsigned int>::max();
        bool stitchStartI = false;
        bool stitchStartJ = false;
        float minDistSquared = std::numeric_limits<float>::max();

        // Loop over all possible pairs of segments
        for (unsigned int i = 0; i < stitchedSegments.size(); ++i)
        {
            const auto segmentI = stitchedSegments.at(i);
            const auto startI = segmentI.front()->GetPositionVector();
            const auto endI = segmentI.back()->GetPositionVector();

            for (unsigned int j = 0; j < i; ++j)
            {
                const auto segmentJ = stitchedSegments.at(j);
                const auto startJ = segmentJ.front()->GetPositionVector();
                const auto endJ = segmentJ.back()->GetPositionVector();

                // Loop over all combinations of start & end points (SS, SE, ES, EE)
                for (const auto &startEndI : {true, false})
                {
                    const auto pointI = startEndI ? startI : endI;
                    for (const auto &startEndJ : {true, false})
                    {
                        const auto pointJ = startEndJ ? startJ : endJ;

                        // Determine if the endpoint distance is below the threshold and the best yet
                        const auto distSquared = pointI.GetDistanceSquared(pointJ);
                        if (distSquared < thresholdSquared && distSquared < minDistSquared)
                        {
                            // Save this match for future reference
                            stitchMade = true;
                            minDistSquared = distSquared;
                            stitchIndexI = i;
                            stitchIndexJ = j;
                            stitchStartI = startEndI;
                            stitchStartJ = startEndJ;
                        }
                    }
                }
            }
        }

        // Make the stitch
        if (stitchMade)
        {
            auto &segmentI = stitchedSegments.at(stitchIndexI);
            auto &segmentJ = stitchedSegments.at(stitchIndexJ);

            // For SS or EE matches, we need to flip the second segment
            const auto shouldFlipJ = (stitchStartI == stitchStartJ);
            if (shouldFlipJ)
                segmentJ.reverse();

            // For SS or SE matches, we should add the second segment in front of the first - for ES or EE we should add it at the back
            const auto shouldPutIFirst = !stitchStartI;
            const auto insertPos = shouldPutIFirst ? segmentI.end() : segmentI.begin();

            // Combine the second segment into the first
            segmentI.insert(insertPos, segmentJ.begin(), segmentJ.end());

            // Remove the second segment
            stitchedSegments.erase(std::next(stitchedSegments.begin(), stitchIndexJ));
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::BuildHitHierarchy(const CartesianVector &vertexPos, const std::vector<CaloHitList> &continuousSegments, CaloHitHierarchyMap &hitHierarchy) const
{
    // Collect all of the hits into a single list
    CaloHitList allHits;
    for (const auto &segment : continuousSegments)
        allHits.insert(allHits.end(), segment.begin(), segment.end());

    if (allHits.empty())
        return;
    
    // Build a separation map with no thresholding
    HitSeparationMap separationMap;
    this->GetHitSeparationMap(allHits, separationMap, std::numeric_limits<float>::max());

    // Make the initial hierarchy map with the segment that is closest to the vertex and mark all others as orphans
    std::vector<CaloHitList> orphanSegments;
    const auto pSeedHit = this->GetClosestHitToVertex(allHits, vertexPos);
    for (const auto &segment : continuousSegments)
    {
        // If the segment doesn't contain the seed hit call it an orphan
        if (std::find(segment.begin(), segment.end(), pSeedHit) == segment.end())
        {
            orphanSegments.push_back(segment);
        }
        else
        {
            this->BuildInitialHitHierarchy(segment, pSeedHit, hitHierarchy);
        }
    }
    
    // Keep making the best possible match until there are no more orphans
    while (!orphanSegments.empty())
    {
        this->MakeBestHierarchyLink(separationMap, orphanSegments, hitHierarchy);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::BuildInitialHitHierarchy(const CaloHitList &segment, const CaloHit *const pSeedHit, CaloHitHierarchyMap &hitHierarchy) const
{
    if (!hitHierarchy.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (segment.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const auto seedHitIter = std::find(segment.begin(), segment.end(), pSeedHit);
    if (seedHitIter == segment.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // If the seed hit is the only one in the input list, then it has no daughters
    if (segment.size() == 1)
    {
        hitHierarchy.emplace(*seedHitIter, std::vector<const CaloHit *>());
        return;
    }
   
    // Move forward through the list from the seed hit
    auto hitIter = seedHitIter;
    while (std::next(hitIter) != segment.end())
    {
        hitHierarchy[*hitIter].push_back(*std::next(hitIter));
        std::advance(hitIter, 1);
    }

    // Move backward through the list from the seed hit
    hitIter = seedHitIter;
    while (hitIter != segment.begin())
    {
        hitHierarchy[*hitIter].push_back(*std::prev(hitIter));
        std::advance(hitIter, -1);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::MakeBestHierarchyLink(const HitSeparationMap &separationMap, std::vector<CaloHitList> &orphanSegments, CaloHitHierarchyMap &hitHierarchy) const
{
    if (orphanSegments.empty())
        return;

    // Extract all of the hits from the hit hierarchy so far
    CaloHitList allHitsInHierarchy;
    for (const auto &entry : hitHierarchy)
    {
        // Parent hits
        const auto pParentHit = entry.first;
        if (std::find(allHitsInHierarchy.begin(), allHitsInHierarchy.end(), pParentHit) == allHitsInHierarchy.end())
            allHitsInHierarchy.push_back(pParentHit);

        // Daughter hits
        for (const auto &pDaughterHit : entry.second)
        {
            if (std::find(allHitsInHierarchy.begin(), allHitsInHierarchy.end(), pDaughterHit) == allHitsInHierarchy.end())
                allHitsInHierarchy.push_back(pDaughterHit);
        }
    }

    if (allHitsInHierarchy.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float minHitDistSquared = std::numeric_limits<float>::max();
    unsigned int nearestSegmentIndex = 0;
    const CaloHit *pNearestHit = nullptr;
    bool matchFront = true;

    for (unsigned int i = 0; i < orphanSegments.size(); i++)
    {
        const auto segment = orphanSegments.at(i);

        // Check the first hit in the segment for a neighbor
        const auto pFrontHit = segment.front();
        const auto pFrontNextHit = *std::next(segment.begin(), segment.size() == 1 ? 0 : 1);

        const CaloHit *pNearestFrontHit = nullptr;
        if (!this->GetNearestNeighbor(separationMap, allHitsInHierarchy, pFrontNextHit, pFrontHit, pNearestFrontHit))
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const auto frontHitDistSquared = pFrontHit->GetPositionVector().GetDistanceSquared(pNearestFrontHit->GetPositionVector());

        // Check the last hit in the segment for a neighbor
        const auto pBackHit = segment.back();
        const auto pBackPrevHit = *std::prev(segment.end(), segment.size() == 1 ? 1 : 2); // ATTN here we go back 1 or 2 since end() points to the position after the back hit

        const CaloHit *pNearestBackHit = nullptr;
        if (!this->GetNearestNeighbor(separationMap, allHitsInHierarchy, pBackPrevHit, pBackHit, pNearestBackHit))
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        
        const auto backHitDistSquared = pBackHit->GetPositionVector().GetDistanceSquared(pNearestBackHit->GetPositionVector());

        // Get the closest of the two hits
        const auto nearestHitDistSquared = std::min(frontHitDistSquared, backHitDistSquared);

        // Check if this is the closest so far
        if (nearestHitDistSquared < minHitDistSquared)
        {
            minHitDistSquared = nearestHitDistSquared;
            nearestSegmentIndex = i;
            matchFront = (frontHitDistSquared < backHitDistSquared);
            pNearestHit = matchFront ? pNearestFrontHit : pNearestBackHit;
        }
    }

    // Make the link
    this->AddSegmentToHitHierarchy(nearestSegmentIndex, matchFront, pNearestHit, orphanSegments, hitHierarchy);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::AddSegmentToHitHierarchy(const unsigned int nearestSegmentIndex, const bool matchFront, const CaloHit *const pNearestHit, std::vector<CaloHitList> &orphanSegments, CaloHitHierarchyMap &hitHierarchy) const
{
    auto &nearestHitDaughters = hitHierarchy[pNearestHit];
    const auto segment = orphanSegments.at(nearestSegmentIndex);

    if (matchFront)
    {
        // Move forward through the list from the seed hit
        nearestHitDaughters.push_back(segment.front());
        auto hitIter = segment.begin();
        while (std::next(hitIter) != segment.end())
        {
            hitHierarchy[*hitIter].push_back(*std::next(hitIter));
            std::advance(hitIter, 1);
        }
    }
    else
    {
        // Move backward through the list from the seed hit
        nearestHitDaughters.push_back(segment.back());
        auto hitIter = segment.end();
        while (hitIter != segment.begin())
        {
            hitHierarchy[*hitIter].push_back(*std::prev(hitIter));
            std::advance(hitIter, -1);
        }
    }

    orphanSegments.erase(std::next(orphanSegments.begin(), nearestSegmentIndex));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::PrintHitHierarchyMap(const CaloHitHierarchyMap &hitHierarchy, const CaloHit *const pSeedHit, const unsigned int depth) const
{
    std::cout << std::string(depth * 4, ' ') << " - " << pSeedHit;

    const auto iter = hitHierarchy.find(pSeedHit);
    if (iter != hitHierarchy.end())
    {
        const auto daughters = iter->second;
        const auto nDaughters = daughters.size();
        std::cout << " (" << nDaughters << ") ";

        if (nDaughters > 1)
        {
            std::cout << " -> ";
            for (const auto pDaughterHit : daughters)
                std::cout << pDaughterHit << "  ";
        }
        std::cout << std::endl;

        for (const auto pDaughterHit : daughters)
        {
            this->PrintHitHierarchyMap(hitHierarchy, pDaughterHit, depth + ((nDaughters == 1) ? 0 : 1));
        }
    }
    else
    {
        std::cout << " (X) " << std::endl << std::endl;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::GetDownstreamHits(const CaloHitHierarchyMap &hitHierarchy, const CaloHit *const pSeedHit, CaloHitList &downstreamHits) const
{
    downstreamHits.push_back(pSeedHit);

    const auto iter = hitHierarchy.find(pSeedHit);
    if (iter != hitHierarchy.end())
    {
        const auto daughters = iter->second;
        for (const auto &pDaughterHit : daughters)
        {
            this->GetDownstreamHits(hitHierarchy, pDaughterHit, downstreamHits);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::GetSplitHits(const std::vector<CaloHitList> &continuousSegments, const HitSeparationMap &separationMap, CaloHitList &splitHits) const
{
    if (continuousSegments.empty())
        return;

    // Collect hits at possible split points: kinks, bifurcations
    CaloHitList collectedHits;
    this->GetBifurcationHits(continuousSegments, separationMap, collectedHits);

    for (const auto &segment : continuousSegments)
        this->GetKinkHits(segment, collectedHits);

    // Add the hits to the output provided they are unique
    for (const auto &pHit : collectedHits)
    {
        if (std::find(splitHits.begin(), splitHits.end(), pHit) == splitHits.end())
            splitHits.push_back(pHit);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::GetKinkHits(const CaloHitList &segment, CaloHitList &kinkHits) const
{
    if (segment.size() < 1 + 2 * m_nSampleHits)
        return;

    // Hits before the current hit
    auto preBegin = segment.begin();
    auto preEnd = std::next(preBegin, m_nSampleHits);
    
    // Current hit
    auto thisHit = std::next(segment.begin(), m_nSampleHits);
        
    // Hits after the current hit
    auto postBegin = std::next(preEnd, 1);
    auto postEnd = std::next(postBegin, m_nSampleHits);

    // Get the hit indices with a kink above the threshold angle
    std::vector<std::pair<unsigned int, float> > hitIndexCosThetaVector;
    for (unsigned int i = m_nSampleHits; i < segment.size() - m_nSampleHits; ++i)
    {
        // Get the hits from the iterators
        const auto thisHitPos = (*thisHit)->GetPositionVector();
        const CaloHitList preHits(preBegin, preEnd);
        const CaloHitList postHits(postBegin, postEnd);

        const auto cosTheta = this->GetKinkAngle(thisHitPos, preHits, postHits);

        if (cosTheta < m_cosAngleThreshold)
            hitIndexCosThetaVector.emplace_back(i, cosTheta);

        preBegin++;
        preEnd++;
        postBegin++;
        postEnd++;
        thisHit++;
    }

    if (hitIndexCosThetaVector.empty())
        return;

    // When we have a bunch of consecutive indices passing the threshold, choose the index with the largest kink angles
    std::vector<std::pair<unsigned int, float> > hitIndexCosThetaBunch = {hitIndexCosThetaVector.front()};
    unsigned int lastIndex = hitIndexCosThetaBunch.front().first;
    for (unsigned int i = 1; i < hitIndexCosThetaVector.size(); ++i)
    {
        const auto entry = hitIndexCosThetaVector.at(i);
        const auto index = entry.first;

        if (index != lastIndex + 1)
        {
            const auto kinkIndex = this->GetIndexWithMaxKinkAngle(hitIndexCosThetaBunch);
            kinkHits.push_back(*(std::next(segment.begin(), kinkIndex)));
            hitIndexCosThetaBunch.clear();
        }

        hitIndexCosThetaBunch.push_back(entry);
        lastIndex = index;
    }

    const auto kinkIndex = this->GetIndexWithMaxKinkAngle(hitIndexCosThetaBunch);
    kinkHits.push_back(*(std::next(segment.begin(), kinkIndex)));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SecondaryInteractionsAlgorithm::GetKinkAngle(const CartesianVector &thisHitPos, const CaloHitList &preHits, const CaloHitList &postHits) const
{
    CartesianVector preCentroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenValues preEigenvalues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors preEigenvectors;
    LArPcaHelper::RunPca(preHits, preCentroid, preEigenvalues, preEigenvectors);
    auto preDir = preEigenvectors.front();

    CartesianVector postCentroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenValues postEigenvalues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors postEigenvectors;
    LArPcaHelper::RunPca(postHits, postCentroid, postEigenvalues, postEigenvectors);
    auto postDir = postEigenvectors.front();

    // Flip the directions so they point along the direction of the "track"
    preDir *= ((thisHitPos - preCentroid).GetUnitVector().GetDotProduct(preDir) > 0 ? 1 : -1);
    postDir *= ((thisHitPos - postCentroid).GetUnitVector().GetDotProduct(postDir) > 0 ? -1 : 1);

    return (preDir).GetCosOpeningAngle(postDir);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
unsigned int SecondaryInteractionsAlgorithm::GetIndexWithMaxKinkAngle(const std::vector<std::pair<unsigned int, float> > &hitIndexCosThetaBunch) const
{
    if (hitIndexCosThetaBunch.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    unsigned int outputIndex = std::numeric_limits<unsigned int>::max();
    float minCosTheta = std::numeric_limits<float>::max();

    for (const auto &entry : hitIndexCosThetaBunch)
    {
        const auto index = entry.first;
        const auto cosTheta = entry.second;

        if (cosTheta > minCosTheta)
            continue;

        outputIndex = index;
        minCosTheta = cosTheta;
    }

    return outputIndex;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::GetBifurcationHits(const std::vector<CaloHitList> &continuousSegments, const HitSeparationMap &separationMap, CaloHitList &bifurcationHits) const
{
    if (continuousSegments.size() < 2)
        return;

    for (unsigned int i = 0; i < continuousSegments.size(); ++i)
    {
        const auto segmentI = continuousSegments.at(i);

        if (segmentI.empty())
            continue;

        // Collect all of the hits in the other segments (ignoring the first and last hits)
        CaloHitList remainingHits;
        for (unsigned int j = 0; j < continuousSegments.size(); ++j)
        {
            if (i == j)
                continue;
        
            const auto segmentJ = continuousSegments.at(j);
            if (segmentJ.size() <= 2)
                continue;
            
            remainingHits.insert(remainingHits.end(), std::next(segmentJ.begin()), std::prev(segmentJ.end()));
        }

        // Check the first hit in the segment for a neighbor
        const auto pFrontHit = segmentI.front();
        const auto pFrontNextHit = *std::next(segmentI.begin(), segmentI.size() == 1 ? 0 : 1);
        
        const CaloHit *pBifurcationHitFront = nullptr;
        if (this->GetNearestNeighbor(separationMap, remainingHits, pFrontNextHit, pFrontHit, pBifurcationHitFront))
            bifurcationHits.push_back(pBifurcationHitFront);
    
        // Check the last hit in the segment for a neighbor
        const auto pBackHit = segmentI.back();
        const auto pBackPrevHit = *std::prev(segmentI.end(), segmentI.size() == 1 ? 1 : 2); // ATTN here we go back 1 or 2 since end() points to the position after the back hit
        
        const CaloHit *pBifurcationHitBack = nullptr;
        if (this->GetNearestNeighbor(separationMap, remainingHits, pBackPrevHit, pBackHit, pBifurcationHitBack))
            bifurcationHits.push_back(pBifurcationHitBack);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
            
void SecondaryInteractionsAlgorithm::Get3DSplitPoints(const ViewToHitsMap &viewToAllHitsMap, const ViewToHitsMap &viewToSplitHitsMap, std::vector<SplitPoint3D> &splitPoints3D) const
{
    const auto &allHitsU = viewToAllHitsMap.at(TPC_VIEW_U);
    const auto &allHitsV = viewToAllHitsMap.at(TPC_VIEW_V);
    const auto &allHitsW = viewToAllHitsMap.at(TPC_VIEW_W);

    auto splitHitsU = viewToSplitHitsMap.at(TPC_VIEW_U);
    auto splitHitsV = viewToSplitHitsMap.at(TPC_VIEW_V);
    auto splitHitsW = viewToSplitHitsMap.at(TPC_VIEW_W);

    while (this->MakeThreeViewMatch(splitHitsU, splitHitsV, splitHitsW, splitPoints3D)) {}
    while (this->MakeTwoViewMatch(allHitsU, allHitsV, allHitsW, splitHitsU, splitHitsV, splitHitsW, splitPoints3D)) {}
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
SecondaryInteractionsAlgorithm::SplitPoint3D SecondaryInteractionsAlgorithm::Select3DSplitPoint(const std::vector<SplitPoint3D> &splitPoints3D, const ViewToHitsMap &viewToAllHitsMap, const ViewToHitHierarchyMap &viewToHitHierarchyMap) const
{
    if (splitPoints3D.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (splitPoints3D.size() == 1)
        return splitPoints3D.front();

    auto bestSplitPoint = splitPoints3D.front();
    float bestScore = -std::numeric_limits<float>::max();

    for (const auto &splitPoint : splitPoints3D)
    {
        const auto score = this->GetSplitPointScore(splitPoint, viewToAllHitsMap, viewToHitHierarchyMap);
        if (score > bestScore)
        {
            bestScore = score;
            bestSplitPoint = splitPoint;
        }
    }

    return bestSplitPoint;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SecondaryInteractionsAlgorithm::GetSplitPointScore(const SplitPoint3D &splitPoint, const ViewToHitsMap &viewToAllHitsMap, const ViewToHitHierarchyMap &viewToHitHierarchyMap) const
{
    float score = 0.f;
    unsigned int nViewsUsed = 0;
    for (const auto &view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        CaloHitList upstreamHits, downstreamHits;
        this->SplitHitListAtPoint(view, viewToAllHitsMap.at(view), viewToHitHierarchyMap.at(view), splitPoint, upstreamHits, downstreamHits);

        const auto nUpstreamHits = upstreamHits.size();
        const auto nDownstreamHits = downstreamHits.size();
        const auto nHitsTotal = nUpstreamHits + nDownstreamHits;

        if (nHitsTotal == 0)
            continue;
    
        // Give more upstream splits a larger score
        const auto upstreamScore = static_cast<float>(nDownstreamHits) / static_cast<float>(nHitsTotal);

        score += upstreamScore;
        nViewsUsed++;
    }

    // Input PFO has no hits at all
    if (nViewsUsed == 0)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Normalize the score to unity
    score *= 1.f / static_cast<float>(nViewsUsed);

    return score;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SecondaryInteractionsAlgorithm::MakeThreeViewMatch(CaloHitList &hitsU, CaloHitList &hitsV, CaloHitList &hitsW, std::vector<SplitPoint3D> &splitPoints3D) const
{
    if (hitsU.empty() || hitsV.empty() || hitsW.empty())
        return false;

    /*
    std::cout << "Looking for 3-view matches for split hits" << std::endl;
    std::cout << "  - NHitsU   : " << hitsU.size() << std::endl;
    std::cout << "  - NHitsV   : " << hitsV.size() << std::endl;
    std::cout << "  - NHitsW   : " << hitsW.size() << std::endl;
    std::cout << "  - NMatches : " << splitPoints3D.size() << std::endl;
    */

    // Collect the viable split points
    std::vector<SplitPoint3D> viableSplitPoints;
    for (const auto &pHitU : hitsU)
    {
        for (const auto &pHitV : hitsV)
        {
            for (const auto &pHitW : hitsW)
            {
                const SplitPoint3D splitPoint(this->GetPandora(), pHitU, pHitV, pHitW);

                /*
                std::cout << "Considering triplet" << std::endl;
                std::cout << "  - U    : " << pHitU << std::endl;
                std::cout << "  - V    : " << pHitV << std::endl;
                std::cout << "  - W    : " << pHitW << std::endl;
                std::cout << "  - dX   : " << splitPoint.m_maxDeltaX << std::endl;
                std::cout << "  - chi2 : " << splitPoint.m_chi2 << std::endl;
                */

                if (splitPoint.m_maxDeltaX > m_maxMatchDeltaX)
                    continue;
                
                //std::cout << "  - Passes dX cut!" << std::endl;
                
                if (splitPoint.m_chi2 > m_maxMatch3ViewChi2)
                    continue;
                
                //std::cout << "  - Passes chi2 cut!" << std::endl;

                viableSplitPoints.push_back(splitPoint);
            }
        }
    }
                
    //std::cout << "Found " << viableSplitPoints.size() << " viable triplets" << std::endl;

    if (viableSplitPoints.empty())
        return false;

    // Get the best split point
    unsigned int bestSplitPointIndex = 0;
    float bestChi2 = std::numeric_limits<float>::max();
    for (unsigned int i = 0; i < viableSplitPoints.size(); ++i)
    {
        const auto chi2 = viableSplitPoints.at(i).m_chi2;
        if (chi2 < bestChi2)
        {
            bestChi2 = chi2;
            bestSplitPointIndex = i;
        }
    }

    // Make the match
    const auto bestSplitPoint = viableSplitPoints.at(bestSplitPointIndex);
    
    /*
    std::cout << "Best triplet is" << std::endl;
    std::cout << "  - U : " << bestSplitPoint.GetHitWithView(TPC_VIEW_U) << std::endl;
    std::cout << "  - V : " << bestSplitPoint.GetHitWithView(TPC_VIEW_V) << std::endl;
    std::cout << "  - W : " << bestSplitPoint.GetHitWithView(TPC_VIEW_W) << std::endl;
    */

    splitPoints3D.push_back(bestSplitPoint);
    hitsU.remove(bestSplitPoint.GetHitWithView(TPC_VIEW_U));
    hitsV.remove(bestSplitPoint.GetHitWithView(TPC_VIEW_V));
    hitsW.remove(bestSplitPoint.GetHitWithView(TPC_VIEW_W));

    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SecondaryInteractionsAlgorithm::MakeTwoViewMatch(const CaloHitList &allHitsU, const CaloHitList &allHitsV, const CaloHitList &allHitsW, CaloHitList &hitsU, CaloHitList &hitsV, CaloHitList &hitsW, std::vector<SplitPoint3D> &splitPoints3D) const
{
    const auto nViewsEmpty = (hitsU.empty() ? 1 : 0) + (hitsV.empty() ? 1 : 0) + (hitsW.empty() ? 1 : 0);
    if (nViewsEmpty > 1)
        return false;

    std::vector<SplitPoint3D> viableSplitPoints;
    this->GetViableTwoViewMatches(hitsU, hitsV, allHitsW, viableSplitPoints);
    this->GetViableTwoViewMatches(hitsU, hitsW, allHitsV, viableSplitPoints);
    this->GetViableTwoViewMatches(hitsV, hitsW, allHitsU, viableSplitPoints);
    
    //std::cout << "Found " << viableSplitPoints.size() << " viable doublets" << std::endl;

    if (viableSplitPoints.empty())
        return false;

    // Get the best split point
    unsigned int bestSplitPointIndex = 0;
    float bestChi2 = std::numeric_limits<float>::max();
    for (unsigned int i = 0; i < viableSplitPoints.size(); ++i)
    {
        const auto chi2 = viableSplitPoints.at(i).m_chi2;
        if (chi2 < bestChi2)
        {
            bestChi2 = chi2;
            bestSplitPointIndex = i;
        }
    }

    // Make the match
    const auto bestSplitPoint = viableSplitPoints.at(bestSplitPointIndex);
    
    /*
    std::cout << "Best doublet is" << std::endl;
    std::cout << "  - A : " << bestSplitPoint.m_pHitA << std::endl;
    std::cout << "  - B : " << bestSplitPoint.m_pHitB << std::endl;
    */

    splitPoints3D.push_back(bestSplitPoint);

    const auto hitTypeA = bestSplitPoint.m_pHitA->GetHitType();
    const auto hitTypeB = bestSplitPoint.m_pHitB->GetHitType();

    if (hitTypeA == TPC_VIEW_U || hitTypeB == TPC_VIEW_U)
        hitsU.remove(bestSplitPoint.GetHitWithView(TPC_VIEW_U));
    
    if (hitTypeA == TPC_VIEW_V || hitTypeB == TPC_VIEW_V)
        hitsV.remove(bestSplitPoint.GetHitWithView(TPC_VIEW_V));
    
    if (hitTypeA == TPC_VIEW_W || hitTypeB == TPC_VIEW_W)
        hitsW.remove(bestSplitPoint.GetHitWithView(TPC_VIEW_W));

    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::GetViableTwoViewMatches(const CaloHitList &hitsA, const CaloHitList &hitsB, const CaloHitList &hitsC, std::vector<SplitPoint3D> &splitPointsAB) const
{
    for (const auto &pHitA : hitsA)
    {
        for (const auto &pHitB : hitsB)
        {
            // Get the remaining hit type
            const auto hitTypeA = pHitA->GetHitType();
            const auto hitTypeB = pHitB->GetHitType();
            HitType hitTypeC = HIT_CUSTOM;

            if ((hitTypeA == TPC_VIEW_U && hitTypeB == TPC_VIEW_V) || (hitTypeA == TPC_VIEW_V && hitTypeB == TPC_VIEW_U))
            {
                hitTypeC = TPC_VIEW_W;
            }
            else if ((hitTypeA == TPC_VIEW_U && hitTypeB == TPC_VIEW_W) || (hitTypeA == TPC_VIEW_W && hitTypeB == TPC_VIEW_U))
            {
                hitTypeC = TPC_VIEW_V;
            }
            else if ((hitTypeA == TPC_VIEW_V && hitTypeB == TPC_VIEW_W) || (hitTypeA == TPC_VIEW_W && hitTypeB == TPC_VIEW_V))
            {
                hitTypeC = TPC_VIEW_U;
            }
            else
            {
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
            }

            const SplitPoint3D splitPoint(this->GetPandora(), pHitA, pHitB);
            
            /*
            std::cout << "Considering doublet" << std::endl;
            std::cout << "  - A    : " << pHitA << " - " << hitTypeA << std::endl;
            std::cout << "  - B    : " << pHitB << " - " << hitTypeB << std::endl;
            std::cout << "  - dX   : " << splitPoint.m_maxDeltaX << std::endl;
            std::cout << "  - chi2 : " << splitPoint.m_chi2 << std::endl;
            */

            if (splitPoint.m_maxDeltaX > m_maxMatchDeltaX)
                continue;

            //std::cout << "  - Passes dX cut!" << std::endl;

            // Get the 3D position as projected into the remaining view
            const auto projectedPosition = LArGeometryHelper::ProjectPosition(this->GetPandora(), splitPoint.m_position3D, hitTypeC);
            const auto isInGap = LArGeometryHelper::IsInGap(this->GetPandora(), projectedPosition, hitTypeC, m_twoViewProjectionThreshold);

            //if (isInGap) std::cout << "  - Projects to gap!" << std::endl;

            bool projectsToHit = false;
            if (!hitsC.empty())
            {
                const auto pNearestHit = this->GetClosestHitToPoint(hitsC, projectedPosition);
                const auto nearestHitDistSquared = projectedPosition.GetDistanceSquared(pNearestHit->GetPositionVector());
                //std::cout << "  - Nearest hit distance (squared) : " << nearestHitDistSquared << std::endl;
                projectsToHit = (nearestHitDistSquared < m_twoViewProjectionThreshold * m_twoViewProjectionThreshold);
            }
            else
            {
                //std::cout << "  - No hits in remaining view!" << std::endl;
            }

            //if (projectsToHit) std::cout << "  - Projects to hit! (not using)" << std::endl;

            if (!isInGap && !projectsToHit)
                continue;
                
            splitPointsAB.push_back(splitPoint);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::SplitHitListAtPoint(const HitType view, const CaloHitList &caloHits, const CaloHitHierarchyMap &hitHierarchyMap, const SplitPoint3D &chosenSplitPoint, CaloHitList &upstreamHits, CaloHitList &downstreamHits) const
{
    // In the case of only one hit, just call it upstream
    if (caloHits.size() <= 1)
    {
        upstreamHits = caloHits;
        downstreamHits.clear();
        return;
    }

    // Get the split hit if possible
    bool hasSplitHit = true;
    const pandora::CaloHit *pSplitHit = nullptr;
    try
    {
        pSplitHit = chosenSplitPoint.GetHitWithView(view);
    }
    catch (const StatusCodeException &) { hasSplitHit = false; }

    if (!hasSplitHit)
    {
        // If the split position projects to a gap, then just find the closest hit 
        const auto projectedPosition = LArGeometryHelper::ProjectPosition(this->GetPandora(), chosenSplitPoint.m_position3D, view);
        pSplitHit = this->GetClosestHitToPoint(caloHits, projectedPosition);
    }

    //// BEGIN DEBUG
    /*
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, "Hits", AUTO));
    for (const auto &pCaloHit : caloHits)
    {
        const auto iter = hitHierarchyMap.find(pCaloHit);
        if (iter == hitHierarchyMap.end())
            continue;

        const auto pos = pCaloHit->GetPositionVector();
        const auto &daughters = iter->second;
        for (const auto &pDaughterHit : daughters)
        {
            const auto daughterPos = pDaughterHit->GetPositionVector();
            
            const auto arrowA = daughterPos + ((daughterPos - pos).GetUnitVector().GetCrossProduct(CartesianVector(0.f, 1.f, 0.f)) * (+1) - (daughterPos - pos).GetUnitVector()) * 0.1;
            const auto arrowB = daughterPos + ((daughterPos - pos).GetUnitVector().GetCrossProduct(CartesianVector(0.f, 1.f, 0.f)) * (-1) - (daughterPos - pos).GetUnitVector()) * 0.1;
            
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pos, &daughterPos, "", RED, 1, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &daughterPos, &arrowA, "", RED, 1, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &daughterPos, &arrowB, "", RED, 1, 1));
        }
    }
    const auto splitHitPos = pSplitHit->GetPositionVector();
    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &splitHitPos, "", BLUE, 2));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    */
    //// END DEBUG
    
    this->GetDownstreamHits(hitHierarchyMap, pSplitHit, downstreamHits);
    for (const auto &pCaloHit : caloHits)
    {
        if (std::find(downstreamHits.begin(), downstreamHits.end(), pCaloHit) == downstreamHits.end())
            upstreamHits.push_back(pCaloHit);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::SplitPfoHits(const ViewToHitsMap &viewToCaloHitsMap, const ViewToHitHierarchyMap &viewToHitHierarchyMap, const SplitPoint3D &chosenSplitPoint, ViewToHitsMap &viewToUpstreamHitsMap, ViewToHitsMap &viewToDownstreamHitsMap) const
{
    // Handle the 2D hits
    for (const auto &view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        auto &upstreamHits = viewToUpstreamHitsMap[view];
        auto &downstreamHits = viewToDownstreamHitsMap[view];
        const auto &caloHits = viewToCaloHitsMap.at(view);
        const auto &hitHierarchyMap = viewToHitHierarchyMap.at(view);

        this->SplitHitListAtPoint(view, caloHits, hitHierarchyMap, chosenSplitPoint, upstreamHits, downstreamHits);
        
        // TODO make sure that the upstream & downstream definitions are consistent between views
    }

    // Handle the 3D hits
    auto &upstreamHits = viewToUpstreamHitsMap[TPC_3D];
    auto &downstreamHits = viewToDownstreamHitsMap[TPC_3D];
    const auto &caloHits = viewToCaloHitsMap.at(TPC_3D);
    
    for (const auto &pHit : caloHits)
    {
        // ATTN the parent of a 3D hit is the 2D hit from which it was produced
        const auto pParent2DHit = static_cast<const CaloHit *>(pHit->GetParentAddress());

        // Determine if the parent hit is in the upstream or downstream list
        const auto view = pParent2DHit->GetHitType();
        const auto &upstreamHitsInView = viewToUpstreamHitsMap.at(view);
        const auto &downstreamHitsInView = viewToDownstreamHitsMap.at(view);
        const auto iterUpstream = std::find(upstreamHitsInView.begin(), upstreamHitsInView.end(), pParent2DHit);
        const auto iterDownstream = std::find(downstreamHitsInView.begin(), downstreamHitsInView.end(), pParent2DHit);
        const auto isUpstream = (iterUpstream != upstreamHitsInView.end());
        const auto isDownstream = (iterDownstream != downstreamHitsInView.end());

        if (isUpstream && isDownstream)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // Add the 3D hit to the same list as it's parent
        if (isUpstream)
            upstreamHits.push_back(pHit);
        else
            downstreamHits.push_back(pHit);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------    
    
void SecondaryInteractionsAlgorithm::SplitPfoClusters(const ParticleFlowObject *const pPfo, const ViewToHitsMap &viewToUpstreamHitsMap, const ViewToHitsMap &viewToDownstreamHitsMap, ClusterList &upstreamClusters, ClusterList &downstreamClusters) const
{
    for (const auto &view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, TPC_3D})
    {
        // Get the original clusters and remove them from the PFO
        ClusterList clusters;
        LArPfoHelper::GetClusters(pPfo, view, clusters);

        if (clusters.empty())
            continue;

        for (const auto &pCluster : clusters)
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfo, pCluster));

        // Fragment the original clusters by splitting them making new upstream and downstream clusters
        std::string clusterListName;
        switch (view)
        {
            case TPC_VIEW_U:
                clusterListName = m_clusterListNameU;
                break;
            case TPC_VIEW_V:
                clusterListName = m_clusterListNameV;
                break;
            case TPC_VIEW_W:
                clusterListName = m_clusterListNameW;
                break;
            case TPC_3D:
                clusterListName = LArPfoHelper::IsTrack(pPfo) ? m_trackClusterListName : m_showerClusterListName;
                break;
            default:
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }

        std::cout << "Splitting cluster in list: " << clusterListName << std::endl;

        const StatusCode listChangeStatusCode(PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));
        if (listChangeStatusCode == STATUS_CODE_NOT_FOUND)
            continue;

        if (listChangeStatusCode != STATUS_CODE_SUCCESS)
            throw StatusCodeException(listChangeStatusCode);

        std::string clusterListToSaveName, clusterListToDeleteName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusters, clusterListToDeleteName, clusterListToSaveName));

        PandoraContentApi::Cluster::Parameters clusterParametersUpstream, clusterParametersDownstream;
        clusterParametersUpstream.m_caloHitList = viewToUpstreamHitsMap.at(view);
        clusterParametersDownstream.m_caloHitList = viewToDownstreamHitsMap.at(view);

        if (!clusterParametersUpstream.m_caloHitList.empty())
        {
            const Cluster *pUpstreamCluster = nullptr;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, clusterParametersUpstream, pUpstreamCluster));
            upstreamClusters.push_back(pUpstreamCluster);
        }

        if (!clusterParametersDownstream.m_caloHitList.empty())
        {
            const Cluster *pDownstreamCluster = nullptr;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, clusterParametersDownstream, pDownstreamCluster));
            downstreamClusters.push_back(pDownstreamCluster);
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::SplitPfo(const ParticleFlowObject *const pPfo, const ClusterList &upstreamClusters, const ClusterList &downstreamClusters, PfoList &pfosToSave, PfoList &pfosToDelete) const
{
    std::cout << "Input PFO" << std::endl;
    std::cout << " - nParents   : " << pPfo->GetParentPfoList().size() << std::endl;
    std::cout << " - nDaughters : " << pPfo->GetDaughterPfoList().size() << std::endl;

    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = pPfo->GetParticleId();
    pfoParameters.m_charge = pPfo->GetCharge();
    pfoParameters.m_mass = pPfo->GetMass();
    pfoParameters.m_energy = pPfo->GetEnergy();
    pfoParameters.m_momentum = pPfo->GetMomentum();
    pfoParameters.m_propertiesToAdd = pPfo->GetPropertiesMap();
    pfoParameters.m_clusterList.clear();
    pfoParameters.m_trackList.clear();
    pfoParameters.m_vertexList.clear();

    // Make the upstream PFO
    const ParticleFlowObject *pUpstreamPfo = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pUpstreamPfo));

    for (const auto &pCluster : upstreamClusters)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pUpstreamPfo, pCluster));
    
    pfosToSave.push_back(pUpstreamPfo);
    
    // Make the downstream PFO
    const ParticleFlowObject *pDownstreamPfo = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pDownstreamPfo));
    
    for (const auto &pCluster : downstreamClusters)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pDownstreamPfo, pCluster));

    pfosToSave.push_back(pDownstreamPfo);

    // Mark the original PFO for deletion
    pfosToDelete.push_back(pPfo);
            
    // Make the parent daughter link
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pUpstreamPfo, pDownstreamPfo));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SecondaryInteractionsAlgorithm::GetDimension(const CaloHitList &hits) const
{
    if (hits.size() <= 1)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float minX = std::numeric_limits<float>::max();
    float minZ = std::numeric_limits<float>::max();
    float maxX = -std::numeric_limits<float>::max();
    float maxZ = -std::numeric_limits<float>::max();

    for (const auto &pHit : hits)
    {
        const auto view = pHit->GetHitType();
        if (view != TPC_VIEW_U && view != TPC_VIEW_V && view != TPC_VIEW_W)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const auto pos = pHit->GetPositionVector();
        const auto x = pos.GetX();
        const auto z = pos.GetZ();

        minX = std::min(x, minX);
        minZ = std::min(z, minZ);
        maxX = std::max(x, maxX);
        maxZ = std::max(z, maxZ);
    }

    // Ensure to floating point accuracy all hits are within the bounds
    minX -= std::numeric_limits<float>::epsilon();
    maxX += std::numeric_limits<float>::epsilon();
    minZ -= std::numeric_limits<float>::epsilon();
    maxZ += std::numeric_limits<float>::epsilon();
   
    const auto deltaX = maxX - minX;
    const auto deltaZ = maxZ - minZ;
    const auto maxBoxSize = std::max(deltaX, deltaZ);

    const auto boxSize = 0.3 * 2;
    const auto nBoxes = this->CountBoxes(hits, boxSize, minX, maxX, minZ, maxZ);

    const auto boxSizeFactor = std::log(maxBoxSize / boxSize);
    const auto nBoxesFactor = std::log(nBoxes);
    const auto D = nBoxesFactor / boxSizeFactor;

    LArFormattingHelper::Table table({"minX", "maxX", "minZ", "maxZ", "", "maxBoxSize", "boxSize", "", "nBoxes", "D"});
    table.AddElement(minX);
    table.AddElement(maxX);
    table.AddElement(minZ);
    table.AddElement(maxZ);
    table.AddElement(maxBoxSize);
    table.AddElement(boxSize);
    table.AddElement(nBoxes);
    table.AddElement(D);
    table.Print();

    return 0.f;    
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int SecondaryInteractionsAlgorithm::CountBoxes(const CaloHitList &hits, const float boxSize, const float minX, const float maxX, const float minZ, const float maxZ) const
{
    const auto nBoxesX = static_cast<unsigned int>(std::ceil((maxX - minX) / boxSize));
    const auto nBoxesZ = static_cast<unsigned int>(std::ceil((maxZ - minZ) / boxSize));

    // Make a "matrix" of size nBoxesX:nBoxesZ with every element initialized to false
    std::vector< std::vector<bool> > counter(nBoxesX, std::vector<bool>(nBoxesZ, false));
    for (const auto &pHit : hits)
    {
        const auto pos = pHit->GetPositionVector();
        const auto x = pos.GetX();
        const auto z = pos.GetZ();

        if (x < minX || x > maxX || z < minZ || z > maxZ)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const auto iX = static_cast<unsigned int>(std::floor((x - minX) / boxSize));
        const auto iZ = static_cast<unsigned int>(std::floor((z - minZ) / boxSize));

        // Set the counter to true in this bin, there is a hit there
        counter.at(iX).at(iZ) = true;
    }

    unsigned int nBoxes = 0;
    for (unsigned int iX = 0; iX < nBoxesX; ++iX)
    {
        for (unsigned int iZ = 0; iZ < nBoxesZ; ++iZ)
        {
            if (counter.at(iX).at(iZ))
                nBoxes++;
        }
    }

    return nBoxes;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

SecondaryInteractionsAlgorithm::SplitPoint3D::SplitPoint3D(const Pandora &pandora, const CaloHit *const pHitA, const CaloHit *const pHitB, const CaloHit *const pHitC) :
    m_pHitA(pHitA),
    m_pHitB(pHitB),
    m_pHitC(pHitC),
    m_hasHitC(true),
    m_position3D(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()),
    m_chi2(-std::numeric_limits<float>::max()),
    m_maxDeltaX(-std::numeric_limits<float>::max())
{
    this->CheckValidHitTypes();
    this->SetMaxDeltaX();
    
    LArGeometryHelper::MergeThreePositions3D(pandora, m_pHitA->GetHitType(), m_pHitB->GetHitType(), m_pHitC->GetHitType(), m_pHitA->GetPositionVector(), m_pHitB->GetPositionVector(), m_pHitC->GetPositionVector(), m_position3D, m_chi2);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SecondaryInteractionsAlgorithm::SplitPoint3D::SplitPoint3D(const Pandora &pandora, const CaloHit *const pHitA, const CaloHit *const pHitB) :
    m_pHitA(pHitA),
    m_pHitB(pHitB),
    m_pHitC(nullptr),
    m_hasHitC(false),
    m_position3D(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()),
    m_chi2(-std::numeric_limits<float>::max()),
    m_maxDeltaX(-std::numeric_limits<float>::max())
{
    this->CheckValidHitTypes();
    this->SetMaxDeltaX();
    
    LArGeometryHelper::MergeTwoPositions3D(pandora, m_pHitA->GetHitType(), m_pHitB->GetHitType(), m_pHitA->GetPositionVector(), m_pHitB->GetPositionVector(), m_position3D, m_chi2);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::SplitPoint3D::CheckValidHitTypes() const
{
    const auto hitTypeA = m_pHitA->GetHitType();
    const auto hitTypeB = m_pHitB->GetHitType();
    const auto hitTypeC = m_hasHitC ? m_pHitC->GetHitType() : HIT_CUSTOM;

    if (hitTypeA == hitTypeB || hitTypeA == hitTypeC || hitTypeB == hitTypeC)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (hitTypeA != TPC_VIEW_U && hitTypeA != TPC_VIEW_V && hitTypeA != TPC_VIEW_W)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    
    if (hitTypeB != TPC_VIEW_U && hitTypeB != TPC_VIEW_V && hitTypeB != TPC_VIEW_W)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    
    if (m_hasHitC && hitTypeC != TPC_VIEW_U && hitTypeC != TPC_VIEW_V && hitTypeC != TPC_VIEW_W)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::SplitPoint3D::SetMaxDeltaX()
{
    const auto hitPosA = m_pHitA->GetPositionVector();
    const auto hitPosB = m_pHitB->GetPositionVector();
    
    const auto deltaXAB = std::abs((hitPosA - hitPosB).GetX());
    m_maxDeltaX = deltaXAB;

    if (m_hasHitC)
    {
        const auto hitPosC = m_pHitC->GetPositionVector();
    
        const auto deltaXAC = std::abs((hitPosA - hitPosC).GetX());
        const auto deltaXBC = std::abs((hitPosB - hitPosC).GetX());

        m_maxDeltaX = std::max(std::max(deltaXAC, deltaXBC), deltaXAB);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *SecondaryInteractionsAlgorithm::SplitPoint3D::GetHitWithView(const HitType &view) const
{
    if (m_pHitA->GetHitType() == view)
        return m_pHitA;
    
    if (m_pHitB->GetHitType() == view)
        return m_pHitB;
    
    if (m_hasHitC)
    {
        if (m_pHitC->GetHitType() == view)
            return m_pHitC;
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

StatusCode SecondaryInteractionsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameU", m_clusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameV", m_clusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameW", m_clusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackClusterListName", m_trackClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerClusterListName", m_showerClusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "IsolatedHitDistance", m_isolatedHitDistance));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "StitchingThresholdDistance", m_stitchingThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TransverseBias", m_transverseBias));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NumberOfSampleHits", m_nSampleHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CosAngleThreshold", m_cosAngleThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MaxMatchDeltaX", m_maxMatchDeltaX));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MaxMatchThreeViewChi2", m_maxMatch3ViewChi2));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TwoViewProjectionThreshold", m_twoViewProjectionThreshold));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
