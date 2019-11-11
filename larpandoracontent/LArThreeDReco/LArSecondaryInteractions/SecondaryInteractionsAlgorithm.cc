/**
 *  @file   larpandoracontent/LArThreeDReco/LArSecondaryInteractions/SecondaryInteractionsAlgorithm.cc
 *
 *  @brief  Implementation of the secondary interactions algorithm class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArSecondaryInteractions/SecondaryInteractionsAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

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
    /* BEGIN DEBUG */
    std::vector<Color> colors = {RED, GREEN, BLUE, MAGENTA, CYAN, VIOLET, PINK};
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    /* END DEBUG */

    // Get the input collections
    PfoList inputPfos;
    this->CollectInputPfos(inputPfos);
    const auto pVertex = this->GetVertex();

    for (const auto &pPfo : inputPfos)
    {
        CaloHitList splitHits;
        CaloHitList allHits;
        for (const auto &view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            // Collect the calo hits
            CaloHitList caloHits;
            LArPfoHelper::GetCaloHits(pPfo, view, caloHits);
            allHits.insert(allHits.end(), caloHits.begin(), caloHits.end());

            HitSeparationMap separationMap;
            this->GetHitSeparationMap(caloHits, separationMap);
            
            // Organise the hits into continuous segments
            std::vector<CaloHitList> segments;
            this->GetContinuousSegments(pVertex, caloHits, separationMap, segments);
            
            // Get the hits that sit at a kink or bifurcation position
            this->GetSplitHits(segments, separationMap, splitHits);
        
            /* BEGIN DEBUG */
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, "AllHits" + std::to_string(view), AUTO));

            for (unsigned int i = 0; i < segments.size(); ++i)
            {
                const auto hits = segments.at(i);

                const auto name = "Segment" + std::to_string(i);
                const auto color = colors.at(i % colors.size());
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, name, color));
            }
            /* END DEBUG */
        }
       
        std::vector<SplitPoint3D> splitPoints3D;
        this->Get3DSplitPoints(allHits, splitHits, splitPoints3D);

        /* BEGIN DEBUG */
        std::cout << "2D split hits" << std::endl;
        for (unsigned int i = 0; i < splitHits.size(); ++i)
        {
            const auto kinkPosition = (*std::next(splitHits.begin(), i))->GetPositionVector();

            const auto name = "Split2D" + std::to_string(i);
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &kinkPosition, name, BLACK, 1));
        
            std::cout << "  - " << name << std::endl;
            std::cout << "  - View : " << (*std::next(splitHits.begin(), i))->GetHitType() << std::endl;
            std::cout << "  - X : " << kinkPosition.GetX() << std::endl;
            std::cout << "  - Z : " << kinkPosition.GetZ() << std::endl;
        }
        
        for (unsigned int i = 0; i < splitPoints3D.size(); ++i)
        {
            const auto splitPosition = splitPoints3D.at(i).m_position3D;
            const auto name = "Split3D" + std::to_string(i);
            const auto color = colors.at(i % colors.size());
            
            const auto splitPositionU = LArGeometryHelper::ProjectPosition(this->GetPandora(), splitPosition, TPC_VIEW_U);
            const auto splitPositionV = LArGeometryHelper::ProjectPosition(this->GetPandora(), splitPosition, TPC_VIEW_V);
            const auto splitPositionW = LArGeometryHelper::ProjectPosition(this->GetPandora(), splitPosition, TPC_VIEW_W);

            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &splitPosition, name, color, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &splitPositionU, name + "U", color, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &splitPositionV, name + "V", color, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &splitPositionW, name + "W", color, 2));
        }

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        /* END DEBUG */
    }

    return STATUS_CODE_SUCCESS;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

const Vertex *SecondaryInteractionsAlgorithm::GetVertex() const
{
    const VertexList *pVertexList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));

    if (pVertexList->size() != 1)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    return pVertexList->front();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::CollectInputPfos(PfoList &inputPfos) const
{
    for (const auto &listName : m_pfoListNames)
    {
        const PfoList *pPfoList = nullptr;

        if (PandoraContentApi::GetList(*this, listName, pPfoList) == STATUS_CODE_SUCCESS)
            inputPfos.insert(inputPfos.end(), pPfoList->begin(), pPfoList->end());
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::GetHitSeparationMap(const CaloHitList &caloHitList, HitSeparationMap &separationMap) const
{
    const auto isolatedHitDistanceSquared = m_isolatedHitDistance * m_isolatedHitDistance;

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

void SecondaryInteractionsAlgorithm::GetContinuousSegments(const Vertex *const pVertex, const CaloHitList &caloHitList, const HitSeparationMap &separationMap, std::vector<CaloHitList> &continuousSegments) const
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
        const auto pSeedHit = this->GetClosestHitToVertex(remainingHits, pVertex);

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

const CaloHit *SecondaryInteractionsAlgorithm::GetClosestHitToVertex(const CaloHitList &caloHitList, const Vertex *const pVertex) const
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
    const auto vertexPosition3D = pVertex->GetPosition();
    const auto vertexPosition = (is2D ? LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition3D, hitType) : vertexPosition3D);

    return this->GetClosestHitToPoint(caloHitList, vertexPosition);
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

void SecondaryInteractionsAlgorithm::GetSplitHits(const std::vector<CaloHitList> &continuousSegments, const HitSeparationMap &separationMap, CaloHitList &splitHits) const
{
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

void SecondaryInteractionsAlgorithm::Get3DSplitPoints(const CaloHitList &allHits, const CaloHitList &splitHits, std::vector<SplitPoint3D> &splitPoints3D) const
{
    CaloHitList allHitsU, allHitsV, allHitsW;
    this->OrganiseHitsByView(allHits, allHitsU, allHitsV, allHitsW);

    CaloHitList splitHitsU, splitHitsV, splitHitsW;
    this->OrganiseHitsByView(splitHits, splitHitsU, splitHitsV, splitHitsW);

    std::cout << "Getting 3D split points" << std::endl;
    std::cout << " - nSplitHitsU : " << splitHitsU.size() << std::endl;
    std::cout << " - nSplitHitsV : " << splitHitsV.size() << std::endl;
    std::cout << " - nSplitHitsW : " << splitHitsW.size() << std::endl;

    while (this->MakeThreeViewMatch(splitHitsU, splitHitsV, splitHitsW, splitPoints3D)) {}
    std::cout << "After 3-view matching found " << splitPoints3D.size() << " split points total" << std::endl;

    while (this->MakeTwoViewMatch(allHitsU, allHitsV, allHitsW, splitHitsU, splitHitsV, splitHitsW, splitPoints3D)) {}
    std::cout << "After 2-view matching found " << splitPoints3D.size() << " split points total" << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::OrganiseHitsByView(const CaloHitList &allHits, CaloHitList &hitsU, CaloHitList &hitsV, CaloHitList &hitsW) const
{
    for (const auto &pHit : allHits)
    {
        switch (pHit->GetHitType())
        {
            case TPC_VIEW_U:
                hitsU.push_back(pHit);
                break;
            case TPC_VIEW_V:
                hitsV.push_back(pHit);
                break;
            case TPC_VIEW_W:
                hitsW.push_back(pHit);
                break;
            default:
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SecondaryInteractionsAlgorithm::MakeThreeViewMatch(CaloHitList &hitsU, CaloHitList &hitsV, CaloHitList &hitsW, std::vector<SplitPoint3D> &splitPoints3D) const
{
    if (hitsU.empty() || hitsV.empty() || hitsW.empty())
        return false;

    std::cout << "Looking for 3-view matches for split hits" << std::endl;
    std::cout << "  - NHitsU   : " << hitsU.size() << std::endl;
    std::cout << "  - NHitsV   : " << hitsV.size() << std::endl;
    std::cout << "  - NHitsW   : " << hitsW.size() << std::endl;
    std::cout << "  - NMatches : " << splitPoints3D.size() << std::endl;

    // Collect the viable split points
    std::vector<SplitPoint3D> viableSplitPoints;
    for (const auto &pHitU : hitsU)
    {
        for (const auto &pHitV : hitsV)
        {
            for (const auto &pHitW : hitsW)
            {
                const SplitPoint3D splitPoint(this->GetPandora(), pHitU, pHitV, pHitW);

                std::cout << "Considering triplet" << std::endl;
                std::cout << "  - U    : " << pHitU << std::endl;
                std::cout << "  - V    : " << pHitV << std::endl;
                std::cout << "  - W    : " << pHitW << std::endl;
                std::cout << "  - dX   : " << splitPoint.m_maxDeltaX << std::endl;
                std::cout << "  - chi2 : " << splitPoint.m_chi2 << std::endl;

                if (splitPoint.m_maxDeltaX > m_maxMatchDeltaX)
                    continue;
                
                std::cout << "  - Passes dX cut!" << std::endl;
                
                if (splitPoint.m_chi2 > m_maxMatch3ViewChi2)
                    continue;
                
                std::cout << "  - Passes chi2 cut!" << std::endl;

                viableSplitPoints.push_back(splitPoint);
            }
        }
    }
                
    std::cout << "Found " << viableSplitPoints.size() << " viable triplets" << std::endl;

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
    
    std::cout << "Best triplet is" << std::endl;
    std::cout << "  - U : " << bestSplitPoint.GetHitWithView(TPC_VIEW_U) << std::endl;
    std::cout << "  - V : " << bestSplitPoint.GetHitWithView(TPC_VIEW_V) << std::endl;
    std::cout << "  - W : " << bestSplitPoint.GetHitWithView(TPC_VIEW_W) << std::endl;

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
    
    std::cout << "Found " << viableSplitPoints.size() << " viable doublets" << std::endl;

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
    
    std::cout << "Best doublet is" << std::endl;
    std::cout << "  - A : " << bestSplitPoint.m_pHitA << std::endl;
    std::cout << "  - B : " << bestSplitPoint.m_pHitB << std::endl;

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
             
            std::cout << "Considering doublet" << std::endl;
            std::cout << "  - A    : " << pHitA << " - " << hitTypeA << std::endl;
            std::cout << "  - B    : " << pHitB << " - " << hitTypeB << std::endl;
            std::cout << "  - dX   : " << splitPoint.m_maxDeltaX << std::endl;
            std::cout << "  - chi2 : " << splitPoint.m_chi2 << std::endl;

            if (splitPoint.m_maxDeltaX > m_maxMatchDeltaX)
                continue;

            std::cout << "  - Passes dX cut!" << std::endl;

            // Get the 3D position as projected into the remaining view
            const auto projectedPosition = LArGeometryHelper::ProjectPosition(this->GetPandora(), splitPoint.m_position3D, hitTypeC);
            const auto isInGap = LArGeometryHelper::IsInGap(this->GetPandora(), projectedPosition, hitTypeC, m_twoViewProjectionThreshold);

            if (isInGap) std::cout << "  - Projects to gap!" << std::endl;

            bool projectsToHit = false;
            if (!hitsC.empty())
            {
                const auto pNearestHit = this->GetClosestHitToPoint(hitsC, projectedPosition);
                const auto nearestHitDistSquared = projectedPosition.GetDistanceSquared(pNearestHit->GetPositionVector());
                std::cout << "  - Nearest hit distance (squared) : " << nearestHitDistSquared << std::endl;
                projectsToHit = (nearestHitDistSquared < m_twoViewProjectionThreshold * m_twoViewProjectionThreshold);
            }
            else
            {
                std::cout << "  - No hits in remaining view!" << std::endl;
            }

            if (projectsToHit) std::cout << "  - Projects to hit! (not using)" << std::endl;

            if (!isInGap && !projectsToHit)
                continue;
                
            splitPointsAB.push_back(splitPoint);
        }
    }
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
