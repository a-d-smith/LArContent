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
    m_cosAngleThreshold(0.9)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

StatusCode SecondaryInteractionsAlgorithm::Run()
{
    // Get the input collections
    PfoList inputPfos;
    this->CollectInputPfos(inputPfos);
    const auto pVertex = this->GetVertex();

    for (const auto &pPfo : inputPfos)
    {
        for (const auto &view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            // Collect the calo hits
            CaloHitList caloHits;
            LArPfoHelper::GetCaloHits(pPfo, view, caloHits);

            HitSeparationMap separationMap;
            this->GetHitSeparationMap(caloHits, separationMap);
            
            // Organise the hits into continuous segments
            std::vector<CaloHitList> segments;
            this->GetContinuousSegments(pVertex, caloHits, separationMap, segments);
            
            // Get the hits that sit at a kink or bifurcation position
            CaloHitList splitHits;
            this->GetSplitHits(segments, separationMap, splitHits);
            
            /* BEGIN DEBUG */
            std::cout << "Found " << segments.size() << " segments and " << splitHits.size() << " kinks in view " << view << std::endl;

            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, "AllHits" + std::to_string(view), AUTO));

            std::vector<Color> colors = {RED, GREEN, BLUE, MAGENTA, CYAN, VIOLET, PINK};
            for (unsigned int i = 0; i < segments.size(); ++i)
            {
                const auto hits = segments.at(i);

                const auto name = "Segment" + std::to_string(i);
                const auto color = colors.at(i % colors.size());
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, name, color));
            }
            
            for (unsigned int i = 0; i < splitHits.size(); ++i)
            {
                const auto kinkPosition = (*std::next(splitHits.begin(), i))->GetPositionVector();

                const auto name = "Kink" + std::to_string(i);
                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &kinkPosition, name, BLACK, 2));
            }
            /* END DEBUG */
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
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

    // Find the closest hit
    float minDistanceSquared = std::numeric_limits<float>::max();
    const CaloHit *pClosestHit = nullptr;
    for (const auto &pHit : caloHitList)
    {
        const auto distSquared = (pHit->GetPositionVector() - vertexPosition).GetMagnitudeSquared();
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

StatusCode SecondaryInteractionsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "IsolatedHitDistance", m_isolatedHitDistance));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "StitchingThresholdDistance", m_stitchingThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TransverseBias", m_transverseBias));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NumberOfSampleHits", m_nSampleHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CosAngleThreshold", m_cosAngleThreshold));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
