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

using namespace pandora;

namespace lar_content
{

SecondaryInteractionsAlgorithm::SecondaryInteractionsAlgorithm() :
    m_isolatedHitDistance(2.f),
    m_transverseBias(2.f),
    m_nSampleHits(3)
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
        CaloHitList caloHitsU, caloHitsV, caloHitsW;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, caloHitsU);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, caloHitsV);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitsW);

        CaloHitList processedHitsU, processedHitsV, processedHitsW;

        this->FilterAndOrderHitsByNearestNeighbor(pVertex, caloHitsU, processedHitsU);
        this->FilterAndOrderHitsByNearestNeighbor(pVertex, caloHitsV, processedHitsV);
        this->FilterAndOrderHitsByNearestNeighbor(pVertex, caloHitsW, processedHitsW);
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

void SecondaryInteractionsAlgorithm::FilterAndOrderHitsByNearestNeighbor(const Vertex *const pVertex, const CaloHitList &caloHitList, pandora::CaloHitList &outputHitList) const
{
    // We don't need to order 0 or 1 hit
    if (caloHitList.size() < 2)
        return;

    HitSeparationMap separationMap;
    this->GetHitSeparationMap(caloHitList, separationMap);

    std::vector<CaloHitList> continuousSegments;
    this->GetContinuousSegments(pVertex, caloHitList, separationMap, continuousSegments);

    std::vector< std::vector<unsigned int> > segmentKinkIndices;
    this->GetKinkIndices(continuousSegments, segmentKinkIndices);

    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "AllHits", AUTO));
    std::cout << "Found " << continuousSegments.size() << " segments" << std::endl;

    std::vector<Color> colors = {RED, GREEN, BLUE, MAGENTA, CYAN, VIOLET, PINK};
    for (unsigned int i = 0; i < continuousSegments.size(); ++i)
    {
        std::cout << "  - Segment " << i << " : " << continuousSegments.at(i).size() << std::endl;
        const auto name = "Segment " + std::to_string(i);
        const auto color = colors.at(i % colors.size());
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &continuousSegments.at(i), name, color));
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &continuousSegments.at(i).front()->GetPositionVector(), name + "Start", color, 4));
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &continuousSegments.at(i).back()->GetPositionVector(), name + "End", color, 2));
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    (void) outputHitList;
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
    CaloHitList remainingHits;
    this->GetNonIsolatedHits(caloHitList, separationMap, remainingHits);

    while (!remainingHits.empty())
    {
        const auto pSeedHit = this->GetClosestHitToVertex(remainingHits, pVertex);

        CaloHitList segment;
        this->CollectHitsInSegment(pSeedHit, separationMap, remainingHits, segment);
        continuousSegments.push_back(segment);
    }
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

void SecondaryInteractionsAlgorithm::GetNonIsolatedHits(const CaloHitList &caloHitList, const HitSeparationMap &separationMap, CaloHitList &nonIsolatedHits) const
{
    for (const auto &pHit : caloHitList)
    {
        // The separation map only contains hits which are within the isolated hit distance threshold
        if (separationMap.find(pHit) != separationMap.end())
            nonIsolatedHits.push_back(pHit);
    }
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
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

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

void SecondaryInteractionsAlgorithm::GetKinkIndices(const std::vector<CaloHitList> &continuousSegments, std::vector< std::vector<unsigned int> > &segmentKinkIndices) const
{
    for (const auto &segment : continuousSegments)
    {
        std::vector<unsigned int> kinkIndices;
        this->GetKinkIndices(segment, kinkIndices);
        segmentKinkIndices.push_back(kinkIndices);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SecondaryInteractionsAlgorithm::GetKinkIndices(const CaloHitList &segment, std::vector<unsigned int> &/*kinkIndices*/) const
{
    if (segment.size() < 1 + 2 * m_nSampleHits)
        return;

    for (unsigned int i = m_nSampleHits; i < segment.size() - m_nSampleHits; ++i)
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &segment, "FullSegment", AUTO));

        // Get the hits before this hit
        auto preBegin = segment.begin();
        std::advance(preBegin, i - m_nSampleHits);
        auto preEnd = preBegin;
        std::advance(preEnd, m_nSampleHits);
        const CaloHitList preHits(preBegin, preEnd);
        
        // Get the hits after this hit
        auto postBegin = segment.begin();
        std::advance(postBegin, i + 1);
        auto postEnd = postBegin;
        std::advance(postEnd, m_nSampleHits);
        const CaloHitList postHits(postBegin, postEnd);
        
        auto thisBegin = segment.begin();
        std::advance(thisBegin, i);
        auto thisEnd = thisBegin;
        std::advance(thisEnd, 1);
        const CaloHitList thisHit(thisBegin, thisBegin);

        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &preHits, "Pre", BLUE));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &thisHit, "Current", RED));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &postHits, "Post", GREEN));

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

StatusCode SecondaryInteractionsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "IsolatedHitDistance", m_isolatedHitDistance));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TransverseBias", m_transverseBias));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_NOT_FOUND, STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NumberOfSampleHits", m_nSampleHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
