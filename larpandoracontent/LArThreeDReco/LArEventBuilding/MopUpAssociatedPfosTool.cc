/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/MopUpAssociatedPfosTool.cc
 *
 *  @brief  Implementation of the end associated pfos tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/MopUpAssociatedPfosTool.h"

using namespace pandora;

namespace lar_content
{

typedef NeutrinoHierarchyAlgorithm::PfoInfo PfoInfo;
typedef NeutrinoHierarchyAlgorithm::PfoInfoMap PfoInfoMap;

MopUpAssociatedPfosTool::MopUpAssociatedPfosTool() :
    m_shouldDisplayMatchingInfo(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpAssociatedPfosTool::Run(NeutrinoHierarchyAlgorithm *const pAlgorithm, const Vertex *const pNeutrinoVertex, PfoInfoMap &pfoInfoMap)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    PfoVector assignedPfos, unassignedPfos;
    pAlgorithm->SeparatePfos(pfoInfoMap, assignedPfos, unassignedPfos);

    if (unassignedPfos.empty())
        return;

    PfoMatchingMap pfoMatchingMap;
    this->PopulatePfoMatchingMap(pfoInfoMap, unassignedPfos, assignedPfos, pNeutrinoVertex, pfoMatchingMap);

    while (!unassignedPfos.empty())
    {
        if (m_shouldDisplayMatchingInfo)
            this->PrintPfoMatchingMap(unassignedPfos, assignedPfos, pfoMatchingMap, pfoInfoMap);

        // Only assign parent / daughter links to already assigned PFOs to avoid circular hierarchies
        const bool wasLinkMade(this->MakeClearLinksToNeutrinoHierarchy(unassignedPfos, assignedPfos, pfoMatchingMap, pfoInfoMap));

        // If no such links were the best choice, then make the next best available link to avoid isolated parts of the hierarchy
        if (!wasLinkMade)
        {
            this->MakeNextBestLinkToNeutrinoHierarchy(unassignedPfos, assignedPfos, pfoMatchingMap, pfoInfoMap);
        }

        unassignedPfos.clear();
        assignedPfos.clear();
        pAlgorithm->SeparatePfos(pfoInfoMap, assignedPfos, unassignedPfos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpAssociatedPfosTool::PopulatePfoMatchingMap(const PfoInfoMap &pfoInfoMap, const PfoVector &unassignedPfos,
    const PfoVector &assignedPfos, const Vertex *const pNeutrinoVertex, PfoMatchingMap &pfoMatchingMap) const
{
    PfoVector allPfos;
    allPfos.insert(allPfos.end(), unassignedPfos.begin(), unassignedPfos.end());
    allPfos.insert(allPfos.end(), assignedPfos.begin(), assignedPfos.end());

    for (const ParticleFlowObject *const pPfo : unassignedPfos)
        pfoMatchingMap.emplace(pPfo, std::make_shared<PfoMatchDetails>(pPfo, pfoInfoMap, pNeutrinoVertex, allPfos));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpAssociatedPfosTool::PrintPfoMatchingMap(const PfoVector &unassignedPfos, const PfoVector &assignedPfos,
    const PfoMatchingMap &pfoMatchingMap, const PfoInfoMap &pfoInfoMap) const
{
    LArFormattingHelper::PrintHeader("Unassigned PFOs");
    LArFormattingHelper::Table pfoTable({"PFO", "", "Vertex score", "", "Best PFO", "Score", "", "Best assigned PFO", "Score"});

    for (const ParticleFlowObject *const pPfo : unassignedPfos)
    {
        const auto &pMatchDetails(pfoMatchingMap.at(pPfo));
        
        pfoTable.AddElement(pPfo);

        bool dummyBool(false);
        float score(-std::numeric_limits<float>::max());

        pMatchDetails->GetVertexMatchScore(score, dummyBool);
        pfoTable.AddElement(score);
        
        const ParticleFlowObject *pBestMatchedPfo(nullptr);
        pMatchDetails->GetBestMatchedPfo(pfoInfoMap, assignedPfos, pBestMatchedPfo, score, dummyBool);
        pfoTable.AddElement(pBestMatchedPfo);
        pfoTable.AddElement(score);
        
        pMatchDetails->GetBestMatchedAssignedPfo(pfoInfoMap, assignedPfos, pBestMatchedPfo, score, dummyBool);
        pfoTable.AddElement(pBestMatchedPfo);
        pfoTable.AddElement(score);
    }

    pfoTable.Print();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MopUpAssociatedPfosTool::MakeClearLinksToNeutrinoHierarchy(const PfoVector &unassignedPfos, const PfoVector &assignedPfos,
    const PfoMatchingMap &pfoMatchingMap, PfoInfoMap &pfoInfoMap) const
{
    bool wasMatchMade(false);

    for (const ParticleFlowObject *const pPfo : unassignedPfos)
    {
        const auto &pMatchDetails(pfoMatchingMap.at(pPfo));
        
        const ParticleFlowObject *pBestMatchedPfo(nullptr);
        float pfoMatchScore(-std::numeric_limits<float>::max());
        bool pfoMatchUseInner(false);
        pMatchDetails->GetBestMatchedPfo(pfoInfoMap, assignedPfos, pBestMatchedPfo, pfoMatchScore, pfoMatchUseInner);

        bool vertexMatchUseInner(false);
        float vertexMatchScore(-std::numeric_limits<float>::max());
        pMatchDetails->GetVertexMatchScore(vertexMatchScore, vertexMatchUseInner);

        if (!pBestMatchedPfo || vertexMatchScore > pfoMatchScore)
        {
            this->MakeVertexLink(pPfo, vertexMatchUseInner, pfoInfoMap);
            wasMatchMade = true;
        } 
        else if (std::find(assignedPfos.begin(), assignedPfos.end(), pBestMatchedPfo) != assignedPfos.end())
        {
            this->MakePfoLink(pPfo, pBestMatchedPfo, pfoMatchUseInner, pfoInfoMap);
            wasMatchMade = true;
        }
    }

    return wasMatchMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpAssociatedPfosTool::MakeNextBestLinkToNeutrinoHierarchy(const PfoVector &unassignedPfos, const PfoVector &assignedPfos,
    const PfoMatchingMap &pfoMatchingMap, PfoInfoMap &pfoInfoMap) const
{
    float bestMatchScore(-std::numeric_limits<float>::max());
    const ParticleFlowObject *pBestChildPfo(nullptr);
    const ParticleFlowObject *pBestParentPfo(nullptr);
    bool useInner(false);

    for (const ParticleFlowObject *const pPfo : unassignedPfos)
    {
        const auto &pMatchDetails(pfoMatchingMap.at(pPfo));
        
        const ParticleFlowObject *pBestMatchedPfo(nullptr);
        float pfoMatchScore(-std::numeric_limits<float>::max());
        bool pfoMatchUseInner(false);
        pMatchDetails->GetBestMatchedAssignedPfo(pfoInfoMap, assignedPfos, pBestMatchedPfo, pfoMatchScore, pfoMatchUseInner);

        bool vertexMatchUseInner(false);
        float vertexMatchScore(-std::numeric_limits<float>::max());
        pMatchDetails->GetVertexMatchScore(vertexMatchScore, vertexMatchUseInner);

        if (vertexMatchScore < bestMatchScore && pfoMatchScore < bestMatchScore)
            continue;

        if (!pBestMatchedPfo || vertexMatchScore >= pfoMatchScore)
        {
            bestMatchScore = vertexMatchScore;
            pBestChildPfo = pPfo;
            pBestParentPfo = nullptr;
            useInner = vertexMatchUseInner;
        }
        else if (pBestMatchedPfo && pfoMatchScore > vertexMatchScore)
        {
            bestMatchScore = pfoMatchScore;
            pBestChildPfo = pPfo;
            pBestParentPfo = pBestMatchedPfo;
            useInner = pfoMatchUseInner;
        }
    }

    if (!pBestChildPfo)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    if (!pBestParentPfo)
    {
        this->MakeVertexLink(pBestChildPfo, useInner, pfoInfoMap);
    }
    else
    {
        this->MakePfoLink(pBestChildPfo, pBestParentPfo, useInner, pfoInfoMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpAssociatedPfosTool::MakeVertexLink(const ParticleFlowObject *const pPfo, const bool useInner, PfoInfoMap &pfoInfoMap) const
{
    PfoInfo *const pPfoInfo(pfoInfoMap.at(pPfo));

    pPfoInfo->SetNeutrinoVertexAssociation(true);
    pPfoInfo->SetInnerLayerAssociation(useInner);

    if (m_shouldDisplayMatchingInfo)
        std::cout << "Matching " << pPfo << " -> Neutrino Vertex" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpAssociatedPfosTool::MakePfoLink(const ParticleFlowObject *const pChildPfo, const ParticleFlowObject *const pParentPfo,
    const bool useInner, PfoInfoMap &pfoInfoMap) const
{
    PfoInfo *const pChildPfoInfo(pfoInfoMap.at(pChildPfo));
    PfoInfo *const pParentPfoInfo(pfoInfoMap.at(pParentPfo));

    pParentPfoInfo->AddDaughterPfo(pChildPfo);
    pChildPfoInfo->SetParentPfo(pParentPfo);
    pChildPfoInfo->SetInnerLayerAssociation(useInner);
    
    if (m_shouldDisplayMatchingInfo)
        std::cout << "Matching " << pChildPfo << " -> " << pParentPfo << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MopUpAssociatedPfosTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DisplayMatchingInfo", m_shouldDisplayMatchingInfo));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MopUpAssociatedPfosTool::PfoMatchDetails::PfoMatchDetails(const ParticleFlowObject *const pPfo,
    const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap, const Vertex *const pNeutrinoVertex, const PfoVector &allPfos) :
    m_pPfo(pPfo)
{
    this->SetVertexMatch(pfoInfoMap, pNeutrinoVertex);

    for (const ParticleFlowObject *const pParentPfo : allPfos)
    {
        if (pParentPfo == pPfo)
            continue;

        this->AddMatch(pfoInfoMap, pParentPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MopUpAssociatedPfosTool::PfoMatchDetails::GetMatchScore(const CartesianVector &daughterPosition,
    const CartesianVector &parentPosition) const
{
    const float magSquared = (daughterPosition - parentPosition).GetMagnitudeSquared();

    if (magSquared <= std::numeric_limits<float>::epsilon())
        return std::numeric_limits<float>::max();

    return 1.f / magSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MopUpAssociatedPfosTool::PfoMatchDetails::SetVertexMatch(const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap,
    const Vertex *const pNeutrinoVertex)
{
    PfoInfo *const pPfoInfo(pfoInfoMap.at(m_pPfo));
    const LArPointingCluster pointingCluster(*(pPfoInfo->GetSlidingFitResult3D()));
    const CartesianVector innerPosition(pointingCluster.GetInnerVertex().GetPosition());
    const CartesianVector outerPosition(pointingCluster.GetOuterVertex().GetPosition());
    const CartesianVector neutrinoPosition(pNeutrinoVertex->GetPosition());

    m_vertexMatchInfo.m_innerScore = this->GetMatchScore(innerPosition, neutrinoPosition);
    m_vertexMatchInfo.m_outerScore = this->GetMatchScore(outerPosition, neutrinoPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------
        
void MopUpAssociatedPfosTool::PfoMatchDetails::AddMatch(const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap,
    const ParticleFlowObject *const pParentPfo)
{
    PfoInfo *const pChildPfoInfo(pfoInfoMap.at(m_pPfo));
    const LArPointingCluster childPointingCluster(*(pChildPfoInfo->GetSlidingFitResult3D()));
    const CartesianVector childInnerPosition(childPointingCluster.GetInnerVertex().GetPosition());
    const CartesianVector childOuterPosition(childPointingCluster.GetOuterVertex().GetPosition());
    
    PfoInfo *const pParentPfoInfo(pfoInfoMap.at(pParentPfo));
    const LArPointingCluster parentPointingCluster(*(pParentPfoInfo->GetSlidingFitResult3D()));
    const CartesianVector parentInnerPosition(parentPointingCluster.GetInnerVertex().GetPosition());
    const CartesianVector parentOuterPosition(parentPointingCluster.GetOuterVertex().GetPosition());

    PfoMatchInfo innerInnerMatchInfo;
    innerInnerMatchInfo.m_pParentPfo = pParentPfo; 
    innerInnerMatchInfo.m_isChildInner = true;
    innerInnerMatchInfo.m_isParentInner = true;
    innerInnerMatchInfo.m_score = this->GetMatchScore(childInnerPosition, parentInnerPosition);
    m_pfoMatchInfoVector.push_back(innerInnerMatchInfo);

    PfoMatchInfo innerOuterMatchInfo;
    innerOuterMatchInfo.m_pParentPfo = pParentPfo; 
    innerOuterMatchInfo.m_isChildInner = true;
    innerOuterMatchInfo.m_isParentInner = false;
    innerOuterMatchInfo.m_score = this->GetMatchScore(childInnerPosition, parentOuterPosition);
    m_pfoMatchInfoVector.push_back(innerOuterMatchInfo);
    
    PfoMatchInfo outerInnerMatchInfo;
    outerInnerMatchInfo.m_pParentPfo = pParentPfo; 
    outerInnerMatchInfo.m_isChildInner = false;
    outerInnerMatchInfo.m_isParentInner = true;
    outerInnerMatchInfo.m_score = this->GetMatchScore(childOuterPosition, parentInnerPosition);
    m_pfoMatchInfoVector.push_back(outerInnerMatchInfo);
    
    PfoMatchInfo outerOuterMatchInfo;
    outerOuterMatchInfo.m_pParentPfo = pParentPfo; 
    outerOuterMatchInfo.m_isChildInner = false;
    outerOuterMatchInfo.m_isParentInner = false;
    outerOuterMatchInfo.m_score = this->GetMatchScore(childOuterPosition, parentOuterPosition);
    m_pfoMatchInfoVector.push_back(outerOuterMatchInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------
        
void MopUpAssociatedPfosTool::PfoMatchDetails::GetBestMatchedAssignedPfo(const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap,
    const pandora::PfoVector &assignedPfos, const pandora::ParticleFlowObject *&pBestMatchedPfo, float &matchScore, bool &useInner) const
{
    pBestMatchedPfo = nullptr;
    matchScore = -std::numeric_limits<float>::max();

    for (const PfoMatchInfo &matchInfo : m_pfoMatchInfoVector)
    {
        if (matchInfo.m_score <= matchScore)
            continue;

        const auto iter(std::find(assignedPfos.begin(), assignedPfos.end(), matchInfo.m_pParentPfo));
        const bool isMatchAlreadyAssigned(iter != assignedPfos.end());

        if (!isMatchAlreadyAssigned)
            continue;

        // ATTN we should only associate to the end points (not the start points) of associated PFOs
        PfoInfo *const pParentPfoInfo(pfoInfoMap.at(matchInfo.m_pParentPfo));
        if (pParentPfoInfo->IsInnerLayerAssociated() == matchInfo.m_isParentInner)
            continue;

        pBestMatchedPfo = matchInfo.m_pParentPfo;
        matchScore = matchInfo.m_score;
        useInner = matchInfo.m_isChildInner; 
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
        
void MopUpAssociatedPfosTool::PfoMatchDetails::GetBestMatchedPfo(const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap,
    const pandora::PfoVector &assignedPfos, const pandora::ParticleFlowObject *&pBestMatchedPfo, float &matchScore, bool &useInner) const
{
    pBestMatchedPfo = nullptr;
    matchScore = -std::numeric_limits<float>::max();

    for (const PfoMatchInfo &matchInfo : m_pfoMatchInfoVector)
    {
        if (matchInfo.m_score <= matchScore)
            continue;

        const auto iter(std::find(assignedPfos.begin(), assignedPfos.end(), matchInfo.m_pParentPfo));
        const bool isMatchAlreadyAssigned(iter != assignedPfos.end());

        if (isMatchAlreadyAssigned)
        {
            // ATTN we should only associate to the end points (not the start points) of associated PFOs, for unassociated PFOs we can
            // link to either the inner or outer vertices
            PfoInfo *const pParentPfoInfo(pfoInfoMap.at(matchInfo.m_pParentPfo));
            if (pParentPfoInfo->IsInnerLayerAssociated() == matchInfo.m_isParentInner)
                continue;
        }

        pBestMatchedPfo = matchInfo.m_pParentPfo;
        matchScore = matchInfo.m_score;
        useInner = matchInfo.m_isChildInner; 
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
        
void MopUpAssociatedPfosTool::PfoMatchDetails::GetVertexMatchScore(float &matchScore, bool &useInner) const
{
    useInner = (m_vertexMatchInfo.m_innerScore > m_vertexMatchInfo.m_outerScore);
    matchScore = useInner ? m_vertexMatchInfo.m_innerScore : m_vertexMatchInfo.m_outerScore;
}

} // namespace lar_content
