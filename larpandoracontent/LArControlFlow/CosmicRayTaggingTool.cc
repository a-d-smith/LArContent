/**
 *  @file   larpandoracontent/LArControlFlow/CosmicRayTaggingTool.cc
 *
 *  @brief  Implementation of the cosmic-ray tagging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/CosmicRayTaggingTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;

namespace lar_content
{

CosmicRayTaggingTool::CosmicRayTaggingTool() :
    m_cutMode("nominal"),
    m_angularUncertainty(5.f),
    m_positionalUncertainty(3.f),
    m_maxAssociationDist(3.f * 18.f),
    m_minimumHits(15),
    m_inTimeMargin(5.f),
    m_inTimeMaxX0(1.f),
    m_marginY(20.f),
    m_marginZ(10.f),
    m_maxNeutrinoCosTheta(0.2f),
    m_minCosmicCosTheta(0.6f),
    m_maxCosmicCurvature(0.04f),
    m_face_Xa(0.f),
    m_face_Xc(0.f),
    m_face_Yb(0.f),
    m_face_Yt(0.f),
    m_face_Zu(0.f),
    m_face_Zd(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CosmicRayTaggingTool::~CosmicRayTaggingTool()
{
    m_pFile->Write();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTaggingTool::Initialize()
{
    if ("nominal" == m_cutMode)
    {
        // Nominal values set in constructor
    }
    else if ("cautious" == m_cutMode)
    {
        m_minCosmicCosTheta = std::numeric_limits<double>::max();
        m_maxCosmicCurvature = -std::numeric_limits<double>::max();
    }
    else if ("aggressive" == m_cutMode)
    {
        m_minCosmicCosTheta = 0.6f;
        m_maxCosmicCurvature = 0.1f;
    }
    else
    {
        std::cout << "CosmicRayTaggingTool - invalid cut mode " << m_cutMode << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    // Setup the output root file
    m_pFile = std::make_shared<TFile>("crTagging.root", "RECREATE");
    m_pTree = std::make_shared<TTree>("events", "events");

    m_pTree->Branch("truth_hasNu", &m_event.truth_hasNu);
    m_pTree->Branch("truth_nuPdgCode", &m_event.truth_nuPdgCode);
    m_pTree->Branch("truth_nuEnergy", &m_event.truth_nuEnergy);
    m_pTree->Branch("truth_nuVertexX", &m_event.truth_nuVertexX);
    m_pTree->Branch("truth_nuVertexY", &m_event.truth_nuVertexY);
    m_pTree->Branch("truth_nuVertexZ", &m_event.truth_nuVertexZ);
    m_pTree->Branch("p_truth_id", &m_event.p_truth_id);
    m_pTree->Branch("p_truth_pdgCode", &m_event.p_truth_pdgCode);
    m_pTree->Branch("p_truth_energy", &m_event.p_truth_energy);
    m_pTree->Branch("p_truth_momentumX", &m_event.p_truth_momentumX);
    m_pTree->Branch("p_truth_momentumY", &m_event.p_truth_momentumY);
    m_pTree->Branch("p_truth_momentumZ", &m_event.p_truth_momentumZ);
    m_pTree->Branch("p_truth_startX", &m_event.p_truth_startX);
    m_pTree->Branch("p_truth_startY", &m_event.p_truth_startY);
    m_pTree->Branch("p_truth_startZ", &m_event.p_truth_startZ);
    m_pTree->Branch("p_truth_endX", &m_event.p_truth_endX);
    m_pTree->Branch("p_truth_endY", &m_event.p_truth_endY);
    m_pTree->Branch("p_truth_endZ", &m_event.p_truth_endZ);
    m_pTree->Branch("p_truth_nHitsU", &m_event.p_truth_nHitsU);
    m_pTree->Branch("p_truth_nHitsV", &m_event.p_truth_nHitsV);
    m_pTree->Branch("p_truth_nHitsW", &m_event.p_truth_nHitsW);
    m_pTree->Branch("p_reco_id", &m_event.p_reco_id);
    m_pTree->Branch("p_reco_nHitsU", &m_event.p_reco_nHitsU);
    m_pTree->Branch("p_reco_nHitsV", &m_event.p_reco_nHitsV);
    m_pTree->Branch("p_reco_nHitsW", &m_event.p_reco_nHitsW);
    m_pTree->Branch("p_reco_canFit", &m_event.p_reco_canFit);
    m_pTree->Branch("p_reco_end1X", &m_event.p_reco_end1X);
    m_pTree->Branch("p_reco_end1Y", &m_event.p_reco_end1Y);
    m_pTree->Branch("p_reco_end1Z", &m_event.p_reco_end1Z);
    m_pTree->Branch("p_reco_end2X", &m_event.p_reco_end2X);
    m_pTree->Branch("p_reco_end2Y", &m_event.p_reco_end2Y);
    m_pTree->Branch("p_reco_end2Z", &m_event.p_reco_end2Z);
    m_pTree->Branch("p_reco_length", &m_event.p_reco_length);
    m_pTree->Branch("p_reco_cosTheta", &m_event.p_reco_cosTheta);
    m_pTree->Branch("p_reco_curvature", &m_event.p_reco_curvature);
    m_pTree->Branch("p_reco_inTime", &m_event.p_reco_inTime);
    m_pTree->Branch("p_reco_isContained", &m_event.p_reco_isContained);
    m_pTree->Branch("p_reco_isTopToBottom", &m_event.p_reco_isTopToBottom);
    m_pTree->Branch("p_reco_inNuSlice", &m_event.p_reco_inNuSlice);
    m_pTree->Branch("p_reco_isCRTagged", &m_event.p_reco_isCRTagged);
    m_pTree->Branch("p_reco_sliceId", &m_event.p_reco_sliceId);
    m_pTree->Branch("p_reco_linkedIds", &m_event.p_reco_linkedIds);
    m_pTree->Branch("p_reco_match_ids", &m_event.p_reco_match_ids);
    m_pTree->Branch("p_reco_match_hitsU", &m_event.p_reco_match_hitsU);
    m_pTree->Branch("p_reco_match_hitsV", &m_event.p_reco_match_hitsV);
    m_pTree->Branch("p_reco_match_hitsW", &m_event.p_reco_match_hitsW);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::FindAmbiguousPfos(const PfoList &parentCosmicRayPfos, PfoList &ambiguousPfos, const MasterAlgorithm *const pAlgorithm)
{
    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // TODO First time only, TODO Refactor with master algorithm
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);

    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    m_face_Xa = parentMinX;
    m_face_Xc = parentMaxX;
    m_face_Yb = parentMinY;
    m_face_Yt = parentMaxY;
    m_face_Zu = parentMinZ;
    m_face_Zd = parentMaxZ;

    PfoToPfoListMap pfoAssociationMap;
    this->GetPfoAssociations(parentCosmicRayPfos, pfoAssociationMap);

    PfoToSliceIdMap pfoToSliceIdMap;
    this->SliceEvent(parentCosmicRayPfos, pfoAssociationMap, pfoToSliceIdMap);

    CRCandidateList candidates;
    this->GetCRCandidates(parentCosmicRayPfos, pfoToSliceIdMap, candidates);

    PfoToBoolMap pfoToInTimeMap;
    this->CheckIfInTime(candidates, pfoToInTimeMap);

    PfoToBoolMap pfoToIsContainedMap;
    this->CheckIfContained(candidates, pfoToIsContainedMap );

    PfoToBoolMap pfoToIsTopToBottomMap;
    this->CheckIfTopToBottom(candidates, pfoToIsTopToBottomMap);

    UIntSet neutrinoSliceSet;
    this->GetNeutrinoSlices(candidates, pfoToInTimeMap, pfoToIsContainedMap, neutrinoSliceSet);

    PfoToBoolMap pfoToIsLikelyCRMuonMap;
    this->TagCRMuons(candidates, pfoToInTimeMap, pfoToIsTopToBottomMap, neutrinoSliceSet, pfoToIsLikelyCRMuonMap);

    //// BEGIN HACK ========================================================================================================================

    // Reset the output event so we can fill it again
    this->ResetOutputEvent();

    // Get the MCParticles
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

    // Find the true neutrino (if it exists)
    MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    // We have a true neutrino
    if (!trueNeutrinos.empty())
    {
        m_event.truth_hasNu = true;

        // If we have multiple neutrinos then find the one with the largest energy (shouldn't happen often - if ever depending on the sample)
        const MCParticle *pSelectedNeutrino = trueNeutrinos.front();
        for (const auto &pNeutrino : trueNeutrinos)
        {
            if (pNeutrino->GetEnergy() > pSelectedNeutrino->GetEnergy())
            {
                pSelectedNeutrino = pNeutrino;
            }
        }

        // Set the neutrino information
        m_event.truth_nuPdgCode = pSelectedNeutrino->GetParticleId();
        m_event.truth_nuEnergy = pSelectedNeutrino->GetEnergy();

        const auto nuVertex = pSelectedNeutrino->GetVertex();
        m_event.truth_nuVertexX = nuVertex.GetX();
        m_event.truth_nuVertexY = nuVertex.GetY();
        m_event.truth_nuVertexZ = nuVertex.GetZ();
    }
    // We don't have a true neutrino
    else
    {
        m_event.truth_hasNu = false;
    }

    // Get the hits
    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, "CaloHitList2D", pCaloHitList));

    // Identify neutrino induced MCParticles, and get mappings to their hits
    // ATTN here we have some requirement that the MCParticles are "reconstructable" - I think there's an implicit assumption that the
    // MCParticle produces at least one hit... so not quite the same as all MCParticles
    LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = 0;
    parameters.m_minHitsForGoodView = 0;
    parameters.m_minPrimaryGoodViews = 0;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);

    // Save the information about the neutrino induced MCParticles
    // First, put them in an ordered container and sort them for reproducibility
    MCParticleVector nuMCParticles;
    LArMonitoringHelper::GetOrderedMCParticleVector({nuMCParticlesToGoodHitsMap}, nuMCParticles);

    // Fill the MCParticle info
    std::unordered_map<const MCParticle *, int> mcParticleIdMap;
    for (size_t id = 0; id < nuMCParticles.size(); ++id)
    {
        const auto &pMCParticle = nuMCParticles.at(id);

        // Make a mapping from the MCParticle to the unique ID for later use
        mcParticleIdMap.emplace(pMCParticle, id);
        m_event.p_truth_id.push_back(id);

        m_event.p_truth_pdgCode.push_back(pMCParticle->GetParticleId());
        m_event.p_truth_energy.push_back(pMCParticle->GetEnergy());

        const auto momentum = pMCParticle->GetMomentum();
        m_event.p_truth_momentumX.push_back(momentum.GetX());
        m_event.p_truth_momentumY.push_back(momentum.GetY());
        m_event.p_truth_momentumZ.push_back(momentum.GetZ());

        const auto start = pMCParticle->GetVertex();
        m_event.p_truth_startX.push_back(start.GetX());
        m_event.p_truth_startY.push_back(start.GetY());
        m_event.p_truth_startZ.push_back(start.GetZ());

        const auto end = pMCParticle->GetEndpoint();
        m_event.p_truth_endX.push_back(end.GetX());
        m_event.p_truth_endY.push_back(end.GetY());
        m_event.p_truth_endZ.push_back(end.GetZ());

        const auto hits = nuMCParticlesToGoodHitsMap.at(pMCParticle);
        m_event.p_truth_nHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, hits));
        m_event.p_truth_nHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, hits));
        m_event.p_truth_nHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, hits));
    }

    // Get the hit sharing maps between Pfos and these primary MCParticles
    LArMCParticleHelper::MCContributionMapVector mcParticlesToGoodHitsMaps({nuMCParticlesToGoodHitsMap});
    LArMCParticleHelper::PfoContributionMap pfoToReconstructable2DHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(parentCosmicRayPfos, mcParticlesToGoodHitsMaps, pfoToReconstructable2DHitsMap);

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcParticleToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToReconstructable2DHitsMap, mcParticlesToGoodHitsMaps, pfoToMCParticleHitSharingMap, mcParticleToPfoHitSharingMap);

    // Get an ordered PFO vector so we can assigned a unique ID to each
    PfoVector orderedPfoVector;
    for (const auto &candidate : candidates)
    {
        orderedPfoVector.push_back(candidate.m_pPfo);
    }

    // Sort by the number of hits
    std::sort(orderedPfoVector.begin(), orderedPfoVector.end(), LArPfoHelper::SortByNHits);

    // Make a mapping from PFO to its ID
    std::unordered_map<const ParticleFlowObject *, int> pfoIdMap;
    for (size_t id = 0; id < orderedPfoVector.size(); ++id)
    {
        pfoIdMap.emplace(orderedPfoVector.at(id), id);
    }

    // Save the information about the reco particles
    for (const auto &candidate : candidates)
    {
        m_event.p_reco_id.push_back(pfoIdMap.at(candidate.m_pPfo));

        // Get the downstream PFOs
        PfoList downstreamPfos;
        LArPfoHelper::GetAllDownstreamPfos(candidate.m_pPfo, downstreamPfos);

        // Get the hits in this PFO (and it's daughters)
        CaloHitList pfoHitsU;
        LArPfoHelper::GetCaloHits(downstreamPfos, TPC_VIEW_U, pfoHitsU);
        LArPfoHelper::GetIsolatedCaloHits(downstreamPfos, TPC_VIEW_U, pfoHitsU);
        m_event.p_reco_nHitsU.push_back(pfoHitsU.size());

        CaloHitList pfoHitsV;
        LArPfoHelper::GetCaloHits(downstreamPfos, TPC_VIEW_V, pfoHitsV);
        LArPfoHelper::GetIsolatedCaloHits(downstreamPfos, TPC_VIEW_V, pfoHitsV);
        m_event.p_reco_nHitsV.push_back(pfoHitsV.size());

        CaloHitList pfoHitsW;
        LArPfoHelper::GetCaloHits(downstreamPfos, TPC_VIEW_W, pfoHitsW);
        LArPfoHelper::GetIsolatedCaloHits(downstreamPfos, TPC_VIEW_W, pfoHitsW);
        m_event.p_reco_nHitsW.push_back(pfoHitsW.size());

        // Get the fitted track information
        m_event.p_reco_canFit.push_back(candidate.m_canFit);
        if (candidate.m_canFit)
        {
            m_event.p_reco_end1X.push_back(candidate.m_endPoint1.GetX());
            m_event.p_reco_end1Y.push_back(candidate.m_endPoint1.GetY());
            m_event.p_reco_end1Z.push_back(candidate.m_endPoint1.GetZ());

            m_event.p_reco_end2X.push_back(candidate.m_endPoint2.GetX());
            m_event.p_reco_end2Y.push_back(candidate.m_endPoint2.GetY());
            m_event.p_reco_end2Z.push_back(candidate.m_endPoint2.GetZ());

            m_event.p_reco_length.push_back(candidate.m_length);
            m_event.p_reco_cosTheta.push_back(candidate.m_theta);
            m_event.p_reco_curvature.push_back(candidate.m_curvature);
        }
        else
        {
            // Store dummy track info
            m_event.p_reco_end1X.push_back(-std::numeric_limits<float>::max());
            m_event.p_reco_end1Y.push_back(-std::numeric_limits<float>::max());
            m_event.p_reco_end1Z.push_back(-std::numeric_limits<float>::max());
            m_event.p_reco_end2X.push_back(-std::numeric_limits<float>::max());
            m_event.p_reco_end2Y.push_back(-std::numeric_limits<float>::max());
            m_event.p_reco_end2Z.push_back(-std::numeric_limits<float>::max());
            m_event.p_reco_length.push_back(-std::numeric_limits<float>::max());
            m_event.p_reco_cosTheta.push_back(-std::numeric_limits<float>::max());
            m_event.p_reco_curvature.push_back(-std::numeric_limits<float>::max());
        }

        // Store the results of each of the "cuts"
        m_event.p_reco_inTime.push_back(pfoToInTimeMap.at(candidate.m_pPfo));
        m_event.p_reco_isContained.push_back(pfoToIsContainedMap.at(candidate.m_pPfo));
        m_event.p_reco_isTopToBottom.push_back(pfoToIsTopToBottomMap.at(candidate.m_pPfo));
        m_event.p_reco_inNuSlice.push_back(neutrinoSliceSet.count(candidate.m_sliceId) != 0);
        m_event.p_reco_isCRTagged.push_back(pfoToIsLikelyCRMuonMap.at(candidate.m_pPfo));

        // Store the associations to other PFOs
        m_event.p_reco_sliceId.push_back(candidate.m_sliceId);

        std::vector<int> linkedIds;
        const auto assocPfosIter = pfoAssociationMap.find(candidate.m_pPfo);
        if (assocPfosIter != pfoAssociationMap.end())
        {
            for (const auto &pLinkedPfo : assocPfosIter->second)
            {
                linkedIds.push_back(pfoIdMap.at(pLinkedPfo));
            }
        }
        m_event.p_reco_linkedIds.push_back(linkedIds);

        // Store the reco-truth matching information
        const auto assocMCIter = pfoToMCParticleHitSharingMap.find(candidate.m_pPfo);

        std::vector<int> matchIds, matchHitsU, matchHitsV, matchHitsW;
        if (assocMCIter != pfoToMCParticleHitSharingMap.end())
        {
            const auto mcParticleHitPairs = assocMCIter->second;
            for (const auto &[pMCParticle, sharedHits] : mcParticleHitPairs)
            {
                matchIds.push_back(mcParticleIdMap.at(pMCParticle));
                matchHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHits));
                matchHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHits));
                matchHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHits));
            }
        }

        m_event.p_reco_match_ids.push_back(matchIds);
        m_event.p_reco_match_hitsU.push_back(matchHitsU);
        m_event.p_reco_match_hitsV.push_back(matchHitsV);
        m_event.p_reco_match_hitsW.push_back(matchHitsW);
    }

    // Fill the tree!
    m_pTree->Fill();

    //// END HACK ==========================================================================================================================

    for (const ParticleFlowObject *const pPfo : parentCosmicRayPfos)
    {
        if (!pfoToIsLikelyCRMuonMap.at(pPfo))
            ambiguousPfos.push_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::ResetOutputEvent()
{
    // Set dummy values
    m_event.truth_hasNu = false;
    m_event.truth_nuPdgCode = -std::numeric_limits<int>::max();
    m_event.truth_nuEnergy = -std::numeric_limits<float>::max();
    m_event.truth_nuVertexX = -std::numeric_limits<float>::max();
    m_event.truth_nuVertexY = -std::numeric_limits<float>::max();
    m_event.truth_nuVertexZ = -std::numeric_limits<float>::max();

    // Clear any vectors
    m_event.p_truth_id.clear();
    m_event.p_truth_pdgCode.clear();
    m_event.p_truth_energy.clear();
    m_event.p_truth_momentumX.clear();
    m_event.p_truth_momentumY.clear();
    m_event.p_truth_momentumZ.clear();
    m_event.p_truth_startX.clear();
    m_event.p_truth_startY.clear();
    m_event.p_truth_startZ.clear();
    m_event.p_truth_endX.clear();
    m_event.p_truth_endY.clear();
    m_event.p_truth_endZ.clear();
    m_event.p_truth_nHitsU.clear();
    m_event.p_truth_nHitsV.clear();
    m_event.p_truth_nHitsW.clear();
    m_event.p_reco_id.clear();
    m_event.p_reco_nHitsU.clear();
    m_event.p_reco_nHitsV.clear();
    m_event.p_reco_nHitsW.clear();
    m_event.p_reco_canFit.clear();
    m_event.p_reco_end1X.clear();
    m_event.p_reco_end1Y.clear();
    m_event.p_reco_end1Z.clear();
    m_event.p_reco_end2X.clear();
    m_event.p_reco_end2Y.clear();
    m_event.p_reco_end2Z.clear();
    m_event.p_reco_length.clear();
    m_event.p_reco_cosTheta.clear();
    m_event.p_reco_curvature.clear();
    m_event.p_reco_inTime.clear();
    m_event.p_reco_isContained.clear();
    m_event.p_reco_isTopToBottom.clear();
    m_event.p_reco_inNuSlice.clear();
    m_event.p_reco_isCRTagged.clear();
    m_event.p_reco_sliceId.clear();
    m_event.p_reco_linkedIds.clear();
    m_event.p_reco_match_ids.clear();
    m_event.p_reco_match_hitsU.clear();
    m_event.p_reco_match_hitsV.clear();
    m_event.p_reco_match_hitsW.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayTaggingTool::GetValid3DCluster(const ParticleFlowObject *const pPfo, const Cluster *&pCluster3D) const
{
    pCluster3D = nullptr;
    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

    // ATTN Only uses first cluster with hits of type TPC_3D
    if (clusters3D.empty() || (clusters3D.front()->GetNCaloHits() < m_minimumHits))
        return false;

    pCluster3D = clusters3D.front();
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::GetPfoAssociations(const PfoList &parentCosmicRayPfos, PfoToPfoListMap &pfoAssociationMap) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pFirstLArTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float layerPitch(pFirstLArTPC->GetWirePitchW());

    PfoToSlidingFitsMap pfoToSlidingFitsMap;

    for (const ParticleFlowObject *const pPfo : parentCosmicRayPfos)
    {
        const pandora::Cluster *pCluster(nullptr);
        if (!this->GetValid3DCluster(pPfo, pCluster) || !pCluster)
            continue;

        (void) pfoToSlidingFitsMap.insert(PfoToSlidingFitsMap::value_type(pPfo, std::make_pair(
            ThreeDSlidingFitResult(pCluster, 5, layerPitch), ThreeDSlidingFitResult(pCluster, 100, layerPitch)))); // TODO Configurable
    }

    for (const ParticleFlowObject *const pPfo1 : parentCosmicRayPfos)
    {
        PfoToSlidingFitsMap::const_iterator iter1(pfoToSlidingFitsMap.find(pPfo1));
        if (pfoToSlidingFitsMap.end() == iter1)
            continue;

        const ThreeDSlidingFitResult &fitPos1(iter1->second.first), &fitDir1(iter1->second.second);

        for (const ParticleFlowObject *const pPfo2 : parentCosmicRayPfos)
        {
            if (pPfo1 == pPfo2)
                continue;

            PfoToSlidingFitsMap::const_iterator iter2(pfoToSlidingFitsMap.find(pPfo2));
            if (pfoToSlidingFitsMap.end() == iter2)
                continue;

            const ThreeDSlidingFitResult &fitPos2(iter2->second.first), &fitDir2(iter2->second.second);

            // TODO Use existing LArPointingClusters and IsEmission/IsNode logic, for consistency
            if (!(this->CheckAssociation(fitPos1.GetGlobalMinLayerPosition(), fitDir1.GetGlobalMinLayerDirection() * -1.f, fitPos2.GetGlobalMinLayerPosition(), fitDir2.GetGlobalMinLayerDirection() * -1.f) ||
                this->CheckAssociation(fitPos1.GetGlobalMinLayerPosition(), fitDir1.GetGlobalMinLayerDirection() * -1.f, fitPos2.GetGlobalMaxLayerPosition(), fitDir2.GetGlobalMaxLayerDirection()) ||
                this->CheckAssociation(fitPos1.GetGlobalMaxLayerPosition(), fitDir1.GetGlobalMaxLayerDirection(), fitPos2.GetGlobalMinLayerPosition(), fitDir2.GetGlobalMinLayerDirection() * -1.f) ||
                this->CheckAssociation(fitPos1.GetGlobalMaxLayerPosition(), fitDir1.GetGlobalMaxLayerDirection(), fitPos2.GetGlobalMaxLayerPosition(), fitDir2.GetGlobalMaxLayerDirection())))
            {
                continue;
            }

            PfoList &pfoList1(pfoAssociationMap[pPfo1]), &pfoList2(pfoAssociationMap[pPfo2]);

            if (pfoList1.end() == std::find(pfoList1.begin(), pfoList1.end(), pPfo2))
                pfoList1.push_back(pPfo2);

            if (pfoList2.end() == std::find(pfoList2.begin(), pfoList2.end(), pPfo1))
                pfoList2.push_back(pPfo1);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayTaggingTool::CheckAssociation(const CartesianVector &endPoint1, const CartesianVector &endDir1, const CartesianVector &endPoint2,
    const CartesianVector &endDir2) const
{
    // TODO This function needs significant tidying and checks (variable names must be self describing, check deltaTheta, etc.)
    const CartesianVector n(endDir1.GetUnitVector());
    const CartesianVector m(endDir2.GetUnitVector());
    const CartesianVector a(endPoint2 - endPoint1);
    const float b(n.GetDotProduct(m));

    // Parallel lines never meet
    if (std::fabs(b - 1.f) < std::numeric_limits<float>::epsilon())
        return false;

    // Distance from endPoint1 along endDir1 to the point of closest approach
    const float lambda((n - m * b).GetDotProduct(a) / (1.f - b * b));

    // Distance from endPoint2 along endDir2 to the point of closest approach
    const float mu((n * b - m).GetDotProduct(a) / (1.f - b * b));

    // Calculate the maximum "vertex uncertainty"
    const float deltaTheta(m_angularUncertainty * M_PI / 180.f);
    const float maxVertexUncertainty(m_maxAssociationDist * std::sin(deltaTheta) + m_positionalUncertainty);

    // Ensure that the distances to the point of closest approch are within the limits
    if ((lambda < -maxVertexUncertainty) || (mu < -maxVertexUncertainty) ||
        (lambda > m_maxAssociationDist + maxVertexUncertainty) || (mu > m_maxAssociationDist + maxVertexUncertainty))
    {
        return false;
    }

    // Ensure point of closest approach is within detector
    const CartesianVector impactPosition((endPoint1 + n * lambda + endPoint2 + m * mu) * 0.5f);

    if ((impactPosition.GetX() < m_face_Xa - maxVertexUncertainty) || (impactPosition.GetX() > m_face_Xc + maxVertexUncertainty) ||
        (impactPosition.GetY() < m_face_Yb - maxVertexUncertainty) || (impactPosition.GetY() > m_face_Yt + maxVertexUncertainty) ||
        (impactPosition.GetZ() < m_face_Zu - maxVertexUncertainty) || (impactPosition.GetZ() > m_face_Zd + maxVertexUncertainty))
    {
        return false;
    }

    // Compare distance of closest approach and max uncertainty in impact parameter
    const float maxImpactDist(std::sin(deltaTheta) * (std::fabs(mu) + std::fabs(lambda)) + m_positionalUncertainty);
    const CartesianVector d(a - n * lambda + m * mu);

    return (d.GetMagnitude() < maxImpactDist);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::SliceEvent(const PfoList &parentCosmicRayPfos, const PfoToPfoListMap &pfoAssociationMap, PfoToSliceIdMap &pfoToSliceIdMap) const
{
    SliceList sliceList;

    for (const ParticleFlowObject *const pPfo : parentCosmicRayPfos)
    {
        bool isAlreadyInSlice(false);

        for (const PfoList &slice : sliceList)
        {
            if (std::find(slice.begin(), slice.end(), pPfo) != slice.end())
            {
                isAlreadyInSlice = true;
                break;
            }
        }

        if (!isAlreadyInSlice)
        {
            sliceList.push_back(PfoList());
            this->FillSlice(pPfo, pfoAssociationMap, sliceList.back());
        }
    }

    unsigned int sliceId(0);
    for (const PfoList &slice : sliceList)
    {
        for (const ParticleFlowObject *const pPfo : slice)
        {
            if (!pfoToSliceIdMap.insert(PfoToSliceIdMap::value_type(pPfo, sliceId)).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }

        ++sliceId;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::FillSlice(const ParticleFlowObject *const pPfo, const PfoToPfoListMap &pfoAssociationMap, PfoList &slice) const
{
    if (std::find(slice.begin(), slice.end(), pPfo) != slice.end())
        return;

    slice.push_back(pPfo);

    PfoToPfoListMap::const_iterator iter(pfoAssociationMap.find(pPfo));

    if (pfoAssociationMap.end() != iter)
    {
        for (const ParticleFlowObject *const pAssociatedPfo : iter->second)
            this->FillSlice(pAssociatedPfo, pfoAssociationMap, slice);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::GetCRCandidates(const PfoList &parentCosmicRayPfos, const PfoToSliceIdMap &pfoToSliceIdMap, CRCandidateList &candidates) const
{
    for (const ParticleFlowObject *const pPfo : parentCosmicRayPfos)
    {
        if (!LArPfoHelper::IsFinalState(pPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        candidates.push_back(CRCandidate(this->GetPandora(), pPfo, pfoToSliceIdMap.at(pPfo)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::CheckIfInTime(const CRCandidateList &candidates, PfoToBoolMap &pfoToInTimeMap) const
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());

    for (const CRCandidate &candidate : candidates)
    {
        // Cosmic-ray muons extending outside of (single) physical volume if given t0 is that of the beam particle
        float minX(std::numeric_limits<float>::max()), maxX(-std::numeric_limits<float>::max());

        if (candidate.m_canFit)
        {
            minX = ((candidate.m_endPoint1.GetX() < candidate.m_endPoint2.GetX()) ? candidate.m_endPoint1.GetX() : candidate.m_endPoint2.GetX());
            maxX = ((candidate.m_endPoint1.GetX() > candidate.m_endPoint2.GetX()) ? candidate.m_endPoint1.GetX() : candidate.m_endPoint2.GetX());
        }
        else
        {
            // Handle any particles with small numbers of 3D hits, for which no 3D sliding fit information is available
            for (const Cluster *const pCluster : candidate.m_pPfo->GetClusterList())
            {
                float clusterMinX(std::numeric_limits<float>::max()), clusterMaxX(-std::numeric_limits<float>::max());
                LArClusterHelper::GetClusterSpanX(pCluster, clusterMinX, clusterMaxX);
                minX = std::min(clusterMinX, minX);
                maxX = std::max(clusterMaxX, maxX);
            }
        }

        bool isInTime((minX > m_face_Xa - m_inTimeMargin) && (maxX < m_face_Xc + m_inTimeMargin));

        // Cosmic-ray muons that have been shifted and stitched across mid plane between volumes
        if (isInTime)
        {
            try
            {
                if (std::fabs(LArPfoHelper::GetVertex(candidate.m_pPfo)->GetX0()) > m_inTimeMaxX0)
                    isInTime = false;
            }
            catch (const StatusCodeException &) {}
        }

        // Cosmic-ray muons extending outside of (any individual) physical volume if given t0 is that of the beam particle
        if (isInTime)
        {
            CaloHitList caloHitList;
            LArPfoHelper::GetCaloHits(candidate.m_pPfo, TPC_VIEW_U, caloHitList);
            LArPfoHelper::GetCaloHits(candidate.m_pPfo, TPC_VIEW_V, caloHitList);
            LArPfoHelper::GetCaloHits(candidate.m_pPfo, TPC_VIEW_W, caloHitList);

            bool isFirstHit(true);
            bool isInSingleVolume(true);
            unsigned int volumeId(std::numeric_limits<unsigned int>::max());

            for (const CaloHit *const pCaloHit : caloHitList)
            {
                const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit*>(pCaloHit));

                if (!pLArCaloHit)
                    continue;

                if (isFirstHit)
                {
                    isFirstHit = false;
                    volumeId = pLArCaloHit->GetLArTPCVolumeId();
                }
                else if (volumeId != pLArCaloHit->GetLArTPCVolumeId())
                {
                    isInSingleVolume = false;
                    break;
                }
            }

            LArTPCMap::const_iterator tpcIter(larTPCMap.find(volumeId));

            if (isInSingleVolume && (larTPCMap.end() != tpcIter))
            {
                const float thisFaceXLow(tpcIter->second->GetCenterX() - 0.5f * tpcIter->second->GetWidthX());
                const float thisFaceXHigh(tpcIter->second->GetCenterX() + 0.5f * tpcIter->second->GetWidthX());

                if (!((minX > thisFaceXLow - m_inTimeMargin) && (maxX < thisFaceXHigh + m_inTimeMargin)))
                    isInTime = false;
            }
        }

        if (!pfoToInTimeMap.insert(PfoToBoolMap::value_type(candidate.m_pPfo, isInTime)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::CheckIfContained(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsContainedMap) const
{
    for (const CRCandidate &candidate : candidates)
    {
        const float upperY((candidate.m_endPoint1.GetY() > candidate.m_endPoint2.GetY()) ? candidate.m_endPoint1.GetY() : candidate.m_endPoint2.GetY());
        const float lowerY((candidate.m_endPoint1.GetY() < candidate.m_endPoint2.GetY()) ? candidate.m_endPoint1.GetY() : candidate.m_endPoint2.GetY());

        const float zAtUpperY((candidate.m_endPoint1.GetY() > candidate.m_endPoint2.GetY()) ? candidate.m_endPoint1.GetZ() : candidate.m_endPoint2.GetZ());
        const float zAtLowerY((candidate.m_endPoint1.GetY() < candidate.m_endPoint2.GetY()) ? candidate.m_endPoint1.GetZ() : candidate.m_endPoint2.GetZ());

        const bool isContained((upperY < m_face_Yt - m_marginY) && (upperY > m_face_Yb + m_marginY) &&
            (lowerY < m_face_Yt - m_marginY) && (lowerY > m_face_Yb + m_marginY) &&
            (zAtUpperY < m_face_Zd - m_marginZ) && (zAtUpperY > m_face_Zu + m_marginZ) &&
            (zAtLowerY < m_face_Zd - m_marginZ) && (zAtLowerY > m_face_Zu + m_marginZ));

        if (!pfoToIsContainedMap.insert(PfoToBoolMap::value_type(candidate.m_pPfo, isContained)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::CheckIfTopToBottom(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsTopToBottomMap) const
{
    for (const CRCandidate &candidate : candidates)
    {
        const float upperY((candidate.m_endPoint1.GetY() > candidate.m_endPoint2.GetY()) ? candidate.m_endPoint1.GetY() : candidate.m_endPoint2.GetY());
        const float lowerY((candidate.m_endPoint1.GetY() < candidate.m_endPoint2.GetY()) ? candidate.m_endPoint1.GetY() : candidate.m_endPoint2.GetY());

        const bool isTopToBottom((upperY > m_face_Yt - m_marginY) && (lowerY < m_face_Yb + m_marginY));

        if (!pfoToIsTopToBottomMap.insert(PfoToBoolMap::value_type(candidate.m_pPfo, isTopToBottom)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::GetNeutrinoSlices(const CRCandidateList &candidates, const PfoToBoolMap &pfoToInTimeMap, const PfoToBoolMap &pfoToIsContainedMap,
    UIntSet &neutrinoSliceSet) const
{
    IntBoolMap sliceIdToIsInTimeMap;

    for (const CRCandidate &candidate : candidates)
    {
        if (sliceIdToIsInTimeMap.find(candidate.m_sliceId) == sliceIdToIsInTimeMap.end())
            sliceIdToIsInTimeMap.insert(std::make_pair(candidate.m_sliceId, true));

        if (!pfoToInTimeMap.at(candidate.m_pPfo))
            sliceIdToIsInTimeMap.at(candidate.m_sliceId) = false;
    }

    for (const CRCandidate &candidate : candidates)
    {
        if (neutrinoSliceSet.count(candidate.m_sliceId))
            continue;

        const bool likelyNeutrino(candidate.m_canFit && sliceIdToIsInTimeMap.at(candidate.m_sliceId) &&
            (candidate.m_theta < m_maxNeutrinoCosTheta || pfoToIsContainedMap.at(candidate.m_pPfo)));

        if (likelyNeutrino)
            (void) neutrinoSliceSet.insert(candidate.m_sliceId);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::TagCRMuons(const CRCandidateList &candidates, const PfoToBoolMap &pfoToInTimeMap, const PfoToBoolMap &pfoToIsTopToBottomMap,
    const UIntSet &neutrinoSliceSet, PfoToBoolMap &pfoToIsLikelyCRMuonMap) const
{
    for (const CRCandidate &candidate : candidates)
    {
        const bool likelyCRMuon(!neutrinoSliceSet.count(candidate.m_sliceId) && (!pfoToInTimeMap.at(candidate.m_pPfo) || (candidate.m_canFit &&
            (pfoToIsTopToBottomMap.at(candidate.m_pPfo) || ((candidate.m_theta > m_minCosmicCosTheta) && (candidate.m_curvature < m_maxCosmicCurvature)))) ));

        if (!pfoToIsLikelyCRMuonMap.insert(PfoToBoolMap::value_type(candidate.m_pPfo, likelyCRMuon)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

CosmicRayTaggingTool::CRCandidate::CRCandidate(const Pandora &pandora, const ParticleFlowObject *const pPfo, const unsigned int sliceId) :
    m_pPfo(pPfo),
    m_sliceId(sliceId),
    m_canFit(false),
    m_endPoint1(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
    m_endPoint2(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
    m_length(std::numeric_limits<float>::max()),
    m_curvature(std::numeric_limits<float>::max()),
    m_theta(std::numeric_limits<float>::max())
{
    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

    if (!clusters3D.empty() && (clusters3D.front()->GetNCaloHits() > 15)) // TODO Configurable
    {
        m_canFit = true;
        const LArTPC *const pFirstLArTPC(pandora.GetGeometry()->GetLArTPCMap().begin()->second);
        const ThreeDSlidingFitResult slidingFitResult(clusters3D.front(), 5, pFirstLArTPC->GetWirePitchW()); // TODO Configurable
        this->CalculateFitVariables(slidingFitResult);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingTool::CRCandidate::CalculateFitVariables(const ThreeDSlidingFitResult &slidingFitResult)
{
    m_endPoint1 = slidingFitResult.GetGlobalMinLayerPosition();
    m_endPoint2 = slidingFitResult.GetGlobalMaxLayerPosition();
    m_length = (m_endPoint2 - m_endPoint1).GetMagnitude();

    if (std::fabs(m_length) > std::numeric_limits<float>::epsilon())
        m_theta = std::fabs(m_endPoint2.GetY() - m_endPoint1.GetY()) / m_length;

    const float layerPitch(slidingFitResult.GetFirstFitResult().GetLayerPitch());

    CartesianPointVector directionList;
    for (int i = slidingFitResult.GetMinLayer(); i < slidingFitResult.GetMaxLayer(); ++i)
    {
        CartesianVector direction(0.f, 0.f, 0.f);
        if (STATUS_CODE_SUCCESS == slidingFitResult.GetGlobalFitDirection(static_cast<float>(i) * layerPitch, direction))
            directionList.push_back(direction);
    }

    CartesianVector meanDirection(0.f, 0.f, 0.f);
    for (const CartesianVector &direction : directionList)
        meanDirection += direction;

    if (!directionList.empty() > 0)
        meanDirection *= 1.f / static_cast<float>(directionList.size());

    m_curvature = 0.f;
    for (const CartesianVector &direction : directionList)
        m_curvature += (direction - meanDirection).GetMagnitude();

    if (!directionList.empty() > 0)
        m_curvature *= 1.f / static_cast<float>(directionList.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTaggingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CutMode", m_cutMode));
    std::transform(m_cutMode.begin(), m_cutMode.end(), m_cutMode.begin(), ::tolower);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AngularUncertainty", m_angularUncertainty ));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PositionalUncertainty", m_positionalUncertainty ));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAssociationDist", m_maxAssociationDist));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HitThreshold", m_minimumHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InTimeMargin", m_inTimeMargin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InTimeMaxX0", m_inTimeMaxX0));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MarginY", m_marginY));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MarginZ", m_marginZ));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxNeutrinoCosTheta", m_maxNeutrinoCosTheta));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosmicCosTheta", m_minCosmicCosTheta));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxCosmicCurvature", m_maxCosmicCurvature));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
