/**
 *  @file   larpandoracontent/LArControlFlow/NeutrinoIdTool.cc
 *
 *  @brief  Implementation of the neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/NeutrinoIdTool.h"

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

NeutrinoIdTool::NeutrinoIdTool() :
    m_useTrainingMode(false),
    m_selectNuanceCode(false),
    m_nuance(-std::numeric_limits<int>::max()),
    m_minPurity(0.9f),
    m_minCompleteness(0.9f),
    m_minProbability(0.0f),
    m_maxNeutrinos(1),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH")
{
    // Setup the output root file
    m_pFile = std::make_shared<TFile>("neutrinoID.root", "RECREATE");
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
    m_pTree->Branch("nHitsU", &m_event.nHitsU);
    m_pTree->Branch("nHitsV", &m_event.nHitsV);
    m_pTree->Branch("nHitsW", &m_event.nHitsW);
    m_pTree->Branch("nNuHitsU", &m_event.nNuHitsU);
    m_pTree->Branch("nNuHitsV", &m_event.nNuHitsV);
    m_pTree->Branch("nNuHitsW", &m_event.nNuHitsW);
    m_pTree->Branch("s_nHitsU", &m_event.s_nHitsU);
    m_pTree->Branch("s_nHitsV", &m_event.s_nHitsV);
    m_pTree->Branch("s_nHitsW", &m_event.s_nHitsW);
    m_pTree->Branch("s_nNuHitsU", &m_event.s_nNuHitsU);
    m_pTree->Branch("s_nNuHitsV", &m_event.s_nNuHitsV);
    m_pTree->Branch("s_nNuHitsW", &m_event.s_nNuHitsW);
    m_pTree->Branch("s_areFeaturesAvailable", &m_event.s_areFeaturesAvailable);
    m_pTree->Branch("s_nuNFinalStatePfos", &m_event.s_nuNFinalStatePfos);
    m_pTree->Branch("s_nuNHitsTotal", &m_event.s_nuNHitsTotal);
    m_pTree->Branch("s_nuVertexY", &m_event.s_nuVertexY);
    m_pTree->Branch("s_nuWeightedDirZ", &m_event.s_nuWeightedDirZ);
    m_pTree->Branch("s_nuNSpacePointsInSphere", &m_event.s_nuNSpacePointsInSphere);
    m_pTree->Branch("s_nuEigenRatioInSphere", &m_event.s_nuEigenRatioInSphere);
    m_pTree->Branch("s_crLongestTrackDirY", &m_event.s_crLongestTrackDirY);
    m_pTree->Branch("s_crLongestTrackDeflection", &m_event.s_crLongestTrackDeflection);
    m_pTree->Branch("s_crFracHitsInLongestTrack", &m_event.s_crFracHitsInLongestTrack);
    m_pTree->Branch("s_nCRHitsMax", &m_event.s_nCRHitsMax);
    m_pTree->Branch("s_score", &m_event.s_score);
}

//------------------------------------------------------------------------------------------------------------------------------------------

NeutrinoIdTool::~NeutrinoIdTool()
{
    m_pFile->Write();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SelectOutputPfos(const Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const unsigned int nSlices(nuSliceHypotheses.size());
    if (nSlices == 0) return;

    SliceFeaturesVector sliceFeaturesVector;
    this->GetSliceFeatures(this, nuSliceHypotheses, crSliceHypotheses, sliceFeaturesVector);

    if (m_useTrainingMode)
    {
        // ATTN in training mode, just return everything as a cosmic-ray
        this->SelectAllPfos(pAlgorithm, crSliceHypotheses, selectedPfos);

        unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());
        if (!this->GetBestMCSliceIndex(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, bestSliceIndex)) return;

        for (unsigned int sliceIndex = 0; sliceIndex < nSlices; ++sliceIndex)
        {
            const SliceFeatures &features(sliceFeaturesVector.at(sliceIndex));
            if (!features.IsFeatureVectorAvailable()) continue;

            LArMvaHelper::MvaFeatureVector featureVector;
            features.GetFeatureVector(featureVector);
            LArMvaHelper::ProduceTrainingExample(m_trainingOutputFile, sliceIndex == bestSliceIndex, featureVector);
        }

        return;
    }

    //// BEGIN HACK ========================================================================================================================

    // Reset the output event so we can fill it again
    this->ResetOutputEvent();

    // Get all hits
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pCaloHitList));

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

    for (size_t id = 0; id < nuMCParticles.size(); ++id)
    {
        const auto &pMCParticle = nuMCParticles.at(id);

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

    // Put the hits in the required format
    CaloHitList caloHitList;
    caloHitList.insert(caloHitList.end(), pCaloHitList->begin(), pCaloHitList->end());

    m_event.nHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, caloHitList);
    m_event.nHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, caloHitList);
    m_event.nHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, caloHitList);

    // Get the total number of neutrino induced hits in the event
    CaloHitList nuHitsTotal;
    this->GetNeutrinoInducedHits(caloHitList, nuHitsTotal);

    m_event.nNuHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, nuHitsTotal);
    m_event.nNuHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, nuHitsTotal);
    m_event.nNuHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, nuHitsTotal);

    // Loop over the slices
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        // Collect up the hits from the CR and neutrino hypotheses
        CaloHitList caloHitsInSlice;
        this->CollectHitsInSlice(nuSliceHypotheses.at(sliceIndex), crSliceHypotheses.at(sliceIndex), caloHitsInSlice);

        // Get the total number of hits in the slice
        m_event.s_nHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, caloHitsInSlice));
        m_event.s_nHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, caloHitsInSlice));
        m_event.s_nHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, caloHitsInSlice));

        // Get the total number of neutrino induced hits in this slice
        CaloHitList nuHitsInSlice;
        this->GetNeutrinoInducedHits(caloHitsInSlice, nuHitsInSlice);

        m_event.s_nNuHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, nuHitsInSlice));
        m_event.s_nNuHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, nuHitsInSlice));
        m_event.s_nNuHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, nuHitsInSlice));

        // Determine if we can access the features for this slice
        const auto sliceFeatures = sliceFeaturesVector.at(sliceIndex);
        const bool areFeaturesAvailable = sliceFeatures.IsFeatureVectorAvailable();
        m_event.s_areFeaturesAvailable.push_back(areFeaturesAvailable);

        if (areFeaturesAvailable)
        {
            // If we can access the features, then get them!
            LArMvaHelper::MvaFeatureVector featureVector;
            sliceFeatures.GetFeatureVector(featureVector);

            // ATTN have to match the order filled!
            if (featureVector.size() != 10)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const float nuNFinalStatePfos = featureVector.at(0).Get();
            const float nuNHitsTotal = featureVector.at(1).Get();
            const float nuVertexY = featureVector.at(2).Get();
            const float nuWeightedDirZ = featureVector.at(3).Get();
            const float nuNSpacePointsInSphere = featureVector.at(4).Get();
            const float nuEigenRatioInSphere = featureVector.at(5).Get();
            const float crLongestTrackDirY = featureVector.at(6).Get();
            const float crLongestTrackDeflection = featureVector.at(7).Get();
            const float crFracHitsInLongestTrack = featureVector.at(8).Get();
            const float nCRHitsMax = featureVector.at(9).Get();

            m_event.s_nuNFinalStatePfos.push_back(nuNFinalStatePfos);
            m_event.s_nuNHitsTotal.push_back(nuNHitsTotal);
            m_event.s_nuVertexY.push_back(nuVertexY);
            m_event.s_nuWeightedDirZ.push_back(nuWeightedDirZ);
            m_event.s_nuNSpacePointsInSphere.push_back(nuNSpacePointsInSphere);
            m_event.s_nuEigenRatioInSphere.push_back(nuEigenRatioInSphere);
            m_event.s_crLongestTrackDirY.push_back(crLongestTrackDirY);
            m_event.s_crLongestTrackDeflection.push_back(crLongestTrackDeflection);
            m_event.s_crFracHitsInLongestTrack.push_back(crFracHitsInLongestTrack);
            m_event.s_nCRHitsMax.push_back(nCRHitsMax);
        }
        else
        {
            // Features aren't available so set dummy values
            m_event.s_nuNFinalStatePfos.push_back(-std::numeric_limits<float>::max());
            m_event.s_nuNHitsTotal.push_back(-std::numeric_limits<float>::max());
            m_event.s_nuVertexY.push_back(-std::numeric_limits<float>::max());
            m_event.s_nuWeightedDirZ.push_back(-std::numeric_limits<float>::max());
            m_event.s_nuNSpacePointsInSphere.push_back(-std::numeric_limits<float>::max());
            m_event.s_nuEigenRatioInSphere.push_back(-std::numeric_limits<float>::max());
            m_event.s_crLongestTrackDirY.push_back(-std::numeric_limits<float>::max());
            m_event.s_crLongestTrackDeflection.push_back(-std::numeric_limits<float>::max());
            m_event.s_crFracHitsInLongestTrack.push_back(-std::numeric_limits<float>::max());
            m_event.s_nCRHitsMax.push_back(-std::numeric_limits<float>::max());
        }

        // Get the neutrino-vs-cosmic score
        const auto score = sliceFeatures.GetNeutrinoProbability(m_supportVectorMachine);
        m_event.s_score.push_back(score);
    }

    // Fill the tree!
    m_pTree->Fill();

    //// END HACK ==========================================================================================================================

    this->SelectPfosByProbability(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, sliceFeaturesVector, selectedPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::ResetOutputEvent()
{
    // Set dummy values
    m_event.truth_hasNu = false;
    m_event.truth_nuPdgCode = -std::numeric_limits<int>::max();
    m_event.truth_nuEnergy = -std::numeric_limits<float>::max();
    m_event.truth_nuVertexX = -std::numeric_limits<float>::max();
    m_event.truth_nuVertexY = -std::numeric_limits<float>::max();
    m_event.truth_nuVertexZ = -std::numeric_limits<float>::max();

    m_event.nHitsU = -std::numeric_limits<int>::max();
    m_event.nHitsV = -std::numeric_limits<int>::max();
    m_event.nHitsW = -std::numeric_limits<int>::max();

    m_event.nNuHitsU = -std::numeric_limits<int>::max();
    m_event.nNuHitsV = -std::numeric_limits<int>::max();
    m_event.nNuHitsW = -std::numeric_limits<int>::max();

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
    m_event.s_nHitsU.clear();
    m_event.s_nHitsV.clear();
    m_event.s_nHitsW.clear();
    m_event.s_nNuHitsU.clear();
    m_event.s_nNuHitsV.clear();
    m_event.s_nNuHitsW.clear();
    m_event.s_areFeaturesAvailable.clear();
    m_event.s_nuNFinalStatePfos.clear();
    m_event.s_nuNHitsTotal.clear();
    m_event.s_nuVertexY.clear();
    m_event.s_nuWeightedDirZ.clear();
    m_event.s_nuNSpacePointsInSphere.clear();
    m_event.s_nuEigenRatioInSphere.clear();
    m_event.s_crLongestTrackDirY.clear();
    m_event.s_crLongestTrackDeflection.clear();
    m_event.s_crFracHitsInLongestTrack.clear();
    m_event.s_nCRHitsMax.clear();
    m_event.s_score.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::CollectHitsInSlice(const PfoList &nuSliceHypothesis, const PfoList &crSliceHypothesis, CaloHitList &caloHitList) const
{
    // Collect up all reconstructed particles in the slice under both hypotheses
    PfoList downstreamPfos;
    LArPfoHelper::GetAllDownstreamPfos(nuSliceHypothesis, downstreamPfos);
    LArPfoHelper::GetAllDownstreamPfos(crSliceHypothesis, downstreamPfos);

    // Collect up the hits in each view
    CaloHitList collectedHits;
    LArPfoHelper::GetCaloHits(downstreamPfos, TPC_VIEW_U, collectedHits);
    LArPfoHelper::GetCaloHits(downstreamPfos, TPC_VIEW_V, collectedHits);
    LArPfoHelper::GetCaloHits(downstreamPfos, TPC_VIEW_W, collectedHits);

    // ATTN the hits from the slice hypotheses are copies of those from the master algorithm, so we need to get their parent. These can be
    // repeated because we reconstruct everything twice, so we have to check for double counting
    for (const auto &pCaloHit : collectedHits)
    {
        const auto pParentHit = static_cast<const CaloHit *const>(pCaloHit->GetParentAddress());

        if (!pParentHit)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        if (std::find(caloHitList.begin(), caloHitList.end(), pParentHit) == caloHitList.end())
            caloHitList.push_back(pParentHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::GetSliceFeatures(const NeutrinoIdTool *const pTool, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, SliceFeaturesVector &sliceFeaturesVector) const
{
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
        sliceFeaturesVector.push_back(SliceFeatures(nuSliceHypotheses.at(sliceIndex), crSliceHypotheses.at(sliceIndex), pTool));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::GetBestMCSliceIndex(const Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, unsigned int &bestSliceIndex) const
{
    unsigned int nHitsInBestSlice(0), nNuHitsInBestSlice(0);

    // Get all hits in all slices to find true number of mc hits
    const CaloHitList *pAllReconstructedCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pAllReconstructedCaloHitList));

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

    // Obtain map: [mc particle -> primary mc particle]
    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);

    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    CaloHitList reconstructableCaloHitList;
    LArMCParticleHelper::PrimaryParameters parameters;
    LArMCParticleHelper::SelectCaloHits(pAllReconstructedCaloHitList, mcToPrimaryMCMap, reconstructableCaloHitList, parameters.m_selectInputHits, parameters.m_maxPhotonPropagation);

    const int nuNHitsTotal(this->CountNeutrinoInducedHits(reconstructableCaloHitList));
    const CaloHitSet reconstructableCaloHitSet(reconstructableCaloHitList.begin(), reconstructableCaloHitList.end());

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        CaloHitList reconstructedCaloHitList;
        this->Collect2DHits(crSliceHypotheses.at(sliceIndex), reconstructedCaloHitList, reconstructableCaloHitSet);

        for (const ParticleFlowObject *const pNeutrino : nuSliceHypotheses.at(sliceIndex))
        {
            const PfoList &nuFinalStates(pNeutrino->GetDaughterPfoList());
            this->Collect2DHits(nuFinalStates, reconstructedCaloHitList, reconstructableCaloHitSet);
        }

        const unsigned int nNuHits(this->CountNeutrinoInducedHits(reconstructedCaloHitList));

        if (nNuHits > nNuHitsInBestSlice)
        {
            nNuHitsInBestSlice = nNuHits;
            nHitsInBestSlice = reconstructedCaloHitList.size();
            bestSliceIndex = sliceIndex;
        }
    }

    // ATTN for events with no neutrino induced hits, default neutrino purity and completeness to zero
    const float purity(nHitsInBestSlice > 0 ? static_cast<float>(nNuHitsInBestSlice) / static_cast<float>(nHitsInBestSlice) : 0.f);
    const float completeness(nuNHitsTotal > 0 ? static_cast<float>(nNuHitsInBestSlice) / static_cast<float>(nuNHitsTotal) : 0.f);

    return this->PassesQualityCuts(pAlgorithm, purity, completeness);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::PassesQualityCuts(const Algorithm *const pAlgorithm, const float purity, const float completeness) const
{
    if (purity < m_minPurity || completeness < m_minCompleteness) return false;
    if (m_selectNuanceCode && (this->GetNuanceCode(pAlgorithm) != m_nuance)) return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::Collect2DHits(const PfoList &pfos, CaloHitList &reconstructedCaloHitList, const CaloHitSet &reconstructableCaloHitSet) const
{
    CaloHitList collectedHits;
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_U, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_V, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_W, collectedHits);

    for (const CaloHit *const pCaloHit : collectedHits)
    {
        if (!reconstructableCaloHitSet.count(pCaloHit))
            continue;

        // Ensure no hits have been double counted
        if (std::find(reconstructedCaloHitList.begin(), reconstructedCaloHitList.end(), pCaloHit) == reconstructedCaloHitList.end())
            reconstructedCaloHitList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int NeutrinoIdTool::CountNeutrinoInducedHits(const CaloHitList &caloHitList) const
{
    CaloHitList neutrinoHits;
    this->GetNeutrinoInducedHits(caloHitList, neutrinoHits);

    return neutrinoHits.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::GetNeutrinoInducedHits(const CaloHitList &caloHitList, CaloHitList &neutrinoHits) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
            if (LArMCParticleHelper::IsNeutrino(LArMCParticleHelper::GetParentMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit))))
            {
                neutrinoHits.push_back(pCaloHit);
            }
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int NeutrinoIdTool::GetNuanceCode(const Algorithm *const pAlgorithm) const
{
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

    MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    if (trueNeutrinos.size() != 1)
    {
        std::cout << "NeutrinoIdTool::GetNuanceCode - Error: number of true neutrinos in event must be exactly one" << std::endl;
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
    }

    return LArMCParticleHelper::GetNuanceCode(trueNeutrinos.front());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SelectAllPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &hypotheses, PfoList &selectedPfos) const
{
    for (const PfoList &pfos : hypotheses)
    {
        for (const ParticleFlowObject *const pPfo : pfos)
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["NuScore"] = -1.f;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
        }

        this->SelectPfos(pfos, selectedPfos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SelectPfosByProbability(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, const SliceFeaturesVector &sliceFeaturesVector, PfoList &selectedPfos) const
{
    // Calculate the probability of each slice that passes the minimum probability cut
    std::vector<UintFloatPair> sliceIndexProbabilityPairs;
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        const float nuProbability(sliceFeaturesVector.at(sliceIndex).GetNeutrinoProbability(m_supportVectorMachine));

        for (const ParticleFlowObject *const pPfo : crSliceHypotheses.at(sliceIndex))
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["NuScore"] = nuProbability;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
        }

        for (const ParticleFlowObject *const pPfo : nuSliceHypotheses.at(sliceIndex))
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["NuScore"] = nuProbability;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
        }

        if (nuProbability < m_minProbability)
        {
            this->SelectPfos(crSliceHypotheses.at(sliceIndex), selectedPfos);
            continue;
        }

        sliceIndexProbabilityPairs.push_back(UintFloatPair(sliceIndex, nuProbability));
    }

    // Sort the slices by probability
    std::sort(sliceIndexProbabilityPairs.begin(), sliceIndexProbabilityPairs.end(), [] (const UintFloatPair &a, const UintFloatPair &b)
    {
        return (a.second > b.second);
    });

    // Select the first m_maxNeutrinos as neutrinos, and the rest as cosmic
    unsigned int nNuSlices(0);
    for (const UintFloatPair &slice : sliceIndexProbabilityPairs)
    {
        if (nNuSlices < m_maxNeutrinos)
        {
            this->SelectPfos(nuSliceHypotheses.at(slice.first), selectedPfos);
            nNuSlices++;
            continue;
        }

        this->SelectPfos(crSliceHypotheses.at(slice.first), selectedPfos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SelectPfos(const PfoList &pfos, PfoList &selectedPfos) const
{
    selectedPfos.insert(selectedPfos.end(), pfos.begin(), pfos.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

// TODO think about how to make this function cleaner when features are more established
NeutrinoIdTool::SliceFeatures::SliceFeatures(const PfoList &nuPfos, const PfoList &crPfos, const NeutrinoIdTool *const pTool) :
    m_isAvailable(false),
    m_pTool(pTool)
{
    try
    {
        const ParticleFlowObject *const pNeutrino(this->GetNeutrino(nuPfos));
        const CartesianVector &nuVertex(LArPfoHelper::GetVertex(pNeutrino)->GetPosition());
        const PfoList &nuFinalStates(pNeutrino->GetDaughterPfoList());

        // Neutrino features
        CartesianVector nuWeightedDirTotal(0.f, 0.f, 0.f);
        unsigned int nuNHitsUsedTotal(0);
        unsigned int nuNHitsTotal(0);
        CartesianPointVector nuAllSpacePoints;
        for (const ParticleFlowObject *const pPfo : nuFinalStates)
        {
            CartesianPointVector spacePoints;
            this->GetSpacePoints(pPfo, spacePoints);

            nuAllSpacePoints.insert(nuAllSpacePoints.end(), spacePoints.begin(), spacePoints.end());
            nuNHitsTotal += spacePoints.size();

            if (spacePoints.size() < 5) continue;

            const CartesianVector dir(this->GetDirectionFromVertex(spacePoints, nuVertex));
            nuWeightedDirTotal += dir * static_cast<float>(spacePoints.size());
            nuNHitsUsedTotal += spacePoints.size();
        }

        if (nuNHitsUsedTotal == 0) return;
        const CartesianVector nuWeightedDir(nuWeightedDirTotal * (1.f / static_cast<float>(nuNHitsUsedTotal)));

        CartesianPointVector pointsInSphere;
        this->GetPointsInSphere(nuAllSpacePoints, nuVertex, 10, pointsInSphere);

        CartesianVector centroid(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        LArPcaHelper::EigenValues eigenValues(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        LArPcaHelper::EigenVectors eigenVectors;
        LArPcaHelper::RunPca(pointsInSphere, centroid, eigenValues, eigenVectors);


        const float nuNFinalStatePfos(static_cast<float>(nuFinalStates.size()));
        const float nuVertexY(nuVertex.GetY());
        const float nuWeightedDirZ(nuWeightedDir.GetZ());
        const float nuNSpacePointsInSphere(static_cast<float>(pointsInSphere.size()));

        if (eigenValues.GetX() <= std::numeric_limits<float>::epsilon()) return;
        const float nuEigenRatioInSphere(eigenValues.GetY() / eigenValues.GetX());

        // Cosmic-ray features
        unsigned int nCRHitsMax(0);
        unsigned int nCRHitsTotal(0);
        float crLongestTrackDirY(std::numeric_limits<float>::max());
        float crLongestTrackDeflection(-std::numeric_limits<float>::max());

        for (const ParticleFlowObject *const pPfo : crPfos)
        {
            CartesianPointVector spacePoints;
            this->GetSpacePoints(pPfo, spacePoints);

            nCRHitsTotal += spacePoints.size();

            if (spacePoints.size() < 5) continue;

            if (spacePoints.size() > nCRHitsMax)
            {
                nCRHitsMax = spacePoints.size();
                const CartesianVector upperDir(this->GetUpperDirection(spacePoints));
                const CartesianVector lowerDir(this->GetLowerDirection(spacePoints));

                crLongestTrackDirY = upperDir.GetY();
                crLongestTrackDeflection = 1.f - upperDir.GetDotProduct(lowerDir * (-1.f));
            }
        }

        if (nCRHitsMax == 0) return;
        if (nCRHitsTotal == 0) return;

        const float crFracHitsInLongestTrack = static_cast<float>(nCRHitsMax)/static_cast<float>(nCRHitsTotal);

        // Push the features to the feature vector
        m_featureVector.push_back(nuNFinalStatePfos);
        m_featureVector.push_back(nuNHitsTotal);
        m_featureVector.push_back(nuVertexY);
        m_featureVector.push_back(nuWeightedDirZ);
        m_featureVector.push_back(nuNSpacePointsInSphere);
        m_featureVector.push_back(nuEigenRatioInSphere);
        m_featureVector.push_back(crLongestTrackDirY);
        m_featureVector.push_back(crLongestTrackDeflection);
        m_featureVector.push_back(crFracHitsInLongestTrack);
        m_featureVector.push_back(nCRHitsMax);

        m_isAvailable = true;
    }
    catch (StatusCodeException &)
    {
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::SliceFeatures::IsFeatureVectorAvailable() const
{
    return m_isAvailable;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SliceFeatures::GetFeatureVector(LArMvaHelper::MvaFeatureVector &featureVector) const
{
    if (!m_isAvailable)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    featureVector.insert(featureVector.end(), m_featureVector.begin(), m_featureVector.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float NeutrinoIdTool::SliceFeatures::GetNeutrinoProbability(const SupportVectorMachine &supportVectorMachine) const
{
    // ATTN if one or more of the features can not be calculated, then default to calling the slice a cosmic ray
    if (!this->IsFeatureVectorAvailable()) return 0.f;

    LArMvaHelper::MvaFeatureVector featureVector;
    this->GetFeatureVector(featureVector);
    return LArMvaHelper::CalculateProbability(supportVectorMachine, featureVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *NeutrinoIdTool::SliceFeatures::GetNeutrino(const PfoList &nuPfos) const
{
    // ATTN we should only ever have one neutrino reconstructed per slice
    if (nuPfos.size() != 1)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    return nuPfos.front();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SliceFeatures::GetSpacePoints(const ParticleFlowObject *const pPfo, CartesianPointVector &spacePoints) const
{
    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

    if (clusters3D.size() > 1)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    if (clusters3D.empty()) return;

    CaloHitList caloHits;
    clusters3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHits);

    for (const CaloHit *const pCaloHit : caloHits)
        spacePoints.push_back(pCaloHit->GetPositionVector());
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector NeutrinoIdTool::SliceFeatures::GetDirectionFromVertex(const CartesianPointVector &spacePoints, const CartesianVector &vertex) const
{
    return this->GetDirection(spacePoints, [&] (const CartesianVector &pointA, const CartesianVector &pointB)
    {
        return ((pointA - vertex).GetMagnitude() < (pointB - vertex).GetMagnitude());
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector NeutrinoIdTool::SliceFeatures::GetUpperDirection(const CartesianPointVector &spacePoints) const
{
    return this->GetDirection(spacePoints, [&] (const CartesianVector &pointA, const CartesianVector &pointB)
    {
        return (pointA.GetY() > pointB.GetY());
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector NeutrinoIdTool::SliceFeatures::GetLowerDirection(const CartesianPointVector &spacePoints) const
{
    return this->GetDirection(spacePoints, [&] (const CartesianVector &pointA, const CartesianVector &pointB)
    {
        return (pointA.GetY() < pointB.GetY());
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector NeutrinoIdTool::SliceFeatures::GetDirection(const CartesianPointVector &spacePoints, std::function<bool(const CartesianVector &pointA, const CartesianVector &pointB)> fShouldChooseA) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pFirstLArTPC(m_pTool->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float layerPitch(pFirstLArTPC->GetWirePitchW());

    const ThreeDSlidingFitResult fit(&spacePoints, 5, layerPitch);
    const CartesianVector endMin(fit.GetGlobalMinLayerPosition());
    const CartesianVector endMax(fit.GetGlobalMaxLayerPosition());
    const CartesianVector dirMin(fit.GetGlobalMinLayerDirection());
    const CartesianVector dirMax(fit.GetGlobalMaxLayerDirection());

    const bool isMinStart(fShouldChooseA(endMin, endMax));
    const CartesianVector startPoint(isMinStart ? endMin : endMax);
    const CartesianVector endPoint(isMinStart ? endMax : endMin);
    const CartesianVector startDir(isMinStart ? dirMin : dirMax);

    const bool shouldFlip((endPoint - startPoint).GetUnitVector().GetDotProduct(startDir) < 0.f);
    return (shouldFlip ? startDir*(-1.f) : startDir);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SliceFeatures::GetPointsInSphere(const CartesianPointVector &spacePoints, const CartesianVector &vertex, const float radius, CartesianPointVector &spacePointsInSphere) const
{
    for (const CartesianVector &point : spacePoints)
    {
        if ((point - vertex).GetMagnitudeSquared() <= radius*radius)
            spacePointsInSphere.push_back(point);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrainingMode", m_useTrainingMode));

    if (m_useTrainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "TrainingOutputFileName", m_trainingOutputFile));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumPurity", m_minPurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumCompleteness", m_minCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectNuanceCode", m_selectNuanceCode));

    if (m_selectNuanceCode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "NuanceCode", m_nuance));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumNeutrinoProbability", m_minProbability));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaximumNeutrinos", m_maxNeutrinos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    if (!m_useTrainingMode)
    {
        std::string svmName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "SvmName", svmName));

        std::string svmFileName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "SvmFileName", svmFileName));

        const std::string fullSvmFileName(LArFileHelper::FindFileInPath(svmFileName, m_filePathEnvironmentVariable));
        m_supportVectorMachine.Initialize(fullSvmFileName, svmName);
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
