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
    m_filePathEnvironmentVariable("FW_SEARCH_PATH"),
    m_useValidationMode(false),
    m_shouldPrintToScreen(true),
    m_shouldVisualizeSlices(false),
    m_shouldWriteToFile(false),
    m_treeName("neutrinoId"),
    m_fileName("NeutrinoIdValidation.root")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

NeutrinoIdTool::~NeutrinoIdTool()
{
    if (m_useValidationMode && m_shouldWriteToFile)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "NeutrinoIdTool: Unable to write tree " << m_treeName << " to file " << m_fileName << std::endl;
        }
    }
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

    // TODO factorise the below
    if (m_useValidationMode)
    {
        unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());
        std::vector<float> purities, completenesses;
        this->GetBestMCSliceIndex(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, bestSliceIndex, purities, completenesses);
    
        std::vector<float> probabilities;
        unsigned int mostProbableSliceIndex(std::numeric_limits<unsigned int>::max());
        float highestProbability(-std::numeric_limits<float>::max());

        for (unsigned int sliceIndex = 0; sliceIndex < nSlices; ++sliceIndex)
        {
            const float nuProbability(sliceFeaturesVector.at(sliceIndex).GetNeutrinoProbability(m_supportVectorMachine));
            probabilities.push_back(nuProbability);

            if (nuProbability > highestProbability)
            {
                highestProbability = nuProbability;
                mostProbableSliceIndex = sliceIndex;
            }
        }
            
        // Get the neutrino
        const MCParticleList *pMCParticleList = nullptr;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

        MCParticleVector trueNeutrinos;
        LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

        if (trueNeutrinos.size() != 1)
        {
            std::cout << "NeutrinoIdTool::GetNuanceCode - Error: number of true neutrinos in event must be exactly one" << std::endl;
            throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
        }

        const MCParticle *const pNeutrino(trueNeutrinos.front());
        const float nuPdg(pNeutrino->GetParticleId());
        const float nuEnergy(pNeutrino->GetEnergy());
        const float nuVtxX(pNeutrino->GetVertex().GetX());
        const float nuVtxY(pNeutrino->GetVertex().GetY());
        const float nuVtxZ(pNeutrino->GetVertex().GetZ());
        const float nuMomX(pNeutrino->GetMomentum().GetX());
        const float nuMomY(pNeutrino->GetMomentum().GetY());
        const float nuMomZ(pNeutrino->GetMomentum().GetZ());
    
        if (m_shouldPrintToScreen)
        {
            std::cout << "True neutrino" << std::endl;
            std::cout << "    PDG      = " << nuPdg << std::endl;
            std::cout << "    Energy   = " << nuEnergy << " GeV" << std::endl;
            std::cout << "    Vertex   = ( " << nuVtxX << ", " << nuVtxY << ", " << nuVtxZ << " ) cm" << std::endl;
            std::cout << "    Momentum = ( " << nuMomX << ", " << nuMomY << ", " << nuMomZ << " ) GeV" << std::endl;
        }

        for (unsigned int sliceIndex = 0; sliceIndex < nSlices; ++sliceIndex)
        {
            const bool isBestSlice(sliceIndex == bestSliceIndex);
            const bool isMostProbable(sliceIndex == mostProbableSliceIndex);
            const bool isCorrect((isBestSlice && isMostProbable) || (!isBestSlice && !isMostProbable));
    
            if (m_shouldPrintToScreen)
            {
                std::cout << "Slice " << sliceIndex << std::endl;
                std::cout << "    Purity        - " << purities.at(sliceIndex) << std::endl;
                std::cout << "    Completeness  - " << completenesses.at(sliceIndex) << std::endl;
                std::cout << "    Probability   - " << probabilities.at(sliceIndex) << std::endl;
                std::cout << "    Best slice    - " << isBestSlice << std::endl;
                std::cout << "    Most probable - " << isMostProbable << std::endl;
                std::cout << "    Correct ID    - " << isCorrect << std::endl;
            }

            if (m_shouldVisualizeSlices)
            {
                PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &nuSliceHypotheses.at(sliceIndex), "Slice " + std::to_string(sliceIndex) + " (Nu)", AUTOID, true, true));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &crSliceHypotheses.at(sliceIndex), "Slice " + std::to_string(sliceIndex) + " (CR)", AUTOID, true,true));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            }

            if (m_shouldWriteToFile)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuPdg", nuPdg));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuEnergy", nuEnergy));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuVtxX", nuVtxX));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuVtxY", nuVtxY));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuVtxZ", nuVtxZ));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuMomX", nuMomX));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuMomY", nuMomY));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuMomZ", nuMomZ));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "purity", purities.at(sliceIndex)));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "completeness", completenesses.at(sliceIndex)));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "score", probabilities.at(sliceIndex)));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isBestNuSlice", static_cast<int>(isBestSlice)));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isMostProbable", static_cast<int>(isMostProbable)));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
            }
        }
    }

    this->SelectPfosByProbability(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, sliceFeaturesVector, selectedPfos);
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
    std::vector<float> purities, completenesses;
    return this->GetBestMCSliceIndex(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, bestSliceIndex, purities, completenesses);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::GetBestMCSliceIndex(const Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, unsigned int &bestSliceIndex, std::vector<float> &purities, std::vector<float> &completenesses) const
{
    if (!purities.empty() || !completenesses.empty())
    {
        std::cout << "NeutrinoIdTool::GetBestMCSliceIndex - Error: input vectors for purity and completeness are not empty" << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

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

    std::vector<unsigned int> nNuHitsInSlice;
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
        const unsigned int nHits(reconstructedCaloHitList.size());
        
        const float purity(nHits > 0 ? static_cast<float>(nNuHits) / static_cast<float>(nHits) : 0.f);
        purities.push_back(purity);
        nNuHitsInSlice.push_back(nNuHits);

        if (nNuHits > nNuHitsInBestSlice)
        {
            nNuHitsInBestSlice = nNuHits;
            nHitsInBestSlice = nHits;
            bestSliceIndex = sliceIndex;
        }
    }

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)                                
    {
        const float completeness(nuNHitsTotal > 0 ? static_cast<float>(nNuHitsInSlice.at(sliceIndex)) / static_cast<float>(nuNHitsTotal) : 0.f);
        completenesses.push_back(completeness);
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
        // ATTN. Here we have to get the parent address of the hits since the PFOs were made in a daughter pandora instance
        const CaloHit *const pParentCaloHit(static_cast<const CaloHit *>(pCaloHit->GetParentAddress()));

        if (!reconstructableCaloHitSet.count(pParentCaloHit))
            continue;

        // Ensure no hits have been double counted
        if (std::find(reconstructedCaloHitList.begin(), reconstructedCaloHitList.end(), pParentCaloHit) == reconstructedCaloHitList.end())
            reconstructedCaloHitList.push_back(pParentCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int NeutrinoIdTool::CountNeutrinoInducedHits(const CaloHitList &caloHitList) const
{
    unsigned int nNuHits(0);
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
            if (LArMCParticleHelper::IsNeutrino(LArMCParticleHelper::GetParentMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit))))
                nNuHits++;
        }
        catch (const StatusCodeException &)
        {
        }
    }

    return nNuHits;
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
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseValidationMode", m_useValidationMode));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintToScreen", m_shouldPrintToScreen));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualizeSlices", m_shouldVisualizeSlices));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToFile", m_shouldWriteToFile));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputFile", m_fileName));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputTree", m_treeName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
