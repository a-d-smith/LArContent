/**
 *  @file   larpandoracontent/LArMonitoring/KinkMonitoringAlgorithm.cc
 * 
 *  @brief  Implementation of the kink monitoring algorithm class 
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/KinkMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArObjects/LArPfoObjects.h"

using namespace pandora;

namespace lar_content
{

KinkMonitoringAlgorithm::KinkMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KinkMonitoringAlgorithm::~KinkMonitoringAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ranges", "ranges.root", "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkMonitoringAlgorithm::Run()
{
    const PfoList *pPfoList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));
   
    PfoList recoNeutrinos;
    LArPfoHelper::GetRecoNeutrinos(pPfoList, recoNeutrinos);

    if (recoNeutrinos.empty())
        return STATUS_CODE_SUCCESS;

    if (recoNeutrinos.size() != 1)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    const auto pNeutrino = recoNeutrinos.front();


    LArFormattingHelper::Table table({"PDG", "E", "P", "decayToMu", "Range"});

    for (const auto &pDaughter : pNeutrino->GetDaughterPfoList())
    {
        const MCParticle *pMCParticle = nullptr;
        try
        {
            pMCParticle = LArMCParticleHelper::GetMainMCParticle(pDaughter);
        }
        catch (const StatusCodeException &)
        {
            continue;
        }

        if (!pMCParticle)
            continue;

        const auto pdg = pMCParticle->GetParticleId();
        const auto energy = pMCParticle->GetEnergy();
        const auto momentum = pMCParticle->GetMomentum().GetMagnitude();

        bool hasMuDaughter = false;
        bool hasNuMuDaughter = false;
        for (const auto &pDaughterMCP : pMCParticle->GetDaughterList())
        {
            if (pDaughterMCP->GetParticleId() == -13)
                hasMuDaughter = true;
            
            if (pDaughterMCP->GetParticleId() == 14)
                hasNuMuDaughter = true;
        }

        const auto decaysToMu = hasMuDaughter && hasNuMuDaughter;

        LArTrackStateVector trackStateVector;
        try
        {
            const auto pVertex = LArPfoHelper::GetVertex(pDaughter);
            const auto pitch = LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W);
            LArPfoHelper::GetSlidingFitTrajectory(pDaughter, pVertex, 5, pitch, trackStateVector);
        }
        catch (const StatusCodeException &)
        {
            continue;
        }

        if (trackStateVector.size() < 2)
            continue;

        float range = 0.f;
        for (unsigned int i = 1; i < trackStateVector.size(); ++i)
            range += (trackStateVector.at(i).GetPosition() - trackStateVector.at(i - 1).GetPosition()).GetMagnitude();

        table.AddElement(pdg);
        table.AddElement(energy);
        table.AddElement(momentum);
        table.AddElement(decaysToMu);
        table.AddElement(range);
            
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ranges", "pdg", pdg));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ranges", "energy", energy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ranges", "momentum", momentum));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ranges", "decaysToMu", decaysToMu ? 1 : 0));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ranges", "range", range));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ranges"));
    }

    table.Print();
    

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkMonitoringAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
