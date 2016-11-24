/**
 *  @file   larpandoracontent/LArMonitoring/SliceValidationAlgorithm.cc
 * 
 *  @brief  Implementation of the slicevalidation algorithm class.
 * 
 *  $Log: $
 */

#include <iomanip>

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArMonitoring/SliceValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SliceValidationAlgorithm::~SliceValidationAlgorithm(){
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_outputTreeName.c_str(), m_outputFileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceValidationAlgorithm::Run()
{

    // Get the input collections 
    // ----------------------------------------------------------------------------------

    // Get the list of MCParticles
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    // Get the list of CaloHits
    const CaloHitList *pCaloHitListW = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListNameW, pCaloHitListW));
    
    const CaloHitList *pCaloHitListU = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListNameU, pCaloHitListU));
    
    const CaloHitList *pCaloHitListV = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListNameV, pCaloHitListV));

    // Get the list of Pfos
    const PfoList *pPfoList = nullptr;
    PfoList inputPfoList((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_pfoListName, pPfoList)) ? PfoList(*pPfoList) : PfoList());
 

    // Run the Algorithm
    // ----------------------------------------------------------------------------------
    
    SimpleEventInfo eventInfo;

    // Get the MCparticle associated with every hit, and check if it is neutrino or cosmic induced
    for (const CaloHit *const pCaloHit : *pCaloHitListW){
        AddHitToEventInfo(pCaloHit, eventInfo, 'W');
    }  
    for (const CaloHit *const pCaloHit : *pCaloHitListU){
        AddHitToEventInfo(pCaloHit, eventInfo, 'U');
    }  
    for (const CaloHit *const pCaloHit : *pCaloHitListV){
        AddHitToEventInfo(pCaloHit, eventInfo, 'V');
    }  

    // Count the number of MCparticles due to neutrinos and due to non-neutrinos (cosmic)
    for (const MCParticle *const pMCParticle : *pMCParticleList){
        if (LArMCParticleHelper::IsNeutrinoInduced(pMCParticle)){
            eventInfo.m_nNeutrinoInducedParticles++;

            // Hitless particle
            Uid id = pMCParticle->GetUid();
            if (eventInfo.m_nNeutrinoHitsWithParticleUid.count(id) == 0){
                eventInfo.m_nNeutrinoHitsWithParticleUid[id] = 0;   
            }
        }
        else{
            eventInfo.m_nCosmicInducedParticles++;

            // Hitless particle
            Uid id = pMCParticle->GetUid();
            if (eventInfo.m_nCosmicHitsWithParticleUid.count(id) == 0){
                eventInfo.m_nCosmicHitsWithParticleUid[id] = 0;   
            }
        }
        if (LArMCParticleHelper::IsNeutrino(pMCParticle)){
            Uid id = pMCParticle->GetUid();
            if (eventInfo.m_nNeutrinoHitsWithNeutrinoUid.count(id) == 0){
                eventInfo.m_nNeutrinoHitsWithNeutrinoUid[id] = 0;
            }        
        }
    }

    // Find the neutrino vertex
    //     This is a bit weird. It turns out that we tend to have mutiple "neutrinos" in a given event
    //     but usually they have no hits associated with them. In about half the events, one (and only
    //     one) will have a number of hits. It is this neutrino that we care about.
    //
    if (eventInfo.m_nNeutrinoInducedHits != 0){
        MCParticleVector trueNeutrinos;
        LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);
        for (const MCParticle *pMCParticle : trueNeutrinos){
            if (eventInfo.m_nNeutrinoHitsWithNeutrinoUid.at(pMCParticle->GetUid()) == eventInfo.m_nNeutrinoInducedHits){
                eventInfo.m_NeutrinoVertex.SetValues(pMCParticle->GetVertex().GetX(),
                                                     pMCParticle->GetVertex().GetY(),
                                                     pMCParticle->GetVertex().GetZ());
                std::cout <<  "(" << eventInfo.m_NeutrinoVertex.GetX();
                std::cout << ", " << eventInfo.m_NeutrinoVertex.GetY();
                std::cout << ", " << eventInfo.m_NeutrinoVertex.GetZ() << ")" << std::endl;
            }
        }
    }

    // Write the output
    // ----------------------------------------------------------------------------------
    WriteOutputToTerminal(eventInfo);   
    WriteOutputToRootFile(eventInfo);   

    m_eventNumber++;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::AddHitToEventInfo(const CaloHit *const hit, SimpleEventInfo &eventInfo, char plane){
    CartesianVector pos = hit->GetPositionVector();

    const MCParticle *bestMCParticle(nullptr);
    if (GetAssociatedParticle(*hit, bestMCParticle)){
        // Neutrino induced hit
        if (LArMCParticleHelper::IsNeutrinoInduced(bestMCParticle)){
            eventInfo.m_nNeutrinoInducedHits++;

            HitPositionPlane hitPos(CartesianVector(pos.GetX(), pos.GetY(), pos.GetZ()), plane);
            eventInfo.m_NeutrinoHitPositionList.push_back(hitPos);

            // Add the particle to the PDG map
            int PDG = bestMCParticle->GetParticleId();
            if (eventInfo.m_nNeutrinoInducedParticlesWithPDG.count(PDG) == 0){
                eventInfo.m_nNeutrinoInducedParticlesWithPDG[PDG] = 1;
            }
            else{
                eventInfo.m_nNeutrinoInducedParticlesWithPDG[PDG]++;
            }

            // Add the hit to the Uid map
            Uid id = bestMCParticle->GetUid();
            if (eventInfo.m_nNeutrinoHitsWithParticleUid.count(id) == 0){
                eventInfo.m_nNeutrinoHitsWithParticleUid[id] = 1;
            }
            else{
                eventInfo.m_nNeutrinoHitsWithParticleUid[id]++;
            }

            // Add the hit to the neutrino Uid map
            id = LArMCParticleHelper::GetParentNeutrino(bestMCParticle)->GetUid();
            if (eventInfo.m_nNeutrinoHitsWithNeutrinoUid.count(id) == 0){
                eventInfo.m_nNeutrinoHitsWithNeutrinoUid[id] = 1;
            }
            else{
                eventInfo.m_nNeutrinoHitsWithNeutrinoUid[id]++;
            }
        }
        else{
            // Cosmic induced hit
            eventInfo.m_nCosmicInducedHits++;

            HitPositionPlane hitPos(CartesianVector(pos.GetX(), pos.GetY(), pos.GetZ()), plane);
            eventInfo.m_CosmicHitPositionList.push_back(hitPos);

            // Add the hit to the Uid map
            Uid id = bestMCParticle->GetUid();
            if (eventInfo.m_nCosmicHitsWithParticleUid.count(id) == 0){
                eventInfo.m_nCosmicHitsWithParticleUid[id] = 1;
            }
            else{
                eventInfo.m_nCosmicHitsWithParticleUid[id]++;
            }
        }
    }
    else{
        // No MCParticle matches the hit!
        eventInfo.m_nBadHits++;

        HitPositionPlane hitPos(CartesianVector(pos.GetX(), pos.GetY(), pos.GetZ()), plane);
        eventInfo.m_BadHitPositionList.push_back(hitPos);
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SliceValidationAlgorithm::GetAssociatedParticle(const CaloHit &hit, const MCParticle* &mcParticle){
    // Get the MCParticleWeight map
    const MCParticleWeightMap &hitMCParticleWeightMap(hit.GetMCParticleWeightMap());

    // Find the MCParticle with the largest weight
    float bestWeight(0.f);

    for (const MCParticleWeightMap::value_type &mapEntry : hitMCParticleWeightMap) {
        if (mapEntry.second > bestWeight) {
            bestWeight = mapEntry.second;
            mcParticle = mapEntry.first;
        }
    }

    if (!mcParticle)
        return false;
    else
        return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::WriteOutputToTerminal(const SimpleEventInfo &eventInfo){
    int totalParticles = eventInfo.m_nNeutrinoInducedParticles + eventInfo.m_nCosmicInducedParticles;
    int totalHits      = eventInfo.m_nNeutrinoInducedHits + eventInfo.m_nCosmicInducedHits + eventInfo.m_nBadHits;

    std::cout << "Total number of MCParticles = " << totalParticles << std:: endl;
    std::cout << "    -> " << eventInfo.m_nNeutrinoInducedParticles << " (" << std::setprecision(4) << (100*((double) eventInfo.m_nNeutrinoInducedParticles) / ((double) totalParticles));
    std::cout << "%) from a neutrino primary" << std::endl;

    std::cout << "    -> " << eventInfo.m_nCosmicInducedParticles   << " (" << std::setprecision(4) << (100*((double) eventInfo.m_nCosmicInducedParticles)   / ((double) totalParticles));
    std::cout << "%) from a cosmic primary" << std::endl;
    

    std::cout << "Total number of Hits = " << totalHits << std:: endl;
    std::cout << "    -> " << eventInfo.m_nNeutrinoInducedHits << " (" << std::setprecision(4) << (100*((double) eventInfo.m_nNeutrinoInducedHits) / ((double) totalHits));
    std::cout << "%) from a neutrino primary" << std::endl;

    typedef PDGtoNumberOfParticlesMap::const_iterator it_type;
    for(it_type iterator = eventInfo.m_nNeutrinoInducedParticlesWithPDG.begin(); iterator != eventInfo.m_nNeutrinoInducedParticlesWithPDG.end(); iterator++) {
        std::cout << "        -> " << iterator->second << " with PDG " << iterator->first << std::endl;
    }

    std::cout << "    -> " << eventInfo.m_nCosmicInducedHits   << " (" << std::setprecision(4) << (100*((double) eventInfo.m_nCosmicInducedHits)   / ((double) totalHits));
    std::cout << "%) from a cosmic primary" << std::endl;
    std::cout << "    -> " << eventInfo.m_nBadHits             << " (" << std::setprecision(4) << (100*((double) eventInfo.m_nBadHits)             / ((double) totalHits));
    std::cout << "%) had no associated MCParticle" << std::endl;

}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::WriteOutputToRootFile(const SimpleEventInfo &eventInfo){
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "eventNumber"              , m_eventNumber                        ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "fileId"                   , m_fileId                             ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "nNeutrinoInducedParticles", eventInfo.m_nNeutrinoInducedParticles));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "nCosmicInducedParticles"  , eventInfo.m_nCosmicInducedParticles  ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "nNeutrinoInducedHits"     , eventInfo.m_nNeutrinoInducedHits     ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "nCosmicInducedHits"       , eventInfo.m_nCosmicInducedHits       ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "nBadHits"                 , eventInfo.m_nBadHits                 ));

    int nNeutrinoProtonHits  = this->GetNumberOfHitsWithPDG(2212, eventInfo.m_nNeutrinoInducedParticlesWithPDG);
    int nNeutrinoMuPlusHits  = this->GetNumberOfHitsWithPDG(-13 , eventInfo.m_nNeutrinoInducedParticlesWithPDG);
    int nNeutrinoMuMinusHits = this->GetNumberOfHitsWithPDG( 13 , eventInfo.m_nNeutrinoInducedParticlesWithPDG);
    int nNeutrinoPiPlusHits  = this->GetNumberOfHitsWithPDG( 211, eventInfo.m_nNeutrinoInducedParticlesWithPDG);
    int nNeutrinoPiMinusHits = this->GetNumberOfHitsWithPDG(-211, eventInfo.m_nNeutrinoInducedParticlesWithPDG);

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "nNeutrinoProtonHits"  , nNeutrinoProtonHits ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "nNeutrinoMuPlusHits"  , nNeutrinoMuPlusHits ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "nNeutrinoMuMinusHits" , nNeutrinoMuMinusHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "nNeutrinoPiPlusHits"  , nNeutrinoPiPlusHits ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "nNeutrinoPiMinusHits" , nNeutrinoPiMinusHits));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "nNeutrinoOtherHits"   , eventInfo.m_nNeutrinoInducedHits
                                                                                                                - nNeutrinoProtonHits
                                                                                                                - nNeutrinoMuPlusHits
                                                                                                                - nNeutrinoMuMinusHits
                                                                                                                - nNeutrinoPiPlusHits
                                                                                                                - nNeutrinoPiMinusHits));

    // Add the positions of the neutrino hits
    FloatVector xPositionsNeutrinoW, zPositionsNeutrinoW;
    FloatVector xPositionsNeutrinoU, zPositionsNeutrinoU;
    FloatVector xPositionsNeutrinoV, zPositionsNeutrinoV;

    for (const HitPositionPlane pHitPositionPlane : eventInfo.m_NeutrinoHitPositionList){
        if (pHitPositionPlane.second == 'W'){
            xPositionsNeutrinoW.push_back(pHitPositionPlane.first.GetX());
            zPositionsNeutrinoW.push_back(pHitPositionPlane.first.GetZ());
        }
        if (pHitPositionPlane.second == 'U'){
            xPositionsNeutrinoU.push_back(pHitPositionPlane.first.GetX());
            zPositionsNeutrinoU.push_back(pHitPositionPlane.first.GetZ());
        }
        if (pHitPositionPlane.second == 'V'){
            xPositionsNeutrinoV.push_back(pHitPositionPlane.first.GetX());
            zPositionsNeutrinoV.push_back(pHitPositionPlane.first.GetZ());
        }
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xPositionsNeutrinoW", &xPositionsNeutrinoW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zPositionsNeutrinoW", &zPositionsNeutrinoW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xPositionsNeutrinoU", &xPositionsNeutrinoU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zPositionsNeutrinoU", &zPositionsNeutrinoU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xPositionsNeutrinoV", &xPositionsNeutrinoV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zPositionsNeutrinoV", &zPositionsNeutrinoV));

    // Add the positions of the cosmic hits
    FloatVector xPositionsCosmicW, zPositionsCosmicW;
    FloatVector xPositionsCosmicU, zPositionsCosmicU;
    FloatVector xPositionsCosmicV, zPositionsCosmicV;

    for (const HitPositionPlane pHitPositionPlane : eventInfo.m_CosmicHitPositionList){
        if (pHitPositionPlane.second == 'W'){
            xPositionsCosmicW.push_back(pHitPositionPlane.first.GetX());
            zPositionsCosmicW.push_back(pHitPositionPlane.first.GetZ());
        }
        if (pHitPositionPlane.second == 'U'){
            xPositionsCosmicU.push_back(pHitPositionPlane.first.GetX());
            zPositionsCosmicU.push_back(pHitPositionPlane.first.GetZ());
        }
        if (pHitPositionPlane.second == 'V'){
            xPositionsCosmicV.push_back(pHitPositionPlane.first.GetX());
            zPositionsCosmicV.push_back(pHitPositionPlane.first.GetZ());
        }
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xPositionsCosmicW", &xPositionsCosmicW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zPositionsCosmicW", &zPositionsCosmicW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xPositionsCosmicU", &xPositionsCosmicU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zPositionsCosmicU", &zPositionsCosmicU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xPositionsCosmicV", &xPositionsCosmicV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zPositionsCosmicV", &zPositionsCosmicV));

    // Add the positions of the bad hits
    FloatVector xPositionsBadW, zPositionsBadW;
    FloatVector xPositionsBadU, zPositionsBadU;
    FloatVector xPositionsBadV, zPositionsBadV;

    for (const HitPositionPlane pHitPositionPlane : eventInfo.m_BadHitPositionList){
        if (pHitPositionPlane.second == 'W'){
            xPositionsBadW.push_back(pHitPositionPlane.first.GetX());
            zPositionsBadW.push_back(pHitPositionPlane.first.GetZ());
        }
        if (pHitPositionPlane.second == 'U'){
            xPositionsBadU.push_back(pHitPositionPlane.first.GetX());
            zPositionsBadU.push_back(pHitPositionPlane.first.GetZ());
        }
        if (pHitPositionPlane.second == 'V'){
            xPositionsBadV.push_back(pHitPositionPlane.first.GetX());
            zPositionsBadV.push_back(pHitPositionPlane.first.GetZ());
        }
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xPositionsBadW", &xPositionsBadW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zPositionsBadW", &zPositionsBadW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xPositionsBadU", &xPositionsBadU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zPositionsBadU", &zPositionsBadU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xPositionsBadV", &xPositionsBadV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zPositionsBadV", &zPositionsBadV));


    // Add the true neutrino vertex
    if (eventInfo.m_NeutrinoVertex.GetX() != std::numeric_limits<float>::max()){
        // A neutrino vertex exists
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xNeutrinoVertexW", LArGeometryHelper::ProjectPosition(this->GetPandora(), eventInfo.m_NeutrinoVertex, TPC_VIEW_W).GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zNeutrinoVertexW", LArGeometryHelper::ProjectPosition(this->GetPandora(), eventInfo.m_NeutrinoVertex, TPC_VIEW_W).GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xNeutrinoVertexU", LArGeometryHelper::ProjectPosition(this->GetPandora(), eventInfo.m_NeutrinoVertex, TPC_VIEW_U).GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zNeutrinoVertexU", LArGeometryHelper::ProjectPosition(this->GetPandora(), eventInfo.m_NeutrinoVertex, TPC_VIEW_U).GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xNeutrinoVertexV", LArGeometryHelper::ProjectPosition(this->GetPandora(), eventInfo.m_NeutrinoVertex, TPC_VIEW_V).GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zNeutrinoVertexV", LArGeometryHelper::ProjectPosition(this->GetPandora(), eventInfo.m_NeutrinoVertex, TPC_VIEW_V).GetZ()));
    }
    else{
        // A neutrino vertex doesn't exist -> use non sensical values
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xNeutrinoVertexW", std::numeric_limits<float>::max()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zNeutrinoVertexW", std::numeric_limits<float>::max()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xNeutrinoVertexU", std::numeric_limits<float>::max()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zNeutrinoVertexU", std::numeric_limits<float>::max()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "xNeutrinoVertexV", std::numeric_limits<float>::max()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName.c_str(), "zNeutrinoVertexV", std::numeric_limits<float>::max()));
    }

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_outputTreeName.c_str()));

}

//------------------------------------------------------------------------------------------------------------------------------------------

int SliceValidationAlgorithm::GetNumberOfHitsWithPDG(const int PDG, const PDGtoNumberOfParticlesMap &map){
    return ((map.count(PDG) == 0) ? 0 : map.at(PDG));
}

//------------------------------------------------------------------------------------------------------------------------------------------
StatusCode SliceValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameW"   , m_caloHitListNameW  ));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameU"   , m_caloHitListNameU  ));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameV"   , m_caloHitListNameV  ));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName" , m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName"        , m_pfoListName       ));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFileName"     , m_outputFileName    ));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTreeName"     , m_outputTreeName    ));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileIdentifier"     , m_fileId            ));

    return STATUS_CODE_SUCCESS;
}

}// namespace lar_content
