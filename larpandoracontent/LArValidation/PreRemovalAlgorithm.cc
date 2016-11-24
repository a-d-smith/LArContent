#include <iomanip>

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArValidation/PreRemovalAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

PreRemovalAlgorithm::PreRemovalAlgorithm(){
    m_eventNumber = 0;
    m_outputPfoTreeName = "HitDataPfo"; 
    m_outputMCParticleTreeName = "HitDataMCParticle"; 
    m_outputHitTreeName = "HitDataHit";
}

PreRemovalAlgorithm::~PreRemovalAlgorithm(){
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_outputPfoTreeName.c_str(), m_outputFileName.c_str(), "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_outputHitTreeName.c_str(), m_outputFileName.c_str(), "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_outputMCParticleTreeName.c_str(), m_outputFileName.c_str(), "UPDATE"));
    std::cout << "File written " << m_outputFileName << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PreRemovalAlgorithm::Run()
{
    // Get the input collections 
    // ----------------------------------------------------------------------------------
    
    // Get the list of MCParticles
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    // Get the list of CaloHits
    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    
    // Get the list of Pfos
    const PfoList *pPfoList = nullptr;
    PfoList inputPfoList((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_pfoListName, pPfoList)) ? PfoList(*pPfoList) : PfoList());

    // Run the Algorithm
    // ----------------------------------------------------------------------------------
   
    // Get the list of all CaloHits and assign them Uids
    std::map<size_t, const CaloHit *const> uidHitMap;
    size_t maxHitId = 0;
    for (const CaloHit *const &pCaloHit : *pCaloHitList){ 
        uidHitMap.insert( std::pair<size_t, const CaloHit *const>(maxHitId, pCaloHit) );
        maxHitId++;
    }

    // Assign each PFO a unique id
    std::map<size_t, const ParticleFlowObject *const> uidPfoMap;
    size_t maxPfoId = 0;
    for (const ParticleFlowObject *const &pPfo : *pPfoList){
        uidPfoMap.insert(std::pair<size_t, const ParticleFlowObject *const>(maxPfoId, pPfo));
        maxPfoId++;
    }

    // A map between hit and Pfo uids
    std::map<size_t, size_t> hitPfoMap;

    // Get information on each Pfo
    std::map<size_t, const ParticleFlowObject *const>::iterator itPfo;
    for (itPfo = uidPfoMap.begin(); itPfo != uidPfoMap.end(); ++itPfo){
        // Uid
        size_t uid = itPfo->first;

        CaloHitList caloHitListW;
        HitType hitTypeW = TPC_VIEW_W;
        LArPfoHelper::GetCaloHits(itPfo->second, hitTypeW, caloHitListW);
       
        CaloHitList caloHitListU;
        HitType hitTypeU = TPC_VIEW_U;
        LArPfoHelper::GetCaloHits(itPfo->second, hitTypeU, caloHitListU);

        CaloHitList caloHitListV;
        HitType hitTypeV = TPC_VIEW_V;
        LArPfoHelper::GetCaloHits(itPfo->second, hitTypeV, caloHitListV);

        CaloHitList caloHitList;
        caloHitList.merge(caloHitListW);
        caloHitList.merge(caloHitListU);
        caloHitList.merge(caloHitListV);

        std::vector<int> *pfoHitList = new std::vector<int>;

        // Get a list of all of the hit uids
        for (const CaloHit *const &pCaloHit : caloHitList){
            // Find the id of this calohit
            size_t hitUid;
            bool foundHitId = false;
            std::map<size_t, const CaloHit *const>::iterator it;
            for (it = uidHitMap.begin(); it != uidHitMap.end(); ++it){
                if (it->second == pCaloHit){
                    hitUid = it->first;
                    foundHitId = true;
                    break;
                }
            }
            if (!foundHitId) {
                std::cerr << "Error: Could not find hit associated with Pfo" << std::endl;
                std::exit(1);
            }
            pfoHitList->push_back((int) hitUid);
            hitPfoMap.insert(std::pair<size_t, size_t>(hitUid, uid));
        }

        // Output this information to a tree
        // --------------------------------------------------
        // Identifier
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputPfoTreeName.c_str(), "FileId"    , m_fileId));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputPfoTreeName.c_str(), "EventId"   , m_eventNumber));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputPfoTreeName.c_str(), "UniqueId"  , (int) uid));
        // Associated Hit List
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputPfoTreeName.c_str(), "HitUidList", pfoHitList));

        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_outputPfoTreeName.c_str()));

        delete pfoHitList;
        
    }


    // A map between hit and mcparticle uids
    std::map<size_t, size_t> hitParticleMap;

    // Get information on each hit
    std::map<size_t, const CaloHit *const>::iterator it;
    for (it = uidHitMap.begin(); it != uidHitMap.end(); ++it){
        // Get the Uid
        size_t uid = it->first;

        // View
        HitType view = it->second->GetHitType();
        int viewInt;
        switch (view){
            case TPC_VIEW_W:
                viewInt = 0;
                break;
            case TPC_VIEW_U:
                viewInt = 1;
                break;
            case TPC_VIEW_V:
                viewInt = 2;
                break;
            default:
                viewInt = -1;
        }

        // Position
        float x = it->second->GetPositionVector().GetX();
        float z = it->second->GetPositionVector().GetZ();

        // Find the associated MCParticle (if exists)
        const MCParticleWeightMap &hitMCParticleWeightMap(it->second->GetMCParticleWeightMap());
        const MCParticle* mcParticle = nullptr;

        // Find the MCParticle with the largest weight
        float bestWeight(0.f);
        for (const MCParticleWeightMap::value_type &mapEntry : hitMCParticleWeightMap) {
            if (mapEntry.second > bestWeight) {
                bestWeight = mapEntry.second;
                mcParticle = mapEntry.first;
            }
        }
        
        // Get the Uid for the associated particle
        size_t mcParticleUid = -1;
        bool isNeutrinoInduced = false;
        if (mcParticle){
            // An associated MCParticle exists -> Not a ghost hit
            mcParticleUid = (size_t)(mcParticle->GetUid());

            if (LArMCParticleHelper::IsNeutrinoInduced(mcParticle)){
                isNeutrinoInduced = true;
            }
        }
        hitParticleMap.insert(std::pair<size_t, size_t>(uid, mcParticleUid));

        // Output this information to a tree
        // --------------------------------------------------
        // Identifier
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputHitTreeName.c_str(), "FileId"           , m_fileId));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputHitTreeName.c_str(), "EventId"          , m_eventNumber));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputHitTreeName.c_str(), "UniqueId"         , (int) uid));
        // Associated particle / pfo
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputHitTreeName.c_str(), "MCParticleId"     , (int) mcParticleUid));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputHitTreeName.c_str(), "IsNeutrinoInduced", (isNeutrinoInduced ? 1 : 0)));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputHitTreeName.c_str(), "PfoId"            , (hitPfoMap.count(uid)==1 ? ((int) hitPfoMap.at(uid)) : -1) ));
        // Positional info
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputHitTreeName.c_str(), "View"             , viewInt));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputHitTreeName.c_str(), "X"                , x));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputHitTreeName.c_str(), "Z"                , z));

        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_outputHitTreeName.c_str()));
    
    }


    // Get a list of all of the MCParticles 
    for (const MCParticle *const &pMCParticle : *pMCParticleList){
        // Identification
        size_t uid = (size_t) (pMCParticle->GetUid());

        // Does this particle come from a neutrino?
        bool isNeutrinoInduced = false;
        if (LArMCParticleHelper::IsNeutrinoInduced(pMCParticle)){
            isNeutrinoInduced = true;
        }
        
        // PDG code
        int pdg = pMCParticle->GetParticleId();

        /*
        // Uid of the primary particle
        size_t primaryUid;
        if (LArMCParticleHelper::IsPrimary(pMCParticle)){
            primaryUid = uid;
        }
        else{
            primaryUid = (size_t) (LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle)->GetUid());
        }
        */

        // List of associated Hits
        std::vector<int> *MCParticleHitList = new std::vector<int>;
        std::map<size_t, size_t>::iterator itMC;
        for (itMC = hitParticleMap.begin(); itMC != hitParticleMap.end(); ++itMC){
            if (itMC->second == uid){
                MCParticleHitList->push_back((int) (itMC->first));
            }
        }

        // Start position
        const CartesianVector StartVtx = pMCParticle->GetVertex();
        float startX = StartVtx.GetX();
        float startY = StartVtx.GetY();
        float startZ = StartVtx.GetZ();

        // End position
        const CartesianVector EndVtx = pMCParticle->GetEndpoint();
        float endX = EndVtx.GetX();
        float endY = EndVtx.GetY();
        float endZ = EndVtx.GetZ();

        // Output this information to a tree
        // --------------------------------------------------
        // Identifier
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "FileId"           , m_fileId));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "EventId"          , m_eventNumber));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "UniqueId"         , (int) uid));
        // Metadata
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "PdgCode"          , pdg));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "IsNeutrinoInduced", (isNeutrinoInduced ? 1 : 0)));
        // Associated Primary Uid
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "PrimaryUid"       , (int) primaryUid));
        // Associated Hit List
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "HitUidList"       , MCParticleHitList));
        // Positional information
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "StartX"           , startX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "StartY"           , startY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "StartZ"           , startZ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "EndX"             , endX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "EndY"             , endY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputMCParticleTreeName.c_str(), "EndZ"             , endZ));

        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_outputMCParticleTreeName.c_str()));

        delete MCParticleHitList;
    }

    m_eventNumber++;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
StatusCode PreRemovalAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName ));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName"   , m_caloHitListName    ));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName"       , m_pfoListName        ));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileId"            , m_fileId             ));

    m_outputFileName = "HitData_PRE_" + std::to_string(m_fileId) + ".root";

    return STATUS_CODE_SUCCESS;
}

}// namespace lar_content
