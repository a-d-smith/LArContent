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
#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"

using namespace pandora;

namespace lar_content
{

SecondaryInteractionsAlgorithm::SecondaryInteractionsAlgorithm() // :
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

StatusCode SecondaryInteractionsAlgorithm::Run()
{
    return STATUS_CODE_SUCCESS;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

StatusCode SecondaryInteractionsAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
