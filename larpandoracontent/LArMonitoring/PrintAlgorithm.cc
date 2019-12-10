/**
 *  @file   larpandoracontent/LArMonitoring/PrintAlgorithm.cc
 * 
 *  @brief  Implementation of the gap monitoring algorithm class 
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/PrintAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"

using namespace pandora;

namespace lar_content
{

StatusCode PrintAlgorithm::Run()
{
    LArFormattingHelper::PrintHeader(m_message);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PrintAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    StringVector message;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "Message", message));

    for (const auto &word : message)
        m_message += word + " ";

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
