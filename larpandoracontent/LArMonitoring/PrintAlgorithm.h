/**
 *  @file   larpandoracontent/LArMonitoring/PrintAlgorithm.h
 * 
 *  @brief  Header file for the print algorithm class
 * 
 *  $Log: $
 */
#ifndef LAR_PRINT_ALGORITHM_H
#define LAR_PRINT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief PrintAlgorithm class
 */
class PrintAlgorithm : public pandora::Algorithm
{
public:
    PrintAlgorithm() {};
private:
    std::string  m_message;  ///< The message to print

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_PRINT_ALGORITHM_H
