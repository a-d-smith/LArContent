/**
 *  @file   larpandoracontent/LArMonitoring/KinkMonitoringAlgorithm.h
 * 
 *  @brief  Header file for the kink monitoring algorithm class
 * 
 *  $Log: $
 */
#ifndef LAR_KINK_MONITORING_ALGORITHM_H
#define LAR_KINK_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief KinkMonitoringAlgorithm class
 */
class KinkMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    KinkMonitoringAlgorithm();
    
    ~KinkMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_KINK_MONITORING_ALGORITHM_H
