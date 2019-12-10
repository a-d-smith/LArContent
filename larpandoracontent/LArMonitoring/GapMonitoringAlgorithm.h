/**
 *  @file   larpandoracontent/LArMonitoring/GapMonitoringAlgorithm.h
 * 
 *  @brief  Header file for the gap monitoring algorithm class
 * 
 *  $Log: $
 */
#ifndef LAR_GAP_MONITORING_ALGORITHM_H
#define LAR_GAP_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief GapMonitoringAlgorithm class
 */
class GapMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    GapMonitoringAlgorithm();
    
    ~GapMonitoringAlgorithm();

private:
    void GetTotalYZExtent(float &minY, float &maxY, float &minZ, float &maxZ) const;
    void GetSampleWidths(float &ySampleWidth, float &zSampleWidth, const float samplesPerWire = 2.f) const;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_GAP_MONITORING_ALGORITHM_H
