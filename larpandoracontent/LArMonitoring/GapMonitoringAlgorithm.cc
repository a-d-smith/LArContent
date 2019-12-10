/**
 *  @file   larpandoracontent/LArMonitoring/GapMonitoringAlgorithm.cc
 * 
 *  @brief  Implementation of the gap monitoring algorithm class 
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/GapMonitoringAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

GapMonitoringAlgorithm::GapMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

GapMonitoringAlgorithm::~GapMonitoringAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "gaps", "gaps.root", "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GapMonitoringAlgorithm::Run()
{
    float ySampleWidth, zSampleWidth;
    this->GetSampleWidths(ySampleWidth, zSampleWidth);
    
    float minY, maxY, minZ, maxZ;
    this->GetTotalYZExtent(minY, maxY, minZ, maxZ);

    const auto nSamplesY = static_cast<int>(std::ceil((maxY - minY) / ySampleWidth));
    const auto nSamplesZ = static_cast<int>(std::ceil((maxZ - minZ) / zSampleWidth));
    const auto nSamplesTot = nSamplesY * nSamplesZ;

    unsigned int nInGapU = 0;
    unsigned int nInGapV = 0;
    unsigned int nInGapW = 0;
    unsigned int nIn0Gap = 0;
    unsigned int nIn1Gap = 0;
    unsigned int nIn2Gap = 0;
    unsigned int nIn3Gap = 0;

    for (unsigned int iY = 0; iY < nSamplesY; ++iY)
    {
        const auto Y = minY + iY * ySampleWidth;
        for (unsigned int iZ = 0; iZ < nSamplesZ; ++iZ)
        {
            const auto Z = minZ + iZ * zSampleWidth;
            const auto point = CartesianVector(0.f, Y, Z);

            const auto isInGapU = LArGeometryHelper::IsInGap3D(this->GetPandora(), point, TPC_VIEW_U);
            const auto isInGapV = LArGeometryHelper::IsInGap3D(this->GetPandora(), point, TPC_VIEW_V);
            const auto isInGapW = LArGeometryHelper::IsInGap3D(this->GetPandora(), point, TPC_VIEW_W);
            
            const auto nUGaps = (isInGapU ? 1 : 0);
            const auto nVGaps = (isInGapV ? 1 : 0);
            const auto nWGaps = (isInGapW ? 1 : 0);
            const auto nGaps = nUGaps + nVGaps + nWGaps;

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "gaps", "Y", Y));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "gaps", "Z", Z));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "gaps", "isInGapU", nUGaps));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "gaps", "isInGapV", nVGaps));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "gaps", "isInGapW", nWGaps));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "gaps", "nGaps", nGaps));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), "gaps"));

            if (isInGapU) nInGapU++;
            if (isInGapV) nInGapV++;
            if (isInGapW) nInGapW++;

            switch (nGaps)
            {
                case 0:
                    nIn0Gap++;
                    break;
                case 1:
                    nIn1Gap++;
                    break;
                case 2:
                    nIn2Gap++;
                    break;
                case 3:
                    nIn3Gap++;
                    break;
                default: break;
            }
        }
    }

    std::cout << "ySampleWidth : " << ySampleWidth << std::endl;
    std::cout << "zSampleWidth : " << zSampleWidth << std::endl;
    std::cout << "nSamplesY    : " << nSamplesY << std::endl;
    std::cout << "nSamplesZ    : " << nSamplesZ << std::endl;
    std::cout << "nSamplesTot  : " << nSamplesTot << std::endl;

    const auto norm = 1.f / static_cast<float>(nSamplesTot);
    std::cout << "nInGapU : " << nInGapU << " - " << (norm * nInGapU) << std::endl;
    std::cout << "nInGapV : " << nInGapV << " - " << (norm * nInGapV) << std::endl;
    std::cout << "nInGapW : " << nInGapW << " - " << (norm * nInGapW) << std::endl;
    std::cout << "nIn0Gap : " << nIn0Gap << " - " << (norm * nIn0Gap) << std::endl;
    std::cout << "nIn1Gap : " << nIn1Gap << " - " << (norm * nIn1Gap) << std::endl;
    std::cout << "nIn2Gap : " << nIn2Gap << " - " << (norm * nIn2Gap) << std::endl;
    std::cout << "nIn3Gap : " << nIn3Gap << " - " << (norm * nIn3Gap) << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GapMonitoringAlgorithm::GetTotalYZExtent(float &minY, float &maxY, float &minZ, float &maxZ) const
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);

    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    minY = parentMinY;
    maxY = parentMaxY;
    minZ = parentMinZ;
    maxZ = parentMaxZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GapMonitoringAlgorithm::GetSampleWidths(float &ySampleWidth, float &zSampleWidth, const float samplesPerWire) const
{
    const auto pitchU = LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U);
    const auto pitchV = LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V);
    const auto pitchW = LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W);

    const auto wireAxisU = LArGeometryHelper::GetWireAxis(this->GetPandora(), TPC_VIEW_U);
    const auto wireAxisV = LArGeometryHelper::GetWireAxis(this->GetPandora(), TPC_VIEW_V);
    const auto wireAxisW = LArGeometryHelper::GetWireAxis(this->GetPandora(), TPC_VIEW_W);
   
    const auto yDir = CartesianVector(0.f, 1.f, 0.f);
    const auto cosThetaUY = yDir.GetDotProduct(wireAxisU);
    const auto cosThetaVY = yDir.GetDotProduct(wireAxisV);
    const auto cosThetaWY = yDir.GetDotProduct(wireAxisW);

    const auto zDir = CartesianVector(0.f, 0.f, 1.f);
    const auto cosThetaUZ = zDir.GetDotProduct(wireAxisU);
    const auto cosThetaVZ = zDir.GetDotProduct(wireAxisV);
    const auto cosThetaWZ = zDir.GetDotProduct(wireAxisW);

    ySampleWidth = 1.f / (samplesPerWire * std::max(std::max(cosThetaUY/pitchU, cosThetaVY/pitchV), cosThetaWY/pitchW));
    ySampleWidth = 1.f / (samplesPerWire * std::max(std::max(cosThetaUY/pitchU, cosThetaVY/pitchV), cosThetaWY/pitchW));
    zSampleWidth = 1.f / (samplesPerWire * std::max(std::max(cosThetaUZ/pitchU, cosThetaVZ/pitchV), cosThetaWZ/pitchW));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GapMonitoringAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
