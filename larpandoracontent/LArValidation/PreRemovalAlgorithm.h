#ifndef PREREMOVAL_ALGORITHM_H
#define PREREMOVAL_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArMonitoring/EventValidationAlgorithm.h"

#ifndef MONITORING
#include "PandoraMonitoringApi.h"
#endif

//#include "larpandoracontent/LArValidation/SimpleObjects.h"

namespace lar_content
{

/**
 *  @brief  PreRemovalAlgorithm class
 */
class PreRemovalAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief Default constructor
     */
    PreRemovalAlgorithm();

    /**
     *  @brief Destructor - output the data to a root file
     *
     */
    ~PreRemovalAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Member variables here
    int m_eventNumber;

    std::string m_mcParticleListName;
    std::string m_caloHitListName;
    std::string m_pfoListName;
    std::string m_outputFileName; 
    std::string m_outputPfoTreeName; 
    std::string m_outputMCParticleTreeName; 
    std::string m_outputHitTreeName; 
    int m_fileId;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *PreRemovalAlgorithm::Factory::CreateAlgorithm() const
{
    return new PreRemovalAlgorithm();
}

} // namespace lar_content

#endif // #ifndef PREREMOVAL_ALGORITHM_H
