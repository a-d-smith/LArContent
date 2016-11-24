#ifndef POSTREMOVAL_ALGORITHM_H
#define POSTREMOVAL_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArMonitoring/EventValidationAlgorithm.h"

#ifndef MONITORING
#include "PandoraMonitoringApi.h"
#endif

//#include "larpandoracontent/LArValidation/SimpleObjects.h"

namespace lar_content
{

/**
 *  @brief  PostRemovalAlgorithm class
 */
class PostRemovalAlgorithm : public pandora::Algorithm
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
    PostRemovalAlgorithm();

    /**
     *  @brief Destructor - output the data to a root file
     *
     */
    ~PostRemovalAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Member variables here
    int m_eventNumber;

    std::string m_mcParticleListName;
    std::string m_caloHitListName;
    std::string m_outputFileName; 
    std::string m_outputMCParticleTreeName; 
    std::string m_outputHitTreeName; 
    int m_fileId;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *PostRemovalAlgorithm::Factory::CreateAlgorithm() const
{
    return new PostRemovalAlgorithm();
}

} // namespace lar_content

#endif // #ifndef POSTREMOVAL_ALGORITHM_H
