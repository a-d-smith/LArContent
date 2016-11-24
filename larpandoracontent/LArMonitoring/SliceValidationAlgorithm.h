/**
 *  @file   WorkshopContent/workshopcontent/Algorithms/SliceValidationAlgorithm.h
 * 
 *  @brief  Header file for the slicevalidation algorithm class.
 * 
 *  $Log: $
 */
#ifndef WORKSHOP_SLICEVALIDATION_ALGORITHM_H
#define WORKSHOP_SLICEVALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArMonitoring/EventValidationAlgorithm.h"

#ifndef MONITORING
#include "PandoraMonitoringApi.h"
#endif

namespace lar_content
{

/**
 *  @brief  SliceValidationAlgorithm class
 */
class SliceValidationAlgorithm : public pandora::Algorithm
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
    SliceValidationAlgorithm() : m_eventNumber(0) {}

    /**
     *  @brief Destructor - output the data to a root file
     *
     */
    ~SliceValidationAlgorithm();

    typedef std::map<int, int> PDGtoNumberOfParticlesMap;
    typedef std::map<pandora::Uid, int> UidtoNumberOfHitsMap;
    typedef std::pair<pandora::CartesianVector, char> HitPositionPlane;
    typedef std::vector<HitPositionPlane> HitPositionList;

    /**
     *  @brief A class to hold the information on a given event
     */
    class SimpleEventInfo
    {
    public:
        int m_event;
        int m_nNeutrinoInducedParticles;
        int m_nCosmicInducedParticles;
        int m_nNeutrinoInducedHits;
        int m_nCosmicInducedHits;
        int m_nBadHits;
        PDGtoNumberOfParticlesMap m_nNeutrinoInducedParticlesWithPDG;
        UidtoNumberOfHitsMap      m_nNeutrinoHitsWithParticleUid;
        UidtoNumberOfHitsMap      m_nCosmicHitsWithParticleUid;
        UidtoNumberOfHitsMap      m_nNeutrinoHitsWithNeutrinoUid;
        HitPositionList           m_NeutrinoHitPositionList;
        HitPositionList           m_CosmicHitPositionList;
        HitPositionList           m_BadHitPositionList;
        pandora::CartesianVector  m_NeutrinoVertex;

        /**
         *  @brief Default constructor
         */
        SimpleEventInfo() : m_event(0), m_nNeutrinoInducedParticles(0), m_nCosmicInducedParticles(0), m_nNeutrinoInducedHits(0), m_nCosmicInducedHits(0), m_nBadHits(0), m_NeutrinoVertex(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),std::numeric_limits<float>::max()) {}
        
        /**
         *  @brief Destructor
         */
        ~SimpleEventInfo() {}
    };

    void AddHitToEventInfo(const pandora::CaloHit *const hit, SimpleEventInfo &eventInfo, char plane);

    /**
     *  @brief Get the MCparticle associated with a given hit. Return false if the MCParticle can not be found
     *
     *  @param hit        the hit for which we want to find the associated MCParticle
     *  @param mcParticle an empty mcParticle pointer which will be set to the associated particle by the function
     */
    bool GetAssociatedParticle(const pandora::CaloHit &hit, const pandora::MCParticle* &mcParticle);

    /**
     *  @brief Output the hit distributions to the terminal
     */
    void WriteOutputToTerminal(const SimpleEventInfo &eventInfo);

    /**
     *  @brief Output the hit distributions to a root file
     */
    void WriteOutputToRootFile(const SimpleEventInfo &eventInfo);

    /**
     *  @brief Get the number of hits with a given PDG from a PDGtoNumberOfParticlesMap
     */
    int GetNumberOfHitsWithPDG(const int PDG, const PDGtoNumberOfParticlesMap &map);

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Member variables here
    std::string m_mcParticleListName; ///< Name of input MC particle list
    std::string m_caloHitListNameU;   ///< Name of input caloHit list U
    std::string m_caloHitListNameV;   ///< Name of input caloHit list V
    std::string m_caloHitListNameW;   ///< Name of input caloHit list W
    std::string m_pfoListName;        ///< Name of input pfo list
    std::string m_outputFileName;     ///< Name of the output file
    std::string m_outputTreeName;     ///< Name of the output tree
    int         m_fileId;             ///< Identifier number to keep track of the file 

    int m_eventNumber;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SliceValidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new SliceValidationAlgorithm();
}

} // namespace workshop_content

#endif // #ifndef WORKSHOP_SLICEVALIDATION_ALGORITHM_H
