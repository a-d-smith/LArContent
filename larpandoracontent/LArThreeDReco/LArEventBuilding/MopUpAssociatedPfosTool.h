/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/MopUpAssociatedPfosTool.h
 *
 *  @brief  Header file for the end associated pfos tool class.
 *
 *  $Log: $
 */
#ifndef LAR_MOP_UP_ASSOCIATED_PFOS_TOOL_H
#define LAR_MOP_UP_ASSOCIATED_PFOS_TOOL_H 1

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoHierarchyAlgorithm.h"

#include <memory>

namespace lar_content
{

/**
 *  @brief  MopUpAssociatedPfosTool class
 */
class MopUpAssociatedPfosTool : public PfoRelationTool
{
public:
    /**
     *  @brief  Default constructor
     */
    MopUpAssociatedPfosTool();

    void Run(NeutrinoHierarchyAlgorithm *const pAlgorithm, const pandora::Vertex *const pNeutrinoVertex, NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap);

private:
    
    /**
     *  @brief  The PFO match details class
     *          Contains details about the matching from a given PFO to other PFOs and the interaction vertex
     */
    class PfoMatchDetails
    {
    public:
        /**
         *  @brief  Default constructor
         *
         *  @param  pPfo the PFO from which the matches are made
         *  @param  pfoInfoMap the map of information about the PFOs, e.g. fits to their 3D hits
         *  @param  pNeutrinoVertex the neutrino vertex position
         *  @param  allPfos the vector of all PFOs to which we should try to match
         */
        PfoMatchDetails(const pandora::ParticleFlowObject *const pPfo, const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap,
            const pandora::Vertex *const pNeutrinoVertex, const pandora::PfoVector &allPfos);

        /**
         *  @brief  Get the best matched PFO that is in the input vector of assigned PFOs
         *
         *  @param  pfoInfoMap the map of information about the PFOs, e.g. fits to their 3D hits
         *  @param  assignedPfos input vector of assigned PFOs
         *  @param  pBestMatchedPfo the output best matched PFO
         *  @param  matchScore the score for this match 
         *  @param  useInner if we should use the inner vertex for child PFO for this match
         */
        void GetBestMatchedAssignedPfo(const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap, const pandora::PfoVector &assignedPfos,
            const pandora::ParticleFlowObject *&pBestMatchedPfo, float &matchScore, bool &useInner) const;

        /**
         *  @brief  Get the best matched PFO
         *
         *  @param  pfoInfoMap the map of information about the PFOs, e.g. fits to their 3D hits
         *  @param  assignedPfos input vector of assigned PFOs
         *  @param  pBestMatchedPfo the output best matched PFO
         *  @param  matchScore the score for this match 
         *  @param  useInner if we should use the inner vertex for child PFO for this match
         */
        void GetBestMatchedPfo(const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap, const pandora::PfoVector &assignedPfos,
            const pandora::ParticleFlowObject *&pBestMatchedPfo, float &matchScore, bool &useInner) const;
        
        /**
         *  @brief  Get the vertex matching score
         *
         *  @param  matchScore the score for the match to the neutrino vertex 
         *  @param  useInner if we should use the inner vertex for child PFO for the match to the neutrino vertex
         */
        void GetVertexMatchScore(float &matchScore, bool &useInner) const;

    private:

        /**
         *  @brief  Get the matching score between two positions
         *
         *  @param  daughterPosition the position of the daughter vertex in question
         *  @param  parentPosition the position of the parent vertex in question
         */
        float GetMatchScore(const pandora::CartesianVector &daughterPosition, const pandora::CartesianVector &parentPosition) const;

        /**
         *  @brief  Calculate the matching to the neutrino vertex
         *
         *  @param  pfoInfoMap the map of information about the PFOs
         *  @param  pNeutrinoVertex the neutrino vertex position
         */
        void SetVertexMatch(const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap, const pandora::Vertex *const pNeutrinoVertex);

        /**
         *  @brief  Calculate the matching to another candidate parent PFO
         *
         *  @param  pfoInfoMap the map of information about the PFOs
         *  @param  pParentPfo the candidate parent PFO
         */
        void AddMatch(const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap, const pandora::ParticleFlowObject *const pParentPfo);

        /**
         *  @brief  The info for a match to a PFO
         */
        class PfoMatchInfo
        {
        public:
            const pandora::ParticleFlowObject *m_pParentPfo;     ///< The parent PFO we are matching to
            bool                               m_isChildInner;   ///< If we are looking at the inner vertex of the child
            bool                               m_isParentInner;  ///< If we are looking at the inner vertex of the candidate parent (not applicable for matches to vertex)
            float                              m_score;          ///< The score for the match
        };

        /**
         *  @brief  The info for a vertex match
         */
        class VertexMatchInfo
        {
        public:
            float  m_innerScore;  ///< The score of the match using the inner child vertex
            float  m_outerScore;  ///< The score of the match using the outer child vertex
        };

        /**
         *  @brief  Vector of match info
         */
        typedef std::vector<PfoMatchInfo> PfoMatchInfoVector;

        const pandora::ParticleFlowObject  *m_pPfo;                ///< The child PFO from which the matches are made
        PfoMatchInfoVector                  m_pfoMatchInfoVector;  ///< The matching info to PFOs
        VertexMatchInfo                     m_vertexMatchInfo;     ///< The match info to the neutrino vertex
    };
      
    /**
     *  @brief  Mapping from PFO to it's match details
     */
    typedef std::unordered_map<const pandora::ParticleFlowObject*, std::shared_ptr<PfoMatchDetails> > PfoMatchingMap;

    /**
     *  @brief  Precalculate the matching info between PFOs for efficiency
     *
     *  @param  pfoInfoMap the map of PFO info
     *  @param  unassignedPfos the current parentless PFOs
     *  @param  assignedPfos the current PFOs with parents
     *  @param  pNeutrinoVertex the neutrino vertex position
     *  @param  pfoMatchingMap the output mapping from unmatched PFO to its matching details
     */
    void PopulatePfoMatchingMap(const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap, const pandora::PfoVector &unassignedPfos,
        const pandora::PfoVector &assignedPfos, const pandora::Vertex *const pNeutrinoVertex, PfoMatchingMap &pfoMatchingMap) const;

    /**
     *  @brief  Print the info from the matching map
     *
     *  @param  unassignedPfos the current parentless PFOs
     *  @param  assignedPfos the current PFOs with parents
     *  @param  pfoMatchingMap the output mapping from unmatched PFO to its matching details
     *  @param  pfoInfoMap the map of PFO info
     */
    void PrintPfoMatchingMap(const pandora::PfoVector &unassignedPfos, const pandora::PfoVector &assignedPfos,
        const PfoMatchingMap &pfoMatchingMap, const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap) const;

    /**
     *  @brief  Make parent / daughter links with unmatched PFOs whose best match is to the neutrino vertex or a currently assigned PFO
     *
     *  @param  pAlgorithm the algorithm calling the tool
     *  @param  unassignedPfos the current parentless PFOs
     *  @param  assignedPfos the current PFOs with parents
     *  @param  pfoMatchingMap the output mapping from unmatched PFO to its matching details
     *  @param  pfoInfoMap the map of PFO info
     *
     *  @return if there was a match made
     */
    bool MakeClearLinksToNeutrinoHierarchy(NeutrinoHierarchyAlgorithm *const pAlgorithm, const pandora::PfoVector &unassignedPfos,
        const pandora::PfoVector &assignedPfos, const PfoMatchingMap &pfoMatchingMap,
        NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap) const;

    /**
     *  @brief  Make parent / daughter links with the unmatched PFO with the next best match to the neutrino vertex or a currently assigned PFO
     *
     *  @param  pAlgorithm the algorithm calling the tool
     *  @param  unassignedPfos the current parentless PFOs
     *  @param  assignedPfos the current PFOs with parents
     *  @param  pfoMatchingMap the output mapping from unmatched PFO to its matching details
     *  @param  pfoInfoMap the map of PFO info
     */
    void MakeNextBestLinkToNeutrinoHierarchy(NeutrinoHierarchyAlgorithm *const pAlgorithm, const pandora::PfoVector &unassignedPfos,
        const pandora::PfoVector &assignedPfos, const PfoMatchingMap &pfoMatchingMap,
        NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap) const;

    /**
     *  @brief  Make a parent / daughter link for the input PFO to the neutrino vertex
     *
     *  @param  pPfo the PFO to link
     *  @param  useInner if the inner vertex of the PFO matches to the vertex (vs. the outer vertex)
     *  @param  score the matching score for the link
     *  @param  pfoInfoMap the map of PFO info
     *  @param  pAlgorithm the algorithm calling the tool
     */
    void MakeVertexLink(const pandora::ParticleFlowObject *const pPfo, const bool useInner, const float score,
        NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap, NeutrinoHierarchyAlgorithm *const pAlgorithm) const;

    /**
     *  @brief  Make a parent / daughter link for the input child PFO to another parent PFO
     *
     *  @param  pChildPfo the PFO to link
     *  @param  pParentPfo the PFO to link pChildPfo to
     *  @param  useInner if the inner vertex of the child PFO matches to the vertex (vs. the outer vertex)
     *  @param  score the matching score for the link
     *  @param  pfoInfoMap the map of PFO info
     *  @param  pAlgorithm the algorithm calling the tool
     */
    void MakePfoLink(const pandora::ParticleFlowObject *const pChildPfo, const pandora::ParticleFlowObject *const pParentPfo,
        const bool useInner, const float score, NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap,
        NeutrinoHierarchyAlgorithm *const pAlgorithm) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_shouldDisplayMatchingInfo;  ///< If we should print matching info to the terminal
};

} // namespace lar_content

#endif // #ifndef LAR_MOP_UP_ASSOCIATED_PFOS_TOOL_H
