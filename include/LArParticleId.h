/**
 *  @file   LArParticleId.h
 * 
 *  @brief  Header file for the lar particle id class.
 * 
 *  $Log: $
 */
#ifndef LAR_PARTICLE_ID_H
#define LAR_PARTICLE_ID_H 1

#include "Pandora/PandoraInternal.h"
#include "Pandora/PandoraInputTypes.h"
#include "Pandora/StatusCodes.h"

#include "Xml/tinyxml.h"

namespace lar
{

/**
 *  @brief  LArParticleId class
 */
class LArParticleId
{
public:
    /**
     *  @brief  TwoDSlidingXZFitResult class
     */
    class TwoDSlidingXZFitResult
    {
    public:
        /**
         *  @brief  class LayerFitResult
         */
        class LayerFitResult
        {
        public:
            /**
             *  @brief  Constructor
             * 
             *  @param  z the z coordinate
             *  @param  fitX the fitted x coordinate
             *  @param  gradient the fitted gradient dx/dz
             *  @param  rms the rms of the fit residuals
             */
            LayerFitResult(const double z, const double fitX, const double gradient, const double rms);

            /**
             *  @brief  Get the z coordinate
             * 
             *  @return the z coordinate
             */
            double GetZ() const;

            /**
             *  @brief  Get the fitted x coordinate
             * 
             *  @return the fitted x coordinate
             */
            double GetFitX() const;

            /**
             *  @brief  Get the fitted gradient dx/dz
             * 
             *  @return the fitted gradient dx/dz
             */
            double GetGradient() const;

            /**
             *  @brief  Get the rms of the fit residuals
             * 
             *  @return the rms of the fit residuals
             */
            double GetRms() const;

            private:
            double                  m_z;                            ///< The z coordinate
            double                  m_fitX;                         ///< The fitted x coordinate
            double                  m_gradient;                     ///< The fitted gradient dx/dz
            double                  m_rms;                          ///< The rms of the fit residuals
        };

        typedef std::map<unsigned int, LayerFitResult> LayerFitResultMap;

        /**
         *  @brief  LayerFitContribution class
         */
        class LayerFitContribution
        {
        public:
            /**
             *  @brief  Constructor
             * 
             *  @param  pCaloHitList address of the calo hit list
             */
            LayerFitContribution(const pandora::CaloHitList *const pCaloHitList);

            /**
             *  @brief  Get the sum x
             * 
             *  @return the sum x
             */
            double GetSumX() const;

            /**
             *  @brief  Get the sum z
             * 
             *  @return the sum z
             */
            double GetSumZ() const;

            /**
             *  @brief  Get the sum x * x
             * 
             *  @return the sum x * x
             */
            double GetSumXX() const;

            /**
             *  @brief  Get the sum z * x
             * 
             *  @return the sum z * x
             */
            double GetSumZX() const;

            /**
             *  @brief  Get the sum z * z
             * 
             *  @return the sum z * z
             */
            double GetSumZZ() const;

            /**
             *  @brief  Get the number of points used
             * 
             *  @return the number of points used
             */
            unsigned int GetNPoints() const;

        private:
            double                  m_sumX;                ///< The sum x
            double                  m_sumZ;                ///< The sum z
            double                  m_sumXX;               ///< The sum x * x
            double                  m_sumZX;               ///< The sum z * x
            double                  m_sumZZ;               ///< The sum z * z
            unsigned int            m_nPoints;             ///< The number of points used
        };

        typedef std::map<unsigned int, LayerFitContribution> LayerFitContributionMap;

        /**
         *  @brief  Get the address of the cluster
         * 
         *  @return the address of the cluster
         */
        const pandora::Cluster *GetCluster() const;

        /**
         *  @brief  Get the layer fit half window
         * 
         *  @return the layer fit half window
         */
        unsigned int GetLayerFitHalfWindow() const;

        /**
         *  @brief  Find the largest scatter in the cluster, if above a threshold value
         * 
         *  @param  largestScatterLayer to receive the layer corresponding to the largest scatter
         * 
         *  @return STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND
         */
        pandora::StatusCode FindLargestScatter(unsigned int &largestScatterLayer) const;

        /**
         *  @brief  Get the sliding fit width
         * 
         *  @return the sliding fit width
         */
        float GetSlidingFitWidth() const;

        /**
         *  @brief  Get the layer fit result map
         * 
         *  @return the layer fit result map
         */
        const LayerFitResultMap &GetLayerFitResultMap() const;

        /**
         *  @brief  Get the layer fit contribution map
         * 
         *  @return the layer fit contribution map
         */
        const LayerFitContributionMap &GetLayerFitContributionMap() const;

    private:
        const pandora::Cluster     *m_pCluster;                 ///< The address of the cluster
        unsigned int                m_layerFitHalfWindow;       ///< The layer fit half window
        LayerFitResultMap           m_layerFitResultMap;        ///< The layer fit result map
        LayerFitContributionMap     m_layerFitContributionMap;  ///< The layer fit contribution map

        friend class LArParticleId;
    };

    typedef TwoDSlidingXZFitResult::LayerFitResultMap LayerFitResultMap;
    typedef TwoDSlidingXZFitResult::LayerFitContributionMap LayerFitContributionMap;

    /**
     *  @brief  Whether a cluster is a candidate electromagnetic shower
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    static bool LArEmShowerId(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Photon identification for use with fine granularity particle flow detectors
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    static bool LArPhotonId(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Electron identification for use with fine granularity particle flow detectors
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    static bool LArElectronId(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Muon identification for use with fine granularity particle flow detectors
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    static bool LArMuonId(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Fraction of occupied layers in cluster
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return float
     */
    static float LArLayerOccupancy(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Fraction of occupied layers in a pair of clusters
     * 
     *  @param  pCluster1 address of the first cluster
     *  @param  pCluster2 address of the second cluster
     * 
     *  @return float
     */
    static float LArLayerOccupancy(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);

    /**
     *  @brief  Measure width of cluster using multiple straight line fits
     * 
     *  @param  pCluster address of the cluster
     * 
     *  @return float
     */
    static float LArTrackWidth(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Perform two dimensional sliding x-z fit
     * 
     *  @param  pCluster address of the cluster
     *  @param  twoDSlidingXZFitResult to receive the fit result
     */
    static void LArTwoDSlidingXZFit(const pandora::Cluster *const pCluster, TwoDSlidingXZFitResult &twoDSlidingXZFitResult);

    /**
     *  @brief  Read the lar particle id settings
     * 
     *  @param  xmlHandle the relevant xml handle
     */
    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    static unsigned int             m_layerFitHalfWindow;       ///< The layer fit half window for sliding 2d x-z fits
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *LArParticleId::TwoDSlidingXZFitResult::GetCluster() const
{
    return m_pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArParticleId::TwoDSlidingXZFitResult::GetLayerFitHalfWindow() const
{
    return m_layerFitHalfWindow;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArParticleId::TwoDSlidingXZFitResult::LayerFitResultMap &LArParticleId::TwoDSlidingXZFitResult::GetLayerFitResultMap() const
{
    return m_layerFitResultMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArParticleId::TwoDSlidingXZFitResult::LayerFitContributionMap &LArParticleId::TwoDSlidingXZFitResult::GetLayerFitContributionMap() const
{
    return m_layerFitContributionMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArParticleId::TwoDSlidingXZFitResult::LayerFitResult::GetZ() const
{
    return m_z;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArParticleId::TwoDSlidingXZFitResult::LayerFitResult::GetFitX() const
{
    return m_fitX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArParticleId::TwoDSlidingXZFitResult::LayerFitResult::GetGradient() const
{
    return m_gradient;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArParticleId::TwoDSlidingXZFitResult::LayerFitResult::GetRms() const
{
    return m_rms;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArParticleId::TwoDSlidingXZFitResult::LayerFitContribution::GetSumX() const
{
    return m_sumX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArParticleId::TwoDSlidingXZFitResult::LayerFitContribution::GetSumZ() const
{
    return m_sumZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArParticleId::TwoDSlidingXZFitResult::LayerFitContribution::GetSumZX() const
{
    return m_sumZX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArParticleId::TwoDSlidingXZFitResult::LayerFitContribution::GetSumZZ() const
{
    return m_sumZZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArParticleId::TwoDSlidingXZFitResult::LayerFitContribution::GetSumXX() const
{
    return m_sumXX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArParticleId::TwoDSlidingXZFitResult::LayerFitContribution::GetNPoints() const
{
    return m_nPoints;
}

} // namespace lar

#endif // #ifndef LAR_PARTICLE_ID_H