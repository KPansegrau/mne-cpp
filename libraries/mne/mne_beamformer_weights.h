//=============================================================================================================
/**
 * @file     mne_beamformer_weights.h
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de>
 * @since    0.1.0
 * @date     January, 2023
 *
 * @section  LICENSE
 *
 * Copyright (C) 2023, Kerstin Pansegrau. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that
 * the following conditions are met:
 *     * Redistributions of source code must retain the above copyright notice, this list of conditions and the
 *       following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
 *       the following disclaimer in the documentation and/or other materials provided with the distribution.
 *     * Neither the name of MNE-CPP authors nor the names of its contributors may be used
 *       to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * @brief     MNEBeamformerWeights class declaration.
 *
 */

#ifndef MNE_BEAMFORMER_WEIGHTS_H
#define MNE_BEAMFORMER_WEIGHTS_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "mne_global.h"
#include "mne_forwardsolution.h"


#include <fiff/fiff_types.h>
#include <fiff/fiff_cov.h>
#include <fiff/fiff_info.h>
#include <fiff/fiff_proj.h>

#include <utils/mnemath.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSharedPointer>
#include <QStringList>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/Core>

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================


//=============================================================================================================
// DEFINE NAMESPACE
//=============================================================================================================

namespace MNELIB
{

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * Beamformer weights.
 *
 * @brief Beamformer weights.
 */
class MNESHARED_EXPORT MNEBeamformerWeights
{

public:
    typedef QSharedPointer<MNEBeamformerWeights> SPtr;            /**< Shared pointer type for MNEBeamformerWeights. */
    typedef QSharedPointer<const MNEBeamformerWeights> ConstSPtr; /**< Const shared pointer type for MNEBeamformerWeights. */

    //=========================================================================================================
    /**
    * Default constructor.
    */
    MNEBeamformerWeights();

    //=========================================================================================================
    /**
     * TODO: edit docu here, constructors analog to inverse operator constructors
     *
     * Constructs beamformer weights.
     *
     * @param[in] info               The measurement info to specify the channels to include. Bad channels in info['bads'] are not used.
     * @param[in] forward            Forward operator.
     * @param[in] p_noise_cov        The noise covariance matrix.
     * @param[in] loose              float in [0, 1]. Value that weights the source variances of the dipole components defining the tangent space of the cortical surfaces.
     * @param[in] depth              float in [0, 1]. Depth weighting coefficients. If None, no depth weighting is performed.
     * @param[in] fixed              Use fixed source orientations normal to the cortical mantle. If True, the loose parameter is ignored.
     * @param[in] limit_depth_chs    If True, use only grad channels in depth weighting (equivalent to MNE C code). If grad chanels aren't present, only mag channels will be used (if no mag, then eeg). If False, use all channels.
     */
    MNEBeamformerWeights(const FIFFLIB::FiffInfo &p_dataInfo,
                            const MNEForwardSolution& p_forward,
                            const FIFFLIB::FiffCov& p_dataCov,
                            const FIFFLIB::FiffCov& p_noiseCov,
                            QString p_sPowMethod = "trace",
                            bool p_bFixedOri = false,
                            bool p_bEstNoisePow = true,
                            bool p_bProjectMom = false,
                            QString p_sWeightnorm = "no",
                            qint32 p_iRegParam = 0,
                            qint32 p_iNAverages = 1);


    //=========================================================================================================
    /**
     * TODO: edit docu for whole class
     *
     * Copy constructor.
     *
     * @param[in] p_MNEBeamformerWeights
     */
    MNEBeamformerWeights(const MNEBeamformerWeights &p_MNEBeamformerWeights);

    //=========================================================================================================
    /**
     * Destroys the MNEBeamformerWeights.
     */
    ~MNEBeamformerWeights();

    //=========================================================================================================
    /**
     * TODO: edit docu, compute Moore Penrose pseudo inverse of matrix
     *
     *
     */


    Eigen::MatrixXd compute_pseudo_inverse(const Eigen::MatrixXd &p_matrix,
                                            double p_dEpsilon = std::numeric_limits<double>::epsilon()) const;


    //=========================================================================================================
    /**
     * TODO: edit docu, invert regularized data covariance matrix
     *
     *
     */


    Eigen::MatrixXd invert_data_cov_mat(const FIFFLIB::FiffCov &p_dataCov);


    //=========================================================================================================
    /**
     * TODO: edit docu
     *
     * Pick desired channels from list of channel names.
     *
     * @param[in] ch_names  - The channel name list to consult.
     * @param[in] include   - Channels to include (if empty, include all available).
     * @param[in] exclude   - Channels to exclude (if empty, do not exclude any).
     *
     */

    QStringList pick_channels(const QStringList &ch_names,
                                const QStringList &include = FIFFLIB::defaultQStringList,
                                const QStringList &exclude = FIFFLIB::defaultQStringList) const;

    //=========================================================================================================
    /**
     * TODO: edit docu
     *
     * Compare two lists of channel names and return the common good channels.
     *
     */

    QStringList compare_ch_names(const QStringList &chNamesA,
                                const QStringList &chNamesB,
                                const QStringList &badsA) const;

    //=========================================================================================================
    /**
     * TODO: edit docu
     *
     * Check that channels in measurement data, forward solution, data covariance matrix and noise covariance matrix are common.
     *
     * @param[in] info      The measurement info.
     * @param[in] forward   The forward solution.
     * @param[in] data_cov  The data covariance matrix.
     * @param[in] noise_cov The noise covariance matrix.
     *
     * @return
     */
    QStringList check_info_bf(const FIFFLIB::FiffInfo &p_info,
                                const MNEForwardSolution &p_forward,
                                const FIFFLIB::FiffCov &p_data_cov,
                                const FIFFLIB::FiffCov &p_noise_cov) const;

    //=========================================================================================================
    /**
     * TODO: edit docu, see .cpp file, check default value for RegParam
     *
     * Assembles the inverse operator.
     *
     * @return the assembled inverse operator.
     */

    MNEBeamformerWeights make_beamformer_weights(const FIFFLIB::FiffInfo &p_dataInfo,
                                                 const MNEForwardSolution &p_forward,
                                                   const FIFFLIB::FiffCov &p_dataCov,
                                                 const FIFFLIB::FiffCov &p_noiseCov,
                                                   QString p_sPowMethod = "trace",
                                                   bool p_bFixedOri = false,
                                                   bool p_bEstNoisePow = true,
                                                   bool p_bProjectMom = false,
                                                 //bool p_bKurtosis = false,
                                                   QString p_sWeightnorm = "no",
                                                   //qint32 &p_iLambda,
                                                   //qint32 &p_iKappa,
                                                   //qint32 &p_iTol
                                                    qint32 p_iRegParam = 0,
                                                 qint32 p_iNAverage = 1
                                                 );


    //=========================================================================================================
    /**
     *
     * Prepare computed beamformer weights for actually computing the inverse solution.
     * TODO: edit docu here
     *
     * @param[in] nave      Number of averages (scales the noise covariance).
     * @param[in] lambda2   The regularization factor.
     * @param[in] dSPM      Compute the noise-normalization factors for dSPM?.
     * @param[in] sLORETA   Compute the noise-normalization factors for sLORETA?.
     *
     * @return the prepared inverse operator.
     */
    MNEBeamformerWeights prepare_beamformer_weights() const;

    //=========================================================================================================


public:

    //TODO: edit docu here

    FIFFLIB::FiffInfoBase info;                     /**< Modified measurement info (contains only common good channels of forward solution, covariance matrices and measurement). */
    Eigen::MatrixXd weights;                        /**< Beamformer weights */
    FIFFLIB::FiffCov::SDPtr data_cov;               /**< Data covariance matrix used to compute the filters. */
    FIFFLIB::FiffCov::SDPtr noise_cov;              /**< Noise covariance matrix used to compute the filters. */
    QString weightNorm;                             /**< Type of applied weight normalization. */
    Eigen::MatrixXd whitener;                       /**< Whitens the data. */
    bool fixedOri;                                  /**< Whether the beamformer was computed for fixed source orientation */
    Eigen::VectorXd optOri;                         /**< If fixedOri, this field contains the estimated optimal orientation the beamformer is computed for (dependent on weightNorm) */
    FIFFLIB::fiff_int_t nsource;                    /**< Number of source points. */
    FIFFLIB::fiff_int_t nchan;                      /**< Number of channels. */
    Eigen::VectorXd sourcePowEst;                   /**< estimates of source power for each source location */
    Eigen::VectorXd noisePowEst;                    /**< Estimates of noise power projected through the filter for each source location */
    QList<FIFFLIB::FiffProj> projs;                 /**< SSP operator. */
    Eigen::MatrixXd proj;                           /**< The projector to apply to the data. */
    MNESourceSpace src;                             /**< Source Space. */
    FIFFLIB::fiff_int_t nave;                       /**< Number of averages used to regularize the solution. Set to 1 on single Epoch by default.*/

    //TODO: check why we need the ssp operator

//TODO we need the source_ori field analog to inverse operator field here


protected:

private:

};

//=============================================================================================================
// INLINE DEFINITIONS
//=============================================================================================================

} //namespace

//HINT: copied from mne_inverse_operator

#ifndef metatype_mnebeamformerweightssptr
#define metatype_mnebeamformerweightssptr
Q_DECLARE_METATYPE(QSharedPointer<MNELIB::MNEBeamformerWeights>); /**< Provides QT META type declaration of the QSharedPointer<MNELIB::MNEBeamformerWeights> type. For signal/slot usage.*/
#endif

#ifndef metatype_mnebeamformerweightss
#define metatype_mnebeamformerweightss
Q_DECLARE_METATYPE(MNELIB::MNEBeamformerWeights); /**< Provides QT META type declaration of the MNELIB::MNEBeamformerWeights type. For signal/slot usage.*/
#endif


#endif // MNE_BEAMFORMER_WEIGHTS_H

