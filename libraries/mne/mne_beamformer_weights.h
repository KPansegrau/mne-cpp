//=============================================================================================================
/**
 * @file     mne_beamformer_weights.h
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de>
 * @since    0.1.9
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
 * MNEBeamformerWeights class.
 *
 * @brief MNEBeamformerWeights class.
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
     *
     * Constructs MNEBeamformerWeights object.
     *
     * @param[in] p_dataInfo        The measurement info to specify the channels to include. Bad channels in info['bads'] are not used.
     * @param[in] p_forward         Forward operator.
     * @param[in] p_dataCov         The data covariance matrix
     * @param[in] p_noiseCov        The noise covariance matrix.
     * @param[in] p_sWeightnorm     The weight normalization method.
     * @param[in] p_bFixedOri       Whether source estimate is computed for fixed dipole orientation or not (free orientation results into vector beamformer). fixed orientation is not debugged/provided yet.
     * @param[in] p_iRegParam       Regularization parameter
     * @param[in] p_iNAverages      Number of averages.
     *
     */

    MNEBeamformerWeights(const FIFFLIB::FiffInfo &p_dataInfo,
                            const MNEForwardSolution& p_forward,
                            const FIFFLIB::FiffCov& p_dataCov,
                            const FIFFLIB::FiffCov& p_noiseCov,
                            QString p_sWeightnorm = "no",
                            bool p_bFixedOri = false,
                            qint32 p_iRegParam = 0,
                            qint32 p_iNAverages = 1);


    //=========================================================================================================
    /**
     *
     * Copy constructor.
     *
     * @param[in] p_MNEBeamformerWeights
     */
    MNEBeamformerWeights(const MNEBeamformerWeights &p_MNEBeamformerWeights);

    //=========================================================================================================
    /**
     * Destroys MNEBeamformerWeights.
     */
    ~MNEBeamformerWeights();

    //=========================================================================================================
    /**
     * Estimates the rank of a covariance matrix.
     * If different channel types are combined, the total rank is the sum of the channel type dependent ranks.
     * (similar to procedure in MNE Python, but performed after regularization (in Covariance Evoked plug-in) here)
     *
     * @param[in] p_infoCov         The info corresponding to the covariance matrix (for separating the different channel types).
     * @param[in] p_exclude         The channels that have been excluded form the covariance matrix before (for excluding them for rank estimation too).
     * @param[in] p_covMat          The covariance matrix for which the rank should be estimated.
     *
     * @return the rank of the covariance matrix.
     *
     *
     */

    qint32 estimateCovarinaceRank(const FIFFLIB::FiffInfo &p_infoCov, const QStringList p_exclude, FIFFLIB::FiffCov p_covMat);

    //=========================================================================================================
    /**
     * Reduces prepared leadfield matrix rank to 2 for each source position
     * This is recommended for MEG input data, which has low sensitivity to radial sources and therefore one orientation often carries a lot of noise.
     */

    //TODO: delete this if rank reduction is not performed

    Eigen::MatrixXd reduceLocLeadfieldRank(const Eigen::MatrixXd &matLocLeadfield);

    //=========================================================================================================
    /**
     * Checks whether cannels of beamformer match with channels of data covariance, noise covariance and input data
     *
     */

    bool check_ch_names(const FIFFLIB::FiffInfo &info) const;

    //=========================================================================================================
    /**
     * Computes the Moore Penrose pseudo inverse by truncated SVD
     * Truncation point is derived from sum of singular values bigger than a tolerance.
     * Tolerance is computed from p_dEpsilon.
     * This method is similar to function pinv in FieldTrips ft_inverse_lcmv
     *
     * @param[in] p_matrix          The matrix to be inverted.
     * @param[in] p_dEpsilon        Machine precision for double numbers used for computation of truncation point.
     *
     * @return the pseudo inverse.
     *
     *
     */

    Eigen::MatrixXd compute_pseudo_inverse(const Eigen::MatrixXd &p_matrix,
                                            double p_dEpsilon = std::numeric_limits<double>::epsilon()) const;


    //=========================================================================================================
    /**
     * Inverts a data covariance matrix.
     * This method computes the pseudo inverse via truncated SVD.
     * The truncation point is derived from the rank of the data covariance matrix.
     *
     * @param[in] p_dataCov     The data covariance matrix to be inverted.
     *
     * @return the pseudo inverse of p_dataCov.
     *
     *
     */

    Eigen::MatrixXd invert_data_cov_mat(const FIFFLIB::FiffCov &p_dataCov);


    //=========================================================================================================
    /**
     *
     * Pick desired channels from list of channel names.
     *
     * @param[in] ch_names  - The channel name list to consult.
     * @param[in] include   - Channels to include (if empty, include all available).
     * @param[in] exclude   - Channels to exclude (if empty, do not exclude any).
     *
     * @return modified channel list with only specified channels
     *
     */

    QStringList pick_channels(const QStringList &ch_names,
                                const QStringList &include = FIFFLIB::defaultQStringList,
                                const QStringList &exclude = FIFFLIB::defaultQStringList) const;

    //=========================================================================================================
    /**
     * Compares two lists of channel names and returns the common good channels.
     *
     * @param[in] chNamesA      First list of channel names.
     * @param[in] chNamesB      Second list of channel names.
     * @param[in] badsA         List of bad channels (from first list).
     *
     * @return list of common good channels.
     *
     */

    QStringList compare_ch_names(const QStringList &chNamesA,
                                const QStringList &chNamesB,
                                const QStringList &badsA) const;

    //=========================================================================================================
    /**
     * Check that channels in measurement data, forward solution, data covariance matrix and noise covariance matrix are good and common.
     *
     * @param[in] p_info      The measurement info.
     * @param[in] p_forward   The forward solution.
     * @param[in] p_data_cov  The data covariance matrix.
     * @param[in] p_noise_cov The noise covariance matrix.
     *
     * @return list of common good channels.
     */
    QStringList check_info_bf(const FIFFLIB::FiffInfo &p_info,
                                const MNEForwardSolution &p_forward,
                                const FIFFLIB::FiffCov &p_data_cov,
                                const FIFFLIB::FiffCov &p_noise_cov) const;


    //=========================================================================================================
    /**
     *
     * Computes the beamformer filter weights. This has Fieldtrips code in ft_inverse_lcmv.m as model.
     *
     * @param[in] p_dataInfo        The measurement info to specify the channels to include. Bad channels in info['bads'] are not used.
     * @param[in] p_forward         Forward operator.
     * @param[in] p_dataCov         The data covariance matrix
     * @param[in] p_noiseCov        The noise covariance matrix.
     * @param[in] p_sWeightnorm     The weight normalization method.
     * @param[in] p_bFixedOri       Whether source estimate is computed for fixed dipole orientation or not (free orientation results into vector beamformer). fixed orientation is not debugged/provided yet.
     * @param[in] p_iRegParam       Regularization parameter
     * @param[in] p_iNAverages      Number of averages.
     *
     * @return the MNEBeamformerWeights object holding the computed beamformer filter weights.
     */

    MNEBeamformerWeights make_beamformer_weights(const FIFFLIB::FiffInfo &p_dataInfo,
                                                 const MNEForwardSolution &p_forward,
                                                   const FIFFLIB::FiffCov &p_dataCov,
                                                 const FIFFLIB::FiffCov &p_noiseCov,
                                                   bool p_bFixedOri = false,
                                                   QString p_sWeightnorm = "no",
                                                    qint32 p_iRegParam = 0,
                                                 qint32 p_iNAverage = 1
                                                 );

    //=========================================================================================================


public:


    FIFFLIB::FiffInfoBase               info;                   /**< Modified measurement info (contains only common good channels of forward solution, covariance matrices and measurement). */
    Eigen::MatrixXd                     weights;                /**< Beamformer weights */
    FIFFLIB::FiffCov::SDPtr             data_cov;               /**< Data covariance matrix used to compute the filters. */
    FIFFLIB::FiffCov::SDPtr             noise_cov;              /**< Noise covariance matrix used to compute the filters. */
    QString                             weightNorm;             /**< Type of applied weight normalization (beamformer type). */
    Eigen::MatrixXd                     whitener;               /**< Prewhitening matrix. */
    bool                                fixedOri;               /**< Whether the beamformer was computed for fixed source orientation */
    Eigen::MatrixXd                     optOri;                 /**< If fixedOri, this field contains the estimated optimal orientation the beamformer is computed for (dependent on weightNorm) */
    FIFFLIB::fiff_int_t                 nsource;                /**< Number of source points. */
    FIFFLIB::fiff_int_t                 nchan;                  /**< Number of channels. */
    QList<FIFFLIB::FiffProj>            projs;                  /**< SSP operator. This is not needed for beamforming but provided for consistency with MNEInverseOperator.*/
    Eigen::MatrixXd                     proj;                   /**< The projector to apply to the data. This is not needed for beamforming but provided for consistency with MNEInverseOperator. */
    MNESourceSpace                      src;                    /**< Source Space. */
    FIFFLIB::fiff_int_t                 nave;                   /**< Number of averages used to regularize the solution. Set to 1 on single Epoch by default.*/


private:

    //=========================================================================================================
    /**
     *
     * Computes the Frobenius norm for each column of a given matrix.
     *
     * @param[in] matrix    Matrix form which the Frobenius norm is to be calculated columnwise.
     *
     * @return the vector of the Frobenius norms of the columns of matrix.
     */

    Eigen::RowVectorXd compute_frobenius_norm(const Eigen::MatrixXd matrix);

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

