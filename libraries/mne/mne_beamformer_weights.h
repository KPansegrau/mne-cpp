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

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSharedPointer>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/Core>

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================



//=============================================================================================================
// DEFINE NAMESPACE Cannot convert result of " Cpp.namespaces('MNEBeamformerWeights')[0]" to string.
//=============================================================================================================

namespace MNELIB
{

//=============================================================================================================
// Cannot convert result of " Cpp.namespaces('MNEBeamformerWeights')[0]" to string. FORWARD DECLARATIONS
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
    MNEBeamformerWeights(const FIFFLIB::FiffInfo &info,
                       const MNEForwardSolution& forward,
                       const FIFFLIB::FiffCov& p_noise_cov,
                       float loose = 0.2f,
                       float depth = 0.8f,
                       bool fixed = false,
                       bool limit_depth_chs = true);

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
     * TODO: edit docu
     *
     * Assembles the inverse operator.
     *
     * @return the assembled inverse operator.
     */
    static MNEBeamformerWeights make_beamformer_weights(const FIFFLIB::FiffInfo &info, MNEForwardSolution forward, const FIFFLIB::FiffCov &p_data_cov, float regularizationFactor, const FIFFLIB::FiffCov &p_noise_cov, float depth);


    //=========================================================================================================
    /**
     * TODO: edit docu
     *
     * @param[in] nave      Number of averages (scales the noise covariance).
     * @param[in] lambda2   The regularization factor.
     *
     */
    MNEBeamformerWeights prepare_beamformer_weights(qint32 nave,
                                                    float lambda2) const;


public:

    FIFFLIB::fiff_int_t nave;                       /**< Number of averages used to regularize the solution. Set to 1 on single Epoch by default.*/
    FIFFLIB::FiffCov::SDPtr noise_cov;              /**< Noise covariance matrix. */




protected:

private:

};

//=============================================================================================================
// INLINE DEFINITIONS
//=============================================================================================================

} //namespace
#endif // MNE_BEAMFORMER_WEIGHTS_H

