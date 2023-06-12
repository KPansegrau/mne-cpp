//=============================================================================================================
/**
 * @file     beamformer.h
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
 * @brief     Beamformer class declaration.
 *
 */

#ifndef BEAMFORMER_H
#define BEAMFORMER_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../inverse_global.h"
#include "../IInverseAlgorithm.h"

#include <mne/mne_beamformer_weights.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSharedPointer>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
// DEFINE NAMESPACE
//=============================================================================================================

namespace INVERSELIB
{

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * Beamformer source estimation algorithm ToDo: Paper references.
 *
 * @brief Beamformer source estimation.
 */
class INVERSESHARED_EXPORT Beamformer : public IInverseAlgorithm
{

public:

    typedef QSharedPointer<Beamformer> SPtr;            /**< Shared pointer type for Beamformer. */
    typedef QSharedPointer<const Beamformer> ConstSPtr; /**< Const shared pointer type for Beamformer. */

    //=========================================================================================================
    // IInverseAlgorithm methods
    //=========================================================================================================


    /**
     * Constructs beamformer inverse algorithm
     * TODO second constructor is used (first can be deleted I guess)
     *
     * @param[in] p_beamformerWeights    The beamformer weights.
     * @param[in] p_fLambda              The regularization factor.
     * @param[in] P_sWeightnorm          Used weight normalization. ("no" | "unitnoisegain" | "arraygain" | "nai")  arraygain and no are the only debugged options.
     *
     * @return the prepared beamformer.
     */

    explicit Beamformer(const MNELIB::MNEBeamformerWeights &p_beamformerWeights, float p_fLambda, const QString p_sWeightnorm);

    //=========================================================================================================
    /**
     *
     * Destructor for beamformer algorithm
     *
     */

    virtual ~Beamformer(){}

    //=========================================================================================================
    /**
     * Calculates source estimate
     *
     * @param[in] p_fiffEvoked   Evoked data.
     * @param[in] pick_normal    If True, rather than pooling the orientations by taking the norm, only the.
     *                           radial component is kept. This is only applied when working with loose orientations.
     *                           This parameter is unused.
     *
     * @return the calculated source estimation.
     */
    virtual MNELIB::MNESourceEstimate calculateInverse(const FIFFLIB::FiffEvoked &p_fiffEvoked, bool pick_normal = false);

    //=========================================================================================================
    /**
     * Calculates source estimate.
     *
     * @param[in] data   Evoked input data.
     * @param[in] tmin          first point in time (from evoked input data)
     * @param[in] tstep         time step between two samples (from sampling frequency)
     * @param[in] pick_normal   If True, rather than pooling the orientations by taking the norm, only the.
     *                           radial component is kept. This is only applied when working with loose orientations.
     *                           This parameter is unused.
     *
     * @return the calculated source estimation.
     */

    virtual MNELIB::MNESourceEstimate calculateInverse(const Eigen::MatrixXd &data, float tmin, float tstep, bool pick_normal = false) const;

    //=========================================================================================================
    /**
     * Get the name of the inverse beamformer algorithm.
     *
     * @return the name of the inverse beamformer algorithm.
     */
    virtual const char* getName() const;

    //=========================================================================================================
    /**
     * Get the source space corresponding to this beamformer.
     *
     * @return the source space corresponding to this beamformer.
     */
    virtual const MNELIB::MNESourceSpace& getSourceSpace() const;

    //=========================================================================================================

    /**
     * TODO: edit docu,
     *
     * Perform the inverse setup: Prepares this beamformer filter weight matrix for use.
     * Parameters are obligatory because of corresponding abstract class but are not used in this virtual method
     *
     * @param[in] nave           Number of averages to use.
     * @param[in] pick_normal    If True, rather than pooling the orientations by taking the norm, only the.
     *                           radial component is kept. This is only applied when working with loose orientations.
     */
    virtual void doInverseSetup(qint32 nave, bool pick_normal);

    //=========================================================================================================
    // further member methods
    //=========================================================================================================

    /**
     * Get the prepared beamformer.
     *
     * @return the prepared beamformer.
     */
    inline MNELIB::MNEBeamformerWeights& getPreparedBeamformer();

    //=========================================================================================================
    /**
     * Get the beamformer filter weights
     *
     * @return the beamformer filter weights.
     */
    inline Eigen::MatrixXd& getBeamformerWeights();

    //=========================================================================================================
    /**
     * Get the beamformer weight normalization method
     *
     * @return the weight normalization method.
     */
    const QString getWeightnorm();

    //=========================================================================================================
    /**
     * Set weight normalization method.
     *
     * @param[in] weightnorm   Weight normalization method to use.
     */
    void setWeightnorm(QString weightnorm);

    //=========================================================================================================
    /**
     * Set regularization factor
     *
     * @param[in] lambda   The regularization factor.
     */
    void setRegularization(float lambda);

private:

    MNELIB::MNEBeamformerWeights    m_beamformerWeights;        /**< The beamformer weights. */
    bool                            m_bBeamformerSetup;         /**< Whether the beamformer weights are set up. */
    MNELIB::MNEBeamformerWeights    m_beamformerWeightsSetup;   /**< The setup beamformer weights. */
    Eigen::MatrixXd                 m_matWTSetup;               /**< The setup beamformer filter weight matrix W^T (the one that is applied to the data). */

    FIFFLIB::FiffInfo               m_dataInfo;                 /**< FiffInfo of the measurement data. */
    FIFFLIB::FiffCov                m_noiseCov;                 /**< The noise covariance matrix. */
    FIFFLIB::FiffCov                m_dataCov;                  /**< The data covariance matrix. */
    MNELIB::MNEForwardSolution      m_forward;                  /**< The forward solution. */
    QString                         m_sWeightnorm;              /**< Selected weight normalization method. */
    float                           m_fLambda;                  /**< Regularization parameter. */

};

//=============================================================================================================
// INLINE DEFINITIONS
//=============================================================================================================

inline Eigen::MatrixXd& Beamformer::getBeamformerWeights()
{
    //HINT: similar to getKernel of minimum norm
    return m_beamformerWeights.weights;
}

//=============================================================================================================

inline MNELIB::MNEBeamformerWeights & Beamformer::getPreparedBeamformer()
{
    //HINT similar to getPreparedInverseOperator
    return m_beamformerWeightsSetup;
}

} //namespace
#endif // BEAMFORMER_H

