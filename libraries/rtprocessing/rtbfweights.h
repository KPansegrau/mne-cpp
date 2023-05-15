//=============================================================================================================
/**
 * @file     rtbfweights.h
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de>;
 *
 * @since    0.1.9
 * @date     February, 2023
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
 * @brief     RtBfWeights class declaration.
 *
 */

#ifndef RTBFWEIGHTS_RTPROCESSING_H
#define RTBFWEIGHTS_RTPROCESSING_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "rtprocessing_global.h"

#include <fiff/fiff_cov.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QThread>
#include <QSharedPointer>

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

namespace FIFFLIB {
    class FiffInfo;
}

namespace MNELIB {
    class MNEForwardSolution;
    class MNEBeamformerWeights;

}

//=============================================================================================================
// DEFINE NAMESPACE RTPROCESSINGLIB
//=============================================================================================================

namespace RTPROCESSINGLIB
{

//=============================================================================================================
// RTPROCESSINGLIB FORWARD DECLARATIONS
//=============================================================================================================

struct RtBfWeightsInput {
    QSharedPointer<FIFFLIB::FiffInfo>           pFiffInfo;
    QSharedPointer<MNELIB::MNEForwardSolution>  pFwd;
    FIFFLIB::FiffCov                            noiseCov;
    FIFFLIB::FiffCov                            dataCov;
    QString                                     sWeightnorm;
};

//=============================================================================================================
/**
 * Real-time beamformer weights worker.
 *
 * @brief Real-time beamformer weights worker.
 */
class RTPROCESINGSHARED_EXPORT RtBfWeightsWorker : public QObject
{
    Q_OBJECT

public:
    //=========================================================================================================
    /**
     * Perform actual beamformer weights creation.
     *
     * @param[in] inputData  Data to estimate the beamformer weights from.
     */
    void doWork(const RtBfWeightsInput &inputData);

signals:
    //=========================================================================================================
    /**
     * Emit this signal whenver new beamformer weights were estimated.
     *
     * @param[in] bfWeights  The final beamformer weights estimation.
     */
    void resultReady(const MNELIB::MNEBeamformerWeights& bfWeights);
};

//=============================================================================================================
/**
 * Real-time beamformer weight computation
 *
 * @brief Real-time beamformer weight computation
 */
class RTPROCESINGSHARED_EXPORT RtBfWeights : public QObject
{
    Q_OBJECT

public:
    typedef QSharedPointer<RtBfWeights> SPtr;             /**< Shared pointer type for RtBfWeights. */
    typedef QSharedPointer<const RtBfWeights> ConstSPtr;  /**< Const shared pointer type for RtBfWeights. */

    //=========================================================================================================
    /**
     * Creates the real-time beamformer weights estimation object
     *
     * @param[in] p_pFiffInfo    Fiff measurement info.
     * @param[in] p_pFwd         Forward solution.
     * @param[in] p_sWeightnorm  The weigth normalization method.
     * @param[in] parent         Parent QObject (optional).
     */
    explicit RtBfWeights(QSharedPointer<FIFFLIB::FiffInfo> &p_pFiffInfo,
                        QSharedPointer<MNELIB::MNEForwardSolution> &p_pFwd,
                         QString &p_sWeightnorm,
                        QObject *parent = 0);

    //=========================================================================================================
    /**
     * Destroys the real-time beamformer weights object.
     */
    ~RtBfWeights();

    //=========================================================================================================
    /**
     * Slot to receive incoming noise and data covariance estimations.
     *
     * @param[in] noiseCov      Noise covariance estimation.
     * @param[in] dataCov       Data covariance estimation.
     */
    void append(const FIFFLIB::FiffCov &noiseCov, const FIFFLIB::FiffCov &dataCov);

    //=========================================================================================================
    /**
     * Slot to receive incoming forward solution.
     *
     * @param[in] pFwd     Forward solution.
     */
    void setFwdSolution(QSharedPointer<MNELIB::MNEForwardSolution> pFwd);

    //=========================================================================================================
    /**
     * Slot to receive incoming weight normalization method.
     *
     * @param[in] weightnorm     The weightnorm method.
     */
    void setWeightnorm(QString weightnorm);


    //=========================================================================================================
    /**
     * Restarts the thread by interrupting its computation queue, quitting, waiting and then starting it again.
     */
    void restart();

    //=========================================================================================================
    /**
     * Stops the thread by interrupting its computation queue, quitting and waiting.
     */
    void stop();

protected:
    //=========================================================================================================
    /**
     * Handles the result
     */
    void handleResults(const MNELIB::MNEBeamformerWeights& bfWeights);

    QSharedPointer<FIFFLIB::FiffInfo>           m_pFiffInfo;        /**< The fiff measurement information. */
    QSharedPointer<MNELIB::MNEForwardSolution>  m_pFwd;             /**< The forward solution. */
    QString                                     m_sWeightnorm;      /**< The weight normalization method. */

    QThread                                     m_workerThread;     /**< The worker thread. */

signals:
    //=========================================================================================================
    /**
     * Signal which is emitted when the beamformer weights are calculated.
     *
     * @param[out] bfWeights  The beamformer weights.
     */
    void bfWeightsCalculated(const MNELIB::MNEBeamformerWeights& bfWeights);

    //=========================================================================================================
    /**
     * Emit this signal whenver the worker should create new beamformer weights.
     *
     * @param[in] inputData  The input data.
     */
    void operate(const RtBfWeightsInput &inputData);
};

//=============================================================================================================
// INLINE DEFINITIONS
//=============================================================================================================
} // NAMESPACE

#endif // RTBFWEIGHTS_RTPROCESSING_H
