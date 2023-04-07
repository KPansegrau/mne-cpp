//=============================================================================================================
/**
 * @file     rtbfweights.cpp
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de>;
 *
 * @since    0.1.0
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
 * @brief     Definition of the RtBfWeights Class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "rtbfweights.h"

#include <fiff/fiff_info.h>

#include <mne/mne_forwardsolution.h>
#include <mne/mne_beamformer_weights.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QDebug>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/Core>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace RTPROCESSINGLIB;
using namespace Eigen;
using namespace MNELIB;
using namespace FIFFLIB;

//=============================================================================================================
// DEFINE MEMBER METHODS RtBfWeightsWorker
//=============================================================================================================

void RtBfWeightsWorker::doWork(const RtBfWeightsInput &inputData)
{
    if(this->thread()->isInterruptionRequested()) {
        return;
    }

    // Restrict forward solution as necessary for MEG
    //TODO: do we need this restriction???, check!
    MNEForwardSolution forwardMeg = inputData.pFwd->pick_types(true, false);

    MNEBeamformerWeights bfWeightsMeg(*inputData.pFiffInfo.data(),
                                        forwardMeg,
                                        inputData.dataCov,
                                        inputData.noiseCov,
                                        inputData.sWeightnorm
                                        );

    //TODO: for debugging, delete later
    qDebug() << "[RtBfWeightsWorker::doWork] bfWeightsMeg.weightnorm: " << bfWeightsMeg.weightNorm;
    qDebug() << "[RtBfWeightsWorker::doWork] bfWeightsMeg dim: " << bfWeightsMeg.weights.rows() << " x " << bfWeightsMeg.weights.cols();
    qDebug() << "[RtBfWeightsWorker::doWork] computed bfWeightsMeg.";

    emit resultReady(bfWeightsMeg);
}

//=============================================================================================================
// DEFINE MEMBER METHODS RtBfWeights
//=============================================================================================================

RtBfWeights::RtBfWeights(FiffInfo::SPtr &p_pFiffInfo,
                            MNEForwardSolution::SPtr &p_pFwd,
                            QString &p_sWeightnorm,
                            QObject *parent)
: QObject(parent)
, m_pFiffInfo(p_pFiffInfo)
, m_pFwd(p_pFwd)
, m_sWeightnorm(p_sWeightnorm)
{
    RtBfWeightsWorker *worker = new RtBfWeightsWorker;
    worker->moveToThread(&m_workerThread);

    connect(&m_workerThread, &QThread::finished,
            worker, &QObject::deleteLater);

    connect(this, &RtBfWeights::operate,
            worker, &RtBfWeightsWorker::doWork);

    connect(worker, &RtBfWeightsWorker::resultReady,
            this, &RtBfWeights::handleResults);

    m_workerThread.start();

    qRegisterMetaType<RtBfWeightsInput>("RtBfWeightsInput");
}

//=============================================================================================================

RtBfWeights::~RtBfWeights()
{
    stop();
}

//=============================================================================================================

void RtBfWeights::append(const FIFFLIB::FiffCov &noiseCov, const FIFFLIB::FiffCov &dataCov)
{
    RtBfWeightsInput inputData;
    inputData.noiseCov = noiseCov;
    inputData.dataCov = dataCov;
    inputData.pFiffInfo = m_pFiffInfo;
    inputData.pFwd = m_pFwd;
    inputData.sWeightnorm = m_sWeightnorm;

    emit operate(inputData);
}

//=============================================================================================================

void RtBfWeights::setFwdSolution(QSharedPointer<MNELIB::MNEForwardSolution> pFwd)
{
    m_pFwd = pFwd;
}

//=============================================================================================================

void RtBfWeights::setWeightnorm(QString weightnorm)
{
    m_sWeightnorm = weightnorm;
}

//=============================================================================================================

void RtBfWeights::handleResults(const MNELIB::MNEBeamformerWeights& bfWeights)
{
    emit bfWeightsCalculated(bfWeights);
}

//=============================================================================================================

void RtBfWeights::restart()
{
    stop();

    RtBfWeightsWorker *worker = new RtBfWeightsWorker;
    worker->moveToThread(&m_workerThread);

    connect(&m_workerThread, &QThread::finished,
            worker, &QObject::deleteLater);

    connect(this, &RtBfWeights::operate,
            worker, &RtBfWeightsWorker::doWork);

    connect(worker, &RtBfWeightsWorker::resultReady,
            this, &RtBfWeights::handleResults);

    m_workerThread.start();
}

//=============================================================================================================

void RtBfWeights::stop()
{
    m_workerThread.requestInterruption();
    m_workerThread.quit();
    m_workerThread.wait();
}
