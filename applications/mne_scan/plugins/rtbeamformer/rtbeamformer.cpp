//=============================================================================================================
/**
 * @file     rtbeamformer.cpp
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de
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
 * @brief    Definition of the RtBeamformer class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "rtbeamformer.h"

#include <fiff/fiff_info.h>

#include <mne/mne_forwardsolution.h>
#include <mne/mne_sourceestimate.h>
#include <mne/mne_epoch_data_list.h>

#include <scMeas/realtimesourceestimate.h>
#include <scMeas/realtimemultisamplearray.h>
#include <scMeas/realtimeevokedset.h>




//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtCore/QtPlugin>
#include <QtConcurrent>
#include <QDebug>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace RTBEAMFORMERPLUGIN;
using namespace SCSHAREDLIB;
using namespace SCMEASLIB;
using namespace UTILSLIB;
using namespace MNELIB;
using namespace FIFFLIB;
using namespace INVERSELIB;
using namespace Eigen;


//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

RtBeamformer::RtBeamformer()
    : m_pCircularMatrixBuffer(CircularBuffer_Matrix_double::SPtr(new CircularBuffer_Matrix_double(40)))
    , m_pCircularEvokedBuffer(CircularBuffer<FIFFLIB::FiffEvoked>::SPtr::create(40))
    , m_bRawInput(false)
    , m_bEvokedInput(false)
    , m_iNumAverages(1)
    , m_sAvrType("3")
{

}

//=============================================================================================================

RtBeamformer::~RtBeamformer()
{
    //HINT:copied from rtmne
    m_future.waitForFinished();

    if(this->isRunning()) {
        stop();
    }
}

//=============================================================================================================

QSharedPointer<AbstractPlugin> RtBeamformer::clone() const
{
    //HINT:copied from rtmne
    QSharedPointer<RtBeamformer> pRtBeamformerClone(new RtBeamformer());
    return pRtBeamformerClone;
}

//=============================================================================================================

void RtBeamformer::init()
{

    //TODO

    // Inits


    //Input
    m_pRTMSAInput = PluginInputData<RealTimeMultiSampleArray>::create(this, "RTBEAMFORMER RTMSA In", "RTBEAMFORMER real-time multi sample array input data");
    connect(m_pRTMSAInput.data(), &PluginInputConnector::notify,
            this, &RtBeamformer::updateRTMSA, Qt::DirectConnection);
    m_inputConnectors.append(m_pRTMSAInput);

    m_pRTESInput = PluginInputData<RealTimeEvokedSet>::create(this, "RTBEAMFORMER RTE In", "RTBEAMFORMER real-time evoked input data");
    connect(m_pRTESInput.data(), &PluginInputConnector::notify,
            this, &RtBeamformer::updateRTE, Qt::DirectConnection);
    m_inputConnectors.append(m_pRTESInput);


    //Output
    m_pRTSEOutput = PluginOutputData<RealTimeSourceEstimate>::create(this, "RTBEAMFORMER Out", "RTBEAMFORMER output data");
    m_outputConnectors.append(m_pRTSEOutput);
    m_pRTSEOutput->measurementData()->setName(this->getName());//Provide name to auto store widget settings


    // Set the annotation and surface data and mri-head transformation

}

//=============================================================================================================

void RtBeamformer::unload()
{
    //HINT: copied form rtcmne
    m_future.waitForFinished();
}

//=============================================================================================================

bool RtBeamformer::start()
{
     //HINT:copied from rtmne
    QThread::start();
    return true;
}

//=============================================================================================================

bool RtBeamformer::stop()
{
     //HINT:copied from rtmne, new clear for new list of data cov channel names
    requestInterruption();
    wait(500);

    m_qListNoiseCovChNames.clear(); //member is renamed to differentiate data and noise cov channel name list
    m_qListDataCovChNames.clear(); //This is new in comparison to rtmne
    m_bEvokedInput = false;
    m_bRawInput = false;
    m_bPluginControlWidgetsInit = false;

    return true;
}

//=============================================================================================================

AbstractPlugin::PluginType RtBeamformer::getType() const
{
     //HINT:copied from rtmne
    return _IAlgorithm;
}

//=============================================================================================================

QString RtBeamformer::getName() const
{
    return "Beamformer Source Localization";
}

//=============================================================================================================

QWidget* RtBeamformer::setupWidget()
{
    //TODO
    //placeholder. This is where the UI will be setup later.
    QWidget* setupWidget = new QWidget();
    return setupWidget;
}

//=============================================================================================================

void RtBeamformer::initPluginControlWidgets()
{

    //TODO: s. notes 20.01.2023 (what can be copied)
}

//=============================================================================================================

bool RtBeamformer::calcFiffInfo()
{
    //HINT:copied from rtmne, some modifications were made


    //HINT: copied
    QMutexLocker locker(&m_qMutex);

    //added && m_qListDataCovChNames.size() > 0 because we have the data covariance matrix as input too
    if(m_qListNoiseCovChNames.size() > 0 && m_qListDataCovChNames.size() > 0 && m_pFiffInfoInput && m_pFiffInfoForward){
        //HINT: debug messages copied, added one for data cov channel names
        qDebug() << "[RtBeamformer::calcFiffInfo] Fiff Infos available";
        qDebug() << "RtcBeamformer::calcFiffInfo - m_qListNoiseCovChNames" << m_qListNoiseCovChNames;
        qDebug() << "RtcBeamformer::calcFiffInfo - m_qListDataCovChNames" << m_qListDataCovChNames;
        qDebug() << "RtcBeamformer::calcFiffInfo - m_pFiffInfoForward->ch_names" << m_pFiffInfoForward->ch_names;
        qDebug() << "RtcBeamformer::calcFiffInfo - m_pFiffInfoInput->ch_names" << m_pFiffInfoInput->ch_names;

        //HINT: copied
        // Align channel names of the forward solution to the incoming averaged (currently acquired) data
        // Find out whether the forward solution depends on only MEG, EEG or both MEG and EEG channels
        QStringList forwardChannelsTypes;
        m_pFiffInfoForward->ch_names.clear();
        int counter = 0;

        //HINT: copied
        //TODO: do we need differentatiation of magnetometer and gradiometers here? (constants does not have FIFFV_MAG_CH etc)
        for(qint32 x = 0; x < m_pFiffInfoForward->chs.size(); ++x) {
            if(forwardChannelsTypes.contains("MEG") && forwardChannelsTypes.contains("EEG"))
                break;

            if(m_pFiffInfoForward->chs[x].kind == FIFFV_MEG_CH && !forwardChannelsTypes.contains("MEG"))
                forwardChannelsTypes<<"MEG";

            if(m_pFiffInfoForward->chs[x].kind == FIFFV_EEG_CH && !forwardChannelsTypes.contains("EEG"))
                forwardChannelsTypes<<"EEG";
        }

        //HINT: copied
        //If only MEG channels are used
        if(forwardChannelsTypes.contains("MEG") && !forwardChannelsTypes.contains("EEG")) {
            for(qint32 x = 0; x < m_pFiffInfoInput->chs.size(); ++x)
            {
                if(m_pFiffInfoInput->chs[x].kind == FIFFV_MEG_CH) {
                    m_pFiffInfoForward->chs[counter].ch_name = m_pFiffInfoInput->chs[x].ch_name;
                    m_pFiffInfoForward->ch_names << m_pFiffInfoInput->chs[x].ch_name;
                    counter++;
                }
            }
        }

        //HINT: copied
        //If only EEG channels are used
        if(!forwardChannelsTypes.contains("MEG") && forwardChannelsTypes.contains("EEG")) {
            for(qint32 x = 0; x < m_pFiffInfoInput->chs.size(); ++x)
            {
                if(m_pFiffInfoInput->chs[x].kind == FIFFV_EEG_CH) {
                    m_pFiffInfoForward->chs[counter].ch_name = m_pFiffInfoInput->chs[x].ch_name;
                    m_pFiffInfoForward->ch_names << m_pFiffInfoInput->chs[x].ch_name;
                    counter++;
                }
            }
        }

        //HINT: copied
        //If both MEG and EEG channels are used
        if(forwardChannelsTypes.contains("MEG") && forwardChannelsTypes.contains("EEG")) {
            qDebug()<<"RtcBeamformer::calcFiffInfo - MEG EEG fwd solution";
            for(qint32 x = 0; x < m_pFiffInfoInput->chs.size(); ++x)
            {
                if(m_pFiffInfoInput->chs[x].kind == FIFFV_MEG_CH || m_pFiffInfoInput->chs[x].kind == FIFFV_EEG_CH) {
                    m_pFiffInfoForward->chs[counter].ch_name = m_pFiffInfoInput->chs[x].ch_name;
                    m_pFiffInfoForward->ch_names << m_pFiffInfoInput->chs[x].ch_name;
                    counter++;
                }
            }
        }

        //HINT Copied, but added data covariance matrix channels
        //Pick only channels which are present in all data structures (covariance, evoked and forward)
        QStringList tmp1_pick_ch_names;
        QStringList tmp2_pick_ch_names; //new, we need this because temporary list should not be overwritten when comparing agains data cov channels
        foreach (const QString &ch, m_pFiffInfoForward->ch_names)
        {
            if(m_pFiffInfoInput->ch_names.contains(ch))
                tmp1_pick_ch_names << ch; //add all channels to pick list that are common for forward and data info
        }
        m_qListPickChannels.clear();

        foreach (const QString &ch, tmp1_pick_ch_names)
        {
            if(m_qListNoiseCovChNames.contains(ch))
                tmp2_pick_ch_names << ch; //add all channels to pick list that are additionally common to noise covariance matrix
        }

        foreach (const QString &ch, tmp1_pick_ch_names) //HINT: this part is new and checks the common channels with data cov matrix
        {
            if(m_qListDataCovChNames.contains(ch))
                m_qListPickChannels << ch; //add all channels to pick list that are additionally common to data covariance matrix
        }

        //HINT: copied until return true
        RowVectorXi sel = m_pFiffInfoInput->pick_channels(m_qListPickChannels);

        qDebug() << "RtBeamformer::calcFiffInfo - m_qListPickChannels.size()" << m_qListPickChannels.size();
        qDebug() << "RtcBeamformer::calcFiffInfo - m_qListPickChannels" << m_qListPickChannels;

        m_pFiffInfo = QSharedPointer<FiffInfo>(new FiffInfo(m_pFiffInfoInput->pick_info(sel)));

        m_pRTSEOutput->measurementData()->setFiffInfo(m_pFiffInfo);

        qDebug() << "RtcBeamformer::calcFiffInfo - m_pFiffInfo" << m_pFiffInfo->ch_names;

        return true;
      }

    return false;

}


//=============================================================================================================

void RtBeamformer::updateRTMSA(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    //HINT: copied from rtcmne

    if(m_pFwd) {
        QSharedPointer<RealTimeMultiSampleArray> pRTMSA = pMeasurement.dynamicCast<RealTimeMultiSampleArray>();

        if(pRTMSA && this->isRunning()) {
            //Fiff Information of the RTMSA
            m_qMutex.lock();
            if(!m_pFiffInfoInput) {
                m_pFiffInfoInput = pRTMSA->info();
                m_iNumAverages = 1;
                m_bRawInput = true;
            }
            m_qMutex.unlock();

            if(!m_bPluginControlWidgetsInit) {
                //TODO: following method still does nothing
                initPluginControlWidgets();
            }

            if(this->isRunning()) {
                // Check for artifacts
                QMap<QString,double> mapReject;
                mapReject.insert("eog", 150e-06);

                for(qint32 i = 0; i < pRTMSA->getMultiSampleArray().size(); ++i) {
                    bool bArtifactDetected = MNEEpochDataList::checkForArtifact(pRTMSA->getMultiSampleArray()[i],
                                                                                *m_pFiffInfoInput,
                                                                                mapReject);

                    if(!bArtifactDetected) {
                        // Please note that we do not need a copy here since this function will block until
                        // the buffer accepts new data again. Hence, the data is not deleted in the actual
                        // Measurement function after it emitted the notify signal.
                        while(!m_pCircularMatrixBuffer->push(pRTMSA->getMultiSampleArray()[i])) {
                            //Do nothing until the circular buffer is ready to accept new data again
                        }
                    } else {
                        qDebug() << "RtBeamformer::updateRTMSA - Reject data block";
                    }
                }
            }
        }
    }
}

//=============================================================================================================

void RtBeamformer::updateRTE(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    //HINT: copied from rtcmne

    if(m_pFwd) {
        if(QSharedPointer<RealTimeEvokedSet> pRTES = pMeasurement.dynamicCast<RealTimeEvokedSet>()) {
            QStringList lResponsibleTriggerTypes = pRTES->getResponsibleTriggerTypes();
            emit responsibleTriggerTypesChanged(lResponsibleTriggerTypes);


            if(!m_bPluginControlWidgetsInit) {
                //TODO: following method does nothing
                initPluginControlWidgets();
            }

            if(!this->isRunning() || !lResponsibleTriggerTypes.contains(m_sAvrType)) {
                return;
            }

            FiffEvokedSet::SPtr pFiffEvokedSet = pRTES->getValue();

            //Fiff Information of the evoked
            if(!m_pFiffInfoInput && pFiffEvokedSet->evoked.size() > 0) {
                QMutexLocker locker(&m_qMutex);

                for(int i = 0; i < pFiffEvokedSet->evoked.size(); ++i) {
                    if(pFiffEvokedSet->evoked.at(i).comment == m_sAvrType) {
                        m_pFiffInfoInput = QSharedPointer<FiffInfo>(new FiffInfo(pFiffEvokedSet->evoked.at(i).info));
                        break;
                    }
                }

                m_bEvokedInput = true;
            }

            if(!m_bPluginControlWidgetsInit) {
                initPluginControlWidgets();
            }

            if(this->isRunning()) {
                for(int i = 0; i < pFiffEvokedSet->evoked.size(); ++i) {
                    if(pFiffEvokedSet->evoked.at(i).comment == m_sAvrType) {
                        // Store current evoked as member so we can dispatch it if the time pick by the user changed
                        m_currentEvoked = pFiffEvokedSet->evoked.at(i).pick_channels(m_qListPickChannels);

                        // Please note that we do not need a copy here since this function will block until
                        // the buffer accepts new data again. Hence, the data is not deleted in the actual
                        // Measurement function after it emitted the notify signal.
                        while(!m_pCircularEvokedBuffer->push(pFiffEvokedSet->evoked.at(i).pick_channels(m_qListPickChannels))) {
                            //Do nothing until the circular buffer is ready to accept new data again
                        }

                            //qDebug()<<"RtcMne::updateRTE - average found type" << m_sAvrType;
                            break;
                        }
                    }
                }
            }
    }
}

//=============================================================================================================


void RtBeamformer::run()
{
     //TODO

    // Wait for fiff info to arrive
    //HINT: copied from rtcmne
    while(!calcFiffInfo()) {
        msleep(200);
    }

    // Init parameters
    //HINT: copied from rtcmne (modification: pBeamformer instead of pMinimumNorm)
    qint32 skip_count = 0;
    FiffEvoked evoked;
    MatrixXd matData;
    MatrixXd matDataResized;
    qint32 j;
    int iTimePointSps = 0;
    int iNumberChannels = 0;
    int iDownSample = 1;
    float tstep;
    float lambda2 = 1.0f / pow(1.0f, 2); //ToDo estimate lambda using covariance
    MNESourceEstimate sourceEstimate;
    bool bEvokedInput = false;
    bool bRawInput = false;
    bool bUpdateBeamformer = false;
    QSharedPointer<INVERSELIB::Beamformer> pBeamformer;
    QStringList lChNamesFiffInfo;
    QStringList lChNamesInvOp;

/*    // Start processing data
    // HINT: copied from rtcmne
    // modifications:
    while(!isInterruptionRequested()) {
        m_qMutex.lock();
        iTimePointSps = m_iTimePointSps;
        bEvokedInput = m_bEvokedInput;
        bRawInput = m_bRawInput;
        iDownSample = m_iDownSample;
        iNumberChannels = m_invOp.noise_cov->names.size();
        tstep = 1.0f / m_pFiffInfoInput->sfreq;
        lChNamesFiffInfo = m_pFiffInfoInput->ch_names;
        lChNamesInvOp = m_invOp.noise_cov->names;
        bUpdateMinimumNorm = m_bUpdateMinimumNorm;
        m_qMutex.unlock();

        if(bUpdateMinimumNorm) {
            m_qMutex.lock();
            pMinimumNorm = MinimumNorm::SPtr(new MinimumNorm(m_invOp, lambda2, m_sMethod));
            m_bUpdateMinimumNorm = false;
            m_qMutex.unlock();

            // Set up the inverse according to the parameters.
            // Use 1 nave here because in case of evoked data as input the minimum norm will always be updated when the source estimate is calculated (see run method).
            pMinimumNorm->doInverseSetup(1,true);
        }

        //Process data from raw data input
        if(bRawInput && pMinimumNorm) {
            if(((skip_count % iDownSample) == 0)) {
                // Get the current raw data
                if(m_pCircularMatrixBuffer->pop(matData)) {
                    //Pick the same channels as in the inverse operator
                    matDataResized.resize(iNumberChannels, matData.cols());

                    for(j = 0; j < iNumberChannels; ++j) {
                        matDataResized.row(j) = matData.row(lChNamesFiffInfo.indexOf(lChNamesInvOp.at(j)));
                    }

                    //TODO: Add picking here. See evoked part as input.
                    sourceEstimate = pMinimumNorm->calculateInverse(matDataResized,
                                                                    0.0f,
                                                                    tstep,
                                                                    true);

                    if(!sourceEstimate.isEmpty()) {
                        if(iTimePointSps < sourceEstimate.data.cols() && iTimePointSps >= 0) {
                            sourceEstimate = sourceEstimate.reduce(iTimePointSps,1);
                            m_pRTSEOutput->measurementData()->setValue(sourceEstimate);
                        } else {
                            m_pRTSEOutput->measurementData()->setValue(sourceEstimate);
                        }
                    }
                }
            } else {
                m_pCircularMatrixBuffer->pop(matData);
            }
        }

        //Process data from averaging input
        if(bEvokedInput && pMinimumNorm) {
            if(m_pCircularEvokedBuffer->pop(evoked)) {
                // Get the current evoked data
                if(((skip_count % iDownSample) == 0)) {
                    sourceEstimate = pMinimumNorm->calculateInverse(evoked);

                    if(!sourceEstimate.isEmpty()) {
                        if(iTimePointSps < sourceEstimate.data.cols() && iTimePointSps >= 0) {
                            sourceEstimate = sourceEstimate.reduce(iTimePointSps,1);
                            m_pRTSEOutput->measurementData()->setValue(sourceEstimate);
                        } else {
                            m_pRTSEOutput->measurementData()->setValue(sourceEstimate);
                        }
                    }
                } else {
                    m_pCircularEvokedBuffer->pop(evoked);
                }
            }
        }

        ++skip_count;
    }

*/

}

//=============================================================================================================

QString RtBeamformer::getBuildInfo()
{
    //HINT: copied from rtcmne
    return QString(RTBEAMFORMERPLUGIN::buildDateTime()) + QString(" - ")  + QString(RTBEAMFORMERPLUGIN::buildHash());

}



