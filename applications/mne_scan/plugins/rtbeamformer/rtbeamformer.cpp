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

#include "FormFiles/rtbeamformersetupwidget.h"

#include "disp/viewers/beamformersettingsview.h"

#include <fs/annotationset.h>
#include <fs/surfaceset.h>

#include <fiff/fiff_info.h>

#include <mne/mne_forwardsolution.h>
#include <mne/mne_sourceestimate.h>
#include <mne/mne_epoch_data_list.h>

#include <inverse/beamformer/beamformer.h>

#include "rtprocessing/rtbfweights.h"

#include <scMeas/realtimesourceestimate.h>
#include <scMeas/realtimemultisamplearray.h>
#include <scMeas/realtimeevokedset.h>
#include <scMeas/realtimeevokedcov.h>
#include <scMeas/realtimefwdsolution.h>

#include <utils/ioutils.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtCore/QtPlugin>
#include <QtConcurrent>
#include <QDebug>
#include <QtWidgets> //TODO check why we need this include here but not in rtcmne

#include <QElapsedTimer>


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
using namespace DISPLIB;
using namespace FSLIB;
using namespace RTPROCESSINGLIB;


//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

RtBeamformer::RtBeamformer()
    : m_pCircularMatrixBuffer(CircularBuffer_Matrix_double::SPtr(new CircularBuffer_Matrix_double(40)))
    , m_pCircularEvokedBuffer(CircularBuffer<FIFFLIB::FiffEvoked>::SPtr::create(40))
    , m_bRawInput(false)
    , m_bEvokedInput(false)
    , m_bUpdateBeamformer(false)
    , m_iNumAverages(1)
    , m_iTimePointSps(0)
    , m_iDownSample(1)
    , m_sWeightnorm("no")
    , m_sAvrType("3")
    , m_sAtlasDir(QCoreApplication::applicationDirPath() + "/MNE-sample-data/subjects/sample/label")
    , m_sSurfaceDir(QCoreApplication::applicationDirPath() + "/MNE-sample-data/subjects/sample/surf")
    , m_fMriHeadTrans(QCoreApplication::applicationDirPath() + "/MNE-sample-data/MEG/sample/all-trans.fif")
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

    qInfo() << "[RtBeamformer::init] Initializing RtBeamformer plugin...";

    // Inits
    //HINT: copied from RtcMne::init()
    m_pAnnotationSet = AnnotationSet::SPtr(new AnnotationSet(m_sAtlasDir+"/lh.aparc.a2009s.annot", m_sAtlasDir+"/rh.aparc.a2009s.annot"));
    m_pSurfaceSet = SurfaceSet::SPtr(new SurfaceSet(m_sSurfaceDir+"/lh.orig", m_sSurfaceDir+"/rh.orig"));
    m_mriHeadTrans = FIFFLIB::FiffCoordTrans(m_fMriHeadTrans);

    //Input
    //HINT: copied from RtcMne::init(), changes are described below

/*    //TODO: this part is for raw data that is not averaged, first draft works only with evoked (averaged) input data
    m_pRTMSAInput = PluginInputData<RealTimeMultiSampleArray>::create(this, "RtBeamformer RTMSA In", "RtBeamformer real-time multi sample array input data");
    connect(m_pRTMSAInput.data(), &PluginInputConnector::notify,
            this, &RtBeamformer::updateRTMSA, Qt::DirectConnection);
    m_inputConnectors.append(m_pRTMSAInput);
*/
    //real-time evoked input
    m_pRTESInput = PluginInputData<RealTimeEvokedSet>::create(this, "RtBeamformer RTE In", "RtBeamformerreal-time evoked input data");
    connect(m_pRTESInput.data(), &PluginInputConnector::notify,
            this, &RtBeamformer::updateRTE, Qt::DirectConnection);
    m_inputConnectors.append(m_pRTESInput);

    //real-time covariance evoked input
    //HINT: this differs from cov input for rtcmne plugin
    m_pRTCEInput = PluginInputData<RealTimeEvokedCov>::create(this,"RtBeamformer RTCE In", "RtBeamformer real-time covariance evoked input data");
    connect(m_pRTCEInput.data(), &PluginInputConnector::notify,
            this, &RtBeamformer::updateRTC, Qt::DirectConnection);
    m_inputConnectors.append(m_pRTCEInput);

    //real-time forward solution input
    m_pRTFSInput = PluginInputData<RealTimeFwdSolution>::create(this, "RtBeamformer RTFS In", "RtBeamformer real-time forward solution input data");
    connect(m_pRTFSInput.data(), &PluginInputConnector::notify,
            this, &RtBeamformer::updateRTFS, Qt::DirectConnection);
    m_inputConnectors.append(m_pRTFSInput);

    //Output
    m_pRTSEOutput = PluginOutputData<RealTimeSourceEstimate>::create(this, "RtBeamformer Out", "RtBeamformer output data");
    m_outputConnectors.append(m_pRTSEOutput);
    m_pRTSEOutput->measurementData()->setName(this->getName());//Provide name to auto store widget settings

    // Set the annotation and surface data and mri-head transformation
    //HINT: copied from RtcMne::init()
    if(m_pAnnotationSet->size() != 0) {
        m_pRTSEOutput->measurementData()->setAnnotSet(m_pAnnotationSet);
    }

    if(m_pSurfaceSet->size() != 0) {
        m_pRTSEOutput->measurementData()->setSurfSet(m_pSurfaceSet);
    }

    if(!m_mriHeadTrans.isEmpty()) {
        m_pRTSEOutput->measurementData()->setMriHeadTrans(m_mriHeadTrans);
    }

    //qDebug() << "[RtBeamformer::init] Finished initializing RtBeamformer plugin.";
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
    m_qListDataCovChNames.clear(); //HINT: This is new in comparison to rtmne
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
    RtBeamformerSetupWidget* setupWidget = new RtBeamformerSetupWidget(this);//widget is later distroyed by CentralWidget - so it has to be created everytime new

    return setupWidget;
}

//=============================================================================================================

void RtBeamformer::initPluginControlWidgets()
{
    qInfo() << "[RtBeamformer::initPluginControlWidgets] Initializing ControlWidgets...";

    //HINT: copied from rtcmne, changed names to beamformer
    QList<QWidget*> plControlWidgets;

    BeamformerSettingsView* pBeamformerSettingsView = new BeamformerSettingsView(QString("MNESCAN/%1").arg(this->getName()));
    connect(this, &RtBeamformer::guiModeChanged,
            pBeamformerSettingsView, &BeamformerSettingsView::setGuiMode);
    pBeamformerSettingsView->setObjectName("group_tab_Settings_Beamformer Source Localization");

    connect(pBeamformerSettingsView, &BeamformerSettingsView::weightnormChanged,
            this, &RtBeamformer::onWeightnormChanged);
    connect(pBeamformerSettingsView, &BeamformerSettingsView::triggerTypeChanged,
            this, &RtBeamformer::onTriggerTypeChanged);
    connect(pBeamformerSettingsView, &BeamformerSettingsView::timePointChanged,
            this, &RtBeamformer::onTimePointValueChanged);
    connect(this, &RtBeamformer::responsibleTriggerTypesChanged,
            pBeamformerSettingsView, &BeamformerSettingsView::setTriggerTypes);

    plControlWidgets.append(pBeamformerSettingsView);

    emit pluginControlWidgetsChanged(plControlWidgets, this->getName());

    m_bPluginControlWidgetsInit = true;

    //qDebug() << "[RtBeamformer::initPluginControlWidgets] Finished initializing ControlWidgets.";
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
 /*       qDebug() << "[RtBeamformer::calcFiffInfo] Fiff Infos available";
        qDebug() << "RtcBeamformer::calcFiffInfo - m_qListNoiseCovChNames" << m_qListNoiseCovChNames;
        qDebug() << "RtcBeamformer::calcFiffInfo - m_qListDataCovChNames" << m_qListDataCovChNames;
        qDebug() << "RtcBeamformer::calcFiffInfo - m_pFiffInfoForward->ch_names" << m_pFiffInfoForward->ch_names;
        qDebug() << "RtcBeamformer::calcFiffInfo - m_pFiffInfoInput->ch_names" << m_pFiffInfoInput->ch_names;
*/
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
            qInfo()<<"RtcBeamformer::calcFiffInfo - MEG fwd solution";
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
            qInfo()<<"RtcBeamformer::calcFiffInfo - EEG fwd solution";
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
            qInfo()<<"RtcBeamformer::calcFiffInfo - MEG EEG fwd solution";
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

        //qDebug() << "[RtBeamformer::calcFiffInfo] m_qListPickChannels.size()" << m_qListPickChannels.size();
        //qDebug() << "[RtBeamformer::calcFiffInfo] m_qListPickChannels" << m_qListPickChannels;

        m_pFiffInfo = QSharedPointer<FiffInfo>(new FiffInfo(m_pFiffInfoInput->pick_info(sel)));

        m_pRTSEOutput->measurementData()->setFiffInfo(m_pFiffInfo);

        qDebug() << "[RtcBeamformer::calcFiffInfo] m_pFiffInfo" << m_pFiffInfo->ch_names;


        return true;
      }

    return false;

}


//=============================================================================================================
/*
void RtBeamformer::updateRTMSA(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    //HINT: copied from rtcmne
    qDebug() << "[RtBeamformer::updateRTMSA] Updating RTMSA...";

    if(m_pFwd) {

        qDebug() << "[RtBeamformer::updateRTMSA] m_pFwd true";

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

    qDebug() << "[RtBeamformer::updateRTMSA] Finished updating RTMSA.";
}
*/
//=============================================================================================================

void RtBeamformer::updateRTE(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    //HINT: copied from rtcmne

    //qDebug() << "[RtBeamformer::updateRTE] Updating RTE....";

    if(m_pFwd) {

        //qDebug() << "[RtBeamformer::updateRTE] m_pFwd true";

        if(QSharedPointer<RealTimeEvokedSet> pRTES = pMeasurement.dynamicCast<RealTimeEvokedSet>()) {
            QStringList lResponsibleTriggerTypes = pRTES->getResponsibleTriggerTypes();
            emit responsibleTriggerTypesChanged(lResponsibleTriggerTypes);


            if(!m_bPluginControlWidgetsInit) {
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

        //qDebug() << "[RtBeamformer::updateRTE] Finished updating RTE.";
}

//=============================================================================================================

void RtBeamformer::updateRTC(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    //TODO
    //HINT: copied from rtcmne::updateRTC, changed realtimecov to realtimeevokedcov


    if(m_pFwd) {

        qInfo() << "[RtBeamformer::updateRTC] Updating RTC...";

        QSharedPointer<RealTimeEvokedCov> pRTC = pMeasurement.dynamicCast<RealTimeEvokedCov>();

        if(pRTC && this->isRunning()) {
            // Init Real-Time inverse estimator
            if(!m_pRtBfWeights && m_pFiffInfo && m_pFwd && !m_sWeightnorm.isEmpty()) {
                m_pRtBfWeights = RtBfWeights::SPtr(new RtBfWeights(m_pFiffInfo, m_pFwd, m_sWeightnorm));
                connect(m_pRtBfWeights.data(), &RtBfWeights::bfWeightsCalculated,
                        this, &RtBeamformer::updateBFWeights);
            }

            //Fiff Information of the covariance
            //HINT: modified to RealTimeEvokedCov
            //TODO: access getValue()->first.names.size() might be wrong, check this
            if(m_qListNoiseCovChNames.size() != pRTC->getValue()->first.names.size()) {
                m_qListNoiseCovChNames = pRTC->getValue()->first.names;
            }
            if(m_qListDataCovChNames.size() != pRTC->getValue()->second.names.size()) {
                m_qListDataCovChNames = pRTC->getValue()->second.names;
            }

            //HINT: added setting of data cov here
            if(this->isRunning() && m_pRtBfWeights){

                uint64_t ms_start_bf_upd_cov = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
                std::cout.precision(64);
                std::cout << "[RtBeamformer::updateRTC] sytem_clock start after new cov: " << ms_start_bf_upd_cov << " milliseconds since the Epoch\n";

                m_pNoiseCov = QSharedPointer<FiffCov>(new FiffCov(pRTC->getValue()->first)); //we only want the first of the pair for the member
                m_pDataCov = QSharedPointer<FiffCov>(new FiffCov(pRTC->getValue()->second));

                //qDebug() << "[RtBeamformer::updateRTC] Noise cov dim: " << m_pNoiseCov->data.rows() << " x " << m_pNoiseCov->data.cols();
                //qDebug() << "[RtBeamformer::updateRTC] Data cov dim: " << m_pDataCov->data.rows() << " x " << m_pDataCov->data.cols();

                m_pRtBfWeights->append(*m_pNoiseCov, *m_pDataCov);

            }
        }
    }

    //qDebug() << "[RtBeamformer::updateRTC] Finished updating RTC.";

}

//=============================================================================================================

void RtBeamformer::updateRTFS(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    //HINT: copied from RtcMne::updateRTFS

    qInfo() << "[RtBeamformer::updateRTFS] Updating RTFS...";

    if(QSharedPointer<RealTimeFwdSolution> pRTFS = pMeasurement.dynamicCast<RealTimeFwdSolution>()) {
        if(pRTFS->isClustered()) {
            m_pFwd = pRTFS->getValue();
            m_pRTSEOutput->measurementData()->setFwdSolution(m_pFwd);

            m_qMutex.lock();
            m_pFiffInfoForward = QSharedPointer<FiffInfoBase>(new FiffInfoBase(m_pFwd->info));
            m_qMutex.unlock();

            // update beamformer weights
            if(this->isRunning() && m_pRtBfWeights) {

                uint64_t ms_start_bf_upd_fwd = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
                std::cout.precision(64);
                std::cout << "[RtBeamformer::updateRTFS] sytem_clock start after new fwd solution: " << ms_start_bf_upd_fwd << " milliseconds since the Epoch\n";

                m_pRtBfWeights->setFwdSolution(m_pFwd);
                m_pRtBfWeights->setWeightnorm(m_sWeightnorm);
                m_pRtBfWeights->append(*m_pNoiseCov, *m_pDataCov);
            }
        } else if(!pRTFS->isClustered()) {
            qWarning() << "[RtBeamformer::updateRTFS] The forward solution has not been clustered yet.";
        }
    }

    //TODO: do we need to set m_bUpdateBeamformer to true after forward solution is updated?

    //qDebug() << "[RtBeamformer::updateRTFS] Finished updating RTFS.";

}

//=============================================================================================================

void RtBeamformer::updateBFWeights(const MNEBeamformerWeights& bfWeights)
{
    //HINT: copied and modified according to rtcmne::updateInvOp()
    QMutexLocker locker(&m_qMutex);

    m_bfWeights = bfWeights;

    m_bUpdateBeamformer = true;

    qDebug() << "[RtBeamformer::updateBFWeights] updated m_bfWeights.";

    uint64_t ms_stop_bf_upd = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    std::cout.precision(64);
    std::cout << "[RtBeamformer::updateBFWeights] sytem_clock stop updateBFWeights: " << ms_stop_bf_upd << " ms since the Epoch\n";

}

//=============================================================================================================

void RtBeamformer::onWeightnormChanged(const QString& weightnorm)
{
    //HINT: copied from rtcmne, changed names
    QMutexLocker locker(&m_qMutex);

    qInfo() << "[RtBeamformer::onWeightnormChanged] Changed weightnorm to " << weightnorm;

    m_sWeightnorm = weightnorm;

    //TODO: add call of beamformer weights constructor similar to update method procedure here

    if(this->isRunning()) {

        if(!m_pFiffInfo || !m_pFwd || !m_pNoiseCov || !m_pDataCov){
                qWarning() << "[RtBeamformer::onWeightnormChanged] Missing FiffInfo and/or foward solution and/or covariance matrices, no computation of beamformer weights possible.";
                return;

        }else if(!m_pRtBfWeights && m_pFiffInfo && m_pFwd && m_pNoiseCov && m_pDataCov) {
            //create m_pRtBfWeights
            m_pRtBfWeights = RtBfWeights::SPtr(new RtBfWeights(m_pFiffInfo, m_pFwd, m_sWeightnorm));
            connect(m_pRtBfWeights.data(), &RtBfWeights::bfWeightsCalculated,
                this, &RtBeamformer::updateBFWeights);

            qDebug() << "[RtBeamformer::onWeightnormChanged] Created new m_pRtBfWeights.";
        }


        if(m_pRtBfWeights && m_pFiffInfo && m_pFwd && m_pNoiseCov && m_pDataCov){


            uint64_t ms_start_bf_upd_weightnorm = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
            std::cout.precision(64);
            std::cout << "[RtBeamformer::onWeightnormChanged] sytem_clock start after weightnorm changed: " << ms_start_bf_upd_weightnorm << " milliseconds since the Epoch\n";

            //update beamformer weights
            m_pRtBfWeights->setFwdSolution(m_pFwd);
            m_pRtBfWeights->setWeightnorm(m_sWeightnorm);
            m_pRtBfWeights->append(*m_pNoiseCov, *m_pDataCov);

            qDebug() << "[RtBeamformer::onWeightnormChanged] Updated m_pRtBfWeights.";


            return;
        }
    }

}

//=============================================================================================================

void RtBeamformer::onTriggerTypeChanged(const QString& triggerType)
{
    //HINT: copied from rtcmne
    m_sAvrType = triggerType;
}

//=============================================================================================================

void RtBeamformer::onTimePointValueChanged(int iTimePointMs)
{
    //HINT: copied from rtcmne
    if(m_pFiffInfoInput && m_pCircularEvokedBuffer) {
        m_qMutex.lock();
        m_iTimePointSps = m_pFiffInfoInput->sfreq * (float)iTimePointMs * 0.001f;
        m_qMutex.unlock();

        if(this->isRunning()) {
            while(!m_pCircularEvokedBuffer->push(m_currentEvoked)) {
                //Do nothing until the circular buffer is ready to accept new data again
            }
        }
    }
}


//=============================================================================================================


void RtBeamformer::run()
{

    qInfo() << "[RtBeamformer::run] Running RtBeamformer plugin....";

    // Wait for fiff info to arrive
    //HINT: copied from rtcmne
    while(!calcFiffInfo()) {
        msleep(200);
    }



    // Init parameters
    //HINT: copied from rtcmne, modifications: names
    qint32 skip_count = 0;
    FiffEvoked evoked;
    int iTimePointSps = 0;
    //int iNumberChannels = 0;
    int iDownSample = 1;
    //float tstep;
    //HINT: in rtcmne this parameter is called lambda2
    float regParam = 1.0f / pow(1.0f, 2); //ToDo estimate lambda using covariance.
    MNESourceEstimate sourceEstimate;
    bool bEvokedInput = false;
    bool bRawInput = false;
    bool bUpdateBeamformer = false; //HINT: changed it from bUpdateMinimumNorm
    QSharedPointer<INVERSELIB::Beamformer> pBeamformer; //HINT: changed MinimumNorm to Beamformer here
    //QStringList lChNamesFiffInfo;
    //QStringList lChNamesBfWeights; //HINT: changed name form lChNamesInvOp




    qInfo() << "[RtBeamformer::run] Start processing data....";


    //Start processing data
    //HINT: copied from rtcmne in parts
    while(!isInterruptionRequested()){

        //HINT: this part is copied from rtcmne, with only name modifications
        m_qMutex.lock();
        iTimePointSps = m_iTimePointSps;
        bEvokedInput = m_bEvokedInput;
        bRawInput = m_bRawInput;
        iDownSample = m_iDownSample;
        //iNumberChannels = m_bfWeights.noise_cov->names.size();
        //tstep = 1.0f / m_pFiffInfoInput->sfreq;
        //lChNamesFiffInfo = m_pFiffInfoInput->ch_names;
        //lChNamesBfWeights = m_bfWeights.noise_cov->names;
        bUpdateBeamformer = m_bUpdateBeamformer;



        m_qMutex.unlock();


        //HINT: this part is copied from rtcmne, changes: names
        if(bUpdateBeamformer) {

            qDebug() << "[RtBeamformer::run] bUpdateBeamformer = true.";

            uint64_t ms_start_bf_prep = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
            std::cout.precision(64);
            std::cout << "[RtBeamformer::run] sytem_clock start after if(bUpdateBeamformer): " << ms_start_bf_prep << " milliseconds since the Epoch\n";

            m_qMutex.lock();
            pBeamformer = Beamformer::SPtr(new Beamformer(m_bfWeights, regParam, m_sWeightnorm));

            //TODO: This constructor combines different cov and bf weights, delete it
            //pBeamformer = Beamformer::SPtr(new Beamformer(m_bfWeights, *m_pFiffInfo, *m_pFwd, *m_pDataCov, *m_pNoiseCov, regParam, m_sWeightnorm));
            qInfo() << "[RtBeamformer::run] Created new pBeamformer with weightnorm method: " << pBeamformer->getWeightnorm();
            m_bUpdateBeamformer = false;
            m_qMutex.unlock();

            // Set up the inverse according to the parameters.
            // in rtcmne::run: doInverseSetup(1,true) Use 1 nave here because in case of evoked data as input the beamformer will always be updated when the source estimate is calculated (see run method).
            pBeamformer->doInverseSetup(1,true); //HINT: this parameters are not used in doInverseSetup
            //qDebug() << "[RtBeamformer::run] Weightnorm after pBeamformer->doInverseSetup(): " << pBeamformer->getWeightnorm();
            //qDebug() << "[RtBeamformer::run] Finished pBeamformer->doInverseSetup().";

            uint64_t ms_stop_bf_prep = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
            std::cout.precision(64);
            std::cout << "[RtBeamformer::run] sytem_clock stop after ->doInverseSetup(1,true): " << ms_stop_bf_prep << " ms since the Epoch\n";

        }

        //Process data from raw data input
        //TODO
        if(bRawInput && pBeamformer) {

            //TODO: this warning is added and should stay here until Todo above is done
            qWarning() << "[RtBeamformer::run] The option for raw data input is not implemented yet. TODO";
            break;

        }


        //Process data from averaging input
        //HINT: copied from rtcmne::run, changes: names
        if(bEvokedInput && pBeamformer) {
            qInfo() << "[RtBeamformer::run] Process data from evoked input.";
            if(m_pCircularEvokedBuffer->pop(evoked)) {
                // Get the current evoked data

                uint64_t ms_start_bf_appl = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
                std::cout.precision(64);
                std::cout << "[RtBeamformer::run] sytem_clock start after ->pop(evoked): " << ms_start_bf_appl << " milliseconds since the Epoch\n";

                qDebug() << "[RtBeamformer::run] sfreq = " << evoked.info.sfreq;


                if(((skip_count % iDownSample) == 0)) {


                    sourceEstimate = pBeamformer->calculateInverse(evoked);
                    qDebug() << "[RtBeamformer::run] Finished computation of source estimate.";

                    if(!sourceEstimate.isEmpty()) {
                        //qDebug() << "RtBeamformer::run] !sourceEstimate.isEmpty()";
                        if(iTimePointSps < sourceEstimate.data.cols() && iTimePointSps >= 0) {
                            sourceEstimate = sourceEstimate.reduce(iTimePointSps,1);


                            //normalize beamformer source estimate by dividing by the spacial maximum
                            //this step ensures that values are displayed in widget
                            sourceEstimate.data /= sourceEstimate.data.maxCoeff();


                            m_pRTSEOutput->measurementData()->setValue(sourceEstimate);

                            uint64_t ms_stop_bf_appl = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
                            std::cout.precision(64);
                            std::cout << "[RtBeamformer::run] sytem_clock stop after ->setValue(sourceEstimate): " << ms_stop_bf_appl << " ms since the Epoch\n";


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

    qDebug() << "[RtBeamformer::run] Finished running RtBeamformer plugin.";

}

//=============================================================================================================

QString RtBeamformer::getBuildInfo()
{
    //HINT: copied from rtcmne
    return QString(RTBEAMFORMERPLUGIN::buildDateTime()) + QString(" - ")  + QString(RTBEAMFORMERPLUGIN::buildHash());

}



