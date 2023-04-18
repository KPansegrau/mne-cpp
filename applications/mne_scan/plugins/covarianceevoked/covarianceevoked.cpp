//=============================================================================================================
/**
 * @file     covarianceevoked.cpp
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
 * @brief    Definition of the CovarianceEvoked class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "covarianceevoked.h"

#include "FormFiles/covarianceevokedsetupwidget.h"

#include <disp/viewers/covarianceevokedsettingsview.h>

#include <scMeas/realtimeevokedcov.h>
#include <scMeas/realtimeevokedset.h>

#include <rtprocessing/rtcov.h>

#include <fiff/fiff_info.h>
#include <fiff/fiff_cov.h>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtCore/QtPlugin>
#include <QDebug>
#include <QtWidgets>


//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace COVARIANCEEVOKEDPLUGIN;
using namespace DISPLIB;
using namespace SCSHAREDLIB;
using namespace SCMEASLIB;
using namespace UTILSLIB;
using namespace RTPROCESSINGLIB;
using namespace FIFFLIB;
using namespace Eigen;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

CovarianceEvoked::CovarianceEvoked()
    : m_pCircularEvokedBuffer(CircularBuffer<FIFFLIB::FiffEvoked>::SPtr::create(40))
    , m_iEstimationSamples(2000)
    , m_sAvrType("3")
    , m_iNumPreStimSamples(-1)
    , m_iNumAverages(-1)
{
    //HINT: copied from Covariance::Covariance()
    // iEstimationSamples is used in Covariance::run() as minimum number of samples for new calculation of covariance (if the number of sample in the data is lower, empty covariance is returned)
}

//=============================================================================================================

CovarianceEvoked::~CovarianceEvoked()
{
        //HINT: copied from Covariance::Covariance()
    if(this->isRunning())
        stop();
}

//=============================================================================================================

QSharedPointer<AbstractPlugin> CovarianceEvoked::clone() const
{
    QSharedPointer<CovarianceEvoked> pCovarianceEvokedClone(new CovarianceEvoked);
    return pCovarianceEvokedClone;
}

//=============================================================================================================

void CovarianceEvoked::init()
{

    //TODO: add those info messages to all methods here for debugging purposes
    qInfo() << "[CovarianceEvoked::init] Initializing CovarianceEvoked plugin...";

    //Load Settings
    //HINT: copied form Covariance::init()
    QSettings settings("MNECPP");
    //TODO: this 5000 is set for RTC-MNE (new cov every 5000 samples, check if this is okay here
    //TODO: maybe two parameters one for naise and one for data cov update
    m_iEstimationSamples = settings.value(QString("MNESCAN/%1/estimationSamples").arg(this->getName()), 5000).toInt();

    //Input
    //HINT: analog to this part in RtcMne::init()
    m_pCovarianceEvokedInput = PluginInputData<RealTimeEvokedSet>::create(this, "CovarianceEvoked In", "CovarianceEvoked real-time evoked input data");
    connect(m_pCovarianceEvokedInput.data(), &PluginInputConnector::notify,
            this, &CovarianceEvoked::updateRTE, Qt::DirectConnection);
    m_inputConnectors.append(m_pCovarianceEvokedInput);

    //Output
    //HINT: new because we need pair of covariances as output
    m_pCovarianceEvokedOutput = PluginOutputData<RealTimeEvokedCov>::create(this,"CovarianceEvoked Out","CovarianceEvoked output data");
    m_pCovarianceEvokedOutput->measurementData()->setName(this->getName());//Provide name to auto store widget settings
    m_outputConnectors.append(m_pCovarianceEvokedOutput);

    //qDebug() << "[CovarianceEvoked::init] Finished initializing CovarianceEvoked plugin.";

}

//=============================================================================================================

void CovarianceEvoked::initPluginControlWidgets()
{

    qInfo() << "[CovarianceEvoked::initPluginControlWidgets] Initializing control widgets...";

    //HINT: copied from Covariance::initPluginControlWidgets(), added the class CovarianceEvokedSettingsView to libraries/disp/viewers
    if(m_pFiffInfoInput) {
        QList<QWidget*> plControlWidgets;

        CovarianceEvokedSettingsView* pCovarianceEvokedWidget = new CovarianceEvokedSettingsView(QString("MNESCAN/%1").arg(this->getName()));
        connect(this, &CovarianceEvoked::guiModeChanged,
                pCovarianceEvokedWidget, &CovarianceEvokedSettingsView::setGuiMode);
        connect(pCovarianceEvokedWidget, &CovarianceEvokedSettingsView::samplesChanged,
                this, &CovarianceEvoked::changeSamples);
        pCovarianceEvokedWidget->setMinSamples(m_pFiffInfoInput->sfreq);
        pCovarianceEvokedWidget->setCurrentSamples(m_iEstimationSamples);
        pCovarianceEvokedWidget->setObjectName("group_Settings");
        plControlWidgets.append(pCovarianceEvokedWidget);

        emit pluginControlWidgetsChanged(plControlWidgets, this->getName());

        m_bPluginControlWidgetsInit = true;
    }

    //qDebug() << "[CovarianceEvoked::initPluginControlWidgets] Finished initializing control widgets.";

}

//=============================================================================================================

void CovarianceEvoked::unload()
{
        //HINT: copied from Covariance::Covariance()
    // Save Settings
    QSettings settings("MNECPP");
    settings.setValue(QString("MNESCAN/%1/estimationSamples").arg(this->getName()), m_iEstimationSamples);

}

//=============================================================================================================

bool CovarianceEvoked::start()
{
        //HINT: copied from Covariance::Covariance()
    // Start thread
    QThread::start();

    return true;
}

//=============================================================================================================

bool CovarianceEvoked::stop()
{
        //HINT: copied from Covariance::Covariance()
    requestInterruption();
    wait(500);

    m_bPluginControlWidgetsInit = false;

    return true;
}

//=============================================================================================================

AbstractPlugin::PluginType CovarianceEvoked::getType() const
{
        //HINT: copied from Covariance::Covariance()
    return _IAlgorithm;
}

//=============================================================================================================

QString CovarianceEvoked::getName() const
{
    return "Covariance Evoked";
}

//=============================================================================================================

void CovarianceEvoked::showCovarianceEvokedWidget()
{
}


//=============================================================================================================

QWidget* CovarianceEvoked::setupWidget()
{
    //HINT: similar to covariance::setupWidget() but new class CovarianceEvokedSetupWidget
    CovarianceEvokedSetupWidget* setupWidget = new CovarianceEvokedSetupWidget(this);//widget is later distroyed by CentralWidget - so it has to be created everytime new
    return setupWidget;
}

//=============================================================================================================

void CovarianceEvoked::updateRTE(SCMEASLIB::Measurement::SPtr pMeasurement)
{
    //HINT: copied from RtcMne::updateRTE, modifications: no setting of m_bEvokedInput, no channel picking, need to get the pre and post stimulative sample number
    //qDebug() << "[CovarianceEvoked::updateRTE] Updating RTE input...";


    if(QSharedPointer<RealTimeEvokedSet> pRTES = pMeasurement.dynamicCast<RealTimeEvokedSet>()) {

        QStringList lResponsibleTriggerTypes = pRTES->getResponsibleTriggerTypes();
        emit responsibleTriggerTypesChanged(lResponsibleTriggerTypes);

        if(!m_bPluginControlWidgetsInit) {
            initPluginControlWidgets();
        }

        if(!this->isRunning() || !lResponsibleTriggerTypes.contains(m_sAvrType)) {
            return;
        }

        //qDebug() << "TEST";

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

            //HINT: we do not need this since the covariance evoked plugin always gets evoked input
            // m_bEvokedInput = true;
        }

        if(!m_bPluginControlWidgetsInit) {
            initPluginControlWidgets();
        }

        //get number of pre stimulative samples
        //HINT: this part is new, we need the number of pre and post stimulative samples for cutting the correct time window for computation of noise and data covariance matrix
        //TODO: check whether this part works
        m_iNumPreStimSamples = pRTES->getNumPreStimSamples();
        qDebug() << "[CovarianceEvoked::updateRTE] m_iNumPreStimSamples = " << m_iNumPreStimSamples;

        //TODOOOOO: get number of epochs that were averaged somehow so we can scale the covariance matrices
 //       m_iNumAverages = pRTES->getNumAverages(); //


        if(this->isRunning()) {
            for(int i = 0; i < pFiffEvokedSet->evoked.size(); ++i) {
                if(pFiffEvokedSet->evoked.at(i).comment == m_sAvrType) {
                    // Store current evoked as member so we can dispatch it if the time pick by the user changed

                    //HINT: no channel picking in this plugin but in beamformer plugin (s. comment below)
                    //_currentEvoked = pFiffEvokedSet->evoked.at(i).pick_channels(m_qListPickChannels);

                    // Please note that we do not need a copy here since this function will block until
                    // the buffer accepts new data again. Hence, the data is not deleted in the actual
                    // Measurement function after it emitted the notify signal.

                    //HINT: we do not need channel picking here since it is perfomed in the rtbeamformer routines when it is known which channels are common to all beamformer inputs
                    //original in RtcMne: while(!m_pCircularEvokedBuffer->push(pFiffEvokedSet->evoked.at(i).pick_channels(m_qListPickChannels))) {
                    while(!m_pCircularEvokedBuffer->push(pFiffEvokedSet->evoked.at(i))) {
                        //Do nothing until the circular buffer is ready to accept new data again
                    }

                        //qDebug()<<"CovarianceEvoked::updateRTE - average found type" << m_sAvrType;
                        break;
                    }
                }
            }
        }

    //qDebug() << "[CovarianceEvoked::updateRTE] Finished updating RTE input.";

}



//=============================================================================================================


void CovarianceEvoked::changeSamples(qint32 samples)
{
    m_iEstimationSamples = samples;
    qDebug() << "[CovarianceEvoked::changeSamples] Changed samples to " << m_iEstimationSamples;
}

//=============================================================================================================

void CovarianceEvoked::onTriggerTypeChanged(const QString& triggerType)
{
    m_sAvrType = triggerType;
    qDebug() << "[CovarianceEvoked::onTriggerTypeChanged] Changed trigger type to " << m_sAvrType;

}

//=============================================================================================================

void CovarianceEvoked::run()
{

    qInfo() << "[CovarianceEvoked::run] Running covariance evoked plugin...";

    // Wait for fiff info
    //HINT: this part is copied from Covariance::run()
    while(true) {
        m_qMutex.lock();
        if(m_pFiffInfoInput) {
            m_qMutex.unlock();
            break;
        }
        m_qMutex.unlock();
        msleep(100);
    }

    //init
    FiffEvoked evokedData;
    MatrixXd matPreStimData;
    MatrixXd matPostStimData;
    FiffCov fiffNoiseCov;
    FiffCov fiffDataCov;
    int iEstimationSamples;
    //TODO: check if its correct to use the same fiff info in both cases
    RTPROCESSINGLIB::RtCov rtNoiseCov(m_pFiffInfoInput);
    RTPROCESSINGLIB::RtCov rtDataCov(m_pFiffInfoInput);


    qInfo() << "[CovarianceEvoked::run] Start processing data.";

    //start processing data
    while(!isInterruptionRequested()) {

        // Get the current data
        if(m_pCircularEvokedBuffer->pop(evokedData)) {
            m_qMutex.lock();
            iEstimationSamples = m_iEstimationSamples;
            m_qMutex.unlock();

            //qDebug() << "[CovarianceEvoked::run] m_iNumPreStimSamples = " << m_iNumPreStimSamples;

            //qDebug() << "[CovarianceEvoked::run] evokedData.data " << evokedData.data.rows() << " x " << evokedData.data.cols();

            //split evoked data into pre and post stimulative part
            matPreStimData = evokedData.data.leftCols(m_iNumPreStimSamples);
            matPostStimData = evokedData.data.rightCols(evokedData.data.cols() - m_iNumPreStimSamples);


            //qDebug() << "[CovarianceEvoked::run] Pre stim matrix " << matPreStimData.rows() << " x " << matPreStimData.cols();
            //qDebug() << "[CovarianceEvoked::run] Post stim matrix " << matPostStimData.rows() << " x " << matPostStimData.cols();


            //TODO: call estimateCovariance for parts of matData
            //TODO: check whether part of evokedData causes issues with rtNoiseCov object because info is not manipulated
            fiffNoiseCov = rtNoiseCov.estimateCovariance(matPreStimData, iEstimationSamples);
            fiffDataCov = rtDataCov.estimateCovariance(matPostStimData, iEstimationSamples);

            //set kind of data covariance manually (during estimateCovariance it is set to FIFFV_MNE_NOISE_COV)
            fiffDataCov.kind = FIFFV_MNE_SIGNAL_COV;

            qDebug() << "[CovarianceEvoked::run] Estimated noise covariance matrix " << fiffNoiseCov.data.rows() << " x " << fiffNoiseCov.data.cols();
            qDebug() << "[CovarianceEvoked::run] Estimated data covariance matrix " << fiffDataCov.data.rows() << " x " << fiffDataCov.data.cols();


            //TODOOOOO decide whether we need this scaling here!
/*            if(m_iNumAverages <= 0)
            {
                qWarning("[CovarianceEvoked::run] The number of averages should be positive. No scaling is performed\n");
            }else{
                qDebug() << "[CovarianceEvoked::run] Scale data and noise covariance matrix...\n";

                //copied from rtcmne (MNEInverseOperator::prepare_inverse_operator) and modified for noise and data cov scaling
                float fScale     = 1/((float)m_iNumAverages); //scaling factor
                noiseCov.data  *= fScale; //scaling data and eigenvalues of noise covariance matrix
                noiseCov.eig   *= fScale;
                dataCov.data *= fScale; //scaling data and eigenvalues of data covariance matrix
                dataCov.eig *= fScale; //TODO: check whether we need to scale the eigenvalues of data covariance matrix (are they used somewhere else? do we need them scaled for consistency)

                qInfo("[MNEBeamformerWeights::make_beamformer_weights] Scaled noise and data covariance with scaling factor = %f (number of averages = %d)\n",fScale,m_iNumAverages);

            }
*/

            if(!fiffNoiseCov.names.isEmpty() && !fiffDataCov.names.isEmpty()) {
                m_pCovarianceEvokedOutput->measurementData()->setValue(fiffNoiseCov,fiffDataCov);
                qDebug() << "[CovarianceEvoked::run] Wrote covariance matrices to plug-in output.";
            }
        }
    }

    qDebug() << "[CovarianceEvoked::run] Finished running covariance evoked plugin.";
}

//=============================================================================================================

QString CovarianceEvoked::getBuildInfo()
{
        //HINT: copied from Covariance::Covariance()
    return QString(COVARIANCEEVOKEDPLUGIN::buildDateTime()) + QString(" - ")  + QString(COVARIANCEEVOKEDPLUGIN::buildHash());
}


