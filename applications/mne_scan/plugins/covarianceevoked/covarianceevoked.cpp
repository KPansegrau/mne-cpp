//=============================================================================================================
/**
 * @file     covarianceevoked.cpp
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

#include <iostream>
#include <fstream>


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
    : m_pCircularEvokedBuffer(CircularBuffer<FIFFLIB::FiffEvoked>::SPtr::create(40)) //same as in rtcmne.cpp: buffer for evoked data input
    , m_iEstimationSamples(2400) //2400 because pre and post stimulative part are 240 samples each
    , m_sAvrType("1") //trigger type is 1 in simulated file, TODO: enable that the trigger type is read from the input evoked data somehow and stored here during updateRTE
    , m_iNumPreStimSamples(-1) //this is set according to the length of pre and post stimulative parts
{
}

//=============================================================================================================

CovarianceEvoked::~CovarianceEvoked()
{
    //HINT: copied from Covariance::Covariance(), completely analog
    if(this->isRunning())
        stop();
}

//=============================================================================================================

QSharedPointer<AbstractPlugin> CovarianceEvoked::clone() const
{
    //HINT: completely analog to Covariance::clone
    QSharedPointer<CovarianceEvoked> pCovarianceEvokedClone(new CovarianceEvoked);
    return pCovarianceEvokedClone;
}

//=============================================================================================================

void CovarianceEvoked::init()
{

    qInfo() << "[CovarianceEvoked::init] Initializing CovarianceEvoked plugin...";

    //Load Settings
    //HINT: copied form Covariance::init()
    QSettings settings("MNECPP");
    m_iEstimationSamples = settings.value(QString("MNESCAN/%1/estimationSamples").arg(this->getName()), 2400).toInt();

    //Input
    //HINT: analog to this part in RtcMne::init()
    m_pCovarianceEvokedInput = PluginInputData<RealTimeEvokedSet>::create(this, "CovarianceEvoked In", "CovarianceEvoked real-time evoked input data");
    connect(m_pCovarianceEvokedInput.data(), &PluginInputConnector::notify,
            this, &CovarianceEvoked::updateRTE, Qt::DirectConnection);
    m_inputConnectors.append(m_pCovarianceEvokedInput);

    //Output
    //HINT: new because we need pair of noise and data covariance as output
    m_pCovarianceEvokedOutput = PluginOutputData<RealTimeEvokedCov>::create(this,"CovarianceEvoked Out","CovarianceEvoked output data");
    m_pCovarianceEvokedOutput->measurementData()->setName(this->getName());//Provide name to auto store widget settings
    m_outputConnectors.append(m_pCovarianceEvokedOutput);

}

//=============================================================================================================

void CovarianceEvoked::initPluginControlWidgets()
{

    //HINT: completely analog to Covariance::initPluginControlWidgets
    qInfo() << "[CovarianceEvoked::initPluginControlWidgets] Initializing control widgets...";

    //HINT: copied from Covariance::initPluginControlWidgets()
    if(m_pFiffInfoInput) {

        //TODO: delete later
//        qDebug() << "[CovarianceEvoked::initPluginControlWidgets] if(m_pFiffInfoInput) true.";

        QList<QWidget*> plControlWidgets;

        CovarianceEvokedSettingsView* pCovarianceEvokedWidget = new CovarianceEvokedSettingsView(QString("MNESCAN/%1").arg(this->getName()));
        connect(this, &CovarianceEvoked::guiModeChanged,
                pCovarianceEvokedWidget, &CovarianceEvokedSettingsView::setGuiMode);
        connect(pCovarianceEvokedWidget, &CovarianceEvokedSettingsView::samplesChanged,
                this, &CovarianceEvoked::changeSamples);

        //HINT: this is new because we need to set the triggert type member when the trigger type of incoming evoked data changes
        connect(this, &CovarianceEvoked::responsibleTriggerTypesChanged,
                this, &CovarianceEvoked::changeTriggerType);


        pCovarianceEvokedWidget->setMinSamples(m_pFiffInfoInput->sfreq);
        pCovarianceEvokedWidget->setCurrentSamples(m_iEstimationSamples);
        pCovarianceEvokedWidget->setObjectName("group_Settings");
        plControlWidgets.append(pCovarianceEvokedWidget);

        emit pluginControlWidgetsChanged(plControlWidgets, this->getName());

        m_bPluginControlWidgetsInit = true;
    }

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
    //HINT: copied from Covariance::start()
    // Start thread
    QThread::start();

    return true;
}

//=============================================================================================================

bool CovarianceEvoked::stop()
{
    //HINT: copied from Covariance::stop()
    requestInterruption();
    wait(500);

    m_bPluginControlWidgetsInit = false;

    return true;
}

//=============================================================================================================

AbstractPlugin::PluginType CovarianceEvoked::getType() const
{
    //HINT: copied from Covariance::getType()
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
    //HINT: copied from RtcMne::updateRTE with some changes
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

            //HINT: we do not need this since the covariance evoked plugin always gets evoked input
            // m_bEvokedInput = true;
        }

        if(!m_bPluginControlWidgetsInit) {
            initPluginControlWidgets();
        }

        //get number of pre stimulative samples (should equal number of poststimulative samples)
        //HINT: this part is new, we need the number of pre and post stimulative samples for cutting the correct time window for computation of noise and data covariance matrix
        m_iNumPreStimSamples = pRTES->getNumPreStimSamples();
//        qDebug() << "[CovarianceEvoked::updateRTE] m_iNumPreStimSamples = " << m_iNumPreStimSamples;

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

}



//=============================================================================================================

void CovarianceEvoked::changeSamples(qint32 samples)
{
    m_iEstimationSamples = samples;
    qDebug() << "[CovarianceEvoked::changeSamples] Changed samples to " << m_iEstimationSamples;
}


//=============================================================================================================

void CovarianceEvoked::changeTriggerType(const QStringList& lTriggerType)
{
    //TODO: this method can only handle the first trigger type. Change it so that full list can be handeled
    m_sAvrType = lTriggerType.at(0);
//    qDebug() << "[CovarianceEvoked::changeTriggerType] Changed trigger type to " << m_sAvrType;

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
    FiffCov fiffNoiseCov;
    FiffCov fiffDataCov;
    m_qMutex.lock();
    int iEstimationSamples;
    m_qMutex.unlock();
    RTPROCESSINGLIB::RtCov rtNoiseCov(m_pFiffInfoInput);
    RTPROCESSINGLIB::RtCov rtDataCov(m_pFiffInfoInput);

    qInfo() << "[CovarianceEvoked::run] Start processing data.";

    std::ofstream timeFileStart;
    std::ofstream timeFileStop;

    std::ofstream testNoiseCovEvokedOut;
    std::ofstream testDataCovEvokedOut;

    testNoiseCovEvokedOut.open("testNoiseCovEvokedOut.txt", std::ofstream::trunc);
    testNoiseCovEvokedOut.close();

    testDataCovEvokedOut.open("testDataCovEvokedOut.txt", std::ofstream::trunc);
    testDataCovEvokedOut.close();

    //start processing data
    while(!isInterruptionRequested()) {

        // Get the current data
        if(m_pCircularEvokedBuffer->pop(evokedData)) {

            //for performance evaluation only
            timeFileStart.open("testTimingCovarianceEvokedStart.txt", std::ios::app);
            uint64_t time_start = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
            timeFileStart <<  time_start << '\n';
            timeFileStart.close();

            m_qMutex.lock();
            iEstimationSamples = m_iEstimationSamples;
            m_qMutex.unlock();

//            qDebug() << "[CovarianceEvoked::run] evokedData.data " << evokedData.data.rows() << " x " << evokedData.data.cols();
            qDebug() << "[CovarianceEvoked::run] iEstimationSamples:  " << iEstimationSamples;
//            qDebug() << "[CovarianceEvoked::run] m_iNumPreStimSamples:  " << m_iNumPreStimSamples;

            //estimate covariance matrices
            fiffNoiseCov = rtNoiseCov.estimateCovariance(evokedData.data.leftCols(m_iNumPreStimSamples), iEstimationSamples);
            fiffDataCov = rtDataCov.estimateCovariance(evokedData.data.rightCols(evokedData.data.cols() - m_iNumPreStimSamples), iEstimationSamples);

            //set kind of data covariance manually (during estimateCovariance it is set to FIFFV_MNE_NOISE_COV)
            fiffDataCov.kind = FIFFV_MNE_SIGNAL_COV;

//            qDebug() << "[CovarianceEvoked::run] Estimated noise covariance matrix " << fiffNoiseCov.data.rows() << " x " << fiffNoiseCov.data.cols();
//            qDebug() << "[CovarianceEvoked::run] Estimated data covariance matrix " << fiffDataCov.data.rows() << " x " << fiffDataCov.data.cols();

            //TODO: only for debugging, delete later
            //this creates a spatial filter when used for beamformer weight computation
//            fiffDataCov.data = MatrixXd::Identity(fiffDataCov.data.rows(),fiffDataCov.data.cols());
//            VectorXd vecSimulatedStdNoise(fiffNoiseCov.data.rows());
//            for(int iChan = 0; iChan < vecSimulatedStdNoise.size()-2; iChan += 3){
//                vecSimulatedStdNoise[iChan] = pow(5e-13,2); //gradiometer sensor noise variance
//                vecSimulatedStdNoise[iChan+1] = pow(5e-13,2);
//                vecSimulatedStdNoise[iChan+2] = pow(20e-15,2); //magnetometer senosor noise variance
//            }
//            MatrixXd matSimulatedNoiseCov = vecSimulatedStdNoise.asDiagonal();
//            fiffNoiseCov.data = matSimulatedNoiseCov;

            //TODO: only for debugging, delete later
//            fiffDataCov.data = fiffNoiseCov.data;
//            fiffDataCov.data = matSimulatedNoiseCov;

            //set diagonal values to zero since this is autokovarianz (Hint from OIPE 2023)
            //TODO: here must be the cause why MNE Scan crashes at the really beginning when pipeline is run
            //TODO: try double loop instead
//            for(int i = 0; i < fiffNoiseCov.data.size(); i++){
//                fiffNoiseCov.data(i,i) = 0.0;
//                fiffDataCov.data(i,i) = 0.0;
//            }

            //TODO: this does not crash, but fiffDataCov.data and fiffNoiseCov.data is empty :/
            for(int iRow = 0; iRow < fiffDataCov.data.rows(); iRow++){
                for(int iCol = 0; iCol < fiffDataCov.data.cols(); iCol++){
                    if(iRow == iCol){
                        fiffDataCov.data(iRow, iCol) = 0.0;
                    }
                }
            }



            //TODO for debugging only, delete later
            testNoiseCovEvokedOut.open("testNoiseCovEvokedOut.txt", std::ios::app);
            for(int iRow = 0; iRow < fiffNoiseCov.data.rows(); iRow++){

                for(int iCol = 0; iCol < fiffNoiseCov.data.cols(); iCol++){

                testNoiseCovEvokedOut << fiffNoiseCov.data(iRow,iCol) << "    ";

                }
                testNoiseCovEvokedOut << '\n' ;

            }
            testNoiseCovEvokedOut << "xxxxxxxxx" << '\n' ;
            testNoiseCovEvokedOut.close();

            //TODO for debugging only, delete later
            testDataCovEvokedOut.open("testDataCovEvokedOut.txt", std::ios::app);
            for(int iRow = 0; iRow < fiffDataCov.data.rows(); iRow++){

                for(int iCol = 0; iCol < fiffDataCov.data.cols(); iCol++){

                testDataCovEvokedOut << fiffDataCov.data(iRow,iCol) << "    ";

                }
                testDataCovEvokedOut << '\n' ;

            }
            testDataCovEvokedOut << "xxxxxxxxx" << '\n' ;
            testDataCovEvokedOut.close();



            if(!fiffNoiseCov.names.isEmpty() && !fiffDataCov.names.isEmpty()) {


//                //TODO: for debugging only, delete later
//                qDebug() << "[CovarianceEvoked::run] fiffDataCov.max = " << fiffDataCov.data.maxCoeff();
//                qDebug() << "[CovarianceEvoked::run] fiffDataCov.mean = " << fiffDataCov.data.mean();
//                qDebug() << "[CovarianceEvoked::run] fiffDataCov.min = " << fiffDataCov.data.minCoeff();
////                [DEBUG] [CovarianceEvoked::run] fiffDataCov.max =  0.00415508
////                [DEBUG] [CovarianceEvoked::run] fiffDataCov.mean =  2.93903e-08
////                [DEBUG] [CovarianceEvoked::run] fiffDataCov.min =  -2.88163e-09
//                FiffCov fiffDataCovConst = fiffDataCov;
//                fiffDataCovConst.data = MatrixXd::Ones(fiffDataCov.data.rows(),fiffDataCov.data.cols()) * 2.93903e-08;



                m_pCovarianceEvokedOutput->measurementData()->setValue(fiffNoiseCov,fiffDataCov);

                // for performance evaluation only
                uint64_t time_stop = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
                timeFileStop.open("testTimingCovarianceEvokedStop.txt", std::ios::app);
                timeFileStop << time_stop << '\n';
                timeFileStop.close();

                qDebug() << "[CovarianceEvoked::run] Wrote new covariance matrices to plug-in output.";

            }else{

                //for performance evaluation only
                timeFileStop.open("testTimingCovarianceEvokedStop.txt", std::ios::app);
                timeFileStop << 0 << '\n';
                timeFileStop.close();
            }

        }
    }

}

//=============================================================================================================

QString CovarianceEvoked::getBuildInfo()
{
    //HINT: copied from Covariance::Covariance()
    return QString(COVARIANCEEVOKEDPLUGIN::buildDateTime()) + QString(" - ")  + QString(COVARIANCEEVOKEDPLUGIN::buildHash());
}


