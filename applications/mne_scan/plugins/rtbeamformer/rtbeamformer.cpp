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


//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

RtBeamformer::RtBeamformer()
    : m_pCircularMatrixBuffer(CircularBuffer_Matrix_double::SPtr(new CircularBuffer_Matrix_double(40)))
    , m_bRawInput(false)
    , m_iNumAverages(1)
{

}

//=============================================================================================================

RtBeamformer::~RtBeamformer()
{
    //TODO
}

//=============================================================================================================

QSharedPointer<AbstractPlugin> RtBeamformer::clone() const
{
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

    //Output
    m_pRTSEOutput = PluginOutputData<RealTimeSourceEstimate>::create(this, "RTBEAMFORMER Out", "RTBEAMFORMER output data");
    m_outputConnectors.append(m_pRTSEOutput);
    m_pRTSEOutput->measurementData()->setName(this->getName());//Provide name to auto store widget settings


    // Set the annotation and surface data and mri-head transformation

}

//=============================================================================================================

void RtBeamformer::unload()
{
    //TODO
}

//=============================================================================================================

bool RtBeamformer::start()
{
    //TODO
    return true;
}

//=============================================================================================================

bool RtBeamformer::stop()
{
    //TODO
    return true;
}

//=============================================================================================================

AbstractPlugin::PluginType RtBeamformer::getType() const
{
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

    //TODO
}

//=============================================================================================================

void RtBeamformer::updateRTMSA(SCMEASLIB::Measurement::SPtr pMeasurement)
{

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
                    //linker errors occur in the following line
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

void RtBeamformer::run()
{
     //TODO
}

//=============================================================================================================

QString RtBeamformer::getBuildInfo()
{
    return QString(RTBEAMFORMERPLUGIN::buildDateTime()) + QString(" - ")  + QString(RTBEAMFORMERPLUGIN::buildHash());

}



