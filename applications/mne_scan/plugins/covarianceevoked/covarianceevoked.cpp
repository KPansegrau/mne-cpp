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
using namespace SCSHAREDLIB;
using namespace UTILSLIB;


//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

CovarianceEvoked::CovarianceEvoked()
: m_iEstimationSamples(2000)
, m_pCircularBuffer(CircularBuffer_Matrix_double::SPtr::create(40))
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
//TODO
}

//=============================================================================================================

void CovarianceEvoked::initPluginControlWidgets()
{
//TODO
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

QWidget* CovarianceEvoked::setupWidget()
{
    //HINT: similar to covariance::setupWidget() but new class CovarianceEvokedSetupWidget
    CovarianceEvokedSetupWidget* setupWidget = new CovarianceEvokedSetupWidget(this);//widget is later distroyed by CentralWidget - so it has to be created everytime new
    return setupWidget;
}


//=============================================================================================================

void CovarianceEvoked::run()
{

}

//=============================================================================================================

QString CovarianceEvoked::getBuildInfo()
{
        //HINT: copied from Covariance::Covariance()
    return QString(COVARIANCEEVOKEDPLUGIN::buildDateTime()) + QString(" - ")  + QString(COVARIANCEEVOKEDPLUGIN::buildHash());
}


