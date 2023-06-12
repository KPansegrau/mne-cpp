//=============================================================================================================
/**
 * @file     covarianceevokedsettingsview.cpp
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de>
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
 * @brief    CovarianceEvokedSettingsView class definition.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "covarianceevokedsettingsview.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QGridLayout>
#include <QSpinBox>
#include <QLabel>
#include <QSettings>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace DISPLIB;

//=============================================================================================================
// DEFINE GLOBAL METHODS
//=============================================================================================================

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

CovarianceEvokedSettingsView::CovarianceEvokedSettingsView(const QString& sSettingsPath,
                                                            QWidget *parent)
: AbstractView(parent)
, m_sSettingsPath(sSettingsPath)
{
    this->setWindowTitle("Covariance Evoked Settings");
    this->setMinimumWidth(330);
    this->setMaximumWidth(330);

    QGridLayout* t_pGridLayout = new QGridLayout;

    QLabel* t_pLabelNumSamples = new QLabel;
    t_pLabelNumSamples->setText("Number of Samples");
    t_pGridLayout->addWidget(t_pLabelNumSamples,0,0,1,1);

    qint32 minSamples = 240;

    m_pSpinBoxNumSamples = new QSpinBox;
    m_pSpinBoxNumSamples->setMinimum(minSamples);
    m_pSpinBoxNumSamples->setMaximum(minSamples*60);
    m_pSpinBoxNumSamples->setSingleStep(minSamples);
    connect(m_pSpinBoxNumSamples, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
            this, &CovarianceEvokedSettingsView::samplesChanged);
    t_pGridLayout->addWidget(m_pSpinBoxNumSamples,0,1,1,1);
    this->setLayout(t_pGridLayout);

    loadSettings();
}

//=============================================================================================================

CovarianceEvokedSettingsView::~CovarianceEvokedSettingsView()
{
    saveSettings();
}

//=============================================================================================================

void CovarianceEvokedSettingsView::setCurrentSamples(int iSamples)
{
    m_pSpinBoxNumSamples->setValue(iSamples);
}

//=============================================================================================================

void CovarianceEvokedSettingsView::setMinSamples(int iSamples)
{
    m_pSpinBoxNumSamples->setMinimum(iSamples);
    m_pSpinBoxNumSamples->setMaximum(iSamples*60);
}

//=============================================================================================================

void CovarianceEvokedSettingsView::saveSettings()
{
    if(m_sSettingsPath.isEmpty()) {
        return;
    }

    // Save Settings
    QSettings settings("MNECPP");
}

//=============================================================================================================

void CovarianceEvokedSettingsView::loadSettings()
{
    if(m_sSettingsPath.isEmpty()) {
        return;
    }

    // Load Settings
    QSettings settings("MNECPP");
}

//=============================================================================================================

void CovarianceEvokedSettingsView::updateGuiMode(GuiMode mode)
{
    switch(mode) {
        case GuiMode::Clinical:
            break;
        default: // default is research mode
            break;
    }
}

//=============================================================================================================

void CovarianceEvokedSettingsView::updateProcessingMode(ProcessingMode mode)
{
    switch(mode) {
        case ProcessingMode::Offline:
            break;
        default: // default is realtime mode
            break;
    }
}

//=============================================================================================================

void CovarianceEvokedSettingsView::clearView()
{

}



