//=============================================================================================================
/**
 * @file     beamformersettingsview.cpp
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
 * @brief    Definition of the BeamformerSettingsView Class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "beamformersettingsview.h"

#include "ui_beamformersettingsview.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSettings>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace DISPLIB;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

BeamformerSettingsView::BeamformerSettingsView(const QString& sSettingsPath,
                                                 QWidget *parent,
                                                 Qt::WindowFlags f)
: AbstractView(parent, f)
, m_pUi(new Ui::BeamformerSettingsViewWidget)
{
    m_sSettingsPath = sSettingsPath;
    m_pUi->setupUi(this);

    connect(m_pUi->m_comboBox_method, static_cast<void (QComboBox::*)(const QString&)>(&QComboBox::currentTextChanged),
            this, &BeamformerSettingsView::onWeightnormChanged);

    connect(m_pUi->m_comboBox_triggerType, static_cast<void (QComboBox::*)(const QString&)>(&QComboBox::currentTextChanged),
            this, &BeamformerSettingsView::onTriggerTypeChanged);

    connect(m_pUi->m_spinBox_timepoint, &QSpinBox::editingFinished,
            this, &BeamformerSettingsView::onTimePointValueChanged);

    this->setWindowTitle("Beamformer Settings");
    this->setMinimumWidth(330);
    this->setMaximumWidth(330);

    loadSettings();
}

//=============================================================================================================

BeamformerSettingsView::~BeamformerSettingsView()
{
    saveSettings();
    delete m_pUi;
}

//=============================================================================================================

void BeamformerSettingsView::setTriggerTypes(const QStringList& lTriggerTypes)
{
    for(const QString &sTriggerType : lTriggerTypes) {
        if(m_pUi->m_comboBox_triggerType->findText(sTriggerType) == -1) {
            m_pUi->m_comboBox_triggerType->addItem(sTriggerType);
        }
    }
}

//=============================================================================================================

void BeamformerSettingsView::saveSettings()
{
    if(m_sSettingsPath.isEmpty()) {
        return;
    }

    // Save Settings
    QSettings settings("MNECPP");
}

//=============================================================================================================

void BeamformerSettingsView::loadSettings()
{
    if(m_sSettingsPath.isEmpty()) {
        return;
    }

    // Load Settings
    QSettings settings("MNECPP");
}

//=============================================================================================================

void BeamformerSettingsView::updateGuiMode(GuiMode mode)
{
    switch(mode) {
        case GuiMode::Clinical:
            break;
        default: // default is research mode
            break;
    }
}

//=============================================================================================================

void BeamformerSettingsView::updateProcessingMode(ProcessingMode mode)
{
    switch(mode) {
        case ProcessingMode::Offline:
            break;
        default: // default is realtime mode
            break;
    }
}

//=============================================================================================================

void BeamformerSettingsView::onWeightnormChanged(const QString& weightnorm)
{
    emit weightnormChanged(weightnorm);
}

//=============================================================================================================

void BeamformerSettingsView::onTriggerTypeChanged(const QString& sTriggerType)
{
    emit triggerTypeChanged(sTriggerType);
}

//=============================================================================================================

void BeamformerSettingsView::onTimePointValueChanged()
{
    emit timePointChanged(m_pUi->m_spinBox_timepoint->value());
}

//=============================================================================================================

void BeamformerSettingsView::clearView()
{

}

