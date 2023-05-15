//=============================================================================================================
/**
 * @file     beamformersettingsview.h
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
 * @brief    Declaration of the BeamformerSettingsView Class.
 *
 */

#ifndef BEAMFORMERSETTINGSVIEW_H
#define BEAMFORMERSETTINGSVIEW_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../disp_global.h"
#include "abstractview.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

namespace Ui {
    class BeamformerSettingsViewWidget;
}

//=============================================================================================================
// DEFINE NAMESPACE DISPLIB
//=============================================================================================================

namespace DISPLIB
{

//=============================================================================================================
// DISPLIB FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * DECLARE CLASS BeamformerSettingsView
 *
 * @brief The BeamformerSettingsView class provides a view to control settings for estiamting functional connectivity
 */
class DISPSHARED_EXPORT BeamformerSettingsView : public AbstractView
{
    Q_OBJECT

public:
    typedef QSharedPointer<BeamformerSettingsView> SPtr;              /**< Shared pointer type for BeamformerSettingsView. */
    typedef QSharedPointer<const BeamformerSettingsView> ConstSPtr;   /**< Const shared pointer type for BeamformerSettingsView. */

    //=========================================================================================================
    /**
     * Constructs a BeamformerSettingsView which is a child of parent.
     *
     * @param[in] parent        parent of widget.
     */
    BeamformerSettingsView(const QString& sSettingsPath = "",
                            QWidget *parent = 0,
                            Qt::WindowFlags f = Qt::Widget);

    //=========================================================================================================
    /**
     * Destroys the BeamformerSettingsView.
     */
    ~BeamformerSettingsView();

    //=========================================================================================================
    /**
     * Destroys the BeamformerSettingsView.
     *
     * @param[in] lTriggerTypes        The new trigger types.
     */
    void setTriggerTypes(const QStringList& lTriggerTypes);

    //=========================================================================================================
    /**
     * Saves all important settings of this view via QSettings.
     */
    void saveSettings();

    //=========================================================================================================
    /**
     * Loads and inits all important settings of this view via QSettings.
     */
    void loadSettings();

    //=========================================================================================================
    /**
     * Clears the view
     */
    void clearView();

protected:
    //=========================================================================================================
    /**
     * Update the views GUI based on the set GuiMode (Clinical=0, Research=1).
     *
     * @param[in] mode     The new mode (Clinical=0, Research=1).
     */
    void updateGuiMode(GuiMode mode);

    //=========================================================================================================
    /**
     * Update the views GUI based on the set ProcessingMode (RealTime=0, Offline=1).
     *
     * @param[in] mode     The new mode (RealTime=0, Offline=1).
     */
    void updateProcessingMode(ProcessingMode mode);

    //=========================================================================================================
    /**
     * Slot called when the weight normalization method changed.
     *
     * @param[in] weightnorm        The new weight normalization method.
     */
    void onWeightnormChanged(const QString& weightnorm);

    //=========================================================================================================
    /**
     * Slot called when the trigger type changed.
     *
     * @param[in] sTriggerType        The new trigger type.
     */
    void onTriggerTypeChanged(const QString& sTriggerType);

    //=========================================================================================================
    /**
     * Slot called when the time point changes.
     */
    void onTimePointValueChanged();

    Ui::BeamformerSettingsViewWidget* m_pUi;

signals:
    //=========================================================================================================
    /**
     * Emit signal whenever the weight normalization method changed.
     *
     * @param[in] method        The new weight normalization method.
     */
    void weightnormChanged(const QString& weightnorm);

    //=========================================================================================================
    /**
     * Emit signal whenever the trigger type changed.
     *
     * @param[in] triggerType        The new trigger type.
     */
    void triggerTypeChanged(const QString& triggerType);

    //=========================================================================================================
    /**
     * Emit signal whenever the time point changed.
     *
     * @param[in] iTimePoint        The new time point.
     */
    void timePointChanged(int iTimePoint);
};
} // NAMESPACE

#endif // BeamformerSettingsView_H
