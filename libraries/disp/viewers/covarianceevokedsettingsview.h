//=============================================================================================================
/**
 * @file     covarianceevokedsettingsview.h
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de>
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
 * @brief     CovarianceEvokedSettingsView class declaration.
 *
 */

#ifndef COVARIANCEEVOKEDSETTINGSVIEW_H
#define COVARIANCEEVOKEDSETTINGSVIEW_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../disp_global.h"
#include "abstractview.h"


//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QWidget>
#include <QSpinBox>
#include <QPair>

#include <QComboBox>
#include <QCheckBox>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

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
 * User GUI control for Covariance Evoked estimation.
 *
 * @brief User GUI control for Covariance Evoked estimation.
 */
class DISPSHARED_EXPORT CovarianceEvokedSettingsView : public AbstractView
{
    Q_OBJECT

public:
    typedef QSharedPointer<CovarianceEvokedSettingsView> SPtr;            /**< Shared pointer type for CovarianceEvokedSettingsView. */
    typedef QSharedPointer<const CovarianceEvokedSettingsView> ConstSPtr; /**< Const shared pointer type for CovarianceEvokedSettingsView. */

    //=========================================================================================================
    /**
    * Constructs a CovarianceEvokedSettingsView object.
    */
    explicit CovarianceEvokedSettingsView(const QString& sSettingsPath = "",
                                          QWidget *parent = 0);

    //=========================================================================================================
    /**
     * Destroys the CovarianceSettingsView.
     */
    ~CovarianceEvokedSettingsView();

    //=========================================================================================================
    /**
     * Set current samples to gather until a new covariance is calculated.
     *
     * @param[in] iSamples     new number samples.
     */
    void setCurrentSamples(int iSamples);

    //=========================================================================================================
    /**
     * Set minimum number of samples to gather until a new covariance is calculated.
     *
     * @param[in] iSamples     new minimum number of samples.
     */
    void setMinSamples(int iSamples);

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

signals:
    void samplesChanged(int iSamples);

private:
    QSpinBox*       m_pSpinBoxNumSamples;
    QString         m_sSettingsPath;            /**< The settings path to store the GUI settings to. */

};

//=============================================================================================================
// INLINE DEFINITIONS
//=============================================================================================================


} //namespace

#endif // COVARIANCEEVOKEDSETTINGSVIEW_H

