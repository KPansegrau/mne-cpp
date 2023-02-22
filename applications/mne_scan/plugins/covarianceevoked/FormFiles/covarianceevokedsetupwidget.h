//=============================================================================================================
/**
 * @file     covarianceevokedsetupwidget.h
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de;
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
 * @brief    Contains the declaration of the CovarianceEvokedSetupWidget class.
 *
 */

#ifndef COVARIANCEEVOKEDSETUPWIDGET_H
#define COVARIANCEEVOKEDSETUPWIDGET_H


//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "ui_covarianceevokedsetupwidget.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtWidgets>

//=============================================================================================================
// DEFINE NAMESPACE COVARIANCEPLUGIN
//=============================================================================================================

namespace COVARIANCEEVOKEDPLUGIN
{

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

class CovarianceEvoked;

//=============================================================================================================
/**
 * DECLARE CLASS CovarianceEvokedSetupWidget
 *
 * @brief The CovarianceEvokedSetupWidget class provides the CovarianceEvokedToolbox configuration window.
 */
class CovarianceEvokedSetupWidget : public QWidget
{
    Q_OBJECT

public:
    //=========================================================================================================
    /**
     * Constructs a CovarianceEvokedSetupWidget which is a child of parent.
     *
     * @param[in] toolbox a pointer to the corresponding CovarianceEvoked toolbox.
     * @param[in] parent pointer to parent widget; If parent is 0, the new CovarianceEvokedSetupWidget becomes a window. If parent is another widget, CovarianceEvokedSetupWidget becomes a child window inside parent. CovarianceEvokedSetupWidget is deleted when its parent is deleted.
     */
    CovarianceEvokedSetupWidget(CovarianceEvoked* toolbox, QWidget *parent = 0);

    //=========================================================================================================
    /**
     * Destroys the CovarianceEvokedSetupWidget.
     * All CovarianceEvokedSetupWidget's children are deleted first. The application exits if CovarianceEvokedSetupWidget is the main widget.
     */
    ~CovarianceEvokedSetupWidget();

private:

    CovarianceEvoked* m_pCovarianceEvoked;        /**< Holds a pointer to corresponding CovarianceEvoked.*/

    Ui::CovarianceEvokedSetupWidgetClass ui;   /**< Holds the user interface for the CovarianceEvokedSetupWidget.*/
};
} // NAMESPACE

#endif // COVARIANCEEVOKEDSETUPWIDGET_H
