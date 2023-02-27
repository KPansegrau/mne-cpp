//=============================================================================================================
/**
 * @file     rtbeamformersetupwidget.h
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
 * @brief    Contains the declaration of the RtBeamformerSetupWidget class.
 *
 */

#ifndef RTBEAMFORMERSETUPWIDGET_H
#define RTBEAMFORMERSETUPWIDGET_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "ui_rtbeamformersetupwidget.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtWidgets>

//=============================================================================================================
// DEFINE NAMESPACE RTCMNEPLUGIN
//=============================================================================================================

namespace RTBEAMFORMERPLUGIN
{

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

class RtBeamformer;

//=============================================================================================================
/**
 * DECLARE CLASS DummySetupWidget
 *
 * @brief The DummySetupWidget class provides the DummyToolbox configuration window.
 */
class RtBeamformerSetupWidget : public QWidget
{
    Q_OBJECT

public:

    //=========================================================================================================
    /**
     * Constructs a RtBeamformerSetupWidget which is a child of parent.
     *
     * @param[in] toolbox a pointer to the corresponding MNEToolbox.
     * @param[in] parent pointer to parent widget; If parent is 0, the new RtBeamformerSetupWidget becomes a window. If parent is another widget, RtBeamformerSetupWidget becomes a child window inside parent. RtBeamformerSetupWidget is deleted when its parent is deleted.
     */
    RtBeamformerSetupWidget(RtBeamformer* toolbox, QWidget *parent = 0);

    //=========================================================================================================
    /**
     * Destroys the RtBeamformerSetupWidget.
     * All RtBeamformerSetupWidget's children are deleted first. The application exits if RtBeamformerSetupWidget is the main widget.
     */
    ~RtBeamformerSetupWidget();

private:

    //=========================================================================================================
    /**
     * Shows atlas selection dialog
     */
    void showAtlasDirDialog();

    //=========================================================================================================
    /**
     * Shows surface selection dialog
     */
    void showSurfaceDirDialog();

    //=========================================================================================================
    /**
     * Shows transformation selection dialog
     */
    void showMriHeadFileDialog();

    RtBeamformer* m_pBeamformer;

    Ui::RtBeamformerSetupWidgetClass ui;   /**< Holds the user interface for the RtBeamformerSetupWidgetClass.*/
};
} // NAMESPACE

#endif // RTBEAMFORMERSETUPWIDGET_H

