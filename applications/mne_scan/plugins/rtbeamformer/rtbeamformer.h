//=============================================================================================================
/**
 * @file     rtbeamformer.h
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de>
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
 * @brief    Contains the declaration of the RtBeamformer class.
 *
 */


#ifndef RTBEAMFORMER_H
#define RTBEAMFORMER_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "rtbeamformer_global.h"

#include <scShared/Plugins/abstractalgorithm.h>

#include "scMeas/realtimemultisamplearray.h"
#include <utils/generics/circularbuffer.h>

#include <fiff/fiff_evoked.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QFuture>
#include <QPointer>
#include <QSharedPointer>

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

namespace SCMEASLIB {
    class RealTimeMultiSampleArray;
    class RealTimeSourceEstimate;
}

namespace MNELIB {
    class MNEForwardSolution;
}

namespace FIFFLIB {
    class FiffInfo;
}


//=============================================================================================================
// DEFINE NAMESPACE RTBEAMFORMERPLUGIN
//=============================================================================================================

namespace RTBEAMFORMERPLUGIN
{

//=============================================================================================================
// RTBEAMFORMERPLUGIN FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * DECLARE CLASS RtBeamformer
 *
 * @brief The RtBeamformer class provides a plugin for estimating source localization in real-time using a beamformer.
 *
 */

class RTBEAMFORMERSHARED_EXPORT RtBeamformer : public SCSHAREDLIB::AbstractAlgorithm
{

    Q_OBJECT
    Q_PLUGIN_METADATA(IID "scsharedlib/1.0" FILE "rtbeamformer.json") //New Qt5 Plugin system replaces Q_EXPORT_PLUGIN2 macro
    // Use the Q_INTERFACES() macro to tell Qt's meta-object system about the interfaces
    Q_INTERFACES(SCSHAREDLIB::AbstractAlgorithm)


public:

    //=========================================================================================================
    /**
     * Constructs a RtBeamformer.
     */
    RtBeamformer();

    //=========================================================================================================
    /**
     * Destroys the RtcMne.
     */
    ~RtBeamformer();

    //=========================================================================================================
    /**
     * AbstractAlgorithm functions
     */

    virtual QSharedPointer<SCSHAREDLIB::AbstractPlugin> clone() const;
    virtual void init();
    virtual void unload();
    virtual bool start();
    virtual bool stop();
    virtual SCSHAREDLIB::AbstractPlugin::PluginType getType() const;
    virtual QString getName() const;
    virtual QWidget* setupWidget();
    virtual QString getBuildInfo();


    //=========================================================================================================
    /**
     * Inits widgets which are used to control this plugin, then emits them in form of a QList.
     */
    void initPluginControlWidgets();

    //=========================================================================================================
    /**
     * Slot to update the real time multi sample array data
     */
    void updateRTMSA(SCMEASLIB::Measurement::SPtr pMeasurement);

    //=========================================================================================================


protected:

    //=========================================================================================================
    /**
     * AbstractAlgorithm functions
     */

    virtual void run();

    //=========================================================================================================

    QSharedPointer<SCSHAREDLIB::PluginInputData<SCMEASLIB::RealTimeMultiSampleArray> >      m_pRTMSAInput;              /**< The RealTimeMultiSampleArray input.*/

    QSharedPointer<UTILSLIB::CircularBuffer_Matrix_double >                                 m_pCircularMatrixBuffer;    /**< Holds incoming RealTimeMultiSampleArray data.*/


    QSharedPointer<SCSHAREDLIB::PluginOutputData<SCMEASLIB::RealTimeSourceEstimate> >       m_pRTSEOutput;              /**< The RealTimeSourceEstimate output.*/

    QSharedPointer<MNELIB::MNEForwardSolution>                                              m_pFwd;                     /**< Forward solution. */




    QSharedPointer<FIFFLIB::FiffInfo>                                                       m_pFiffInfo;                /**< Fiff information. */
    QSharedPointer<FIFFLIB::FiffInfo>                                                       m_pFiffInfoInput;           /**< Fiff information of the evoked. */

    bool                            m_bRawInput;                /**< Flag whether a raw data input was received. */


    QMutex                          m_qMutex;                   /**< The mutex ensuring thread safety. */

    qint32                          m_iNumAverages;             /**< The number of trials/averages to store. */


};

}   // namespace

#endif // RTBEAMFORMER_H
