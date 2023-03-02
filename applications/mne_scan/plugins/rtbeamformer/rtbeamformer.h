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

#include <utils/generics/circularbuffer.h>

#include <fiff/fiff_evoked.h>

#include <mne/mne_beamformer_weights.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QFuture>
#include <QPointer>
#include <QSharedPointer>
#include <QFile>

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

namespace DISPLIB {
    class BeamformerSettingsView;
}

namespace MNELIB {
    class MNEForwardSolution;
    class MNEBeamformerWeights;
}

namespace FIFFLIB {
    class FiffInfo;
    class FiffInfoBase;
}

namespace INVERSELIB {
    class Beamformer;
}

namespace RTPROCESSINGLIB {
    class RtBfWeights;
}

namespace FSLIB {
    class AnnotationSet;
    class SurfaceSet;
}

namespace SCMEASLIB {
    class RealTimeMultiSampleArray;
    class RealTimeEvokedSet;
    class RealTimeEvokedCov;
    class RealTimeFwdSolution;
    class RealTimeSourceEstimate;
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

    friend class RtBeamformerSetupWidget;


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

    //=========================================================================================================
    /**
     * Inits widgets which are used to control this plugin, then emits them in form of a QList.
     */
    void initPluginControlWidgets();


    virtual void unload();
    virtual bool start();
    virtual bool stop();
    virtual SCSHAREDLIB::AbstractPlugin::PluginType getType() const;
    virtual QString getName() const;
    virtual QWidget* setupWidget();
    virtual QString getBuildInfo();



    //=========================================================================================================
    /**
     * Slot called when the fiff info is to be calculated.
     */
    bool calcFiffInfo();

    //=========================================================================================================
    /**
     * Slot to update the real time multi sample array data
     */
    void updateRTMSA(SCMEASLIB::Measurement::SPtr pMeasurement);

    //=========================================================================================================

    /**
     * Slot to update the fiff evoked
     *
     * @param[in] pMeasurement   The evoked to be appended.
     */
    void updateRTE(SCMEASLIB::Measurement::SPtr pMeasurement);

    //=========================================================================================================
    /**
     * Slot to update the fiff covariance
     */
    void updateRTCE(SCMEASLIB::Measurement::SPtr pMeasurement);

    //=========================================================================================================

    /**
     * Slot to update the real time forward solution
     */
    void updateRTFS(SCMEASLIB::Measurement::SPtr pMeasurement);

    //=========================================================================================================
    /**
     * Slot to update the beamformer weights
     *
     * @param[in] invOp    The beamformer weights to update.
     */
    void updateBFWeights(const MNELIB::MNEBeamformerWeights& bfWeights);

    //=========================================================================================================


protected:
    //=========================================================================================================
    /**
     * Slot called when the weight normalization method changed.
     *
     * @param[in] method        The new weight normalization.
     */
    void onWeightnormChanged(const QString &weightnorm);

    //=========================================================================================================
    /**
     * Slot called when the trigger type changed.
     *
     * @param[in] triggerType        The new trigger type.
     */
    void onTriggerTypeChanged(const QString& triggerType);

    //=========================================================================================================
    /**
     * Slot called when the time point changes.
     *
     * @param[in] iTimePointMs        The new time point in ms.
     */
    void onTimePointValueChanged(int iTimePointMs);

    //=========================================================================================================
    /**
     * AbstractAlgorithm functions
     */

    virtual void run();

    //=========================================================================================================

    QSharedPointer<SCSHAREDLIB::PluginInputData<SCMEASLIB::RealTimeMultiSampleArray> >      m_pRTMSAInput;              /**< The RealTimeMultiSampleArray input.*/
    QSharedPointer<SCSHAREDLIB::PluginInputData<SCMEASLIB::RealTimeEvokedSet> >             m_pRTESInput;               /**< The RealTimeEvoked input.*/
    QSharedPointer<SCSHAREDLIB::PluginInputData<SCMEASLIB::RealTimeEvokedCov> >             m_pRTCEInput;                /**< The RealTimeEvokedCov input.*/
    QSharedPointer<SCSHAREDLIB::PluginInputData<SCMEASLIB::RealTimeFwdSolution> >           m_pRTFSInput;               /**< The RealTimeFwdSolution input.*/
    QSharedPointer<SCSHAREDLIB::PluginOutputData<SCMEASLIB::RealTimeSourceEstimate> >       m_pRTSEOutput;              /**< The RealTimeSourceEstimate output.*/


    QSharedPointer<UTILSLIB::CircularBuffer_Matrix_double >                                 m_pCircularMatrixBuffer;    /**< Holds incoming RealTimeMultiSampleArray data.*/
    QSharedPointer<UTILSLIB::CircularBuffer<FIFFLIB::FiffEvoked> >                          m_pCircularEvokedBuffer;    /**< Holds incoming RealTimeMultiSampleArray data.*/

    QSharedPointer<RTPROCESSINGLIB::RtBfWeights>                                            m_pRtBfWeights;                 /**< Real-time beamformer weights. */



    QSharedPointer<MNELIB::MNEForwardSolution>                                              m_pFwd;                     /**< Forward solution. */
    QSharedPointer<FIFFLIB::FiffCov>                                                        m_pNoiseCov;                     /**< Noise Covariance Matrix. */
    QSharedPointer<FIFFLIB::FiffCov>                                                        m_pDataCov;                     /**< Data Covariance Matrix. */



    QSharedPointer<FIFFLIB::FiffInfoBase>                                                   m_pFiffInfoForward;         /**< Fiff information of the forward solution. */
    QSharedPointer<FIFFLIB::FiffInfo>                                                       m_pFiffInfo;                /**< Fiff information. */
    QSharedPointer<FIFFLIB::FiffInfo>                                                       m_pFiffInfoInput;           /**< Fiff information of the evoked. */

    QSharedPointer<FSLIB::AnnotationSet>                                                    m_pAnnotationSet;           /**< Annotation set. */
    QSharedPointer<FSLIB::SurfaceSet>                                                       m_pSurfaceSet;              /**< Surface set. */


    bool                            m_bRawInput;                /**< Flag whether a raw data input was received. */
    bool                            m_bEvokedInput;             /**< Flag whether an evoked input was received. */
    bool                            m_bUpdateBeamformer;       /**< Flag whether to update the beamformer object. */


    QMutex                          m_qMutex;                   /**< The mutex ensuring thread safety. */
    QFuture<void>                   m_future;                   /**< The future monitoring the clustering. */

    FIFFLIB::FiffEvoked             m_currentEvoked;
    FIFFLIB::FiffCoordTrans         m_mriHeadTrans;             /**< the Mri Head transformation. */

    qint32                          m_iNumAverages;             /**< The number of trials/averages to store. */
    qint32                          m_iTimePointSps;            /**< The time point to pick from the data in samples. */
    qint32                          m_iDownSample;              /**< Down sample factor. */


    QString                         m_sWeightnorm;               /**< The method for Weight normalization: "no" | "unitnoisegain" | "arraygain" | "nai". */

    QString                         m_sAvrType;                 /**< The average type. */
    QString                         m_sAtlasDir;                /**< File to Atlas. */
    QString                         m_sSurfaceDir;              /**< File to Surface. */
    QFile                           m_fMriHeadTrans;            /**< The Head - Mri transformation. */


    QStringList                     m_qListNoiseCovChNames;          /**< Noise Covariance channel names. */
    QStringList                     m_qListDataCovChNames;         /**< Data Covariance channel names. */  //HINT: new field here
    QStringList                     m_qListPickChannels;        /**< Channels to pick. */

    MNELIB::MNEBeamformerWeights      m_bfWeights;                    /**< The beamformer weights. */


signals:
    void responsibleTriggerTypesChanged(const QStringList& lResponsibleTriggerTypes);


};

}   // namespace

#endif // RTBEAMFORMER_H
