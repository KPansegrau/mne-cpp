//=============================================================================================================
/**
 * @file     covarianceevoked.h
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
 * @brief    Contains the declaration of the CovarianceEvoked class.
 *
 */

#ifndef COVARIANCEEVOKED_H
#define COVARIANCEEVOKED_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "covarianceevoked_global.h"

#include <scShared/Plugins/abstractalgorithm.h>
#include <utils/generics/circularbuffer.h>

#include <fiff/fiff_evoked.h>


//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QFuture>
#include <QPointer>
#include <QSharedPointer>

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

namespace FIFFLIB {
    class FiffCov;
    class FiffInfo;
}

namespace RTPROCESSINGLIB {
    class RtCov;
}

namespace SCMEASLIB {
    class RealTimeEvokedSet;
    class RealTimeEvokedCov;
}

//=============================================================================================================
// DEFINE NAMESPACE COVARIANCEPLUGIN
//=============================================================================================================

namespace COVARIANCEEVOKEDPLUGIN
{

//=============================================================================================================
// COVARIANCEEVOKEDPLUGIN FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * DECLARE CLASS CovarianceEvoked
 *
 * @brief The CovariancEvoked class provides a CovarianceEvoked algorithm structure.
 */
class COVARIANCEEVOKEDSHARED_EXPORT CovarianceEvoked : public SCSHAREDLIB::AbstractAlgorithm
{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "scsharedlib/1.0" FILE "covarianceevoked.json") //New Qt5 Plugin system replaces Q_EXPORT_PLUGIN2 macro
    // Use the Q_INTERFACES() macro to tell Qt's meta-object system about the interfaces
    Q_INTERFACES(SCSHAREDLIB::AbstractAlgorithm)

    friend class CovarianceEvoked;

public:
    //=========================================================================================================
    /**
     * Constructs a CovarianceEvoked.
     */
    CovarianceEvoked();

    //=========================================================================================================
    /**
     * Destroys the Covariance.
     */
    ~CovarianceEvoked();

    //=========================================================================================================
    /**
     * Initialise input and output connectors.
     */
    virtual void init();

    //=========================================================================================================
    /**
     * Inits widgets which are used to control this plugin, then emits them in form of a QList.
     */
    void initPluginControlWidgets();

    //=========================================================================================================
    /**
     * Is called when plugin is detached of the stage. Can be used to safe settings.
     */
    virtual void unload();

    //=========================================================================================================

    //AbstractAlgorithm methods

    virtual QSharedPointer<SCSHAREDLIB::AbstractPlugin> clone() const;

    virtual bool start();
    virtual bool stop();

    virtual SCSHAREDLIB::AbstractPlugin::PluginType getType() const;
    virtual QString getName() const;

    virtual QWidget* setupWidget();

    virtual QString getBuildInfo();

    void updateRTE(SCMEASLIB::Measurement::SPtr pMeasurement);

    //TODO: maybe we can delete this method since it is not used (check and delete)
    void showCovarianceEvokedWidget();

    void changeSamples(qint32 samples);


protected:

    //=========================================================================================================
    /**
     * Slot called when the trigger type changed.
     *
     * @param[in] triggerType        The new trigger type.
     */
    void onTriggerTypeChanged(const QString& triggerType);

    //=========================================================================================================


    virtual void run();

private:
    //HINT: this variables are copied from rtcmne.h
    //TODO: check whether we need all
    QSharedPointer<SCSHAREDLIB::PluginInputData<SCMEASLIB::RealTimeEvokedSet> >             m_pCovarianceEvokedInput;               /**< The CovarianceEvoked input.*/
    QSharedPointer<SCSHAREDLIB::PluginOutputData<SCMEASLIB::RealTimeEvokedCov> >            m_pCovarianceEvokedOutput;     /**< The CovarianceEvoked output.*/

    QSharedPointer<FIFFLIB::FiffInfo>                                                       m_pFiffInfoInput;           /**< Fiff information of the evoked. */

    QSharedPointer<UTILSLIB::CircularBuffer<FIFFLIB::FiffEvoked> > m_pCircularEvokedBuffer;    /**< Holds incoming RealTimeMultiSampleArray data.*/

    QMutex                          m_qMutex;                   /**< The mutex ensuring thread safety. */


    qint32      m_iEstimationSamples; //number of samples used for estimation of covariance matrices

    QString                         m_sAvrType;                 /**< The average type. */

    quint32 m_iNumPreStimSamples;   /**< The number of pre stimulative samples. */

signals:
    void responsibleTriggerTypesChanged(const QStringList& lResponsibleTriggerTypes);

};

} // NAMESPACE

#endif // COVARIANCEEVOKED_H
