//=============================================================================================================
/**
 * @file     realtimeevokedcov.h
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de>;
 *
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
 * @brief    Contains the declaration of the RealTimeEvokedCov class.
 *
 */
#ifndef REALTIMEEVOKEDCOV_H
#define REALTIMEEVOKEDCOV_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "scmeas_global.h"
#include "measurement.h"
#include "realtimesamplearraychinfo.h"

#include <fiff/fiff_cov.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSharedPointer>
#include <QVector>
#include <QList>
#include <QColor>
#include <QMutex>
#include <QMutexLocker>

//=============================================================================================================
// DEFINE NAMESPACE SCMEASLIB
//=============================================================================================================

namespace SCMEASLIB
{

//=========================================================================================================
/**
 * DECLARE CLASS RealTimeEvokedCov
 *
 * @brief The RealTimeEvokedCov class provides a container for real-time covariance estimations from evoked input data.
 */
class SCMEASSHARED_EXPORT RealTimeEvokedCov : public Measurement
{
    Q_OBJECT

public:
    typedef QSharedPointer<RealTimeEvokedCov> SPtr;               /**< Shared pointer type for RealTimeEvokedCov. */
    typedef QSharedPointer<const RealTimeEvokedCov> ConstSPtr;    /**< Const shared pointer type for RealTimeEvokedCov. */

    //=========================================================================================================
    /**
     * Constructs a RealTimeEvokedCov.
     */
    explicit RealTimeEvokedCov(QObject *parent = 0);

    //=========================================================================================================
    /**
     * Destroys the RealTimeEvokedCov.
     */
    virtual ~RealTimeEvokedCov();

    //=========================================================================================================
    /**
     * Set the fiff info
     *
     * @param[in] pFiffInfo     the new fiff info.
     */
    void setFiffInfo(QSharedPointer<FIFFLIB::FiffInfo> pFiffInfo);

    //=========================================================================================================
    /**
     * Get the fiff info
     *
     * @return     the current fiff info.
     */
    QSharedPointer<FIFFLIB::FiffInfo> getFiffInfo();

    //=========================================================================================================
    /**
     * New pair of covariances to distribute
     *
     * @param[in] noiseCov     the noise covariance which should be distributed.
     * @param[in] dataCov      the data covariance which should be distributed.
     */
    virtual void setValue(const FIFFLIB::FiffCov& noiseCov, const FIFFLIB::FiffCov& dataCov);

    //=========================================================================================================
    /**
     * Returns the current pair of covariance matrices.
     * This method is inherited by Measurement.
     *
     * @return the last attached pair of covariance matrices.
     */
    virtual QSharedPointer<QPair<FIFFLIB::FiffCov,FIFFLIB::FiffCov> >& getValue();

    //=========================================================================================================
    /**
     * Returns whether RealTimeEvokedCov contains values
     *
     * @return whether RealTimeEvokedCov contains values.
     */
    inline bool isInitialized() const;

private:

    mutable QMutex                                              m_qMutex;           /**< Mutex to ensure thread safety. */

    QSharedPointer<QPair<FIFFLIB::FiffCov,FIFFLIB::FiffCov> >   m_pFiffCovPair;     /**< Pair of Noise Covariance and Data Covariance data set. */

    FIFFLIB::FiffInfo::SPtr                                     m_pFiffInfo;        /**< The Fiff Info. */

    bool                                                        m_bInitialized;     /**< If values are stored.*/
};

//=============================================================================================================
// INLINE DEFINITIONS
//=============================================================================================================

inline bool RealTimeEvokedCov::isInitialized() const
{
    QMutexLocker locker(&m_qMutex);
    return m_bInitialized;
}
} // NAMESPACE

Q_DECLARE_METATYPE(SCMEASLIB::RealTimeEvokedCov::SPtr)

#endif // REALTIMEEVOKEDCOV_H

