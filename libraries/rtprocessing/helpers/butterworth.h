//=============================================================================================================

//TODO edit documentation here

/**
 * @file     butterworth.h
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de>
 * @since    0.1.9
 * @date     April, 2022
 *
 * @section  LICENSE
 *
 * Copyright (C) 2022, Kerstin Pansegrau. All rights reserved.
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
 * @brief    Declaration of the Butterworth class
 *
 */


#ifndef BUTTERWORTH_H
#define BUTTERWORTH_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../rtprocessing_global.h"

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/Core>

//=============================================================================================================
// DEFINE NAMESPACE RTPROCESSINGLIB
//=============================================================================================================

namespace RTPROCESSINGLIB
{

//=============================================================================================================
/**
 * Creates an IIR butterworth filter.
 *
 * @brief
 */
class RTPROCESINGSHARED_EXPORT Butterworth
{
public:
    enum TPassType {LPF, HPF, BPF, NOTCH };

    //=========================================================================================================
    /**
     * Constructs a default Butterworth object.
     *
     */
    Butterworth();

    //=========================================================================================================
    /**
     * Constructs a Butterworth object.
     *
     * @param[in] iType             Type of the filter: LPF, HPF, BPF, NOTCH (from enum).
     * @param[in] iOrder            The filter order.
     * @param[in] dCenterFreq       The center frequency of the filter in Hz, not normed.
     * @param[in] dBandwidth        The bandwidth of the filter in Hz, not normed.
     * @param[in] dSFreq            The sampling frequency in Hz.
     *
     *
     **/
    Butterworth(int iType,
                int iOrder,
                double dCenterFreq,
                double dBandwidth,
                double dSFreq);

    //=========================================================================================================
    /**
     * Creates an analog lowpass Butterworth filter prototype with filter order iOrder and prototype gain of 1 (stored in m_dPrototypeGain).
     * The complex poles of the prototype are stored as conjugate complex pairs in m_vecPrototypePoles
     * with regard to ascending absolute real part.
     * The prototype has no zeros.
     *
     *
     * @param[in]   iOrder          The filter order.
     *
     *
     * References:  [1] MATLAB function butterap(filterOrder)
     */
    void createAnalogLowpassPrototype(int iOrder);

    //=========================================================================================================
    /**
     * Converts the analog LP prototype (cutoff = 1 rad/s) into an analog LP with desired cutoff frequency.
     * The poles of the analog LP are stored in m_vecAnalogPoles, analog LP has no zeros.
     * The new analog gain is stored in m_dAnalogGain.
     *
     */

    void convertPrototype2Lowpass();

    //=========================================================================================================
    /**
     * Converts the analog LP prototype (cutoff = 1 rad/s) into an analog HP with desired cutoff frequency.
     * The poles of the analog HP are stored in m_vecAnalogPoles, resulting analog zeros are stored in m_vecAnalogZeros.
     * The new analog gain is stored in m_dAnalogGain.
     *
     */

    void convertPrototype2Highpass();

    //=========================================================================================================
    /**
     * Converts the analog LP prototype (cutoff = 1 rad/s) into an analog BP with desired cutoff frequency.
     * The poles of the analog BP are stored in m_vecAnalogPoles, resulting analog zeros are stored in m_vecAnalogZeros.
     * The new analog gain is stored in m_dAnalogGain.
     */

    void convertPrototype2Bandpass();

    //=========================================================================================================
    /**
     * Converts the analog LP prototype (cutoff = 1 rad/s) into an analog BS with desired cutoff frequency.
     * The poles of the analog BS (NOTCH) are stored in m_vecAnalogPoles, resulting analog zeros are stored in m_vecAnalogZeros.
     * The new analog gain is stored in m_dAnalogGain.
     */

    void convertPrototype2Bandstop();

    //=========================================================================================================
    /**
     * Converts the analog filter into a digital one by using the bilinear transfrom.
     * Digital poles, zeros and gain are stored in m_vecDigitalPoles, m_vecDigitalZeros and m_dDigitalGain, respectively.
     *
     * References:  [1] MATLAB function bilinear() https://ch.mathworks.com/help/signal/ref/bilinear.html
     *
     */

    void bilinearTransform();

    //=========================================================================================================
    /**
     * Converts the digital filter transfer function H[z] in pole zero gain form into second order biquads.
     *
     *
     */

    void splitIntoBiquads();

    //=========================================================================================================
    /**
     * Calculates the Butterworth Recursion Coefficients.
     */
    bool calculateButterworthCoeffs();

    //=========================================================================================================


    Eigen::RowVectorXcd getPrototypePoles() const;
    double getPrototypeGain() const;

    Eigen::RowVectorXcd getAnalogPoles() const;
    Eigen::RowVectorXcd getAnalogZeros() const;
    double getAnalogGain() const;

    Eigen::RowVectorXcd getDigitalPoles() const;
    Eigen::RowVectorXcd getDigitalZeros() const;
    double getDigitalGain() const;

    Eigen::RowVectorXd getRecCoeffsA() const;
    Eigen::RowVectorXd getRecCoeffsB() const;

    //=========================================================================================================

    Eigen::RowVectorXd  m_vecRecCoeffsA;    /**< the recursion coefficients of the denominator of H[z] (a's) in second order section form. */
    Eigen::RowVectorXd  m_vecRecCoeffsB;    /**< the recursion coefficients of the numerator of H[z] (b's) in second order section form. */


    //=========================================================================================================
private:
    //=========================================================================================================

    int     m_iFilterType;                  /**< represents the type of the filter instance.*/
    int     m_iFilterOrder;                 /**< represents the order of the filter instance. */
    double  m_dBandwidth;                   /**< the bandwidth of the filter in Hz, not normed. */
    double  m_dCenterFreq;                  /**< the center frequency of the filter in Hz, not normed. */

    double  m_dLowCutoff;                   /**< the lower cutoff frequency in Hz, not normed. */
    double  m_dHighCutoff;                  /**< the higher cutoff frequency in Hz, not normed. */
    double  m_dSFreq;                       /**< the sampling frequency in Hz. */

    double m_dOmegaLowPrewarped;            /**< the prewarped lower cutoff frequency in rad/s */
    double m_dOmegaHighPrewarped;           /**< the prewarped higher cutoff frequency in rad/s */
    double m_dPrototypeGain;                /**< the gain of the analog LP prototype. */
    double m_dAnalogGain;                   /**< the gain of the analog filter after type conversion. */
    double m_dDigitalGain;                  /**< the gain of the digital filter after bilinear transform. */

    Eigen::RowVectorXcd m_vecPrototypePoles;/**< the poles of the analog LP prototype. */
    Eigen::RowVectorXcd m_vecAnalogPoles;   /**< the poles of the analog filter. */
    Eigen::RowVectorXcd m_vecAnalogZeros;   /**< the zeros of the analog filter. */
    Eigen::RowVectorXcd m_vecDigitalPoles;  /**< the poles of the digital filter. */
    Eigen::RowVectorXcd m_vecDigitalZeros;  /**< the zeros of the digital filter. */

    int m_iNumStages;                       /**< the number of necessary stages (second order sections (+ one first order section for odd filter orders)) for the digital filter. For odd filter orders the first order section is counted as well. */



};
} // NAMESPACE RTPROCESSINGLIB

#endif // BUTTERWORTH_H
