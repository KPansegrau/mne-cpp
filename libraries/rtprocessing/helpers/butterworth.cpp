//=============================================================================================================
/**
 * @file     butterworth.cpp
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
 * @brief    Definition of the Butterworth class
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================


#include "butterworth.h"

#include <fiff/fiff_info.h>

#include <complex>
#include <iostream>
#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>


//=============================================================================================================
// QT INCLUDES
//=============================================================================================================


//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/Core>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace RTPROCESSINGLIB;
using namespace Eigen;

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

Butterworth::Butterworth(int iType,
                         int iOrder,
                         double dCenterFreq,
                         double dBandwidth,
                         double dSFreq)

: m_iFilterType(iType)
, m_iFilterOrder(iOrder)
, m_dBandwidth(dBandwidth) //in Hz, not normed
, m_dCenterFreq(dCenterFreq) //in Hz, not normed
, m_dSFreq(dSFreq) //in Hz

{

    calculateButterworthCoeffs();
}

//=============================================================================================================

void Butterworth::createAnalogLowpassPrototype(int iOrder)
{

    RowVectorXcd vecPoles = RowVectorXcd::Zero(iOrder);

    //for loop calculates all poles of the analog prototype according to Euler's equation
    //prototype poles are ordered by descending absolut value of imaginary part
    int iCount{1};
    for(int i = 0,iIndex = 0; i < (iOrder/2); i++, iIndex+=2){
        double dPhi = iCount * (M_PI / (2*iOrder)) + (M_PI / 2);
        vecPoles[iIndex] = std::complex<double> (cos(dPhi), sin(dPhi));
        vecPoles[iIndex+1] = std::complex<double> (cos(dPhi), -sin(dPhi));
        iCount += 2;
    }
    if((iOrder%2) != 0){ //for odd filter orders add one pole on real axis
        vecPoles[iOrder-1] = std::complex<double> (-1,0);
    }

    m_vecPrototypePoles = vecPoles;

    //check for stability: all poles of the prototype have to be in the left half of s-plane
    int iNumPoles = m_vecPrototypePoles.size();
    for(int i = 0; i < iNumPoles; i++){
        if(m_vecPrototypePoles[i].real() > 0){
            qCritical() << "[Butterworth::createAnalogLowpassPrototype] Analog Prototype not stable due to poles outside of left half of s-plane. ";
            break;
        }
    }

    // for Debugging
    std::cout.precision(17);
    qDebug() << "[Butterworth::createAnalogLowpassPrototype] Prototype poles:";
    for(int i{0}; i < iOrder; i++){
       std::cout << "Prototype pole #" << i << " : " << m_vecPrototypePoles[i] << "\n"; //for debugging print prototype poles
    }

    //calculate prototype gain (should be 1)
    std::complex<double> dPoleProduct(1.0,0.0);
    for(int i{0}; i < iNumPoles; i++){
        dPoleProduct *= (-1.0 * m_vecPrototypePoles[i]);
    }
    m_dPrototypeGain = dPoleProduct.real();

    //for Debugging
    std::cout << "[Butterworth::createAnalogLowpassPrototype] Prototype gain: " << m_dPrototypeGain << '\n';
}

//=============================================================================================================

void Butterworth::convertPrototype2Lowpass()
{
    int iNumPoles = m_vecPrototypePoles.size();

    //for debugging: critical frequency (cutoff) of new analog LP is m_dOmegaLowPrewarped
    qDebug() << "[Butterworth::convertPrototype2Lowpass] Analog LP cutoff frequency in rad/s: " << m_dOmegaLowPrewarped;

    //calculate new analog LP poles (elementwise operation)
    //critical frequency (cutoff) of new analog LP is m_dOmegaLowPrewarped
    m_vecAnalogPoles = m_dOmegaLowPrewarped * m_vecPrototypePoles.array();


    //for Debugging
    qDebug() << "[Butterworth::convertPrototype2Lowpass] Analog LP poles:";
    for(int i{0}; i < iNumPoles; i++){
       std::cout << "Analog LP pole #" << i << " : " << m_vecAnalogPoles[i] << "\n";
    }

    //check for stability: all poles of the analog LP have to be in the left half of s-plane
    for(int i = 0; i < iNumPoles; i++){
        if(m_vecAnalogPoles[i].real() > 0){
            qCritical() << "[Butterworth::convertPrototype2Lowpass] Analog LP not stable due to poles outside of left half of s-plane. ";
            break;
        }
    }

    //Butterworth LP has only zeros at w = infinty so m_vecAnalogZeros remains empty

    //calculate new analog gain
    m_dAnalogGain = m_dPrototypeGain * pow(m_dOmegaLowPrewarped, iNumPoles);

    std::cout << "[Butterworth::convertPrototype2Lowpass] Analog LP gain: " << m_dAnalogGain << "\n";

}

//=============================================================================================================

void Butterworth::convertPrototype2Highpass()
{
    int iNumPoles = m_vecPrototypePoles.size();

    //for debugging: critical frequency (cutoff) of new analog HP is m_dOmegaHighPrewarped
    qDebug() << "[Butterworth::convertPrototype2Highpass] Analog HP cutoff frequency in rad/s: " << m_dOmegaHighPrewarped;

    //calculate new analog HP poles (elementwise operation)
    m_vecAnalogPoles = m_dOmegaHighPrewarped / m_vecPrototypePoles.array();


    //for Debugging
    qDebug() << "[Butterworth::convertPrototype2Highpass] Analog HP poles:";
    for(int i{0}; i < iNumPoles; i++){
       std::cout << "Analog HP pole #" << i << " : " << m_vecAnalogPoles[i] << "\n";
    }

    //check for stability: all poles of the analog HP have to be in the left half of s-plane
    for(int i = 0; i < iNumPoles; i++){
        if(m_vecAnalogPoles[i].real() > 0){
            qCritical() << "[Butterworth::convertPrototype2Highpass] Analog HP not stable due to poles outside of left half of s-plane. ";
            break;
        }
    }

    //calculate new analog HP zeros --> each prototype pole causes a zero at origin
    m_vecAnalogZeros = RowVectorXcd::Zero(iNumPoles);
    //for Debugging
    qDebug() << "[Butterworth::convertPrototype2Highpass] Analog HP zeros:";
    for(int i{0}; i < iNumPoles; i++){
       std::cout << "Analog HP zero #" << i << " : " << m_vecAnalogZeros[i] << "\n";
    }

    //calculate new analog gain
    std::complex<double> dPoleProduct(1.0,0.0);
    for(int i{0}; i < iNumPoles; i++){
        dPoleProduct *= (-1.0 / m_vecPrototypePoles[i]);
    }
    m_dAnalogGain = m_dPrototypeGain * dPoleProduct.real();
    qDebug() << "[Butterworth::convertPrototype2Highpass] Analog HP gain: " << m_dAnalogGain;

}

//=============================================================================================================

void Butterworth::convertPrototype2Bandpass()
{
    int iNumPoles = m_vecPrototypePoles.size();

    //calculate bandwidth and center frequency from prewarped frequencies (result is in rad/s)
    double dBandwidthPrewarped = m_dOmegaHighPrewarped - m_dOmegaLowPrewarped;
    double dCenterFreqPrewarped = sqrt(m_dOmegaHighPrewarped * m_dOmegaLowPrewarped);

    //for debugging: bandwidth and center frequency of the new analog BP
    std::cout.precision(17);
    std::cout << "[Butterworth::convertPrototype2Bandpass] Analog BP bandwidth in rad/s (after prewarping): " << dBandwidthPrewarped << "\n";
    std::cout << "[Butterworth::convertPrototype2Bandpass] Analog BP center frequency in rad/s (after prewarping): " << dCenterFreqPrewarped << "\n";

    //calculate new analog BP poles --> each prototype pole causes a pair of new poles
    m_vecAnalogPoles = RowVectorXcd::Zero(iNumPoles*2);
    for(int iPole{0},iCount{0}; iPole < iNumPoles; iPole++, iCount+=2){
        m_vecAnalogPoles[iCount] = 0.5 * (m_vecPrototypePoles[iPole] * dBandwidthPrewarped + sqrt(pow(m_vecPrototypePoles[iPole], 2) * pow(dBandwidthPrewarped, 2) - (4 * pow(dCenterFreqPrewarped,2))));
        m_vecAnalogPoles[iCount + 1] = 0.5 * (m_vecPrototypePoles[iPole] * dBandwidthPrewarped - sqrt(pow(m_vecPrototypePoles[iPole], 2) * pow(dBandwidthPrewarped, 2) - (4 * pow(dCenterFreqPrewarped,2))));
    }

    //we need the poles ordered by descending absolute value of imag part (so that conjugate complex pairs are stored in subsequent vector elements, descending order because this matches Matlabs order)
    std::sort(m_vecAnalogPoles.data(),m_vecAnalogPoles.data()+m_vecAnalogPoles.size(),
            [&](std::complex<double> a, std::complex<double> b){if(std::abs(a.imag()) == std::abs(b.imag())){ return a.imag() > b.imag();} return (std::abs(a.imag()) > std::abs(b.imag()));});


    //for Debugging
    qDebug() << "[Butterworth::convertPrototype2Bandpass] Analog BP poles:";
    for(int i{0}; i < m_vecAnalogPoles.size(); i++){
       std::cout << "Analog BP pole #" << i << " : " << m_vecAnalogPoles[i] << "\n";
    }

    //check for stability: all poles of the analog BP have to be in the left half of s-plane
    for(int i = 0; i < iNumPoles; i++){
        if(m_vecAnalogPoles[i].real() > 0){
            qCritical() << "[Butterworth::convertPrototype2Bandpass] Analog BP not stable due to poles outside of left half of s-plane. ";
            break;
        }
    }

    //calculate new analog BP zeros --> each prototype pole causes a zero at origin
    m_vecAnalogZeros = RowVectorXcd::Zero(iNumPoles);
    //for Debugging
    qDebug() << "[Butterworth::convertPrototype2Bandpass] Analog BP zeros:";
    for(int i{0}; i < iNumPoles; i++){
       std::cout << "Analog BP zero #" << i << " : " << m_vecAnalogZeros[i] << "\n";
    }

    //calculate new analog gain
    m_dAnalogGain = m_dPrototypeGain * pow(dBandwidthPrewarped, iNumPoles);
    std::cout << "[Butterworth::convertPrototype2Bandpass] Analog BP gain: " << m_dAnalogGain << "\n";

}

//=============================================================================================================

void Butterworth::convertPrototype2Bandstop()
{
    int iNumPoles = m_vecPrototypePoles.size();

    //calculate bandwidth and center frequency from prewarped frequencies (result is in rad/s)
    double dBandwidthPrewarped = m_dOmegaHighPrewarped - m_dOmegaLowPrewarped;
    double dCenterFreqPrewarped = sqrt(m_dOmegaHighPrewarped * m_dOmegaLowPrewarped);

    //for debugging: bandwidth and center frequency of the new analog BP
    std::cout.precision(17);
    std::cout << "[Butterworth::convertPrototype2Bandpass] Analog BS bandwidth in rad/s (after prewarping): " << dBandwidthPrewarped << "\n";
    std::cout << "[Butterworth::convertPrototype2Bandpass] Analog BS center frequency in rad/s (after prewarping): " << dCenterFreqPrewarped << "\n";

    //calculate new analog BS poles --> each prototype pole causes a pair of new poles
    m_vecAnalogPoles = RowVectorXcd::Zero(iNumPoles*2);
    for(int iPole{0}, iCount{0}; iPole < iNumPoles; iPole++, iCount+=2){
        m_vecAnalogPoles[iCount] = -0.5 * (-1 * dBandwidthPrewarped / m_vecPrototypePoles[iPole] + sqrt(pow(dBandwidthPrewarped, 2) / pow(m_vecPrototypePoles[iPole], 2) - (4 * pow(dCenterFreqPrewarped,2))));
        m_vecAnalogPoles[iCount + 1] = -0.5 * (-1 * dBandwidthPrewarped / m_vecPrototypePoles[iPole] - sqrt(pow(dBandwidthPrewarped, 2) / pow(m_vecPrototypePoles[iPole], 2) - (4 * pow(dCenterFreqPrewarped,2))));
    }

    //we need the poles ordered by descending absolute value of imag part (so that conjugate complex pairs are stored in subsequent vector elements)
    //this sorting does not match Matlabs order of analog poles but descending order is closer than ascending order
    std::sort(m_vecAnalogPoles.data(),m_vecAnalogPoles.data()+m_vecAnalogPoles.size(),
            [&](std::complex<double> a, std::complex<double> b){if(std::abs(a.imag()) == std::abs(b.imag())){ return a.imag() > b.imag();} return (std::abs(a.imag()) > std::abs(b.imag()));});

    //for Debugging
    qDebug() << "[Butterworth::convertPrototype2Bandstop] Analog BS poles:";
    for(int i{0}; i < m_vecAnalogPoles.size(); i++){
       std::cout << "Analog BS pole #" << i << " : " << m_vecAnalogPoles[i] << "\n";
    }

    //calculate new analog BS zeros --> each prototype pole causes a pair of imaginary zeros
    m_vecAnalogZeros = RowVectorXcd::Zero(iNumPoles*2);
    for(int iPole{0}, iCount{0}; iPole < iNumPoles; iPole++, iCount+=2){
        m_vecAnalogZeros[iCount] = std::complex<double> (0.0, dCenterFreqPrewarped);
        m_vecAnalogZeros[iCount + 1] = conj(m_vecAnalogZeros[iCount]);
    }


    //for Debugging
    qDebug() << "[Butterworth::convertPrototype2Bandstop] Analog BS zeros:";
    for(int i{0}; i < m_vecAnalogZeros.size(); i++){
       std::cout << "Analog BS zero #" << i << " : " << m_vecAnalogZeros[i] << "\n";
    }

    //check for stability: all poles of the analog BS have to be in the left half of s-plane
    for(int i = 0; i < iNumPoles; i++){
        if(m_vecAnalogPoles[i].real() > 0){
            qCritical() << "[Butterworth::convertPrototype2Bandstop] Analog BS not stable due to poles outside of left half of s-plane. ";
            break;
        }
    }

    //calculate new analog gain
    std::complex<double> dPoleProduct(1.0,0.0);
    for(int i{0}; i < iNumPoles; i++){
        dPoleProduct *= (-1.0 / m_vecPrototypePoles[i]);
    }
    m_dAnalogGain = m_dPrototypeGain * dPoleProduct.real();
    qDebug() << "[Butterworth::convertPrototype2Bandstop] Analog BS gain: " << m_dAnalogGain;

}

//=============================================================================================================

void Butterworth::bilinearTransform()
{   //use bilinear transform with s = K * (z - 1)/(z + 1) with K = 2fs

    double dPrewarpFactorK = 2 * m_dSFreq;

    //transform s-plane poles into z-plane poles
    int iNumPoles = m_vecAnalogPoles.size();
    m_vecDigitalPoles = RowVectorXcd::Zero(iNumPoles);
    for(int iPole{0}; iPole < iNumPoles; iPole++){
        m_vecDigitalPoles[iPole] = (1.0 + (m_vecAnalogPoles[iPole]/dPrewarpFactorK)) / (1.0 - (m_vecAnalogPoles[iPole]/dPrewarpFactorK));
    }

    //for Debugging
    std::cout.precision(17);
    qDebug() << "[Butterworth::bilinearTransform] Digital poles after bilinear transform:";
    for(int i{0}; i < iNumPoles; i++){
       std::cout << "Digital pole #" << i << " : " << m_vecDigitalPoles[i] << "\n";
    }

    //check for stability: all poles of the digital filter have to be inside the unit circle of z-plane
    for(int i = 0; i < iNumPoles; i++){
        if((m_vecDigitalPoles[i].real() > 1) && (m_vecDigitalPoles[i].imag() > 1) ){
            qCritical() << "[Butterworth::bilinearTransform] Digital filter not stable due to poles outside of z-plane unit circle.";
            break;
        }
    }

    //transform s-plane zeros into z-plane zeros
    //we need a second for loop because number of analog poles and zeros differs in case of LP and BP
    int iNumZeros = m_vecAnalogZeros.size();
    m_vecDigitalZeros = RowVectorXcd::Zero(std::max(iNumZeros,iNumPoles));
    for(int iZero{0}; iZero < iNumZeros; iZero++){
        m_vecDigitalZeros[iZero] = (1.0 + (m_vecAnalogZeros[iZero]/dPrewarpFactorK)) / (1.0 - (m_vecAnalogZeros[iZero]/dPrewarpFactorK));
    }


    //add further zeros at (-1,0) so that number of zeros equals number of poles (numerator and denominator with same order)
    if(iNumZeros < iNumPoles){
        for(int iZero{iNumZeros}; iZero < iNumPoles; iZero++){
            m_vecDigitalZeros[iZero] = std::complex<double> (-1.0,0.0);
        }
        iNumZeros = m_vecDigitalZeros.size();
    }


    //for Debugging
    qDebug() << "[Butterworth::bilinearTransform] Digital zeros after bilinear transform:";
    for(int i{0}; i < iNumZeros; i++){
       std::cout << "Digital zero #" << i << " : " << m_vecDigitalZeros[i] << "\n";
    }

    //check for stability: all zeros of the digital filter have to be inside the unit circle of z-plane
    for(int i = 0; i < iNumZeros; i++){
        if((m_vecDigitalZeros[i].real() > 1) && (m_vecDigitalZeros[i].imag() > 1) ){
            qCritical() << "[Butterworth::bilinearTransform] Digital filter not stable due to zeros outside of z-plane unit circle.";
            break;
        }
    }

    //transform s-plane gain into z-plane gain
    // we must use the numbers of analog zeros and poles here
    std::complex<double> dPoleProduct(1.0,0.0);
    for(int i{0}; i < iNumPoles; i++){
        dPoleProduct *= (dPrewarpFactorK - m_vecAnalogPoles[i]);
    }

    std::complex<double> dZeroProduct(1.0,0.0);
    for(int i{0}; i < (m_vecAnalogZeros.size()); i++){
        dZeroProduct *= (dPrewarpFactorK - m_vecAnalogZeros[i]);
    }

    std::complex<double> dTmpGain = m_dAnalogGain * dZeroProduct / dPoleProduct;
    m_dDigitalGain = dTmpGain.real();

    qDebug() << "[Butterworth::bilinearTransform] Digital gain after bilinear transform: " << m_dDigitalGain;


}


//=============================================================================================================

void Butterworth::splitIntoBiquads()
{

    //numbers of digital poles and zeros before splitting into stages
    int iNumPoles = m_vecDigitalPoles.size();
    int iNumZeros = m_vecDigitalZeros.size();

    //for Debugging
    qDebug() << "[Butterworth::splitIntoBiquads] Number of digital poles: " << iNumPoles;
    qDebug() << "[Butterworth::splitIntoBiquads] Number of digital zeros: " << iNumZeros;

    //get number of necessary biquads
    bool bDifferentLastStage = ((iNumPoles % 2) != 0);
    qDebug() << "[Butterworth::splitIntoBiquads] bDifferentLastStage: " << bDifferentLastStage;
    if(bDifferentLastStage){
        // odd number of poles and zeros --> last stage with only one pole and one zero (last stage is first order)
        m_iNumStages = (iNumPoles / 2) + 1;
    }else{
        //even number of poles and zeros --> all stages with two poles and two zeros (all stages are biquads)
        m_iNumStages = iNumPoles / 2;
    }
    qDebug() << "[Butterworth::splitIntoBiquads] Number of necessary Stages: " << m_iNumStages;

    //TODO: implement grouping of poles and zeros (decide which poles and zeros are combined for each biquad) to minimize errors due to numberical precision issues

    int iNumCoeffs = m_iNumStages * 3;
    qDebug() << "[Butterworth::splitIntoBiquads] Number of recursion coefficients a and b each: " << iNumCoeffs;


    //calculate biquad coeffs and store in m_vecRecCoeffsA, and m_vecRecCoeffsB
    //first three entries of m_vecRecCoeffsA are a0, a1 and a2 of the first biquad, next three entries contain a0, a1, a2 of the second biquad etc.
    m_vecRecCoeffsA = RowVectorXd::Zero(iNumCoeffs);
    m_vecRecCoeffsB = RowVectorXd::Zero(iNumCoeffs);
    for(int iCoeff{0}, i{0}; iCoeff < iNumCoeffs; iCoeff += 3, i += 2){

        if((iCoeff < iNumCoeffs-3) || (!bDifferentLastStage)){

        //calculate a's from poles
        m_vecRecCoeffsA[iCoeff] = 1.0; //a0
        m_vecRecCoeffsA[iCoeff + 1] = (- (m_vecDigitalPoles[i] + m_vecDigitalPoles[i + 1])).real(); //a1
        m_vecRecCoeffsA[iCoeff + 2] = (m_vecDigitalPoles[i] * m_vecDigitalPoles[i + 1]).real(); //a2

        //calculate b's from zeros
        m_vecRecCoeffsB[iCoeff] = 1.0; //b0
        m_vecRecCoeffsB[iCoeff + 1] = (- (m_vecDigitalZeros[i] + m_vecDigitalZeros[i + 1])).real(); //b1
        m_vecRecCoeffsB[iCoeff + 2] = (m_vecDigitalZeros[i] * m_vecDigitalZeros[i +1]).real(); //b2
        }

        //if number of poles is odd, calculate coefficients for last stage with single pole and zero differently
        if((iCoeff == (iNumCoeffs - 3)) && bDifferentLastStage){
            //calculate a's from the single pole
            m_vecRecCoeffsA[iCoeff] = 1.0; //a0
            m_vecRecCoeffsA[iCoeff + 1] = - m_vecDigitalPoles[i].real(); //a1
            m_vecRecCoeffsA[iCoeff + 2] = 0.0; //a2

            //calculate b's from the single zero
            m_vecRecCoeffsB[iCoeff] = 1.0; //b0
            m_vecRecCoeffsB[iCoeff + 1] = (- m_vecDigitalZeros[i]).real(); //b1
            m_vecRecCoeffsB[iCoeff + 2] = 0.0; //b2
        }
   }


    //apply overall gain to the first biquad by multiplying the coefficients b0, b1, b2 of the first biquad by the digital gain
    for(int i{0}; i < 3; i++){
        m_vecRecCoeffsB[i] *= m_dDigitalGain;
    }

    //for Debugging
    qDebug() << "[Butterworth::splitIntoBiquads] Recursion coefficients a0, a1, a2 for each biquad:";
    for(int i{0}; i < iNumCoeffs; i++){
       std::cout << "A Coefficient #" << i << " : " << m_vecRecCoeffsA[i] << "\n";
    }
    qDebug() << "[Butterworth::splitIntoBiquads] Recursion coefficients b0, b1, b2 for each biquad:";
    for(int i{0}; i < iNumCoeffs; i++){
       std::cout << "B Coefficient #" << i << " : " << m_vecRecCoeffsB[i] << "\n";
    }


}

//=============================================================================================================

bool Butterworth::calculateButterworthCoeffs()
{

    //set cutoff frequencies depending on selected filter type
    // for LP and HP one of the cutoff frequencies have to be zero
    //center frequency and bandwidth are input parameters in Hz (not normed)
    switch(m_iFilterType) {
        case 0:
            m_dHighCutoff = 0;
            m_dLowCutoff = m_dCenterFreq;
        break;

        case 1:
            m_dHighCutoff = m_dCenterFreq;
            m_dLowCutoff = 0;
        break;

        case 2:
            m_dHighCutoff = m_dCenterFreq + m_dBandwidth/2;
            m_dLowCutoff = m_dCenterFreq - m_dBandwidth/2;
        case 3:
            m_dHighCutoff = m_dCenterFreq + m_dBandwidth/2;
            m_dLowCutoff = m_dCenterFreq - m_dBandwidth/2;
        break;
    }

    //for Debugging
    std::cout.precision(17);
    std::cout << "[Butterworth::calculateButterworthCoeffs] Lower cutoff frequency in Hz: " << m_dLowCutoff << "\n";
    std::cout << "[Butterworth::calculateButterworthCoeffs] Higher cutoff frequency in Hz: " << m_dHighCutoff << "\n";
    std::cout << "[Butterworth::calculateButterworthCoeffs] Sampling frequency in Hz: " << m_dSFreq << "\n";

    //convert cutoff frequencies from Hz to rad/s and perform frequency prewarping
    //reason: bilinear transform results in non linear compression especially at high frequencies which can be compensated by this prewarping step
    m_dOmegaLowPrewarped = 2 * m_dSFreq * tan(m_dLowCutoff * M_PI / m_dSFreq);
    m_dOmegaHighPrewarped = 2 * m_dSFreq * tan(m_dHighCutoff * M_PI / m_dSFreq);
    std::cout << "[Butterworth::calculateButterworthCoeffs] Prewarped lower cutoff frequency in rad/s: from " << m_dLowCutoff*2*M_PI << " to " << m_dOmegaLowPrewarped << "\n";
    std::cout << "[Butterworth::calculateButterworthCoeffs] Prewarped higher cutoff frequency in rad/s: from " << m_dHighCutoff*2*M_PI << " to " << m_dOmegaHighPrewarped << "\n";

    //create an anlog lowpass prototype of desired order (its cutoff frequency in rad/s is one)
    //the prototype has iOrder poles, no zeros and gain of one
    createAnalogLowpassPrototype(m_iFilterOrder);

    //convert prototype into desired filter type and adjust cutoff frequencies
    //new poles and zeros are in m_vecAnalogPoles and m_vecAnalogZeros
    //m_dAnalogGain is set to the adjusted gain after conversion
    switch(m_iFilterType){
        case 0: {   //desired type is LP  (iOrder poles, no zeros)
            convertPrototype2Lowpass();
            break;
        }
        case 1: {   //desired type is HP (iOrder poles, iOrder zeros)
            convertPrototype2Highpass();
            break;
        }
        case 2: {   //desired type is BP (2*iOrder poles, iOrder zeros)
            convertPrototype2Bandpass();
            break;
        }
        case 3: {   //desired type is NOTCH (2*iOrder poles, 2*iOrder zeros)
            convertPrototype2Bandstop();
            break;
        }
        default: {  //desired type is unknown
            qCritical() << "[Butterworth::calculateButterworthCoeffs] Desired filter type is unknown.";
            return false;
        }

    }

    //perform bilinear transform to get a digital filter from the analog one
    bilinearTransform();

    //convert digital poles, zeros and gain into second order sections (biquads)
    splitIntoBiquads();

    return true;

}

//=============================================================================================================

RowVectorXcd Butterworth::getPrototypePoles() const
{
    return m_vecPrototypePoles;
}

//=============================================================================================================

double Butterworth::getPrototypeGain() const
{
    return m_dPrototypeGain;
}

//=============================================================================================================

RowVectorXcd Butterworth::getAnalogPoles() const
{
    return m_vecAnalogPoles;
}

//=============================================================================================================

RowVectorXcd Butterworth::getAnalogZeros() const
{
    return m_vecAnalogZeros;
}

//=============================================================================================================

double Butterworth::getAnalogGain() const
{
    return m_dAnalogGain;
}

//=============================================================================================================

RowVectorXcd Butterworth::getDigitalPoles() const
{
    return m_vecDigitalPoles;
}


//=============================================================================================================

RowVectorXcd Butterworth::getDigitalZeros() const
{
    return m_vecDigitalZeros;
}

//=============================================================================================================

double Butterworth::getDigitalGain() const
{
    return m_dDigitalGain;
}

//=============================================================================================================

Eigen::RowVectorXd Butterworth::getRecCoeffsA() const
{
    return m_vecRecCoeffsA;
}

//=============================================================================================================
