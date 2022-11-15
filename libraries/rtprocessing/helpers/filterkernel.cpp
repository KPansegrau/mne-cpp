
// TODO edit docu here

//=============================================================================================================
/**
 * @file     filterkernel.cpp
 * @author   Lorenz Esch <lesch@mgh.harvard.edu>;
 *           Ruben Doerfel <Ruben.Doerfel@tu-ilmenau.de>;
 *           Christoph Dinh <chdinh@nmr.mgh.harvard.edu>
 * @since    0.1.0
 * @date     February, 2014
 *
 * @section  LICENSE
 *
 * Copyright (C) 2014, Lorenz Esch, Ruben Doerfel, Christoph Dinh. All rights reserved.
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
 * @brief    Contains all FilterKernel.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "filterkernel.h"

#include <utils/mnemath.h>

#include "parksmcclellan.h"
#include "cosinefilter.h"
#include "butterworth.h"

#include <iostream>

#define _USE_MATH_DEFINES
#include <math.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QDebug>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/SparseCore>
//#ifndef EIGEN_FFTW_DEFAULT
//#define EIGEN_FFTW_DEFAULT
//#endif
#include <unsupported/Eigen/FFT>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace RTPROCESSINGLIB;
using namespace Eigen;
using namespace UTILSLIB;



//=============================================================================================================
// INIT STATIC MEMBERS
//=============================================================================================================

QVector<RTPROCESSINGLIB::FilterParameter> FilterKernel::m_designMethods ({
    FilterParameter(QString("Cosine"), QString("A cosine filter")),
    FilterParameter(QString("Tschebyscheff"), QString("A tschebyscheff filter")),
    FilterParameter(QString("Butterworth"), QString("An IIR butterworth filter"))
//    FilterParameter(QString("External"), QString("An external filter"))
});
QVector<RTPROCESSINGLIB::FilterParameter> FilterKernel::m_filterTypes ({
    FilterParameter(QString("LPF"), QString("An LPF filter")),
    FilterParameter(QString("HPF"), QString("An HPF filter")),
    FilterParameter(QString("BPF"), QString("A BPF filter")),
    FilterParameter(QString("NOTCH"), QString("A NOTCH filter")),
    FilterParameter(QString("UNKNOWN"), QString("An UNKNOWN filter"))
});

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

FilterKernel::FilterKernel()
: m_sFreq(1000)
, m_dCenterFreq(0.5)
, m_dBandwidth(0.1)
, m_dParksWidth(0.1)
, m_dLowpassFreq(40)
, m_dHighpassFreq(4)
, m_iFilterOrder(80)
, m_iDesignMethod(m_designMethods.indexOf(FilterParameter("Cosine")))
, m_iFilterType(m_filterTypes.indexOf(FilterParameter("BPF")))
, m_sFilterName("Unknown")
, m_sFilterShortDescription("")
{
    designFilter();
}

//=============================================================================================================

FilterKernel::FilterKernel(const QString& sFilterName,
                           int iFilterType,
                           int iOrder,
                           double dCenterfreq,
                           double dBandwidth,
                           double dParkswidth,
                           double dSFreq,
                           int iDesignMethod)
: m_sFreq(dSFreq)
, m_dCenterFreq(dCenterfreq)
, m_dBandwidth(dBandwidth)
, m_dParksWidth(dParkswidth)
, m_iFilterOrder(iOrder)
, m_iDesignMethod(iDesignMethod)
, m_iFilterType(iFilterType)
, m_sFilterName(sFilterName)
, m_sFilterShortDescription()
{
    if(m_iDesignMethod != 2){ //for FIR filters
        if(iOrder < 9) {
           qWarning() << "[FilterKernel::FilterKernel] Less than 9 taps were provided. Setting number of taps to 9.";
        }
    }

    designFilter();
}

//=============================================================================================================

void FilterKernel::prepareFilter(int iDataSize)
{
    int iFftLength, exp;

    iFftLength = iDataSize + m_vecCoeff.cols();
    exp = ceil(MNEMath::log2(iFftLength));
    iFftLength = pow(2, exp);

    // Transform coefficients anew if needed
    if(m_vecCoeff.cols() != (iFftLength/2+1)) {
        fftTransformCoeffs(iFftLength);
    }
}

//=============================================================================================================

RowVectorXd FilterKernel::applyConvFilter(const RowVectorXd& vecData,
                                          bool bKeepOverhead) const
{
    //Do zero padding or mirroring depending on user input
    RowVectorXd vecDataZeroPad = RowVectorXd::Zero(2*m_vecCoeff.cols() + vecData.cols());
    RowVectorXd vecFilteredTime = RowVectorXd::Zero(2*m_vecCoeff.cols() + vecData.cols());

    vecDataZeroPad.segment(m_vecCoeff.cols(), vecData.cols()) = vecData;

    //Do the convolution
    for(int i = m_vecCoeff.cols(); i < vecFilteredTime.cols(); i++) {
        vecFilteredTime(i-m_vecCoeff.cols()) = vecDataZeroPad.segment(i-m_vecCoeff.cols(),m_vecCoeff.cols()) * m_vecCoeff.transpose();
    }

    //Return filtered data
    if(!bKeepOverhead) {
        return vecFilteredTime.segment(m_vecCoeff.cols()/2, vecData.cols());
    }

    return vecFilteredTime.head(vecData.cols()+m_vecCoeff.cols());
}

//=============================================================================================================

void FilterKernel::applyFftFilter(RowVectorXd& vecData,
                                  bool bKeepOverhead)
{
    #ifdef EIGEN_FFTW_DEFAULT
    fftw_make_planner_thread_safe();
    #endif

    // Make sure we always have the correct FFT length for the given input data and filter overlap
    int iFftLength = vecData.cols() + m_vecCoeff.cols();
    int exp = ceil(MNEMath::log2(iFftLength));
    iFftLength = pow(2, exp);

    // Transform coefficients anew if needed
    if(m_vecFftCoeff.cols() != (iFftLength/2+1)) {
        fftTransformCoeffs(iFftLength);
    }

    //generate fft object
    Eigen::FFT<double> fft;
    fft.SetFlag(fft.HalfSpectrum);

    // Zero padd if necessary. Please note: The zero padding in Eigen's FFT is only working for column vectors -> We have to zero pad manually here
    int iOriginalSize = vecData.cols();
    if (vecData.cols() < iFftLength) {
        int iResidual = iFftLength - vecData.cols();
        vecData.conservativeResize(iFftLength);
        vecData.tail(iResidual).setZero();
    }

    //fft-transform data sequence
    RowVectorXcd vecFreqData;
    fft.fwd(vecFreqData, vecData, iFftLength);

    //perform frequency-domain filtering
    vecFreqData = m_vecFftCoeff.array() * vecFreqData.array();

    //inverse-FFT
    fft.inv(vecData, vecFreqData);

    //Return filtered data
    if(!bKeepOverhead) {
        vecData = vecData.segment(m_vecCoeff.cols()/2, iOriginalSize).eval();
    } else {
        vecData = vecData.head(iOriginalSize + m_vecCoeff.cols()).eval();
    }
}

//=============================================================================================================

void FilterKernel::applyIirFilter(Eigen::RowVectorXd& vecData)
{

    int iNumBiquads = m_vecRecCoeffsA.cols() / 3; //we have three a's (and b's) per biquad
    int iNumSamples = vecData.size(); //number of samples

    //this filer application implementation is adapted from I. Orifes filter application in his GitHub repository "Butterworth-Filter-Design" (see BiquadChain Class)

    //initialize delay variables and pointers
    //this part is from Orifes BiquadChain::allocate and BiquadChain::reset
    RowVectorXd _yn = RowVectorXd::Zero(iNumBiquads);
    RowVectorXd _yn1 = RowVectorXd::Zero(iNumBiquads);
    RowVectorXd _yn2 = RowVectorXd::Zero(iNumBiquads);
    double _xn1 = 0;
    double _xn2 = 0;
    double* yn = &_yn[0];
    double* yn1 = &_yn1[0];
    double* yn2 = &_yn2[0];


    //this part is adapted from Orifes BiquadChain::processBiquad
    for (int iSample{0}; iSample < iNumSamples; iSample++){

        //current input sample
        double xn = vecData[iSample];

        //iteration #0 as starting point
        yn[0] = m_vecRecCoeffsB[0] * xn + m_vecRecCoeffsB[1] * _xn1 + m_vecRecCoeffsB[2] * _xn2
                - m_vecRecCoeffsA[1] * yn1[0] - m_vecRecCoeffsA[2] * yn2[0];

        //iterations #1 to iNumBiquad-1
        for(int i = 1; i < iNumBiquads; i++){
            int iCoeffStartIdx = i * 3;
            yn[i] = m_vecRecCoeffsB[iCoeffStartIdx] * yn[i - 1]
                    + m_vecRecCoeffsB[iCoeffStartIdx+1] * yn1[i - 1]
                    + m_vecRecCoeffsB[iCoeffStartIdx+2] * yn2[i - 1]
                    - m_vecRecCoeffsA[iCoeffStartIdx+1] * yn1[i]
                    - m_vecRecCoeffsA[iCoeffStartIdx+2] * yn2[i];
        }

        // Shift delay elements so that they can be used for the next sample
        for(int i = 0; i < iNumBiquads; i++){
            yn2[i] = yn1[i];
            yn1[i] = yn[i];
        }
        _xn2 = _xn1;
        _xn1 = xn;

        //write result
        vecData[iSample] = yn[iNumBiquads - 1];
    }

}

//=============================================================================================================

QString FilterKernel::getName() const
{
    return m_sFilterName;
}

//=============================================================================================================

void FilterKernel::setName(const QString& sFilterName)
{
    m_sFilterName = sFilterName;
}

//=============================================================================================================

double FilterKernel::getSamplingFrequency() const
{
    return m_sFreq;
}

//=============================================================================================================

void FilterKernel::setSamplingFrequency(double dSFreq)
{
    m_sFreq = dSFreq;
}

//=============================================================================================================

int FilterKernel::getFilterOrder() const
{
    return m_iFilterOrder;
}

//=============================================================================================================

void FilterKernel::setFilterOrder(int iOrder)
{
    m_iFilterOrder = iOrder;
}

//=============================================================================================================

double FilterKernel::getCenterFrequency() const
{
    return m_dCenterFreq;
}

//=============================================================================================================

void FilterKernel::setCenterFrequency(double dCenterFreq)
{
    m_dCenterFreq = dCenterFreq;
}

//=============================================================================================================

double FilterKernel::getBandwidth() const
{
    return m_dBandwidth;
}

//=============================================================================================================

void FilterKernel::setBandwidth(double dBandwidth)
{
    m_dBandwidth = dBandwidth;
}

//=============================================================================================================

double FilterKernel::getParksWidth() const
{
    return m_dParksWidth;
}

//=============================================================================================================

void FilterKernel::setParksWidth(double dParksWidth)
{
    m_dParksWidth = dParksWidth;
}

//=============================================================================================================

double FilterKernel::getHighpassFreq() const
{
    return m_dHighpassFreq;
}

//=============================================================================================================

void FilterKernel::setHighpassFreq(double dHighpassFreq)
{
    m_dHighpassFreq = dHighpassFreq;
}

//=============================================================================================================

double FilterKernel::getLowpassFreq() const
{
    return m_dLowpassFreq;
}

//=============================================================================================================

void FilterKernel::setLowpassFreq(double dLowpassFreq)
{
    m_dLowpassFreq = dLowpassFreq;
}

//=============================================================================================================

Eigen::RowVectorXd FilterKernel::getCoefficients() const
{
    return m_vecCoeff;
}

//=============================================================================================================

void FilterKernel::setCoefficients(const Eigen::RowVectorXd& vecCoeff)
{
    m_vecCoeff = vecCoeff;
}

//=============================================================================================================

Eigen::RowVectorXcd FilterKernel::getFftCoefficients() const
{
    return m_vecFftCoeff;
}

//=============================================================================================================

void FilterKernel::setFftCoefficients(const Eigen::RowVectorXcd& vecFftCoeff)
{
    m_vecFftCoeff = vecFftCoeff;
}

//=============================================================================================================

Eigen::RowVectorXd FilterKernel::getRecCoeffsA() const
{
    return m_vecRecCoeffsA;
}

//=============================================================================================================

Eigen::RowVectorXd FilterKernel::getRecCoeffsB() const
{
    return m_vecRecCoeffsB;
}

//=============================================================================================================


bool FilterKernel::fftTransformCoeffs(int iFftLength)
{
    #ifdef EIGEN_FFTW_DEFAULT
        fftw_make_planner_thread_safe();
    #endif

    if(m_vecCoeff.cols() > iFftLength) {
        std::cout <<"[FilterKernel::fftTransformCoeffs] The number of filter taps is bigger than the FFT length."<< std::endl;
        return false;
    }

    //generate fft object
    Eigen::FFT<double> fft;
    fft.SetFlag(fft.HalfSpectrum);

    // Zero padd if necessary. Please note: The zero padding in Eigen's FFT is only working for column vectors -> We have to zero pad manually here
    RowVectorXd vecInputFft;
    if (m_vecCoeff.cols() < iFftLength) {
        vecInputFft.setZero(iFftLength);
        vecInputFft.block(0,0,1,m_vecCoeff.cols()) = m_vecCoeff;
    } else {
        vecInputFft = m_vecCoeff;
    }

    //fft-transform filter coeffs
    RowVectorXcd vecFreqData;
    fft.fwd(vecFreqData, vecInputFft, iFftLength);
    m_vecFftCoeff = vecFreqData;;

    return true;
}

//=============================================================================================================

void FilterKernel::designFilter()
{

    // Make sure we only use a minimum needed FFT size
    int iFftLength = m_iFilterOrder;
    int exp = ceil(MNEMath::log2(iFftLength));
    iFftLength = pow(2, exp);

    switch(m_iDesignMethod) {
        case 1: {
            ParksMcClellan filter(m_iFilterOrder,
                                  m_dCenterFreq,
                                  m_dBandwidth,
                                  m_dParksWidth,
                                  static_cast<ParksMcClellan::TPassType>(m_iFilterType));
            m_vecCoeff = filter.FirCoeff;

            //fft-transform m_vecCoeff in order to be able to perform frequency-domain filtering
            fftTransformCoeffs(iFftLength);

            break;
        }

        case 0: {
            CosineFilter filtercos;

            switch(m_iFilterType) {
                case 0:
                    filtercos = CosineFilter(iFftLength,
                                             (m_dCenterFreq)*(m_sFreq/2.),
                                             m_dParksWidth*(m_sFreq/2),
                                             (m_dCenterFreq)*(m_sFreq/2),
                                             m_dParksWidth*(m_sFreq/2),
                                             m_sFreq,
                                             static_cast<CosineFilter::TPassType>(m_iFilterType));

                    break;

                case 1:
                    filtercos = CosineFilter(iFftLength,
                                             (m_dCenterFreq)*(m_sFreq/2),
                                             m_dParksWidth*(m_sFreq/2),
                                             (m_dCenterFreq)*(m_sFreq/2),
                                             m_dParksWidth*(m_sFreq/2),
                                             m_sFreq,
                                             static_cast<CosineFilter::TPassType>(m_iFilterType));

                    break;

                case 2:
                    filtercos = CosineFilter(iFftLength,
                                             (m_dCenterFreq + m_dBandwidth/2)*(m_sFreq/2),
                                             m_dParksWidth*(m_sFreq/2),
                                             (m_dCenterFreq - m_dBandwidth/2)*(m_sFreq/2),
                                             m_dParksWidth*(m_sFreq/2),
                                             m_sFreq,
                                             static_cast<CosineFilter::TPassType>(m_iFilterType));

                    break;
            }

            //This filter is designed in the frequency domain, hence the time domain impulse response need to be shortend by the users dependent number of taps
            m_vecCoeff.resize(m_iFilterOrder);

            m_vecCoeff.head(m_iFilterOrder/2) = filtercos.m_vecCoeff.tail(m_iFilterOrder/2);
            m_vecCoeff.tail(m_iFilterOrder/2) = filtercos.m_vecCoeff.head(m_iFilterOrder/2);

            //Now generate the fft version of the shortened impulse response
            fftTransformCoeffs(iFftLength);

            break;            
        }

        case 2:{


            //create new IIR Butterworth filter
            Butterworth butterfilt(m_iFilterType,
                                   m_iFilterOrder,
                                   m_dCenterFreq*(m_sFreq/2),
                                   m_dBandwidth*(m_sFreq/2),
                                   m_sFreq);

            //get calculated recursion coefficients of the new Butterworth filter
            //coefficients are in second order section form (a0,a1,a2 of first biquad then a0,a1,a2 of second ...)
            m_vecRecCoeffsA = butterfilt.m_vecRecCoeffsA;
            m_vecRecCoeffsB = butterfilt.m_vecRecCoeffsB;

        }
    }

    if(m_iDesignMethod != 2){ //only true for FIR filter design methods (cosine and parksmcclellan)
        switch(m_iFilterType) {
            case 0:
                m_dLowpassFreq = 0;
                m_dHighpassFreq = m_dCenterFreq*(m_sFreq/2);
            break;

            case 1:
                m_dLowpassFreq = m_dCenterFreq*(m_sFreq/2);
                m_dHighpassFreq = 0;
            break;

            case 2:
                m_dLowpassFreq = (m_dCenterFreq + m_dBandwidth/2)*(m_sFreq/2);
                m_dHighpassFreq = (m_dCenterFreq - m_dBandwidth/2)*(m_sFreq/2);
            break;
        }
    }

    getShortDescription();
}


//=============================================================================================================

QString FilterKernel::getShortDescription() const
{
    QString description(m_designMethods.at(m_iDesignMethod).getName() + "  -  " + \
                                QString::number(m_dHighpassFreq,'g',4) + "Hz to " + QString::number(m_dLowpassFreq,'g',4) + "Hz  -  " \
                                "Ord: " + QString::number(m_iFilterOrder));
    return description;
}

//=============================================================================================================

RTPROCESSINGLIB::FilterParameter FilterKernel::getDesignMethod() const
{
    if(m_iDesignMethod < 0){
        return m_designMethods.at(0);
    }
    return m_designMethods.at(m_iDesignMethod);
}

//=============================================================================================================

RTPROCESSINGLIB::FilterParameter FilterKernel::getFilterType() const
{
    if(m_iFilterType < 0){
        return m_filterTypes.at(0);
    }
    return m_filterTypes.at(m_iFilterType);
}

//=============================================================================================================

void FilterKernel::setDesignMethod(int iDesignMethod)
{
    if(iDesignMethod < 0){
        m_iDesignMethod = 0;
    } else {
        m_iDesignMethod = iDesignMethod;
    }
}

//=============================================================================================================

void FilterKernel::setFilterType(int iFilterType)
{
    if(iFilterType < 0){
        m_iFilterType = 0;
    } else {
        m_iFilterType = iFilterType;
    }
}

//=============================================================================================================

FilterParameter::FilterParameter()
:FilterParameter("Unknown", "")
{
}

//=============================================================================================================

FilterParameter::FilterParameter(QString sName)
:FilterParameter(sName,"")
{  
}

//=============================================================================================================

FilterParameter::FilterParameter(QString sName,
                           QString sDescription)
: m_sName(sName)
, m_sDescription(sDescription)
{

}

//=============================================================================================================

QString FilterParameter::getName() const
{
    return m_sName;
}

//=============================================================================================================
