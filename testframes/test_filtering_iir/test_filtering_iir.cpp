//=============================================================================================================
/**
 * @file     test_filtering_iir.cpp
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de>
 * @since    0.1.0
 * @date     June, 2022
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
 * @brief     Test for Butterworth filter creation and filtering.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include <utils/generics/applicationlogger.h>

#include <iostream>
#include <vector>
#include <complex>
#include <math.h>


#include <fiff/fiff.h>
#include <rtprocessing/helpers/butterworth.h>
#include <rtprocessing/helpers/filterkernel.h>
#include <rtprocessing/filter.h>
#include <utils/ioutils.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QtCore/QCoreApplication>
#include <QFile>
#include <QCommandLineParser>
#include <QtTest>

//=============================================================================================================
// Eigen
//=============================================================================================================

#include <Eigen/Dense>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace FIFFLIB;
using namespace UTILSLIB;
using namespace RTPROCESSINGLIB;
using namespace Eigen;

//=============================================================================================================



/**
 * DECLARE CLASS TestFilteringIir
 *
 * @brief The TestFilteringIir class provides Butterworth filter creation and filter application tests.
 *
 */
class TestFilteringIir: public QObject
{
    Q_OBJECT

public:
    TestFilteringIir();

private slots:
    void initTestCase();

    //test Butterworth filter kernel cration
    void comparePrototypePoles();
    void comparePrototypeGain();

    void compareAnalogPoles();
    void compareAnalogZeros();
    void compareAnalogGain();

    void compareDigitalPoles();
    void compareDigitalZeros();
    void compareDigitalGain();

    //test Butterworth offline filtering
    void compareData();
    void compareTimes();

    void cleanupTestCase();

private:
    // declare your thresholds, variables and error values here
    double dEpsilon;   
    int iOrder;

    RowVectorXcd m_vecPrototypePoles;
    RowVectorXcd m_vecRefPrototypePoles;
    double m_dPrototypeGain;
    double m_dRefPrototypeGain;

    RowVectorXcd m_vecAnalogPoles;
    RowVectorXcd m_vecRefAnalogPoles;
    RowVectorXcd m_vecAnalogZeros;
    RowVectorXcd m_vecRefAnalogZeros;
    double m_dAnalogGain;
    double m_dRefAnalogGain;

    RowVectorXcd m_vecDigitalPoles;
    RowVectorXcd m_vecRefDigitalPoles;
    RowVectorXcd m_vecDigitalZeros;
    RowVectorXcd m_vecRefDigitalZeros;
    double m_dDigitalGain;
    double m_dRefDigitalGain;

    MatrixXd mFirstInData;
    MatrixXd mFirstInTimes;
    MatrixXd mFirstFiltered;

    MatrixXd mRefInData;
    MatrixXd mRefInTimes;
    MatrixXd mRefFiltered;

};

//=============================================================================================================

TestFilteringIir::TestFilteringIir()
    : dEpsilon(DBL_EPSILON)
{
}

//=============================================================================================================

void TestFilteringIir::initTestCase()
{
    qInstallMessageHandler(UTILSLIB::ApplicationLogger::customLogWriter);
    qDebug() << "Epsilon" << dEpsilon;

    QFile t_fileIn(QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/MEG/sample/sample_audvis_trunc_raw_no_time_offset.fif");
    QFile t_fileOut(QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/MEG/sample/sample_audvis_trunc_raw_butter_filt_out.fif");

/*    // Filter in Python is created with following function: mne.filter.design_mne_c_filter(raw.info['sfreq'], 5, 10, 1, 1)
    // This will create a filter with with 8193 elements/taps/Order. In order to be concise with the MNE-CPP implementation
    // the filter is cut to the Order used in mne-cpp (1024, see below).//
    // The actual filtering was performed with the function: mne.filter._overlap_add_filter(dataIn, filter_python, phase = 'linear')
*/


    //TODO: edit docu here for used matlab filter for reference file
//    QFile t_fileRef(QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/Result/ref_sample_audvis_trunc_raw_butter_BP_filt.fif");
//    QFile t_fileRef(QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/Result/ref_sample_audvis_trunc_raw_butter_BS_filt.fif");
    QFile t_fileRef(QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/Result/ref_sample_audvis_trunc_raw_butter_LP_filt.fif");
//    QFile t_fileRef(QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/Result/ref_sample_audvis_trunc_raw_butter_HP_filt.fif");


    // Make sure test folder exists
    QFileInfo t_fileOutInfo(t_fileOut);
    QDir().mkdir(t_fileOutInfo.path());

    // Setup for reading the raw data
    FiffRawData rawFirstInRaw;
    rawFirstInRaw = FiffRawData(t_fileIn);

    //*********************************************************************************************************
    // Creation of Butterworth Filter Kernel
    //*********************************************************************************************************

    printf(">>>>>>>>>>>>>>>>>>>>>>>>> Creation of Butterworth Filter Kernel >>>>>>>>>>>>>>>>>>>>>>>>>\n");

    // fifth order butterworth filter with center frequency 10, bandwith 10 and specified type (LPF, HPF, BPF or NOTCH)

    //for lowpass
    int iFilterType = FilterKernel::m_filterTypes.indexOf(FilterParameter("LPF"));
    iOrder = 16; //filter order 16 for unit test
    double dCenterFreq = 40; //center frequency in Hz (cutoff frequency for LP)
    double dBandwidth = 0; //bandwidth in Hz (not needed for LP)



/*    //for highpass
    int iFilterType = FilterKernel::m_filterTypes.indexOf(FilterParameter("HPF"));
    iOrder = 6; //filter order 6
    double dCenterFreq = 4; //center frequency in Hz (cutoff frequency for HP)
    double dBandwidth = 0; //bandwidth in Hz (not needed for HP)
*/


/*    //init filter parameters for bandpass
    int iFilterType = FilterKernel::m_filterTypes.indexOf(FilterParameter("BPF"));
    iOrder = 5; //filter order 5
    double dCenterFreq = 10; //center frequency in Hz
    double dBandwidth = 10; //bandwidth in Hz
*/

/*    //for bandstop
    int iFilterType = FilterKernel::m_filterTypes.indexOf(FilterParameter("NOTCH"));
    iOrder = 10; //filter order 10 for unit test
    double dCenterFreq = 60; //center frequency in Hz
    double dBandwidth = 4; //bandwidth in Hz
*/


    double dSFreq = rawFirstInRaw.info.sfreq; //get sampling frequency of data that is filtered later (dSFreq = 300.3074951171875)

    //normalize and denormalize parameters (according to the procedure during filterData, so that filter parameters are consistent with filter parameters generated during filter application test by using filterData)
    dCenterFreq = (dCenterFreq / (dSFreq/2.0)) * (dSFreq/2.0);
    dBandwidth = (dBandwidth / (dSFreq/2.0)) * (dSFreq/2.0);


    //create Butterworth filter kernel
    Butterworth butterfilt(iFilterType,
                           iOrder,
                           dCenterFreq,
                           dBandwidth,
                           dSFreq);

    //init test variables with filter kernel results
    m_vecPrototypePoles = butterfilt.getPrototypePoles();
    m_dPrototypeGain = butterfilt.getPrototypeGain();

    m_vecAnalogPoles = butterfilt.getAnalogPoles();
    m_vecAnalogZeros = butterfilt.getAnalogZeros();
    m_dAnalogGain = butterfilt.getAnalogGain();

    m_vecDigitalPoles = butterfilt.getDigitalPoles();
    m_vecDigitalZeros = butterfilt.getDigitalZeros();
    m_dDigitalGain = butterfilt.getDigitalGain();

    printf(">>>>>>>>>>>>>>>>>>>>>>>>> Creation of Butterworth Filter Kernel Finished >>>>>>>>>>>>>>>>>>>>>>>>>\n");


    //*********************************************************************************************************
    // Read Matlab Results for Butterworth Kernel Creation As Reference
    //*********************************************************************************************************


    printf(">>>>>>>>>>>>>>>>>>>>>>>>> Read Matlab Results For Butterworth Kernel Creation As Reference >>>>>>>>>>>>>>>>>>>>>>>>>\n");

    //TODO: edit documentation here
    //The txt file contains the prototype poles calcuated with Matlab's function buttap(5) (5th order prototype) in the first five rows
    //Each row represents one pole, where the fist column contais the real part, the second column the imaginary part


    //set reference file name depending on filter type
    QString refKernelFileName;
    switch(iFilterType){
        case 0:{ //lowpass filter
            refKernelFileName = QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/Result/butterworthFilter/ref_test_butter_LP_creation.txt";
            break;
        }
        case 1:{ //highpass filter
            refKernelFileName = QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/Result/butterworthFilter/ref_test_butter_HP_creation.txt";
            break;
        }
        case 2:{ //bandpass filter
            refKernelFileName = QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/Result/butterworthFilter/ref_test_butter_BP_creation.txt";
            break;
        }
        case 3:{ //bandstop filter
            refKernelFileName = QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/Result/butterworthFilter/ref_test_butter_BS_creation.txt";
            break;
        }
    }

    Eigen::MatrixXd matRefData;
    UTILSLIB::IOUtils::read_eigen_matrix(matRefData, refKernelFileName);

    //counter for tracking the rownumber to of matRefData to read from
    int iRefDataRowIdx{0};

    //PROTOTYPE

    //connect real and imaginary parts of reference prototype poles to one complex number for each pole
    m_vecRefPrototypePoles = RowVectorXcd::Zero(iOrder);
    for(int iPole{iRefDataRowIdx}; iPole < m_vecRefPrototypePoles.size(); iPole++){
        m_vecRefPrototypePoles[iPole] = std::complex<double> (matRefData(iPole,1), matRefData(iPole,2));
    }
    iRefDataRowIdx += m_vecRefPrototypePoles.size(); //update row counter for reference data matrix

    //read prototype gain as reference
    m_dRefPrototypeGain = matRefData(iRefDataRowIdx,1);
    iRefDataRowIdx++;

    //ANALOG FILTER

    //read analog poles and zeros as reference (number of poles and zeros depends on filter type)
    int iNumAnalogPoles;
    int iNumAnalogZeros;
    switch(iFilterType){
        case 0:{ //lowpass filter
            iNumAnalogPoles = m_vecRefPrototypePoles.size();
            iNumAnalogZeros = 0;
            break;
        }
        case 1:{ //highpass filter
            iNumAnalogPoles = m_vecRefPrototypePoles.size();
            iNumAnalogZeros = iNumAnalogPoles;
            break;
        }
        case 2:{ //bandpass filter
            iNumAnalogPoles = 2 * m_vecRefPrototypePoles.size();
            iNumAnalogZeros = m_vecRefPrototypePoles.size();
            break;
        }
        case 3:{ //bandstop filter
            iNumAnalogPoles = 2 * m_vecRefPrototypePoles.size();
            iNumAnalogZeros = iNumAnalogPoles;
            break;
        }
    }

    //read analog poles as reference
    m_vecRefAnalogPoles = RowVectorXcd::Zero(iNumAnalogPoles);
    for(int iPole{iRefDataRowIdx}, iCount{0}; iPole < (iRefDataRowIdx + iNumAnalogPoles); iPole++, iCount++){
        m_vecRefAnalogPoles[iCount] = std::complex<double> (matRefData(iPole,1), matRefData(iPole,2));
    }
    iRefDataRowIdx += m_vecRefAnalogPoles.size();

    //read analog zeros as reference
    m_vecRefAnalogZeros = RowVectorXcd::Zero(iNumAnalogZeros);
    for(int iZero{iRefDataRowIdx}, iCount{0}; iZero < (iRefDataRowIdx + iNumAnalogZeros); iZero++, iCount++){
        m_vecRefAnalogZeros[iCount] = std::complex<double> (matRefData(iZero,1), matRefData(iZero,2));
    }
    iRefDataRowIdx += m_vecRefAnalogZeros.size();

    //read analog gain as reference
    m_dRefAnalogGain = matRefData(iRefDataRowIdx,1);
    iRefDataRowIdx++;

    //DIGITAL FILTER

    int iNumDigitalPoles = iNumAnalogPoles;
    int iNumDigitalZeros = iNumAnalogPoles;

    //read digital poles as reference
    m_vecRefDigitalPoles = RowVectorXcd::Zero(iNumDigitalPoles);
    for(int iPole{iRefDataRowIdx}, iCount{0}; iPole < (iRefDataRowIdx + iNumDigitalPoles); iPole++, iCount++){
        m_vecRefDigitalPoles[iCount] = std::complex<double> (matRefData(iPole,1), matRefData(iPole,2));
    }
    iRefDataRowIdx += m_vecRefDigitalPoles.size();

    //read digital zeros as reference
    m_vecRefDigitalZeros = RowVectorXcd::Zero(iNumDigitalZeros);
    for(int iZero{iRefDataRowIdx}, iCount{0}; iZero < (iRefDataRowIdx + iNumDigitalZeros); iZero++, iCount++){
        m_vecRefDigitalZeros[iCount] = std::complex<double> (matRefData(iZero,1), matRefData(iZero,2));
    }
    iRefDataRowIdx += m_vecRefDigitalZeros.size();

    //read digital gain as reference
    m_dRefDigitalGain = matRefData(iRefDataRowIdx,1);
    iRefDataRowIdx++;


    printf(">>>>>>>>>>>>>>>>>>>>>>>>> Read Matlab Results For Butterworth Kernel Creation As Reference Finished>>>>>>>>>>>>>>>>>>>>>>>>>\n");

    //*********************************************************************************************************
    // Application of Butterworth Filter Kernel
    //*********************************************************************************************************

    printf(">>>>>>>>>>>>>>>>>>>>>>>>> Application of Butterworth Filter Kernel (Read, Filter and Write) >>>>>>>>>>>>>>>>>>>>>>>>>\n");

    //TODO: code copied parts from test_filtering.cpp

    //TODO: copied
    //Only filter MEG channels (eeg channels were set to true in test_filtering, changed it)
    RowVectorXi vPicks = rawFirstInRaw.info.pick_types(true, false, false);
    RowVectorXd vCals;
    FiffStream::SPtr outfid = FiffStream::start_writing_raw(t_fileOut, rawFirstInRaw.info, vCals);

    //TODO:copied
    // Set up the reading parameters
    // To read the whole file at once set
    fiff_int_t from = rawFirstInRaw.first_samp;
    fiff_int_t to = rawFirstInRaw.last_samp;

    MatrixXd mDataFiltered;

    //TODO: copied
    // Reading
    if(!rawFirstInRaw.read_raw_segment(mFirstInData, mFirstInTimes, from, to)) {
        printf("error during read_raw_segment\n");
    }


    //TODO: copied
    // Filtering
    printf("Filtering...");

    //TODO: dTransition is ignored for IIR filter kernel; its value is arbitrary set to one in order to use existing filterData
    //only forward filtering
    mFirstFiltered = RTPROCESSINGLIB::filterData(mFirstInData,
                                                 iFilterType,
                                                 dCenterFreq,
                                                 dBandwidth,
                                                 1,
                                                 dSFreq,
                                                 iOrder,
                                                 RTPROCESSINGLIB::FilterKernel::m_designMethods.indexOf(FilterParameter("Butterworth")),
                                                 vPicks,
                                                 true,
                                                 false,
                                                 false);


    //test twopass filtering
/*    mFirstFiltered = RTPROCESSINGLIB::filterData(mFirstInData,
                                                 iFilterType,
                                                 dCenterFreq,
                                                 dBandwidth,
                                                 1,
                                                 dSFreq,
                                                 iOrder,
                                                 RTPROCESSINGLIB::FilterKernel::m_designMethods.indexOf(FilterParameter("Butterworth")),
                                                 vPicks,
                                                 true,
                                                 false,
                                                 true);

*/


    //TODO:copied
    printf("[done]\n");

        //TODO:copied
    // Writing
    printf("Writing...");
    outfid->write_int(FIFF_FIRST_SAMPLE, &from);
    outfid->write_raw_buffer(mFirstFiltered,vCals);
    printf("[done]\n");

        //TODO:copied
    outfid->finish_writing_raw();

        //TODO:copied
    // Read filtered data from the filtered output file to check if read and write is working correctly
    FiffRawData rawSecondInRaw;
    rawSecondInRaw = FiffRawData(t_fileOut);

        //TODO:copied
    // Reading
    if (!rawSecondInRaw.read_raw_segment(mFirstFiltered,mFirstInTimes,from,to,vPicks)) {
        printf("error during read_raw_segment\n");
    }

    printf("<<<<<<<<<<<<<<<<<<<<<<<<< Application of Butterworth Filter Kernel (Read, Filter and Write) Finished <<<<<<<<<<<<<<<<<<<<<<<<<\n");

    //*********************************************************************************************************
    // Read MATLAB Results For Filter Application As Reference
    //*********************************************************************************************************

    printf(">>>>>>>>>>>>>>>>>>>>>>>>> Read Matlab Results For Filter Application As Reference >>>>>>>>>>>>>>>>>>>>>>>>>\n");

    FiffRawData ref_in_raw;
    ref_in_raw = FiffRawData(t_fileRef);

    //TODO: from and to should be equal to ref_from and ref_to but they are not due modification during fieldtrip2fif in Matlab
    //TODO: delete this part when fieldtrip modification is fixed
    fiff_int_t ref_from = ref_in_raw.first_samp;
    fiff_int_t ref_to = ref_in_raw.last_samp;

    // Reading
    if (!ref_in_raw.read_raw_segment(mRefFiltered,mRefInTimes,ref_from,ref_to,vPicks)) {
        printf("error during read_raw_segment\n");
    }

    printf(">>>>>>>>>>>>>>>>>>>>>>>>> Read Matlab Results For Filter Application As Reference Finished>>>>>>>>>>>>>>>>>>>>>>>>>\n");

}

//=============================================================================================================

void TestFilteringIir::comparePrototypePoles()
{
    // compare prototype poles to reference prototype poles
    Eigen::RowVectorXcd vecPrototypePoleDiffComplex = m_vecPrototypePoles.array() - m_vecRefPrototypePoles.array();
    Eigen::RowVectorXd vecPrototypePoleDiff = vecPrototypePoleDiffComplex.array().abs();
    qDebug() << "[TestFilteringIir::comparePrototypePoles] Maximum of absolute differences = " << vecPrototypePoleDiff.maxCoeff();
    QVERIFY( vecPrototypePoleDiff.maxCoeff() < dEpsilon );

}

//=============================================================================================================

void TestFilteringIir::comparePrototypeGain()
{
     //compare prototype gain to reference prototype gain
    double dPrototypeGainDiff = abs(m_dPrototypeGain - m_dRefPrototypeGain);
    qDebug() << "[TestFilteringIir::comparePrototypeGain] Absolute difference = " << dPrototypeGainDiff;
    QVERIFY( dPrototypeGainDiff < dEpsilon );

}



//=============================================================================================================

void TestFilteringIir::compareAnalogPoles()
{
    // compare prototype poles to reference prototype poles
    Eigen::RowVectorXcd vecAnalogPoleDiffComplex = m_vecAnalogPoles.array() - m_vecRefAnalogPoles.array();
    Eigen::RowVectorXd vecAnalogPoleDiff = vecAnalogPoleDiffComplex.array().abs();
    qDebug() << "[TestFilteringIir::compareAnalogPoles] Maximum of absolute differences = " << vecAnalogPoleDiff.maxCoeff();
    QVERIFY( vecAnalogPoleDiff.maxCoeff() < dEpsilon );

}

//=============================================================================================================

void TestFilteringIir::compareAnalogZeros()
{

    if(m_vecAnalogZeros.size() != 0){ //this unit test is not for LP filters (they have no analog zeros)
        // compare prototype poles to reference prototype poles
        Eigen::RowVectorXcd vecAnalogZeroDiffComplex = m_vecAnalogZeros.array() - m_vecRefAnalogZeros.array();
        Eigen::RowVectorXd vecAnalogZeroDiff = vecAnalogZeroDiffComplex.array().abs();
        qDebug() << "[TestFilteringIir::compareAnalogZeros] Maximum of absolute differences = " << vecAnalogZeroDiff.maxCoeff();
        QVERIFY( vecAnalogZeroDiff.maxCoeff() < dEpsilon );
    }else{
        qDebug() << "[TestFilteringIir::compareAnalogZeros] Your LP filter has no analog zeros.";
    }

}

//=============================================================================================================

void TestFilteringIir::compareAnalogGain()
{
     //compare prototype gain to reference prototype gain
    double dAnalogGainDiff = abs(m_dAnalogGain - m_dRefAnalogGain);
    qDebug() << "[TestFilteringIir::compareAnalogGain] Absolute difference = " << dAnalogGainDiff;
    QVERIFY( dAnalogGainDiff < dEpsilon );

}

//=============================================================================================================

void TestFilteringIir::compareDigitalPoles()
{
    // compare prototype poles to reference prototype poles
    Eigen::RowVectorXcd vecDigitalPoleDiffComplex = m_vecDigitalPoles.array() - m_vecRefDigitalPoles.array();
    Eigen::RowVectorXd vecDigitalPoleDiff = vecDigitalPoleDiffComplex.array().abs();
    qDebug() << "[TestFilteringIir::compareDigitalPoles] Maximum of absolute differences = " << vecDigitalPoleDiff.maxCoeff();
    QVERIFY( vecDigitalPoleDiff.maxCoeff() < dEpsilon );

}

//=============================================================================================================

void TestFilteringIir::compareDigitalZeros()
{
    // compare prototype poles to reference prototype poles
    Eigen::RowVectorXcd vecDigitalZeroDiffComplex = m_vecDigitalZeros.array() - m_vecRefDigitalZeros.array();
    Eigen::RowVectorXd vecDigitalZeroDiff = vecDigitalZeroDiffComplex.array().abs();
    qDebug() << "[TestFilteringIir::compareDigitalZeros] Maximum of absolute differences = " << vecDigitalZeroDiff.maxCoeff();
    QVERIFY( vecDigitalZeroDiff.maxCoeff() < dEpsilon );

}

//=============================================================================================================

void TestFilteringIir::compareDigitalGain()
{
    //compare prototype gain to reference prototype gain
    double dDigitalGainDiff = abs(m_dDigitalGain - m_dRefDigitalGain);
    qDebug() << "[TestFilteringIir::compareDigitalGain] Absolute difference = " << dDigitalGainDiff;
    QVERIFY( dDigitalGainDiff < dEpsilon );

}

//=============================================================================================================

void TestFilteringIir::compareData()
{
    //compare filtered data
    MatrixXd mDataDiff = mFirstFiltered - mRefFiltered;
    mDataDiff = mDataDiff.cwiseAbs();
    qDebug() << "[TestFilteringIir::compareData] Maximum of absolute differences = " << mDataDiff.maxCoeff();
    QVERIFY( mDataDiff.maxCoeff() < dEpsilon );

}

//=============================================================================================================

void TestFilteringIir::compareTimes()
{
    MatrixXd mTimesDiff = mFirstInTimes - mRefInTimes;
    mTimesDiff = mTimesDiff.cwiseAbs();
    qDebug() << "[TestFilteringIir::compareTimes] Maximum of absolute differences = " << mTimesDiff.maxCoeff();
    QVERIFY( mTimesDiff.maxCoeff() < dEpsilon );
}

void TestFilteringIir::cleanupTestCase()
{
/*  //TODO: enable this in case all unit tests pass
    QFile t_fileOut(QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/MEG/sample/sample_audvis_trunc_raw_butter_filt_out.fif");
    t_fileOut.remove();
*/
}

//=============================================================================================================
// MAIN
//=============================================================================================================

QTEST_GUILESS_MAIN(TestFilteringIir)
#include "test_filtering_iir.moc"

