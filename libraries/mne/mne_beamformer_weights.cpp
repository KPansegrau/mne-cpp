//=============================================================================================================
/**
 * @file     mne_beamformer_weights.cpp
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de>
 * @since    0.1.9
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
 * @brief    MNEBeamformerWeights class definition.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "mne_beamformer_weights.h"

#include <math.h>

//TODO: delete later (only for debugging)
#include <fstream>
#include <utils/ioutils.h>


//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QFuture>
#include <QtConcurrent>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

#include <Eigen/Dense>
#include <Eigen/SVD>

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace MNELIB;
using namespace FIFFLIB;
using namespace Eigen;
using namespace UTILSLIB;

//=============================================================================================================
// DEFINE GLOBAL METHODS
//=============================================================================================================

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

//HINT this constructor is analog to the one for inverse operator class

MNEBeamformerWeights::MNEBeamformerWeights()
//: weights(defaultMatrixXd)
    : data_cov(new FiffCov)
, noise_cov(new FiffCov)
, fixedOri(false)
, nsource(-1)
, nchan(-1)
, nave(-1)
{

    qRegisterMetaType<QSharedPointer<MNELIB::MNEBeamformerWeights> >("QSharedPointer<MNELIB::MNEBeamformerWeights>");
    qRegisterMetaType<MNELIB::MNEBeamformerWeights>("MNELIB::MNEBeamformerWeights");
}



//=============================================================================================================


MNEBeamformerWeights::MNEBeamformerWeights(const FiffInfo &p_dataInfo,
                                           const MNEForwardSolution& p_forward,
                                           const FiffCov& p_dataCov,
                                           const FiffCov& p_noiseCov,
                                           QString p_sWeightnorm,
                                           bool p_bFixedOri,                                         
                                           qint32 p_iRegParam,
                                           qint32 p_iNAverages)
{

    *this = MNEBeamformerWeights::make_beamformer_weights(p_dataInfo,p_forward,p_dataCov,p_noiseCov,p_bFixedOri,p_sWeightnorm,p_iRegParam,p_iNAverages);

   qRegisterMetaType<QSharedPointer<MNELIB::MNEBeamformerWeights> >("QSharedPointer<MNELIB::MNEBeamformerWeights>");
   qRegisterMetaType<MNELIB::MNEBeamformerWeights>("MNELIB::MNEBeamformerWeights");

}


//=============================================================================================================

//HINT this constructor is analog to constructor in inverse operator class

MNEBeamformerWeights::MNEBeamformerWeights(const MNEBeamformerWeights &p_MNEBeamformerWeights)
    : info(p_MNEBeamformerWeights.info)
    , weights(p_MNEBeamformerWeights.weights)
    , data_cov(p_MNEBeamformerWeights.data_cov)
    , noise_cov(p_MNEBeamformerWeights.noise_cov)
    , weightNorm(p_MNEBeamformerWeights.weightNorm)
    , whitener(p_MNEBeamformerWeights.whitener)
    , fixedOri(p_MNEBeamformerWeights.fixedOri)
    , optOri(p_MNEBeamformerWeights.optOri)
    , nsource(p_MNEBeamformerWeights.nsource)
    , nchan(p_MNEBeamformerWeights.nchan)
    , projs(p_MNEBeamformerWeights.projs)
    , proj(p_MNEBeamformerWeights.proj)
    , src(p_MNEBeamformerWeights.src)
    , nave(p_MNEBeamformerWeights.nave)

{

    qRegisterMetaType<QSharedPointer<MNELIB::MNEBeamformerWeights> >("QSharedPointer<MNELIB::MNEBeamformerWeights>");
    qRegisterMetaType<MNELIB::MNEBeamformerWeights>("MNELIB::MNEBeamformerWeights");
}



//=============================================================================================================

MNEBeamformerWeights::~MNEBeamformerWeights()

{
}

//=============================================================================================================

MatrixXd MNEBeamformerWeights::reduceLocLeadfieldRank(const MatrixXd &matLocLeadfield){

    //decompose matrix with svd
    JacobiSVD<MatrixXd> svdLF(matLocLeadfield, ComputeThinU | ComputeThinV);
    VectorXd vecSingLocLF = svdLF.singularValues(); //sorted in decreasing order
    MatrixXd matULocLF = svdLF.matrixU();
    MatrixXd matVLocLF = svdLF.matrixV();

//    std::cout << "[MNEBeamformerWeights::reduceLocLeadfieldRank] vecSingLocLF: " << vecSingLocLF;

    //set smallest singular value to zero
    vecSingLocLF[vecSingLocLF.size() -1] = 0;

//    std::cout << "[MNEBeamformerWeights::reduceLocLeadfieldRank] vecSingLocLF: " << vecSingLocLF;

//    qDebug() << "[MNEBeamformerWeights::reduceLocLeadfieldRank] U dim : " << matULocLF.rows() << " x " << matULocLF.cols();
//    qDebug() << "[MNEBeamformerWeights::reduceLocLeadfieldRank] V dim : " << matVLocLF.rows() << " x " << matVLocLF.cols();
//    qDebug() << "[MNEBeamformerWeights::reduceLocLeadfieldRank] sing.asDiagonal() dim : " << vecSingLocLF.asDiagonal().rows() << " x " << vecSingLocLF.asDiagonal().cols();


    //backproject to rank reduced local leadfield matrix
    MatrixXd matRankReducedLF = matULocLF * vecSingLocLF.asDiagonal() * matVLocLF.adjoint();

    qDebug() << "[MNEBeamformerWeights::reduceLocLeadfieldRank] matRankReducedLF dim : " << matRankReducedLF.rows() << " x " << matRankReducedLF.cols();

    return matRankReducedLF;

}


//=============================================================================================================

bool MNEBeamformerWeights::check_ch_names(const FiffInfo &info) const
{
    //HINT: copied from MNEInverseOperator::check_ch_names and modified for beamformer purposes

    QStringList bf_ch_names = this->info.ch_names;


    bool t_bContainsNoiseCov = true;
    if(this->info.ch_names.size() != this->noise_cov->names.size())
        t_bContainsNoiseCov = false;
    else
    {
        for(qint32 i = 0; i < this->noise_cov->names.size(); ++i)
        {
            if(bf_ch_names[i] != this->noise_cov->names[i])
            {
                t_bContainsNoiseCov = false;
                break;
            }
        }
    }

    if(!t_bContainsNoiseCov)
    {
        qCritical("[MNEBeamformerWeights::check_ch_names] Channels of beamformer do not match noise covariance channels.");
        return false;
    }


    bool t_bContainsDataCov = true;
    if(this->info.ch_names.size() != this->data_cov->names.size())
        t_bContainsDataCov = false;
    else
    {
        for(qint32 i = 0; i < this->data_cov->names.size(); ++i)
        {
            if(bf_ch_names[i] != this->data_cov->names[i])
            {
                t_bContainsDataCov = false;
                break;
            }
        }
    }

    if(!t_bContainsDataCov)
    {
        qCritical("[MNEBeamformerWeights::check_ch_names] Channels of beamformer do not match data covariance channels.");
        return false;
    }




    QStringList data_ch_names = info.ch_names;

    //TODO: for debugging delete later
 //   qDebug() << "[MNEBeamformerWeights::check_ch_names] data ch names size: " << data_ch_names.size();


    QStringList missing_ch_names;
    for(qint32 i = 0; i < bf_ch_names.size(); ++i)
        if(!data_ch_names.contains(bf_ch_names[i]))
            missing_ch_names.append(bf_ch_names[i]);

    qint32 n_missing = missing_ch_names.size();

    if(n_missing > 0)
    {
        qCritical() << "[MNEBeamformerWeights::check_ch_names] " << n_missing << "channels in beamformer are not present in the data (" << missing_ch_names << ")";
        return false;
    }

 //   qDebug() << " [MNEBeamformerWeights::check_ch_names] Finished checking names.";

    return true;
}



//=============================================================================================================


RowVectorXd MNEBeamformerWeights::compute_frobenius_norm(const MatrixXd matrix){

    RowVectorXd vecFrobNorm = RowVectorXd::Zero(matrix.cols());


    for(int iCol = 0; iCol < matrix.cols(); iCol++){
        double dSumSquares = 0;
        for(int iRow = 0; iRow < matrix.rows(); iRow++){
            dSumSquares += (matrix(iRow,iCol) * matrix(iRow,iCol));
        }
        vecFrobNorm(iCol) = sqrt(dSumSquares);
    }

    return vecFrobNorm;


}


//=============================================================================================================

MatrixXd MNEBeamformerWeights::compute_pseudo_inverse(const MatrixXd &p_matrix, double p_dEpsilon) const
{

    //HINT: similar to function pinv in ft_inverse_lcmv (same as Matlab but default tolerance is twice as high, but without handling of tolerance value as input)


    //dimensions of input matrix
    qint32 nrows = p_matrix.rows();
    qint32 ncols = p_matrix.cols();


    if(nrows != ncols){
        qWarning() << "[MNEBeamformerWeights::compute_pseudo_inverse] Pseudo inverse can only be computed for squared matrices. Return default Matrix.";
        return defaultMatrixXd;
    }


    JacobiSVD<MatrixXd> svd(p_matrix, ComputeFullU | ComputeFullV);
    VectorXd vecSing = svd.singularValues(); //sorted in decreasing order
    MatrixXd matU = svd.matrixU();
    MatrixXd matV = svd.matrixV();

//    std::cout << "[MNEBeamformerWeights::compute_pseudo_inverse] vecSing: " << vecSing << '\n';

    double tol = 10.0 * nrows * vecSing[0] * p_dEpsilon;
    double sumSing = 0;
    for(int i = 0; i < vecSing.size(); i++){
        if(vecSing[i] > tol){
            sumSing += vecSing(i);
        }
    }

//    std::cout << "[MNEBeamformerWeights::compute_pseudo_inverse] tol: " << tol << '\n';
//    std::cout << "[MNEBeamformerWeights::compute_pseudo_inverse] sumSing: " << sumSing << '\n';

    VectorXd vecInvSing = VectorXd::Zero(vecSing.size());
    if(sumSing > 0){
        for(int i = 0; (i < vecSing.size()) && (i < sumSing); i++){
            vecInvSing[i] = 1.0 / vecSing[i];
        }
    }else{
        qWarning() << "[MNEBeamformerWeights::compute_pseudo_inverse] Sum of singular values equals zero. Return default matrix.";
        return defaultMatrixXd;
    }

    //compute the pseudo inverse by recompose svd results
    MatrixXd matInv = matV * vecInvSing.asDiagonal() *  matU.adjoint();


//    //TODO: for debugging only
//    std::cout.precision(17);
//    MatrixXd testMatProd = (p_matrix * matInv * p_matrix);
//    MatrixXd testPinvCond = testMatProd - p_matrix;
//    std::cout << "[MNEBeamformerWeights::compute_pseudo_inverse] Test inv constraint (A A+ A) -A  = 0: " << testPinvCond << "\n";
//    std::cout << "[MNEBeamformerWeights::compute_pseudo_inverse] Test inv diff min: " << testPinvCond.minCoeff() << "\n";
//    std::cout << "[MNEBeamformerWeights::compute_pseudo_inverse] Test inv diff max: " << testPinvCond.maxCoeff() << "\n";


//    qDebug() << "[MNEBeamformerWeights::compute_pseudo_inverse] matInv dim: " << matInv.rows() << " x " << matInv.cols();
//    qDebug() << "[MNEBeamformerWeights::compute_pseudo_inverse] Finished computation of pseudo inverse.";



    return matInv;


}



//=============================================================================================================

MatrixXd MNEBeamformerWeights::invert_data_cov_mat(const FiffCov &p_dataCov)
{
    //invert data covariance matrix
    //HINT: implemented own inversion method similar to mnepy _reg_pinv without regularization since this is already done,

//    qDebug() << "[MNEBeamformerWeights::invert_data_cov_mat] Start inverting data covariance matrix...";


    //decompose matrix with svd
    JacobiSVD<MatrixXd> svd(p_dataCov.data, ComputeFullU | ComputeFullV);
    VectorXd vecSing = svd.singularValues(); //sorted in decreasing order
    MatrixXd matU = svd.matrixU();
    MatrixXd matV = svd.matrixV();


    //invert only non-zero singular values from index 0 to rank-1
    qint32 rank = MNEMath::rank(p_dataCov.data);
    qDebug() << "[MNEBeamformerWeights::invert_data_cov_mat] rank p_dataCov.data = " << rank;

    VectorXd vecSingInv = VectorXd::Zero(vecSing.size());
    for(int i = 0; i < rank; i++){
        if(vecSing[i] > 0){
            vecSingInv[i] = 1 / vecSing[i];
        }else{
            vecSingInv[i] = 0;
        }
    }

    //compute the pseudo inverse by recompose svd results
    //HINT: this should be moore penrose pseude inversion as described in Wikipedia
    MatrixXd tmpMatSingInv = vecSingInv.asDiagonal();
    MatrixXd matCovInv = matV * tmpMatSingInv.transpose() *  matU.adjoint();



//    //TODO: for debugging only
//    std::cout.precision(17);
//    MatrixXd testPinvCond = (p_dataCov.data * matCovInv * p_dataCov.data) - p_dataCov.data;
//    std::cout << "[MNEBeamformerWeights::invert_data_cov_mat] Test inv diff min: " << testPinvCond.minCoeff() << "\n";
//    std::cout << "[MNEBeamformerWeights::invert_data_cov_mat] Test inv diff max: " << testPinvCond.maxCoeff() << "\n";



//    qDebug() << "[MNEBeamformerWeights::invert_data_cov_mat] Finished inversion of data covariance matrix.";

    return matCovInv;

}



//=============================================================================================================


QStringList MNEBeamformerWeights::pick_channels(const QStringList &ch_names,
                                                   const QStringList &include,
                                                   const QStringList &exclude) const
{

    QStringList qSelectedNames;
    for(qint32 i = 0; i < ch_names.size(); i++)
    {
        if((include.size() == 0 || include.contains(ch_names[i])) && !exclude.contains(ch_names[i])){
            qSelectedNames.append(ch_names[i]);
        }
    }

    return qSelectedNames;
}



//=============================================================================================================

QStringList MNEBeamformerWeights::compare_ch_names(const QStringList &chNamesA,
                                                   const QStringList &chNamesB,
                                                   const QStringList &bads) const
{

    //HINT: similar to compare_channel_names form mnepy
    //HINT: code adapted from check_ch_names of inverse operator
    QStringList ch_names; //list of common channels

    for(qint32 i = 0; i < chNamesA.size(); ++i){
        if((!bads.contains(chNamesA[i])) && chNamesB.contains(chNamesA[i])){
            ch_names.append(chNamesA[i]);
        }
    }

    return ch_names;
}


//=============================================================================================================

QStringList MNEBeamformerWeights::check_info_bf(const FiffInfo &p_info,
                                                 const MNEForwardSolution &p_forward,
                                                 const FiffCov &p_data_cov,
                                                 const FiffCov &p_noise_cov) const

{

    //UNDER CONSTRUCTION

    QStringList qCommonChNames; //list of common and good channels

//    QStringList data_names = p_info.ch_names;
//    QStringList fwd_names = p_forward.info.ch_names;
//    QStringList dataCov_names = p_data_cov.names;
//    QStringList noiseCov_names = p_noise_cov.names;


   //idee: fwd als maßstab nehmen, weil da schon bads raussortiert sind,
    //absolut identisch zu _check_info_inv in check.py implementieren

        //compare data to fwd (return all channels in data that are not bad and are present in fwd in ch_names)
    //handle channels from forward model and info (similar to _compare_ch_names (in file check.py))
    qCommonChNames = compare_ch_names(p_info.ch_names,p_forward.info.ch_names,p_info.bads);

//    std::cout << "[MNEBeamformerWeights::check_info_bf] p_info.bads: " << &p_info.bads;

    //make sure no meg reference channels are left
    QStringList qRefCh; // list of reference channels that should be excluded
    fiff_int_t kind;
    for(qint32 i = 0; i < p_info.chs.size(); ++i){
        kind = p_info.chs[i].kind;
        if(kind == FIFFV_REF_MEG_CH){
            qRefCh.append(p_info.ch_names[i]);
        }
    }
    // exclude reference channels from list of common good channels, dont include some channels explicitely
//    qCommonChNames = pick_channels(qCommonChNames,defaultQStringList,qRefCh);
   //exclude all reference channels form list of common good channels
    QStringList commonNoRefs;
    for(qint32 i = 0; i < qCommonChNames.size(); ++i){
        if(!qRefCh.contains(qCommonChNames[i])){
            commonNoRefs.append(qCommonChNames[i]);
        }
    }


    //inform if bad channels of data cov and noise cov do not align with bad channels in data and that all bad channels are to be excluded
    if(p_info.bads != p_data_cov.bads)
    {
        qWarning("List of bad channels of data and data covariance matrix do not match, excluding bad channels from both.");
    }
    if(p_info.bads != p_noise_cov.bads)
    {
        qWarning("List of bad channels of data and noise covariance matrix do not match, excluding bad channels from both.");
    }


    //OLD code here:

//    //compare channels of measurement data, forward solution and both covariance matrices and return list of common good channel indices
//    //HINT: Steps analog to _check_info_inv of mnepy

//    QStringList qCommonChNames; //list of common and good channels



//    //make sure no meg reference channels are left
//    QStringList qRefCh; // list of reference channels that should be excluded
//    fiff_int_t kind;
//    for(qint32 i = 0; i < p_info.chs.size(); ++i){
//        kind = p_info.chs[i].kind;
//        if(kind == FIFFV_REF_MEG_CH){
//            qRefCh.append(p_info.ch_names[i]);
//        }
//    }
//    // exclude reference channels from list of common good channels, dont include some channels explicitely
//    qCommonChNames = pick_channels(qCommonChNames,defaultQStringList,qRefCh);

//    //inform if bad channels of data cov and noise cov do not align with bad channels in data and that all bad channels are to be excluded
//    if(p_info.bads != p_data_cov.bads)
//    {
//        qWarning("List of bad channels of data and data covariance matrix do not match, excluding bad channels from both.");
//    }
//    if(p_info.bads != p_noise_cov.bads)
//    {
//        qWarning("List of bad channels of data and noise covariance matrix do not match, excluding bad channels from both.");
//    }

//    //compare current list of common good channels to data covariance channels (return all channels in current list that are present in data cov and are not bad)
//    qCommonChNames = compare_ch_names(qCommonChNames,p_data_cov.names,p_data_cov.bads);

//    //compare current list of common good channels to noise covariance channels (return all channels in current list that are present in noise cov and are not bad)
//    qCommonChNames = compare_ch_names(qCommonChNames,p_noise_cov.names,p_data_cov.bads);



    return qCommonChNames;
}



//=============================================================================================================

MNEBeamformerWeights MNEBeamformerWeights::make_beamformer_weights(const FiffInfo &p_dataInfo, //the info of the measurement data (no input in fieldtrip)
                                                                   const MNEForwardSolution &p_forward, //the forward solution
                                                                   const FiffCov &p_dataCov, //the data covariance matrix
                                                                   const FiffCov &p_noiseCov, //the noise covariance matrix
                                                                   bool p_bFixedOri,
                                                                   QString p_sWeightnorm, //normalize the beamformer weights ('no' (=unitgain,default), 'unitnoisegain', 'arraygain' or 'nai')
                                                                   qint32 p_iRegParam, //the regularization parameter called lambda in fieldtrip
                                                                   qint32 p_iNAverage //number of averages if beamformer should be applied to evoked data
                                                                   ) {
    //HINT: resembling the FieldTrip code in ft_inverse_lcmv.m




    //TODO: this is only for debugging, delete this later
//    std::ofstream fileDataCovSingIn;
//    std::ofstream fileDataCovSingWhitened;


    qInfo("[MNEBeamformerWeights::make_beamformer_weights] Start calculation of beamformer weights...\n");

    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] p_sWeightnorm: " << p_sWeightnorm;
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] p_bFixedOri: " << p_bFixedOri;
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] p_iRegParam: " << p_iRegParam;
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] p_iNAverage: " << p_iNAverage;

    //prepare output
    MNEBeamformerWeights p_MNEBeamformerWeights;

    //use this logical flags instead of string comparisons
    bool bNormNo = (p_sWeightnorm == "no");
    bool bNormUnitNoise = (p_sWeightnorm == "unitnoisegain");
    bool bNormArray = (p_sWeightnorm == "arraygain");
    bool bNormNai = (p_sWeightnorm == "nai");

    //TODO only for debugging
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] bNormNo: " << bNormNo;
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] bNormUnitNoise: " << bNormUnitNoise;
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] bNormArray: " << bNormArray;
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] bNormNai: " << bNormNai;

    //init local variables
    FiffInfo dataInfo = p_dataInfo;
    MNEForwardSolution forward = p_forward;
    FiffCov dataCov = p_dataCov;
    FiffCov noiseCov = p_noiseCov;

    //TODO: this is only for debugging, delete later (write singular values to file)
//    JacobiSVD<MatrixXd> svdDataCovIn(dataCov.data, ComputeFullU | ComputeFullV);
//    VectorXd vecSingDataCovIn = svdDataCovIn.singularValues(); //sorted in decreasing order
//    fileDataCovSingIn.open("testDataCovSingIn.txt", std::ios::app);
//    for(int i = 0; i < vecSingDataCovIn.size(); i++){
//        fileDataCovSingIn << vecSingDataCovIn[i] << "    ";
//    }
//    fileDataCovSingIn << '\n' << "xxxxxxxxx" << '\n' ;
//    fileDataCovSingIn.close();


//    //TODO: only for debugging, delete later
//    //this creates a spatial filter
//    dataCov.data = MatrixXd::Identity(dataCov.data.rows(),dataCov.data.cols());
//    VectorXd vecSimulatedStdNoise(noiseCov.data.rows());
//    for(int iChan = 0; iChan < vecSimulatedStdNoise.size()-2; iChan += 3){
//        vecSimulatedStdNoise[iChan] = 5e-13; //gradiometer sensor noise std
//        vecSimulatedStdNoise[iChan+1] = 5e-13;
//        vecSimulatedStdNoise[iChan+2] = 20e-15; //magnetometer senosor noise std
//    }
//    MatrixXd matSimulatedNoiseCov = vecSimulatedStdNoise.asDiagonal();
//    noiseCov.data = matSimulatedNoiseCov;


    //TODO only for debugging, delete later
//    noiseCov.data = MatrixXd::Identity(p_noiseCov.data.rows(),p_noiseCov.data.cols());

//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] forward.sol dim: " << forward.sol->nrow << " x " << forward.sol->ncol;
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] dataCov dim: " << dataCov.data.rows() << " x " << dataCov.data.cols();
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] noiseCov dim: " << noiseCov.data.rows() << " x " << noiseCov.data.cols();


    //check invalid parameter combinations (partly copied from make_inverse_operator mnecpp)
    bool is_fwd_fixed_ori = forward.isFixedOrient();
    if(is_fwd_fixed_ori && !p_bFixedOri)
    {
        qWarning("[MNEBeamformerWeights::make_beamformer_weights] Setting p_bFixedOri parameter = true, because the given forward operator has fixed orientation and can only be used to make a fixed-orientation beamformer weights.");
        p_bFixedOri = true;
    }

    //TODO
    //ensure that p_sWeightnorm option no or arraygain because all other options are not debugged yet
    if(!bNormNo && !bNormArray){
        qWarning("[MNEBeamformerWeights::make_beamformer_weights] Weightnorm option not debugged. TODO\n");
    }

    //ensure that p_sWeightnorm option is valid, if not set to default (no weight normalization)
    if(!bNormNo && !bNormUnitNoise && !bNormArray && !bNormNai){
        qWarning("[MNEBeamformerWeights::make_beamformer_weights] Invalid weightnorm option. Set p_sWeightnorm = 'no' (default, baisc unit-gain LCMV beamformer).");
        p_sWeightnorm = "no";
    }

    //ensure that number of averages is positive
    if(p_iNAverage <= 0)
    {
        qWarning("[MNEBeamformerWeights::make_beamformer_weights] The number of averages should be positive. Returned default MNEBeamformerWeights.\n");
        return p_MNEBeamformerWeights;
    }


    qDebug("[MNEBeamformerWeights::make_beamformer_weights] Prepare forward solution...");

    //prepare forward operator: we use the prepare_forward method form MNEForwardSolution
    //HINT: these lines are copied from make_inverse_operator (MNEInverseOperator)
    FiffInfo leadfieldInfo;
    MatrixXd matLeadfield;
    FiffCov outNoiseCov;
    MatrixXd matWhitener;
    qint32 n_nzero;

    //TODO: only for debugging, delete later
//    noiseCov.data = MatrixXd::Identity(noiseCov.data.rows(),noiseCov.data.cols());

    //existing prepare_forward method: selects common channel, prepares noise cov matrix, calculates whitener from noise cov, gain matrix
    //TODO: the computed whitener seems to be the PCA whitener although pca is set to false! Check this.
    forward.prepare_forward(dataInfo, noiseCov, false, leadfieldInfo, matLeadfield, outNoiseCov, matWhitener, n_nzero);


    //TODO: forward rank reduction, if it does not work, delete this here (this makes WT computation very slow in debug mode!!!, check if suited for RT)
//    for(int iPos = 0; iPos < forward.nsource; iPos++){
//        matLeadfield.block(0,iPos*3,matLeadfield.rows(),3) = reduceLocLeadfieldRank(matLeadfield.block(0,iPos*3,matLeadfield.rows(),3));
//    }

//    //TODO: only for debugging. delete later   (this writes prepared forward solution for debugging in Matlab)
//    std::ofstream testFilePreparedForward;
//    testFilePreparedForward.open("testPreparedForward.txt", std::ios::app);
//    for(int iRow = 0; iRow < matLeadfield.rows(); iRow++){

//        for(int iCol = 0; iCol < matLeadfield.cols(); iCol++){

//        testFilePreparedForward << matLeadfield(iRow,iCol) << "    ";

//        }
//        testFilePreparedForward << '\n' ;

//    }
//    testFilePreparedForward << "xxxxxxxxx" << '\n' ;
//    testFilePreparedForward.close();


//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] outNoiseCov.eigvec.adjoint() dim: " << outNoiseCov.eigvec.adjoint().rows() << " x " << outNoiseCov.eigvec.adjoint().cols();
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matWhitener dim: " << matWhitener.rows() << " x " << matWhitener.cols();

//    //TODO: only for debugging, delete later
//    matWhitener = MatrixXd::Identity(matLeadfield.rows(),matLeadfield.rows());

    //TODO: only for debugging, delete later if it does not solve the jumping solution problem
    //Together with the computations in prepare_forward, this creates matWhitener as described in Wester2022 (zca case)
//    matWhitener = outNoiseCov.eigvec * matWhitener;

    //TODO: this is the python non pca whitener (combined with computations in prepare_forward)
//    matWhitener = outNoiseCov.eigvec.adjoint() * matWhitener; //rows of eigvec are eigenvectors

//    //TODO: only for debugging, delete later
    matWhitener = MatrixXd::Identity(matLeadfield.rows(),matLeadfield.rows());


    //Whiten the forward solution with whitener matrix
    matLeadfield = matWhitener * matLeadfield;


//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matLeadfield dim: " << matLeadfield.rows() << " x " << matLeadfield.cols();

    qDebug("[MNEBeamformerWeights::make_beamformer_weights] Finished preparation of forward solution.");


   // ensure that measurement data channels match forward model, noise covariance matrix and data covariance matrix channels (from mnepy make inverse operator, first step there)
//    QStringList lCommonChanNames = p_MNEBeamformerWeights.check_info_bf(dataInfo, forward, dataCov, noiseCov);

    //TODO: only for debugging, delet it if this does not work
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matLeadfield dim: " << matLeadfield.rows() << " x " << matLeadfield.cols();
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] leadfieldInfo.ch_names.size(): " << leadfieldInfo.ch_names.size();
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] outNoiseCov dim: " << outNoiseCov.data.rows() << " x " << outNoiseCov.data.cols();
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matWhitener dim: " << matWhitener.rows() << " x " << matWhitener.cols();


    //check if leadfield and outNoiseCov have the same channel structure
    //TODO: only for debugging, can be deleted afterwards
    qint32 iNumChDifferences = 0;
    for(int i = 0; i < leadfieldInfo.ch_names.size(); ++i){
        if(leadfieldInfo.ch_names[i] != outNoiseCov.names[i]){
            iNumChDifferences++;
        }
    }
    if(iNumChDifferences > 0){
        qWarning() << "[MNEBeamformerWeights::make_beamformer_weights] iNumChDifferences = " << iNumChDifferences;

    }

//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] leadfieldInfo.ch_names : " << leadfieldInfo.ch_names;
//    qDebug() << " XXXXXXXXXXXXXXXXXXXXX";
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] outNoiseCov.names : " << outNoiseCov.names;
//    qDebug() << " XXXXXXXXXXXXXXXXXXXXX";

    QStringList lCommonChanNames = leadfieldInfo.ch_names;
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] Number of common channels = " << lCommonChanNames.size();



    qDebug("[MNEBeamformerWeights::make_beamformer_weights] Prepare data covariance matrix...");

    //TODO: the following part is not debugged yet!
    //estimate noise level in the covariance matrix by the smallest non-zero singular value, always needed for NAI weight normalization#
    qint32 iNoiseLevelDataCov; //only for nai BF relevant
    if(bNormNai){

        qWarning() << "[MNEBeamformerWeights::make_beamformer_weights] Weight normalization method nai not debugged yet.";

        if(MNEMath::rank(dataCov.data) > dataCov.data.size()){
            qWarning("[MNEBeamformerWeights::make_beamformer_weights] Cannot compute a noise subspace with a full-rank covariance matrix and no regularization.");

        }else{
            JacobiSVD<MatrixXd> svd(dataCov.data);
            VectorXd p_sing = svd.singularValues();
            iNoiseLevelDataCov = p_sing[MNEMath::rank(dataCov.data)];
            iNoiseLevelDataCov = std::max(iNoiseLevelDataCov, p_iRegParam);
        }
    }


    //forward solution (+whitener and noise cov matrix) is restricted to sensor channels only (MEG -1, MEG+EE or only EEG) during prepare_forward
    //so we need to reduce the dataCov to the same channels
    FiffCov dataCovPicked =  dataCov.pick_channels(lCommonChanNames);


    //TODO: this is only for debugging, delete later
    //The txt files contais the data covariance matrix estimated from simulated data without noise
//    QString idealDataCovFileName = QCoreApplication::applicationDirPath() + "/mne-cpp-test-data/ideal_noisefree_datacov_simulated_data_isi900_stim340_raw.txt";
//    Eigen::MatrixXd idealDataCov;
//    UTILSLIB::IOUtils::read_eigen_matrix(idealDataCov, idealDataCovFileName);
//    dataCovPicked.data = idealDataCov;

    //TODO: delete this, it makes a spatial filter
//    dataCovPicked.data = matLeadfield * matLeadfield.transpose();


    //whiten data cov mat (after picking because whitener is from picked forward solution) before regularization
//    dataCovPicked.data = matWhitener * dataCovPicked.data;
    //TODO: this is data cov whitening as performed in MNE Python (delete if it does not work)
    dataCovPicked.data = matWhitener * dataCovPicked.data * matWhitener.adjoint();
    dataCovPicked.data = (dataCovPicked.data + dataCovPicked.data.adjoint())/2.0;


//    //TODO: this is only for debugging, delete later
    //TODO: this is only for debugging, delete later (write singular values to file)
//    JacobiSVD<MatrixXd> svdDataCovWhitened(dataCovPicked.data, ComputeFullU | ComputeFullV);
//    VectorXd vecSingDataCovWhitened = svdDataCovWhitened.singularValues(); //sorted in decreasing order
//    fileDataCovSingWhitened.open("testDataCovSingWhitened.txt", std::ios::app);
//    for(int i = 0; i < vecSingDataCovWhitened.size(); i++){
//        fileDataCovSingWhitened << vecSingDataCovWhitened[i] << "    ";
//    }
//    fileDataCovSingWhitened << '\n' << "xxxxxxxxx" << '\n' ;
//    fileDataCovSingWhitened.close();


    //invert data covariance matrix
    MatrixXd matDataCovInv = invert_data_cov_mat(dataCovPicked);
    //check dimension of inverted data cov mat
    if(matDataCovInv.size() != dataCovPicked.data.size()){
        qWarning("[MNEBeamformerWeights::make_beamformer_weights] Size of inverted data covariance matrix does not match size of data covariance matrix.");
    }


    //compute square of inverse data covaricance matrix (needed for unitnoisegain and nai constraint)
    MatrixXd matDataCovInvSquared = matDataCovInv * matDataCovInv;

    qDebug("[MNEBeamformerWeights::make_beamformer_weights] Finished preparation of data covariance matrix.");


    //start scanning the source grid
    qDebug("[MNEBeamformerWeights::make_beamformer_weights] Start scanning the grid of source positions and calculating virtual sensors...");

    MatrixXd matOptimalOri = MatrixXd::Zero(dataInfo.nchan, matLeadfield.cols()/3); //matrix where optimal orientations for fixed orientation beamformers can be stored
    MatrixXd matWT = MatrixXd::Zero(matLeadfield.cols(),matLeadfield.rows()); //matrix of virtual sensors (nsources*3 x nchannels)

    for(int iPos = 0; iPos < forward.nsource; iPos++){ //iterate through all source positions in source grid

        //HINT: print iPos +1 because the forward.nsource value is showed in GUI so counting from 0 to forward.nsource-1 in command window might be irritating
        //internal iteration form 0 to forward.nsources-1
        qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] Scanning source position " << iPos+1 << " of " << forward.nsource;

        //extract local leadfield for source location iPos from leadfield matrix
        //nsource = ncol/3 -> 3 columns per dipole in leadfield matrix
        MatrixXd matLocLeadfield = matLeadfield.block(0,3*iPos,matLeadfield.rows(),3); //extract a block including all rows and only the dipole at iPos (3 columns from idx iPos on)

        //TODO: this is rank reduction of forward solution (recommended for MEG by Westner et al. 2022)
        //TODO: maybe delete this if it does not work, if it works put this in a separate method
        matLocLeadfield = reduceLocLeadfieldRank(matLocLeadfield);
        qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matLocLeadfield dim after rank reduction: " << matLocLeadfield.rows() << " x " << matLocLeadfield.cols();



        //if fixed orientation: calculate optimal orientation for beamformer with respect to weight normalization type
        //HINT: copied from fieldtrip ft_inverse_lcmv
        if(p_bFixedOri){

            qWarning() << "[MNEBeamformerWeights::make_beamformer_weights] p_bFixedOri true. TODO: this is not debugged yet.";

            if(bNormUnitNoise || bNormNai){
                //optimal orientation calculation for unit-noise gain beamformer; also applies nai weightnorm constraint
                //based on equation 4.47 from Sekihara and Nagarajan 2008
                //the following is a reformulation of the generalized eigenproblem

                MatrixXd matTmp = compute_pseudo_inverse((matLocLeadfield.transpose() * matDataCovInvSquared * matLocLeadfield)) * matLocLeadfield.transpose() * matDataCovInv * matLocLeadfield;
                JacobiSVD<MatrixXd> svd(matTmp); //, ComputeThinU ); //TODO:check whether thin is correct here (full in case of square matrix for inversion) thin full defines the dimension of V and U, Check it in whole code
                VectorXd vecSing = svd.singularValues(); //singular values sorted in descending order
                MatrixXd matU = svd.matrixU(); //right singular vectors (TODO?: does this equals the right eigenvectors of EVD?)
                VectorXd vecUnitNoiseOri = matU.block(0,0,matU.rows(),1); //get eigenvector corresponding to biggest eigenvalue
                matLocLeadfield *= vecUnitNoiseOri; //adjust local leadfield to optimal orientation
                matOptimalOri.block(0,iPos,vecUnitNoiseOri.rows(),1) = vecUnitNoiseOri; //store optimal orientation for source position iPos
                //TODO: in fieldtrip the ratio between largest and second larges eigenvalue ist stored for each iPos too, do we need that? if not then we can delete calc of vecSing here (check this for all weightnorm options here)

            }else if(bNormArray){

                //qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] bNormArray true";

                //based on equation 4.44 from Sekihara and Nagarajan 2008
                MatrixXd matTmp = compute_pseudo_inverse((matLocLeadfield.transpose() * matDataCovInv * matLocLeadfield)) * matLocLeadfield.transpose() * matLocLeadfield;
                JacobiSVD<MatrixXd> svd(matTmp); //, ComputeThinU ); //TODO:check whether thin is correct here (full in case of square matrix for inversion) thin full defines the dimension of V and U, Check it in whole code
                VectorXd vecSing = svd.singularValues(); //singular values sorted in descending order
                MatrixXd matU = svd.matrixU(); //right singular vectors (TODO?: does this equals the right eigenvectors of EVD?)
                VectorXd vecArrayGainOri = matU.block(0,0,matU.rows(),1); //get eigenvector corresponding to biggest eigenvalue
                matLocLeadfield *= vecArrayGainOri; //adjust local leadfield to optimal orientation
                matOptimalOri.block(0,iPos,vecArrayGainOri.rows(),1) = vecArrayGainOri; //store optimal orientation for source position iPos
                //TODO: in fieldtrip the ratio between largest and second larges eigenvalue ist stored for each iPos too, do we need that? if not then we can delete calc of vecSing here (check this for all weightnorm options here)

            }else{

                //qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] no weight normalization";


                //compute leadfield for optimal dipole orientation that maximizes beamformer output
                //subsequently the leadfield for only that dipole orientation will be used for the final filter computation
                //in this step the filter computation is not necessary, use the quick way to compute the voxel level covariance (cf. van Veen 1997)
                MatrixXd matTmp = compute_pseudo_inverse((matLocLeadfield.transpose() * matDataCovInv * matLocLeadfield));
                JacobiSVD<MatrixXd> svd(matTmp.real()); //, ComputeThinU ); //TODO:check whether thin is correct here (full in case of square matrix for inversion) thin full defines the dimension of V and U, Check it in whole code
                VectorXd vecSing = svd.singularValues(); //singular values sorted in descending order
                MatrixXd matU = svd.matrixU(); //right singular vectors
                VectorXd vecMaxPowerOri = matU.block(0,0,matU.rows(),1); //get sigular vector corresponding to biggest singular value
                matLocLeadfield *= vecMaxPowerOri; //adjust local leadfield to optimal orientation
                matOptimalOri.block(0,iPos,vecMaxPowerOri.rows(),1) = vecMaxPowerOri; //store optimal orientation for source position iPos
                //TODO: in fieldtrip the ratio between largest and second larges singular value ist stored for each iPos too, do we need that? if not then we can delete calc of vecSing here (check this for all weightnorm options here)

            }//if weightnorm

//            qDebug("[MNEBeamformerWeights::make_beamformer_weights] Calculated optimal orientation.");

        }//if p_bFixedOri


        //construct virtual sensor for source location iPos and add it to W with respect to weight normalization type
        //HINT: copied from fieldtrip ft_inverse_lcmv

        MatrixXd matLocVirtSens = MatrixXd::Zero(3,matLeadfield.rows()); //matrix (3 x nchannels) of local virtual sensor

        if(bNormNai){

            qWarning("[MNEBeamformerWeights::make_beamformer_weights] NAI weight normalization is debugged yet. TODO");
            //qDebug() << "[MNEBeamformerweights::make_beamformer_weights] bNormNai = true";

            // Van Veen's Neural Activity Index
            //HINT: different calculation in Fieldtrip but this is easier to understand; different calculation in Knösche 2022 (check original paper van veen)
            //the scaling term in the denominator is sqrt of projected noise, as per eqn. 2.67 of Sekihara & Nagarajan 2008

            // Van Veen's Neural Activity Index
            // below equation is equivalent to following:
            // filt = pinv(lf' * invC * lf) * lf' * invC;
            // filt = filt/sqrt(noise*filt*filt');
            // the scaling term in the denominator is sqrt of projected noise, as per eqn. 2.67 of Sekihara & Nagarajan 2008 (S&N)

            if(p_bFixedOri){

                matLocVirtSens = compute_pseudo_inverse((iNoiseLevelDataCov * matLocLeadfield.transpose() * matDataCovInvSquared * matLocLeadfield).array().sqrt()) * matLocLeadfield.transpose() * matDataCovInv; // based on S&N eqn. 4.08
                matWT.block(iPos*3,0,3,matLocVirtSens.cols()) = matLocVirtSens; //store local virtual sensor with dimension (3 x nchannels) in filter weights matrix matW
                // iPos * 3 is row index where to write the virtual sensor (containing of 3 row vectors) for source position iPos into matW
            }else{
                qWarning("[MNEBeamformerWeights::make_beamformer_weights] Vector version of NAI weight normalization is not implemented TODO");
            }

        }else if(bNormUnitNoise){

            qWarning("[MNEBeamformerWeights::make_beamformer_weights] Unit noise gain weight normalization is debugged yet. TODO");

            // filt*filt' = I
            // Unit-noise gain minimum variance (aka Borgiotti-Kaplan) beamformer
            // below equation is equivalent to following:
            // filt = pinv(lf' * invC * lf) * lf' * invC;
            // filt = filt/sqrt(filt*filt');
            if(p_bFixedOri){
                matLocVirtSens = compute_pseudo_inverse((matLocLeadfield.transpose() * matDataCovInvSquared * matLocLeadfield).array().sqrt()) * matLocLeadfield.transpose(); // S&N eqn. 4.15
                matWT.block(iPos*3,0,3,matLocVirtSens.cols()) = matLocVirtSens; //local virtual sensor with dimension (3 x nchannels)
            }else{
                // compute the matrix that is used for scaling of the filter's rows, as per eqn. 4.83
                MatrixXd matDenom = compute_pseudo_inverse(matLocLeadfield.transpose() * matDataCovInv * matLocLeadfield);
                MatrixXd matGamma = matDenom * (matLocLeadfield.transpose() * matDataCovInvSquared * matLocLeadfield) * matDenom;
                // compute the spatial filter, as per eqn. 4.85
                matLocVirtSens = (1/(matGamma.diagonal()).array().sqrt()).matrix().asDiagonal() * matDenom * matLocLeadfield.transpose() * matDataCovInv; //see notes 07.02.2023
                matWT.block(iPos*3,0,3,matLocVirtSens.cols()) = matLocVirtSens; //local virtual sensor with dimension (3 x nchannels)
            }

        }else if(bNormArray){


            if(p_bFixedOri){
                qWarning("[MNEBeamformerWeights::make_beamformer_weights] Scalar version of array gain beamformer is not implemented TODO");


            }else{

                // filt*lf = ||lf||, applies to scalar leadfield, and to one of the possibilities of the vector version, eqn. 4.75
                //TODO: Frobenius norm is better s. Wester 2022 so norm() is used which returns Frobenius norm of matrix
                //normalize local leadfield to Frobenius Norm (s. Westner 2022, why Frobenius is better than L2 as used in Fieldtrip)
                MatrixXd matLocLeadfieldNorm = MatrixXd::Zero(matLocLeadfield.rows(),matLocLeadfield.cols());
                RowVectorXd vecFrobNorm = compute_frobenius_norm(matLocLeadfield);

                for(int i = 0; i < vecFrobNorm.size(); i++){
                    matLocLeadfieldNorm.col(i) = matLocLeadfield.col(i) / vecFrobNorm(i);
                }

                matLocVirtSens = compute_pseudo_inverse(matLocLeadfieldNorm.transpose() * matDataCovInv * matLocLeadfieldNorm) * matLocLeadfieldNorm.transpose() * matDataCovInv; // S&N eqn. 4.09 (scalar version), and eqn. 4.75 (vector version)
                matWT.block(iPos*3,0,3,matLocVirtSens.cols()) = matLocVirtSens; //store local virtual sensor with dimension (3 x nchannels) in filter weights matrix matWT

//                //test array gain beamformer condition
//                MatrixXd matTestBFCond = matLocVirtSens * matLocLeadfieldNorm;
//                std::cout.precision(18);
//                std::cout << "[MNEBeamformerWeights::make_beamformer_weights] matTestBFCond arraygain should be identity matrix: " << matTestBFCond << '\n';
            }

        }else{ //unitgain BF without weight normalization


            // this is the 'standard' unit gain constraint spatial filter: filt*lf=I, applies both to vector and scalar leadfields
            //HINT: compute_pseudo_inverse should work correctly

            //TODO only for debugging, delete later
//            matLocVirtSens = compute_pseudo_inverse(matLocLeadfield.transpose() * matLocLeadfield) * matLocLeadfield.transpose();




            matLocVirtSens = compute_pseudo_inverse(matLocLeadfield.transpose() * matDataCovInv * matLocLeadfield) * matLocLeadfield.transpose() * matDataCovInv;
            matWT.block(iPos*3,0,3,matLocVirtSens.cols()) = matLocVirtSens; //store local virtual sensor with dimension (3 x nchannels) in filter weights matrix matWT


//            qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matLocVirtSens dim: " << matLocVirtSens.rows() << " x " << matLocVirtSens.cols();

//            MatrixXd matTestBFCond = matLocVirtSens * matLocLeadfield;
//            std::cout.precision(18);
//            std::cout << "[MNEBeamformerWeights::make_beamformer_weights] matTestBFCond should be identity matrix: " << matTestBFCond << '\n';

        }

    }//for iPos (scanning the source grid)

    qDebug("[MNEBeamformerWeights::make_beamformer_weights] Finished scanning the grid of source positions and calculating virtual sensors for each position.");

//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matWT mean: " << matWT.mean();

//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matWT dim: " << matWT.rows() << " x " << matWT.cols();
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] outNoiseCov dim: " << outNoiseCov.data.rows() << " x " << outNoiseCov.data.cols();
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] dataCovPicked.data dim: " << dataCovPicked.data.rows() << " x " << dataCovPicked.data.cols();

//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] forward.nsources: " << forward.nsource;
//    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] lCommonChanNames.size(): " << lCommonChanNames.size();

    //for source reconstruction accuracy evaluation
//    std::ofstream fileMakeBeamformerWeights;
//    fileMakeBeamformerWeights.open("testfileMakeBeamformerWeights.txt", std::ios::app);
//    fileMakeBeamformerWeights <<  "  " << '\n';
//    fileMakeBeamformerWeights.close();

    //TODO: only for debugging delete later
    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] p_MNEBeamformerWeights.weightNorm: " << p_MNEBeamformerWeights.weightNorm;



    //store filter weights and additional info
    p_MNEBeamformerWeights.info = leadfieldInfo;
    p_MNEBeamformerWeights.weights = matWT;
//        p_MNEBeamformerWeights.weights = MatrixXd::Ones(matWT.rows(),matWT.cols());
    p_MNEBeamformerWeights.data_cov = FiffCov::SDPtr(new FiffCov(dataCovPicked));
    p_MNEBeamformerWeights.noise_cov = FiffCov::SDPtr(new FiffCov(outNoiseCov));
    p_MNEBeamformerWeights.weightNorm = p_sWeightnorm;
    p_MNEBeamformerWeights.whitener = matWhitener;
    p_MNEBeamformerWeights.fixedOri = p_bFixedOri;
    p_MNEBeamformerWeights.optOri = matOptimalOri;
    p_MNEBeamformerWeights.nsource = forward.nsource;
    p_MNEBeamformerWeights.nchan = lCommonChanNames.size();
    p_MNEBeamformerWeights.projs = dataInfo.projs;
    p_MNEBeamformerWeights.src = forward.src;
    p_MNEBeamformerWeights.nave = p_iNAverage;

    qDebug("[MNEBeamformerWeights::make_beamformer_weights] Finished calculation of beamformer weights.");

    return p_MNEBeamformerWeights;

}
