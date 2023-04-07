//=============================================================================================================
/**
 * @file     mne_beamformer_weights.cpp
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
 * @brief    MNEBeamformerWeights class definition.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "mne_beamformer_weights.h"

#include <math.h>


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
: weights(defaultMatrixXd)
, data_cov(new FiffCov)
, noise_cov(new FiffCov)
, whitener(defaultMatrixXd)
, fixedOri(false)
, optOri(defaultMatrixXd)
, nsource(-1)
, nchan(-1)
, sourcePowEst(defaultVectorXd)
, noisePowEst(defaultVectorXd)
, nave(-1)
{
    qDebug() << "[MNEBeamformerWeights::MNEBeamformerWeights] TEST 1";
    qRegisterMetaType<QSharedPointer<MNELIB::MNEBeamformerWeights> >("QSharedPointer<MNELIB::MNEBeamformerWeights>");
    qRegisterMetaType<MNELIB::MNEBeamformerWeights>("MNELIB::MNEBeamformerWeights");
}



//=============================================================================================================


MNEBeamformerWeights::MNEBeamformerWeights(const FiffInfo &p_dataInfo,
                                           const MNEForwardSolution& p_forward,
                                           const FiffCov& p_dataCov,
                                           const FiffCov& p_noiseCov,
                                           QString p_sWeightnorm,
                                           QString p_sPowMethod,
                                           bool p_bFixedOri,
                                           bool p_bEstNoisePow,
                                           bool p_bProjectMom,                                           
                                           qint32 p_iRegParam,
                                           qint32 p_iNAverages)
{

    qDebug() << "[MNEBeamformerWeights::MNEBeamformerWeights] TEST 2";

    *this = MNEBeamformerWeights::make_beamformer_weights(p_dataInfo,p_forward,p_dataCov,p_noiseCov,p_bFixedOri,p_bEstNoisePow,p_bProjectMom,p_sWeightnorm,p_iRegParam,p_iNAverages);
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
    , sourcePowEst(p_MNEBeamformerWeights.sourcePowEst)
    , noisePowEst(p_MNEBeamformerWeights.noisePowEst)
    , projs(p_MNEBeamformerWeights.projs)
    , proj(p_MNEBeamformerWeights.proj)
    , src(p_MNEBeamformerWeights.src)
    , nave(p_MNEBeamformerWeights.nave)

{
        qDebug() << "[MNEBeamformerWeights::MNEBeamformerWeights] TEST 3";

    qRegisterMetaType<QSharedPointer<MNELIB::MNEBeamformerWeights> >("QSharedPointer<MNELIB::MNEBeamformerWeights>");
    qRegisterMetaType<MNELIB::MNEBeamformerWeights>("MNELIB::MNEBeamformerWeights");

}



//=============================================================================================================

MNEBeamformerWeights::~MNEBeamformerWeights()

{
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
    //TODO maybe move this to mne math in the end
    //HINT: similar to function pinv in ft_inverse_lcmv (same as Matlab but default tolerance is twice as high, but without handling of tolerance value as input
    //HINT: used this version for calculation with Eigen from GitHub https://gist.github.com/pshriwise/67c2ae78e5db3831da38390a8b2a209f, adapted it a bit

    qDebug() << "[MNEBeamformerWeights::compute_pseudo_inverse] Start computing pseudo inverse...";

    //dimensions of input matrix
    qint32 nrows = p_matrix.rows();
    qint32 ncols = p_matrix.cols();

    //qDebug() << "[MNEBeamformerWeights::compute_pseudo_inverse] nrows = " << nrows;
    //ug() << "[MNEBeamformerWeights::compute_pseudo_inverse] ncols = " << ncols;

    if(nrows != ncols){
        qWarning() << "[MNEBeamformerWeights::compute_pseudo_inverse] Pseudo inverse can only be computed for squared matrices. TODO.";
    }

    JacobiSVD<MatrixXd> svd(p_matrix ,ComputeFullU | ComputeFullV);


    //calculate tolerance
    double tolerance = 10 * std::max(ncols, nrows) * svd.singularValues().array().abs()(0) * p_dEpsilon;

    // return inverse (U is real in case p_matrix is real -> adjoint of U = U^T; changed this here)
    //MatrixXd matInv = svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
    MatrixXd matInv = svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().transpose();


    qDebug() << "[MNEBeamformerWeights::compute_pseudo_inverse] matInv dim: " << matInv.rows() << " x " << matInv.cols();
    qDebug() << "[MNEBeamformerWeights::compute_pseudo_inverse] Finished computation of pseudo inverse.";

    return matInv;


/*    //original code from GitHub
 *
 *    template<typename _Matrix_Type_>
    _Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon())
    {

        Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeFullU | Eigen::ComputeFullV);
            // For a non-square matrix
            // Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
        double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
        return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
    }
*/



}



//=============================================================================================================

MatrixXd MNEBeamformerWeights::invert_data_cov_mat(const FiffCov &p_dataCov)
{
    //invert data covariance matrix
    //HINT: implemented own inversion method similar to mnepy _reg_pinv without regularization since this is already done,

    qDebug() << "[MNEBeamformerWeights::invert_data_cov_mat] Start inverting data covariance matrix...";


    /* mnepy docu of _reg_pinv in numerics.py
    Compute a regularized pseudoinverse of Hermitian matrices.
    Regularization is performed by adding a constant value to each diagonal
    element of the matrix before inversion. This is known as "diagonal
    loading". The loading factor is computed as ``reg * np.trace(x) / len(x)``.

    The pseudo-inverse is computed through SVD decomposition and inverting the
    singular values. When the matrix is rank deficient, some singular values
    will be close to zero and will not be used during the inversion. The number
    of singular values to use can either be manually specified or automatically
    estimated.
    */

    //TODO?: ensure that cov mat is square?, all cov mat are square
    //TODO?: ensure that input matrix is Hermitian (symmetric), all covariance matrices of real valued random variables are hermitian
    //maybe dont check these things because they are true by definition

//    qDebug() << "[MNEBeamformerWeights::invert_data_cov_mat] p_dataCov.data dim: " << p_dataCov.data.rows() << " x " << p_dataCov.data.cols();


    //decompose matrix with svd
    JacobiSVD<MatrixXd> svd(p_dataCov.data, ComputeFullU | ComputeFullV); //TODO: check whether we need full U and V here (full/thin defines dimension of U and V)
    VectorXd vecSing = svd.singularValues(); //sorted in decreasing order
    MatrixXd matU = svd.matrixU();
    MatrixXd matV = svd.matrixV();

/*    qDebug() << "[MNEBeamformerWeights::invert_data_cov_mat] vecSing size: " << vecSing.size();
    qDebug() << "[MNEBeamformerWeights::invert_data_cov_mat] matU dim: " << matU.rows() << " x " << matU.cols();
    qDebug() << "[MNEBeamformerWeights::invert_data_cov_mat] matV dim: " << matV.rows() << " x " << matV.cols();
*/

    //invert only non-zero singular values (invert only singular values with index
    qint32 rank = MNEMath::rank(p_dataCov.data); //dont know whether this rank calculation is too easy here (no differentiation between channel types)
    qDebug() << "[MNEBeamformerWeights::invert_data_cov_mat] rank p_dataCov.data = " << rank;
    VectorXd vecSingInv = VectorXd::Zero(vecSing.size());
    for(int i = 0; i < (p_dataCov.dim - rank); i++){
        vecSingInv[i] = 1 / vecSing[i];
    }

    //qDebug() << "[MNEBeamformerWeights::invert_data_cov_mat] vecSingInv size: " << vecSingInv.size();


    //compute the pseudo inverse by recompose svd results
    MatrixXd matCovInv = matV * vecSingInv.asDiagonal() *  matU;

    qDebug() << "[MNEBeamformerWeights::invert_data_cov_mat] Finished inversion of data covariance matrix.";

    return matCovInv;

}



//=============================================================================================================


QStringList MNEBeamformerWeights::pick_channels(const QStringList &ch_names,
                                                   const QStringList &include,
                                                   const QStringList &exclude) const
{
    //TODO maybe move this to fiff library in the end
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
    //compare channels of measurement data, forward solution and both covariance matrices and return row vector of common good channel indices
    //HINT: Steps analog to _check_info_inv of mnepy

    QStringList qCommonChNames; //list of common and good channels

    //compare data to fwd (return all channels in data that are not bad and are present in fwd in ch_names)
    qCommonChNames = compare_ch_names(p_info.ch_names,p_forward.info.ch_names,p_info.bads);

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
    qCommonChNames = pick_channels(qCommonChNames,defaultQStringList,qRefCh);

    //inform if bad channels of data cov and noise cov do not align with bad channels in data and that all bad channels are to be excluded
    if(p_info.bads != p_data_cov.bads)
    {
        qWarning("List of bad channels of data and data covariance matrix do not match, excluding bad channels from both.");
    }
    if(p_info.bads != p_noise_cov.bads)
    {
        qWarning("List of bad channels of data and noise covariance matrix do not match, excluding bad channels from both.");
    }

    //compare current list of common good channels to data covariance channels (return all channels in current list that are present in data cov and are not bad)
    qCommonChNames = compare_ch_names(qCommonChNames,p_data_cov.names,p_data_cov.bads);

    //compare current list of common good channels to noise covariance channels (return all channels in current list that are present in noise cov and are not bad)
    qCommonChNames = compare_ch_names(qCommonChNames,p_noise_cov.names,p_data_cov.bads);



    return qCommonChNames;
}



//=============================================================================================================

MNEBeamformerWeights MNEBeamformerWeights::make_beamformer_weights(const FiffInfo &p_dataInfo, //the info of the measurement data (no input in fieldtrip)
                                                                   const MNEForwardSolution &p_forward, //the forward solution
                                                                   const FiffCov &p_dataCov, //the data covariance matrix
                                                                   const FiffCov &p_noiseCov, //can be 'trace' (default s. Van Veen 1997) or 'lambda1'
                                                                   bool p_bFixedOri, // true for fixed orientation, false for free orientation (default=false)
                                                                   bool p_bEstNoisePow, // estimate noise power projected through filter (default=true)
                                                                   bool p_bProjectMom, // true: project the dipole moment time course on the direction of maximal power (default=false), TODO: check what this is for, when do we need it and why?
                                                                   QString p_sWeightnorm, //normalize the beamformer weights ('no' (=unitgain,default), 'unitnoisegain', 'arraygain' or 'nai')
                                                                   qint32 p_iRegParam, //the regularization parameter //called lambda in fieldtrip
                                                                   qint32 p_iNAverage //number of averages if beamformer should be applied to evoked data
                                                                   ) {
    //HINT: resembling the FieldTrip code in ft_inverse_lcmv.m
    //because mnepy code is highly modular and too complex for first lcmv implementation
    //decided on FieldTrip implementation because it is easier to understand (maybe less overload code)
    // TODO: add to outlook: for better symmetry to exsisting code in mnecpp, restructure code according to mnepy implementation

    //TODO?:split this method into prepare_beamformer_input and make beamformer weights

    //TODO: check whether mnepy includes some useful steps we need here

    //TODO: option for rank reduction (s. Westner 2022) since for MEG it is recommende to reduce the rank of the forward field

    //TODO: check all equations from Sekihara Nagarajan 2008 (maybe in Westner 2022 and Knösche 2022 so no need of book)


    qInfo("[MNEBeamformerWeights::make_beamformer_weights] Start calculation of beamformer weights...\n");


    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] p_sWeightnorm: " << p_sWeightnorm;

    //prepare output
    MNEBeamformerWeights p_MNEBeamformerWeights;

    //use this logical flags instead of string comparisons (idea from fieldtrip ft_inverse_lcmv, saves time during scanning loop)
    bool bNormNo = (p_sWeightnorm == "no");
    bool bNormUnitNoise = (p_sWeightnorm == "unitnoisegain");
    bool bNormArray = (p_sWeightnorm == "arraygain");
    bool bNormNai = (p_sWeightnorm == "nai");

    //init local variables (because parameters need to be const)
    FiffInfo dataInfo = p_dataInfo;
    MNEForwardSolution forward = p_forward;
    FiffCov dataCov = p_dataCov;
    FiffCov noiseCov = p_noiseCov;


    //check invalid parameter combinations (partly copied from make_inverse_operator mnecpp)
    bool is_fixed_ori = forward.isFixedOrient();
    std::cout << "[MNEBeamformerWeights::make_beamformer_weights] TODO: do surf_ori check\n"; //TODO: this warning was copied from make_inverse_operator, check whether we need it here too
    if(is_fixed_ori && !p_bFixedOri)
    {
        qWarning("[MNEBeamformerWeights::make_beamformer_weights] Setting p_bFixedOri parameter = true. Because the given forward operator has fixed orientation and can only be used to make a fixed-orientation beamformer weights.\n");
        p_bFixedOri = true;
    }

    //TODO
    //ensure that p_sWeightnorm option no because all other options are not debugged yet
    if(!bNormNo){
        qWarning("[MNEBeamformerWeights::make_beamformer_weights] Weightnorm option not debugged. TODO\n");
        //p_sWeightnorm = "no";
    }

    //ensure that p_sWeightnorm option is valid, if not set to default (no weight normalization)
    if(!bNormNo && !bNormUnitNoise && !bNormArray && !bNormNai){
        qWarning("[MNEBeamformerWeights::make_beamformer_weights] Invalid weightnorm option. Set p_sWeightnorm = 'no' (default) and compute a unit-gain beamformer without weight normalization.\n");
        p_sWeightnorm = "no";
    }

    //ensure that number of averages is positive
    if(p_iNAverage <= 0)
    {
        qWarning("[MNEBeamformerWeights::make_beamformer_weights] The number of averages should be positive. Returned default MNEBeamformerWeights.\n");
        return p_MNEBeamformerWeights;
    }


    //TODO: do we need to ensure that there is only one channel type? check whether fieldtrip code is only for eeg or meg but not mixed channels

    qInfo("[MNEBeamformerWeights::make_beamformer_weights] Prepare measurement info...\n");

    // ensure that measurement data channels match forward model, noise covariance matrix and data covariance matrix channels (from mnepy make inverse operator, first step there)
    //TODO: channel selection stuff is performerd in miminumnorms calculate Inverse, move it?
    //TODO: calcFiffInfo should ensure that we have the same channels in all four instances, maybe leave this as double check or delete it and trust the calcFiff Implementation? (double check would be better I guess but it means a bit code overload)
    // return the channels that are common to all four objects
    //from mnepy
    QStringList lCommonChanNames = p_MNEBeamformerWeights.check_info_bf(dataInfo, forward, dataCov, noiseCov);
    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] Number of common channels = " << lCommonChanNames.size();

    //restrict data info to common channels
    RowVectorXi picksCommonChan; //indices of channels that are common to data, forward, and both covariance matrices
    picksCommonChan = FiffInfoBase::pick_channels(lCommonChanNames);
    dataInfo.pick_info(picksCommonChan);

    qInfo("[MNEBeamformerWeights::make_beamformer_weights] Finished preparation of measurement info.\n");

    //scaling with number of averages (necessary if W should be applied to averaged data later on, covariance matrices need to be scaled here because they are fundamental for W computation (e.g. whitener depends on it))
    //HINT: copied from prepare_inverse_operator and adapted for the purposes here
    //TODO: check whether we need this scaling her
    if(p_iNAverage != 1){

        qInfo("[MNEBeamformerWeights::make_beamformer_weights] Scale data and noise covariance matrix...\n");

        float fScale     = 1/((float)p_iNAverage); //scaling factor
        noiseCov.data  *= fScale; //scaling data and eigenvalues of noise covariance matrix
        noiseCov.eig   *= fScale;
        dataCov.data *= fScale; //scaling data and eigenvalues of data covariance matrix
        dataCov.eig *= fScale; //TODO: check whether we need to scale the eigenvalues of data covariance matrix (are they used somewhere else? do we need them scaled for consistency)

        qInfo("[MNEBeamformerWeights::make_beamformer_weights] Scaled noise and data covariance with scaling factor = %f (number of averages = %d)\n",fScale,p_iNAverage);

    }

    qInfo("[MNEBeamformerWeights::make_beamformer_weights] Prepare forward solution...\n");

    //prepare forward operator: we use the prepare_forward method form MNEForwardSolution
    //existing prepare_forward method: selects common channel, prepares noise cov matrix, calculates whitener from noise cov, gain matrix
    //HINT: these lines are copied from make_inverse_operator
    FiffInfo leadfieldInfo;
    MatrixXd matLeadfield;
    FiffCov outNoiseCov;
    MatrixXd matWhitener; //necessary for whitening of data covariance matrix and forward solution
    qint32 n_nzero; //total rank of ??? noise covariance matrix? leadfield? whitener? all of them? (used in make_inverse_operator to scale the source covariance matrix) maybe we dont need that value here
    forward.prepare_forward(dataInfo, noiseCov, false, leadfieldInfo, matLeadfield, outNoiseCov, matWhitener, n_nzero);

    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matLeadfield dim: " << matLeadfield.rows() << " x " << matLeadfield.cols();
    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] outNoiseCov dim: " <<  outNoiseCov.data.rows() << " x " << outNoiseCov.data.cols();
    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matWhitener dim: " << matWhitener.rows() << " x " << matWhitener.cols();

    //Whiten the forward solution with whitener
    matLeadfield = matWhitener * matLeadfield;

    qInfo("[MNEBeamformerWeights::make_beamformer_weights] Finished preparation of forward solution.\n");


    qInfo("[MNEBeamformerWeights::make_beamformer_weights] Prepare data covariance matrix...\n");

    //estimate noise level in the covariance matrix by the smallest non-zero singular value, always needed for NAI weight normalization
    //TODO: maybe we can do this step below? e.g. below regularization does not work because else statement uses cov data which is modified during regularization
    qint32 iNoiseLevelDataCov; //only for nai and unitnoisegain BF relevant
    if(bNormNai || p_bEstNoisePow){
        // MNEMath::rank method does not include that there are different sensor type with different amplitude ranges
        // but fieldtrip uses matlabs rank function which does not include this problem either (matlab svd returns singular values in descending order)
        if(MNEMath::rank(dataCov.data) > dataCov.data.size()){ //TODO:check whether this is implemented correct. c++ size() == python len()?, rank caluclated independently from channel types?
            //TODO:  add if no regularization for this warning
            qWarning("[MNEBeamformerWeights::make_beamformer_weights] Cannot compute a noise subspace with a full-rank covariance matrix and no regularization.\n");
            //TODO:end if no regularization
            //TODO: else set noise level to regularization parameter
        }else{
            JacobiSVD<MatrixXd> svd(dataCov.data);
            VectorXd p_sing = svd.singularValues();
            //TODO: maybe use the eig from FiffCov instance here instead of computing the singular values
            iNoiseLevelDataCov = p_sing[MNEMath::rank(dataCov.data)]; //denoted as noise in fieldtrip,
            iNoiseLevelDataCov = std::max(iNoiseLevelDataCov, p_iRegParam);
        }
    } //end if nai or unitnoisegain

    //prepare data covariance matrix for inversion
    //restrict data covariance matrix to common channels
    //TODO: forward solution (+whitener and noise cov matrix)is restricted to sensor channels only (MEG -1 = 305, MEG+EE or only EEG) during prepare_forward
    //so we need to reduce the dataCov to the same channels
    //old: dataCov.pick_channels(lCommonChanNames);
    FiffCov dataCovPicked =  dataCov.pick_channels(outNoiseCov.names);


    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] dataCovPicked.data dim: " << dataCovPicked.data.rows() << " x " << dataCovPicked.data.cols();

    //whiten data cov mat (after picking because whitener is from picked forward solution) before regularization
    dataCovPicked.data = matWhitener * dataCovPicked.data;


    //HINT: fieldTrip does regularization of data cov here
    //not necessary, because regularization is performerd during creation of data covariance matrix


    //invert data covariance matrix
    MatrixXd matDataCovInv = p_MNEBeamformerWeights.invert_data_cov_mat(dataCovPicked);
    //check dimension of inverted data cov mat
    if(matDataCovInv.size() != dataCovPicked.data.size()){
        qWarning("[MNEBeamformerWeights::make_beamformer_weights] Size of inverted data covariance matrix does not match size of data covariance matrix.\n");
    }


    //compute square of inverse data covaricance matrix (needed for unitnoisegain and nai constraint)
    MatrixXd matDataCovInvSquared = matDataCovInv * matDataCovInv;


    qInfo("[MNEBeamformerWeights::make_beamformer_weights] Finished preparation of data covariance matrix.\n");


    //TODO: maybe from check_info_bf up to here in prepare_beamformer_input?


    //start scanning the source grid (LCMV BF uses local leadfield for computation of filter weights)
    qInfo("[MNEBeamformerWeights::make_beamformer_weights] Start scanning the grid of source positions and calculating virtual sensors...\n");

    //matrix where optimal orientations for fixed orientation beamformers can be stored
    //TODO make this a FiffNamedMatrix to preserve channel names (we need new list of column names then)
    MatrixXd matOptimalOri = MatrixXd::Zero(dataInfo.nchan, matLeadfield.cols()/3);
    MatrixXd matWT = MatrixXd::Zero(matLeadfield.cols(),matLeadfield.rows()); //matrix of virtual sensors (nsources*3 x nchannels)
    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matWT definition dim: " << matWT.rows() << " x " << matWT.cols();

    VectorXd vecSourcePow = VectorXd::Zero(matLeadfield.cols()/3,1); //TODO check whether dimension is correct here; vector of activity strength (power) for each source position
    VectorXd vecNoisePow = VectorXd::Zero(matLeadfield.cols()/3,1); //vector of estimated noise power that is projected through the filter for each source position

    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] vecSourcePow size: " << vecSourcePow.size();

    for(int iPos = 0; iPos < forward.nsource; iPos++){ //iterate through all source positions in source grid

        //HINT: print iPos +1 because the forward.nsource value is showed in GUI so counting from 0 to forward.nsource-1 in command window might be irritating
        //internal iteration form 0 to forward.nsources-1
        printf("Scanning source position %i/%i\n", iPos+1, forward.nsource); //count the source positions to show progress

        //extract local leadfield for source location iPos from leadfield matrix
        //nsource = ncol/3 -> 3 columns per dipole in leadfield matrix
        //qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] iPos = " << iPos;
        //qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] forward.nchan = " << forward.nchan;

        //TODO: check if number 3 should be a variable (does it depend on fixed, free orientation?)
        MatrixXd matLocLeadfield = matLeadfield.block(0,iPos,matLeadfield.rows(),3); //extract a block including all rows and only the dipole at iPos (3 columns from idx iPos on)

        qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matLocLeadfield dim: " << matLocLeadfield.rows() << " x " << matLocLeadfield.cols();

        //if fixed orientation: calculate optimal orientation for beamformer with respect to weight normalization type
        //HINT: copied from fieldtrip ft_inverse_lcmv
        //TODO: check all these equations form the book

        qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] p_bFixedOri = " << p_bFixedOri;

        if(p_bFixedOri){

            qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] p_bFixedOri true. TODO: this is not debugged yet.";

            if(bNormUnitNoise || bNormNai){ //HINT: dont use switch case here because it does not work for string cases
                //optimal orientation calculation for unit-noise gain beamformer; also applies nai weightnorm constraint
                //based on equation 4.47 from Sekihara and Nagarajan 2008
                //the following is a reformulation of the generalized eigenproblem

                qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] bNormUnitNoise || bNormNai true";

                MatrixXd matTmp = compute_pseudo_inverse((matLocLeadfield.transpose() * matDataCovInvSquared * matLocLeadfield)) * matLocLeadfield.transpose() * matDataCovInv * matLocLeadfield;
                JacobiSVD<MatrixXd> svd(matTmp); //, ComputeThinU ); //TODO:check whether thin is correct here (full in case of square matrix for inversion) thin full defines the dimension of V and U, Check it in whole code
                VectorXd vecSing = svd.singularValues(); //singular values sorted in descending order
                MatrixXd matU = svd.matrixU(); //right singular vectors (TODO?: does this equals the right eigenvectors of EVD?)
                VectorXd vecUnitNoiseOri = matU.block(0,0,matU.rows(),1); //get eigenvector corresponding to biggest eigenvalue
                matLocLeadfield *= vecUnitNoiseOri; //adjust local leadfield to optimal orientation
                matOptimalOri.block(0,iPos,vecUnitNoiseOri.rows(),1) = vecUnitNoiseOri; //store optimal orientation for source position iPos
                //TODO: in fieldtrip the ratio between largest and second larges eigenvalue ist stored for each iPos too, do we need that? if not then we can delete calc of vecSing here (check this for all weightnorm options here)

            }else if(bNormArray){

                qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] bNormArray true";


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

                qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] no weight normalization";


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

            qInfo("[MNEBeamformerWeights::make_beamformer_weights] Calculated optimal orientation.\n");

        }//if p_bFixedOri


        //construct virtual sensor for source location iPos and add it to W with respect to weight normalization type
        //HINT: copied from fieldtrip ft_inverse_lcmv

        MatrixXd matLocVirtSens; //matrix (3 x nchannels) of local virtual sensor

        if(bNormNai){

            qDebug() << "[MNEBeamformerweights::make_beamformer_weights] bNormNai = true";

            // Van Veen's Neural Activity Index
            //HINT: different calculation in Fieldtrip but this is easier to understand; different calculation in Knösche 2022 (check original paper van veen)
            //the scaling term in the denominator is sqrt of projected noise, as per eqn. 2.67 of Sekihara & Nagarajan 2008

            // Van Veen's Neural Activity Index
            // below equation is equivalent to following:
            // filt = pinv(lf' * invC * lf) * lf' * invC;
            // filt = filt/sqrt(noise*filt*filt');
            // the scaling term in the denominator is sqrt of projected noise, as per eqn. 2.67 of Sekihara & Nagarajan 2008 (S&N)

            if(p_bFixedOri){
                //TODO: check whether this cwise sqrt operation works correct here (was coincindence that I tried this version of sqrt application)
                matLocVirtSens = compute_pseudo_inverse((iNoiseLevelDataCov * matLocLeadfield.transpose() * matDataCovInvSquared * matLocLeadfield).array().sqrt()) * matLocLeadfield.transpose() * matDataCovInv; // based on S&N eqn. 4.08
                matWT.block(iPos*3,0,3,matLocVirtSens.cols()) = matLocVirtSens; //store local virtual sensor with dimension (3 x nchannels) in filter weights matrix matW
                // iPos * 3 is row index where to write the virtual sensor (containing of 3 row vectors) for source position iPos into matW
            }else{
                qWarning("[MNEBeamformerWeights::make_beamformer_weights] Vector version of NAI weight normalization is not implemented TODO");
                //TODO. maybe quit calculations here
            }

        }else if(bNormUnitNoise){
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

            qDebug() << "[MNEBeamformerweights::make_beamformer_weights] bNormArray = true";

            // filt*lf = ||lf||, applies to scalar leadfield, and to one of the possibilities of the vector version, eqn. 4.75
            //TODO: Frobenius norm is better s. Wester 2022 so norm() is used which returns Frobenius norm of matrix
            //normalize local leadfield to Frobenius Norm (s. Westner 2022, why Frobenius is better than L2 as used in Fieldtrip)
            MatrixXd matLocLeadfieldNorm = MatrixXd::Zero(matLocLeadfield.rows(),matLocLeadfield.cols());
            RowVectorXd vecFrobNorm = compute_frobenius_norm(matLocLeadfield);

            for(int i = 0; i < vecFrobNorm.size(); i++){
                matLocLeadfieldNorm.col(i) = matLocLeadfield.col(i) / vecFrobNorm(i);
            }

            std::cout.precision(18);
            std::cout << "[MNEBeamformerWeights::compute_forbenius_norm] vecFrobNorm = " << vecFrobNorm << '\n';

            //L2 norm for normalization (rotation variant, so Frobenius norm is prefered)
            //MatrixXd matLocLeadfieldNorm = matLocLeadfield / matLocLeadfield.norm();
            //MatrixXd matLocLeadfieldNorm = matLocLeadfield / fFrobNorm;

            matLocVirtSens = compute_pseudo_inverse(matLocLeadfieldNorm.transpose() * matDataCovInv * matLocLeadfieldNorm) * matLocLeadfieldNorm.transpose() * matDataCovInv; // S&N eqn. 4.09 (scalar version), and eqn. 4.75 (vector version)
            matWT.block(iPos*3,0,3,matLocVirtSens.cols()) = matLocVirtSens; //store local virtual sensor with dimension (3 x nchannels) in filter weights matrix matWT

            //test array gain beamformer condition
            MatrixXd matTestBFCond = matLocVirtSens * matLocLeadfieldNorm;
            std::cout.precision(18);
            std::cout << "[MNEBeamformerWeights::make_beamformer_weights] matTestBFCond arraygain: " << matTestBFCond << '\n';

        }else{ //unitgain BF without weight normalization



            // this is the 'standard' unit gain constraint spatial filter: filt*lf=I, applies both to vector and scalar leadfields
            //HINT: compute_pseudo_inverse should work correctly
            matLocVirtSens = compute_pseudo_inverse(matLocLeadfield.transpose() * matDataCovInv * matLocLeadfield) * matLocLeadfield.transpose() * matDataCovInv;
            matWT.block(iPos*3,0,3,matLocVirtSens.cols()) = matLocVirtSens; //store local virtual sensor with dimension (3 x nchannels) in filter weights matrix matWT

            MatrixXd matTestBFCond = matLocVirtSens * matLocLeadfield;
            std::cout.precision(18);
            std::cout << "[MNEBeamformerWeights::make_beamformer_weights] matTestBFCond no weight normalization: " << matTestBFCond << '\n';

        }


        qInfo("[MNEBeamformerWeights::make_beamformer_weights] Constructed local virtual sensor.\n");


        qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matLocVirtSens dim: " << matLocVirtSens.rows() << " x " << matLocVirtSens.cols();
        qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] dataCovPicked.data dim: " << dataCovPicked.data.rows() << " x " << dataCovPicked.data.cols();
        qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matLocVirtSens.transpose() dim: " << matLocVirtSens.transpose().rows() << " x " << matLocVirtSens.transpose().cols();
        qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matWT dim: " << matWT.rows() << " x " << matWT.cols();


        //calculate source cov matrix (necessary for projecting the dipole moment timecourse on the direction of maximal power)
        //matrix of local source covariance
        MatrixXd matSourceCov = matLocVirtSens * dataCovPicked.data * matLocVirtSens.transpose();


        //if projectmom: project the dipole moment timecourse on the direction of maximal power (store it)
        //TODO: understand when and why we need this
        if(p_bProjectMom){
            JacobiSVD<MatrixXd> svd(matSourceCov);
            MatrixXd matU = svd.matrixU();
            VectorXd vecMom = matU.block(0,0,matU.rows(),1); //first right singular vector is dominant dipole direction
            matLocVirtSens = vecMom.transpose() * matLocVirtSens; //TODO: find out why we modify  virtual sensor here
        }

        //HINT: in contrast to fieldtrip no option for power estimation, we use the trace method
        vecSourcePow[iPos] = matSourceCov.trace();


        //if project noise: estimate noise power projected through the filter
        //HINT: in contrast to fieldtrip we use the trace option for source power estimation
        if(p_bEstNoisePow){
            vecNoisePow[iPos] = iNoiseLevelDataCov * (matLocVirtSens * matLocVirtSens.transpose()).trace();
        }

    }//for iPos (scanning the source grid)

    qInfo("[MNEBeamformerWeights::make_beamformer_weights] Finished scanning the grid of source positions and calculating virtual sensors for each position.");

    //std::cout << "[MNEBeamformerWeights::make_beamformer_weights] matWT: " << matWT;
    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matWT max: " << matWT.maxCoeff();
    qDebug() << "[MNEBeamformerWeights::make_beamformer_weights] matWT min: " << matWT.minCoeff();


    //store filter weights and additional info
    //HINT: partly copied from make_inverse_operator
    p_MNEBeamformerWeights.info = dataInfo; //TODO: this info is wrong here since prepare forward solution modifies number of used channels -> reduce this info to this channels too?
    p_MNEBeamformerWeights.weights = matWT;
    p_MNEBeamformerWeights.data_cov = FiffCov::SDPtr(new FiffCov(dataCovPicked)); //TODO check whether this creation of FiffCov matrix works correct here
    p_MNEBeamformerWeights.noise_cov = FiffCov::SDPtr(new FiffCov(outNoiseCov));
    p_MNEBeamformerWeights.weightNorm = p_sWeightnorm;
    p_MNEBeamformerWeights.whitener = matWhitener;
    p_MNEBeamformerWeights.fixedOri = p_bFixedOri;
    p_MNEBeamformerWeights.optOri = matOptimalOri;
    p_MNEBeamformerWeights.nsource = forward.nsource;
    p_MNEBeamformerWeights.nchan = forward.nchan;
    p_MNEBeamformerWeights.sourcePowEst = vecSourcePow;
    p_MNEBeamformerWeights.noisePowEst = vecNoisePow;
    p_MNEBeamformerWeights.projs = dataInfo.projs;
    p_MNEBeamformerWeights.src = forward.src;
    p_MNEBeamformerWeights.nave = p_iNAverage;


    qInfo("[MNEBeamformerWeights::make_beamformer_weights] Finished calculation of beamformer weights.");

    return p_MNEBeamformerWeights;

}

//=============================================================================================================

MNEBeamformerWeights MNEBeamformerWeights::prepare_beamformer_weights() const
{

    // TODO: this method is implemented because the inverse operator routine includes one similar method,
    // this inverse operator prepare method constructs the inverse operator by making one, (we need this for beamformer weights)
    // does some scaling stuff dependent on the number of averages (might be necessary for lcmv too, but not sure -> first draft only for raw data),
    // creates a regularized inverter (i think we dont need that here),
    // an ssp projector (ssp for noise reduction might be helpful here (mnepy does include proj too),
    // creates a whitener (we dont need that because whitener has to be constructed prior W computation in make_beamformer_weights)
    // and creates noise normalization factors (not necessary for lcmv)



    //construct MNEBeamformerWeights
    //HINT: analog to prepare_inverse_operator
    qInfo("[MNEBeamformerWeights::prepare_beamformer_weights] Preparing the beamformer weights for use...\n");
    MNEBeamformerWeights p_MNEBeamformerWeights(*this);





    //TODO for debugging
    qDebug() << "[MNEBeamformerWeights::prepare_beamformer_weights] p_MNEBeamformerWeights.weights dim: " << p_MNEBeamformerWeights.weights.rows() << " x " << p_MNEBeamformerWeights.weights.cols();


    //
    //   Create the projection operator
    // HINT: copied from prepare_inverse_operator, we dont need it here probably
    //TODO: this projector is not needed
    //
/*    qint32 ncomp = FiffProj::make_projector(p_MNEBeamformerWeights.projs,p_MNEBeamformerWeights.noise_cov->names, p_MNEBeamformerWeights.proj);
    if (ncomp > 0)
        printf("\t MNEBeamformerWeights::prepare_beamformer_weights Created an SSP operator (subspace dimension = %d)\n",ncomp);
*/

    return p_MNEBeamformerWeights;

}

//=============================================================================================================


MNEBeamformerWeights MNEBeamformerWeights::prepare_beamformer_weights(const FiffInfo &p_dataInfo, const MNEForwardSolution &p_forward, const FiffCov &p_dataCov, const FiffCov &p_noiseCov, QString p_sWeightnorm) const
{

    // TODO: this method is implemented because the inverse operator routine includes one similar method,
    // this inverse operator prepare method constructs the inverse operator by making one, (we need this for beamformer weights)
    // does some scaling stuff dependent on the number of averages (might be necessary for lcmv too, but not sure -> first draft only for raw data),
    // creates a regularized inverter (i think we dont need that here),
    // an ssp projector (ssp for noise reduction might be helpful here (mnepy does include proj too),
    // creates a whitener (we dont need that because whitener has to be constructed prior W computation in make_beamformer_weights)
    // and creates noise normalization factors (not necessary for lcmv)



    //construct MNEBeamformerWeights
    //HINT: analog to prepare_inverse_operator
    qInfo("[MNEBeamformerWeights::prepare_beamformer_weights] Preparing the beamformer weights for use...\n");

    //This constructor creates MNEBeamformerWeights object
    MNEBeamformerWeights p_MNEBeamformerWeights(p_dataInfo,
                                                p_forward,
                                                p_dataCov,
                                                p_noiseCov,
                                                p_sWeightnorm);



    //TODO for debugging
    qDebug() << "[MNEBeamformerWeights::prepare_beamformer_weights] p_MNEBeamformerWeights.weights dim: " << p_MNEBeamformerWeights.weights.rows() << " x " << p_MNEBeamformerWeights.weights.cols();


    //
    //   Create the projection operator
    // HINT: copied from prepare_inverse_operator, we dont need it here probably
    //TODO: this projector is not needed
    //
/*    qint32 ncomp = FiffProj::make_projector(p_MNEBeamformerWeights.projs,p_MNEBeamformerWeights.noise_cov->names, p_MNEBeamformerWeights.proj);
    if (ncomp > 0)
        printf("\t MNEBeamformerWeights::prepare_beamformer_weights Created an SSP operator (subspace dimension = %d)\n",ncomp);
*/

    return p_MNEBeamformerWeights;

}

//=============================================================================================================




