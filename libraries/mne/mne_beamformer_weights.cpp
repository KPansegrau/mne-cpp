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
{
    qRegisterMetaType<QSharedPointer<MNELIB::MNEBeamformerWeights> >("QSharedPointer<MNELIB::MNEBeamformerWeights>");
    qRegisterMetaType<MNELIB::MNEBeamformerWeights>("MNELIB::MNEBeamformerWeights");
}



//=============================================================================================================


MNEBeamformerWeights::MNEBeamformerWeights(FiffInfo &p_dataInfo,
                                           MNEForwardSolution &p_forward,
                                           FiffCov &p_dataCov,
                                           const FiffCov &p_noiseCov,
                                           QString p_sPowMethod,
                                           bool p_bFixedOri,
                                           bool p_bEstNoisePow,
                                           bool p_bProjectMom,
                                           QString p_sWeightnorm,
                                           qint32 p_iRegParam)
{
    *this = MNEBeamformerWeights::make_beamformer_weights(p_dataInfo,p_forward,p_dataCov,p_noiseCov,p_sPowMethod,p_bFixedOri,p_bEstNoisePow,p_bProjectMom,p_sWeightnorm,p_iRegParam);
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

{
    qRegisterMetaType<QSharedPointer<MNELIB::MNEBeamformerWeights> >("QSharedPointer<MNELIB::MNEBeamformerWeights>");
    qRegisterMetaType<MNELIB::MNEBeamformerWeights>("MNELIB::MNEBeamformerWeights");

}



//=============================================================================================================

MNEBeamformerWeights::~MNEBeamformerWeights()

{
}



//=============================================================================================================

MatrixXd MNEBeamformerWeights::compute_pseudo_inverse(const MatrixXd &p_matrix, double p_dEpsilon) const
{
    //TODO maybe move this to mne math in the end
    //HINT: similar to function pinv in ft_inverse_lcmv (same as Matlab but default tolerance is twice as high, but without handling of tolerance value as input
    //HINT: used this version for calculation with Eigen from GitHub https://gist.github.com/pshriwise/67c2ae78e5db3831da38390a8b2a209f, adapted it a bit



    //dimensions of input matrix
    qint32 nrows = p_matrix.rows();
    qint32 ncols = p_matrix.cols();

    JacobiSVD<MatrixXd> svd;
    if(nrows == ncols){//squared matrix
        JacobiSVD<MatrixXd> svd(p_matrix ,ComputeFullU | ComputeFullV);
    }else{
        // For a non-square matrix
        JacobiSVD<MatrixXd> svd(p_matrix ,ComputeThinU | ComputeThinV);
    }

    //calculate tolerance
    double tolerance = p_dEpsilon * std::max(ncols, nrows) *svd.singularValues().array().abs()(0);

    // return inverse
    MatrixXd matInv = svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
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

    //decompose matrix with svd
    JacobiSVD<MatrixXd> svd(p_dataCov.data, ComputeFullU | ComputeFullV); //TODO: check whether we need full U and V here (full/thin defines dimension of U and V)
    VectorXd vecSing = svd.singularValues(); //sorted in decreasing order
    MatrixXd matU = svd.matrixU();
    MatrixXd matV = svd.matrixV();

    //invert only non-zero singular values (invert only singular values with index
    qint32 rank = MNEMath::rank(p_dataCov.data); //dont know whether this rank calculation is too easy here (no differentiation between channel types)
    VectorXd vecSingInv = VectorXd::Zero(vecSing.size());
    for(int i = 0; i < (p_dataCov.dim - rank); i++){
        vecSingInv[i] = 1 / vecSing[i];
    }

    //compute the pseudo inverse by recompose svd results
    MatrixXd matCovInv = matV * vecSingInv *  matU;
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

MNEBeamformerWeights MNEBeamformerWeights::make_beamformer_weights(//const MatrixXd &p_matData, //the measurement data matrix (input in fieldtrip, we dont need this because no subspace projection procedure in this first implementation)
                                                                   FiffInfo &p_dataInfo, //the info of the measurement data (no input in fieldtrip)
                                                                   MNEForwardSolution &p_forward, //the forward solution
                                                                   FiffCov &p_dataCov, //the data covariance matrix
                                                                   const FiffCov &p_noiseCov, //the noise covariance matrix (no input in fieldtrip but we need it for whitening)
                                                                   QString p_sPowMethod, //can be 'trace' (default s. Van Veen 1997) or 'lambda1'
                                                                   bool p_bFixedOri, // true for fixed orientation, false for free orientation (default=false)
                                                                   bool p_bEstNoisePow, // estimate noise power projected through filter (default=true)
                                                                   bool p_bProjectMom, // true: project the dipole moment time course on the direction of maximal power (default=false), TODO: check what this is for, when do we need it and why?
                                                                   //bool p_bKurtosis, //compute kurtosis of dipole timeseries (s. goodnotes zeigt ob time series spiked an betrachteter Position, s. Hall 2018)
                                                                   QString p_sWeightnorm, //normalize the beamformer weights ('no' (=unitgain,default), 'unitnoisegain', 'arraygain' or 'nai')
                                                                   qint32 p_iRegParam //the regularization parameter //called lambda in fieldtrip
                                                                   //qint32 &p_iKappa,
                                                                   //qint32 &p_iTol
                                                                   ) {
    //HINT: resemling the FieldTrip code in ft_inverse_lcmv.m
    //because mnepy code is highly modular and too complex for first lcmv implementation
    //decided on FieldTrip implementation because it is easier to understand (maybe less overload code)
    // TODO: add to outlook: for better symmetry to exsisting code in mnecpp, restructure code according to mnepy implementation

    //TODO?:split this method into prepare_beamformer_input and make beamformer weights

    //TODO: check whether mnepy includes some useful steps we need here

    //TODO: option for rank reduction (s. Westner 2022) since for MEG it is recommende to reduce the rank of the forward field

    //TODO: check all equations from Sekihara Nagarajan 2008 (maybe in Westner 2022 and Knösche 2022 so no need of book)

    //TODO: include number of averages for the evoked input data

    qInfo("MNEBeamformerWeights::make_beamformer_weights Start calculation of beamformer weights...\n");


    //prepare output
    MNEBeamformerWeights p_MNEBeamformerWeights;

    //use this logical flags instead of string comparisons (idea from fieldtrip ft_inverse_lcmv, saves time during scanning loop)
    bool bPowTrace = (p_sPowMethod == "trace");
    bool bPowLambda1 = (p_sPowMethod == "lambda1");
    bool bNormNo = (p_sWeightnorm == "no");
    bool bNormUnitNoise = (p_sWeightnorm == "unitnoisegain");
    bool bNormArray = (p_sWeightnorm == "arraygain");
    bool bNormNai = (p_sWeightnorm == "nai");

    //check invalid parameter combinations (partly copied from make_inverse_operator mnecpp)
    bool is_fixed_ori = p_forward.isFixedOrient();
    std::cout << "ToDo MNEBeamformerWeights::make_beamformer_weights: do surf_ori check\n"; //TODO: this warning was copied from make_inverse_operator, check whether we need it here too
    if(is_fixed_ori && !p_bFixedOri)
    {
        qWarning("MNEBeamformerWeights::make_beamformer_weights Warning: Setting p_bFixedOri parameter = true. Because the given forward operator has fixed orientation and can only be used to make a fixed-orientation beamformer weights.\n");
        p_bFixedOri = true;
    }

    //ensure that p_sPowMethod option is validm, if not set to default (trace).
    if(!bPowTrace && !bPowLambda1){
        qWarning("MNEBeamformerWeights::make_beamformer_weights Warning: Invalid option for p_sPowMethod. Setting p_sPowMethod = 'trace' (default).\n");
        p_sPowMethod = "trace";
    }

    //ensure that p_sWeightnorm option is valid, if not set to default (no weight normalization)
    if(!bNormNo && !bNormUnitNoise && !bNormArray && !bNormNai){
        qWarning("MNEBeamformerWeights::make_beamformer_weights Warning: Invalid weightnorm option. Set p_sWeightnorm = 'no' (default) and compute a unit-gain beamformer without weight normalization.\n");
        p_sWeightnorm = "no";
    }


    //TODO: do we need to ensure that there is only one channel type? check whether fieldtrip code is only for eeg or meg but not mixed channels

    qInfo("MNEBeamformerWeights::make_beamformer_weights Prepare measurement info...\n");

    // ensure that measurement data channels match forward model, noise covariance matrix and data covariance matrix channels (from mnepy make inverse operator, first step there)
    //TODO: channel selection stuff is performerd in miminumnorms calculate Inverse, move it?
    // return the channels that are common to all four objects
    //from mnepy
    QStringList lCommonChanNames = p_MNEBeamformerWeights.check_info_bf(p_dataInfo, p_forward, p_dataCov, p_noiseCov);
    //restrict data info to common channels
    RowVectorXi picksCommonChan; //indices of channels that are common to data, forward, and both covariance matrices
    picksCommonChan = FiffInfoBase::pick_channels(lCommonChanNames);
    p_dataInfo.pick_info(picksCommonChan);

    qInfo("MNEBeamformerWeights::make_beamformer_weights Finished preparation of measurement info.\n");

    qInfo("MNEBeamformerWeights::make_beamformer_weights Prepare forward solution...\n");

    //prepare forward operator: we use the prepare_forward method form MNEForwardSolution
    //existing prepare_forward method: selects common channel, prepares noise cov matrix, calculates whitener from noise cov, gain matrix
    //HINT: these lines are copied from make_inverse_operator
    FiffInfo leadfieldInfo;
    MatrixXd matLeadfield;
    FiffCov outNoiseCov;
    MatrixXd matWhitener; //necessary for whitening of data covariance matrix and forward solution
    qint32 n_nzero; //total rank of ??? noise covariance matrix? leadfield? whitener? all of them? (used in make_inverse_operator to scale the source covariance matrix) maybe we dont need that value here
    p_forward.prepare_forward(p_dataInfo, p_noiseCov, false, leadfieldInfo, matLeadfield, outNoiseCov, matWhitener, n_nzero);

    //Whiten the forward solution with whitener
    matLeadfield = matWhitener * matLeadfield;

    qInfo("MNEBeamformerWeights::make_beamformer_weights Finished preparation of forward solution.\n");


    qInfo("MNEBeamformerWeights::make_beamformer_weights Prepare data covariance matrix...\n");

    //estimate noise level in the covariance matrix by the smallest non-zero singular value, always needed for NAI weight normalization
    //TODO: maybe we can do this step below? e.g. below regularization does not work because else statement uses cov data which is modified during regularization
    qint32 iNoiseLevelDataCov; //only for nai and unitnoisegain BF relevant
    if(bNormNai || p_bEstNoisePow){
        // MNEMath::rank method does not include that there are different sensor type with different amplitude ranges
        // but fieldtrip uses matlabs rank function which does not include this problem either (matlab svd returns singular values in descending order)
        if(MNEMath::rank(p_dataCov.data) > p_dataCov.data.size()){ //TODO:check whether this is implemented correct. c++ size() == python len()?, rank caluclated independently from channel types?
            //TODO:  add if no regularization for this warning
            qWarning("MNEBeamformerWeights::make_beamformer_weights Warning: Cannot compute a noise subspace with a full-rank covariance matrix and no regularization.\n");
            //TODO:end if no regularization
            //TODO: else set noise level to regularization parameter
        }else{
            JacobiSVD<MatrixXd> svd(p_dataCov.data);
            VectorXd p_sing = svd.singularValues();
            //TODO: maybe use the eig from FiffCov instance here instead of computing the singular values
            iNoiseLevelDataCov = p_sing[MNEMath::rank(p_dataCov.data)]; //denoted as noise in fieldtrip,
            iNoiseLevelDataCov = std::max(iNoiseLevelDataCov, p_iRegParam);
        }
    } //end if nai or unitnoisegain


    //prepare data covariance matrix for inversion
    //restrict data covariance matrix to common channels
    p_dataCov.pick_channels(lCommonChanNames);

    //TODO:whiten data cov mat (after picking because whitener is from picked forward solution) before regularization
    p_dataCov.data = matWhitener * p_dataCov.data;


    //regularize data covariance matrix by adding values to diagonal entries
    //TODO: regularize takes further parameters that are added according to channel type, use them here
    p_dataCov.regularize(p_dataInfo);

    //invert data covariance matrix
    MatrixXd matDataCovInv = p_MNEBeamformerWeights.invert_data_cov_mat(p_dataCov);
    //check dimension of inverted data cov mat
    if(matDataCovInv.size() != p_dataCov.data.size()){
        qWarning("MNEBeamformerWeights::make_beamformer_weights Warning: Size of inverted data covariance matrix does not match size of data covariance matrix.\n");
    }

    //compute square of inverse data covaricance matrix (needed for unitnoisegain and nai constraint
    MatrixXd matDataCovInvSquared = matDataCovInv * matDataCovInv;

    qInfo("MNEBeamformerWeights::make_beamformer_weights Finished preparation of data covariance matrix.\n");




    //TODO: maybe from check_info_bf up to here in prepare_beamformer_input?


    //start scanning the source grid (LCMV BF uses local leadfield for computation of filter weights)
    qInfo("MNEBeamformerWeights::make_beamformer_weights Start scanning the grid of source positions and calculating virtual sensors for each position...\n");

    //matrix where optimal orientations for fixed orientation beamformers can be stored
    //TODO make this a FiffNamedMatrix to preserve channel names (we need new list of column names then)
    //TODO this matrix has to be stored in MNEBeamformerWeights
    MatrixXd matOptimalOri = MatrixXd::Zero(p_dataInfo.nchan, matLeadfield.cols()/3);
    MatrixXd matW = MatrixXd::Zero(matLeadfield.cols(),matLeadfield.rows()); //matrix of virtual sensors (nsources*3 x nchannels)
    VectorXd vecSourcePow = VectorXd::Zero(matLeadfield.cols()/3,1); //TODO check whether dimension is correct here; vector of activity strength (power) for each source position
    VectorXd vecNoisePow = VectorXd::Zero(matLeadfield.cols()/3,1); //vector of estimated noise power that is projected through the filter for each source position

    for(int iPos = 0; iPos < p_forward.nsource; iPos++){ //iterate through all source positions in source grid

        printf("Scanning source position %i/%i\n", iPos, p_forward.nsource); //count the source positions to show progress

        //extract local leadfield for source location iPos from leadfield matrix
        //nsource = ncol/3 -> 3 columns per dipole in leadfield matrix
        MatrixXd matLocLeadfield = matLeadfield.block(0,iPos,p_forward.nchan,3); //extract a block including all rows and only the dipole at iPos (3 columns from idx iPos on)

        //if fixed orientation: calculate optimal orientation for beamformer with respect to weight normalization type
        //HINT: copied from fieldtrip ft_inverse_lcmv
        //TODO: check all these equations form the book
        if(p_bFixedOri){
            if(bNormUnitNoise || bNormNai){ //HINT: dont use switch case here because it does not work for string cases
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

            qInfo("MNEBeamformerWeights::make_beamformer_weights Calculated optimal orientation.\n");

        }//if p_bFixedOri


        //construct virtual sensor for source location iPos and add it to W with respect to weight normalization type
        //HINT: copied from fieldtrip ft_inverse_lcmv

        MatrixXd matLocVirtSens; //matrix (3 x nchannels) of local virtual sensor

        if(bNormNai){
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
                matW.block(iPos*3,0,3,matLocVirtSens.cols()) = matLocVirtSens; //store local virtual sensor with dimension (3 x nchannels) in filter weights matrix matW
                // iPos * 3 is row index where to write the virtual sensor (containing of 3 row vectors) for source position iPos into matW
            }else{
                qWarning("MNEBeamformerWeights::make_beamformer_weights Vector version of NAI weight normalization is not implemented.\n");
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
                matW.block(iPos*3,0,3,matLocVirtSens.cols()) = matLocVirtSens; //local virtual sensor with dimension (3 x nchannels)
            }else{
                // compute the matrix that is used for scaling of the filter's rows, as per eqn. 4.83
                MatrixXd matDenom = compute_pseudo_inverse(matLocLeadfield.transpose() * matDataCovInv * matLocLeadfield);
                MatrixXd matGamma = matDenom * (matLocLeadfield.transpose() * matDataCovInvSquared * matLocLeadfield) * matDenom;
                // compute the spatial filter, as per eqn. 4.85
                matLocVirtSens = (1/(matGamma.diagonal()).array().sqrt()).matrix().asDiagonal() * matDenom * matLocLeadfield.transpose() * matDataCovInv; //see notes 07.02.2023
                matW.block(iPos*3,0,3,matLocVirtSens.cols()) = matLocVirtSens; //local virtual sensor with dimension (3 x nchannels)
            }

        }else if(bNormArray){
            // filt*lf = ||lf||, applies to scalar leadfield, and to one of the possibilities of the vector version, eqn. 4.75
            //TODO: Frobenius norm is better s. Wester 2022 so norm() is used which returns Frobenius norm of matrix
            //normalize local leadfield to Frobenius Norm (s. Westner 2022, why Frobenius is better than L2 as used in Fieldtrip)
            MatrixXd matLocLeadfieldNorm = matLocLeadfield / matLocLeadfield.norm();
            matLocVirtSens = compute_pseudo_inverse(matLocLeadfieldNorm.transpose() * matDataCovInv * matLocLeadfieldNorm) * matLocLeadfieldNorm.transpose() * matDataCovInv; // S&N eqn. 4.09 (scalar version), and eqn. 4.75 (vector version)

        }else{ //unitgain BF without weight normalization
            // this is the 'standard' unit gain constraint spatial filter: filt*lf=I, applies both to vector and scalar leadfields
            matLocVirtSens = compute_pseudo_inverse(matLocLeadfield.transpose() * matDataCovInv * matLocLeadfield) * matLocLeadfield.transpose() * matDataCovInv;
        }


        qInfo("MNEBeamformerWeights::make_beamformer_weights Constructed local virtual sensor.\n");


        //calculate source cov matrix (TODO: add why we need it)
        //matrix of local source covariance
        MatrixXd matSourceCov = matLocVirtSens * p_dataCov.data * matLocVirtSens.transpose();


        //if projectmom: project the dipole moment timecourse on the direction of maximal power (store it)
        //TODO: understand when and why we need this
        if(p_bProjectMom){
            JacobiSVD<MatrixXd> svd(matSourceCov);
            MatrixXd matU = svd.matrixU();
            VectorXd vecMom = matU.block(0,0,matU.rows(),1); //first right singular vector is dominant dipole direction
            matLocVirtSens = vecMom.transpose() * matLocVirtSens; //TODO: find out why we modify  virtual sensor here
        }

        //if powlambda1 elseif powtrace calculate lambda1 or trace of source covariance matrix
        if(bPowLambda1){
            // determine the largest singular value, which corresponds to the power along the dominant direction
            JacobiSVD<MatrixXd> svd(matSourceCov);
            VectorXd vecSingVals = svd.singularValues();
            vecSourcePow[iPos] = vecSingVals[0];
        }else{ //bPowTrace (default)
            vecSourcePow[iPos] = matSourceCov.trace();
        }


        //if keepmom and we have data: estimate dipole moment at current position (filt*data; TODO: this in apply lcmv method)

        //if compute kurtosis and we have data: caluclate kurtosis of dipole time series (moment computed before; kurt(filt*data), TODO: this in apply lcmv method)

        //if project noise: estimate noise power projected through the filter
        if(p_bEstNoisePow){
            if(bPowLambda1){
                JacobiSVD<MatrixXd> svd(matLocVirtSens * matLocVirtSens.transpose());
                VectorXd vecSingVals = svd.singularValues();
                vecNoisePow[iPos] = iNoiseLevelDataCov * vecSingVals[0];
            }else{//bPowTrace is true (default)
                vecNoisePow[iPos] = iNoiseLevelDataCov * (matLocVirtSens * matLocVirtSens.transpose()).trace();
            }
        }

    }//for iPos (scanning the source grid)

    qInfo("MNEBeamformerWeights::make_beamformer_weights Finished scanning the grid of source positions and calculating virtual sensors for each position.");


    //store filter weights and additional info
    //HINT: partly copied from make_inverse_operator
    p_MNEBeamformerWeights.info = p_dataInfo;
    p_MNEBeamformerWeights.weights = matW;
    p_MNEBeamformerWeights.data_cov = FiffCov::SDPtr(new FiffCov(p_dataCov)); //TODO check whether this creation of FiffCov matrix works correct here
    p_MNEBeamformerWeights.noise_cov = FiffCov::SDPtr(new FiffCov(outNoiseCov));
    p_MNEBeamformerWeights.weightNorm = p_sWeightnorm;
    p_MNEBeamformerWeights.whitener = matWhitener;
    p_MNEBeamformerWeights.fixedOri = p_bFixedOri;
    p_MNEBeamformerWeights.optOri = matOptimalOri;
    p_MNEBeamformerWeights.nsource = p_forward.nsource;
    p_MNEBeamformerWeights.nchan = p_forward.nchan;
    p_MNEBeamformerWeights.sourcePowEst = vecSourcePow;
    p_MNEBeamformerWeights.noisePowEst = vecNoisePow;


    qInfo("MNEBeamformerWeights::make_beamformer_weights Finished calculation of beamformer weights.");

    return p_MNEBeamformerWeights;

}







