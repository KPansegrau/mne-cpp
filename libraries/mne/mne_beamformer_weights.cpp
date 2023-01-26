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

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QFuture>
#include <QtConcurrent>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace MNELIB;
using namespace FIFFLIB;
using namespace Eigen;

//=============================================================================================================
// DEFINE GLOBAL METHODS
//=============================================================================================================

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

MNEBeamformerWeights::MNEBeamformerWeights()
{
}

//=============================================================================================================

MNEBeamformerWeights::MNEBeamformerWeights(const MNEBeamformerWeights &p_MNEBeamformerWeights)
{
}

//=============================================================================================================

MNEBeamformerWeights::~MNEBeamformerWeights()
{
}

//=============================================================================================================

QStringList MNEBeamformerWeights::compare_ch_names(const QStringList &chNamesA,
                                                   const QStringList &chNamesB) const
{

    // HINT: adapted from inverse operator check_ch_names
    // HINt: similar to compare_channel_names form mnepy but without handling of bads
    QStringList ch_names; //list of common channels

    for(qint32 i = 0; i < chNamesA.size(); ++i){
        if(chNamesB.contains(chNamesA[i])){
            ch_names.append(chNamesA[i]);
        }
    }

    return ch_names;
}


//=============================================================================================================

QStringList MNEBeamformerWeights::get_common_channels(const FiffInfo &info,
                                                 const MNEForwardSolution &forward,
                                                 const FiffCov &p_data_cov,
                                                 const FiffCov &p_noise_cov) const

{
    //TODO: compare channels and return QStringList of common channels
    //HINT: Steps analog to _check_info_inv of mnepy

    QStringList qCommonChNames;

    QStringList qDataChNames = info.ch_names;
    QStringList qFwdChNames = forward.info.ch_names;
    QStringList qDataCovChNames = p_data_cov.names;
    QStringList qNoiseCovChNames = p_noise_cov.names;

    //compare data to fwd (return all channels in data that are not bad and are present in fwd in ch_names)
    qCommonChNames = compare_ch_names(qDataChNames,qFwdChNames);

    //TODO start here next time

    //see notes 26.01.2023 use pick_info from fiff_info_base


    //make sure no reference channels are left
    //handle if data cov and noise cov have no bads -> set bads to data bad channels
    //if data cov not none: compare channels common to fwd that are in data cov and not are not bad, update ch_names
    //if noise cov not none: same procedure as for data cov

    return qCommonChNames;
}


//=============================================================================================================

MNEBeamformerWeights MNEBeamformerWeights::make_beamformer_weights(const FiffInfo &info,
                                                                   MNEForwardSolution forward,
                                                                   const FiffCov &p_data_cov,
                                                                   float regularizationFactor,
                                                                   const FiffCov &p_noise_cov,
                                                                   //label
                                                                   //pick_ori
                                                                   //rank
                                                                   //weight_norm
                                                                   bool reduce_rank,
                                                                   //depth
                                                                   //inversion
                                                                   float depth)
{
    //HINT: this method should resemble make_inverse_operator (mne cpp) and make_lcmv (mne py)

    //bool is_fixed_ori = forward.isFixedOrient();
    MNEBeamformerWeights p_MNEBeamformerWeights; //HINT: named similar to p_MNEInverseOperator, dont understand why p_ yet

    //Step 1: ensure that measurement data channels match forward model, noise covariance matrix and data covariance matrix channels
    // return the channels that are common to all four objects
    //compare p_data_cov.names p_noise_cov.names and info.ch_names




    //Step 2: restrict info to the common channels in forward model, noise and data covariance matrix

    //Step 3: compute rank of data and noise covariance matrix, error if they do not match, set rank to data cov rank

    //Step 4: check depth (in mnepy) but not needed here because depth 0 (=NONE)

    //Step 5: if inversion single (dont know if i need this, maybe only implement one method for inversion in case there is no inversion method yet

    //Step 6: call prepare_beamformer_input, get channels present in prepared variables

    //Step 7: use data cov channels that are present in prepared variables, get data cov matrix and integrate rank values,

    //Step 8: call _compute_beamformer which computes the weights W and max_power_orientation

    //Step 9: store source type and subject

    //Step 10: is computed BF scalar or vector

    //Step 11: create Beamformer object with calculated W and parameters (maybe not here but in beamformer class (inverse lib))

    //Step X: check which sensor types are present
    //HINT: adapted from make_inverse_operator (l. 862)
    //TODO: changes: differentiate between mag and grad (instead of only meg) because whitening is also necessary if two different meg sensor types are present



    //set members



    //HINT: make_inverse_operator
/*    bool is_fixed_ori = forward.isFixedOrient();
    MNEInverseOperator p_MNEInverseOperator;

    std::cout << "ToDo MNEInverseOperator::make_inverse_operator: do surf_ori check" << std::endl;

    //Check parameters
    if(fixed && loose > 0)
    {
        qWarning("Warning: When invoking make_inverse_operator with fixed = true, the loose parameter is ignored.\n");
        loose = 0.0f;
    }

    if(is_fixed_ori && !fixed)
    {
        qWarning("Warning: Setting fixed parameter = true. Because the given forward operator has fixed orientation and can only be used to make a fixed-orientation inverse operator.\n");
        fixed = true;
    }

    if(forward.source_ori == -1 && loose > 0)
    {
        qCritical("Error: Forward solution is not oriented in surface coordinates. loose parameter should be 0 not %f.\n", loose);
        return p_MNEInverseOperator;
    }

    if(loose < 0 || loose > 1)
    {
        qWarning("Warning: Loose value should be in interval [0,1] not %f.\n", loose);
        loose = loose > 1 ? 1 : 0;
        printf("Setting loose to %f.\n", loose);
    }

    if(depth < 0 || depth > 1)
    {
        qWarning("Warning: Depth value should be in interval [0,1] not %f.\n", depth);
        depth = depth > 1 ? 1 : 0;
        printf("Setting depth to %f.\n", depth);
    }

    //
    // 1. Read the bad channels
    // 2. Read the necessary data from the forward solution matrix file
    // 3. Load the projection data
    // 4. Load the sensor noise covariance matrix and attach it to the forward
    //
    FiffInfo gain_info;
    MatrixXd gain;
    MatrixXd whitener;
    qint32 n_nzero;
    FiffCov p_outNoiseCov;
    forward.prepare_forward(info, p_noise_cov, false, gain_info, gain, p_outNoiseCov, whitener, n_nzero);

    //
    // 5. Compose the depth weight matrix
    //
    FiffCov::SDPtr p_depth_prior;
    MatrixXd patch_areas;
    if(depth > 0)
    {
        std::cout << "ToDo: patch_areas" << std::endl;
//        patch_areas = forward.get('patch_areas', None)
        p_depth_prior = FiffCov::SDPtr(new FiffCov(MNEForwardSolution::compute_depth_prior(gain, gain_info, is_fixed_ori, depth, 10.0, patch_areas, limit_depth_chs)));
    }
    else
    {
        p_depth_prior->data = MatrixXd::Ones(gain.cols(), gain.cols());
        p_depth_prior->kind = FIFFV_MNE_DEPTH_PRIOR_COV;
        p_depth_prior->diag = true;
        p_depth_prior->dim = gain.cols();
        p_depth_prior->nfree = 1;
    }

    // Deal with fixed orientation forward / inverse
    if(fixed)
    {
        if(depth < 0 || depth > 1)
        {
            // Convert the depth prior into a fixed-orientation one
            printf("\tToDo: Picked elements from a free-orientation depth-weighting prior into the fixed-orientation one.\n");
        }
        if(!is_fixed_ori)
        {
            // Convert to the fixed orientation forward solution now
            qint32 count = 0;
            for(qint32 i = 2; i < p_depth_prior->data.rows(); i+=3)
            {
                p_depth_prior->data.row(count) = p_depth_prior->data.row(i);
                ++count;
            }
            p_depth_prior->data.conservativeResize(count, 1);

//            forward = deepcopy(forward)
            forward.to_fixed_ori();
            is_fixed_ori = forward.isFixedOrient();
            forward.prepare_forward(info, p_outNoiseCov, false, gain_info, gain, p_outNoiseCov, whitener, n_nzero);
        }
    }
    printf("\tComputing inverse operator with %d channels.\n", gain_info.ch_names.size());

    //
    // 6. Compose the source covariance matrix
    //
    printf("\tCreating the source covariance matrix\n");
    FiffCov::SDPtr p_source_cov = p_depth_prior;

    // apply loose orientations
    FiffCov::SDPtr p_orient_prior;
    if(!is_fixed_ori)
    {
        p_orient_prior = FiffCov::SDPtr(new FiffCov(forward.compute_orient_prior(loose)));
        p_source_cov->data.array() *= p_orient_prior->data.array();
    }

    // 7. Apply fMRI weighting (not done)

    //
    // 8. Apply the linear projection to the forward solution
    // 9. Apply whitening to the forward computation matrix
    //
    printf("\tWhitening the forward solution.\n");
    gain = whitener*gain;

    // 10. Exclude the source space points within the labels (not done)

    //
    // 11. Do appropriate source weighting to the forward computation matrix
    //

    // Adjusting Source Covariance matrix to make trace of G*R*G' equal
    // to number of sensors.
    printf("\tAdjusting source covariance matrix.\n");
    RowVectorXd source_std = p_source_cov->data.array().sqrt().transpose();

    for(qint32 i = 0; i < gain.rows(); ++i)
        gain.row(i) = gain.row(i).array() * source_std.array();

    double trace_GRGT = (gain * gain.transpose()).trace();//pow(gain.norm(), 2);
    double scaling_source_cov = (double)n_nzero / trace_GRGT;

    p_source_cov->data.array() *= scaling_source_cov;

    gain.array() *= sqrt(scaling_source_cov);

    // now np.trace(np.dot(gain, gain.T)) == n_nzero
    // logger.info(np.trace(np.dot(gain, gain.T)), n_nzero)

    //
    // 12. Decompose the combined matrix
    //
    printf("Computing SVD of whitened and weighted lead field matrix.\n");
    JacobiSVD<MatrixXd> svd(gain, ComputeThinU | ComputeThinV);
    std::cout << "ToDo Sorting Necessary?" << std::endl;
    VectorXd p_sing = svd.singularValues();
    MatrixXd t_U = svd.matrixU();
    MNEMath::sort<double>(p_sing, t_U);
    FiffNamedMatrix::SDPtr p_eigen_fields = FiffNamedMatrix::SDPtr(new FiffNamedMatrix( svd.matrixU().cols(),
                                                                                        svd.matrixU().rows(),
                                                                                        defaultQStringList,
                                                                                        gain_info.ch_names,
                                                                                        t_U.transpose() ));

    p_sing = svd.singularValues();
    MatrixXd t_V = svd.matrixV();
    MNEMath::sort<double>(p_sing, t_V);
    FiffNamedMatrix::SDPtr p_eigen_leads = FiffNamedMatrix::SDPtr(new FiffNamedMatrix( svd.matrixV().rows(),
                                                                                       svd.matrixV().cols(),
                                                                                       defaultQStringList,
                                                                                       defaultQStringList,
                                                                                       t_V ));
    printf("\tlargest singular value = %f\n", p_sing.maxCoeff());
    printf("\tscaling factor to adjust the trace = %f\n", trace_GRGT);

    qint32 p_nave = 1.0;

    // Handle methods
    bool has_meg = false;
    bool has_eeg = false;

    RowVectorXd ch_idx(info.chs.size());
    qint32 count = 0;
    for(qint32 i = 0; i < info.chs.size(); ++i)
    {
        if(gain_info.ch_names.contains(info.chs[i].ch_name))
        {
            ch_idx[count] = i;
            ++count;
        }
    }
    ch_idx.conservativeResize(count);

    for(qint32 i = 0; i < ch_idx.size(); ++i)
    {
        QString ch_type = info.channel_type(ch_idx[i]);
        if (ch_type == "eeg")
            has_eeg = true;
        if ((ch_type == "mag") || (ch_type == "grad"))
            has_meg = true;
    }

    qint32 p_iMethods;

    if(has_eeg && has_meg)
        p_iMethods = FIFFV_MNE_MEG_EEG;
    else if(has_meg)
        p_iMethods = FIFFV_MNE_MEG;
    else
        p_iMethods = FIFFV_MNE_EEG;

    // We set this for consistency with mne C code written inverses
    if(depth == 0)
        p_depth_prior = FiffCov::SDPtr();

    p_MNEInverseOperator.eigen_fields = p_eigen_fields;
    p_MNEInverseOperator.eigen_leads = p_eigen_leads;
    p_MNEInverseOperator.sing = p_sing;
    p_MNEInverseOperator.nave = p_nave;
    p_MNEInverseOperator.depth_prior = p_depth_prior;
    p_MNEInverseOperator.source_cov = p_source_cov;
    p_MNEInverseOperator.noise_cov = FiffCov::SDPtr(new FiffCov(p_outNoiseCov));
    p_MNEInverseOperator.orient_prior = p_orient_prior;
    p_MNEInverseOperator.projs = info.projs;
    p_MNEInverseOperator.eigen_leads_weighted = false;
    p_MNEInverseOperator.source_ori = forward.source_ori;
    p_MNEInverseOperator.mri_head_t = forward.mri_head_t;
    p_MNEInverseOperator.methods = p_iMethods;
    p_MNEInverseOperator.nsource = forward.nsource;
    p_MNEInverseOperator.coord_frame = forward.coord_frame;
    p_MNEInverseOperator.source_nn = forward.source_nn;
    p_MNEInverseOperator.src = forward.src;
    p_MNEInverseOperator.info = forward.info;
    p_MNEInverseOperator.info.bads = info.bads;

    return p_MNEInverseOperator;
*/


}

//=============================================================================================================

MNEBeamformerWeights MNEBeamformerWeights::prepare_beamformer_weights(qint32 nave,
                                                                      float lambda2) const
{
    //HINT: this method should resemple prepare_inverse_operator (mne cpp) and prepare_beamformer_input (mne py)
    //HINT: prepare_inverse_operator is called in doInverseSetup of MinimumNorm class





    //HINT: prepare_inverse_operator
/*    if(nave <= 0)
    {
        printf("The number of averages should be positive\n");
        return MNEInverseOperator();
    }
    printf("Preparing the inverse operator for use...\n");
    MNEInverseOperator inv(*this);
    //
    //   Scale some of the stuff
    //
    float scale     = ((float)inv.nave)/((float)nave);
    inv.noise_cov->data  *= scale;
    inv.noise_cov->eig   *= scale;
    inv.source_cov->data *= scale;
    //
    if (inv.eigen_leads_weighted)
        inv.eigen_leads->data *= sqrt(scale);
    //
    printf("\tScaled noise and source covariance from nave = %d to nave = %d\n",inv.nave,nave);
    inv.nave = nave;
    //
    //   Create the diagonal matrix for computing the regularized inverse
    //
    VectorXd tmp = inv.sing.cwiseProduct(inv.sing) + VectorXd::Constant(inv.sing.size(), lambda2);
//    if(inv.reginv)
//        delete inv.reginv;
    inv.reginv = VectorXd(inv.sing.cwiseQuotient(tmp));
    printf("\tCreated the regularized inverter\n");
    //
    //   Create the projection operator
    //

    qint32 ncomp = FiffProj::make_projector(inv.projs, inv.noise_cov->names, inv.proj);
    if (ncomp > 0)
        printf("\tCreated an SSP operator (subspace dimension = %d)\n",ncomp);

    //
    //   Create the whitener
    //
//    if(inv.whitener)
//        delete inv.whitener;
    inv.whitener = MatrixXd::Zero(inv.noise_cov->dim, inv.noise_cov->dim);

    qint32 nnzero, k;
    if (inv.noise_cov->diag == 0)
    {
        //
        //   Omit the zeroes due to projection
        //
        nnzero = 0;

        for (k = ncomp; k < inv.noise_cov->dim; ++k)
        {
            if (inv.noise_cov->eig[k] > 0)
            {
                inv.whitener(k,k) = 1.0/sqrt(inv.noise_cov->eig[k]);
                ++nnzero;
            }
        }
        //
        //   Rows of eigvec are the eigenvectors
        //
        inv.whitener *= inv.noise_cov->eigvec;
        printf("\tCreated the whitener using a full noise covariance matrix (%d small eigenvalues omitted)\n", inv.noise_cov->dim - nnzero);
    }
    else
    {
        //
        //   No need to omit the zeroes due to projection
        //
        for (k = 0; k < inv.noise_cov->dim; ++k)
            inv.whitener(k,k) = 1.0/sqrt(inv.noise_cov->data(k,0));

        printf("\tCreated the whitener using a diagonal noise covariance matrix (%d small eigenvalues discarded)\n",ncomp);
    }
    //
    //   Finally, compute the noise-normalization factors
    //
    if (dSPM || sLORETA)
    {
        VectorXd noise_norm = VectorXd::Zero(inv.eigen_leads->nrow);
        VectorXd noise_weight;
        if (dSPM)
        {
           printf("\tComputing noise-normalization factors (dSPM)...");
           noise_weight = VectorXd(inv.reginv);
        }
        else
        {
           printf("\tComputing noise-normalization factors (sLORETA)...");
           VectorXd tmp = (VectorXd::Constant(inv.sing.size(), 1) + inv.sing.cwiseProduct(inv.sing)/lambda2);
           noise_weight = inv.reginv.cwiseProduct(tmp.cwiseSqrt());
        }
        VectorXd one;
        if (inv.eigen_leads_weighted)
        {
           for (k = 0; k < inv.eigen_leads->nrow; ++k)
           {
              one = inv.eigen_leads->data.block(k,0,1,inv.eigen_leads->data.cols()).cwiseProduct(noise_weight);
              noise_norm[k] = sqrt(one.dot(one));
           }
        }
        else
        {
//            qDebug() << 32;
            double c;
            for (k = 0; k < inv.eigen_leads->nrow; ++k)
            {
//                qDebug() << 321;
                c = sqrt(inv.source_cov->data(k,0));
//                qDebug() << 322;
//                qDebug() << "inv.eigen_leads->data" << inv.eigen_leads->data.rows() << "x" << inv.eigen_leads->data.cols();
//                qDebug() << "noise_weight" << noise_weight.rows() << "x" << noise_weight.cols();
                one = c*(inv.eigen_leads->data.row(k).transpose()).cwiseProduct(noise_weight);//ToDo eigenleads data -> pointer
                noise_norm[k] = sqrt(one.dot(one));
//                qDebug() << 324;
            }
        }

//        qDebug() << 4;

        //
        //   Compute the final result
        //
        VectorXd noise_norm_new;
        if (inv.source_ori == FIFFV_MNE_FREE_ORI)
        {
            //
            //   The three-component case is a little bit more involved
            //   The variances at three consequtive entries must be squeared and
            //   added together
            //
            //   Even in this case return only one noise-normalization factor
            //   per source location
            //
            VectorXd* t = MNEMath::combine_xyz(noise_norm.transpose());
            noise_norm_new = t->cwiseSqrt();//double otherwise values are getting too small
            delete t;
            //
            //   This would replicate the same value on three consequtive
            //   entries
            //
            //   noise_norm = kron(sqrt(mne_combine_xyz(noise_norm)),ones(3,1));
        }
        VectorXd vOnes = VectorXd::Ones(noise_norm_new.size());
        VectorXd tmp = vOnes.cwiseQuotient(noise_norm_new.cwiseAbs());
//        if(inv.noisenorm)
//            delete inv.noisenorm;

        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        tripletList.reserve(noise_norm_new.size());
        for(qint32 i = 0; i < noise_norm_new.size(); ++i)
            tripletList.push_back(T(i, i, tmp[i]));

        inv.noisenorm = SparseMatrix<double>(noise_norm_new.size(),noise_norm_new.size());
        inv.noisenorm.setFromTriplets(tripletList.begin(), tripletList.end());

        printf("[done]\n");
    }
    else
    {
//        if(inv.noisenorm)
//            delete inv.noisenorm;
        inv.noisenorm = SparseMatrix<double>();
    }

    return inv;
*/
}





