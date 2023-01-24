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

MNEBeamformerWeights MNEBeamformerWeights::make_beamformer_weights()
{
}

//=============================================================================================================

MNEBeamformerWeights MNEBeamformerWeights::prepare_beamformer_weights(qint32 nave,
                                                                      float lambda2) const
{
    if(nave <= 0)
    {
        printf("The number of averages should be positive\n");
        return MNEBeamformerWeights();
    }
    printf("Preparing the inverse operator for use...\n");
    MNEBeamformerWeights bfw(*this);
    //
    //   Scale some of the stuff
    //
    float scale     = ((float)bfw.nave)/((float)nave);
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

}





