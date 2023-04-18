//=============================================================================================================
/**
 * @file     beamformer.cpp
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
 * @brief    Beamformer class definition.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "beamformer.h"

#include <mne/mne_sourceestimate.h>
#include <fiff/fiff_evoked.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================



//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace INVERSELIB;
using namespace FIFFLIB;
using namespace MNELIB;
using namespace Eigen;
using namespace UTILSLIB;

//=============================================================================================================
// DEFINE GLOBAL METHODS
//=============================================================================================================

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

Beamformer::Beamformer(const MNEBeamformerWeights &p_beamformerWeights, float p_fLambda, const QString p_sWeightnorm)
    : m_beamformerWeights(p_beamformerWeights),
      m_bBeamformerSetup(false)
{

    qDebug() << "[Beamformer::Beamformer] Created Beamformer object.";


    //HINT: this constructor is similar to the ones in MinimumNorm, but lambda is p_lambda here and method is used for weight normalization option
    //TODO: check whether we really need the regularization parameter here
    this->setRegularization(p_fLambda);
    this->setWeightnorm(p_sWeightnorm);

}

//=============================================================================================================

Beamformer::Beamformer(const MNEBeamformerWeights &p_beamformerWeights, const FiffInfo &p_dataInfo, const MNEForwardSolution &p_forward, const FiffCov &p_dataCov, const FiffCov &p_noiseCov, float p_fLambda, const QString p_sWeightnorm)
    : m_beamformerWeights(p_beamformerWeights)
    , m_bBeamformerSetup(false)
    , m_dataInfo(p_dataInfo)
    , m_noiseCov(p_noiseCov)
    , m_dataCov(p_dataCov)
    , m_forward(p_forward)

{

    qDebug() << "[Beamformer::Beamformer] Created Beamformer object. It is not set up yet.";


    //HINT: this constructor is similar to the ones in MinimumNorm, but lambda is p_lambda here and method is used for weight normalization option
    //TODO: check whether we really need the regularization parameter here
    this->setRegularization(p_fLambda);
    this->setWeightnorm(p_sWeightnorm);

}

//=============================================================================================================

MNESourceEstimate Beamformer::calculateInverse(const FiffEvoked &p_fiffEvoked, bool pick_normal)
{
    Q_UNUSED(pick_normal);

    //HINT: copied from MinimumNorm::calculateInverse for evoked input
    //HINT: modifications: names

    //
    //   Set up the inverse according to the parameters
    //
    qint32 nave = p_fiffEvoked.nave;
    qDebug() << "[Beamformer::calculateInverse] nave = " << nave;
    qDebug() << "[Beamformer::calculateInverse] pick_normal = " << pick_normal;

    //check whether channels are compatible
    if(!m_beamformerWeights.check_ch_names(p_fiffEvoked.info)) {
        qWarning("[Beamformer::calculateInverse] Channel name check failed.");
        return MNESourceEstimate();
    }

    //TODO: this doInverseSetup method does nothing
    //doInverseSetup(nave, pick_normal);

    if(!m_bBeamformerSetup)
    {
        qWarning("[Beamformer::calculateInverse] Beamformer not setup -> call doInverseSetup first! Return default MNESourceEstimate");
        return MNESourceEstimate();
    }

    //
    //   Pick the correct channels from the data
    //
    FiffEvoked t_fiffEvoked = p_fiffEvoked.pick_channels(m_beamformerWeightsSetup.noise_cov->names);

    qInfo("[Beamformer::calculateInverse] Picked %d channels from the data\n",t_fiffEvoked.info.nchan);

    //Results
    float tmin = p_fiffEvoked.times[0];
    float tstep = 1/t_fiffEvoked.info.sfreq;

    return calculateInverse(t_fiffEvoked.data, tmin, tstep, pick_normal);
}

//=============================================================================================================

MNESourceEstimate Beamformer::calculateInverse(const MatrixXd &data, float tmin, float tstep, bool pick_normal) const
{


    //TODO: this parameter is unused in this function body but necessary in signature for correct overloading (idea from rtc music implementation of calculateInverse)
    Q_UNUSED(pick_normal);

    //HINT: these ifs are copied from minimumnorm method calculateInverse and slightly adapted
    if(!m_bBeamformerSetup)
    {
        qWarning("[Beamformer::calculateInverse] Beamformer not setup -> call doInverseSetup first! Return default MNESourceEstimate");
        return MNESourceEstimate();
    }

    qDebug() << "[Beamformer::calculateInverse] dim m_matWTransposed: " << m_matWTSetup.rows() << " x " << m_matWTSetup.cols();
    qDebug() << "[Beamformer::calculateInverse] dim data: " << data.rows() << " x " << data.cols();


    if(m_matWTSetup.cols() != data.rows()) {
        qWarning() << "[Beamformer::calculateInverse] Dimension mismatch between m_matWTransposed.cols() and data.rows() -" << m_matWTSetup.cols() << "and" << data.rows() << ". Return default MNESourceEstimate.\n";
        return MNESourceEstimate();
    }

    //TODO: add some other output options (source power as estimated source activity in vecSourcePow of MNEBeamformer vs beamformer filter output)
    //these options should be user user adjustable similar to sLoreta etc methods




    //apply beamformer filter matrix to data to get filter output
    //output matrix has dimension (3*nsource x ntimes)
    MatrixXd sol = m_matWTSetup * data; //filter output
    std::cout << "[Beamformer::calculateInverse] Filter output dim: sol " << sol.rows() << " x " << sol.cols() << std::endl;


    //copied from MinimumNorm::calculateInverse
    if (!m_beamformerWeightsSetup.fixedOri && pick_normal == false)
    {
        printf("[Beamformer::calculateInverse] Combining the current components...\n");

        MatrixXd sol1(sol.rows()/3,sol.cols());
        for(qint32 i = 0; i < sol.cols(); ++i)
        {
            VectorXd* tmp = MNEMath::combine_xyz(sol.col(i));
            sol1.block(0,i,sol.rows()/3,1) = tmp->cwiseSqrt();
            delete tmp;
        }
        sol.resize(sol1.rows(),sol1.cols());
        sol = sol1;
    }




    qDebug() << "[Beamformer::calculateInverse] Data dim: " << data.rows() << " x " << data.cols();
    qDebug() << "[Beamformer::calculateInverse] m_matWTSetup dim: " << m_matWTSetup.rows() << " x " << m_matWTSetup.cols();
    qDebug() << "[Beamformer::calculateInverse] sol dim: " << sol.rows() << " x " << sol.cols();

    qDebug() << "[Beamformer::calculateInverse] data.maxCoeff() = " << data.maxCoeff();
    qDebug() << "[Beamformer::calculateInverse] data.minCoeff() = " << data.minCoeff();
    qDebug() << "[Beamformer::calculateInverse] data.mean() = " << data.mean();

    qDebug() << "[Beamformer::calculateInverse] m_matWTSetup.maxCoeff() = " << m_matWTSetup.maxCoeff();
    qDebug() << "[Beamformer::calculateInverse] m_matWTSetup.minCoeff() = " << m_matWTSetup.minCoeff();
    qDebug() << "[Beamformer::calculateInverse] m_matWTSetup.mean() = " << m_matWTSetup.mean();

    qDebug() << "[Beamformer::calculateInverse] sol.maxCoeff() = " << sol.maxCoeff();
    qDebug() << "[Beamformer::calculateInverse] sol.minCoeff() = " << sol.minCoeff();
    qDebug() << "[Beamformer::calculateInverse] sol.mean() = " << sol.mean();

    printf("[Beamformer::calculateInverse] Source estimate (beamformer filter output) computed. \n");


//    //Scale beamformer output so that values are displayed in the source estimation widget
//    //scaling values were determined empirically
//    double dScale = 1;
//    if(m_sWeightnorm == "no"){
//        //TODO fix this value
//        dScale = 50000;

//    }else if(m_sWeightnorm == "arraygain"){
//        //TODO: this value has to be set
//        dScale = 1e6;
//    }else{
//        qWarning() << "[Beamformer::calculateInverse] Unknown weight normalization method or method not debugged. No scaling of inverse solution. ";
//    }

//    sol *= dScale;
//    qDebug() << "[Beamformer::calculateInverse] Scaled inverse solution with " << dScale;



//    qDebug() << "[Beamformer::calculateInverse] sol.maxCoeff() after scaling = " << sol.maxCoeff();
//    qDebug() << "[Beamformer::calculateInverse] sol.minCoeff() after scaling= " << sol.minCoeff();
//    qDebug() << "[Beamformer::calculateInverse] sol.mean() after scaling= " << sol.mean();

    //Results
    //HINT: copied from calculateInverse of minimumnorm
    VectorXi p_vecVertices(m_beamformerWeightsSetup.src[0].vertno.size() + m_beamformerWeightsSetup.src[1].vertno.size());
    p_vecVertices << m_beamformerWeightsSetup.src[0].vertno, m_beamformerWeightsSetup.src[1].vertno;

    //for Debugging
    //return MNESourceEstimate(sol*(1e18), p_vecVertices, tmin, tstep);

    return MNESourceEstimate(sol, p_vecVertices, tmin, tstep);

}



//=============================================================================================================

void Beamformer::doInverseSetup(qint32 nave, bool pick_normal) //parameters are not used in this class
{
    Q_UNUSED(nave);
    Q_UNUSED(pick_normal);
    //HINT: copied from minimumnorm, adapted to beamformer weights preparation

    //TODO: check whether we need this pick_normal option for BF (does it make sense for BF?)

    //
    //   Set up the beamformer weights
    //
    //m_beamformerWeightsSetup = m_beamformerWeights.prepare_beamformer_weights();
    //TODOOO: start here next time?
/*    m_beamformerWeightsSetup = m_beamformerWeights.prepare_beamformer_weights(m_dataInfo, m_forward, m_dataCov, m_noiseCov, m_sWeightnorm);
    qInfo("[Beamformer::doInverseSetup] Prepared the beamformer weights.");

    //TODO: this is only for debugging
    m_matWTransposed = m_beamformerWeightsSetup.weights;
    std::cout << "[Beamformer::doInverseSetup] W^T dim: " << m_matWTransposed.rows() << " x " << m_matWTransposed.cols() << std::endl;
*/

    //prepare the beamformer filter matrix W^T for use by multiplying it by the whitener matrix
    //this preparation step is performed in order to apply the whitener to the input data too without having to access the whitener matrix from within the calculateInverse method called in the run method of the beamformer plug-in
    m_beamformerWeightsSetup = m_beamformerWeights;
    m_beamformerWeightsSetup.weights *= m_beamformerWeights.whitener;
    m_bBeamformerSetup = true;
    qDebug() << "[Beamformer::doInverseSetup] set m_bBeamformerSetup = true.";

    m_matWTSetup = m_beamformerWeightsSetup.weights;

    qDebug() << "[Beamformer::doInverseSetup] Set m_matWTSetup with dim: " << m_matWTSetup.rows() << " x " << m_matWTSetup.cols();



}

//=============================================================================================================

const char* Beamformer::getName() const
{
    return "Beamformer Source Estimate";
}

//=============================================================================================================

const MNESourceSpace& Beamformer::getSourceSpace() const
{
    return m_beamformerWeights.src;
}

//=============================================================================================================

const QString Beamformer::getWeightnorm()
{
    return m_sWeightnorm;
}

//=============================================================================================================

void Beamformer::setWeightnorm(QString weightnorm)
{
    if(weightnorm.compare("no") == 0)
        m_sWeightnorm = QString("no");
    else if(weightnorm.compare("unitnoisegain") == 0)
        m_sWeightnorm = QString("unitnoisegain");
    else if(weightnorm.compare("arraygain") == 0)
        m_sWeightnorm = QString("arraygain");
    else if(weightnorm.compare("nai") == 0)
        m_sWeightnorm = QString("nai");
    else
    {
        qWarning("[Beamformer::setWeightnorm] Weight normalization method not recognized! Activating no weight normalization.");
        m_sWeightnorm = QString("no");

    }

    printf("\t [Beamformer::setWeightnorm] Set weight normalization method to %s.\n", m_sWeightnorm.toUtf8().constData());
}

//=============================================================================================================

void Beamformer::setRegularization(float lambda)
{
    m_fLambda = lambda;
}
