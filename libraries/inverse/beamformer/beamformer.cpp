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
    //HINT: modifications: names, no channel checking

    //
    //   Set up the inverse according to the parameters
    //
    qint32 nave = p_fiffEvoked.nave;

    //HINT: this part not necessary since channel checking was already performed during computation of beamformer weights
/*    if(!m_beamformerWeights.check_ch_names(p_fiffEvoked.info)) {
        qWarning("Channel name check failed.");
        return MNESourceEstimate();
    }
*/
    doInverseSetup(nave,pick_normal); //HINT: this parameters are not needed in doInverseSetup (but have to be function parameters since this is a virtual function)

    //
    //   Pick the correct channels from the data
    //
    FiffEvoked t_fiffEvoked = p_fiffEvoked.pick_channels(m_beamformerWeightsSetup.noise_cov->names);

    printf("Picked %d channels from the data\n",t_fiffEvoked.info.nchan);

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
        qWarning("Beamformer::calculateInverse - Beamformer not setup -> call doInverseSetup first! Return default MNESourceEstimate");
        return MNESourceEstimate();
    }

    if(m_W_transposed.cols() != data.rows()) {
        qWarning() << "Beamformer::calculateInverse - Dimension mismatch between m_W_transposed.cols() and data.rows() -" << m_W_transposed.cols() << "and" << data.rows() << ". Return default MNESourceEstimate.\n";
        return MNESourceEstimate();
    }

    //TODO: add some other output options (source power as estimated source activity in vecSourcePow of MNEBeamformer vs beamformer filter output)
    //these options should be user user adjustable similar to sLoreta etc methods
    //maybe in do Inverse setup and then if statement here (if filter output sol = Wt*data, if activity strengh sol = vecSourcePow etc)

    //TODO: we need to apply the SSP projector to the data first, do it here!
    printf("Beamformer::calculateInverse TODO: Apply SSP projector to data...\n");

    //apply beamformer filter matrix to data to get filter output
    //output matrix has dimension (3*nsource x ntimes)
    MatrixXd sol = m_W_transposed * data; //filter output
    std::cout << "Beamformer::calculateInverse: Filter output dimension: sol " << sol.rows() << " x " << sol.cols() << std::endl;

    printf("Beamformer::calculateInverse Source estimate (beamformer filter output) computed. \n");

    //Results
    //HINT: copied from calculateInverse of minimumnorm
    VectorXi p_vecVertices(m_beamformerWeightsSetup.src[0].vertno.size() + m_beamformerWeightsSetup.src[1].vertno.size());
    p_vecVertices << m_beamformerWeightsSetup.src[0].vertno, m_beamformerWeightsSetup.src[1].vertno;


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
    m_beamformerWeightsSetup = m_beamformerWeights.prepare_beamformer_weights();
    qInfo("Beamformer::doInverseSetup Prepared the beamformer weights.");

    //TODO: do we need a method assemble beamformer weights here similar to assemble kernel? (this should handle pick normal?)

    m_W_transposed = m_beamformerWeightsSetup.weights;
    std::cout << "Beamformer::doInverseSetup: W^T " << m_W_transposed.rows() << " x " << m_W_transposed.cols() << std::endl;

    m_bBeamformerSetup = true;
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
        qWarning("Weight normalization method not recognized! Activating no weight normalization.");
        m_sWeightnorm = QString("no");

    }

    printf("\tSet weight normalization method to %s.\n", m_sWeightnorm.toUtf8().constData());
}

//=============================================================================================================

void Beamformer::setRegularization(float lambda)
{
    m_fLambda = lambda;
}
