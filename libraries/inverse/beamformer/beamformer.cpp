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

//=============================================================================================================
// DEFINE GLOBAL METHODS
//=============================================================================================================

//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

Beamformer::Beamformer(const MNEInverseOperator &p_inverseOperator, float lambda)
    : m_inverseOperator(p_inverseOperator)
    , inverseSetup(false)
{
    this->setRegularization(lambda);
}

//=============================================================================================================

MNESourceEstimate Beamformer::calculateInverse(const FiffEvoked &p_fiffEvoked, bool pick_normal)
{
    //HINT: copied from minimumnorm class

    //
    //   Set up the inverse according to the parameters
    //
    qint32 nave = p_fiffEvoked.nave;

    if(!m_inverseOperator.check_ch_names(p_fiffEvoked.info)) {
        qWarning("Channel name check failed.");
        return MNESourceEstimate();
    }

    doInverseSetup(nave,pick_normal);

    //
    //   Pick the correct channels from the data
    //
    FiffEvoked t_fiffEvoked = p_fiffEvoked.pick_channels(inv.noise_cov->names);

    printf("Picked %d channels from the data\n",t_fiffEvoked.info.nchan);

    //Results
    float tmin = p_fiffEvoked.times[0];
    float tstep = 1/t_fiffEvoked.info.sfreq;

    return calculateInverse(t_fiffEvoked.data, tmin, tstep, pick_normal);
}
//=============================================================================================================

void Beamformer::doInverseSetup(qint32 nave, bool pick_normal)
{

    //HINT: copied from minimumnorm

    //
    //   Set up the inverse according to the parameters
    //
    inv = m_inverseOperator.prepare_inverse_operator(nave, m_fLambda, m_bdSPM, m_bsLORETA);

    printf("Computing inverse...\n");
    inv.assemble_kernel(label, m_sMethod, pick_normal, K, noise_norm, vertno);

    std::cout << "K " << K.rows() << " x " << K.cols() << std::endl;

    inverseSetup = true;
}
