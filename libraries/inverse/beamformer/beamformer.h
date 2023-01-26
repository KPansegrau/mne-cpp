//=============================================================================================================
/**
 * @file     beamformer.h
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
 * @brief     Beamformer class declaration.
 *
 */

#ifndef BEAMFORMER_H
#define BEAMFORMER_H

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "../inverse_global.h"
#include "../IInverseAlgorithm.h"

#include <mne/mne_beamformer_weights.h>

//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QSharedPointer>

//=============================================================================================================
// EIGEN INCLUDES
//=============================================================================================================

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
// DEFINE NAMESPACE
//=============================================================================================================

namespace INVERSELIB
{

//=============================================================================================================
// FORWARD DECLARATIONS
//=============================================================================================================

//=============================================================================================================
/**
 * TODO: edit docu
 *
 * Beamformer source estimation algorithm ToDo: Paper references.
 *
 * @brief Beamformer source estimation.
 */
class INVERSESHARED_EXPORT Beamformer : public IInverseAlgorithm
{

public:

    typedef QSharedPointer<Beamformer> SPtr;            /**< Shared pointer type for Beamformer. */
    typedef QSharedPointer<const Beamformer> ConstSPtr; /**< Const shared pointer type for Beamformer. */

    //=========================================================================================================
    // IInverseAlgorithm methods
    //=========================================================================================================


    /**
     * TODO: edit docu
     *
     * Constructs beamformer inverse algorithm
     *
     *
     */

    explicit Beamformer(const MNEInverseOperator &p_inverseOperator, float lambda);

    //=========================================================================================================
    /**
     *
     *
     * Destructor for beamformer algorithm
     *
     */

    virtual ~Beamformer(){}

    //=========================================================================================================
    /**
     * TODO: edit docu
     *
     * @param[in] p_fiffEvoked   Evoked data.
     * @param[in] pick_normal    If True, rather than pooling the orientations by taking the norm, only the.
     *                           radial component is kept. This is only applied when working with loose orientations.
     *
     * @return the calculated source estimation.
     */
    virtual MNELIB::MNESourceEstimate calculateInverse(const FIFFLIB::FiffEvoked &p_fiffEvoked, bool pick_normal = false);

    //=========================================================================================================
    /**
     * TODO: edit docu
     *
     * @param[in] p_fiffEvoked   Evoked data.
     * @param[in] pick_normal    If True, rather than pooling the orientations by taking the norm, only the.
     *                           radial component is kept. This is only applied when working with loose orientations.
     *
     * @return the calculated source estimation.
     */

    virtual MNELIB::MNESourceEstimate calculateInverse(const Eigen::MatrixXd &data, float tmin, float tstep, bool pick_normal = false) const;

    //=========================================================================================================
    /**
     * TODO: edit docu
     *
     * Get the name of the inverse beamformer algorithm.
     *
     * @return the name of the inverse beamformer algorithm.
     */
    virtual const char* getName() const;

    //=========================================================================================================
    /**
     * TODO: edit docu
     *
     * Get the source space corresponding to this inverse operator.
     *
     * @return the source space corresponding to this inverse operator.
     */
    virtual const MNELIB::MNESourceSpace& getSourceSpace() const;

    //=========================================================================================================
    // further member methods
    //=========================================================================================================

    /**
     * TODO: edit docu
     *
     * Perform the inverse setup: Prepares this inverse operator and assembles the kernel.
     *
     * @param[in] nave           Number of averages to use.
     * @param[in] pick_normal    If True, rather than pooling the orientations by taking the norm, only the.
     *                           radial component is kept. This is only applied when working with loose orientations.
     */
    virtual void doInverseSetup(qint32 nave, bool pick_normal = false);


protected:

private:

};

//=============================================================================================================
// INLINE DEFINITIONS
//=============================================================================================================

} //namespace
#endif // BEAMFORMER_H

