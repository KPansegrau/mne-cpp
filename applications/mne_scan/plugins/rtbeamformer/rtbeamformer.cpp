//=============================================================================================================
/**
 * @file     rtbeamformer.cpp
 * @author   Kerstin Pansegrau <kerstin.pansegrau@tu-ilmenau.de
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
 * @brief    Definition of the RtBeamformer class.
 *
 */

//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "rtbeamformer.h"




//=============================================================================================================
// QT INCLUDES
//=============================================================================================================



//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace RTBEAMFORMERPLUGIN;
using namespace SCSHAREDLIB;


//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

RtBeamformer::RtBeamformer()
{

}

//=============================================================================================================

RtBeamformer::~RtBeamformer()
{

}

//=============================================================================================================

QSharedPointer<AbstractPlugin> RtBeamformer::clone() const
{
    QSharedPointer<RtBeamformer> pRtBeamformerClone(new RtBeamformer());
    return pRtBeamformerClone;
}

//=============================================================================================================

void RtBeamformer::init()
{

}

//=============================================================================================================

void RtBeamformer::unload()
{

}

//=============================================================================================================

bool RtBeamformer::start()
{
    return true;
}

//=============================================================================================================

bool RtBeamformer::stop()
{
    return true;
}

//=============================================================================================================

AbstractPlugin::PluginType RtBeamformer::getType() const
{
    return _IAlgorithm;
}

//=============================================================================================================

QString RtBeamformer::getName() const
{
    return "Beamformer Source Localization";
}

//=============================================================================================================

QWidget* RtBeamformer::setupWidget()
{
    //placeholder. This is where the UI will be setup later
    QWidget* setupWidget = new QWidget();
    return setupWidget;
}

//=============================================================================================================

void RtBeamformer::initPluginControlWidgets()
{

}

//=============================================================================================================

void RtBeamformer::run()
{

}

//=============================================================================================================

QString RtBeamformer::getBuildInfo()
{
    return QString(RTBEAMFORMERPLUGIN::buildDateTime()) + QString(" - ")  + QString(RTBEAMFORMERPLUGIN::buildHash());

}



