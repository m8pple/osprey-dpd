/* **********************************************************************
Copyright 2020  Dr. J. C. Shillcock and Prof. Dr. R. Lipowsky, Director at the Max Planck Institute (MPI) of Colloids and Interfaces; Head of Department Theory and Bio-Systems.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************** */
// IModifySimStateIntegration.cpp: implementation of the IModifySimStateIntegration class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "SimDefs.h"
#include "IModifySimStateIntegration.h" 

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
// Standard constructor

IModifySimStateIntegration::IModifySimStateIntegration(CSimState& simState) : IModifySimState(simState)
{

}

IModifySimStateIntegration::~IModifySimStateIntegration()
{

}

// Function to change the integration time step in response to 
// a command. It passes the command on to the CSimState.

void IModifySimStateIntegration::SetTimeStep(double dt)
{
	m_rSimState.SetTimeStep(dt);
}

double IModifySimStateIntegration::GetTimeStep() const
{
	return m_rSimState.GetTimeStep();
}

// Function to extend the total simulation time in response to
// a command. It passed the command on to the CSimState.
//
// The argument specifies the NEW total simulation time.

void IModifySimStateIntegration::SetTotalTime(long newTime)
{
	m_rSimState.SetTotalTime(newTime);
}

// Function to toggle the DPD bead interactions on and off.

bool IModifySimStateIntegration::ToggleDPDBeadForces()
{
	return m_rSimState.ToggleDPDBeadForces();
}

// Function to toggle the DPD bead conservative interactions on and off.

bool IModifySimStateIntegration::ToggleDPDBeadConservativeForces()
{
	return m_rSimState.ToggleDPDBeadConservativeForces();
}

// Function to toggle the DPD bead random and dissipative interactions on and off.

bool IModifySimStateIntegration::ToggleDPDBeadThermostat()
{
	return m_rSimState.ToggleDPDBeadThermostat();
}
