/* **********************************************************************
Copyright 2020  Dr. J. C. Shillcock and Prof. Dr. R. Lipowsky, Director at the Max Planck Institute (MPI) of Colloids and Interfaces; Head of Department Theory and Bio-Systems.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************** */
// ccOptimisePolymerOrderingForPDPD.cpp: implementation of the ccOptimisePolymerOrderingForPDPD class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "SimDefs.h"
#include "ccOptimisePolymerOrderingForPDPD.h"
#include "ISimCmd.h"
#include "InputData.h"

#include "AbstractBead.h"
#include "Bead.h"
#include "SimBox.h"
#include "ISimBox.h"
#include "Polymer.h"
#include "Bond.h"
#include "BondPair.h"

#include <cmath>
#include <cassert>

#include <unordered_map>

//////////////////////////////////////////////////////////////////////
// Global members
//////////////////////////////////////////////////////////////////////

// Static member variable containing the identifier for this command. 
// The static member function GetType() is invoked by the xxCommandObject 
// to compare the type read from the control data file with each
// xxCommand-derived class so that it can create the appropriate object 
// to hold the command data.

long ccOptimisePolymerOrderingForPDPD::GetArgumentTotal()
{
	return 0;
}

const zString ccOptimisePolymerOrderingForPDPD::m_Type = "OptimisePolymerOrderingForPDPD";

const zString ccOptimisePolymerOrderingForPDPD::GetType()
{
	return m_Type;
}

// We use an anonymous namespace to wrap the call to the factory object
// so that it is not accessible from outside this file. The identifying
// string for the command is stored in the m_Type static member variable.
//
// Note that the Create() function is not a member function of the
// command class but a global function hidden in the namespace.

namespace private_ns_ccOptimisePolymerOrderingForPDPD
{
	xxCommand* Create(long executionTime) {return new ccOptimisePolymerOrderingForPDPD(executionTime);}

	const zString id = ccOptimisePolymerOrderingForPDPD::GetType();

	const bool bRegistered = acfCommandFactory::Instance()->Register(id, Create);
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ccOptimisePolymerOrderingForPDPD::ccOptimisePolymerOrderingForPDPD(long executionTime) : xxCommand(executionTime)
{
}

ccOptimisePolymerOrderingForPDPD::ccOptimisePolymerOrderingForPDPD(const ccOptimisePolymerOrderingForPDPD& oldCommand) : xxCommand(oldCommand)
{
}


ccOptimisePolymerOrderingForPDPD::~ccOptimisePolymerOrderingForPDPD()
{

}

// Member functions to write/read the data specific to the command.
// Note that spaces are written here not in the xxCommand function.
//
// Arguments
// *********


zOutStream& ccOptimisePolymerOrderingForPDPD::put(zOutStream& os) const
{
#if EnableXMLCommands == SimXMLEnabled

	// XML output
	putXMLStartTags(os);
	putXMLEndTags(os);

#elif EnableXMLCommands == SimXMLDisabled

	// ASCII output 
	putASCIIStartTags(os);
	putASCIIEndTags(os);

#endif

	return os;
}

zInStream& ccOptimisePolymerOrderingForPDPD::get(zInStream& is)
{

	if(!is.good())
		SetCommandValid(false);

	return is;
}


// Implementation of the command that is sent by the SimBox to each xxCommand
// object to see if it is the right time for it to carry out its operation.
// We return a boolean so that the SimBox can see if the command executed or not
// as this may be useful for considering several commands.



bool ccOptimisePolymerOrderingForPDPD::Execute(long simTime, ISimCmd* const pISimCmd) const
{
	if(simTime != GetExecutionTime()){
		return false;
	}

	CSimBox *box=dynamic_cast<CSimBox*>(pISimCmd);
	if(!box){
		fprintf(stderr, "Export only works with basic sim box.");
		exit(1);
	}

	fprintf(stderr, "ccOptimisePolymerOrderingForPDPD\n");

	auto beads=box->GetSimBoxBeads();

	std::vector<CPolymer*> polymers;
	std::vector<CPolymer*> monomers;

	for(CPolymer *p : box->GetAllPolymers()) {
		if(p->GetSize()==1){
			monomers.push_back(p);
		}else{
			polymers.push_back(p);
		}
	}

	fprintf(stderr, "Found %u monomers and %u polymers.\n", (unsigned)monomers.size(), (unsigned)polymers.size());

	for(unsigned i=0; i<polymers.size(); i++){
		polymers[i]->SetId(i+1);
	}
	for(unsigned i=0; i<monomers.size(); i++){
		monomers[i]->SetId(i+1+polymers.size());
	}


	return true;
}

const zString ccOptimisePolymerOrderingForPDPD::GetCommandType() const
{
	return m_Type;
}

// Function to return a pointer to a copy of the current command.

const xxCommand* ccOptimisePolymerOrderingForPDPD::GetCommand() const
{
	return new ccOptimisePolymerOrderingForPDPD(*this);
}

// Function to check that the command is valid.
// TODO

bool ccOptimisePolymerOrderingForPDPD::IsDataValid(const CInputData &riData) const
{
	return true;
}

