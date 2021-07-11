/* **********************************************************************
Copyright 2020  Dr. J. C. Shillcock and Prof. Dr. R. Lipowsky, Director at the Max Planck Institute (MPI) of Colloids and Interfaces; Head of Department Theory and Bio-Systems.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************** */
// ccUnchargeBeadType.cpp: implementation of the ccUnchargeBeadType class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "SimDefs.h"
#include "ccUnchargeBeadType.h"
#include "ISimCmd.h"
#include "InputData.h"
#include "Builder.h"		// Needed for pi

//////////////////////////////////////////////////////////////////////
// Global members
//////////////////////////////////////////////////////////////////////

// Static member variable containing the identifier for this command. 
// The static member function GetType() is invoked by the xxCommandObject 
// to compare the type read from the control data file with each
// xxCommand-derived class so that it can create the appropriate object 
// to hold the command data.

const zString ccUnchargeBeadType::m_Type = "UnchargeBeadType";

const zString ccUnchargeBeadType::GetType()
{
	return m_Type;
}

// We use an anonymous namespace to wrap the call to the factory object
// so that it is not accessible from outside this file. The identifying
// string for the command is stored in the m_Type static member variable.
//
// Note that the Create() function is not a member function of the
// command class but a global function hidden in the namespace.

namespace private_ns_ccUnchargeBeadType
{
	xxCommand* Create(long executionTime) {return new ccUnchargeBeadType(executionTime);}

	const zString id = ccUnchargeBeadType::GetType();

	const bool bRegistered = acfCommandFactory::Instance()->Register(id, Create);
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ccUnchargeBeadType::ccUnchargeBeadType(long executionTime) : xxCommand(executionTime),
									m_Name("")
{
}

ccUnchargeBeadType::ccUnchargeBeadType(const ccUnchargeBeadType& oldCommand) : xxCommand(oldCommand),
									 m_Name(oldCommand.m_Name)
{
}

ccUnchargeBeadType::~ccUnchargeBeadType()
{
}

// Member functions to read/write the data specific to the command.
//
// Arguments
// *********
//
//	Name				Bead type name of charged species

zOutStream& ccUnchargeBeadType::put(zOutStream& os) const
{
#if EnableXMLCommands == SimXMLEnabled

	// XML output
	putXMLStartTags(os);
	os << "<Name>" << m_Name << "</Name>" << zEndl;
	putXMLEndTags(os);

#elif EnableXMLCommands == SimXMLDisabled

	// ASCII output 
	putASCIIStartTags(os);
	os << m_Name;
	putASCIIEndTags(os);

#endif

	return os;
}

zInStream& ccUnchargeBeadType::get(zInStream& is)
{
	// We check that the name corresponds to a valid bead type in 
	// IsDataValid() as the CInputData object has to first create the map 
	// holding the valid bead type names.

	is >> m_Name;

	if(!is.good())
	   SetCommandValid(false);


	return is;
}

// Non-static function to return the type of the command

const zString ccUnchargeBeadType::GetCommandType() const
{
	return m_Type;
}

// Function to return a pointer to a copy of the current command.

const xxCommand* ccUnchargeBeadType::GetCommand() const
{
	return new ccUnchargeBeadType(*this);
}


// Implementation of the command that is sent by the SimBox to each xxCommand
// object to see if it is the right time for it to carry out its operation.
// We return a boolean so that the SimBox can see if the command executed or not
// as this may be useful for considering several commands. 

bool ccUnchargeBeadType::Execute(long simTime, ISimCmd* const pISimCmd) const
{
	if(simTime == GetExecutionTime())
	{
		pISimCmd->UnchargeBeadType(this);
		return true;
	}
	else
		return false;
}

// Function to check that the bead type name is valid. Note that 
// we have already validated the force's strength and range above.

bool ccUnchargeBeadType::IsDataValid(const CInputData& riData) const
{
	// Check bead type name is unique and consists only of allowed characters.

	if(!riData.IsExternalNameValid(m_Name))
		return ErrorTrace("Invalid Bead type name");

	return true;
}
