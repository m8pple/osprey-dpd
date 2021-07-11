/* **********************************************************************
Copyright 2020  Dr. J. C. Shillcock and Prof. Dr. R. Lipowsky, Director at the Max Planck Institute (MPI) of Colloids and Interfaces; Head of Department Theory and Bio-Systems.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************** */
// ccRemovefActinBond.cpp: implementation of the ccRemovefActinBond class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "SimDefs.h"
#include "ccRemovefActinBond.h"
#include "ISimCmd.h"
#include "InputData.h"

//////////////////////////////////////////////////////////////////////
// Global members
//////////////////////////////////////////////////////////////////////

// Static member variable containing the identifier for this command. 
// The static member function GetType() is invoked by the xxCommandObject 
// to compare the type read from the control data file with each
// xxCommand-derived class so that it can create the appropriate object 
// to hold the command data.

const zString ccRemovefActinBond::m_Type = "RemovefActinBond";

const zString ccRemovefActinBond::GetType()
{
	return m_Type;
}

// We use an anonymous namespace to wrap the call to the factory object
// so that it is not accessible from outside this file. The identifying
// string for the command is stored in the m_Type static member variable.
//
// Note that the Create() function is not a member function of the
// command class but a global function hidden in the namespace.

namespace private_ns_ccRemovefActinBond
{
	xxCommand* Create(long executionTime) {return new ccRemovefActinBond(executionTime);}

	const zString id = ccRemovefActinBond::GetType();

#if EnableACNCommand == SimCommandEnabled
	const bool bRegistered = acfCommandFactory::Instance()->Register(id, Create);
#endif
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ccRemovefActinBond::ccRemovefActinBond(long executionTime) : xxCommand(executionTime),
															 m_HeadMonomerName(""),
															 m_TailMonomerName("")
{
}

ccRemovefActinBond::ccRemovefActinBond(const ccRemovefActinBond& oldCommand) : xxCommand(oldCommand),
											m_HeadMonomerName(oldCommand.m_HeadMonomerName),
											m_TailMonomerName(oldCommand.m_TailMonomerName)
{
}


ccRemovefActinBond::~ccRemovefActinBond()
{
}

// Member functions to read/write the data specific to the command.
//
// Arguments
// *********
//
//	m_HeadMonomerName		- Name of polymer type that is bound to this active bond
//	m_TailMonomerName		- Name of polymer type that is wrapped by this active bond

zOutStream& ccRemovefActinBond::put(zOutStream& os) const
{
#if EnableACNCommand == SimCommandEnabled
#if EnableXMLCommands == SimXMLEnabled

	// XML output
	putXMLStartTags(os);
	os << "<HeadMonomerName>" << m_HeadMonomerName << "</HeadMonomerName>" << zEndl;
	os << "<TailMonomerName>" << m_TailMonomerName << "</TailMonomerName>" << zEndl;
	putXMLEndTags(os);

#elif EnableXMLCommands == SimXMLDisabled

	// ASCII output 
	putASCIIStartTags(os);
	os << " " << m_HeadMonomerName << " " << m_TailMonomerName;
	putASCIIEndTags(os);

#endif
#endif

	return os;
}

zInStream& ccRemovefActinBond::get(zInStream& is)
{
#if EnableACNCommand == SimCommandEnabled
	is >> m_HeadMonomerName >> m_TailMonomerName;

	if(!is.good())
	   SetCommandValid(false);
#endif

	return is;
}

// Non-static function to return the type of the command

const zString ccRemovefActinBond::GetCommandType() const
{
	return m_Type;
}

// Function to return a pointer to a copy of the current command.

const xxCommand* ccRemovefActinBond::GetCommand() const
{
	return new ccRemovefActinBond(*this);
}

// Implementation of the command that is sent by the SimBox to each xxCommand
// object to see if it is the right time for it to carry out its operation.
// We return a boolean so that the SimBox can see if the command executed or not
// as this may be useful for considering several commands. 

bool ccRemovefActinBond::Execute(long simTime, ISimCmd* const pISimCmd) const
{
	if(simTime == GetExecutionTime())
	{
#if EnableACNCommand == SimCommandEnabled
		pISimCmd->RemovefActinBond(this);
#endif
		return true;
	}
	else
		return false;
}

// Function to check that the data defining the active bond type is valid.
// We only have to check that the polymer name specified as the monomers
// that are wrapped by the active bonds is valid. Because we are removing
// the active bond type, we don't need to read or validate the bond force data.

bool ccRemovefActinBond::IsDataValid(const CInputData& riData) const
{
	if(!riData.IsExternalNameValid(m_HeadMonomerName))
		return ErrorTrace("Invalid active bond head monomer name");
	else if(!riData.IsPolymerinMap(m_HeadMonomerName))
		return ErrorTrace("Active bond head monomer type does not exist");
	else if(!riData.IsExternalNameValid(m_TailMonomerName))
		return ErrorTrace("Invalid active bond tail monomer name");
	else if(!riData.IsPolymerinMap(m_TailMonomerName))
		return ErrorTrace("Active bond tail monomer type does not exist");

	return true;
}
