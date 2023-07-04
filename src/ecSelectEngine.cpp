/* **********************************************************************
Copyright 2020  Dr. J. C. Shillcock and Prof. Dr. R. Lipowsky, Director at the Max Planck Institute (MPI) of Colloids and Interfaces; Head of Department Theory and Bio-Systems.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************** */
// ecSelectEngine.cpp: implementation of the ecSelectEngine class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "SimDefs.h"
#include "ecSelectEngine.h"
#include "ISimCmd.h"
#include "InputData.h"

#include "IIntegrationEngine.h"

//////////////////////////////////////////////////////////////////////
// Global members
//////////////////////////////////////////////////////////////////////

// Static member variable containing the identifier for this command. 
// The static member function GetType() is invoked by the xxCommandObject 
// to compare the type read from the control data file with each
// xxCommand-derived class so that it can create the appropriate object 
// to hold the command data.

const zString ecSelectEngine::m_Type = "SelectSimEngine";

const zString ecSelectEngine::GetType()
{
	return m_Type;
}

// We use an anonymous namespace to wrap the call to the factory object
// so that it is not accessible from outside this file. The identifying
// string for the command is stored in the m_Type static member variable.
//
// Note that the Create() function is not a member function of the
// command class but a global function hidden in the namespace.

namespace
{
	xxCommand* Create(long executionTime) {
		return new ecSelectEngine(executionTime);
	}

	const zString id = ecSelectEngine::GetType();

	const bool bRegistered = acfCommandFactory::Instance()->Register(id, Create);
}


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ecSelectEngine::ecSelectEngine(long executionTime) : xxCommand(executionTime)
{
}

// Constructor for use when executing the command internally. We use a flag 
// to indicate if the command's execution should be logged or not.

ecSelectEngine::ecSelectEngine(long executionTime, bool bLog) : xxCommand(executionTime, bLog)
{
}

ecSelectEngine::ecSelectEngine(const ecSelectEngine& oldCommand)
	: xxCommand(oldCommand)
	, m_engine_name(oldCommand.m_engine_name)
{
}

ecSelectEngine::~ecSelectEngine()
{
}

// Member functions to read/write the data specific to the command.
//
// Arguments
// *********
//
//

long ecSelectEngine::GetCommandArgumentTotal() const
{
	return 1;
}

zOutStream& ecSelectEngine::put(zOutStream& os) const
{
#if EnableXMLCommands == SimXMLEnabled

	// XML output
	putXMLStartTags(os);
	putXMLEndTags(os);

#elif EnableXMLCommands == SimXMLDisabled

	// ASCII output 
	putASCIIStartTags(os);
	os << "    EngineName " << m_engine_name << zEndl;
	putASCIIEndTags(os);

#endif

	return os;
}

zInStream& ecSelectEngine::get(zInStream& is)
{
	return is >> m_engine_name;
}

// Non-static function to return the type of the command

const zString ecSelectEngine::GetCommandType() const
{
	return m_Type;
}

// Function to return a pointer to a copy of the current command.

const xxCommand* ecSelectEngine::GetCommand() const
{
	return new ecSelectEngine(*this);
}


// Implementation of the command that is sent by the SimBox to each xxCommand
// object to see if it is the right time for it to carry out its operation.
// We return a boolean so that the SimBox can see if the command executed or not
// as this may be useful for considering several commands. 

bool ecSelectEngine::Execute(long simTime, ISimCmd* const pISimCmd) const
{
	if(simTime == GetExecutionTime())
	{
		auto inst = IIntegrationEngineFactory::CreateEngineInstanceByName(m_engine_name);
		if(!inst){
			ErrorTrace("Couldn't find engine called "+m_engine_name);
		}

		IIntegrationEngine::SetGlobalEngine(inst );
		return true;
	}
	else
		return false;
}

// Function to check that the command data are valid: no data are required for this command.

bool ecSelectEngine::IsDataValid(const CInputData& riData) const
{
	if(IIntegrationEngineFactory::IsKnownEngine(m_engine_name)){
		return true;
	}
	ErrorTrace("SelectSimEngine : Unknown SimEngine "+m_engine_name);
	return false;
}
