/* **********************************************************************
Copyright 2020  Dr. J. C. Shillcock and Prof. Dr. R. Lipowsky, Director at the Max Planck Institute (MPI) of Colloids and Interfaces; Head of Department Theory and Bio-Systems.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************** */
// pmSetDPDBeadConsIntByType.cpp: implementation of the pmSetDPDBeadConsIntByType class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "SimDefs.h"
#include "pmSetDPDBeadConsIntByType.h"
#include "ccSetDPDBeadConsIntByType.h"

//////////////////////////////////////////////////////////////////////
// Global members
//////////////////////////////////////////////////////////////////////

// Static member variable containing the identifier for this message. 

const zString pmSetDPDBeadConsIntByType::m_Type = "SetDPDBeadConsIntByType";

const zString pmSetDPDBeadConsIntByType::GetType()
{
	return m_Type;
}

// We use an anonymous namespace to wrap the call to the factory object
// so that it is not accessible from outside this file. The identifying
// string for the command is stored in the m_Type static member variable.
//
// Note that the Create() function is not a member function of the
// message class but a global function hidden in the namespace.

namespace private_ns_pmSetDPDBeadConsIntByType
{
	mpmMessage* Create() {return new pmSetDPDBeadConsIntByType();}

	const zString id = pmSetDPDBeadConsIntByType::GetType();

	const bool bRegistered = acfParallelMessageFactory::Instance()->Register(id, Create);
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

pmSetDPDBeadConsIntByType::pmSetDPDBeadConsIntByType() : mpmCommandWrapper(), 
											     m_FirstType(-1),
											     m_SecondType(-1),
											     m_ConsInt(0.0)
{
}

pmSetDPDBeadConsIntByType::pmSetDPDBeadConsIntByType(const pmSetDPDBeadConsIntByType& oldMessage) : mpmCommandWrapper(oldMessage),
											     m_FirstType(oldMessage.m_FirstType),
											     m_SecondType(oldMessage.m_SecondType),
											     m_ConsInt(oldMessage.m_ConsInt)
{
}

pmSetDPDBeadConsIntByType::~pmSetDPDBeadConsIntByType()
{
}

// Member functions to read/write the data specific to the message.
//
// Arguments
// *********
//

zOutStream& pmSetDPDBeadConsIntByType::put(zOutStream& os) const
{
#if EnableXMLParallelMessagess == SimXMLEnabled

	// XML output
	putXMLStartTags(os);
	os << "<MessageId>"    << GetId()   << "</MessageId>" << zEndl;
	putXMLEndTags(os);

#elif EnableXMLParallelMessagess == SimXMLDisabled

	// ASCII output 
	putASCIIStartTags(os);
	os << GetId();
	putASCIIEndTags(os);

#endif

	return os;
}

zInStream& pmSetDPDBeadConsIntByType::get(zInStream& is)
{
	return is;
}

// Non-static function to return the type of the message

const zString pmSetDPDBeadConsIntByType::GetMessageType() const
{
	return m_Type;
}

// Function to return a pointer to a copy of the current message.

const mpmMessage* pmSetDPDBeadConsIntByType::GetMessage() const
{
	return new pmSetDPDBeadConsIntByType(*this);
}

// ****************************************
// Function to copy all the required data for the command into the message.
// Because each command message wraps a known command we can extract the
// data using dynamic cast.

void pmSetDPDBeadConsIntByType::SetMessageData(const xxCommand* const pCommand)
{
	const ccSetDPDBeadConsIntByType* const pCmd = dynamic_cast<const ccSetDPDBeadConsIntByType*>(pCommand);

    if(pCmd)
    {
        m_FirstType  = pCmd->GetFirstType();
        m_SecondType = pCmd->GetSecondType();
		m_ConsInt    = pCmd->GetConsInt();
    }
}

// ****************************************
// Pure virtual function to check that the data are valid prior to sending a message.
// Note that we allow any real value for the conservative interaction parameter.

bool pmSetDPDBeadConsIntByType::Validate()
{
    bool bSuccess = true;

    if(m_FirstType < 0 || m_SecondType < 0)
    {
        bSuccess = false;
    }

    return bSuccess;
}

// Wrapper for the MPI function:
// MPI_Send(message, length, MPI_CHAR, m_ProcId, tag, MPI_COMM_WORLD);

void pmSetDPDBeadConsIntByType::SendAllP()
{
#if SimMPS == SimulationEnabled

    // The data consist of 1 integer and two doubles, so we pack them
    // sequentially in the buffer and send them to all the other processors: 
	// note that we assume that the sending processor is P0, so we start the 
	// loop at processor rank 1.

    char buffer[100];
    int  position;
	position = 0;

    MPI_Pack(&m_FirstType,   1, MPI_LONG,   buffer, 100, &position, MPI_COMM_WORLD);
    MPI_Pack(&m_SecondType,  1, MPI_LONG,   buffer, 100, &position, MPI_COMM_WORLD);
    MPI_Pack(&m_ConsInt,     1, MPI_DOUBLE, buffer, 100, &position, MPI_COMM_WORLD);

    for(m_ProcId=1; m_ProcId<GetWorld(); m_ProcId++)
    {
        MPI_Send(buffer, 100, MPI_PACKED, m_ProcId, 0, MPI_COMM_WORLD);
   }

#endif
}

// Wrapper for the Mpi function:
// MPI_Recv(message, length, MPI_CHAR, m_ProcId, tag, MPI_COMM_WORLD);

// We retrieve the data from the message and fill up this instance's member variables. 
// We do not check for validity here as we assume that the sending object has 
// already done it. Note that this function does not propagate the data out into 
// the code, it only fills up the receiving message instance. Each messaging class 
// is responsible for providing further functions to pass the data to the rest of the code.

void pmSetDPDBeadConsIntByType::Receive()
{
#if SimMPS == SimulationEnabled

    char buffer[100];
    int  position;

    int tag = 0;
    MPI_Status status;

    MPI_Recv(buffer, 100, MPI_PACKED, 0, tag, MPI_COMM_WORLD, &status);

    // Unpack the data into member variables

    position = 0;
    MPI_Unpack(buffer, 100, &position, &m_FirstType,   1, MPI_LONG,   MPI_COMM_WORLD);
    MPI_Unpack(buffer, 100, &position, &m_SecondType,  1, MPI_LONG,   MPI_COMM_WORLD);
    MPI_Unpack(buffer, 100, &position, &m_ConsInt,     1, MPI_DOUBLE, MPI_COMM_WORLD);

    // Now instantiate the concrete command class using the internal constructor that
    // takes the command's data as arguments
	
    m_pCommand = new ccSetDPDBeadConsIntByType(m_SimTime, m_FirstType, m_SecondType, m_ConsInt);

#endif
}

// Stub function for the overloaded Receive function to avoid compiler warnings on the Hermit1.hww.de platform.

void pmSetDPDBeadConsIntByType::Receive(long procId)
{
#if SimMPS == SimulationEnabled

#endif
}

// PVF to assemble the disparate data items into a struct suitable for sending to 
// the MPI messaging functions. All the data must be in the calling object, and the
// data type is not usable outside this class.
// The SimMPS flag must be defined for this function to work, otherwise it does nothing.

void pmSetDPDBeadConsIntByType::BuildDerivedType(MsgDataTypePtr pMsgDataType)
{
}

