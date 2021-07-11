/* **********************************************************************
Copyright 2020  Dr. J. C. Shillcock and Prof. Dr. R. Lipowsky, Director at the Max Planck Institute (MPI) of Colloids and Interfaces; Head of Department Theory and Bio-Systems.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************** */
// ccAddBeadsInRow.cpp: implementation of the ccAddBeadsInRow class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "SimDefs.h"
#include "SimAlgorithmFlags.h"
#include "ccAddBeadsInRow.h"
#include "ccSelectBeadsInRow.h"
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

const zString ccAddBeadsInRow::m_Type = "AddBeadsInRow";

const zString ccAddBeadsInRow::GetType()
{
	return m_Type;
}

// We use an anonymous namespace to wrap the call to the factory object
// so that it is not accessible from outside this file. The identifying
// string for the command is stored in the m_Type static member variable.
//
// Note that the Create() function is not a member function of the
// command class but a global function hidden in the namespace.

namespace private_ns_ccAddBeadsInRow
{
	xxCommand* Create(long executionTime) {return new ccAddBeadsInRow(executionTime);}

	const zString id = ccAddBeadsInRow::GetType();

	const bool bRegistered = acfCommandFactory::Instance()->Register(id, Create);
}

//////////////////////////////////////////////////////////////////////
// Global functions
//////////////////////////////////////////////////////////////////////

bool operator==(const ccAddBeadsInRow& a, const ccAddBeadsInRow& b)
{
	if(dynamic_cast<const ccSelection&>(a) == dynamic_cast<const ccSelection&>(b) &&
		a.m_SliceIndex == b.m_SliceIndex &&
		a.m_RowIndex   == b.m_RowIndex)
		return true;
	else
		return false;
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ccAddBeadsInRow::ccAddBeadsInRow(long executionTime) : ccSelection(executionTime),
														 m_SliceIndex(0),
														 m_RowIndex(0)
{
}

ccAddBeadsInRow::ccAddBeadsInRow(const ccAddBeadsInRow& oldCommand) : ccSelection(oldCommand),
																	   m_SliceIndex(oldCommand.m_SliceIndex),
																	   m_RowIndex(oldCommand.m_RowIndex)
{
}


ccAddBeadsInRow::~ccAddBeadsInRow()
{
}

// Member functions to write/read the data specific to the command.
// Note that spaces are written here not in the ccSelection or xxCommand
// classes. Also note that we have to pass the command to put or get up 
// to the ccSelection base class so that its data can be read/written 
// before we get the data for this class.
//
// Arguments
// *********
//
//	m_SliceIndex	Identifier showing the slice
//	m_RowIndex		Identifier showing the location of a row in the slice
//

zOutStream& ccAddBeadsInRow::put(zOutStream& os) const
{
#if EnableXMLCommands == SimXMLEnabled

	// XML output - write the start tags first then write the base class
	// data before writing data in this class

	putXMLStartTags(os);
	ccSelection::put(os);
	os << "<SliceIndex>" << m_SliceIndex << "</SliceIndex>" << zEndl;
	os << "<RowIndex>"   << m_RowIndex   << "</RowIndex>"   << zEndl;
	putXMLEndTags(os);

#elif EnableXMLCommands == SimXMLDisabled

	// ASCII output 
	putASCIIStartTags(os);
	ccSelection::put(os);
	os << m_SliceIndex << " " << m_RowIndex;
	putASCIIEndTags(os);

#endif

	return os;
}

zInStream& ccAddBeadsInRow::get(zInStream& is)
{
	// Read the base class data first

	ccSelection::get(is);

	// Now read the location of the slice and the row within the slice. 
	// Note that the CNT cell structure, and hence the data read here, 
	// is zero-based so the slice and row indices can both be zero.
	// The normal read in by ccSelection cannot though.

	is >> m_SliceIndex >> m_RowIndex;

	if(!is.good())
		SetCommandValid(false);

	return is;
}


// Implementation of the command that is sent by the SimBox to each xxCommand
// object to see if it is the right time for it to carry out its operation.
// We return a boolean so that the SimBox can see if the command executed or not
// as this may be useful for considering several commands.
//

bool ccAddBeadsInRow::Execute(long simTime, ISimCmd* const pISimCmd) const
{
	if(simTime == GetExecutionTime())
	{
		pISimCmd->AddBeadsInRow(this);
		return true;
	}
	else
		return false;
}

const zString ccAddBeadsInRow::GetCommandType() const
{
	return m_Type;
}

// Function to return a pointer to a copy of the current command.

const xxCommand* ccAddBeadsInRow::GetCommand() const
{
	return new ccAddBeadsInRow(*this);
}


// Function to check that the data required by the AddBeadsInRow command
// is valid. We cannot check the generic data by a call to the base class's
// IsDataValid() function because we don't want to add target names to the
// list but to check that they already exist.
// 
// Check that the slice and row indices specified lie within the SimBox 
// in the direction specified by its normal.
//

bool ccAddBeadsInRow::IsDataValid(const CInputData &riData) const
{
	// Check target name exists and consists only of allowed characters.

	if( !riData.IsExternalNameValid(GetLabel()) )
		return ErrorTrace("Invalid Add force target name");
	else if( !riData.IsForceTargetPresent(GetLabel()) )
		return ErrorTrace("Add force target does not exist");

	// Check that the bead names exist in the CInputData's bead names map

	const StringSequence targetBeadNames = GetNames();

	for(cStringSequenceIterator iBeadName=targetBeadNames.begin(); iBeadName != targetBeadNames.end(); iBeadName++)
	{
		if(!riData.IsExternalNameValid(*iBeadName))
			return ErrorTrace("Invalid Add command bead name");
		else if(!riData.IsBeadinMap(*iBeadName))
			return ErrorTrace("Add command bead does not exist");
	}

	// Check the slice and row indices lie within the SimBox

	if(GetX() == 1 && (GetSliceIndex() < 0 || GetSliceIndex() >= riData.GetCNTXCellNo() ||
			           GetRowIndex()   < 0 || GetRowIndex()   >= riData.GetCNTZCellNo()) )
	{
		return ErrorTrace("Selection command Z Row outside SimBox");
	}
	if(GetY() == 1 && (GetSliceIndex() < 0 || GetSliceIndex() >= riData.GetCNTYCellNo() ||
			           GetRowIndex()   < 0 || GetRowIndex()   >= riData.GetCNTXCellNo()) )
	{
		return ErrorTrace("Selection command X Row outside SimBox");
	}
	if(GetZ() == 1 && (GetSliceIndex() < 0 || GetSliceIndex() >= riData.GetCNTZCellNo() ||
			           GetRowIndex()   < 0 || GetRowIndex()   >= riData.GetCNTYCellNo()) )
	{
		return ErrorTrace("Selection command Y Row outside SimBox");
	}

	// The following tests see if the command is a duplicate of a previous one.
	// Because there are several classes all inheriting from ccSelection this
	// is a non-trivial exercise.
	//
	// Find the Select command's target corresponding to this Add command. We only 
	// check commands that precede the current one as the Selectxxx command must
	// occur before any Addxxx command. We use the command id's to search.

	cCommandIterator iCmd = riData.GetCommands().begin();

	while((*iCmd)->GetId() < GetId())
	{
		// Find the Select command corresponding to the current Add command

		const ccSelectBeadsInRow* const pTarget = dynamic_cast<const ccSelectBeadsInRow*>(*iCmd);

		if(pTarget && pTarget->GetLabel() == GetLabel())
		{
			// Now we have the target check equality between the ccSelection 
			// part of the current Addxxx and the Selectxxx or Addxxx commands

			const ccSelection* const pS1 = dynamic_cast<const ccSelection*>(*iCmd);
			const ccSelection* const pS2 = dynamic_cast<const ccSelection*>(this);

			if(*pS1 == *pS2)
			{
				// Found the target now check the slice and row indices

				if(pTarget->GetSliceIndex() == GetSliceIndex() &&
					pTarget->GetRowIndex()  == GetRowIndex())
				{
					return ErrorTrace("Error: Add command duplicates Select/Add command");
				}
			}
			else if(pS1->GetX() == pS2->GetX() && pS1->GetY() == pS2->GetY() && 
				    pS1->GetZ() == pS2->GetZ() && 
					pTarget->GetSliceIndex() == GetSliceIndex() &&
					pTarget->GetRowIndex()   == GetRowIndex())
			{
				// If the normals and slice, row indices are the same check that command 
				// bead lists do not overlap. This means that no entry in the current 
				// Addxxx command can occur in the iCmd list.

				bool bDistinctNames = true;

				for(long unsigned int i=0; i<GetNames().size(); i++)
				{
					zString name = GetNames().at(i);

					const StringSequence targetBeads = pTarget->GetNames();

					if(std::find(targetBeads.begin(), targetBeads.end(), name) != targetBeads.end())
						bDistinctNames = false;
				}
				if(!bDistinctNames)
					return ErrorTrace("Error: Beads overlap with previous Select/Add command");;
			}
		}
		iCmd++;
	}

	return true;
}
