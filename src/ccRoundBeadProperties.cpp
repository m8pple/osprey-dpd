/* **********************************************************************
Copyright 2020  Dr. J. C. Shillcock and Prof. Dr. R. Lipowsky, Director at the Max Planck Institute (MPI) of Colloids and Interfaces; Head of Department Theory and Bio-Systems.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************** */
// ccRoundBeadProperties.cpp: implementation of the ccRoundBeadProperties class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "SimDefs.h"
#include "ccRoundBeadProperties.h"
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

long ccRoundBeadProperties::GetArgumentTotal()
{
	return 3;
}

const zString ccRoundBeadProperties::m_Type = "RoundBeadProperties";

const zString ccRoundBeadProperties::GetType()
{
	return m_Type;
}

// We use an anonymous namespace to wrap the call to the factory object
// so that it is not accessible from outside this file. The identifying
// string for the command is stored in the m_Type static member variable.
//
// Note that the Create() function is not a member function of the
// command class but a global function hidden in the namespace.

namespace private_ns_ccRoundBeadProperties
{
	xxCommand* Create(long executionTime) {return new ccRoundBeadProperties(executionTime);}

	const zString id = ccRoundBeadProperties::GetType();

	const bool bRegistered = acfCommandFactory::Instance()->Register(id, Create);
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ccRoundBeadProperties::ccRoundBeadProperties(long executionTime) : xxCommand(executionTime),
										   m_base(10),
										   m_posFracDigits(4),
										   m_vecMantDigits(4)
{
}

ccRoundBeadProperties::ccRoundBeadProperties(long executionTime, int base, int posFracDigits, int vecMantDigits) : xxCommand(executionTime),
											m_base(base),
											m_posFracDigits(posFracDigits),
											m_vecMantDigits(vecMantDigits)
{
}

ccRoundBeadProperties::ccRoundBeadProperties(const ccRoundBeadProperties& oldCommand) : xxCommand(oldCommand),
													m_base(oldCommand.m_base),
													m_posFracDigits(oldCommand.m_posFracDigits),
													m_vecMantDigits(oldCommand.m_vecMantDigits)
{
}


ccRoundBeadProperties::~ccRoundBeadProperties()
{

}

// Member functions to write/read the data specific to the command.
// Note that spaces are written here not in the xxCommand function.
//
// Arguments
// *********
//
// m_Path	- Where to write the file to
//

zOutStream& ccRoundBeadProperties::put(zOutStream& os) const
{
#if EnableXMLCommands == SimXMLEnabled

	// XML output
	putXMLStartTags(os);
	os << "<Base>" << m_base << "</Base>" << zEndl;
	os << "<PosFracDigits>" << m_posFracDigits << "</PosFracDigits>" << zEndl;
	os << "<VecMantDigits>" << m_base << "</VecMantDigits>" << zEndl;
	putXMLEndTags(os);

#elif EnableXMLCommands == SimXMLDisabled

	// ASCII output 
	putASCIIStartTags(os);
	os << m_base << " "<<m_posFracDigits<<" "<<m_vecMantDigits;
	putASCIIEndTags(os);

#endif

	return os;
}

zInStream& ccRoundBeadProperties::get(zInStream& is)
{
	is >> m_base;
	is >> m_posFracDigits;
	is >> m_vecMantDigits;

	if(!is.good())
		SetCommandValid(false);

	return is;
}


// Implementation of the command that is sent by the SimBox to each xxCommand
// object to see if it is the right time for it to carry out its operation.
// We return a boolean so that the SimBox can see if the command executed or not
// as this may be useful for considering several commands.


// https://stackoverflow.com/a/29583280
static double frexp10(double arg, int * exp)
{
   *exp = (arg == 0) ? 0 : 1 + (int)std::floor(std::log10(std::fabs(arg) ) );
   return arg * std::pow(10 , -(*exp));    
}


bool ccRoundBeadProperties::Execute(long simTime, ISimCmd* const pISimCmd) const
{
	if(simTime != GetExecutionTime()){
		return false;
	}

	CSimBox *box=dynamic_cast<CSimBox*>(pISimCmd);
	if(!box){
		fprintf(stderr, "Export only works with basic sim box.");
		exit(1);
	}

	auto beads=box->GetSimBoxBeads();

	double pos_scale_dec=0, pos_scale_dec_inv=0;
	int pos_digits_bin=m_posFracDigits;
	double mant_scale_dec=0, mant_scale_dec_inv=0;
	int mant_digits_bin=m_vecMantDigits;

	if(m_base==10){
		pos_scale_dec=pow(10, m_posFracDigits);
		pos_scale_dec_inv=1.0/pos_scale_dec;
		mant_scale_dec=pow(10, m_vecMantDigits);
		mant_scale_dec_inv=1.0/mant_scale_dec;
	}

	auto round_pos = [&](double x) -> double 
	{
		if(m_base==2){
			x=ldexp( floor(ldexp( x, pos_digits_bin )), -pos_digits_bin );
		}else if(m_base==10){
			x=floor( x * pos_scale_dec ) * pos_scale_dec_inv;
		}else{
			assert(0);
		}
		return x;
	};

	auto round_mant = [&](double x) -> double
	{
		if(x!=0){
			int exponent;
			if(m_base==2){
				double f=frexp(x, &exponent);
				f=ldexp( round( ldexp(f, mant_digits_bin) ), -mant_digits_bin);
				x=ldexp(f, exponent);
			}else{
				// TODO : This is extremely slow.
				double f=frexp10(x, &exponent);
				f=round( f * mant_scale_dec ) * mant_scale_dec_inv;
				x=std::pow(10, exponent) * f;
			}
			
		}
		return x;
	};

	for(CAbstractBead *b : beads)
	{
		b->SetXPos( round_pos(b->GetXPos()) );
		b->SetYPos( round_pos(b->GetYPos()) );
		b->SetZPos( round_pos(b->GetZPos()) );

		b->SetXMom( round_mant(b->GetXMom()) );
		b->SetYMom( round_mant(b->GetYMom()) );
		b->SetZMom( round_mant(b->GetZMom()) );

		b->SetXForce(round_mant(b->GetXForce()) );
		b->SetYForce( round_mant(b->GetYForce()) );
		b->SetZForce( round_mant(b->GetZForce()) );
	}

	return true;
}

const zString ccRoundBeadProperties::GetCommandType() const
{
	return m_Type;
}

// Function to return a pointer to a copy of the current command.

const xxCommand* ccRoundBeadProperties::GetCommand() const
{
	return new ccRoundBeadProperties(*this);
}

// Function to check that the command is valid.
// TODO

bool ccRoundBeadProperties::IsDataValid(const CInputData &riData) const
{
	if(m_base!=2 && m_base!=10){
		return false;
	}
	if(m_posFracDigits<2) return false;
	if(m_vecMantDigits<2) return false;

	return true;
}

