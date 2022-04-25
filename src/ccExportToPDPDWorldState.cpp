/* **********************************************************************
Copyright 2020  Dr. J. C. Shillcock and Prof. Dr. R. Lipowsky, Director at the Max Planck Institute (MPI) of Colloids and Interfaces; Head of Department Theory and Bio-Systems.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************** */
// ccExportToPDPDWorldState.cpp: implementation of the ccExportToPDPDWorldState class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "SimDefs.h"
#include "ccExportToPDPDWorldState.h"
#include "ISimCmd.h"
#include "InputData.h"

#include "AbstractBead.h"
#include "Bead.h"
#include "SimBox.h"
#include "ISimBox.h"
#include "Polymer.h"
#include "Bond.h"
#include "BondPair.h"

// TODO - make this portable. Turn it off for non-GNU
#include <ext/stdio_filebuf.h>

#include <cmath>

#include <unordered_map>

//////////////////////////////////////////////////////////////////////
// Global members
//////////////////////////////////////////////////////////////////////

// Static member variable containing the identifier for this command. 
// The static member function GetType() is invoked by the xxCommandObject 
// to compare the type read from the control data file with each
// xxCommand-derived class so that it can create the appropriate object 
// to hold the command data.

long ccExportToPDPDWorldState::GetArgumentTotal()
{
	return 1; // Just the path
}

const zString ccExportToPDPDWorldState::m_Type = "ExportToPDPDWorldState";

const zString ccExportToPDPDWorldState::GetType()
{
	return m_Type;
}

// We use an anonymous namespace to wrap the call to the factory object
// so that it is not accessible from outside this file. The identifying
// string for the command is stored in the m_Type static member variable.
//
// Note that the Create() function is not a member function of the
// command class but a global function hidden in the namespace.

namespace private_ns_ccExportToPDPDWorldState
{
	xxCommand* Create(long executionTime) {return new ccExportToPDPDWorldState(executionTime);}

	const zString id = ccExportToPDPDWorldState::GetType();

	const bool bRegistered = acfCommandFactory::Instance()->Register(id, Create);
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ccExportToPDPDWorldState::ccExportToPDPDWorldState(long executionTime) : xxCommand(executionTime),
										   m_Path("pdpd-world-state.txt")
{
}

ccExportToPDPDWorldState::ccExportToPDPDWorldState(long executionTime, zString path) : xxCommand(executionTime),
											m_Path(path)
{
}

ccExportToPDPDWorldState::ccExportToPDPDWorldState(const ccExportToPDPDWorldState& oldCommand) : xxCommand(oldCommand),
													m_Path(oldCommand.m_Path)
{
}


ccExportToPDPDWorldState::~ccExportToPDPDWorldState()
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

zOutStream& ccExportToPDPDWorldState::put(zOutStream& os) const
{
#if EnableXMLCommands == SimXMLEnabled

	// XML output
	putXMLStartTags(os);
	os << "<ExportFilePath>" << m_Path << "</ExportFilePath>" << zEndl;
	putXMLEndTags(os);

#elif EnableXMLCommands == SimXMLDisabled

	// ASCII output 
	putASCIIStartTags(os);
	os << m_Path;
	putASCIIEndTags(os);

#endif

	return os;
}

zInStream& ccExportToPDPDWorldState::get(zInStream& is)
{
	is >> m_Path;

	if(!is.good())
		SetCommandValid(false);

	return is;
}


// Implementation of the command that is sent by the SimBox to each xxCommand
// object to see if it is the right time for it to carry out its operation.
// We return a boolean so that the SimBox can see if the command executed or not
// as this may be useful for considering several commands.


struct polymer_template_t
{
	struct bead_info_t
	{
		unsigned offset;
		unsigned type;
	};

	struct bond_info_t
	{
		unsigned offset;
		unsigned head_offset;
		unsigned tail_offset;
		double kappa;
		double r0;
	};

	struct bond_pair_info_t
	{
		unsigned offset;
		unsigned head_offset;
		unsigned tail_offset;
		double kappa;
		double theta0;
	};


	unsigned id;
	std::string name;
	std::vector<bead_info_t> beads;
	std::vector<bond_info_t> bonds;
	std::vector<bond_pair_info_t> bond_pairs;
};

polymer_template_t extract_polymer_template(const CSimBox *box, const CPolymer *p)
{
	polymer_template_t res;
	const auto &beads=p->GetBeads();
	const auto &bonds=p->GetBonds();
	const auto &bond_pairs=p->GetBondPairs();

	res.id=p->GetType();
	res.name=box->GetPolymerNameFromType(p->GetType());

	std::unordered_map<const CAbstractBead*,unsigned> bead_ptr_to_offset;
	for(const auto *b : beads){
		typename polymer_template_t::bead_info_t t;
		t.offset=res.beads.size();
		t.type=b->GetType();

		res.beads.push_back(t);
		auto it=bead_ptr_to_offset.insert({b, t.offset});
		if(!it.second){
			fprintf(stderr, "Corrupt polymer : bead appears twice.");
			exit(1);
		}
	}

	std::unordered_map<const CBond*,unsigned> bond_ptr_to_offset;
	for(const auto *b : bonds){
		typename polymer_template_t::bond_info_t t;

		t.offset=res.bonds.size();
		t.head_offset=bead_ptr_to_offset.at(b->GetHead());
		t.tail_offset=bead_ptr_to_offset.at(b->GetTail());
		t.r0=b->GetUnStrLength();
		t.kappa=b->GetSprConst();
		res.bonds.push_back(t);
		
		auto it=bond_ptr_to_offset.insert({b, t.offset});
		if(!it.second){
			fprintf(stderr, "Corrupt polymer : bond appears twice.");
			exit(1);
		}
	}


	std::unordered_map<const CBondPair*,unsigned> bond_pair_ptr_to_offset;
	for(const auto *bp : bond_pairs){
		typename polymer_template_t::bond_pair_info_t t;

		t.offset=res.bond_pairs.size();
		t.head_offset=bond_ptr_to_offset.at(bp->GetFirst());
		t.tail_offset=bond_ptr_to_offset.at(bp->GetSecond());
		t.theta0=bp->GetPhi0();
		t.kappa=bp->GetModulus();
		res.bond_pairs.push_back(t);

		auto it=bond_pair_ptr_to_offset.insert({bp, t.offset});
		if(!it.second){
			fprintf(stderr, "Corrupt polymer : bond pair appears twice.");
			exit(1);
		}
	}

	return res;
}

bool ccExportToPDPDWorldState::Execute(long simTime, ISimCmd* const pISimCmd) const
{
	bool res;
	if(m_Path.size()>3 && m_Path.substr(m_Path.size()-3)==".gz"){
		std::string cmd="gzip -9 -c > "+m_Path;
        FILE *f=popen(cmd.c_str(), "w");
        if(!f){
            throw std::runtime_error("Error when spawning gunzip command '"+cmd+"'");
        }
        
        __gnu_cxx::stdio_filebuf<char> buf(f, std::ios_base::out);
        std::ostream fs(&buf);

        res=ExecuteImpl(simTime, pISimCmd, fs);
		fs.flush();

        pclose(f);
	}else{
		std::ofstream dst(m_Path.c_str());
		if(!dst.is_open()){
			fprintf(stderr, "ccExportToPDPDWorldState::Execute - Couldn't open destination file '%s'", m_Path.c_str());
			exit(1);
		}
		res=ExecuteImpl(simTime, pISimCmd, dst);
	}
	return res;
}

bool ccExportToPDPDWorldState::ExecuteImpl(long simTime, ISimCmd* const pISimCmd, std::ostream &dst) const
{
	if(simTime != GetExecutionTime()){
		return false;
	}

	CSimBox *box=dynamic_cast<CSimBox*>(pISimCmd);
	if(!box){
		fprintf(stderr, "Export only works with basic sim box.");
		exit(1);
	}

	fprintf(stderr, "ccExportToPDPDWorldState to %s\n", m_Path.c_str());

	auto bead_types=box->GetBeadTypes();
	auto polymers=box->GetAllPolymers();
	auto polymer_types=box->GetPolymerTypes();
	auto beads=box->GetSimBoxBeads();

	

	// TODO: Make this selectable?
	// We print in very high precision. This is to avoid quantisation due to printing
	// having a big effect (it is very noticeable for the default 6 figures). 
	// If you want to reduce precision, use ccRoundBeadProperties to round
	// the velocities and forces off first.
	dst.precision(14);

	dst<<"WorldState v0 "<<bead_types.size()<<" "<<polymer_types.size()<<" "<<beads.size()<<" "<<polymers.size()<<"\n";

	dst<<"T "<<box->GetCurrentTime()<<" "<<box->GetTimeStepSize()<<"\n";
	dst<<"Lambda "<<box->GetISimBox()->GetLambda()<<"\n";
    dst<<"Origin "<<box->GetSimBoxXOrigin()<<" "<<box->GetSimBoxYOrigin()<<" "<<box->GetSimBoxZOrigin()<<"\n";
	dst<<"Box "<<box->GetSimBoxXOrigin()+box->GetSimSpaceXLength()<<" "<<box->GetSimBoxYOrigin()+box->GetSimSpaceYLength()<<" "<<box->GetSimBoxZOrigin()+box->GetSimSpaceZLength()<<"\n";
	dst<<"Seed "<<std::abs(pISimCmd->GetISimBox()->GetRNGSeed())<<"\n"; // TODO
	for(int i=0; i<bead_types.size(); i++){
		dst<<"ConservativeStrength "<<i;
		for(int j=0; j<bead_types.size(); j++){
			dst<<" "<<box->GetISimBox()->GetDPDBeadConsInt(i, j);
		}
		dst<<"\n";
	}
	for(int i=0; i<bead_types.size(); i++){
		dst<<"DissipativeStrength "<<i;
		for(int j=0; j<bead_types.size(); j++){
			dst<<" "<<box->GetISimBox()->GetDPDBeadDissInt(i, j);
		}
		dst<<"\n";
	}

	for(int i=0; i<bead_types.size(); i++){
		const auto &bt=bead_types[i];
		// The id always seems to be -1 for all bead types?
		/*if(bt->GetId()!=i){
			fprintf(stderr, "Bead type order is incorrect : type %u at index %u.\n", bt->GetId(), i);
			//exit(1);
		}*/
		dst<<"BeadType "<<i<<" "<<box->GetBeadNameFromType(i)<<" "<<bt->GetRadius()<<"\n";
	}

	std::vector<polymer_template_t> templates;
	for(const CPolymer *p : polymer_types){
		templates.push_back(extract_polymer_template(box, p));
		const auto &pt = templates.back();
		dst<<"PolymerType "<<pt.id<<" "<<pt.name<<" "<<pt.beads.size()<<" "<<pt.bonds.size()<<" "<<pt.bond_pairs.size()<<"\n";
		dst<<"  BeadTypeIndices ";
		for(const auto &b : pt.beads){
			dst<<" "<<b.type;
		}
		dst<<"\n";
		for(const auto &b : pt.bonds){
			dst<<"  Bond "<<b.head_offset<<" "<<b.tail_offset<<" "<<b.kappa<<" "<<b.r0<<"\n";
		}
		for(const auto &bp : pt.bond_pairs){
			dst<<"  BondPair "<<bp.tail_offset<<" "<<bp.head_offset<<" "<<bp.kappa<<" "<<bp.theta0<<"\n";
		}
		dst<<"\n";
	}

	auto require = [](bool cond, const char *msg)
	{
		if(!cond){
			fprintf(stderr, "ccExportToPDPDWorldState::Execute - polymer instance does not match template : %s.", msg);
			exit(1);
		}
	};

	std::unordered_map<const CAbstractBead*,unsigned> bead_ptr_to_offset;
	std::unordered_map<const CBond*,unsigned> bond_ptr_to_offset;
	for(const CPolymer *p : polymers){
		// Grab the map
		const auto &pt = templates.at(p->GetType());

		bead_ptr_to_offset.clear();
		unsigned off=0;
		for(const auto *b : p->GetBeads()){
			dst<<"B "<<b->GetId()-1<<" "<<p->GetId()-1<<" "<<p->GetType()<<" "<<off;
			dst<<" "<<b->GetXPos()<<" "<<b->GetYPos()<<" "<<b->GetZPos();
			dst<<" "<<b->GetXMom()<<" "<<b->GetYMom()<<" "<<b->GetZMom();
			dst<<" "<<b->GetXForce()<<" "<<b->GetYForce()<<" "<<b->GetZForce();
			dst<<"\n";

			require(b->GetType() == pt.beads[off].type, "Bead type does not match bead type in template");
			require(bead_ptr_to_offset.insert({b,off}).second, "Bead exists twice in polymer instance.");

			++off;
		}

		bond_ptr_to_offset.clear();
		off=0;
		for(const CBond *b : p->GetBonds()){
			require(b->GetUnStrLength()==pt.bonds[off].r0, "Bond length does not match.");
			require(b->GetSprConst()==pt.bonds[off].kappa, "Bond constant does not match");
			require( bead_ptr_to_offset.at(b->GetHead()) == pt.bonds[off].head_offset, "Head of bond." );
			require( bead_ptr_to_offset.at(b->GetTail()) == pt.bonds[off].tail_offset, "Tail of bond." );

			require( bond_ptr_to_offset.insert({b,off}).second, "Bond appears twice.");

			++off;
		}

		off=0;
		for(const CBondPair *bp : p->GetBondPairs()){
			require(bp->GetPhi0() == pt.bond_pairs[off].theta0, "BondPair angle.");
			require(bp->GetModulus() == pt.bond_pairs[off].kappa, "BondPair constant.");
			require( bond_ptr_to_offset.at(bp->GetFirst()) == pt.bond_pairs[off].head_offset, "BondPair head.");
			require( bond_ptr_to_offset.at(bp->GetSecond()) == pt.bond_pairs[off].tail_offset, "BondPair tail.");
			++off;
		}
	}

	return true;
}

const zString ccExportToPDPDWorldState::GetCommandType() const
{
	return m_Type;
}

// Function to return a pointer to a copy of the current command.

const xxCommand* ccExportToPDPDWorldState::GetCommand() const
{
	return new ccExportToPDPDWorldState(*this);
}

// Function to check that the command is valid.
// TODO

bool ccExportToPDPDWorldState::IsDataValid(const CInputData &riData) const
{
	return true;
}

