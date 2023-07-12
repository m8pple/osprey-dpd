#include <vector>
#include <atomic>
#include <vector>
#include <mutex>
#include "DebugAssert.h"
#include <cstdlib>
#include <set>
#include <cstring>
#include <unordered_map>

#include "StdAfx.h"
#include "SimDefs.h"
#include "xxBase.h"
#include "AbstractBead.h"
#include "CNTCell.h"

#include "IIntegrationEngine.h"

#include "SimBox.h"
#include "mpsSimBox.h"
#include "aeActiveSimBox.h"

#include "ISimBox.h"

#include "StateLogger.h"

#include "ParallelContext.h"

class SimEngineFast
    : public IIntegrationEngine
{
public:
    SimEngineFast()
    {}

    std::string Name() const override
    {
        return "SimEngineFast";
    }

    virtual bool IsParallel() const override
    { return false; }

    bool IsProductionReady() const override
    { return true; }

    int GetEstimatedMerit() const override
    { return 10; }  

    run_result Run(ISimBox *box, bool modified, unsigned start_sim_time, unsigned num_steps) override
    {
        CSimBox *mbox=const_cast<CSimBox*>(box->GetSimBox());

        if(modified){
            auto tmp=CanSupport(box);
            if(tmp.status!=Supported){
                return {tmp.status, tmp.reason, 0};
            }
        }
     
        ImportGlobals(mbox);
     
        if(0){
            EvolveFast(mbox, start_sim_time, num_steps);
        }else{
            for(unsigned i=0; i<num_steps; i++){
                EvolveFast(mbox, start_sim_time+i, 1);
            }
        }

        return {Supported, {}, num_steps};
    }

private:
    void ImportGlobals(CSimBox *mbox)
    {
        DEBUG_ASSERT(CanSupport(mbox->GetISimBox()).status==Supported);

        m_dt=CCNTCell::m_dt;
        m_halfdt=CCNTCell::m_halfdt;
        m_halfdt2=CCNTCell::m_halfdt2;
        m_SimBoxXLength=CCNTCell::m_SimBoxXLength;
        m_SimBoxYLength=CCNTCell::m_SimBoxYLength;
        m_SimBoxZLength=CCNTCell::m_SimBoxZLength;
        m_CNTXCellNo=CCNTCell::m_CNTXCellNo;
        m_CNTYCellNo=CCNTCell::m_CNTYCellNo;
        m_CNTZCellNo=CCNTCell::m_CNTZCellNo;
    }

    double m_dt;
    double m_halfdt;
    double m_halfdt2;
    double m_SimBoxXLength;
    double m_SimBoxYLength;
    double m_SimBoxZLength;
    unsigned m_CNTXCellNo, m_CNTYCellNo, m_CNTZCellNo;

    void PrefetchHint(const CCNTCell &cell)
	{
		for(unsigned i=0; i<cell.m_lBeads.size(); i++){
			__builtin_prefetch(&cell.m_lBeads[i]->m_Type);
		}
	}

    //#pragma GCC push_options
    //#pragma GCC optimize("fast-math")
    void UpdateForceFast(CCNTCell *cell)
    {
        // Forwarding to CCNTCell members
        const bool m_bExternal=cell->m_bExternal;
        auto &m_lBeads=cell->m_lBeads;
        auto &m_aIntNNCells=cell->m_aIntNNCells;

        /////////////////////////////////////////////////////////////
        // Begin the adapted body from CCNTCell::UpdateForce

        mpsSimBox::GlobalCellCounter++;  // increment the counter for intra-cell force calculations
        
        long  localCellCellCounter = 0;

        if(cell->m_lBeads.empty()){
            localCellCellCounter = 13;
            mpsSimBox::GlobalCellCellIntCounter += localCellCellCounter;
            return;
        }

        PrefetchHint( *cell->m_aIntNNCells[0] );

        double rng_scale =  CCNTCell::m_invrootdt;

        // DPD and MD equations of motion

        AbstractBeadVectorIterator iterBead1;
        AbstractBeadVectorIterator iterBead2;

        double dx[3], dv[3], dx_dv[3], newForce[3];
        double dr, dr2;

        double gammap, rdotv, wr, wr2;
        double conForce, dissForce, randForce, normTotalForce;

        const double min_r2 = 0.000000001 * 0.000000001;

        /*
        for( iterBead1=m_lBeads.begin(); iterBead1!=m_lBeads.end(); iterBead1++ )
        {
            // First add interactions between beads in the current cell. Note that
            // we don't have to check the PBCs here and we perform a reverse loop
            // over the neighbouring beads until the iterators are equal. Because you can't
            // compare a forward and reverse iterator we compare the bead ids for
            // the terminating condition.
            for( riterBead2=m_lBeads.rbegin(); (*riterBead2)->m_id!=(*iterBead1)->m_id; ++riterBead2 )
            {
        */

        // It's better to avoid the data-dependent loop termination as it leads to
        // more branch mis-predicts. Better to loop based on size.
        int numLocal=m_lBeads.size();
        iterBead1=m_lBeads.begin();
        for(int ii1=0; ii1 < numLocal-1; ii1++, iterBead1++){
            iterBead2 = std::next(iterBead1);
            for(int ii2=ii1+1; ii2 < m_lBeads.size(); ii2++, iterBead2++){

                double dx2[3];
                for(int d=0; d<3; d++){
                    dx[d] = ((*iterBead1)->m_Pos[d] - (*iterBead2)->m_Pos[d]);
                    dx2[d] = dx[d] * dx[d];
                }

                dr2 = dx2[0] + dx2[1] + dx2[2];
                
                if(StateLogger::IsEnabled()){
                    int id1=(*iterBead1)->GetId()-1, id2=(*iterBead2)->GetId()-1;
                    StateLogger::LogBeadPairRefl("dpd_dx", id1, id2, dx);
                    StateLogger::LogBeadPairRefl("dpd_dr2", id1, id2, dr2);
                }
        
        

    // Calculate the interactions between the two beads for each simulation type.
    // For the DPD interactions we use the flag UseDPDBeadRadii to determine whether
    // the interaction radius is bead-specific or not. Note that when we use the
    // interaction radius we have to compare the actual distance between beads
    // not its square!

                if( dr2 < 1.0 && dr2 > min_r2 )
                {		
                    for(int d=0; d<3; d++){
                        dv[d] = ((*iterBead1)->m_Mom[d] - (*iterBead2)->m_Mom[d]);
                        dx_dv[d] = dx[d] * dv[d];
                    }

                    dr = sqrt(dr2);
                    wr = (1.0 - dr);
                    wr2 = wr*wr;
                    double inv_dr = 1.0/dr;

    // Conservative force magnitude

                    conForce  = CCNTCell::m_vvConsInt[(*iterBead1)->GetType()][(*iterBead2)->GetType()]*wr;

    // Dissipative and random force magnitudes. Note dr factor in newForce calculation

                    rdotv		= (dx_dv[0] + dx_dv[1] + dx_dv[2]) * inv_dr;
                    gammap		= CCNTCell::m_vvDissInt[(*iterBead1)->GetType()][(*iterBead2)->GetType()]*wr2;

                    dissForce	= -gammap*rdotv;			
                    double randNum    = CCNTCell::RandUniformBetweenBeads(*iterBead1, *iterBead2);	
                    randForce	= sqrt(gammap) * rng_scale * randNum;
                    normTotalForce = (conForce + dissForce + randForce) * inv_dr;

                    if(StateLogger::IsEnabled()){
						int id1=(*iterBead1)->GetId()-1, id2=(*iterBead2)->GetId()-1;
						StateLogger::LogBeadPairRefl("dpd_conForce", id1, id2, conForce);
						StateLogger::LogBeadPairRefl("dpd_randNum", id1, id2, randNum);
                        StateLogger::LogBeadPairRefl("dpd_randForce", id1, id2, randForce);
                        StateLogger::LogBeadPairRefl("dpd_dissForce", id1, id2, dissForce);
						StateLogger::LogBeadPairRefl("dpd_newForce", id1, id2, newForce);
					}

                    for(int d=0; d<3; d++){
                        newForce[d] = normTotalForce * dx[d];
                        (*iterBead1)->m_Force[d] += newForce[d];
                        (*iterBead2)->m_Force[d] -= newForce[d];
                    }
                }
            }
        }

        // This is pulled in from SetMom
        // We do it here as (hopefully) the beads will now all
        // be sat in the cache.
        for(unsigned i=0; i<m_lBeads.size(); i++){
            m_lBeads[i]->SetMovable();
        }

        // Next add in interactions with beads in neighbouring cells taking the
        // PBCs into account and the presence of a wall. The PBCs are only applied
        // if both the current CNT cell and the neighbouring one are external.

        for( int i=0; i<13; i++ )
        {	
            if(i<12){
                PrefetchHint(*m_aIntNNCells[i+1]);
            }

            bool both_external = m_bExternal && m_aIntNNCells[i]->IsExternal();

            iterBead2=m_aIntNNCells[i]->m_lBeads.begin();
            for(int jj1=0; jj1<m_aIntNNCells[i]->m_lBeads.size(); jj1++, iterBead2++){

                iterBead1=m_lBeads.begin();
                for(int jj2=0; jj2<numLocal; jj2++, iterBead1++){
                    localCellCellCounter++;  // Increment the local cell-cell inteaction counter

                    for(int d=0; d<3; d++){
                        dx[d] = ((*iterBead1)->m_Pos[d] - (*iterBead2)->m_Pos[d]);
                    }

                    if( both_external )
                    {
                        if( dx[0] > CCNTCell::m_HalfSimBoxXLength )
                            dx[0] = dx[0] - CCNTCell::m_SimBoxXLength;
                        else if( dx[0] < -CCNTCell::m_HalfSimBoxXLength )
                            dx[0] = dx[0] + CCNTCell::m_SimBoxXLength;

                        if( dx[1] > CCNTCell::m_HalfSimBoxYLength )
                            dx[1] = dx[1] - CCNTCell::m_SimBoxYLength;
                        else if( dx[1] < -CCNTCell::m_HalfSimBoxYLength )
                            dx[1] = dx[1] + CCNTCell::m_SimBoxYLength;

                        if( dx[2] > CCNTCell::m_HalfSimBoxZLength )
                            dx[2] = dx[2] - CCNTCell::m_SimBoxZLength;
                        else if( dx[2] < -CCNTCell::m_HalfSimBoxZLength )
                            dx[2] = dx[2] + CCNTCell::m_SimBoxZLength;

                    }

                    dr2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

                    if(StateLogger::IsEnabled()){
                        int id1=(*iterBead1)->GetId()-1, id2=(*iterBead2)->GetId()-1;
                        StateLogger::LogBeadPairRefl("dpd_dx", id1, id2, dx);
                        StateLogger::LogBeadPairRefl("dpd_dr2", id1, id2, dr2);
                    }

                    if( dr2 < 1.0 && dr2 > min_r2)
                    {		
                        for(int d=0; d<3; d++){
                            dv[d] = ((*iterBead1)->m_Mom[d] - (*iterBead2)->m_Mom[d]);
                        }

                        dr = sqrt(dr2);
                    
                        wr = (1.0 - dr);
                        wr2 = wr*wr;
                        double inv_dr = 1.0/dr;

                        conForce	= CCNTCell::m_vvConsInt[(*iterBead1)->GetType()][(*iterBead2)->GetType()]*wr;

                        rdotv		= (dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2]) * inv_dr;
                        gammap		= CCNTCell::m_vvDissInt[(*iterBead1)->GetType()][(*iterBead2)->GetType()]*wr2;

                        dissForce	= -gammap*rdotv;		
                        double randNum = CCNTCell::RandUniformBetweenBeads(*iterBead1, *iterBead2);
                        randForce	= sqrt(gammap) * rng_scale * randNum;
                        normTotalForce = (conForce + dissForce + randForce) * inv_dr;

                        for(int d=0; d<3; d++){
                            newForce[d] = normTotalForce * dx[d];
                            (*iterBead1)->m_Force[d] += newForce[d];
                            (*iterBead2)->m_Force[d] -= newForce[d];
                        }

                        if(StateLogger::IsEnabled()){
                            int id1=(*iterBead1)->GetId()-1, id2=(*iterBead2)->GetId()-1;
                            StateLogger::LogBeadPairRefl("dpd_conForce", id1, id2, conForce);
                            StateLogger::LogBeadPairRefl("dpd_randNum", id1, id2, randNum);
                            StateLogger::LogBeadPairRefl("dpd_randForce", id1, id2, randForce);
                            StateLogger::LogBeadPairRefl("dpd_dissForce", id1, id2, dissForce);
                            StateLogger::LogBeadPairRefl("dpd_newForce", id1, id2, newForce);
                        }
                    }
                }
            }
        }

    // Divide the local cell-cell counter by the number of beads in this cell and add the result to the global cell-cell counter
        
        DEBUG_ASSERT(m_lBeads.size() > 0);
        localCellCellCounter /= m_lBeads.size();
        
        mpsSimBox::GlobalCellCellIntCounter += localCellCellCounter;
    }
    //#pragma GCC pop_options

    void UpdateMomThenPosFastV2(CSimBox *mbox, CCNTCell *cell)
    {
        // Forwarding to CCNTCell members
        const bool m_bExternal=cell->m_bExternal;
        auto &m_lBeads=cell->m_lBeads;
        auto &m_aIntNNCells=cell->m_aIntNNCells;

        /////////////////////////////////////////////////////////////
        // Begin the adapted body from CCNTCell::UpdateForce


        #ifndef NDEBUG
        if(CCNTCell::m_lambda!=0.5){
            cell->ErrorTrace("Attempt to call UpdatePosFast with lambda!=0.5");
            exit(1);
        }
        #endif

        double dx[3];

        int index=m_lBeads.size()-1;

        while(index >= 0){
            DEBUG_ASSERT(!m_lBeads.empty() && index < m_lBeads.size() );

            CAbstractBead *bead=m_lBeads[index];

            // Only allow bead to move if its IsMovable flag is true. This allows
            // us to indicate when a bead has already crossed a cell boundary and
            // should not be moved again in this timestep.

            if(!bead->GetMovable()){
                index--;
                continue;
            }
            
            /*for(int d=0; d<3; d++){
                DEBUG_ASSERT( m_BLCoord[d] <= bead->m_Pos[d] );
                DEBUG_ASSERT( m_TRCoord[d] >= bead->m_Pos[d] );
            }*/

            bead->SetNotMovable();

            // DPD and MD simulations
            // Store current values of position, velocity and force for later use

            // We do not write to m_oldPos, m_oldMom, or m_oldForce, as it is not needed in fast path (causes extra memory traffic)
    #ifndef NDEBUG
            bead->m_oldPos[0] = nanf("");
            bead->m_oldMom[0] = nanf("");
            bead->m_oldForce[0] = nanf("");
    #endif		

            // Apply the mom from the end of previous step (loop skewed)
            bead->m_Mom[0] = bead->m_Mom[0] + m_halfdt * bead->m_Force[0] ;
            bead->m_Mom[1] = bead->m_Mom[1] + m_halfdt* bead->m_Force[1];
            bead->m_Mom[2] = bead->m_Mom[2] + m_halfdt * bead->m_Force[2];

            // Update position coordinates

            dx[0] = m_dt*bead->m_Mom[0] + m_halfdt2*bead->m_Force[0];
            dx[1] = m_dt*bead->m_Mom[1] + m_halfdt2*bead->m_Force[1];
            dx[2] = m_dt*bead->m_Mom[2] + m_halfdt2*bead->m_Force[2];

            bead->m_Pos[0] += dx[0];
            bead->m_Pos[1] += dx[1];
            bead->m_Pos[2] += dx[2];

            // Update the unPBC coordinates for use in calculating bond lengths
            // where we don't want to have to check for beads at opposite side
            // of the simulation box.

            bead->m_unPBCPos[0] += dx[0];
            bead->m_unPBCPos[1] += dx[1];
            bead->m_unPBCPos[2] += dx[2];

            // We do not update m_dPos
    #ifndef NDEBUG
            bead->m_dPos[0] = nanf("");
    #endif

            // Update intermediate velocity

            // We know that lambda=0.5 for fast path
            DEBUG_ASSERT(CCNTCell::m_lamdt == CCNTCell::m_halfdt);
            bead->m_Mom[0] = bead->m_Mom[0] + m_halfdt*bead->m_Force[0];
            bead->m_Mom[1] = bead->m_Mom[1] + m_halfdt*bead->m_Force[1];
            bead->m_Mom[2] = bead->m_Mom[2] + m_halfdt*bead->m_Force[2];

            // Zero current force on beads so that UpdateForce() just has to form
            // a sum of all bead-bead interactions

            // Force counting is not tracked in fast mode
    #ifndef NDEBUG
            bead->m_ForceCounter = -1000;
    #endif

            bead->m_Force[0] = 0.0;
            bead->m_Force[1] = 0.0;
            bead->m_Force[2] = 0.0;

            /* 
            Idea here is to make it a branchless as possible. Expected case is that beads
            mostly stay in the same cell, so all boundary conditions need to be evaluated.
            */
            bool moved=0;
            for(int d=0; d<3; d++){
                moved |= (bead->m_Pos[d] > cell->m_TRCoord[d]) | (bead->m_Pos[d] < cell->m_BLCoord[d]);
            }
            if(!moved){
                for(int d=0; d<3; d++){
                    DEBUG_ASSERT( cell->m_BLCoord[d] <= bead->m_Pos[d] );
                    DEBUG_ASSERT( cell->m_TRCoord[d] >= bead->m_Pos[d] );
                }

                // If the bead did not change cells increment the
                // iterator by hand.
                index--;
                continue;
            }

            // Enforce wrapping. Unconditional on all dims, as movement is relatively rare
            // Use while loop as it is safer for super-fast beads and small boxes, and costs little (hopefully)
            while(bead->m_Pos[0] < 0)	bead->m_Pos[0] += m_SimBoxXLength;
            while(bead->m_Pos[1] < 0)	bead->m_Pos[1] += m_SimBoxYLength;
            while(bead->m_Pos[2] < 0)	bead->m_Pos[2] += m_SimBoxZLength;
            while(bead->m_Pos[0] >= m_SimBoxXLength)	bead->m_Pos[0] -= m_SimBoxXLength;
            while(bead->m_Pos[1] >= m_SimBoxYLength)	bead->m_Pos[1] -= m_SimBoxYLength;
            while(bead->m_Pos[2] >= m_SimBoxZLength)	bead->m_Pos[2] -= m_SimBoxZLength;

            int ix=floor(bead->m_Pos[0]), iy=floor(bead->m_Pos[1]), iz=floor(bead->m_Pos[2]);
            int cell_index = m_CNTXCellNo*(m_CNTYCellNo*iz+iy) + ix;

            m_lBeads.erase( m_lBeads.begin() + index );

            mbox->AddBeadToCNTCell(cell_index, bead);
            index--;
        }
    }

    void UpdateMomFastReverse(CSimBox *mbox, CCNTCell *cell)
    {   
        // Forwarding to CCNTCell members
        auto &m_lBeads=cell->m_lBeads;

        /////////////////////////////////////////////////////////////
        // Begin the adapted body from CCNTCell::UpdateMom

        for( AbstractBeadVectorIterator iterBead=m_lBeads.begin(); iterBead!=m_lBeads.end(); iterBead++ )
        {
            if((*iterBead)->SetMovable())	// flag ignored by immovable beads
            {
                // We have already updated m_Mom to m_oldMom + m_halfdt * m_oldForce
                (*iterBead)->m_Mom[0] = (*iterBead)->m_Mom[0] -= m_halfdt * (*iterBead)->m_Force[0] ;
                (*iterBead)->m_Mom[1] = (*iterBead)->m_Mom[1] -= m_halfdt* (*iterBead)->m_Force[1];
                (*iterBead)->m_Mom[2] = (*iterBead)->m_Mom[2] -= m_halfdt * (*iterBead)->m_Force[2];		

    #ifndef NDEBUG
                (*iterBead)->m_AngMom[0]	= nanf(""); // angular momentum not updated for fast path
    #endif
            }
        }
    }

    void UpdateMomFast(CSimBox *mbox, CCNTCell *cell)
    {
        // Forwarding to CCNTCell members
        auto &m_lBeads=cell->m_lBeads;

        /////////////////////////////////////////////////////////////
        // Begin the adapted body from CCNTCell::UpdateMom
        
        for( AbstractBeadVectorIterator iterBead=m_lBeads.begin(); iterBead!=m_lBeads.end(); iterBead++ )
        {
            // Note that moveable calculations have been moved into the force update loop,
            // so here we only read it.
            if((*iterBead)->GetMovable())	// flag ignored by immovable beads
            {
                // We have already updated m_Mom to m_oldMom + m_halfdt * m_oldForce
                (*iterBead)->m_Mom[0] = (*iterBead)->m_Mom[0] + m_halfdt * (*iterBead)->m_Force[0] ;
                (*iterBead)->m_Mom[1] = (*iterBead)->m_Mom[1] + m_halfdt* (*iterBead)->m_Force[1];
                (*iterBead)->m_Mom[2] = (*iterBead)->m_Mom[2] + m_halfdt * (*iterBead)->m_Force[2];		

    #ifndef NDEBUG
                (*iterBead)->m_AngMom[0]	= nanf(""); // angular momentum not updated for fast path
    #endif

                if(StateLogger::IsEnabled()){
                    StateLogger::LogBead("x_next", (*iterBead)->GetId()-1, (*iterBead)->m_Pos);
                    StateLogger::LogBead("v_next", (*iterBead)->GetId()-1, (*iterBead)->m_Mom);
                    StateLogger::LogBead("f_next", (*iterBead)->GetId()-1, (*iterBead)->m_Force);
                }
            }
        }
    }

    void AddForceFast(CSimBox *mbox, CBond *bond)
    {
        DEBUG_ASSERT( EnableParallelSimBox != SimMPSEnabled);
        // Serial code uses the unPBC coordinates to ensure bonds that span the SimBox boundaries have the correct length

        auto &m_pHead=bond->m_pHead;
        auto &m_pTail=bond->m_pTail;
        auto &dx = bond->m_dx;
        auto &dy = bond->m_dy;
        auto &dz = bond->m_dz;
        auto &Length=bond->m_Length;
        
        dx = m_pHead->GetunPBCXPos() - m_pTail->GetunPBCXPos();
        dy = m_pHead->GetunPBCYPos() - m_pTail->GetunPBCYPos();
        dz = m_pHead->GetunPBCZPos() - m_pTail->GetunPBCZPos();

        Length = sqrt(dx*dx + dy*dy + dz*dz);
        
        // The overall minus sign is used to swap the order of m_Length - m_UnStrLen

        double SprConst=bond->m_SprConst;
        double UnStrLen=bond->m_UnStrLen;
        double fx = SprConst*(UnStrLen - Length)*dx/Length;
        double fy = SprConst*(UnStrLen - Length)*dy/Length;
        double fz = SprConst*(UnStrLen - Length)*dz/Length;

        if(StateLogger::IsEnabled()){
            unsigned id1=m_pHead->GetId()-1, id2=m_pTail->GetId()-1;
            StateLogger::LogBeadPair("bond_dx", id1, id2, {dx, dy, dz});
            StateLogger::LogBeadPair("bond_r", id1, id2, Length);
            StateLogger::LogBeadPair("bond_f", id1, id2, {fx, fy, fz});
        }

        m_pHead->m_Force[0]+= fx;
        m_pHead->m_Force[1]+= fy;
        m_pHead->m_Force[2]+= fz;

        m_pTail->m_Force[0]-= fx;
        m_pTail->m_Force[1]-= fy;
        m_pTail->m_Force[2]-= fz;

        // No stresses calculated
    #ifndef NDEBUG
        m_pHead->m_Stress[0] = nanf("");
    #endif
    }

    void AddBondForcesFast(CSimBox *box)
    {
        auto &m_vAllPolymers=box->m_vAllPolymers;

        for(PolymerVectorIterator iterPoly=m_vAllPolymers.begin(); iterPoly!=m_vAllPolymers.end(); iterPoly++)
        {
            auto &poly=**iterPoly;

            for(BondVectorIterator iterBond=poly.m_vBonds.begin(); iterBond!=poly.m_vBonds.end(); iterBond++)
            {
                AddForceFast(box, (*iterBond));
            }
        }
    }


    void AddForceFast(CSimBox *mbox, CBondPair *bp)
    {
        #ifndef NDEBUG
        bp->m_P2=nanf("");
        bp->m_PhiDiff=nanf("");
        bp->m_PhiDiffSq=nanf("");
        bp->m_PhiDiff=nanf("");
        bp->m_PhiDiffSq=nanf("");
        bp->m_BeadXForce[0]=nanf("");
        bp->m_BeadYForce[0]=nanf("");
        bp->m_BeadZForce[0]=nanf("");
        #endif

        auto &m_pFirst=bp->m_pFirst;
        auto &m_pSecond=bp->m_pSecond;

        double FirstLength   = m_pFirst->m_Length;
        double SecondLength  = m_pSecond->m_Length;

        const double magProduct = FirstLength*SecondLength;

        // Calculate the angle between the bonds: if either bond length is zero
        // then we set the forces to zero.  Note that the bond vectors are defined
        // in the direction: head - tail.

        if(magProduct > 0.0001)
        {
            double first[3];
            first[0] = m_pFirst->m_dx;
            first[1] = m_pFirst->m_dy;
            first[2] = m_pFirst->m_dz;

            double second[3];
            second[0] = m_pSecond->m_dx;
            second[1] = m_pSecond->m_dy;
            second[2] = m_pSecond->m_dz;

            const double b1MagSq		= FirstLength*FirstLength;
            const double b2MagSq		= SecondLength*SecondLength;
            const double b1Dotb2		= first[0]*second[0] + first[1]*second[1] + first[2]*second[2];
            const double b1b2Overb1Sq	= b1Dotb2/b1MagSq;
            const double b1b2Overb2Sq	= b1Dotb2/b2MagSq;
            const double cosPhiSq		= b1b2Overb1Sq*b1b2Overb2Sq;

            double forceMag = 0.0;

            // Check that the bond angle is not exactly 90 deg but allow the cosine to be < 0

            if(fabs(b1Dotb2) > 0.000001)
            {			
                #ifndef NDEBUG
                bp->m_Prefactor=nanf("");
                #endif
                
                // Add the restoring force depending on whether there is a preferred angle
                // for the bond pair or not

                if(bp->m_Phi0 > 0.0)
                {
                    double Prefactor = sqrt(1.0/cosPhiSq - 1.0);
                    forceMag = bp->m_Modulus*(bp->m_CosPhi0 - bp->m_SinPhi0/Prefactor)/magProduct;
                }
                else
                {
                    forceMag = bp->m_Modulus/magProduct;
                }
            }
            else
            {
                forceMag	= bp->m_Modulus/magProduct;
            }


            double BeadXForce[3], BeadYForce[3], BeadZForce[3];

            BeadXForce[0] = forceMag*(b1b2Overb1Sq*first[0] - second[0]);
            BeadYForce[0] = forceMag*(b1b2Overb1Sq*first[1] - second[1]);
            BeadZForce[0] = forceMag*(b1b2Overb1Sq*first[2] - second[2]);

            m_pFirst->m_pTail->m_Force[0] += BeadXForce[0];
            m_pFirst->m_pTail->m_Force[1] += BeadYForce[0];
            m_pFirst->m_pTail->m_Force[2] += BeadZForce[0];

            BeadXForce[2] = forceMag*(first[0] - b1b2Overb2Sq*second[0]);
            BeadYForce[2] = forceMag*(first[1] - b1b2Overb2Sq*second[1]);
            BeadZForce[2] = forceMag*(first[2] - b1b2Overb2Sq*second[2]);

            m_pSecond->m_pHead->m_Force[0] += BeadXForce[2];
            m_pSecond->m_pHead->m_Force[1] += BeadYForce[2];
            m_pSecond->m_pHead->m_Force[2] += BeadZForce[2];
            
            BeadXForce[1] = -BeadXForce[0] - BeadXForce[2];
            BeadYForce[1] = -BeadYForce[0] - BeadYForce[2];
            BeadZForce[1] = -BeadZForce[0] - BeadZForce[2];

            m_pFirst->m_pHead->m_Force[0] += BeadXForce[1];
            m_pFirst->m_pHead->m_Force[1] += BeadYForce[1];
            m_pFirst->m_pHead->m_Force[2] += BeadZForce[1];

            if(StateLogger::IsEnabled()){
                unsigned idH=m_pFirst->m_pTail->GetId()-1, idM=m_pFirst->m_pHead->GetId()-1, idT=m_pSecond->m_pHead->GetId()-1;
                //DEBUG_ASSERT(idH!=idM && idM!=idT && idT!=idH);
                //DEBUG_ASSERT( idM == m_pSecond->m_pTail->GetId()-1 );
                StateLogger::LogBeadTriple("bondpair_dx1", idH, idM, idT, {m_pFirst->m_dx, m_pFirst->m_dy, m_pFirst->m_dz});
                StateLogger::LogBeadTriple("bondpair_dx2", idH, idM, idT, {m_pSecond->m_dx, m_pSecond->m_dy, m_pSecond->m_dz});

                StateLogger::LogBeadTriple("bondpair_fH", idH, idM, idT, {BeadXForce[0], BeadYForce[0], BeadZForce[0] });
                StateLogger::LogBeadTriple("bondpair_fM", idH, idM, idT, {BeadXForce[1], BeadYForce[1], BeadZForce[1] });
                StateLogger::LogBeadTriple("bondpair_fT", idH, idM, idT, {BeadXForce[2], BeadYForce[2], BeadZForce[2] });
            }
        }
    }

    void AddBondPairForcesFast(CSimBox *box)
    {
        auto &m_vAllPolymers=box->m_vAllPolymers;

        for(PolymerVectorIterator iterPoly=m_vAllPolymers.begin(); iterPoly!=m_vAllPolymers.end(); iterPoly++)
        {
            auto &poly=**iterPoly;

            for(BondPairVectorIterator iterBP=poly.m_vBondPairs.begin(); iterBP!=poly.m_vBondPairs.end(); iterBP++)
            {
                AddForceFast(box, (*iterBP));
            }
        }
    }

    void EvolveFast(CSimBox *box, unsigned sim_start_time, unsigned nSteps)
    {
        auto &m_vCNTCells=box->m_vCNTCells;

        CNTCellIterator iterCell;  // used in all three loops below

        for(iterCell=m_vCNTCells.begin(); iterCell!=m_vCNTCells.end(); iterCell++)
        {
            UpdateMomFastReverse(box, (*iterCell));
        }

        for(unsigned i=0; i<nSteps; i++){
            StateLogger::BeginStep(sim_start_time+i);

            CCNTCell::PreCalculateDPDForces(box->GetRNGSeed(), sim_start_time+i);

            iterCell = m_vCNTCells.begin();
            PrefetchHint(**iterCell);

            while(iterCell != m_vCNTCells.end()){
                auto curr=*iterCell;

                ++iterCell;
                if(iterCell!=m_vCNTCells.end()){
                    PrefetchHint(**iterCell);
                }

                UpdateMomThenPosFastV2(box, curr);
            }

            // Next calculate the forces between all pairs of beads in NN CNT cells
            // that can potentially interact. No monitor accumulations are performed.

            if(!IsParallel()){
                iterCell = m_vCNTCells.begin();
                PrefetchHint(**iterCell);

                while(iterCell != m_vCNTCells.end()){
                    auto curr=*iterCell;

                    ++iterCell;
                    if(iterCell!=m_vCNTCells.end()){
                        PrefetchHint(**iterCell);
                    }

                    UpdateForceFast(curr);
                }
            }else{
                range_3d box{
                  {0, m_CNTXCellNo},
                  {0, m_CNTYCellNo},
                  {0, m_CNTZCellNo}  
                };

                ParallelContext::GetDefaultContext().ParForWithSafeHalo(
                    box,
                    2, // Means minimum parallelism comes with dimension of 2*2*2 = 8. Too fine-grain?
                    [&](const range_3d &local){
                        for(auto p : local){
                            unsigned index = box.to_linear_index(p);
                            auto cell=m_vCNTCells[index];
                            DEBUG_ASSERT(cell->GetBLXIndex()==p.x);
                            DEBUG_ASSERT(cell->GetBLYIndex()==p.y);
                            DEBUG_ASSERT(cell->GetBLZIndex()==p.z);
                            UpdateForceFast(cell);
                        }
                    }
                );
            }

            CCNTCell::PostCalculateDPDForces();

            if(StateLogger::IsEnabled()){
                for(auto b : box->GetAllBeadsInCNTCells()){
                    StateLogger::LogBead("dpd_f_total", b->GetId()-1, b->m_Force);
                }
            }

            // Add in the forces between bonded beads and the stiff bond force. Note that
            // AddBondPairForces() must be called after AddBondForces() because it relies
            // on the bond lengths having already been calculated in CBond::AddForce().

            AddBondForcesFast(box);
            AddBondPairForcesFast(box);
        }

        for(iterCell=m_vCNTCells.begin(); iterCell!=m_vCNTCells.end(); iterCell++)
        {
            // Any logging here appears for the final time-step
            UpdateMomFast(box, (*iterCell));
        } 
    }
};


class SimEngineFastPar
    : public SimEngineFast
{
public:

    std::string Name() const override
    {
        return "SimEngineFastPar";
    }

    virtual bool IsParallel() const override
    { return true; }
};

static bool reg_SimEngineFast = SimEngineBase<SimEngineFast>::Register();
static bool reg_SimEngineFastPar = SimEngineBase<SimEngineFastPar>::Register();
