#include <vector>
#include <atomic>
#include <vector>
#include <mutex>
#include <cassert>
#include <cstdlib>
#include <set>
#include <cstring>
#include <unordered_map>

#include "StdAfx.h"
#include "SimDefs.h"
#include "xxBase.h"
#include "AbstractBead.h"
#include "CNTCell.h"

#include "ISimEngine.h"

#include "SimBox.h"
#include "mpsSimBox.h"

#include "ISimBox.h"


struct SimEngineFast
    : public ISimEngine
{
public:
    std::string Name() const override
    {
        return "SimEngineFast";
    }

    std::string CanSupport(const ISimBox *box) const override
    {
        return {};
    }

    bool IsParallel() const override
    { return false; }

    void Run(ISimBox *box, bool modified, unsigned num_steps) override
    {
        CSimBox *mbox=const_cast<CSimBox*>(box->GetSimBox());
        mbox->EvolveFast(num_steps);
    }

private:
    //#pragma GCC push_options
    //#pragma GCC optimize("fast-math")
    void UpdateForceFast(CCNTCell *cell)
    {
        #if (SimDimension!=3)
        ErrorTrace("Compile-time conditions for UpdateForceFast not met - dimensions.");
        exit(1);
        #elif (EnableStressTensorSphere == SimMiscEnabled)
        ErrorTrace("Compile-time conditions for UpdateForceFast not met - EnableStressTensorSphere.");
        exit(1);
        #elif (EnableParallelSimBox != SimMPSDisabled)
        ErrorTrace("Compile-time conditions for UpdateForceFast not met - parallel simb box.");
        exit(1);
        #elif defined(UseDPDBeadRadii)
        ErrorTrace("Compile-time conditions for UpdateForceFast not met - bead radii.");
        exit(1);
        #endif

        mpsSimBox::GlobalCellCounter++;  // increment the counter for intra-cell force calculations
        
        long  localCellCellCounter = 0;

        if(cell->m_lBeads.empty()){
            localCellCellCounter = 13;
            mpsSimBox::GlobalCellCellIntCounter += localCellCellCounter;
            return;
        }

        cell->m_aIntNNCells[0]->PrefetchHint();

        double rng_scale =  CCNTCell::m_invrootdt * CCNTCell::m_Inv2Power32;

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
        int numLocal=cell->m_lBeads.size();
        iterBead1=cell->m_lBeads.begin();
        for(int ii1=0; ii1 < numLocal-1; ii1++, iterBead1++){
            iterBead2 = std::next(iterBead1);
            for(int ii2=ii1+1; ii2 < cell->m_lBeads.size(); ii2++, iterBead2++){

                double dx2[3];
                for(int d=0; d<3; d++){
                    dx[d] = ((*iterBead1)->m_Pos[d] - (*iterBead2)->m_Pos[d]);
                    dx2[d] = dx[d] * dx[d];
                }
        
                dr2 = dx2[0] + dx2[1] + dx2[2];

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

                    conForce  = cell->m_vvConsInt[(*iterBead1)->GetType()][(*iterBead2)->GetType()]*wr;

    // Dissipative and random force magnitudes. Note dr factor in newForce calculation

                    rdotv		= (dx_dv[0] + dx_dv[1] + dx_dv[2]) * inv_dr;
                    gammap		= cell->m_vvDissInt[(*iterBead1)->GetType()][(*iterBead2)->GetType()]*wr2;

                    dissForce	= -gammap*rdotv;				
                    randForce	= sqrt(gammap)*cell->RandUnifScaled(rng_scale);
                    normTotalForce = (conForce + dissForce + randForce) * inv_dr;

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
                m_aIntNNCells[i+1]->PrefetchHint();
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

                    if( dr2 < 1.0 && dr2 > min_r2)
                    {		
                        for(int d=0; d<3; d++){
                            dv[d] = ((*iterBead1)->m_Mom[d] - (*iterBead2)->m_Mom[d]);
                        }

                        dr = sqrt(dr2);
                    
                        wr = (1.0 - dr);
                        wr2 = wr*wr;
                        double inv_dr = 1.0/dr;

                        conForce	= m_vvConsInt[(*iterBead1)->GetType()][(*iterBead2)->GetType()]*wr;

                        rdotv		= (dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2]) * inv_dr;
                        gammap		= m_vvDissInt[(*iterBead1)->GetType()][(*iterBead2)->GetType()]*wr2;

                        dissForce	= -gammap*rdotv;		
                        randForce	= sqrt(gammap)*RandUnifScaled(rng_scale);
                        normTotalForce = (conForce + dissForce + randForce) * inv_dr;

                        for(int d=0; d<3; d++){
                            newForce[d] = normTotalForce * dx[d];
                            (*iterBead1)->m_Force[d] += newForce[d];
                            (*iterBead2)->m_Force[d] -= newForce[d];
                        }
                    }
                }
            }
        }

    // Divide the local cell-cell counter by the number of beads in this cell and add the result to the global cell-cell counter
        
        assert(m_lBeads.size() > 0);
        localCellCellCounter /= m_lBeads.size();
        
        mpsSimBox::GlobalCellCellIntCounter += localCellCellCounter;
    }
    //#pragma GCC pop_options

    #if 0
    void CCNTCell::UpdatePosFast()
    {
        #ifndef NDEBUG
        if(m_lambda!=0.5){
            ErrorTrace("Attempt to call UpdatePosFast with lambda!=0.5");
            exit(1);
        }
        #endif

        enum DirIndex : char{
            UTR = 26,
            DTR = 8,
            _TR = 17,
            UBR = 20,
            DBR = 2,
            _BR = 11,
            U_R = 23,
            D_R = 5,
            __R = 14,
            UTL = 24,
            DTL = 6,
            _TL = 15,
            UBL = 18,
            DBL = 0,
            _BL = 9,
            U_L = 21,
            D_L = 3,
            __L = 12,
            UT_ = 25,
            DT_ = 7,
            _T_ = 16,
            UB_ = 19,
            DB_ = 1,
            _B_ = 10,
            U__ = 22,
            D__ = 4,
            ___ = -1
        };

        double dx[3];

        unsigned index=0;

        while(index < m_lBeads.size()){
            assert(!m_lBeads.empty());

            CAbstractBead *bead=m_lBeads[index];

            // Only allow bead to move if its IsMovable flag is true. This allows
            // us to indicate when a bead has already crossed a cell boundary and
            // should not be moved again in this timestep.

            if(!bead->GetMovable()){
                index++;
                continue;
            }
            
            /*for(int d=0; d<3; d++){
                assert( m_BLCoord[d] <= bead->m_Pos[d] );
                assert( m_TRCoord[d] >= bead->m_Pos[d] );
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
            assert(m_lamdt == m_halfdt);
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
            int deltaParts[3];
            int moved=0;
            for(int d=0; d<3; d++){
                deltaParts[d] = (bead->m_Pos[d] > m_TRCoord[d]) - (bead->m_Pos[d] < m_BLCoord[d]);
                moved |= deltaParts[d];
            }
            if(!moved){
                for(int d=0; d<3; d++){
                    assert( m_BLCoord[d] <= bead->m_Pos[d] );
                    assert( m_TRCoord[d] >= bead->m_Pos[d] );
                }

                // If the bead did not change cells increment the
                // iterator by hand.
                index++;
                continue;
            }

            // delta[0] : 1 = R, -1 = L
            // delta[1] : 1 = T, -1 = B
            // delta[2] : 1 = U, -1 = D
            //int deltaFull=1 + deltaParts[0] + 3 + 3*deltaParts[1] + 9 + 9*deltaParts[2];
            int deltaFull=(1+3+9) + deltaParts[0] + 3*deltaParts[1] + 9*deltaParts[2];
            static const DirIndex indexMapping[27] = {
                DBL, DB_, DBR,
                D_L, D__, D_R,
                DTL, DT_, DTR,

                _BL, _B_, _BR,
                __L, ___, __R,
                _TL, _T_, _TR,

                UBL, UB_, UBR,
                U_L, U__, U_R,
                UTL, UT_, UTR,
            };

            DirIndex directDir = indexMapping[deltaFull];

            CCNTCell *destCell = m_aNNCells[directDir];
            if(m_bExternal && destCell->IsExternal()){
                bead->m_Pos[0] += m_SimBoxXLength * ( (bead->m_Pos[0] < 0) - (bead->m_Pos[0] >= m_SimBoxXLength) );
                bead->m_Pos[1] += m_SimBoxYLength * ( (bead->m_Pos[1] < 0) - (bead->m_Pos[1] >= m_SimBoxYLength) );
                bead->m_Pos[2] += m_SimBoxZLength * ( (bead->m_Pos[2] < 0) - (bead->m_Pos[2] >= m_SimBoxZLength) );
            }

            m_aNNCells[directDir]->m_lBeads.push_back(bead);
            m_lBeads[index] = m_lBeads.back();
            m_lBeads.pop_back();

        //	assert(dir == directDir);

            /*for(int i=0; i<27; i++){
                CCNTCell *c=m_aNNCells[i];
                for(auto lbead : c->m_lBeads){
                    for(int d=0; d<3; d++){
                        assert( c->m_BLCoord[d] <= lbead->m_Pos[d] );
                        assert( c->m_TRCoord[d] >= lbead->m_Pos[d] );
                    }
                }
            }*/
        }

        /*for(int i=0; i<27; i++){
            CCNTCell *c=m_aNNCells[i];
            for(auto bead : c->m_lBeads){
                for(int d=0; d<3; d++){
                    assert( c->m_BLCoord[d] <= bead->m_Pos[d] );
                    assert( c->m_TRCoord[d] >= bead->m_Pos[d] );
                }
            }
        }*/
    }
    #endif

    void CCNTCell::UpdateMomThenPosFastV2()
    {
        #ifndef NDEBUG
        if(m_lambda!=0.5){
            ErrorTrace("Attempt to call UpdatePosFast with lambda!=0.5");
            exit(1);
        }
        #endif

        double dx[3];

        int index=m_lBeads.size()-1;

        while(index >= 0){
            assert(!m_lBeads.empty());

            CAbstractBead *bead=m_lBeads[index];

            // Only allow bead to move if its IsMovable flag is true. This allows
            // us to indicate when a bead has already crossed a cell boundary and
            // should not be moved again in this timestep.

            if(!bead->GetMovable()){
                index++;
                continue;
            }
            
            /*for(int d=0; d<3; d++){
                assert( m_BLCoord[d] <= bead->m_Pos[d] );
                assert( m_TRCoord[d] >= bead->m_Pos[d] );
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
            assert(m_lamdt == m_halfdt);
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
                moved |= (bead->m_Pos[d] > m_TRCoord[d]) - (bead->m_Pos[d] < m_BLCoord[d]);
            }
            if(!moved){
                for(int d=0; d<3; d++){
                    assert( m_BLCoord[d] <= bead->m_Pos[d] );
                    assert( m_TRCoord[d] >= bead->m_Pos[d] );
                }

                // If the bead did not change cells increment the
                // iterator by hand.
                index++;
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

            assert( fmod(m_SimBoxXLength, 1) == 0);
            assert( fmod(m_SimBoxYLength, 1) == 0);
            assert( fmod(m_SimBoxZLength, 1) == 0);
            assert( m_CNTXCellWidth == 1 && m_CNTYCellWidth==1 && m_CNTZCellWidth==1);

            int ix=floor(bead->m_Pos[0]), iy=floor(bead->m_Pos[1]), iz=floor(bead->m_Pos[2]);
            int cell_index = m_CNTXCellNo*(m_CNTYCellNo*iz+iy) + ix; 

            if(index+1 < m_lBeads.size()){
                std::swap(m_lBeads[index], m_lBeads.back());
            }
            m_lBeads.pop_back();

            m_pISimBox->GetSimBox()->AddBeadToCNTCell(cell_index, bead);
        }
    }


    void UpdateMomThenPosFast()
    {
        #ifndef NDEBUG
        if(m_lambda!=0.5){
            ErrorTrace("Attempt to call UpdateMomThenPosFast with lambda!=0.5");
            exit(1);
        }
        #endif

        enum DirIndex : char{
            UTR = 26,
            DTR = 8,
            _TR = 17,
            UBR = 20,
            DBR = 2,
            _BR = 11,
            U_R = 23,
            D_R = 5,
            __R = 14,
            UTL = 24,
            DTL = 6,
            _TL = 15,
            UBL = 18,
            DBL = 0,
            _BL = 9,
            U_L = 21,
            D_L = 3,
            __L = 12,
            UT_ = 25,
            DT_ = 7,
            _T_ = 16,
            UB_ = 19,
            DB_ = 1,
            _B_ = 10,
            U__ = 22,
            D__ = 4,
            ___ = -1
        };

        double dx[3];

        unsigned index=0;

        while(index < m_lBeads.size()){
            assert(!m_lBeads.empty());

            CAbstractBead *bead=m_lBeads[index];

            // Only allow bead to move if its IsMovable flag is true. This allows
            // us to indicate when a bead has already crossed a cell boundary and
            // should not be moved again in this timestep.

            if(!bead->GetMovable()){
                index++;
                continue;
            }
            
            /*for(int d=0; d<3; d++){
                assert( m_BLCoord[d] <= bead->m_Pos[d] );
                assert( m_TRCoord[d] >= bead->m_Pos[d] );
            }*/

            bead->SetNotMovable();

            // Apply the mom from the end of previous step (loop skewed)
            bead->m_Mom[0] = bead->m_Mom[0] + m_halfdt * bead->m_Force[0] ;
            bead->m_Mom[1] = bead->m_Mom[1] + m_halfdt* bead->m_Force[1];
            bead->m_Mom[2] = bead->m_Mom[2] + m_halfdt * bead->m_Force[2];


            // We do not write to m_oldPos, m_oldMom, or m_oldForce, as it is not needed in fast path (causes extra memory traffic)
    #ifndef NDEBUG
            bead->m_oldPos[0] = nanf("");
            bead->m_oldMom[0] = nanf("");
            bead->m_oldForce[0] = nanf("");
    #endif		

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
            assert(m_lamdt == m_halfdt);
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
            int deltaParts[3];
            int moved=0;
            for(int d=0; d<3; d++){
                deltaParts[d] = (bead->m_Pos[d] > m_TRCoord[d]) - (bead->m_Pos[d] < m_BLCoord[d]);
                moved |= deltaParts[d];
            }
            if(!moved){
                for(int d=0; d<3; d++){
                    assert( m_BLCoord[d] <= bead->m_Pos[d] );
                    assert( m_TRCoord[d] >= bead->m_Pos[d] );
                }

                // If the bead did not change cells increment the
                // iterator by hand.
                index++;
                continue;
            }

            // delta[0] : 1 = R, -1 = L
            // delta[1] : 1 = T, -1 = B
            // delta[2] : 1 = U, -1 = D
            //int deltaFull=1 + deltaParts[0] + 3 + 3*deltaParts[1] + 9 + 9*deltaParts[2];
            int deltaFull=(1+3+9) + deltaParts[0] + 3*deltaParts[1] + 9*deltaParts[2];
            static const DirIndex indexMapping[27] = {
                DBL, DB_, DBR,
                D_L, D__, D_R,
                DTL, DT_, DTR,

                _BL, _B_, _BR,
                __L, ___, __R,
                _TL, _T_, _TR,

                UBL, UB_, UBR,
                U_L, U__, U_R,
                UTL, UT_, UTR,
            };

            DirIndex directDir = indexMapping[deltaFull];

            CCNTCell *destCell = m_aNNCells[directDir];
            if(m_bExternal && destCell->IsExternal()){
                bead->m_Pos[0] += m_SimBoxXLength * ( (bead->m_Pos[0] < 0) - (bead->m_Pos[0] >= m_SimBoxXLength) );
                bead->m_Pos[1] += m_SimBoxYLength * ( (bead->m_Pos[1] < 0) - (bead->m_Pos[1] >= m_SimBoxYLength) );
                bead->m_Pos[2] += m_SimBoxZLength * ( (bead->m_Pos[2] < 0) - (bead->m_Pos[2] >= m_SimBoxZLength) );
            }

            m_aNNCells[directDir]->m_lBeads.push_back(bead);
            m_lBeads[index] = m_lBeads.back();
            m_lBeads.pop_back();

        //	assert(dir == directDir);

            /*for(int i=0; i<27; i++){
                CCNTCell *c=m_aNNCells[i];
                for(auto lbead : c->m_lBeads){
                    for(int d=0; d<3; d++){
                        assert( c->m_BLCoord[d] <= lbead->m_Pos[d] );
                        assert( c->m_TRCoord[d] >= lbead->m_Pos[d] );
                    }
                }
            }*/
        }

        /*for(int i=0; i<27; i++){
            CCNTCell *c=m_aNNCells[i];
            for(auto bead : c->m_lBeads){
                for(int d=0; d<3; d++){
                    assert( c->m_BLCoord[d] <= bead->m_Pos[d] );
                    assert( c->m_TRCoord[d] >= bead->m_Pos[d] );
                }
            }
        }*/
    }


    void UpdateMomFastReverse(CCNTCell *cell)
    {
        #if SimIdentifier == BD
        ErrorTrace("Attempt to call UpdateMomFast with DB");
        exit(1);
        #endif
        #ifndef NDEBUG
        if(CCNTCell::m_lambda!=0.5){
            ErrorTrace("Attempt to call UpdateMomFast with lambda!=0.5");
            exit(1);
        }
        #endif
        
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

    void CCNTCell::UpdateMomFast()
    {
        #if SimIdentifier == BD
        ErrorTrace("Attempt to call UpdateMomFast with DB");
        exit(1);
        #endif
        #ifndef NDEBUG
        if(m_lambda!=0.5){
            ErrorTrace("Attempt to call UpdateMomFast with lambda!=0.5");
            exit(1);
        }
        #endif
        
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
            }
        }
    }
};

static bool reg_SimEngineFast = SimEngineBase<SimEngineFast>::Register();