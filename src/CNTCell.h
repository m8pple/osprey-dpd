// CNTCell.h: interface for the CCNTCell class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CNTCELL_H__E36E3E80_32EC_11D3_820E_0060088AD300__INCLUDED_)
#define AFX_CNTCELL_H__E36E3E80_32EC_11D3_820E_0060088AD300__INCLUDED_



// Forward declarations

class CMonitor;
class ISimBox;
class mpsBorder;


#include "AbstractCell.h"
#include "AbstractBead.h"

class CCNTCell : public CAbstractCell  
{
	// friend functions of the CMonitor class needed to calculate observables
	// and time series data for output. And the CExternalCNTCell class needs
    // to copy data from the original CNT cells.

	friend class CMonitor;
    friend class CExternalCNTCell;

	// This provides faster versions of the functions in CCNTCell wrapped into a SimEngine
	// It directly accesses members of CCNTCell, but is kept seperate.
	friend class SimEngineFast;

	// ****************************************
	// Construction/Destruction: base class has protected constructor
public:

	CCNTCell();
	CCNTCell(long index, long coords[3], bool flag);

    // Copy constructor
    CCNTCell(const CCNTCell& oldCell);

	virtual ~CCNTCell();


	// ****************************************
	// Global functions, static member functions and variables
public:

	static double GetRandomNo();
	static double GetExponentialRandomNo();
	static double GetGaussRandomNo();
	static double Randf();
	static double RandUniformBetweenBeads(uint32_t bead_id1, uint32_t bead_id2); // Returns random force in [-0.5,0.5) between two beads
	static double RandUniformBetweenBeads(const CAbstractBead *bead1, const CAbstractBead *bead2); // Returns random force in [-0.5,0.5) between two beads
	static double Gasdev();
	static double Expdev();
	static double GetLambda();
	static double GetKt();

	static const zString GetRandomString();
	static const zString GetRandomString(const zString prefix, const zString separator);
	static const zString GetRandomString(const zString prefix, const zString separator, long counter);
	static const zString GetRandomAlphabeticString(long size);

    // ****************************************
	// Functions to initialise constants used in the simulation

	static void SetMonitorPointer(CMonitor* pMon);

	static void SetISimBoxPointer(ISimBox* pISimBox);

	static void SetRNGSeed(long idum);

	static void SetSimBoxLengths(long nx, long ny, long nz, double cntlx, double cntly, double cntlz);

	static void SetTimeStepConstants(double dt, double lambda, double cutoffradius, double kT);

	static void SetBDBeadStructure(const zArray2dDouble* pvvConsInt, const zArray2dDouble* pvvDissInt);

    // Standard DPD conservative and dissipative forces
	static void SetDPDBeadStructure(const zArray2dDouble* pvvConsInt, const zArray2dDouble* pvvDissInt);

    // Overloaded function to include the density-dependent force parameters
	static void SetDPDBeadStructure(const zArray2dDouble* pvvConsInt, const zArray2dDouble* pvvDissInt, 
                                    const zArray2dDouble* pvvLGInt);

	static void SetMDBeadStructure(const zArray2dDouble* pvvLJDepth, const zArray2dDouble* pvvLJRange,
								   const zArray2dDouble* pvvSCDepth, const zArray2dDouble* pvvSCRange);

	// ****************************************
	// Command handler functions that implement commands directed to the CSimBox to 
	// modify the evolution of the simulation.
	// 
	// Change the interactions of a DPD bead type dynamically

	static void SetDPDBeadConsInt(long firstType, long secondType, double newValue);
	static void SetDPDBeadDissInt(long firstType, long secondType, double newValue);
	static void SetDPDBeadLGInt(long firstType, long secondType, double newValue);

	// Allow inspection of current dynamic values CCNTCell is using
	static double GetDPDBeadConsInt(long firstType, long secondType);
	static double GetDPDBeadDissInt(long firstType, long secondType);


	// Change the integration step size and dependent parameters

	static void ChangeTimeStepConstants(const double dt);

	// Toggle the DPD bead forces to zero 

	static void ToggleDPDBeadConservativeForces(bool bZero);
	static void ToggleDPDBeadForces(const double dt);
	static void ToggleDPDBeadThermostat(bool bZero);

	// Add a new DPD bead type and initialise its interactions

	static void AddDPDBeadType(long oldType);

	// Called before calculating DPD forces in each time step.
	// Used to do any CNTCell global setup
	static void PreCalculateDPDForces(long global_seed, unsigned simTime);

	// Called after calculating DPD forces, just in case there is any tear-down to do.
	static void PostCalculateDPDForces();

	// ****************************************
	// PVFs that must be overridden by all derived classes
public:

    // VF to calculate the forces between beads in two adjacent cells that
    // are both internal to one processor

	virtual void UpdateForceBetweenCells(bool bExternal, CAbstractBead* const pBead);


	// ****************************************
	// Public access functions
public:

	inline bool IsExternal() const {return m_bExternal;}

	// Function to return a pointer to a neighbouring CNT Cell

	inline CCNTCell*  GetNNCell(short int i) const {return m_aNNCells[i];}
	inline CCNTCell*  GetIntNNCell(short int i) const {return m_aIntNNCells[i];}

	// Main functions used by the integration scheme to evolve the simulation.

	void UpdateForce();
	void UpdateMom();
	void UpdatePos();

    // Parallel versions of the updating functions

	void UpdateForceP();

    // Alternative force calculation to include the DPD density-dependent force.

	void UpdateLGForce();
	void UpdateLGDensity();

	// Function to calculate the kinetic and potential energy of the cell including
	// both interactions within the cell and with its immediate neighbours

	void UpdateTotalEnergy(double* const pKinetic, double* const pPotential) const;

	double GetPotentialEnergy(CAbstractBead* pBead) const;	// Used for MC relaxation
	long CellBeadTotal() const;
	void AddBeadtoCell(CAbstractBead* pBead);
	void RemoveBeadFromCell(CAbstractBead* const pBead);
	void RemoveAllBeadsFromCell();
	void SetNNCellIndex(long index, CCNTCell* pCell);
	void SetIntNNCellIndex(long index, CCNTCell* pCell);
	bool CheckBeadsinCell();

	/* 	Hook to allow different RNGs to be inserted into CCNTCell. Used both for testing
		purposes, and to make parts of CCNTCell thread safe.
		If this method is not called, then the original LCG RNG is used.
		- CustomRngProc : returns a number in the range (-0.5,+0.5) for the given bead pair
		- CustomRNGBeginTimeStep : performs any initialisation for a given time-step, returning an opaque state handle
		- CustomRNGEndTimeStep : cleans up resources associated with time-step
	*/
	static void SetCustomRNGProc(
		float (*CustomRNGProc)(uintptr_t state, uint32_t bead_id1, uint32_t bead_id2) = 0,
		uintptr_t (*CustomRNGBeginTimeStep)(uint64_t global_seed, uint64_t step_index) = 0,
		void (*CustomRNGEndTimeStep)(uintptr_t state) = 0,
		bool isThreadSafe = false,    	// If true, then multiple threads can access at the same time
		bool isRepeatable = false 		// If true, then a sequential thread calling (bead_id1,bead_id2) will always get the same value
		
	);
	static bool IsRandUniformBetweenBeadsThreadSafe(); // Returns true if the RandUniformBetweenBeads is thread-safe (though it may still be non-deterministic)
	static bool IsRandUniformBetweenBeadsRepeatable(); // Returns true if the RandUniformBetweenBeads is repeatable (though it may not be thread safe)

	// ****************************************
	// Protected local functions
protected:


	// ****************************************
	// Implementation


	// ****************************************
	// Private functions
private:

	void CheckPBCDrift(const CAbstractBead *bg);

	double GetExternalRandomNumber();  // Helper function to RNG tests

    static uint32_t lcg(uint64_t &state);  // Internal helper function for RNG

	// ****************************************
	// Data members
private:
	
    static bool   m_bReadFileOnce;           // Flag showing if external random number file should be read once only
    static long   m_NextRNIndex;             // Index of next random number to use
    static zDoubleVector m_RandomNumbers;    // Externally-generated random numbers for testing
    static long   m_StringSize;              // Default size of randomly-generated strings
    static const zString m_StringSeparator;  // Separator character used in random strings
    static const zString m_AlphabetChars;    // Set of characters that can be in the first position of a type name

	static double m_SimBoxXLength;
	static double m_SimBoxYLength;
	static double m_SimBoxZLength;
	static double m_HalfSimBoxXLength;
	static double m_HalfSimBoxYLength;
	static double m_HalfSimBoxZLength;

	static long m_CNTXCellNo;
	static long m_CNTYCellNo;
	static long m_CNTZCellNo;
	static double m_CNTXCellWidth;
	static double m_CNTYCellWidth;
	static double m_CNTZCellWidth;

	static double m_dt;
	static double m_lambda;
	static double m_cutoffradius;		// Potential cut-off radius for both DPD and MD
	static double m_coradius2;			// Square of cut-off radius
	static double m_kT;
	static double m_halfdt;
	static double m_halfdt2;
	static double m_invrootdt;
	static double m_lamdt;
    static double m_lgnorm;            // Constant part of DPD density-dependent force prefactor

    static double m_dtoverkt;          // Prefactor of the BD force term: not including diffusion constant
    static double m_dispmag;           // Prefactor of the BD displacement term: not including diffusion constant

    static uint64_t m_RNGSeed;   // 64-bit seed for the lcg RNG
    static long double m_2Power32;         // 2**32
    static long double m_Inv2Power32;      // Inverse of 2**32

	// A hook point allowing a custom RNG to be inserted at runtime.
	// It is only active if m_CustomRNGProc is non-null. When null
	// it will be skipped and should be branch-predicted very accurately.
	static float (*m_CustomRNGProc)(uintptr_t state, uint32_t bead_id1, uint32_t bead_id2);
	static uintptr_t m_CustomRNGState;
	// Called before the forces are calculated in the time step. 
	static uintptr_t (*m_CustomRNGBeginTimeStep)(uint64_t global_seed, uint64_t step_index);
	// Called after all forces are calculate in the time step, just in case state points at something dynamically allocated.
	static void (*m_CustomRNGEndTimeStep)(uintptr_t state);
	static bool m_CustomRNGIsThreadSafe;
	static bool m_CustomRNGIsRepeatable;
	
	static CMonitor* m_pMonitor;	// Pointer to CMonitor to allow on-the-fly analysis
	static ISimBox*  m_pISimBox;	// Pointer to ISimBox to allow on-the-fly analysis

	// Bead internal structure data including bead-bead interaction matrices
	// for all simulation types: DPD, MD. We use static pointers so that we
	// can initialize them but then copy the vectors into local storage in the 
	// constructor to avoid having to dereference the pointers all the time.

	static const zArray2dDouble* m_pvvConsInt;	// DPD
	static const zArray2dDouble* m_pvvDissInt;
	static const zArray2dDouble* m_pvvLGInt;

	static const zArray2dDouble* m_pvvLJDepth;	// MD
	static const zArray2dDouble* m_pvvLJRange;
	static const zArray2dDouble* m_pvvSCDepth;
	static const zArray2dDouble* m_pvvSCRange;

	// Static arrays to hold the above data without requiring a dereference

	static zArray2dDouble m_vvConsInt;		// DPD
	static zArray2dDouble m_vvDissInt;
	static zArray2dDouble m_vvLGInt;

	static zArray2dDouble m_vvConsIntBackup;	// DPD zero force array
	static zArray2dDouble m_vvDissIntBackup;
	static zArray2dDouble m_vvLGIntBackup;

	static zArray2dDouble m_vvLJDepth;	// MD
	static zArray2dDouble m_vvLJRange;
	static zArray2dDouble m_vvSCDepth;
	static zArray2dDouble m_vvSCRange;

	static zArray2dDouble m_vvLJDelta;	// Shift in LJ potential 
	static zArray2dDouble m_vvLJSlope;   // Slope of shifted LJ potential
	static zArray2dDouble m_vvSCDelta;	// Ditto for SC potential
	static zArray2dDouble m_vvSCSlope;

	// Local data members

	bool m_bExternal;

    CCNTCell* m_aNNCells[27];		// Allow for both 2d and 3d
    CCNTCell* m_aIntNNCells[13];
};

#endif // !defined(AFX_CNTCELL_H__E36E3E80_32EC_11D3_820E_0060088AD300__INCLUDED_)
