#if !defined(AFX_ccOptimisePolymerOrderingForPDPD_H__B75BC544_E47F_11D3_BF23_004095086186__INCLUDED_)
#define AFX_ccOptimisePolymerOrderingForPDPD_H__B75BC544_E47F_11D3_BF23_004095086186__INCLUDED_


#include "xxCommand.h"

/* In POETS DPD we pack bead type and polymer identity into 32-bits. As a result
it is important to try to give multi-bead polymers low indices, with monomers
(e.g. water) at high indices.
*/
class ccOptimisePolymerOrderingForPDPD : public xxCommand  
{
public:
	ccOptimisePolymerOrderingForPDPD(long executionTime);
	ccOptimisePolymerOrderingForPDPD(const ccOptimisePolymerOrderingForPDPD& oldCommand);
	virtual ~ccOptimisePolymerOrderingForPDPD();

	// Access functions to the command's data

	zOutStream& put(zOutStream& os) const;
	zInStream&  get(zInStream& is);

	virtual bool Execute(long simTime, ISimCmd* const pISimCmd) const;

	static const zString GetType();	// Return the type of command
	static long  GetArgumentTotal(); // Return the number of arguments 

	virtual const xxCommand* GetCommand() const;

	virtual bool IsDataValid(const CInputData& riData) const;

protected:
	virtual const zString GetCommandType() const;

private:
	static const zString m_Type;	// Identifier used in control data file for command
};

#endif

