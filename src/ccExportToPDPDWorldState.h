// ccExportToPDPDWorldState.h: interface for the ccExportToPDPDWorldState class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ccExportToPDPDWorldState_H__B75BC544_E47F_11D3_BF23_004095086186__INCLUDED_)
#define AFX_ccExportToPDPDWorldState_H__B75BC544_E47F_11D3_BF23_004095086186__INCLUDED_


#include "xxCommand.h"

class ccExportToPDPDWorldState : public xxCommand  
{
public:
	ccExportToPDPDWorldState(long executionTime);
	ccExportToPDPDWorldState(long executionTime, zString path);
	ccExportToPDPDWorldState(const ccExportToPDPDWorldState& oldCommand);
	virtual ~ccExportToPDPDWorldState();

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
	zString m_Path; // Name of the file to export to

};

#endif // !defined(AFX_ccExportToPDPDWorldState_H__B75BC544_E47F_11D3_BF23_004095086186__INCLUDED_)
