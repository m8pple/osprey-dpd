#if !defined(AFX_ccRoundBeadProperties_H__B75BC544_E47F_11D3_BF23_004095086186__INCLUDED_)
#define AFX_ccRoundBeadProperties_H__B75BC544_E47F_11D3_BF23_004095086186__INCLUDED_


#include "xxCommand.h"

class ccRoundBeadProperties : public xxCommand  
{
public:
	ccRoundBeadProperties(long executionTime);
	ccRoundBeadProperties(long executionTime, int base, int posFracDigits, int velMantDigits);
	ccRoundBeadProperties(const ccRoundBeadProperties& oldCommand);
	virtual ~ccRoundBeadProperties();

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
	int m_base;
	int m_posFracDigits;
	int m_vecMantDigits;
};

#endif

