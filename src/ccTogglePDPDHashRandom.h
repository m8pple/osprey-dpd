// ccToggleForceLogging.h: interface for the ccToggleForceLogging class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(ccTogglePDPDHashRandom_h)
#define ccTogglePDPDHashRandom_h


#include "xxCommand.h"

class ccTogglePDPDHashRandom : public xxCommand  
{
public:
	ccTogglePDPDHashRandom(long executionTime);
	ccTogglePDPDHashRandom(const ccTogglePDPDHashRandom& oldCommand);

	virtual ~ccTogglePDPDHashRandom();


	// xxCommand functions
	
	zOutStream& put(zOutStream& os) const;
	zInStream&  get(zInStream& is);

	// The following pure virtual functions must be provided by all derived classes
	// so that they may have data read into them given only an xxCommand pointer,
	// respond to the SimBox's request to execute and return the name of the command.

	virtual bool Execute(long simTime, ISimCmd* const pISimCmd) const;

	static const zString GetType();	// Return the type of command

	virtual const xxCommand* GetCommand() const;

	virtual bool IsDataValid(const CInputData& riData) const;

	// Access functions to command's data

	// Local data access functions


protected:
	virtual const zString GetCommandType() const;

private:
	static const zString m_Type;	// Identifier used in control data file for command


};

#endif 

