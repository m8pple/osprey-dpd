// BeadChargeWrapper.h: interface for the CBeadChargeWrapper class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BEADCHARGEWRAPPER_H__9662103E_8034_4119_8CB7_21488FA8D8C2__INCLUDED_)
#define AFX_BEADCHARGEWRAPPER_H__9662103E_8034_4119_8CB7_21488FA8D8C2__INCLUDED_


// Include files to give access to the underlying CAbstractBead object

#include "AbstractBead.h"

#include "xxBase.h"

class CBeadChargeWrapper  
{
	// The SimBox must be a friend so that it can create the charged beads
	// and tell them to calculate their forces. The constructor is private 
	// so that this wrapper class cannot be instantiated elsewhere, and the
	// destructor is not virtual because we have no base classes.
	// Other command classes that need to modify charged beads must also be
	// declared as friends here.

	friend class CSimBox;
	friend class ccChargeBeadByTypeImpl;
	friend class ccUnchargeBeadByTypeImpl;

public:
	~CBeadChargeWrapper();


	// Forwarding functions to the contained CAbstractBead object

	inline long	  GetId()			const {return m_pBead->GetId();}
	inline long	  GetType()			const {return m_pBead->GetType();}
	inline bool	  GetMovable()		const {return m_pBead->GetMovable();}
	inline double GetXPos()			const {return m_pBead->GetXPos();}
	inline double GetYPos()			const {return m_pBead->GetYPos();}
	inline double GetZPos()			const {return m_pBead->GetZPos();}
	inline double GetXForce()		const {return m_pBead->GetXForce();}
	inline double GetYForce()		const {return m_pBead->GetYForce();}
	inline double GetZForce()		const {return m_pBead->GetZForce();}

private:
	CBeadChargeWrapper(CAbstractBead* const pBead, double strength, double range, double kappa);

	void AddForce(CBeadChargeWrapper* const pChargedBead2, double dx[3]);


private:
	CAbstractBead* const m_pBead;		// Contained bead pointer

	double m_Strength;					// Strength of screened charge force
	double m_Range;						// Range of force
	double m_Kappa;						// Inverse range of force

};

typedef xxBaselist<CBeadChargeWrapper*>::iterator				ChargedBeadListIterator;
typedef xxBaselist<CBeadChargeWrapper*>::const_iterator			cChargedBeadListIterator;
typedef xxBaselist<CBeadChargeWrapper*>::reverse_iterator		rChargedBeadListIterator;
typedef xxBaselist<CBeadChargeWrapper*>::const_reverse_iterator	crChargedBeadListIterator;

typedef xxBasevector<CBeadChargeWrapper*>::iterator					ChargedBeadVectorIterator;
typedef xxBasevector<CBeadChargeWrapper*>::const_iterator			cChargedBeadVectorIterator;
typedef xxBasevector<CBeadChargeWrapper*>::reverse_iterator			rChargedBeadVectorIterator;
typedef xxBasevector<CBeadChargeWrapper*>::const_reverse_iterator	crChargedBeadVectorIterator;


#endif // !defined(AFX_BEADCHARGEWRAPPER_H__9662103E_8034_4119_8CB7_21488FA8D8C2__INCLUDED_)
