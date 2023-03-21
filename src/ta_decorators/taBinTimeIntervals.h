// taBinTimeIntervals.h: interface for the taBinTimeIntervals class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TABINTIMEINTERVALS_H__DDA01256_8832_45EC_8067_19F09E9F5182__INCLUDED_)
#define AFX_TABINTIMEINTERVALS_H__DDA01256_8832_45EC_8067_19F09E9F5182__INCLUDED_



#include "taHistogramDecorator.h"

class taBinTimeIntervals : public taHistogramDecorator
{
public:
	// ****************************************
	// Construction/Destruction - protected constructor declared below
public:

    taBinTimeIntervals(const zString label, CCommandTargetNode* const pDec, bool bBinTotal);

	virtual ~taBinTimeIntervals();


    // ****************************************
	// Global functions, static member functions and variables
public:

    static const zString GetType();		// return the target's type

private:

	static const zString m_Type;

	// ****************************************
	// PVFs that must be implemented by all instantiated derived classes 
public:

    const zString GetTargetType() const;    // return the target's type

    // Functions to add data points and normalise the histogram.

    virtual void AddDataPoint(double x);
    virtual void Normalise(double norm);
     
    // Implementation of the ISerialiseInclusiveRestartState interface
    // to allow this class to read/write data that can be modified
    // for restarts.
    
    virtual zInStream& Read(zInStream& is);
    virtual zOutStream& Write(zOutStream& is) const;

    // ****************************************
	// Public access functions
public:
   
    // Functions to initialize the histogram calculation, add data and normalize it
    // when required. We provide two overloaded functions to reset the histogram.
    // The first requires the number of bins and calculates their width from the data,
    // the second requires the bin width and calculates the number of bins from the data.

	void ResetFixedBinTotal(long bins);
	void ResetFixedBinWidth(double width);

	// ****************************************
	// Protected local functions
protected:



	// ****************************************
	// Implementation


	// ****************************************
	// Private functions
private:


	// ****************************************
	// Data members

protected:

private:

    bool  m_bFixedBinTotal;   // Flag showing how the histogram is to be initialised
    double m_OldValue;    // Previous value of time needed to calculate the interval

    long m_BinTotal;      // Local variable to allow the size of the array to be changed
};


#endif // !defined(AFX_TABINTIMEINTERVALS_H__DDA01256_8832_45EC_8067_19F09E9F5182__INCLUDED_)
