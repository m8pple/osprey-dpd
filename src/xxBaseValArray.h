// xxBaseValArray.h: interface for the xxBase class.
//
// Extracts ValArray specialisations from xxBase.h
// 
//////////////////////////////////////////////////////////////////////

#ifndef AFX_XXBASE_VALARRAY_H
#define AFX_XXBASE_VALARRAY_H

#include "xxBase.h"

// ****************************************
#if xxBasePlatform == xxBaseCONSOLE	// PC console application

	#include <valarray>

	// The following typedef uses valarrays where available or vectors where not

	typedef std::valarray<size_t>				zLongValArray;
	typedef std::valarray<double>				zDoubleValArray;

// ****************************************
#elif xxBasePlatform == xxBaseXCMAC	// Mac OS X/XCode compiler console app

    #include <valarray>

	// The following typedef uses valarrays where available or vectors where not

	typedef std::valarray<size_t>				zLongValArray;
	typedef std::valarray<double>				zDoubleValArray;

// ****************************************
#elif xxBasePlatform == xxBaseCRAYHERMIT	// hermit1.hww.de Cray

	#include <valarray>

	// The following typedef uses valarrays where available or vectors where not

	typedef std::valarray<size_t>				zLongValArray;
	typedef std::valarray<double>				zDoubleValArray;

// ****************************************
#elif xxBasePlatform == xxBaseCW55MAC	// Mac OS X/Code Warrior compiler

	#include <valarray>

	// The following typedef uses valarrays where available or vectors where not

	typedef std::valarray<size_t>				zLongValArray;
	typedef std::valarray<double>				zDoubleValArray;

// ****************************************
#elif xxBasePlatform == xxBaseDECALPHA	// Dec alpha machines

	#include <valarray>

	// The following typedef uses valarrays where available or vectors where not

	typedef std::valarray<size_t>			zLongValArray;
	typedef std::valarray<double>			zDoubleValArray;


// ****************************************
#elif xxBasePlatform == xxBaseSGICC	// SGI platform with CC compiler

	#include <valarray>

	// The following typedef uses valarrays where available or vectors where not

	typedef std::vector<long>					zLongValArray;
	typedef std::vector<double>					zDoubleValArray;


// ****************************************
#elif xxBasePlatform == xxBaseCRAYJ90		// Cray platform with CC compiler

	#include <vector>

	// The following typedef uses valarrays where available or vectors where not

	typedef vector<long>						zLongValArray;
	typedef vector<double>						zDoubleValArray;

// ****************************************
#elif xxBasePlatform == xxBaseBORLAND6	// Borland 6 compiler in Windows application

	#include <valarray>

	// The following typedef uses valarrays where available or vectors where not

	typedef std::valarray<size_t>				zLongValArray;
	typedef std::valarray<double>				zDoubleValArray;

#endif

#endif


