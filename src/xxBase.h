// xxBase.h: interface for the xxBase class.
//
// This class provides error message handling facilities to all derived classes.
//
// It also defines types that are used throughout the program to represent
// collections of beads, bonds, polymers, commands, etc. To enable this it 
// includes all the STL header files needed by its derived classes. Although this
// makes available container classes that are not needed by all the derived
// classes, it ensures that the typedefs are uniformly defined and used.
// 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_XXBASE_H__83457702_3C59_11D3_820E_0060088AD300__INCLUDED_)
#define AFX_XXBASE_H__83457702_3C59_11D3_820E_0060088AD300__INCLUDED_


#include "SimACNFlags.h"

#include <stddef.h>

// Forward declarations

#if EnableShadowSimBox == SimACNEnabled
class aeActiveBond;
class aeActiveCellNetwork;
class aeActivePolymer;
class aevActiveEvent;
class IACNAccessControl;
class IModifyActiveCellNetwork;
#endif

class aaScalarProfile;
class aaStressTensorPoint;
class aeCNTCell;
class CAbstractBead;
class CAnalysis;
class CAnalysisData;
class CAnalysisObject;
class CAnalysisTool;
class CBead;
class CBeadType;
class CBeadChargeWrapper;
class CBond;
class CBondType;
class CBondPair;
class CBondPairType;
class CCell;
class CCellProfile;
class CCellProfileSet;
class CCNTCell;
class CCommandTargetNode;
class CCurrentState;
class CDensityFieldState;
class CDensityState;
class CForceTarget;
class CGridObservable;
class CInitialStateData;
class CNanoparticle;
class CObservable;
class CPolymer;
class CPolymerType;
class CPolymerCrossLink;
class CRaft;
class CRow;
class CRowProfile;
class CRowProfileSet;
class CSlice;
class CStressGridCell;
class CTimeSeriesData;
class IActiveSimBox;
class IModifyProcess;
class IRegionAnalysis;
class mpuExtendedBond;
class mpuExtendedBondPair;
class mpuExtendedNanoparticle;
class mpuExtendedNPBond;
class mpuExtendedPolymer;
class mpuGhostBead;
class mpmMessage;
class pmSendGhostBeadCoords;
class taEventSourceDecorator;
class taEventAnalysisDecorator;
class tguArgumentType;
class tguCommandGroup;
class xxCommand;
class xxCommandObject;
class xxEvent;
class xxEventObject;
class xxMessage;
class xxProcess;
class xxProcessObject;


// **********************************************************************
// The following statements are platform-dependent but to avoid including the
// SimDefs.h header file here (where it would be included in every header file)
// we manually change the platform name. To distinguish the platform from the
// one in SimDefs.h I prefix it with xxBase.
// **********************************************************************

#define xxBaseMFC			1
#define xxBaseCONSOLE		2
#define xxBaseDECALPHA		3
#define xxBaseSGICC			4
#define xxBaseCRAYJ90		5
#define xxBaseBORLAND6		6
#define xxBaseI7XEON		7
#define xxBaseI7ITANIUM		8
#define xxBaseGCC 			9
#define xxBaseCW55MAC		10
#define xxBaseXCMAC			11
#define xxBaseCRAYHERMIT    12
#define xxBaseNEWPLATFORM2 	13
#define xxBaseNEWPLATFORM3 	14
#define xxBaseNEWPLATFORM4 	15
#define xxBaseNEWPLATFORM5 	16
#define xxBaseNEWPLATFORM6 	17
#define xxBaseNEWPLATFORM7 	18
#define xxBaseNEWPLATFORM8 	19
#define xxBaseNEWPLATFORM9 	20
#define xxBaseNEWPLATFORM10 21

// Select appropriate definitions for the current platform.

#define xxBasePlatform xxBaseXCMAC


// ****************************************
#if xxBasePlatform == xxBaseCONSOLE	// PC console application

	#include <string>
	#include <list>
	#include <vector>
	#include <map>
	#include <iostream>
	#include <fstream>
	#include <sstream>
	#include <iomanip>
	
	// **********************************************************************
	// Define typedefs to allow each platform to access the STL containers
	// and other standard libraries using std or not as required.
	// Typedefs beginning with z are used throughout the code. But macros
	// that begin with xxBase are only used in this header file to simplify
	// the replacement of the std name.

	typedef std::istream				zInStream;
	typedef std::ostream				zOutStream;
	typedef std::ifstream				zInFileStream;
	typedef std::ofstream				zOutFileStream;
	typedef std::fstream				zFileStream;
	typedef std::istringstream			zInStringStream;
	typedef std::ostringstream			zOutStringStream;
	typedef std::ios					zIos;

	// Define macros to wrap the std namespace name for different platforms

	#define zEndl						std::endl
	#define zFlush						std::flush
	#define zLeft						std::left
	#define xxBaselist					std::list
	#define xxBasemap					std::map
	#define xxBasemultimap				std::multimap
	#define xxBasestring				std::basic_string<char>

	// **********************************************************************
	// class xxBasevector
	// ******************
	//
	// Local vector class to provide the member at() functions for those platforms that
	// do not have it in their C++ language. Both non-const and const versions are
	// needed.
	// I have to provide a default constructor with no arguments and one which sets each
	// element of the vector to a specified initial value to accomodate the various
	// cases in the code. In particular, see Builder.cpp where the initial value form is
	// required. A copy constructor is also needed, but I hope that the assignment and 
	// comparison operators will default to the base class ones.

	template <class T>
	class xxBasevector : public std::vector<T>
	{
	public:
		xxBasevector() : std::vector<T>(){}
		xxBasevector(size_t num, const T& val=T()) : std::vector<T>(num,val){}
		xxBasevector(const xxBasevector& oldVec) : std::vector<T>(oldVec){}
		~xxBasevector(){}

//		inline T& at(size_t i)             {return (*(begin() + i)); }
//		inline const T& at(size_t i) const {return (*(begin() + i)); }
	};

	// The following pair functions have to take different forms on the Cray platform

	typedef std::pair<long, long>					zPairLL;
	typedef std::pair<long, double>					zPairLD;
	typedef std::pair<long, xxBasestring>			zPairLS;
	typedef std::pair<xxBasestring, long>			zPairSL;
	typedef std::pair<xxBasestring, double>			zPairSD;
	typedef std::pair<xxBasestring, xxBasestring>	zPairSS;
	typedef std::pair<long, CBeadType*>				zPairLongBeadType;
	typedef std::pair<long, CBondType*>				zPairLongBondType;
	typedef std::pair<long, CBondPairType*>			zPairLongBondPairType;
	typedef std::pair<long, CPolymerType*>			zPairLongPolymerType;
	typedef std::pair<long, CAbstractBead*>			zPairLongBead;
	typedef std::pair<long, CBond*>					zPairLongBond;
	typedef std::pair<long, CBondPair*>				zPairLongBondPair;
	typedef std::pair<long, CPolymer*>				zPairLongPolymer;
	typedef std::pair<long, mpuGhostBead*>	        zPairLongGhostBead;
	typedef std::pair<long, mpuExtendedBond*>	    zPairLongExtendedBond;
	typedef std::pair<long, mpuExtendedNPBond*>	    zPairLongExtendedNPBond;
	typedef std::pair<long, mpuExtendedPolymer*>	zPairLongExtendedPolymer;
	typedef std::pair<long, CCommandTargetNode*>	zPairLongTarget;
	typedef std::pair<long, CNanoparticle*>	        zPairLongNanoparticle;

	typedef std::pair<xxBasestring, tguArgumentType*>			    zPairSVAT;
	typedef std::pair<xxBasestring, CCommandTargetNode*>			zPairST;
	typedef std::pair<xxBasestring, IModifyProcess*>				zPairSMP;
	typedef std::pair<xxBasestring, taEventSourceDecorator*>		zPairSESD;
	typedef std::pair<xxBasestring, taEventAnalysisDecorator*>		zPairSEAD;

#if EnableShadowSimBox == SimACNEnabled
	typedef std::pair<xxBasestring, aeActiveCellNetwork*>		    zPairSACN;
	typedef std::pair<xxBasestring, IModifyActiveCellNetwork*>		zPairSMACN;
	typedef std::pair<xxBasestring, IACNAccessControl*>				zPairSAACN;

    typedef aeActiveCellNetwork* (*CreateACNCallback)(IActiveSimBox* const pShadow, long xn, long yn, long zn,
					                                  double xw, double yw, double zw);
	typedef aevActiveEvent*		(*CreateActiveEventCallback)();
#endif

// ****************************************
#elif xxBasePlatform == xxBaseXCMAC	// Mac OS X/XCode compiler console app

    #include <cstdint>     // Needed for uint32_t and uint64_t types
    #include <string.h>    // Needed for parallel compilation using mpicxx
    #include <string>
	#include <list>
	#include <vector>
	#include <map>
	#include <iostream>
	#include <fstream>
	#include <sstream>
	#include <iomanip>
	
	// **********************************************************************
	// Define typedefs to allow each platform to access the STL containers
	// and other standard libraries using std or not as required.
	// Typedefs beginning with z are used throughout the code. But macros
	// that begin with xxBase are only used in this header file to simplify
	// the replacement of the std name.

	typedef std::istream				zInStream;
	typedef std::ostream				zOutStream;
	typedef std::ifstream				zInFileStream;
	typedef std::ofstream				zOutFileStream;
	typedef std::fstream				zFileStream;
	typedef std::istringstream			zInStringStream;
	typedef std::ostringstream			zOutStringStream;
	typedef std::ios					zIos;

	// Define macros to wrap the std namespace name for different platforms

	#define zEndl						std::endl
	#define zFlush						std::flush
	#define zLeft						std::left
	#define xxBaselist					std::list
	#define xxBasemap					std::map
	#define xxBasemultimap				std::multimap
	#define xxBasestring				std::basic_string<char>

	// **********************************************************************
	// class xxBasevector
	// ******************
	//
	// Local vector class to provide the member at() functions for those platforms that
	// do not have it in their C++ language. Both non-const and const versions are
	// needed.
	// I have to provide a default constructor with no arguments and one which sets each
	// element of the vector to a specified initial value to accomodate the various
	// cases in the code. In particular, see Builder.cpp where the initial value form is
	// required. A copy constructor is also needed, but I hope that the assignment and 
	// comparison operators will default to the base class ones.

	template <class T>
	class xxBasevector : public std::vector<T>
	{
	public:
		xxBasevector() : std::vector<T>(){}
		xxBasevector(size_t num, const T& val=T()) : std::vector<T>(num,val){}
		xxBasevector(const xxBasevector& oldVec) : std::vector<T>(oldVec){}
		~xxBasevector(){}

//		inline T& at(size_t i)             {return (*(begin() + i)); }
//		inline const T& at(size_t i) const {return (*(begin() + i)); }
	};

	// The following pair functions have to take different forms on the Cray platform

	typedef std::pair<long, long>					zPairLL;
	typedef std::pair<long, double>					zPairLD;
	typedef std::pair<long, xxBasestring>			zPairLS;
	typedef std::pair<xxBasestring, long>			zPairSL;
	typedef std::pair<xxBasestring, double>			zPairSD;
	typedef std::pair<xxBasestring, xxBasestring>	zPairSS;
	typedef std::pair<long, CBeadType*>				zPairLongBeadType;
	typedef std::pair<long, CBondType*>				zPairLongBondType;
	typedef std::pair<long, CBondPairType*>			zPairLongBondPairType;
	typedef std::pair<long, CPolymerType*>			zPairLongPolymerType;
	typedef std::pair<long, CAbstractBead*>			zPairLongBead;
	typedef std::pair<long, CBond*>					zPairLongBond;
	typedef std::pair<long, CBondPair*>				zPairLongBondPair;
	typedef std::pair<long, CPolymer*>				zPairLongPolymer;
	typedef std::pair<long, mpuGhostBead*>	        zPairLongGhostBead;
	typedef std::pair<long, mpuExtendedBond*>	    zPairLongExtendedBond;
	typedef std::pair<long, mpuExtendedNPBond*>	    zPairLongExtendedNPBond;
	typedef std::pair<long, mpuExtendedPolymer*>	zPairLongExtendedPolymer;
	typedef std::pair<long, CCommandTargetNode*>	zPairLongTarget;
	typedef std::pair<long, CNanoparticle*>	        zPairLongNanoparticle;

	typedef std::pair<xxBasestring, tguArgumentType*>			    zPairSVAT;
	typedef std::pair<xxBasestring, CCommandTargetNode*>			zPairST;
	typedef std::pair<xxBasestring, IModifyProcess*>				zPairSMP;
	typedef std::pair<xxBasestring, taEventSourceDecorator*>		zPairSESD;
	typedef std::pair<xxBasestring, taEventAnalysisDecorator*>		zPairSEAD;

#if EnableShadowSimBox == SimACNEnabled
	typedef std::pair<xxBasestring, aeActiveCellNetwork*>		    zPairSACN;
	typedef std::pair<xxBasestring, IModifyActiveCellNetwork*>		zPairSMACN;
	typedef std::pair<xxBasestring, IACNAccessControl*>				zPairSAACN;
	typedef aeActiveCellNetwork* (*CreateACNCallback)(IActiveSimBox* const pShadow, long xn, long yn, long zn,
					                                  double xw, double yw, double zw);
	typedef aevActiveEvent*		(*CreateActiveEventCallback)();
#endif

// ****************************************
#elif xxBasePlatform == xxBaseCRAYHERMIT	// hermit1.hww.de Cray

	#include <string.h>
	#include <list>
	#include <vector>
	#include <map>
	#include <iostream>
	#include <fstream>
	#include <sstream>
	#include <iomanip>
	
	// **********************************************************************
	// Define typedefs to allow each platform to access the STL containers
	// and other standard libraries using std or not as required.
	// Typedefs beginning with z are used throughout the code. But macros
	// that begin with xxBase are only used in this header file to simplify
	// the replacement of the std name.

	typedef std::istream				zInStream;
	typedef std::ostream				zOutStream;
	typedef std::ifstream				zInFileStream;
	typedef std::ofstream				zOutFileStream;
	typedef std::fstream				zFileStream;
	typedef std::istringstream			zInStringStream;
	typedef std::ostringstream			zOutStringStream;
	typedef std::ios					zIos;

	// Define macros to wrap the std namespace name for different platforms

	#define zEndl						std::endl
	#define zFlush						std::flush
	#define zLeft						std::left
	#define xxBaselist					std::list
	#define xxBasemap					std::map
	#define xxBasemultimap				std::multimap
	#define xxBasestring				std::basic_string<char>

	// **********************************************************************
	// class xxBasevector
	// ******************
	//
	// Local vector class to provide the member at() functions for those platforms that
	// do not have it in their C++ language. Both non-const and const versions are
	// needed.
	// I have to provide a default constructor with no arguments and one which sets each
	// element of the vector to a specified initial value to accomodate the various
	// cases in the code. In particular, see Builder.cpp where the initial value form is
	// required. A copy constructor is also needed, but I hope that the assignment and 
	// comparison operators will default to the base class ones.

	template <class T>
	class xxBasevector : public std::vector<T>
	{
	public:
		xxBasevector() : std::vector<T>(){}
		xxBasevector(size_t num, const T& val=T()) : std::vector<T>(num,val){}
		xxBasevector(const xxBasevector& oldVec) : std::vector<T>(oldVec){}
		~xxBasevector(){}

//		inline T& at(size_t i)             {return (*(begin() + i)); }
//		inline const T& at(size_t i) const {return (*(begin() + i)); }
	};

	// The following pair functions have to take different forms on the Cray platform

	typedef std::pair<long, long>					zPairLL;
	typedef std::pair<long, double>					zPairLD;
	typedef std::pair<long, xxBasestring>			zPairLS;
	typedef std::pair<xxBasestring, long>			zPairSL;
	typedef std::pair<xxBasestring, double>			zPairSD;
	typedef std::pair<xxBasestring, xxBasestring>	zPairSS;
	typedef std::pair<long, CBeadType*>				zPairLongBeadType;
	typedef std::pair<long, CBondType*>				zPairLongBondType;
	typedef std::pair<long, CBondPairType*>			zPairLongBondPairType;
	typedef std::pair<long, CPolymerType*>			zPairLongPolymerType;
	typedef std::pair<long, CAbstractBead*>			zPairLongBead;
	typedef std::pair<long, CBond*>					zPairLongBond;
	typedef std::pair<long, CBondPair*>				zPairLongBondPair;
	typedef std::pair<long, CPolymer*>				zPairLongPolymer;
	typedef std::pair<long, mpuGhostBead*>	        zPairLongGhostBead;
	typedef std::pair<long, mpuExtendedBond*>	    zPairLongExtendedBond;
	typedef std::pair<long, mpuExtendedNPBond*>	    zPairLongExtendedNPBond;
	typedef std::pair<long, mpuExtendedPolymer*>	zPairLongExtendedPolymer;
	typedef std::pair<long, CCommandTargetNode*>	zPairLongTarget;
	typedef std::pair<long, CNanoparticle*>	        zPairLongNanoparticle;

	typedef std::pair<xxBasestring, tguArgumentType*>			    zPairSVAT;
	typedef std::pair<xxBasestring, CCommandTargetNode*>			zPairST;
	typedef std::pair<xxBasestring, IModifyProcess*>				zPairSMP;
	typedef std::pair<xxBasestring, taEventSourceDecorator*>		zPairSESD;
	typedef std::pair<xxBasestring, taEventAnalysisDecorator*>		zPairSEAD;

#if EnableShadowSimBox == SimACNEnabled
	typedef std::pair<xxBasestring, aeActiveCellNetwork*>		    zPairSACN;
	typedef std::pair<xxBasestring, IModifyActiveCellNetwork*>		zPairSMACN;
	typedef std::pair<xxBasestring, IACNAccessControl*>				zPairSAACN;
	typedef aeActiveCellNetwork* (*CreateACNCallback)(IActiveSimBox* const pShadow, long xn, long yn, long zn,
					                                  double xw, double yw, double zw);
	typedef aevActiveEvent*		(*CreateActiveEventCallback)();
#endif

// ****************************************
#elif xxBasePlatform == xxBaseCW55MAC	// Mac OS X/Code Warrior compiler

	#include <string>
	#include <list>
	#include <vector>
	#include <map>
	#include <iostream>
	#include <fstream>
	#include <sstream>
	#include <iomanip>
	
	// **********************************************************************
	// Define typedefs to allow each platform to access the STL containers
	// and other standard libraries using std or not as required.
	// Typedefs beginning with z are used throughout the code. But macros
	// that begin with xxBase are only used in this header file to simplify
	// the replacement of the std name.

	typedef std::istream				zInStream;
	typedef std::ostream				zOutStream;
	typedef std::ifstream				zInFileStream;
	typedef std::ofstream				zOutFileStream;
	typedef std::fstream				zFileStream;
	typedef std::istringstream			zInStringStream;
	typedef std::ostringstream			zOutStringStream;
	typedef std::ios					zIos;

	// Define macros to wrap the std namespace name for different platforms

	#define zEndl						std::endl
	#define zFlush						std::flush
	#define zLeft						std::left
	#define xxBaselist					std::list
	#define xxBasemap					std::map
	#define xxBasemultimap				std::multimap
	#define xxBasestring				std::basic_string<char>

	// **********************************************************************
	// class xxBasevector
	// ******************
	//
	// Local vector class to provide the member at() functions for those platforms that
	// do not have it in their C++ language. Both non-const and const versions are
	// needed.
	// I have to provide a default constructor with no arguments and one which sets each
	// element of the vector to a specified initial value to accomodate the various
	// cases in the code. In particular, see Builder.cpp where the initial value form is
	// required. A copy constructor is also needed, but I hope that the assignment and 
	// comparison operators will default to the base class ones.

	template <class T>
	class xxBasevector : public std::vector<T>
	{
	public:
		xxBasevector() : std::vector<T>(){}
		xxBasevector(size_t num, const T& val=T()) : std::vector<T>(num,val){}
		xxBasevector(const xxBasevector& oldVec) : std::vector<T>(oldVec){}
		~xxBasevector(){}

//		inline T& at(size_t i)             {return (*(begin() + i)); }
//		inline const T& at(size_t i) const {return (*(begin() + i)); }
	};

	// The following pair functions have to take different forms on the Cray platform

	typedef std::pair<long, long>					zPairLL;
	typedef std::pair<long, double>					zPairLD;
	typedef std::pair<long, xxBasestring>			zPairLS;
	typedef std::pair<xxBasestring, long>			zPairSL;
	typedef std::pair<xxBasestring, double>			zPairSD;
	typedef std::pair<xxBasestring, xxBasestring>	zPairSS;
	typedef std::pair<long, CBeadType*>				zPairLongBeadType;
	typedef std::pair<long, CBondType*>				zPairLongBondType;
	typedef std::pair<long, CBondPairType*>			zPairLongBondPairType;
	typedef std::pair<long, CPolymerType*>			zPairLongPolymerType;
	typedef std::pair<long, CAbstractBead*>			zPairLongBead;
	typedef std::pair<long, CBond*>					zPairLongBond;
	typedef std::pair<long, CBondPair*>				zPairLongBondPair;
	typedef std::pair<long, CPolymer*>				zPairLongPolymer;
	typedef std::pair<long, mpuGhostBead*>	        zPairLongGhostBead;
	typedef std::pair<long, mpuExtendedBond*>	    zPairLongExtendedBond;
	typedef std::pair<long, mpuExtendedNPBond*>	    zPairLongExtendedNPBond;
	typedef std::pair<long, mpuExtendedPolymer*>	zPairLongExtendedPolymer;
	typedef std::pair<long, CCommandTargetNode*>	zPairLongTarget;
	typedef std::pair<long, CNanoparticle*>	        zPairLongNanoparticle;

	typedef std::pair<xxBasestring, tguArgumentType*>			    zPairSVAT;
	typedef std::pair<xxBasestring, CCommandTargetNode*>			zPairST;
	typedef std::pair<xxBasestring, IModifyProcess*>				zPairSMP;
	typedef std::pair<xxBasestring, taEventSourceDecorator*>		zPairSESD;
	typedef std::pair<xxBasestring, taEventAnalysisDecorator*>		zPairSEAD;

#if EnableShadowSimBox == SimACNEnabled
	typedef std::pair<xxBasestring, aeActiveCellNetwork*>		    zPairSACN;
	typedef std::pair<xxBasestring, IModifyActiveCellNetwork*>		zPairSMACN;
	typedef std::pair<xxBasestring, IACNAccessControl*>				zPairSAACN;
	typedef aeActiveCellNetwork* (*CreateACNCallback)(IActiveSimBox* const pShadow, long xn, long yn, long zn,
					                                  double xw, double yw, double zw);
	typedef aevActiveEvent*		(*CreateActiveEventCallback)();
#endif

// ****************************************
#elif xxBasePlatform == xxBaseDECALPHA	// Dec alpha machines

	#include <string>
	#include <list>
	#include <vector>
	#include <map>
	#include <iostream>
	#include <fstream>
	#include <sstream>
	#include <iomanip>
	
	typedef std::istream				zInStream;
	typedef std::ostream				zOutStream;
	typedef std::ifstream				zInFileStream;
	typedef std::ofstream				zOutFileStream;
	typedef std::fstream				zFileStream;
	typedef std::istringstream			zInStringStream;
	typedef std::ostringstream			zOutStringStream;
	typedef std::ios					zIos;

	#define zEndl						std::endl
	#define zFlush						std::flush
	#define zLeft						std::left
	#define xxBaselist					std::list
	#define xxBasemap					std::map
	#define xxBasemultimap				std::multimap

	template <class T>
	class xxBasevector : public std::vector<T>
	{
	public:
		xxBasevector() : std::vector<T>(){}
		xxBasevector(size_t num, const T& val=T()) : std::vector<T>(num,val){}
		xxBasevector(const xxBasevector& oldVec) : std::vector<T>(oldVec){}
		~xxBasevector(){}
	};

	class xxBasestring : public std::basic_string<char>
	{
	public:
		xxBasestring() : std::basic_string<char>(){}
		xxBasestring(std::string bStr) : std::basic_string<char>(bStr){}
		xxBasestring(const char* cStr) : std::basic_string<char>(cStr){}
		xxBasestring(const xxBasestring& oldStr) : std::basic_string<char>(oldStr){}
		~xxBasestring(){}
	};

	// The following pair functions have to take different forms on the Cray platform

	typedef pair<long, long>						zPairLL;
	typedef pair<long, double>						zPairLD;
	typedef pair<long, xxBasestring>				zPairLS;
	typedef pair<xxBasestring, long>				zPairSL;
	typedef pair<xxBasestring, double>			    zPairSD;
	typedef pair<xxBasestring, xxBasestring>		zPairSS;
	typedef pair<long, CBeadType*>					zPairLongBeadType;
	typedef pair<long, CBondType*>					zPairLongBondType;
	typedef pair<long, CBondPairType*>				zPairLongBondPairType;
	typedef pair<long, CPolymerType*>				zPairLongPolymerType;
	typedef pair<long, CAbstractBead*>				zPairLongBead;
	typedef pair<long, CBond*>						zPairLongBond;
	typedef pair<long, CBondPair*>					zPairLongBondPair;
	typedef pair<long, CPolymer*>					zPairLongPolymer;
	typedef pair<long, mpuGhostBead*>	            zPairLongGhostBead;
	typedef pair<long, mpuExtendedBond*>	        zPairLongExtendedBond;
	typedef pair<long, mpuExtendedNPBond*>          zPairLongExtendedNPBond;
	typedef pair<long, mpuExtendedPolymer*>	        zPairLongExtendedPolymer;
	typedef pair<long, CCommandTargetNode*>			zPairLongTarget;
	typedef pair<long, CNanoparticle*>	            zPairLongNanoparticle;

	typedef pair<xxBasestring, tguArgumentType*>			    zPairSVAT;
	typedef pair<xxBasestring, CCommandTargetNode*>				zPairST;
	typedef pair<xxBasestring, IModifyProcess*>					zPairSMP;
	typedef pair<xxBasestring, taEventSourceDecorator*>	        zPairSESD;
	typedef pair<xxBasestring, taEventAnalysisDecorator*>       zPairSEAD;

#if EnableShadowSimBox == SimACNEnabled
	typedef pair<xxBasestring, aeActiveCellNetwork*>		    zPairSACN;
	typedef pair<xxBasestring, IModifyActiveCellNetwork*>		zPairSMACN;
	typedef pair<xxBasestring, IACNAccessControl*>				zPairSAACN;
	typedef aeActiveCellNetwork* (*CreateACNCallback)(IActiveSimBox* const pShadow, long xn, long yn, long zn,
					                                  double xw, double yw, double zw);
	typedef aevActiveEvent*		(*CreateActiveEventCallback)();
#endif


// ****************************************
#elif xxBasePlatform == xxBaseSGICC	// SGI platform with CC compiler

	#include <utility>
	#include <vector.h>
	#include <string>
	#include <list>
	#include <map>
	#include <iostream.h>
	#include <fstream.h>
	#include <strstream.h>
	#include <iomanip.h>

	typedef istream					zInStream;
	typedef ostream					zOutStream;
	typedef ifstream				zInFileStream;
	typedef ofstream				zOutFileStream;
	typedef fstream					zFileStream;
	typedef istrstream				zInStringStream;
	typedef ostrstream				zOutStringStream;
	typedef ios						zIos;

	#define zEndl						"\n"
	#define zFlush						" "
	#define zLeft						" "
	#define xxBaselist					std::list
	#define xxBasemap					std::map
	#define xxBasemultimap				std::multimap

	// The vector class is not in the std namespace, so don't prefix with std::.

	template <class T>
	class xxBasevector : public vector<T>
	{
	public:
		xxBasevector() : vector<T>(){}
		xxBasevector(size_t num, const T& val=T()) : vector<T>(num,val){}
		xxBasevector(const xxBasevector& oldVec) : vector<T>(oldVec){}
		~xxBasevector(){}

		inline T& at(size_t i)             {return (*(begin() + i)); }
		inline const T& at(size_t i) const {return (*(begin() + i)); }
	};

	// The basic_string<> class is in the std namespace, so we use the prefix std::

	class xxBasestring : public std::basic_string<char>
	{
	public:
		xxBasestring() : std::basic_string<char>(){}
		xxBasestring(std::basic_string<char> bStr) : std::basic_string<char>(bStr){}
		xxBasestring(const char* cStr) : std::basic_string<char>(cStr){}
		xxBasestring(const xxBasestring& oldStr) : std::basic_string<char>(oldStr){}
		~xxBasestring(){}

		friend inline bool operator!=(const xxBasestring& L, const xxBasestring& R) {return (!(L == R));}
		friend inline bool operator!=(const char *L, const xxBasestring& R)			{return (!(L == R));}
		friend inline bool operator!=(const xxBasestring& L, const char *R)			{return (!(L == R));}
	};

	// The following pair functions have to take different forms on the Cray platform

	typedef pair<long, long>						zPairLL;
	typedef pair<long, double>						zPairLD;
	typedef pair<long, xxBasestring>				zPairLS;
	typedef pair<xxBasestring, long>				zPairSL;
	typedef pair<xxBasestring, double>			    zPairSD;
	typedef pair<xxBasestring, xxBasestring>		zPairSS;
	typedef pair<long, CBeadType*>					zPairLongBeadType;
	typedef pair<long, CBondType*>					zPairLongBondType;
	typedef pair<long, CBondPairType*>				zPairLongBondPairType;
	typedef pair<long, CPolymerType*>				zPairLongPolymerType;
	typedef pair<long, CAbstractBead*>				zPairLongBead;
	typedef pair<long, CBond*>						zPairLongBond;
	typedef pair<long, CBondPair*>					zPairLongBondPair;
	typedef pair<long, CPolymer*>					zPairLongPolymer;
	typedef pair<long, mpuGhostBead*>	            zPairLongGhostBead;
	typedef pair<long, mpuExtendedBond*>	        zPairLongExtendedBond;
	typedef pair<long, mpuExtendedNPBond*>	        zPairLongExtendedNPBond;
	typedef pair<long, mpuExtendedPolymer*>	        zPairLongExtendedPolymer;
	typedef pair<long, CCommandTargetNode*>			zPairLongTarget;
	typedef pair<long, CNanoparticle*>	            zPairLongNanoparticle;

	typedef pair<xxBasestring, tguArgumentType*>			    zPairSVAT;
	typedef pair<xxBasestring, CCommandTargetNode*>				zPairST;
	typedef pair<xxBasestring, IModifyProcess*>					zPairSMP;
	typedef pair<xxBasestring, taEventSourceDecorator*>	        zPairSESD;
	typedef pair<xxBasestring, taEventAnalysisDecorator*>	    zPairSEAD;

#if EnableShadowSimBox == SimACNEnabled
	typedef pair<xxBasestring, aeActiveCellNetwork*>		    zPairSACN;
	typedef pair<xxBasestring, IModifyActiveCellNetwork*>		zPairSMACN;
	typedef pair<xxBasestring, IACNAccessControl*>				zPairSAACN;
	typedef aeActiveCellNetwork* (*CreateACNCallback)(IActiveSimBox* const pShadow, long xn, long yn, long zn,
					                                  double xw, double yw, double zw);
	typedef aevActiveEvent*		(*CreateActiveEventCallback)();
#endif


// ****************************************
#elif xxBasePlatform == xxBaseCRAYJ90		// Cray platform with CC compiler

	#include <stl_config.h>
	#include <mstring.h>
	#include <list.h>
	#include <vector.h>
	#include <map.h>
	#include <pair.h>
	#include <iostream.h>
	#include <fstream.h>
	#include <strstream.h>
	#include <iomanip.h>

	typedef istream					zInStream;
	typedef ostream					zOutStream;
	typedef ifstream				zInFileStream;
	typedef ofstream				zOutFileStream;
	typedef fstream					zFileStream;
	typedef istrstream				zInStringStream;
	typedef ostrstream				zOutStringStream;
	typedef ios						zIos;


	#define zEndl					"\n"
	#define zFlush					" "
	#define zLeft					" "
	#define xxBaselist				list
	#define xxBasemap				map
	#define xxBasemultimap			multimap


	template <class T>
	class xxBasevector : public vector<T>
	{
	public:
		xxBasevector() : vector<T>(){}
		xxBasevector(size_t num, const T& val=T()) : vector<T>(num,val){}
		xxBasevector(const xxBasevector& oldVec) : vector<T>(oldVec){}
		~xxBasevector(){}

		inline T& at(size_t i)             {return (*(begin() + i)); }
		inline const T& at(size_t i) const {return (*(begin() + i)); }
	};


	class xxBasestring : public basic_string<char>
	{
	public:
		xxBasestring() : basic_string<char>(){}
		xxBasestring(basic_string<char> bStr) : basic_string<char>(bStr){}
		xxBasestring(const char* cStr) : basic_string<char>(cStr){}
		xxBasestring(const xxBasestring& oldStr) : basic_string<char>(oldStr){}
		~xxBasestring(){}

		inline const char& at(size_t i) const {return (*(begin()+i));}

//		friend inline bool operator!=(const xxBasestring& L, const xxBasestring& R) {return (!(L == R));}
//		friend inline bool operator!=(const char *L, const xxBasestring& R)		    {return (!(L == R));}
//		friend inline bool operator!=(const xxBasestring& L, const char *R)		    {return (!(L == R));}

	};


	// The following pair functions have to take different forms on the Cray platform

	typedef pair<const long, long>					zPairLL;
	typedef pair<long, double>						zPairLD;
	typedef pair<const long, xxBasestring>			zPairLS;
	typedef pair<const xxBasestring, long>			zPairSL;
	typedef pair<xxBasestring, double>			    zPairSD;
	typedef pair<xxBasestring, xxBasestring>		zPairSS;
	typedef pair<const long, CBeadType*>			zPairLongBeadType;
	typedef pair<const long, CBondType*>			zPairLongBondType;
	typedef pair<const long, CBondPairType*>		zPairLongBondPairType;
	typedef pair<const long, CPolymerType*>			zPairLongPolymerType;
	typedef pair<const long, CAbstractBead*>		zPairLongBead;
	typedef pair<const long, CBond*>				zPairLongBond;
	typedef pair<const long, CBondPair*>			zPairLongBondPair;
	typedef pair<const long, CPolymer*>				zPairLongPolymer;
	typedef pair<long, mpuGhostBead*>	            zPairLongGhostBead;
	typedef pair<long, mpuExtendedBond*>	        zPairLongExtendedBond;
	typedef pair<long, mpuExtendedNPBond*>	        zPairLongExtendedNPBond;
	typedef pair<long, mpuExtendedPolymer*>	        zPairLongExtendedPolymer;
	typedef pair<long, CCommandTargetNode*>			zPairLongTarget;
	typedef pair<long, CNanoparticle*>	            zPairLongNanoparticle;

	typedef pair<xxBasestring, tguArgumentType*>			    zPairSVAT;
	typedef pair<xxBasestring, CCommandTargetNode*>				zPairST;
	typedef pair<xxBasestring, IModifyProcess*>					zPairSMP;
	typedef pair<xxBasestring, taEventSourceDecorator*>	        zPairSESD;
	typedef pair<xxBasestring, taEventAnalysisDecorator*>	    zPairSEAD;

#if EnableShadowSimBox == SimACNEnabled
	typedef pair<xxBasestring, aeActiveCellNetwork*>		    zPairSACN;
	typedef pair<xxBasestring, IModifyActiveCellNetwork*>		zPairSMACN;
	typedef pair<xxBasestring, IACNAccessControl*>				zPairSAACN;
	typedef aeActiveCellNetwork* (*CreateACNCallback)(IActiveSimBox* const pShadow, long xn, long yn, long zn,
					                                  double xw, double yw, double zw);
	typedef aevActiveEvent*		(*CreateActiveEventCallback)();
#endif

// ****************************************
#elif xxBasePlatform == xxBaseBORLAND6	// Borland 6 compiler in Windows application

	#include <string>
	#include <list>
	#include <vector>
	#include <map>
	#include <iostream>
	#include <fstream>
	#include <sstream>
	#include <iomanip>
	
	// **********************************************************************
	// Define typedefs to allow each platform to access the STL containers
	// and other standard libraries using std or not as required.
	// Typedefs beginning with z are used throughout the code. But macros
	// that begin with xxBase are only used in this header file to simplify
	// the replacement of the std name.

	typedef std::istream				zInStream;
	typedef std::ostream				zOutStream;
	typedef std::ifstream				zInFileStream;
	typedef std::ofstream				zOutFileStream;
	typedef std::fstream				zFileStream;
	typedef std::istringstream			zInStringStream;
	typedef std::ostringstream			zOutStringStream;
	typedef std::ios					zIos;

	// Define macros to wrap the std namespace names for different platforms

	#define zEndl						std::endl
	#define zFlush						std::flush
	#define zLeft						std::left
	#define xxBaselist					std::list
	#define xxBasemap					std::map
	#define xxBasemultimap				std::multimap
	#define xxBasestring			    std::basic_string<char>
    #define xxBasevector                std::vector

	// The following pair functions have to take different forms on the Cray platform

	typedef std::pair<long, long>					zPairLL;
	typedef std::pair<long, double>					zPairLD;
	typedef std::pair<long, xxBasestring>		    zPairLS;
	typedef std::pair<xxBasestring, long>		    zPairSL;
	typedef std::pair<xxBasestring, double>			zPairSD;
	typedef std::pair<xxBasestring, xxBasestring>	zPairSS;
	typedef std::pair<long, CBeadType*>				zPairLongBeadType;
	typedef std::pair<long, CBondType*>				zPairLongBondType;
	typedef std::pair<long, CBondPairType*>		    zPairLongBondPairType;
	typedef std::pair<long, CPolymerType*>		    zPairLongPolymerType;
	typedef std::pair<long, CAbstractBead*>	        zPairLongBead;
	typedef std::pair<long, CBond*>					zPairLongBond;
	typedef std::pair<long, CBondPair*>				zPairLongBondPair;
	typedef std::pair<long, CPolymer*>				zPairLongPolymer;
	typedef std::pair<long, mpuGhostBead*>	        zPairLongGhostBead;
	typedef std::pair<long, mpuExtendedBond*>	    zPairLongExtendedBond;
	typedef std::pair<long, mpuExtendedNPBond*>	    zPairLongExtendedNPBond;
	typedef std::pair<long, mpuExtendedPolymer*>	zPairLongExtendedPolymer;
	typedef std::pair<long, CCommandTargetNode*>	zPairLongTarget;
	typedef pair<long, CNanoparticle*>	            zPairLongNanoparticle;

	typedef std::pair<xxBasestring, tguArgumentType*>			    zPairSVAT;
	typedef std::pair<xxBasestring, CCommandTargetNode*>		    zPairST;
	typedef std::pair<xxBasestring, IModifyProcess*>			    zPairSMP;
	typedef std::pair<xxBasestring, taEventSourceDecorator*>	    zPairSESD;
	typedef std::pair<xxBasestring, taEventAnalysisDecorator*>	    zPairSEAD;

#if EnableShadowSimBox == SimACNEnabled
	typedef std::pair<xxBasestring, aeActiveCellNetwork*>		    zPairSACN;
	typedef std::pair<xxBasestring, IModifyActiveCellNetwork*>		zPairSMACN;
	typedef std::pair<xxBasestring, IACNAccessControl*>				zPairSAACN;
	typedef aeActiveCellNetwork* (*CreateACNCallback)(IActiveSimBox* const pShadow, long xn, long yn, long zn,
					                                  double xw, double yw, double zw);
	typedef aevActiveEvent*		(*CreateActiveEventCallback)();
#endif


#endif

// Typedef for the factory classes that are used to create the commands,
// events, processes, etc.
// These are called: acfCommandFactory, etc.

typedef CAnalysisData*		   (*CreateAnalysisCallback)();
typedef CInitialStateData* (*CreateInitialStateCallback)();
typedef mpmMessage*		(*CreateParallelMessageCallback)();
typedef xxCommand*				(*CreateCommandCallback)(long executionTime);
typedef xxEvent*				  (*CreateEventCallback)();
typedef xxProcess*				(*CreateProcessCallback)();


// **********************************************************************
// Typedefs to be propagated to the derived classes and others by including this 
// header file. Types whose names begin with z are used to name containers of
// the standard types (long, double, etc) while others have descriptive names.

// STL string containers. The Cray platform only has a const_iterator.

typedef xxBasestring									zString;
typedef xxBasestring::const_iterator					czStringIterator;

#if xxBasePlatform != xxBaseCRAYJ90
	typedef xxBasestring::iterator						zStringIterator;
#endif


// STL vector containers

typedef xxBasevector<bool>								zBoolVector;
typedef xxBasevector<bool>::iterator					zBoolVectorIterator;
typedef xxBasevector<bool>::const_iterator				czBoolVectorIterator;

typedef xxBasevector<long>								zLongVector;
typedef xxBasevector<long>::iterator					zLongVectorIterator;
typedef xxBasevector<long>::const_iterator				czLongVectorIterator;

typedef xxBasevector<double>							zDoubleVector;
typedef xxBasevector<double>::iterator					zDoubleVectorIterator;
typedef xxBasevector<double>::const_iterator			czDoubleVectorIterator;

typedef	xxBasevector<xxBasevector<long> >				zArray2dLong;
typedef	xxBasevector<xxBasevector<double> >				zArray2dDouble;
typedef	xxBasevector<xxBasevector<zString> >	        zArray2dString;


typedef xxBasevector<zString>							StringSequence;
typedef xxBasevector<zString>::iterator					StringSequenceIterator;
typedef xxBasevector<zString>::const_iterator			cStringSequenceIterator;

#if EnableShadowSimBox == SimACNEnabled
typedef xxBasevector<aeCNTCell*>						aeCNTCellVector;
typedef xxBasevector<aeCNTCell*>::iterator				aeCNTCellIterator;
typedef xxBasevector<aeCNTCell*>::const_iterator		caeCNTCellIterator;
#endif

typedef xxBasevector<CBeadType*>						BeadTypeVector;
typedef xxBasevector<CAbstractBead*>					AbstractBeadVector;
typedef xxBasevector<CBead*>							        BeadVector;
typedef xxBasevector<CBeadChargeWrapper*>							ChargedBeadVector;
typedef xxBasevector<CBond*>							        BondVector;
typedef xxBasevector<CBondPair*>						        BondPairVector;
typedef xxBasevector<CPolymer*>							        PolymerVector;
typedef xxBasevector<CPolymerCrossLink*>					    PolymerCrossLinkVector;

typedef xxBasevector<mpuGhostBead*>					            GhostBeadVector;
typedef xxBasevector<mpuExtendedBond*>						    ExtendedBondVector;
typedef xxBasevector<mpuExtendedBondPair*>					    ExtendedBondPairVector;
typedef xxBasevector<mpuExtendedNPBond*>					    ExtendedNPBondVector;
typedef xxBasevector<mpuExtendedNanoparticle*>					ExtendedNanoparticleVector;

typedef xxBasevector<mpuExtendedPolymer*>					    ExtendedPolymerVector;
typedef xxBasevector<pmSendGhostBeadCoords*>					GhostBeadCoordsMsgVector;
typedef xxBasevector<CNanoparticle*>							NanoparticleSequence;

typedef xxBasevector<CCNTCell*>							        CNTCellVector;
typedef xxBasevector<CStressGridCell*>						    StressCellVector;
typedef xxBasevector<CObservable*>						        ObservableSequence;
typedef xxBasevector<CObservable*>*						        ObservableSequencePointer;
typedef xxBasevector<CGridObservable*>					        GridObservableSequence;
typedef xxBasevector<CCurrentState*>					        CurrentStateSequence;
typedef xxBasevector<CDensityState*>					        DensityStateSequence;
typedef xxBasevector<CTimeSeriesData*>							TimeSeriesSequence;
typedef xxBasevector<CAnalysis*>								AggregateSequence;
typedef xxBasevector<CAnalysisData*>							AggregateDataSequence;
typedef xxBasevector<CAnalysisTool*>							ToolSequence;
typedef xxBasevector<IRegionAnalysis*>							RegionSequence;

#if EnableShadowSimBox == SimACNEnabled
typedef xxBasevector<aeActiveBond*>							    ActiveBondSequence;
typedef xxBasevector<aeActiveBond*>::iterator					ActiveBondIterator;
typedef xxBasevector<aeActiveBond*>::const_iterator			    cActiveBondIterator;
typedef xxBasevector<aeActiveBond*>::reverse_iterator			rActiveBondIterator;
typedef xxBasevector<aeActiveBond*>::const_reverse_iterator	    crActiveBondIterator;

typedef xxBasevector<aeActivePolymer*>							ActivePolymerSequence;
typedef xxBasevector<aeActivePolymer*>::iterator			    ActivePolymerIterator;
typedef xxBasevector<aeActivePolymer*>::const_iterator			cActivePolymerIterator;
typedef xxBasevector<aeActivePolymer*>::reverse_iterator	    rActivePolymerIterator;
typedef xxBasevector<aeActivePolymer*>::const_reverse_iterator	crActivePolymerIterator;

typedef xxBasevector<aevActiveEvent*>							ActiveEventSequence;
typedef xxBasevector<aevActiveEvent*>::iterator					ActiveEventIterator;
typedef xxBasevector<aevActiveEvent*>::const_iterator			cActiveEventIterator;
typedef xxBasevector<aevActiveEvent*>::reverse_iterator			rActiveEventIterator;
typedef xxBasevector<aevActiveEvent*>::const_reverse_iterator	crActiveEventIterator;
#endif

typedef xxBasevector<xxProcess*>								ProcessSequence;
typedef xxBasevector<xxEvent*>									EventSequence;
typedef xxBasevector<xxMessage*>								MessageSequence;
typedef xxBasevector<const CAnalysisObject*>					AnalysisObjectSequence;
typedef xxBasevector<const xxProcessObject*>					ProcessObjectSequence;
typedef xxBasevector<const xxEventObject*>						EventObjectSequence;

typedef xxBasevector<const xxCommandObject*>					CommandObjectSequence;
typedef xxBasevector<const xxCommandObject*>::iterator			CommandObjectIterator;
typedef xxBasevector<const xxCommandObject*>::const_iterator	cCommandObjectIterator;

typedef xxBasevector<CSlice*>							SliceVector;
typedef xxBasevector<CRow*>								RowVector;
typedef xxBasevector<CCell*>							CellVector;
typedef xxBasevector<CRowProfile*>						RowProfileVector;
typedef xxBasevector<CCellProfile*>						CellProfileVector;
typedef xxBasevector<CCellProfileSet*>					CellProfileSetVector;

typedef xxBasevector<aaScalarProfile*>					ScalarProfileSequence;
typedef xxBasevector<aaStressTensorPoint*>				StressVector;

typedef xxBasevector<CRaft*>							RaftSequence;

typedef xxBasevector<tguArgumentType*>					       tguArgumentSequence;

typedef xxBasemap<zString, CreateAnalysisCallback>							StringAnalysisMap;
typedef xxBasemap<zString, CreateInitialStateCallback>						StringInitialStateMap;
typedef xxBasemap<zString, CreateCommandCallback>							StringCommandMap;
typedef xxBasemap<zString, CreateEventCallback>								StringEventMap;
typedef xxBasemap<zString, CreateParallelMessageCallback>					StringParallelMessageMap;
typedef xxBasemap<zString, CreateProcessCallback>							StringProcessMap;


// STL list containers

typedef xxBaselist<CAbstractBead*>								BeadList;
typedef xxBaselist<CBond*>										BondList;
typedef xxBaselist<CBondPair*>									BondPairList;
typedef xxBaselist<CPolymer*>									PolymerList;
typedef xxBaselist<CNanoparticle*>								NanoparticleList;
typedef xxBaselist<mpuGhostBead*>								GhostBeadList;
typedef xxBaselist<mpuExtendedBond*>							ExtendedBondList;
typedef xxBaselist<mpuExtendedBondPair*>						    ExtendedBondPairList;
typedef xxBaselist<mpuExtendedNanoparticle*>							ExtendedNanoparticleList;
typedef xxBaselist<mpuExtendedNPBond*>							    ExtendedNPBondList;
typedef xxBaselist<mpuExtendedPolymer*>							    ExtendedPolymerList;
typedef xxBaselist<pmSendGhostBeadCoords*>                          GhostBeadCoordsMsgList;
typedef xxBaselist<CBeadChargeWrapper*>							ChargedBeadList;
typedef xxBaselist<CCommandTargetNode*>							CommandTargetSequence;
typedef xxBaselist<CDensityFieldState*>							DensityFieldSequence;


typedef xxBaselist<CForceTarget*>								ForceTargetSequence;
typedef xxBaselist<const xxCommand*>							CommandSequence;
typedef xxBaselist<tguCommandGroup*>						    CommandGroupSequence;

#if EnableShadowSimBox == SimACNEnabled
typedef xxBaselist<aeActiveBond*>								ActiveBondList;
typedef xxBaselist<aeActiveBond*>::iterator						ActiveBondListIterator;
typedef xxBaselist<aeActiveBond*>::const_iterator				cActiveBondListIterator;
typedef xxBaselist<aeActiveBond*>::reverse_iterator				rActiveBondListIterator;
typedef xxBaselist<aeActiveBond*>::const_reverse_iterator		crActiveBondListIterator;

typedef xxBaselist<aeActivePolymer*>							ActivePolymerList;
typedef xxBaselist<aeActivePolymer*>::iterator					ActivePolymerListIterator;
typedef xxBaselist<aeActivePolymer*>::const_iterator			cActivePolymerListIterator;
typedef xxBaselist<aeActivePolymer*>::reverse_iterator			rActivePolymerListIterator;
typedef xxBaselist<aeActivePolymer*>::const_reverse_iterator	crActivePolymerListIterator;

typedef xxBaselist<aevActiveEvent*>								ActiveEventList;
typedef xxBaselist<aevActiveEvent*>::iterator					ActiveEventListIterator;
typedef xxBaselist<aevActiveEvent*>::const_iterator				cActiveEventListIterator;
typedef xxBaselist<aevActiveEvent*>::reverse_iterator			rActiveEventListIterator;
typedef xxBaselist<aevActiveEvent*>::const_reverse_iterator		crActiveEventListIterator;
#endif

// STL map containers

typedef xxBasemap<long,long>							LongLongMap;
typedef xxBasemap<long,long>::iterator					LongLongIterator;
typedef xxBasemap<long,long>::const_iterator			cLongLongIterator;

typedef xxBasemap<long,double>							LongDoubleMap;
typedef xxBasemap<long,double>::iterator				LongDoubleIterator;
typedef xxBasemap<long,double>::const_iterator			cLongDoubleIterator;

typedef xxBasemap<long,zString>							LongStringMap;
typedef xxBasemap<long,zString>::iterator				LongStringIterator;
typedef xxBasemap<long,zString>::const_iterator			cLongStringIterator;

typedef xxBasemap<zString,long>							StringLongMap;
typedef xxBasemap<zString,long>::iterator				StringLongIterator;
typedef xxBasemap<zString,long>::const_iterator			cStringLongIterator;

typedef xxBasemap<zString,double>						StringDoubleMap;
typedef xxBasemap<zString,double>::iterator			    StringDoubleIterator;
typedef xxBasemap<zString,double>::const_iterator		cStringDoubleIterator;

typedef xxBasemap<zString,zString>						StringStringMap;
typedef xxBasemap<zString,zString>::iterator			StringStringIterator;
typedef xxBasemap<zString,zString>::const_iterator		cStringStringIterator;
typedef xxBasemap<long,CBeadType*>						LongBeadTypeMap;
typedef xxBasemap<long,CBondType*>						LongBondTypeMap;
typedef xxBasemap<long,CBondPairType*>					LongBondPairTypeMap;
typedef xxBasemap<long,CPolymerType*>					LongPolymerTypeMap;
typedef xxBasemap<long,CAbstractBead*>                  LongBeadMap;
typedef xxBasemap<long,CBond*>					        LongBondMap;
typedef xxBasemap<long,CBondPair*>					    LongBondPairMap;
typedef xxBasemap<long,CPolymer*>					     LongPolymerMap;
typedef xxBasemap<long,CNanoparticle*>					 LongNanoparticleMap;


typedef xxBasemap<long,mpuGhostBead*>				                        LongGhostBeadMap;
typedef xxBasemap<long,mpuExtendedBond*>				                    LongExtendedBondMap;
typedef xxBasemap<long,mpuExtendedNPBond*>				                    LongExtendedNPBondMap;
typedef xxBasemap<long,mpuExtendedNanoparticle*>				            LongExtendedNanoparticleMap;
typedef xxBasemap<long,mpuExtendedPolymer*>				                    LongExtendedPolymerMap;
typedef xxBasemap<long,CCommandTargetNode*>									LongTargetMap;
typedef xxBasemap<zString,CCommandTargetNode*>								StringTargetMap;
typedef xxBasemap<zString,tguArgumentType*>		                            StringArgumentTypeMap;



#if EnableShadowSimBox == SimACNEnabled
typedef xxBasemap<zString, CreateACNCallback>						        StringACNMap;
typedef xxBasemap<zString, CreateACNCallback>::iterator				        StringACNIterator;
typedef xxBasemap<zString, CreateACNCallback>::const_iterator		        cStringACNIterator;

typedef xxBasemap<zString, CreateActiveEventCallback>						StringActiveEventMap;
typedef xxBasemap<zString, CreateActiveEventCallback>::iterator				StringActiveEventIterator;
typedef xxBasemap<zString, CreateActiveEventCallback>::const_iterator		cStringActiveEventIterator;
#endif

typedef xxBasemap<zString, IModifyProcess*>									StringModifiableProcessMap;

#if EnableShadowSimBox == SimACNEnabled
typedef xxBasemap<zString, aeActiveCellNetwork*>							StringActiveACNMap;
typedef xxBasemap<zString, aeActiveCellNetwork*>::iterator					StringActiveACNIterator;
typedef xxBasemap<zString, aeActiveCellNetwork*>::const_iterator			cStringActiveACNIterator;

typedef xxBasemap<zString, IModifyActiveCellNetwork*>						StringModifiableACNMap;
typedef xxBasemap<zString, IModifyActiveCellNetwork*>::iterator				StringModifiableACNIterator;
typedef xxBasemap<zString, IModifyActiveCellNetwork*>::const_iterator		cStringModifiableACNIterator;

typedef xxBasemap<zString, IACNAccessControl*>								StringAccessibleACNMap;
typedef xxBasemap<zString, IACNAccessControl*>::iterator					StringAccessibleACNIterator;
typedef xxBasemap<zString, IACNAccessControl*>::const_iterator			    cStringAccessibleACNIterator;
#endif

// STL multimap containers

typedef xxBasemultimap<zString,long>			                             StringLongMMap;
//typedef xxBasemultimap<zString,long>::iterator					         StringLongMMIterator;
//typedef xxBasemultimap<zString,long>::const_iterator			             cStringLongMMIterator;
typedef xxBasemultimap<zString,double>			                             StringDoubleMMap;
//typedef xxBasemultimap<zString,double>::iterator					         StringDoubleMMIterator;
//typedef xxBasemultimap<zString,double>::const_iterator			         cStringDoubleMMIterator;
typedef xxBasemultimap<zString,zString>			                             StringStringMMap;
typedef xxBasemultimap<zString,zString>::iterator					         StringStringMMIterator;
typedef xxBasemultimap<zString,zString>::const_iterator			             cStringStringMMIterator;

typedef xxBasemultimap<long,long>			                                 LongLongMMap;
typedef xxBasemultimap<long,long>::iterator					                 LongLongMMIterator;
typedef xxBasemultimap<long,long>::const_iterator			                 cLongLongMMIterator;
typedef xxBasemultimap<long,double>			                                 LongDoubleMMap;
//typedef xxBasemultimap<long,double>::iterator					             LongDoubleMMIterator;
//typedef xxBasemultimap<long,double>::const_iterator			             cLongDoubleMMIterator;
typedef xxBasemultimap<long,zString>			                             LongStringMMap;
typedef xxBasemultimap<long,zString>::iterator					             LongStringMMIterator;
typedef xxBasemultimap<long,zString>::const_iterator			             cLongStringMMIterator;


typedef xxBasemultimap<long,CAbstractBead*>			                         LongBeadMMap;
typedef xxBasemultimap<long,CBond*>					                         LongBondMMap;
typedef xxBasemultimap<long,CBondPair*>				                         LongBondPairMMap;
typedef xxBasemultimap<long,CPolymer*>				                         LongPolymerMMap;
typedef xxBasemultimap<zString,taEventSourceDecorator*>	                     StringEventSourceMMap;
typedef xxBasemultimap<zString,taEventAnalysisDecorator*>	                 StringEventAnalysisMMap;

typedef xxBasemultimap<long,mpuGhostBead*>				                     LongGhostBeadMMap;
typedef xxBasemultimap<long,mpuExtendedBond*>				                 LongExtendedBondMMap;
typedef xxBasemultimap<long,mpuExtendedNPBond*>				                 LongExtendedNPBondMMap;




// **********************************************************************

class xxBase  
{
public:
	xxBase();				
	xxBase(const xxBase& oldBase);				
	virtual ~xxBase();


public:			// static constants used in initialisation routines

	static const double m_globalPI;
	static const double m_globalTwoPI;
	static const double m_globalFourPI;
	static const double m_globalPIByTwo;

    static const zString GetFilePrefix();

    static const zString GetAAPrefix();
    static const zString GetASPrefix();
    static const zString GetCDFPrefix();
    static const zString GetCHPrefix();
    static const zString GetCQPrefix();
    static const zString GetCSPrefix();
    static const zString GetDPPrefix();
    static const zString GetDSPrefix();
	static const zString GetEADPrefix();
    static const zString GetHSPrefix();
    static const zString GetISPrefix();
    static const zString GetLSPrefix();
    static const zString GetOSPrefix();
    static const zString GetPSPrefix();
    static const zString GetPSDFPrefix();
    static const zString GetRAPrefix();
    static const zString GetRSPrefix();
	static const zString GetTADSPrefix();

public:

	zString ToString(unsigned long x) const;
	zString ToString(long x) const;
	zString ToString(double x) const;

protected:
	virtual bool ErrorTrace(zString errStr) const;
	virtual void TraceEndl() const;
	virtual void Trace(zString msgStr) const;
	virtual void TraceStringNoEndl(zString msgStr) const;
	virtual void TraceIntNoEndl(long var1) const;
	virtual void TraceInt(zString msgStr, long var1) const;
	virtual void TraceInt2(zString msgStr, long var1, long var2) const;
	virtual void TraceInt3(zString msgStr, long var1, long var2, long var3) const;
	virtual void TraceDouble( zString msgStr, double var1) const;
	virtual void TraceDoubleNoEndl(double var1) const;
	virtual void TraceDouble2(zString msgStr, double var1, double var2) const;
	virtual void TraceVector(zString msgStr, double var1, double var2, double var3) const;


};

#endif // !defined(AFX_XXBASE_H__83457702_3C59_11D3_820E_0060088AD300__INCLUDED_)


