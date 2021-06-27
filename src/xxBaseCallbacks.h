// xxBaseCallbacks.h: interface for the xxBase Callbacks
//
// This defines types specifically related to callabacks
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_XXBASE_CALLBACKS_H)
#define AFX_XXBASE_CALLBACKS_

#include "xxBase.h"


typedef xxBasemap<zString, CreateAnalysisCallback>::iterator				StringAnalysisIterator;
typedef xxBasemap<zString, CreateAnalysisCallback>::const_iterator			cStringAnalysisIterator;

typedef xxBasemap<zString, CreateInitialStateCallback>::iterator			StringInitialStateIterator;
typedef xxBasemap<zString, CreateInitialStateCallback>::const_iterator		cStringInitialStateIterator;

typedef xxBasemap<zString, CreateCommandCallback>::iterator					StringCommandIterator;
typedef xxBasemap<zString, CreateCommandCallback>::const_iterator			cStringCommandIterator;

typedef xxBasemap<zString, CreateEventCallback>::iterator					StringEventIterator;
typedef xxBasemap<zString, CreateEventCallback>::const_iterator				cStringEventIterator;

typedef xxBasemap<zString, CreateParallelMessageCallback>::iterator			StringParallelMessageIterator;
typedef xxBasemap<zString, CreateParallelMessageCallback>::const_iterator	cStringParallelMessageIterator;

typedef xxBasemap<zString, CreateProcessCallback>::iterator					StringProcessIterator;
typedef xxBasemap<zString, CreateProcessCallback>::const_iterator			cStringProcessIterator;

#endif



