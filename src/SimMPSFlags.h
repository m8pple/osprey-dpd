// Header file to hold flags and options relating to the parallel code.
//
// This file must be included in all .h and .cpp files that refer to any of 
// the parallelised classes. Note that because it is used inside
// other header files, it should NOT contain any further #include statements.
//
// **********************************************************************
//
// Version History
// ***************
//
//	22/10/07    Baseline version. 
//              General flags relating to the parallel code are in SimDefs.h
//              where they are included in every compilation unit. This 
//              header file contains parallel-aware flags that only need to be
//              included in a limited number of files and, particularly, 
//              header files that need to distinguish between the standard code
//              and the parallelised version.
//  23/10/07    I added a flag (SimMPSSystemShape) showing the shape of the system. 
//              It can take one of the following values for the corresponding number
//              of processors:
//
//              SimMPSLinearShape = 1,   1, LZ
//              SimMPSPlanarShape = LX, LY, 1
//              SimMPSCubicShape  = LX, LY, LZ
//  30/10/07    I added more flags to allow different code features to be
//              toggled between serial and parallel implementations.
//           
// **********************************************************************

#define SimMPSEnabled	1
#define SimMPSDisabled	2

#ifndef SimMPS
#define SimMPS          SimMPSDisabled
#endif

#if SimMPS==SimMPSEnabled
	#define EnableParallelACN           SimMPSEnabled
	#define EnableParallelAggregates    SimMPSEnabled
	#define EnableParallelCommands      SimMPSEnabled
	#define EnableParallelEvents        SimMPSEnabled
	#define EnableParallelExperiment    SimMPSEnabled
	#define EnableParallelMonitor       SimMPSEnabled
	#define EnableParallelProcesses     SimMPSEnabled
	#define EnableParallelRestart       SimMPSEnabled
	#define EnableParallelSimBox        SimMPSEnabled
	#define EnableParallelTargets       SimMPSEnabled
#else	
	#define EnableParallelACN           SimMPSDisabled
	#define EnableParallelAggregates    SimMPSDisabled
	#define EnableParallelCommands      SimMPSDisabled
	#define EnableParallelEvents        SimMPSDisabled
	#define EnableParallelExperiment    SimMPSDisabled
	#define EnableParallelMonitor       SimMPSDisabled
	#define EnableParallelProcesses     SimMPSDisabled
	#define EnableParallelRestart       SimMPSDisabled
	#define EnableParallelSimBox        SimMPSDisabled
	#define EnableParallelTargets       SimMPSDisabled
#endif

