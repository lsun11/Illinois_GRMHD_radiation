 /*@@
   @file      ExternalVariables.h
   @date      Wed May 22 02:32:10 2002
   @author    Ian Hawke
   @desc 
   The header file containing the local variables used across routines.
   These are the arrays containing GF indexes for all types of variables,
   and the number of each type of variable currently in use (the 
   parameters only give the maximum possible number).
   No function prototypes are defined in this file, so we do not protect
   it with an ifdef so that we can do inclusion within multiple routines
   in the same file.
   @enddesc 
   @version   $Header: /cactusdevcvs/CactusBase/MoL/src/ExternalVariables.h,v 1.4 2006/08/01 12:57:45 hawke Exp $
 @@*/


extern CCTK_INT *EvolvedVariableIndex;
extern CCTK_INT *RHSVariableIndex;
extern CCTK_INT *ConstrainedVariableIndex;
extern CCTK_INT *SandRVariableIndex;


extern CCTK_INT MoLNumEvolvedVariables;
extern CCTK_INT MoLNumConstrainedVariables;
extern CCTK_INT MoLNumSandRVariables;



extern CCTK_INT *EvolvedComplexVariableIndex;
extern CCTK_INT *RHSComplexVariableIndex;
extern CCTK_INT *ConstrainedComplexVariableIndex;
extern CCTK_INT *SandRComplexVariableIndex;


extern CCTK_INT MoLNumEvolvedComplexVariables;
extern CCTK_INT MoLNumConstrainedComplexVariables;
extern CCTK_INT MoLNumSandRComplexVariables;



extern CCTK_INT *EvolvedArrayVariableIndex;
extern CCTK_INT *RHSArrayVariableIndex;
extern CCTK_INT *ConstrainedArrayVariableIndex;
extern CCTK_INT *SandRArrayVariableIndex;


extern CCTK_INT MoLNumEvolvedArrayVariables;
extern CCTK_INT MoLNumConstrainedArrayVariables;
extern CCTK_INT MoLNumSandRArrayVariables;



extern CCTK_INT *EvolvedComplexArrayVariableIndex;
extern CCTK_INT *RHSComplexArrayVariableIndex;
extern CCTK_INT *ConstrainedComplexArrayVariableIndex;
extern CCTK_INT *SandRComplexArrayVariableIndex;


extern CCTK_INT MoLNumEvolvedComplexArrayVariables;
extern CCTK_INT MoLNumConstrainedComplexArrayVariables;
extern CCTK_INT MoLNumSandRComplexArrayVariables;


extern CCTK_INT ScheduleStatus;


extern CCTK_REAL *ArrayScratchSpace;
extern CCTK_INT *ArrayScratchSizes;
extern CCTK_INT CurrentArrayScratchSize;
