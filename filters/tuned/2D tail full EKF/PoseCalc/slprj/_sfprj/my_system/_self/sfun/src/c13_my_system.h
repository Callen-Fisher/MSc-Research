#ifndef __c13_my_system_h__
#define __c13_my_system_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_c13_ResolvedFunctionInfo
#define typedef_c13_ResolvedFunctionInfo

typedef struct {
  const char * context;
  const char * name;
  const char * dominantType;
  const char * resolved;
  uint32_T fileTimeLo;
  uint32_T fileTimeHi;
  uint32_T mFileTimeLo;
  uint32_T mFileTimeHi;
} c13_ResolvedFunctionInfo;

#endif                                 /*typedef_c13_ResolvedFunctionInfo*/

#ifndef typedef_SFc13_my_systemInstanceStruct
#define typedef_SFc13_my_systemInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c13_sfEvent;
  boolean_T c13_isStable;
  boolean_T c13_doneDoubleBufferReInit;
  uint8_T c13_is_active_c13_my_system;
  real_T c13_sampleTime;
  real_T c13_l1;
  real_T c13_l2;
  real_T c13_Cd1;
  real_T c13_Cd2;
  real_T c13_m1;
  real_T c13_m2;
  real_T c13_g;
  real_T c13_J1;
  real_T c13_J2;
  real_T c13_magVecX;
  real_T c13_magVecY;
  real_T c13_magVecZ;
  real_T c13_A1;
  real_T c13_A2;
  real_T c13_rho;
  real_T c13_initTh1;
  real_T c13_initPh1;
  real_T c13_initTh2;
  real_T c13_initPh2;
  real_T c13_states[12];
  boolean_T c13_states_not_empty;
  real_T c13_ddth2T;
  boolean_T c13_ddth2T_not_empty;
  real_T c13_ddph2T;
  boolean_T c13_ddph2T_not_empty;
  real_T c13_ddph1T;
  boolean_T c13_ddph1T_not_empty;
  real_T c13_covP[144];
  boolean_T c13_covP_not_empty;
} SFc13_my_systemInstanceStruct;

#endif                                 /*typedef_SFc13_my_systemInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c13_my_system_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c13_my_system_get_check_sum(mxArray *plhs[]);
extern void c13_my_system_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
