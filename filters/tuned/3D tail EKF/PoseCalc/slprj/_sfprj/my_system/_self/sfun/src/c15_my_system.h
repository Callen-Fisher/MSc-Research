#ifndef __c15_my_system_h__
#define __c15_my_system_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_c15_ResolvedFunctionInfo
#define typedef_c15_ResolvedFunctionInfo

typedef struct {
  const char * context;
  const char * name;
  const char * dominantType;
  const char * resolved;
  uint32_T fileTimeLo;
  uint32_T fileTimeHi;
  uint32_T mFileTimeLo;
  uint32_T mFileTimeHi;
} c15_ResolvedFunctionInfo;

#endif                                 /*typedef_c15_ResolvedFunctionInfo*/

#ifndef typedef_SFc15_my_systemInstanceStruct
#define typedef_SFc15_my_systemInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c15_sfEvent;
  boolean_T c15_isStable;
  boolean_T c15_doneDoubleBufferReInit;
  uint8_T c15_is_active_c15_my_system;
  real_T c15_sampleTime;
  real_T c15_states[3];
  boolean_T c15_states_not_empty;
} SFc15_my_systemInstanceStruct;

#endif                                 /*typedef_SFc15_my_systemInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c15_my_system_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c15_my_system_get_check_sum(mxArray *plhs[]);
extern void c15_my_system_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
