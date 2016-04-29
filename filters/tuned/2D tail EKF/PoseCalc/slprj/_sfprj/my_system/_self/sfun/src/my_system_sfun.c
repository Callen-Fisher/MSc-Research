/* Include files */

#include "my_system_sfun.h"
#include "my_system_sfun_debug_macros.h"
#include "c13_my_system.h"
#include "c14_my_system.h"
#include "c15_my_system.h"
#include "c16_my_system.h"
#include "c17_my_system.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _my_systemMachineNumber_;
real_T _sfTime_;

/* Function Declarations */

/* Function Definitions */
void my_system_initializer(void)
{
}

void my_system_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_my_system_method_dispatcher(SimStruct *simstructPtr, unsigned
  int chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==13) {
    c13_my_system_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==14) {
    c14_my_system_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==15) {
    c15_my_system_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==16) {
    c16_my_system_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==17) {
    c17_my_system_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

unsigned int sf_my_system_process_check_sum_call( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[20];
  if (nrhs<1 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the checksum */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"sf_get_check_sum"))
    return 0;
  plhs[0] = mxCreateDoubleMatrix( 1,4,mxREAL);
  if (nrhs>1 && mxIsChar(prhs[1])) {
    mxGetString(prhs[1], commandName,sizeof(commandName)/sizeof(char));
    commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
    if (!strcmp(commandName,"machine")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2468738197U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(587713135U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3012334349U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1562132118U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1465757251U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2244233252U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1649039863U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3073978159U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 13:
        {
          extern void sf_c13_my_system_get_check_sum(mxArray *plhs[]);
          sf_c13_my_system_get_check_sum(plhs);
          break;
        }

       case 14:
        {
          extern void sf_c14_my_system_get_check_sum(mxArray *plhs[]);
          sf_c14_my_system_get_check_sum(plhs);
          break;
        }

       case 15:
        {
          extern void sf_c15_my_system_get_check_sum(mxArray *plhs[]);
          sf_c15_my_system_get_check_sum(plhs);
          break;
        }

       case 16:
        {
          extern void sf_c16_my_system_get_check_sum(mxArray *plhs[]);
          sf_c16_my_system_get_check_sum(plhs);
          break;
        }

       case 17:
        {
          extern void sf_c17_my_system_get_check_sum(mxArray *plhs[]);
          sf_c17_my_system_get_check_sum(plhs);
          break;
        }

       default:
        ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0.0);
      }
    } else if (!strcmp(commandName,"target")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3564696471U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(678668628U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1090454852U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3896867807U);
    } else {
      return 0;
    }
  } else {
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(464494296U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1967889976U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3985102111U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1242855064U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_my_system_autoinheritance_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[32];
  char aiChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the autoinheritance_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_autoinheritance_info"))
    return 0;
  mxGetString(prhs[2], aiChksum,sizeof(aiChksum)/sizeof(char));
  aiChksum[(sizeof(aiChksum)/sizeof(char)-1)] = '\0';

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 13:
      {
        if (strcmp(aiChksum, "Oe29uawHUkZxW5LrsG2kTD") == 0) {
          extern mxArray *sf_c13_my_system_get_autoinheritance_info(void);
          plhs[0] = sf_c13_my_system_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 14:
      {
        if (strcmp(aiChksum, "J5raknxOvQuOHBNJbqIN2F") == 0) {
          extern mxArray *sf_c14_my_system_get_autoinheritance_info(void);
          plhs[0] = sf_c14_my_system_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 15:
      {
        if (strcmp(aiChksum, "U0WLiAL691ISwMrbvrJdKH") == 0) {
          extern mxArray *sf_c15_my_system_get_autoinheritance_info(void);
          plhs[0] = sf_c15_my_system_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 16:
      {
        if (strcmp(aiChksum, "U0WLiAL691ISwMrbvrJdKH") == 0) {
          extern mxArray *sf_c16_my_system_get_autoinheritance_info(void);
          plhs[0] = sf_c16_my_system_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 17:
      {
        if (strcmp(aiChksum, "mHwpEqptu722zW1gN1IX1B") == 0) {
          extern mxArray *sf_c17_my_system_get_autoinheritance_info(void);
          plhs[0] = sf_c17_my_system_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_my_system_get_eml_resolved_functions_info( int nlhs, mxArray *
  plhs[], int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[64];
  if (nrhs<2 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the get_eml_resolved_functions_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_eml_resolved_functions_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 13:
      {
        extern const mxArray *sf_c13_my_system_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c13_my_system_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 14:
      {
        extern const mxArray *sf_c14_my_system_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c14_my_system_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 15:
      {
        extern const mxArray *sf_c15_my_system_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c15_my_system_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 16:
      {
        extern const mxArray *sf_c16_my_system_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c16_my_system_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 17:
      {
        extern const mxArray *sf_c17_my_system_get_eml_resolved_functions_info
          (void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c17_my_system_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_my_system_third_party_uses_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the third_party_uses_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_third_party_uses_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 13:
      {
        if (strcmp(tpChksum, "VDsMTmgojjvP2vpsILfmHF") == 0) {
          extern mxArray *sf_c13_my_system_third_party_uses_info(void);
          plhs[0] = sf_c13_my_system_third_party_uses_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "ttnNd5JRgtHsnwXsxYXZ7") == 0) {
          extern mxArray *sf_c14_my_system_third_party_uses_info(void);
          plhs[0] = sf_c14_my_system_third_party_uses_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "3YRugGLz9U4vcTbxtI8cKG") == 0) {
          extern mxArray *sf_c15_my_system_third_party_uses_info(void);
          plhs[0] = sf_c15_my_system_third_party_uses_info();
          break;
        }
      }

     case 16:
      {
        if (strcmp(tpChksum, "3YRugGLz9U4vcTbxtI8cKG") == 0) {
          extern mxArray *sf_c16_my_system_third_party_uses_info(void);
          plhs[0] = sf_c16_my_system_third_party_uses_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "2ffaOmsPSm4ZfOuEwEkUBE") == 0) {
          extern mxArray *sf_c17_my_system_third_party_uses_info(void);
          plhs[0] = sf_c17_my_system_third_party_uses_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

void my_system_debug_initialize(struct SfDebugInstanceStruct* debugInstance)
{
  _my_systemMachineNumber_ = sf_debug_initialize_machine(debugInstance,
    "my_system","sfun",0,5,0,0,0);
  sf_debug_set_machine_event_thresholds(debugInstance,_my_systemMachineNumber_,0,
    0);
  sf_debug_set_machine_data_thresholds(debugInstance,_my_systemMachineNumber_,0);
}

void my_system_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_my_system_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("my_system",
      "my_system");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_my_system_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}
