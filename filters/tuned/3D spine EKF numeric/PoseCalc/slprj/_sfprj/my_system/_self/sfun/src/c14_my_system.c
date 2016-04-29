/* Include files */

#include <stddef.h>
#include "blas.h"
#include "my_system_sfun.h"
#include "c14_my_system.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "my_system_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c14_debug_family_names[20] = { "accNoise", "gyroNoise",
  "magNoise1", "magNoise2", "magNoise3", "magNoise4", "maxA", "val1", "val2",
  "val3", "val4", "val5", "val6", "nargin", "nargout", "sampleTime", "J1", "J2",
  "Q", "R" };

/* Function Declarations */
static void initialize_c14_my_system(SFc14_my_systemInstanceStruct
  *chartInstance);
static void initialize_params_c14_my_system(SFc14_my_systemInstanceStruct
  *chartInstance);
static void enable_c14_my_system(SFc14_my_systemInstanceStruct *chartInstance);
static void disable_c14_my_system(SFc14_my_systemInstanceStruct *chartInstance);
static void c14_update_debugger_state_c14_my_system
  (SFc14_my_systemInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c14_my_system(SFc14_my_systemInstanceStruct *
  chartInstance);
static void set_sim_state_c14_my_system(SFc14_my_systemInstanceStruct
  *chartInstance, const mxArray *c14_st);
static void finalize_c14_my_system(SFc14_my_systemInstanceStruct *chartInstance);
static void sf_c14_my_system(SFc14_my_systemInstanceStruct *chartInstance);
static void initSimStructsc14_my_system(SFc14_my_systemInstanceStruct
  *chartInstance);
static void registerMessagesc14_my_system(SFc14_my_systemInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c14_machineNumber, uint32_T
  c14_chartNumber);
static const mxArray *c14_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData);
static void c14_emlrt_marshallIn(SFc14_my_systemInstanceStruct *chartInstance,
  const mxArray *c14_R, const char_T *c14_identifier, real_T c14_y[144]);
static void c14_b_emlrt_marshallIn(SFc14_my_systemInstanceStruct *chartInstance,
  const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId, real_T c14_y[144]);
static void c14_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData);
static const mxArray *c14_b_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData);
static real_T c14_c_emlrt_marshallIn(SFc14_my_systemInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId);
static void c14_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData);
static real_T c14_mpower(SFc14_my_systemInstanceStruct *chartInstance, real_T
  c14_a);
static void c14_eml_scalar_eg(SFc14_my_systemInstanceStruct *chartInstance);
static const mxArray *c14_c_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData);
static int32_T c14_d_emlrt_marshallIn(SFc14_my_systemInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId);
static void c14_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData);
static uint8_T c14_e_emlrt_marshallIn(SFc14_my_systemInstanceStruct
  *chartInstance, const mxArray *c14_b_is_active_c14_my_system, const char_T
  *c14_identifier);
static uint8_T c14_f_emlrt_marshallIn(SFc14_my_systemInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId);
static void init_dsm_address_info(SFc14_my_systemInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c14_my_system(SFc14_my_systemInstanceStruct
  *chartInstance)
{
  chartInstance->c14_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c14_is_active_c14_my_system = 0U;
}

static void initialize_params_c14_my_system(SFc14_my_systemInstanceStruct
  *chartInstance)
{
  real_T c14_d0;
  real_T c14_d1;
  real_T c14_d2;
  sf_set_error_prefix_string(
    "Error evaluating data 'sampleTime' in the parent workspace.\n");
  sf_mex_import_named("sampleTime", sf_mex_get_sfun_param(chartInstance->S, 2, 0),
                      &c14_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c14_sampleTime = c14_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'J1' in the parent workspace.\n");
  sf_mex_import_named("J1", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c14_d1, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c14_J1 = c14_d1;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'J2' in the parent workspace.\n");
  sf_mex_import_named("J2", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      &c14_d2, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c14_J2 = c14_d2;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c14_my_system(SFc14_my_systemInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c14_my_system(SFc14_my_systemInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c14_update_debugger_state_c14_my_system
  (SFc14_my_systemInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c14_my_system(SFc14_my_systemInstanceStruct *
  chartInstance)
{
  const mxArray *c14_st;
  const mxArray *c14_y = NULL;
  int32_T c14_i0;
  real_T c14_u[144];
  const mxArray *c14_b_y = NULL;
  int32_T c14_i1;
  real_T c14_b_u[144];
  const mxArray *c14_c_y = NULL;
  uint8_T c14_hoistedGlobal;
  uint8_T c14_c_u;
  const mxArray *c14_d_y = NULL;
  real_T (*c14_R)[144];
  real_T (*c14_Q)[144];
  c14_R = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 2);
  c14_Q = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 1);
  c14_st = NULL;
  c14_st = NULL;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_createcellarray(3), FALSE);
  for (c14_i0 = 0; c14_i0 < 144; c14_i0++) {
    c14_u[c14_i0] = (*c14_Q)[c14_i0];
  }

  c14_b_y = NULL;
  sf_mex_assign(&c14_b_y, sf_mex_create("y", c14_u, 0, 0U, 1U, 0U, 2, 12, 12),
                FALSE);
  sf_mex_setcell(c14_y, 0, c14_b_y);
  for (c14_i1 = 0; c14_i1 < 144; c14_i1++) {
    c14_b_u[c14_i1] = (*c14_R)[c14_i1];
  }

  c14_c_y = NULL;
  sf_mex_assign(&c14_c_y, sf_mex_create("y", c14_b_u, 0, 0U, 1U, 0U, 2, 12, 12),
                FALSE);
  sf_mex_setcell(c14_y, 1, c14_c_y);
  c14_hoistedGlobal = chartInstance->c14_is_active_c14_my_system;
  c14_c_u = c14_hoistedGlobal;
  c14_d_y = NULL;
  sf_mex_assign(&c14_d_y, sf_mex_create("y", &c14_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c14_y, 2, c14_d_y);
  sf_mex_assign(&c14_st, c14_y, FALSE);
  return c14_st;
}

static void set_sim_state_c14_my_system(SFc14_my_systemInstanceStruct
  *chartInstance, const mxArray *c14_st)
{
  const mxArray *c14_u;
  real_T c14_dv0[144];
  int32_T c14_i2;
  real_T c14_dv1[144];
  int32_T c14_i3;
  real_T (*c14_Q)[144];
  real_T (*c14_R)[144];
  c14_R = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 2);
  c14_Q = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c14_doneDoubleBufferReInit = TRUE;
  c14_u = sf_mex_dup(c14_st);
  c14_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c14_u, 0)), "Q",
                       c14_dv0);
  for (c14_i2 = 0; c14_i2 < 144; c14_i2++) {
    (*c14_Q)[c14_i2] = c14_dv0[c14_i2];
  }

  c14_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c14_u, 1)), "R",
                       c14_dv1);
  for (c14_i3 = 0; c14_i3 < 144; c14_i3++) {
    (*c14_R)[c14_i3] = c14_dv1[c14_i3];
  }

  chartInstance->c14_is_active_c14_my_system = c14_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c14_u, 2)),
     "is_active_c14_my_system");
  sf_mex_destroy(&c14_u);
  c14_update_debugger_state_c14_my_system(chartInstance);
  sf_mex_destroy(&c14_st);
}

static void finalize_c14_my_system(SFc14_my_systemInstanceStruct *chartInstance)
{
}

static void sf_c14_my_system(SFc14_my_systemInstanceStruct *chartInstance)
{
  int32_T c14_i4;
  int32_T c14_i5;
  real_T c14_hoistedGlobal;
  real_T c14_b_hoistedGlobal;
  real_T c14_c_hoistedGlobal;
  real_T c14_b_sampleTime;
  real_T c14_b_J1;
  real_T c14_b_J2;
  uint32_T c14_debug_family_var_map[20];
  real_T c14_accNoise;
  real_T c14_gyroNoise;
  real_T c14_magNoise1;
  real_T c14_magNoise2;
  real_T c14_magNoise3;
  real_T c14_magNoise4;
  real_T c14_maxA;
  real_T c14_val1;
  real_T c14_val2;
  real_T c14_val3;
  real_T c14_val4;
  real_T c14_val5;
  real_T c14_val6;
  real_T c14_nargin = 3.0;
  real_T c14_nargout = 2.0;
  real_T c14_Q[144];
  real_T c14_R[144];
  real_T c14_b;
  real_T c14_y;
  real_T c14_a;
  real_T c14_b_b;
  real_T c14_c_b;
  real_T c14_b_y;
  real_T c14_b_a;
  real_T c14_d_b;
  int32_T c14_i6;
  int32_T c14_i7;
  int32_T c14_i8;
  int32_T c14_i9;
  int32_T c14_i10;
  int32_T c14_i11;
  int32_T c14_i12;
  int32_T c14_i13;
  int32_T c14_i14;
  int32_T c14_i15;
  real_T (*c14_b_Q)[144];
  real_T (*c14_b_R)[144];
  c14_b_R = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 2);
  c14_b_Q = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 1);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, chartInstance->c14_sfEvent);
  for (c14_i4 = 0; c14_i4 < 144; c14_i4++) {
    _SFD_DATA_RANGE_CHECK((*c14_b_Q)[c14_i4], 0U);
  }

  for (c14_i5 = 0; c14_i5 < 144; c14_i5++) {
    _SFD_DATA_RANGE_CHECK((*c14_b_R)[c14_i5], 1U);
  }

  _SFD_DATA_RANGE_CHECK(chartInstance->c14_sampleTime, 2U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c14_J1, 3U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c14_J2, 4U);
  chartInstance->c14_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c14_sfEvent);
  c14_hoistedGlobal = chartInstance->c14_sampleTime;
  c14_b_hoistedGlobal = chartInstance->c14_J1;
  c14_c_hoistedGlobal = chartInstance->c14_J2;
  c14_b_sampleTime = c14_hoistedGlobal;
  c14_b_J1 = c14_b_hoistedGlobal;
  c14_b_J2 = c14_c_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 20U, 20U, c14_debug_family_names,
    c14_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_accNoise, 0U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_gyroNoise, 1U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_magNoise1, 2U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_magNoise2, 3U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_magNoise3, 4U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_magNoise4, 5U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c14_maxA, 6U, c14_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_val1, 7U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_val2, 8U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_val3, 9U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_val4, 10U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_val5, 11U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_val6, 12U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_nargin, 13U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_nargout, 14U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_b_sampleTime, 15U,
    c14_b_sf_marshallOut, c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_b_J1, 16U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_b_J2, 17U, c14_b_sf_marshallOut,
    c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c14_Q, 18U, c14_sf_marshallOut,
    c14_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c14_R, 19U, c14_sf_marshallOut,
    c14_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 4);
  c14_accNoise = 0.000213858;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 5);
  c14_gyroNoise = 1.8992E-5;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 6);
  c14_magNoise1 = 0.0029;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 7);
  c14_magNoise2 = 0.0029;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 8);
  c14_magNoise3 = 0.0029;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 9);
  c14_magNoise4 = 0.0029;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 11);
  c14_maxA = 20.0;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 12);
  c14_b = c14_b_J1;
  c14_y = 20.0 * c14_b;
  c14_a = c14_mpower(chartInstance, c14_y);
  c14_b_b = c14_mpower(chartInstance, c14_b_sampleTime);
  c14_val1 = c14_a * c14_b_b;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 13);
  c14_c_b = c14_b_J2;
  c14_b_y = 20.0 * c14_c_b;
  c14_b_a = c14_mpower(chartInstance, c14_b_y);
  c14_d_b = c14_mpower(chartInstance, c14_b_sampleTime);
  c14_val2 = c14_b_a * c14_d_b;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 15);
  c14_val3 = 2.05E-5;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 16);
  c14_val4 = 0.00205;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 17);
  c14_val5 = 2.05E-5;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 18);
  c14_val6 = 0.00205;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 20);
  c14_Q[0] = c14_val1;
  c14_Q[12] = 0.0;
  c14_Q[24] = 0.0;
  c14_Q[36] = 0.0;
  c14_Q[48] = 0.0;
  c14_Q[60] = 0.0;
  c14_Q[72] = 0.0;
  c14_Q[84] = 0.0;
  c14_Q[96] = 0.0;
  c14_Q[108] = 0.0;
  c14_Q[120] = 0.0;
  c14_Q[132] = 0.0;
  c14_Q[1] = 0.0;
  c14_Q[13] = c14_val1;
  c14_Q[25] = 0.0;
  c14_Q[37] = 0.0;
  c14_Q[49] = 0.0;
  c14_Q[61] = 0.0;
  c14_Q[73] = 0.0;
  c14_Q[85] = 0.0;
  c14_Q[97] = 0.0;
  c14_Q[109] = 0.0;
  c14_Q[121] = 0.0;
  c14_Q[133] = 0.0;
  c14_Q[2] = 0.0;
  c14_Q[14] = 0.0;
  c14_Q[26] = c14_val2;
  c14_Q[38] = 0.0;
  c14_Q[50] = 0.0;
  c14_Q[62] = 0.0;
  c14_Q[74] = 0.0;
  c14_Q[86] = 0.0;
  c14_Q[98] = 0.0;
  c14_Q[110] = 0.0;
  c14_Q[122] = 0.0;
  c14_Q[134] = 0.0;
  c14_Q[3] = 0.0;
  c14_Q[15] = 0.0;
  c14_Q[27] = 0.0;
  c14_Q[39] = c14_val2;
  c14_Q[51] = 0.0;
  c14_Q[63] = 0.0;
  c14_Q[75] = 0.0;
  c14_Q[87] = 0.0;
  c14_Q[99] = 0.0;
  c14_Q[111] = 0.0;
  c14_Q[123] = 0.0;
  c14_Q[135] = 0.0;
  c14_Q[4] = 0.0;
  c14_Q[16] = 0.0;
  c14_Q[28] = 0.0;
  c14_Q[40] = 0.0;
  c14_Q[52] = c14_val3;
  c14_Q[64] = 0.0;
  c14_Q[76] = 0.0;
  c14_Q[88] = 0.0;
  c14_Q[100] = 0.0;
  c14_Q[112] = 0.0;
  c14_Q[124] = 0.0;
  c14_Q[136] = 0.0;
  c14_Q[5] = 0.0;
  c14_Q[17] = 0.0;
  c14_Q[29] = 0.0;
  c14_Q[41] = 0.0;
  c14_Q[53] = 0.0;
  c14_Q[65] = c14_val4;
  c14_Q[77] = 0.0;
  c14_Q[89] = 0.0;
  c14_Q[101] = 0.0;
  c14_Q[113] = 0.0;
  c14_Q[125] = 0.0;
  c14_Q[137] = 0.0;
  c14_Q[6] = 0.0;
  c14_Q[18] = 0.0;
  c14_Q[30] = 0.0;
  c14_Q[42] = 0.0;
  c14_Q[54] = 0.0;
  c14_Q[66] = 0.0;
  c14_Q[78] = c14_val5;
  c14_Q[90] = 0.0;
  c14_Q[102] = 0.0;
  c14_Q[114] = 0.0;
  c14_Q[126] = 0.0;
  c14_Q[138] = 0.0;
  c14_Q[7] = 0.0;
  c14_Q[19] = 0.0;
  c14_Q[31] = 0.0;
  c14_Q[43] = 0.0;
  c14_Q[55] = 0.0;
  c14_Q[67] = 0.0;
  c14_Q[79] = 0.0;
  c14_Q[91] = c14_val6;
  c14_Q[103] = 0.0;
  c14_Q[115] = 0.0;
  c14_Q[127] = 0.0;
  c14_Q[139] = 0.0;
  c14_i6 = 0;
  for (c14_i7 = 0; c14_i7 < 12; c14_i7++) {
    c14_Q[c14_i6 + 8] = 0.0;
    c14_i6 += 12;
  }

  c14_i8 = 0;
  for (c14_i9 = 0; c14_i9 < 12; c14_i9++) {
    c14_Q[c14_i8 + 9] = 0.0;
    c14_i8 += 12;
  }

  c14_i10 = 0;
  for (c14_i11 = 0; c14_i11 < 12; c14_i11++) {
    c14_Q[c14_i10 + 10] = 0.0;
    c14_i10 += 12;
  }

  c14_i12 = 0;
  for (c14_i13 = 0; c14_i13 < 12; c14_i13++) {
    c14_Q[c14_i12 + 11] = 0.0;
    c14_i12 += 12;
  }

  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 34);
  c14_R[0] = c14_accNoise;
  c14_R[12] = 0.0;
  c14_R[24] = 0.0;
  c14_R[36] = 0.0;
  c14_R[48] = 0.0;
  c14_R[60] = 0.0;
  c14_R[72] = 0.0;
  c14_R[84] = 0.0;
  c14_R[96] = 0.0;
  c14_R[108] = 0.0;
  c14_R[120] = 0.0;
  c14_R[132] = 0.0;
  c14_R[1] = 0.0;
  c14_R[13] = c14_accNoise;
  c14_R[25] = 0.0;
  c14_R[37] = 0.0;
  c14_R[49] = 0.0;
  c14_R[61] = 0.0;
  c14_R[73] = 0.0;
  c14_R[85] = 0.0;
  c14_R[97] = 0.0;
  c14_R[109] = 0.0;
  c14_R[121] = 0.0;
  c14_R[133] = 0.0;
  c14_R[2] = 0.0;
  c14_R[14] = 0.0;
  c14_R[26] = c14_gyroNoise;
  c14_R[38] = 0.0;
  c14_R[50] = 0.0;
  c14_R[62] = 0.0;
  c14_R[74] = 0.0;
  c14_R[86] = 0.0;
  c14_R[98] = 0.0;
  c14_R[110] = 0.0;
  c14_R[122] = 0.0;
  c14_R[134] = 0.0;
  c14_R[3] = 0.0;
  c14_R[15] = 0.0;
  c14_R[27] = 0.0;
  c14_R[39] = c14_gyroNoise;
  c14_R[51] = 0.0;
  c14_R[63] = 0.0;
  c14_R[75] = 0.0;
  c14_R[87] = 0.0;
  c14_R[99] = 0.0;
  c14_R[111] = 0.0;
  c14_R[123] = 0.0;
  c14_R[135] = 0.0;
  c14_R[4] = 0.0;
  c14_R[16] = 0.0;
  c14_R[28] = 0.0;
  c14_R[40] = 0.0;
  c14_R[52] = c14_magNoise1;
  c14_R[64] = 0.0;
  c14_R[76] = 0.0;
  c14_R[88] = 0.0;
  c14_R[100] = 0.0;
  c14_R[112] = 0.0;
  c14_R[124] = 0.0;
  c14_R[136] = 0.0;
  c14_R[5] = 0.0;
  c14_R[17] = 0.0;
  c14_R[29] = 0.0;
  c14_R[41] = 0.0;
  c14_R[53] = 0.0;
  c14_R[65] = c14_magNoise2;
  c14_R[77] = 0.0;
  c14_R[89] = 0.0;
  c14_R[101] = 0.0;
  c14_R[113] = 0.0;
  c14_R[125] = 0.0;
  c14_R[137] = 0.0;
  c14_R[6] = 0.0;
  c14_R[18] = 0.0;
  c14_R[30] = 0.0;
  c14_R[42] = 0.0;
  c14_R[54] = 0.0;
  c14_R[66] = 0.0;
  c14_R[78] = c14_accNoise;
  c14_R[90] = 0.0;
  c14_R[102] = 0.0;
  c14_R[114] = 0.0;
  c14_R[126] = 0.0;
  c14_R[138] = 0.0;
  c14_R[7] = 0.0;
  c14_R[19] = 0.0;
  c14_R[31] = 0.0;
  c14_R[43] = 0.0;
  c14_R[55] = 0.0;
  c14_R[67] = 0.0;
  c14_R[79] = 0.0;
  c14_R[91] = c14_accNoise;
  c14_R[103] = 0.0;
  c14_R[115] = 0.0;
  c14_R[127] = 0.0;
  c14_R[139] = 0.0;
  c14_R[8] = 0.0;
  c14_R[20] = 0.0;
  c14_R[32] = 0.0;
  c14_R[44] = 0.0;
  c14_R[56] = 0.0;
  c14_R[68] = 0.0;
  c14_R[80] = 0.0;
  c14_R[92] = 0.0;
  c14_R[104] = c14_gyroNoise;
  c14_R[116] = 0.0;
  c14_R[128] = 0.0;
  c14_R[140] = 0.0;
  c14_R[9] = 0.0;
  c14_R[21] = 0.0;
  c14_R[33] = 0.0;
  c14_R[45] = 0.0;
  c14_R[57] = 0.0;
  c14_R[69] = 0.0;
  c14_R[81] = 0.0;
  c14_R[93] = 0.0;
  c14_R[105] = 0.0;
  c14_R[117] = c14_gyroNoise;
  c14_R[129] = 0.0;
  c14_R[141] = 0.0;
  c14_R[10] = 0.0;
  c14_R[22] = 0.0;
  c14_R[34] = 0.0;
  c14_R[46] = 0.0;
  c14_R[58] = 0.0;
  c14_R[70] = 0.0;
  c14_R[82] = 0.0;
  c14_R[94] = 0.0;
  c14_R[106] = 0.0;
  c14_R[118] = 0.0;
  c14_R[130] = c14_magNoise3;
  c14_R[142] = 0.0;
  c14_R[11] = 0.0;
  c14_R[23] = 0.0;
  c14_R[35] = 0.0;
  c14_R[47] = 0.0;
  c14_R[59] = 0.0;
  c14_R[71] = 0.0;
  c14_R[83] = 0.0;
  c14_R[95] = 0.0;
  c14_R[107] = 0.0;
  c14_R[119] = 0.0;
  c14_R[131] = 0.0;
  c14_R[143] = c14_magNoise4;
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, -34);
  _SFD_SYMBOL_SCOPE_POP();
  for (c14_i14 = 0; c14_i14 < 144; c14_i14++) {
    (*c14_b_Q)[c14_i14] = c14_Q[c14_i14];
  }

  for (c14_i15 = 0; c14_i15 < 144; c14_i15++) {
    (*c14_b_R)[c14_i15] = c14_R[c14_i15];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c14_sfEvent);
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_my_systemMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void initSimStructsc14_my_system(SFc14_my_systemInstanceStruct
  *chartInstance)
{
}

static void registerMessagesc14_my_system(SFc14_my_systemInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c14_machineNumber, uint32_T
  c14_chartNumber)
{
}

static const mxArray *c14_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData)
{
  const mxArray *c14_mxArrayOutData = NULL;
  int32_T c14_i16;
  int32_T c14_i17;
  int32_T c14_i18;
  real_T c14_b_inData[144];
  int32_T c14_i19;
  int32_T c14_i20;
  int32_T c14_i21;
  real_T c14_u[144];
  const mxArray *c14_y = NULL;
  SFc14_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc14_my_systemInstanceStruct *)chartInstanceVoid;
  c14_mxArrayOutData = NULL;
  c14_i16 = 0;
  for (c14_i17 = 0; c14_i17 < 12; c14_i17++) {
    for (c14_i18 = 0; c14_i18 < 12; c14_i18++) {
      c14_b_inData[c14_i18 + c14_i16] = (*(real_T (*)[144])c14_inData)[c14_i18 +
        c14_i16];
    }

    c14_i16 += 12;
  }

  c14_i19 = 0;
  for (c14_i20 = 0; c14_i20 < 12; c14_i20++) {
    for (c14_i21 = 0; c14_i21 < 12; c14_i21++) {
      c14_u[c14_i21 + c14_i19] = c14_b_inData[c14_i21 + c14_i19];
    }

    c14_i19 += 12;
  }

  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", c14_u, 0, 0U, 1U, 0U, 2, 12, 12),
                FALSE);
  sf_mex_assign(&c14_mxArrayOutData, c14_y, FALSE);
  return c14_mxArrayOutData;
}

static void c14_emlrt_marshallIn(SFc14_my_systemInstanceStruct *chartInstance,
  const mxArray *c14_R, const char_T *c14_identifier, real_T c14_y[144])
{
  emlrtMsgIdentifier c14_thisId;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_R), &c14_thisId, c14_y);
  sf_mex_destroy(&c14_R);
}

static void c14_b_emlrt_marshallIn(SFc14_my_systemInstanceStruct *chartInstance,
  const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId, real_T c14_y[144])
{
  real_T c14_dv2[144];
  int32_T c14_i22;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), c14_dv2, 1, 0, 0U, 1, 0U, 2, 12,
                12);
  for (c14_i22 = 0; c14_i22 < 144; c14_i22++) {
    c14_y[c14_i22] = c14_dv2[c14_i22];
  }

  sf_mex_destroy(&c14_u);
}

static void c14_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData)
{
  const mxArray *c14_R;
  const char_T *c14_identifier;
  emlrtMsgIdentifier c14_thisId;
  real_T c14_y[144];
  int32_T c14_i23;
  int32_T c14_i24;
  int32_T c14_i25;
  SFc14_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc14_my_systemInstanceStruct *)chartInstanceVoid;
  c14_R = sf_mex_dup(c14_mxArrayInData);
  c14_identifier = c14_varName;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_R), &c14_thisId, c14_y);
  sf_mex_destroy(&c14_R);
  c14_i23 = 0;
  for (c14_i24 = 0; c14_i24 < 12; c14_i24++) {
    for (c14_i25 = 0; c14_i25 < 12; c14_i25++) {
      (*(real_T (*)[144])c14_outData)[c14_i25 + c14_i23] = c14_y[c14_i25 +
        c14_i23];
    }

    c14_i23 += 12;
  }

  sf_mex_destroy(&c14_mxArrayInData);
}

static const mxArray *c14_b_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData)
{
  const mxArray *c14_mxArrayOutData = NULL;
  real_T c14_u;
  const mxArray *c14_y = NULL;
  SFc14_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc14_my_systemInstanceStruct *)chartInstanceVoid;
  c14_mxArrayOutData = NULL;
  c14_u = *(real_T *)c14_inData;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", &c14_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c14_mxArrayOutData, c14_y, FALSE);
  return c14_mxArrayOutData;
}

static real_T c14_c_emlrt_marshallIn(SFc14_my_systemInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId)
{
  real_T c14_y;
  real_T c14_d3;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), &c14_d3, 1, 0, 0U, 0, 0U, 0);
  c14_y = c14_d3;
  sf_mex_destroy(&c14_u);
  return c14_y;
}

static void c14_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData)
{
  const mxArray *c14_b_J2;
  const char_T *c14_identifier;
  emlrtMsgIdentifier c14_thisId;
  real_T c14_y;
  SFc14_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc14_my_systemInstanceStruct *)chartInstanceVoid;
  c14_b_J2 = sf_mex_dup(c14_mxArrayInData);
  c14_identifier = c14_varName;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_y = c14_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_b_J2),
    &c14_thisId);
  sf_mex_destroy(&c14_b_J2);
  *(real_T *)c14_outData = c14_y;
  sf_mex_destroy(&c14_mxArrayInData);
}

const mxArray *sf_c14_my_system_get_eml_resolved_functions_info(void)
{
  const mxArray *c14_nameCaptureInfo;
  c14_ResolvedFunctionInfo c14_info[9];
  c14_ResolvedFunctionInfo (*c14_b_info)[9];
  const mxArray *c14_m0 = NULL;
  int32_T c14_i26;
  c14_ResolvedFunctionInfo *c14_r0;
  c14_nameCaptureInfo = NULL;
  c14_nameCaptureInfo = NULL;
  c14_b_info = (c14_ResolvedFunctionInfo (*)[9])c14_info;
  (*c14_b_info)[0].context = "";
  (*c14_b_info)[0].name = "mtimes";
  (*c14_b_info)[0].dominantType = "double";
  (*c14_b_info)[0].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  (*c14_b_info)[0].fileTimeLo = 1289516092U;
  (*c14_b_info)[0].fileTimeHi = 0U;
  (*c14_b_info)[0].mFileTimeLo = 0U;
  (*c14_b_info)[0].mFileTimeHi = 0U;
  (*c14_b_info)[1].context = "";
  (*c14_b_info)[1].name = "mpower";
  (*c14_b_info)[1].dominantType = "double";
  (*c14_b_info)[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  (*c14_b_info)[1].fileTimeLo = 1286818842U;
  (*c14_b_info)[1].fileTimeHi = 0U;
  (*c14_b_info)[1].mFileTimeLo = 0U;
  (*c14_b_info)[1].mFileTimeHi = 0U;
  (*c14_b_info)[2].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  (*c14_b_info)[2].name = "power";
  (*c14_b_info)[2].dominantType = "double";
  (*c14_b_info)[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m";
  (*c14_b_info)[2].fileTimeLo = 1348191930U;
  (*c14_b_info)[2].fileTimeHi = 0U;
  (*c14_b_info)[2].mFileTimeLo = 0U;
  (*c14_b_info)[2].mFileTimeHi = 0U;
  (*c14_b_info)[3].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  (*c14_b_info)[3].name = "eml_scalar_eg";
  (*c14_b_info)[3].dominantType = "double";
  (*c14_b_info)[3].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  (*c14_b_info)[3].fileTimeLo = 1286818796U;
  (*c14_b_info)[3].fileTimeHi = 0U;
  (*c14_b_info)[3].mFileTimeLo = 0U;
  (*c14_b_info)[3].mFileTimeHi = 0U;
  (*c14_b_info)[4].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  (*c14_b_info)[4].name = "eml_scalexp_alloc";
  (*c14_b_info)[4].dominantType = "double";
  (*c14_b_info)[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  (*c14_b_info)[4].fileTimeLo = 1352421260U;
  (*c14_b_info)[4].fileTimeHi = 0U;
  (*c14_b_info)[4].mFileTimeLo = 0U;
  (*c14_b_info)[4].mFileTimeHi = 0U;
  (*c14_b_info)[5].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  (*c14_b_info)[5].name = "floor";
  (*c14_b_info)[5].dominantType = "double";
  (*c14_b_info)[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  (*c14_b_info)[5].fileTimeLo = 1343830380U;
  (*c14_b_info)[5].fileTimeHi = 0U;
  (*c14_b_info)[5].mFileTimeLo = 0U;
  (*c14_b_info)[5].mFileTimeHi = 0U;
  (*c14_b_info)[6].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  (*c14_b_info)[6].name = "eml_scalar_floor";
  (*c14_b_info)[6].dominantType = "double";
  (*c14_b_info)[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m";
  (*c14_b_info)[6].fileTimeLo = 1286818726U;
  (*c14_b_info)[6].fileTimeHi = 0U;
  (*c14_b_info)[6].mFileTimeLo = 0U;
  (*c14_b_info)[6].mFileTimeHi = 0U;
  (*c14_b_info)[7].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power";
  (*c14_b_info)[7].name = "eml_scalar_eg";
  (*c14_b_info)[7].dominantType = "double";
  (*c14_b_info)[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  (*c14_b_info)[7].fileTimeLo = 1286818796U;
  (*c14_b_info)[7].fileTimeHi = 0U;
  (*c14_b_info)[7].mFileTimeLo = 0U;
  (*c14_b_info)[7].mFileTimeHi = 0U;
  (*c14_b_info)[8].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power";
  (*c14_b_info)[8].name = "mtimes";
  (*c14_b_info)[8].dominantType = "double";
  (*c14_b_info)[8].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  (*c14_b_info)[8].fileTimeLo = 1289516092U;
  (*c14_b_info)[8].fileTimeHi = 0U;
  (*c14_b_info)[8].mFileTimeLo = 0U;
  (*c14_b_info)[8].mFileTimeHi = 0U;
  sf_mex_assign(&c14_m0, sf_mex_createstruct("nameCaptureInfo", 1, 9), FALSE);
  for (c14_i26 = 0; c14_i26 < 9; c14_i26++) {
    c14_r0 = &c14_info[c14_i26];
    sf_mex_addfield(c14_m0, sf_mex_create("nameCaptureInfo", c14_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c14_r0->context)), "context", "nameCaptureInfo",
                    c14_i26);
    sf_mex_addfield(c14_m0, sf_mex_create("nameCaptureInfo", c14_r0->name, 15,
      0U, 0U, 0U, 2, 1, strlen(c14_r0->name)), "name", "nameCaptureInfo",
                    c14_i26);
    sf_mex_addfield(c14_m0, sf_mex_create("nameCaptureInfo",
      c14_r0->dominantType, 15, 0U, 0U, 0U, 2, 1, strlen(c14_r0->dominantType)),
                    "dominantType", "nameCaptureInfo", c14_i26);
    sf_mex_addfield(c14_m0, sf_mex_create("nameCaptureInfo", c14_r0->resolved,
      15, 0U, 0U, 0U, 2, 1, strlen(c14_r0->resolved)), "resolved",
                    "nameCaptureInfo", c14_i26);
    sf_mex_addfield(c14_m0, sf_mex_create("nameCaptureInfo", &c14_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c14_i26);
    sf_mex_addfield(c14_m0, sf_mex_create("nameCaptureInfo", &c14_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c14_i26);
    sf_mex_addfield(c14_m0, sf_mex_create("nameCaptureInfo",
      &c14_r0->mFileTimeLo, 7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo",
                    c14_i26);
    sf_mex_addfield(c14_m0, sf_mex_create("nameCaptureInfo",
      &c14_r0->mFileTimeHi, 7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo",
                    c14_i26);
  }

  sf_mex_assign(&c14_nameCaptureInfo, c14_m0, FALSE);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c14_nameCaptureInfo);
  return c14_nameCaptureInfo;
}

static real_T c14_mpower(SFc14_my_systemInstanceStruct *chartInstance, real_T
  c14_a)
{
  real_T c14_b_a;
  real_T c14_c_a;
  real_T c14_ak;
  real_T c14_d_a;
  real_T c14_e_a;
  real_T c14_b;
  c14_b_a = c14_a;
  c14_c_a = c14_b_a;
  c14_eml_scalar_eg(chartInstance);
  c14_ak = c14_c_a;
  c14_d_a = c14_ak;
  c14_eml_scalar_eg(chartInstance);
  c14_e_a = c14_d_a;
  c14_b = c14_d_a;
  return c14_e_a * c14_b;
}

static void c14_eml_scalar_eg(SFc14_my_systemInstanceStruct *chartInstance)
{
}

static const mxArray *c14_c_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData)
{
  const mxArray *c14_mxArrayOutData = NULL;
  int32_T c14_u;
  const mxArray *c14_y = NULL;
  SFc14_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc14_my_systemInstanceStruct *)chartInstanceVoid;
  c14_mxArrayOutData = NULL;
  c14_u = *(int32_T *)c14_inData;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", &c14_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c14_mxArrayOutData, c14_y, FALSE);
  return c14_mxArrayOutData;
}

static int32_T c14_d_emlrt_marshallIn(SFc14_my_systemInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId)
{
  int32_T c14_y;
  int32_T c14_i27;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), &c14_i27, 1, 6, 0U, 0, 0U, 0);
  c14_y = c14_i27;
  sf_mex_destroy(&c14_u);
  return c14_y;
}

static void c14_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData)
{
  const mxArray *c14_b_sfEvent;
  const char_T *c14_identifier;
  emlrtMsgIdentifier c14_thisId;
  int32_T c14_y;
  SFc14_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc14_my_systemInstanceStruct *)chartInstanceVoid;
  c14_b_sfEvent = sf_mex_dup(c14_mxArrayInData);
  c14_identifier = c14_varName;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_y = c14_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_b_sfEvent),
    &c14_thisId);
  sf_mex_destroy(&c14_b_sfEvent);
  *(int32_T *)c14_outData = c14_y;
  sf_mex_destroy(&c14_mxArrayInData);
}

static uint8_T c14_e_emlrt_marshallIn(SFc14_my_systemInstanceStruct
  *chartInstance, const mxArray *c14_b_is_active_c14_my_system, const char_T
  *c14_identifier)
{
  uint8_T c14_y;
  emlrtMsgIdentifier c14_thisId;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_y = c14_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c14_b_is_active_c14_my_system), &c14_thisId);
  sf_mex_destroy(&c14_b_is_active_c14_my_system);
  return c14_y;
}

static uint8_T c14_f_emlrt_marshallIn(SFc14_my_systemInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId)
{
  uint8_T c14_y;
  uint8_T c14_u0;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), &c14_u0, 1, 3, 0U, 0, 0U, 0);
  c14_y = c14_u0;
  sf_mex_destroy(&c14_u);
  return c14_y;
}

static void init_dsm_address_info(SFc14_my_systemInstanceStruct *chartInstance)
{
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c14_my_system_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2843527626U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2778276455U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(4066182034U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1603789138U);
}

mxArray *sf_c14_my_system_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("J5raknxOvQuOHBNJbqIN2F");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(12);
      pr[1] = (double)(12);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(12);
      pr[1] = (double)(12);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c14_my_system_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

static const mxArray *sf_get_sim_state_info_c14_my_system(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[5],T\"Q\",},{M[1],M[6],T\"R\",},{M[8],M[0],T\"is_active_c14_my_system\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c14_my_system_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc14_my_systemInstanceStruct *chartInstance;
    chartInstance = (SFc14_my_systemInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _my_systemMachineNumber_,
           14,
           1,
           1,
           5,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           ssGetPath(S),
           (void *)S);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          init_script_number_translation(_my_systemMachineNumber_,
            chartInstance->chartNumber);
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_my_systemMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _my_systemMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,2,0,1,"Q");
          _SFD_SET_DATA_PROPS(1,2,0,1,"R");
          _SFD_SET_DATA_PROPS(2,10,0,0,"sampleTime");
          _SFD_SET_DATA_PROPS(3,10,0,0,"J1");
          _SFD_SET_DATA_PROPS(4,10,0,0,"J2");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,2654);
        _SFD_TRANS_COV_WTS(0,0,0,1,0);
        if (chartAlreadyPresent==0) {
          _SFD_TRANS_COV_MAPS(0,
                              0,NULL,NULL,
                              0,NULL,NULL,
                              1,NULL,NULL,
                              0,NULL,NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 12;
          dimVector[1]= 12;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c14_sf_marshallOut,(MexInFcnForType)
            c14_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 12;
          dimVector[1]= 12;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c14_sf_marshallOut,(MexInFcnForType)
            c14_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c14_b_sf_marshallOut,(MexInFcnForType)
          c14_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c14_b_sf_marshallOut,(MexInFcnForType)
          c14_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c14_b_sf_marshallOut,(MexInFcnForType)
          c14_b_sf_marshallIn);

        {
          real_T (*c14_Q)[144];
          real_T (*c14_R)[144];
          c14_R = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 2);
          c14_Q = (real_T (*)[144])ssGetOutputPortSignal(chartInstance->S, 1);
          _SFD_SET_DATA_VALUE_PTR(0U, *c14_Q);
          _SFD_SET_DATA_VALUE_PTR(1U, *c14_R);
          _SFD_SET_DATA_VALUE_PTR(2U, &chartInstance->c14_sampleTime);
          _SFD_SET_DATA_VALUE_PTR(3U, &chartInstance->c14_J1);
          _SFD_SET_DATA_VALUE_PTR(4U, &chartInstance->c14_J2);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _my_systemMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "ttnNd5JRgtHsnwXsxYXZ7";
}

static void sf_opaque_initialize_c14_my_system(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc14_my_systemInstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c14_my_system((SFc14_my_systemInstanceStruct*)
    chartInstanceVar);
  initialize_c14_my_system((SFc14_my_systemInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c14_my_system(void *chartInstanceVar)
{
  enable_c14_my_system((SFc14_my_systemInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c14_my_system(void *chartInstanceVar)
{
  disable_c14_my_system((SFc14_my_systemInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c14_my_system(void *chartInstanceVar)
{
  sf_c14_my_system((SFc14_my_systemInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c14_my_system(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c14_my_system
    ((SFc14_my_systemInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c14_my_system();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c14_my_system(SimStruct* S, const mxArray *
  st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c14_my_system();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c14_my_system((SFc14_my_systemInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c14_my_system(SimStruct* S)
{
  return sf_internal_get_sim_state_c14_my_system(S);
}

static void sf_opaque_set_sim_state_c14_my_system(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c14_my_system(S, st);
}

static void sf_opaque_terminate_c14_my_system(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc14_my_systemInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_my_system_optimization_info();
    }

    finalize_c14_my_system((SFc14_my_systemInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc14_my_system((SFc14_my_systemInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c14_my_system(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c14_my_system((SFc14_my_systemInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c14_my_system(SimStruct *S)
{
  /* Actual parameters from chart:
     J1 J2 sampleTime
   */
  const char_T *rtParamNames[] = { "J1", "J2", "sampleTime" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for J1*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);

  /* registration for J2*/
  ssRegDlgParamAsRunTimeParam(S, 1, 1, rtParamNames[1], SS_DOUBLE);

  /* registration for sampleTime*/
  ssRegDlgParamAsRunTimeParam(S, 2, 2, rtParamNames[2], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_my_system_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      14);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,14,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,14,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,14);
    if (chartIsInlinable) {
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,14,2);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=2; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,14);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2253267467U));
  ssSetChecksum1(S,(2139347840U));
  ssSetChecksum2(S,(3723217707U));
  ssSetChecksum3(S,(3890106939U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c14_my_system(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c14_my_system(SimStruct *S)
{
  SFc14_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc14_my_systemInstanceStruct *)utMalloc(sizeof
    (SFc14_my_systemInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc14_my_systemInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c14_my_system;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c14_my_system;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c14_my_system;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c14_my_system;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c14_my_system;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c14_my_system;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c14_my_system;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c14_my_system;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c14_my_system;
  chartInstance->chartInfo.mdlStart = mdlStart_c14_my_system;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c14_my_system;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->S = S;
  ssSetUserData(S,(void *)(&(chartInstance->chartInfo)));/* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c14_my_system_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c14_my_system(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c14_my_system(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c14_my_system(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c14_my_system_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
