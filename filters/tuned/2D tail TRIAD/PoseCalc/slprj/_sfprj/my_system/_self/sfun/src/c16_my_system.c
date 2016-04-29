/* Include files */

#include <stddef.h>
#include "blas.h"
#include "my_system_sfun.h"
#include "c16_my_system.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "my_system_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c16_debug_family_names[6] = { "nargin", "nargout", "u",
  "sampleTime", "y", "states" };

/* Function Declarations */
static void initialize_c16_my_system(SFc16_my_systemInstanceStruct
  *chartInstance);
static void initialize_params_c16_my_system(SFc16_my_systemInstanceStruct
  *chartInstance);
static void enable_c16_my_system(SFc16_my_systemInstanceStruct *chartInstance);
static void disable_c16_my_system(SFc16_my_systemInstanceStruct *chartInstance);
static void c16_update_debugger_state_c16_my_system
  (SFc16_my_systemInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c16_my_system(SFc16_my_systemInstanceStruct *
  chartInstance);
static void set_sim_state_c16_my_system(SFc16_my_systemInstanceStruct
  *chartInstance, const mxArray *c16_st);
static void finalize_c16_my_system(SFc16_my_systemInstanceStruct *chartInstance);
static void sf_c16_my_system(SFc16_my_systemInstanceStruct *chartInstance);
static void initSimStructsc16_my_system(SFc16_my_systemInstanceStruct
  *chartInstance);
static void registerMessagesc16_my_system(SFc16_my_systemInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c16_machineNumber, uint32_T
  c16_chartNumber);
static const mxArray *c16_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData);
static void c16_emlrt_marshallIn(SFc16_my_systemInstanceStruct *chartInstance,
  const mxArray *c16_b_states, const char_T *c16_identifier, real_T c16_y[3]);
static void c16_b_emlrt_marshallIn(SFc16_my_systemInstanceStruct *chartInstance,
  const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId, real_T c16_y[3]);
static void c16_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData);
static const mxArray *c16_b_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData);
static void c16_c_emlrt_marshallIn(SFc16_my_systemInstanceStruct *chartInstance,
  const mxArray *c16_y, const char_T *c16_identifier, real_T c16_b_y[3]);
static void c16_d_emlrt_marshallIn(SFc16_my_systemInstanceStruct *chartInstance,
  const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId, real_T c16_y[3]);
static void c16_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData);
static const mxArray *c16_c_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData);
static real_T c16_e_emlrt_marshallIn(SFc16_my_systemInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId);
static void c16_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData);
static const mxArray *c16_d_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData);
static int32_T c16_f_emlrt_marshallIn(SFc16_my_systemInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId);
static void c16_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData);
static uint8_T c16_g_emlrt_marshallIn(SFc16_my_systemInstanceStruct
  *chartInstance, const mxArray *c16_b_is_active_c16_my_system, const char_T
  *c16_identifier);
static uint8_T c16_h_emlrt_marshallIn(SFc16_my_systemInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId);
static void init_dsm_address_info(SFc16_my_systemInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c16_my_system(SFc16_my_systemInstanceStruct
  *chartInstance)
{
  chartInstance->c16_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c16_states_not_empty = FALSE;
  chartInstance->c16_is_active_c16_my_system = 0U;
}

static void initialize_params_c16_my_system(SFc16_my_systemInstanceStruct
  *chartInstance)
{
  real_T c16_d0;
  sf_set_error_prefix_string(
    "Error evaluating data 'sampleTime' in the parent workspace.\n");
  sf_mex_import_named("sampleTime", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c16_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c16_sampleTime = c16_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c16_my_system(SFc16_my_systemInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c16_my_system(SFc16_my_systemInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c16_update_debugger_state_c16_my_system
  (SFc16_my_systemInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c16_my_system(SFc16_my_systemInstanceStruct *
  chartInstance)
{
  const mxArray *c16_st;
  const mxArray *c16_y = NULL;
  int32_T c16_i0;
  real_T c16_u[3];
  const mxArray *c16_b_y = NULL;
  int32_T c16_i1;
  real_T c16_b_u[3];
  const mxArray *c16_c_y = NULL;
  uint8_T c16_hoistedGlobal;
  uint8_T c16_c_u;
  const mxArray *c16_d_y = NULL;
  real_T (*c16_e_y)[3];
  c16_e_y = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c16_st = NULL;
  c16_st = NULL;
  c16_y = NULL;
  sf_mex_assign(&c16_y, sf_mex_createcellarray(3), FALSE);
  for (c16_i0 = 0; c16_i0 < 3; c16_i0++) {
    c16_u[c16_i0] = (*c16_e_y)[c16_i0];
  }

  c16_b_y = NULL;
  sf_mex_assign(&c16_b_y, sf_mex_create("y", c16_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_setcell(c16_y, 0, c16_b_y);
  for (c16_i1 = 0; c16_i1 < 3; c16_i1++) {
    c16_b_u[c16_i1] = chartInstance->c16_states[c16_i1];
  }

  c16_c_y = NULL;
  if (!chartInstance->c16_states_not_empty) {
    sf_mex_assign(&c16_c_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c16_c_y, sf_mex_create("y", c16_b_u, 0, 0U, 1U, 0U, 1, 3),
                  FALSE);
  }

  sf_mex_setcell(c16_y, 1, c16_c_y);
  c16_hoistedGlobal = chartInstance->c16_is_active_c16_my_system;
  c16_c_u = c16_hoistedGlobal;
  c16_d_y = NULL;
  sf_mex_assign(&c16_d_y, sf_mex_create("y", &c16_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c16_y, 2, c16_d_y);
  sf_mex_assign(&c16_st, c16_y, FALSE);
  return c16_st;
}

static void set_sim_state_c16_my_system(SFc16_my_systemInstanceStruct
  *chartInstance, const mxArray *c16_st)
{
  const mxArray *c16_u;
  real_T c16_dv0[3];
  int32_T c16_i2;
  real_T c16_dv1[3];
  int32_T c16_i3;
  real_T (*c16_y)[3];
  c16_y = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c16_doneDoubleBufferReInit = TRUE;
  c16_u = sf_mex_dup(c16_st);
  c16_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c16_u, 0)),
    "y", c16_dv0);
  for (c16_i2 = 0; c16_i2 < 3; c16_i2++) {
    (*c16_y)[c16_i2] = c16_dv0[c16_i2];
  }

  c16_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c16_u, 1)),
                       "states", c16_dv1);
  for (c16_i3 = 0; c16_i3 < 3; c16_i3++) {
    chartInstance->c16_states[c16_i3] = c16_dv1[c16_i3];
  }

  chartInstance->c16_is_active_c16_my_system = c16_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c16_u, 2)),
     "is_active_c16_my_system");
  sf_mex_destroy(&c16_u);
  c16_update_debugger_state_c16_my_system(chartInstance);
  sf_mex_destroy(&c16_st);
}

static void finalize_c16_my_system(SFc16_my_systemInstanceStruct *chartInstance)
{
}

static void sf_c16_my_system(SFc16_my_systemInstanceStruct *chartInstance)
{
  int32_T c16_i4;
  int32_T c16_i5;
  real_T c16_hoistedGlobal;
  int32_T c16_i6;
  real_T c16_u[3];
  real_T c16_b_sampleTime;
  uint32_T c16_debug_family_var_map[6];
  real_T c16_nargin = 2.0;
  real_T c16_nargout = 1.0;
  real_T c16_y[3];
  int32_T c16_i7;
  int32_T c16_i8;
  real_T c16_b_hoistedGlobal[3];
  int32_T c16_i9;
  real_T c16_a[3];
  real_T c16_b;
  int32_T c16_i10;
  int32_T c16_i11;
  int32_T c16_i12;
  int32_T c16_i13;
  real_T (*c16_b_y)[3];
  real_T (*c16_b_u)[3];
  c16_b_y = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c16_b_u = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 3U, chartInstance->c16_sfEvent);
  for (c16_i4 = 0; c16_i4 < 3; c16_i4++) {
    _SFD_DATA_RANGE_CHECK((*c16_b_u)[c16_i4], 0U);
  }

  for (c16_i5 = 0; c16_i5 < 3; c16_i5++) {
    _SFD_DATA_RANGE_CHECK((*c16_b_y)[c16_i5], 1U);
  }

  _SFD_DATA_RANGE_CHECK(chartInstance->c16_sampleTime, 2U);
  chartInstance->c16_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 3U, chartInstance->c16_sfEvent);
  c16_hoistedGlobal = chartInstance->c16_sampleTime;
  for (c16_i6 = 0; c16_i6 < 3; c16_i6++) {
    c16_u[c16_i6] = (*c16_b_u)[c16_i6];
  }

  c16_b_sampleTime = c16_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c16_debug_family_names,
    c16_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_nargin, 0U, c16_c_sf_marshallOut,
    c16_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_nargout, 1U, c16_c_sf_marshallOut,
    c16_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c16_u, 2U, c16_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_b_sampleTime, 3U,
    c16_c_sf_marshallOut, c16_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c16_y, 4U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(chartInstance->c16_states, 5U,
    c16_sf_marshallOut, c16_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, 3);
  _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, 4);
  if (CV_EML_IF(0, 1, 0, !chartInstance->c16_states_not_empty)) {
    _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, 5);
    for (c16_i7 = 0; c16_i7 < 3; c16_i7++) {
      chartInstance->c16_states[c16_i7] = 0.0;
    }

    chartInstance->c16_states_not_empty = TRUE;
  }

  _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, 7);
  for (c16_i8 = 0; c16_i8 < 3; c16_i8++) {
    c16_b_hoistedGlobal[c16_i8] = chartInstance->c16_states[c16_i8];
  }

  for (c16_i9 = 0; c16_i9 < 3; c16_i9++) {
    c16_a[c16_i9] = c16_u[c16_i9];
  }

  c16_b = c16_b_sampleTime;
  for (c16_i10 = 0; c16_i10 < 3; c16_i10++) {
    c16_a[c16_i10] *= c16_b;
  }

  for (c16_i11 = 0; c16_i11 < 3; c16_i11++) {
    c16_y[c16_i11] = c16_b_hoistedGlobal[c16_i11] + c16_a[c16_i11];
  }

  _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, 8);
  for (c16_i12 = 0; c16_i12 < 3; c16_i12++) {
    chartInstance->c16_states[c16_i12] = c16_y[c16_i12];
  }

  _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, -8);
  _SFD_SYMBOL_SCOPE_POP();
  for (c16_i13 = 0; c16_i13 < 3; c16_i13++) {
    (*c16_b_y)[c16_i13] = c16_y[c16_i13];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 3U, chartInstance->c16_sfEvent);
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_my_systemMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void initSimStructsc16_my_system(SFc16_my_systemInstanceStruct
  *chartInstance)
{
}

static void registerMessagesc16_my_system(SFc16_my_systemInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c16_machineNumber, uint32_T
  c16_chartNumber)
{
}

static const mxArray *c16_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData)
{
  const mxArray *c16_mxArrayOutData = NULL;
  int32_T c16_i14;
  real_T c16_b_inData[3];
  int32_T c16_i15;
  real_T c16_u[3];
  const mxArray *c16_y = NULL;
  SFc16_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc16_my_systemInstanceStruct *)chartInstanceVoid;
  c16_mxArrayOutData = NULL;
  for (c16_i14 = 0; c16_i14 < 3; c16_i14++) {
    c16_b_inData[c16_i14] = (*(real_T (*)[3])c16_inData)[c16_i14];
  }

  for (c16_i15 = 0; c16_i15 < 3; c16_i15++) {
    c16_u[c16_i15] = c16_b_inData[c16_i15];
  }

  c16_y = NULL;
  if (!chartInstance->c16_states_not_empty) {
    sf_mex_assign(&c16_y, sf_mex_create("y", NULL, 0, 0U, 1U, 0U, 2, 0, 0),
                  FALSE);
  } else {
    sf_mex_assign(&c16_y, sf_mex_create("y", c16_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  }

  sf_mex_assign(&c16_mxArrayOutData, c16_y, FALSE);
  return c16_mxArrayOutData;
}

static void c16_emlrt_marshallIn(SFc16_my_systemInstanceStruct *chartInstance,
  const mxArray *c16_b_states, const char_T *c16_identifier, real_T c16_y[3])
{
  emlrtMsgIdentifier c16_thisId;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c16_b_states), &c16_thisId,
    c16_y);
  sf_mex_destroy(&c16_b_states);
}

static void c16_b_emlrt_marshallIn(SFc16_my_systemInstanceStruct *chartInstance,
  const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId, real_T c16_y[3])
{
  real_T c16_dv2[3];
  int32_T c16_i16;
  if (mxIsEmpty(c16_u)) {
    chartInstance->c16_states_not_empty = FALSE;
  } else {
    chartInstance->c16_states_not_empty = TRUE;
    sf_mex_import(c16_parentId, sf_mex_dup(c16_u), c16_dv2, 1, 0, 0U, 1, 0U, 1,
                  3);
    for (c16_i16 = 0; c16_i16 < 3; c16_i16++) {
      c16_y[c16_i16] = c16_dv2[c16_i16];
    }
  }

  sf_mex_destroy(&c16_u);
}

static void c16_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData)
{
  const mxArray *c16_b_states;
  const char_T *c16_identifier;
  emlrtMsgIdentifier c16_thisId;
  real_T c16_y[3];
  int32_T c16_i17;
  SFc16_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc16_my_systemInstanceStruct *)chartInstanceVoid;
  c16_b_states = sf_mex_dup(c16_mxArrayInData);
  c16_identifier = c16_varName;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c16_b_states), &c16_thisId,
    c16_y);
  sf_mex_destroy(&c16_b_states);
  for (c16_i17 = 0; c16_i17 < 3; c16_i17++) {
    (*(real_T (*)[3])c16_outData)[c16_i17] = c16_y[c16_i17];
  }

  sf_mex_destroy(&c16_mxArrayInData);
}

static const mxArray *c16_b_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData)
{
  const mxArray *c16_mxArrayOutData = NULL;
  int32_T c16_i18;
  real_T c16_b_inData[3];
  int32_T c16_i19;
  real_T c16_u[3];
  const mxArray *c16_y = NULL;
  SFc16_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc16_my_systemInstanceStruct *)chartInstanceVoid;
  c16_mxArrayOutData = NULL;
  for (c16_i18 = 0; c16_i18 < 3; c16_i18++) {
    c16_b_inData[c16_i18] = (*(real_T (*)[3])c16_inData)[c16_i18];
  }

  for (c16_i19 = 0; c16_i19 < 3; c16_i19++) {
    c16_u[c16_i19] = c16_b_inData[c16_i19];
  }

  c16_y = NULL;
  sf_mex_assign(&c16_y, sf_mex_create("y", c16_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c16_mxArrayOutData, c16_y, FALSE);
  return c16_mxArrayOutData;
}

static void c16_c_emlrt_marshallIn(SFc16_my_systemInstanceStruct *chartInstance,
  const mxArray *c16_y, const char_T *c16_identifier, real_T c16_b_y[3])
{
  emlrtMsgIdentifier c16_thisId;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c16_y), &c16_thisId, c16_b_y);
  sf_mex_destroy(&c16_y);
}

static void c16_d_emlrt_marshallIn(SFc16_my_systemInstanceStruct *chartInstance,
  const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId, real_T c16_y[3])
{
  real_T c16_dv3[3];
  int32_T c16_i20;
  sf_mex_import(c16_parentId, sf_mex_dup(c16_u), c16_dv3, 1, 0, 0U, 1, 0U, 1, 3);
  for (c16_i20 = 0; c16_i20 < 3; c16_i20++) {
    c16_y[c16_i20] = c16_dv3[c16_i20];
  }

  sf_mex_destroy(&c16_u);
}

static void c16_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData)
{
  const mxArray *c16_y;
  const char_T *c16_identifier;
  emlrtMsgIdentifier c16_thisId;
  real_T c16_b_y[3];
  int32_T c16_i21;
  SFc16_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc16_my_systemInstanceStruct *)chartInstanceVoid;
  c16_y = sf_mex_dup(c16_mxArrayInData);
  c16_identifier = c16_varName;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c16_y), &c16_thisId, c16_b_y);
  sf_mex_destroy(&c16_y);
  for (c16_i21 = 0; c16_i21 < 3; c16_i21++) {
    (*(real_T (*)[3])c16_outData)[c16_i21] = c16_b_y[c16_i21];
  }

  sf_mex_destroy(&c16_mxArrayInData);
}

static const mxArray *c16_c_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData)
{
  const mxArray *c16_mxArrayOutData = NULL;
  real_T c16_u;
  const mxArray *c16_y = NULL;
  SFc16_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc16_my_systemInstanceStruct *)chartInstanceVoid;
  c16_mxArrayOutData = NULL;
  c16_u = *(real_T *)c16_inData;
  c16_y = NULL;
  sf_mex_assign(&c16_y, sf_mex_create("y", &c16_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c16_mxArrayOutData, c16_y, FALSE);
  return c16_mxArrayOutData;
}

static real_T c16_e_emlrt_marshallIn(SFc16_my_systemInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId)
{
  real_T c16_y;
  real_T c16_d1;
  sf_mex_import(c16_parentId, sf_mex_dup(c16_u), &c16_d1, 1, 0, 0U, 0, 0U, 0);
  c16_y = c16_d1;
  sf_mex_destroy(&c16_u);
  return c16_y;
}

static void c16_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData)
{
  const mxArray *c16_b_sampleTime;
  const char_T *c16_identifier;
  emlrtMsgIdentifier c16_thisId;
  real_T c16_y;
  SFc16_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc16_my_systemInstanceStruct *)chartInstanceVoid;
  c16_b_sampleTime = sf_mex_dup(c16_mxArrayInData);
  c16_identifier = c16_varName;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_y = c16_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c16_b_sampleTime),
    &c16_thisId);
  sf_mex_destroy(&c16_b_sampleTime);
  *(real_T *)c16_outData = c16_y;
  sf_mex_destroy(&c16_mxArrayInData);
}

const mxArray *sf_c16_my_system_get_eml_resolved_functions_info(void)
{
  const mxArray *c16_nameCaptureInfo;
  c16_ResolvedFunctionInfo c16_info[1];
  c16_ResolvedFunctionInfo (*c16_b_info)[1];
  const mxArray *c16_m0 = NULL;
  int32_T c16_i22;
  c16_ResolvedFunctionInfo *c16_r0;
  c16_nameCaptureInfo = NULL;
  c16_nameCaptureInfo = NULL;
  c16_b_info = (c16_ResolvedFunctionInfo (*)[1])c16_info;
  (*c16_b_info)[0].context = "";
  (*c16_b_info)[0].name = "mtimes";
  (*c16_b_info)[0].dominantType = "double";
  (*c16_b_info)[0].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  (*c16_b_info)[0].fileTimeLo = 1289516092U;
  (*c16_b_info)[0].fileTimeHi = 0U;
  (*c16_b_info)[0].mFileTimeLo = 0U;
  (*c16_b_info)[0].mFileTimeHi = 0U;
  sf_mex_assign(&c16_m0, sf_mex_createstruct("nameCaptureInfo", 1, 1), FALSE);
  for (c16_i22 = 0; c16_i22 < 1; c16_i22++) {
    c16_r0 = &c16_info[c16_i22];
    sf_mex_addfield(c16_m0, sf_mex_create("nameCaptureInfo", c16_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c16_r0->context)), "context", "nameCaptureInfo",
                    c16_i22);
    sf_mex_addfield(c16_m0, sf_mex_create("nameCaptureInfo", c16_r0->name, 15,
      0U, 0U, 0U, 2, 1, strlen(c16_r0->name)), "name", "nameCaptureInfo",
                    c16_i22);
    sf_mex_addfield(c16_m0, sf_mex_create("nameCaptureInfo",
      c16_r0->dominantType, 15, 0U, 0U, 0U, 2, 1, strlen(c16_r0->dominantType)),
                    "dominantType", "nameCaptureInfo", c16_i22);
    sf_mex_addfield(c16_m0, sf_mex_create("nameCaptureInfo", c16_r0->resolved,
      15, 0U, 0U, 0U, 2, 1, strlen(c16_r0->resolved)), "resolved",
                    "nameCaptureInfo", c16_i22);
    sf_mex_addfield(c16_m0, sf_mex_create("nameCaptureInfo", &c16_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c16_i22);
    sf_mex_addfield(c16_m0, sf_mex_create("nameCaptureInfo", &c16_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c16_i22);
    sf_mex_addfield(c16_m0, sf_mex_create("nameCaptureInfo",
      &c16_r0->mFileTimeLo, 7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo",
                    c16_i22);
    sf_mex_addfield(c16_m0, sf_mex_create("nameCaptureInfo",
      &c16_r0->mFileTimeHi, 7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo",
                    c16_i22);
  }

  sf_mex_assign(&c16_nameCaptureInfo, c16_m0, FALSE);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c16_nameCaptureInfo);
  return c16_nameCaptureInfo;
}

static const mxArray *c16_d_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData)
{
  const mxArray *c16_mxArrayOutData = NULL;
  int32_T c16_u;
  const mxArray *c16_y = NULL;
  SFc16_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc16_my_systemInstanceStruct *)chartInstanceVoid;
  c16_mxArrayOutData = NULL;
  c16_u = *(int32_T *)c16_inData;
  c16_y = NULL;
  sf_mex_assign(&c16_y, sf_mex_create("y", &c16_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c16_mxArrayOutData, c16_y, FALSE);
  return c16_mxArrayOutData;
}

static int32_T c16_f_emlrt_marshallIn(SFc16_my_systemInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId)
{
  int32_T c16_y;
  int32_T c16_i23;
  sf_mex_import(c16_parentId, sf_mex_dup(c16_u), &c16_i23, 1, 6, 0U, 0, 0U, 0);
  c16_y = c16_i23;
  sf_mex_destroy(&c16_u);
  return c16_y;
}

static void c16_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData)
{
  const mxArray *c16_b_sfEvent;
  const char_T *c16_identifier;
  emlrtMsgIdentifier c16_thisId;
  int32_T c16_y;
  SFc16_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc16_my_systemInstanceStruct *)chartInstanceVoid;
  c16_b_sfEvent = sf_mex_dup(c16_mxArrayInData);
  c16_identifier = c16_varName;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_y = c16_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c16_b_sfEvent),
    &c16_thisId);
  sf_mex_destroy(&c16_b_sfEvent);
  *(int32_T *)c16_outData = c16_y;
  sf_mex_destroy(&c16_mxArrayInData);
}

static uint8_T c16_g_emlrt_marshallIn(SFc16_my_systemInstanceStruct
  *chartInstance, const mxArray *c16_b_is_active_c16_my_system, const char_T
  *c16_identifier)
{
  uint8_T c16_y;
  emlrtMsgIdentifier c16_thisId;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_y = c16_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c16_b_is_active_c16_my_system), &c16_thisId);
  sf_mex_destroy(&c16_b_is_active_c16_my_system);
  return c16_y;
}

static uint8_T c16_h_emlrt_marshallIn(SFc16_my_systemInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId)
{
  uint8_T c16_y;
  uint8_T c16_u0;
  sf_mex_import(c16_parentId, sf_mex_dup(c16_u), &c16_u0, 1, 3, 0U, 0, 0U, 0);
  c16_y = c16_u0;
  sf_mex_destroy(&c16_u);
  return c16_y;
}

static void init_dsm_address_info(SFc16_my_systemInstanceStruct *chartInstance)
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

void sf_c16_my_system_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(267966996U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2980175585U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3222015561U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(365266718U);
}

mxArray *sf_c16_my_system_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("U0WLiAL691ISwMrbvrJdKH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

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
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c16_my_system_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

static const mxArray *sf_get_sim_state_info_c16_my_system(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[5],T\"y\",},{M[4],M[0],T\"states\",S'l','i','p'{{M1x2[52 58],M[0],}}},{M[8],M[0],T\"is_active_c16_my_system\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c16_my_system_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc16_my_systemInstanceStruct *chartInstance;
    chartInstance = (SFc16_my_systemInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _my_systemMachineNumber_,
           16,
           1,
           1,
           3,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"u");
          _SFD_SET_DATA_PROPS(1,2,0,1,"y");
          _SFD_SET_DATA_PROPS(2,10,0,0,"sampleTime");
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
        _SFD_CV_INIT_EML(0,1,1,1,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,137);
        _SFD_CV_INIT_EML_IF(0,1,0,59,77,-1,102);
        _SFD_TRANS_COV_WTS(0,0,0,1,0);
        if (chartAlreadyPresent==0) {
          _SFD_TRANS_COV_MAPS(0,
                              0,NULL,NULL,
                              0,NULL,NULL,
                              1,NULL,NULL,
                              0,NULL,NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)
            c16_b_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_c_sf_marshallOut,(MexInFcnForType)
          c16_c_sf_marshallIn);

        {
          real_T (*c16_u)[3];
          real_T (*c16_y)[3];
          c16_y = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
          c16_u = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c16_u);
          _SFD_SET_DATA_VALUE_PTR(1U, *c16_y);
          _SFD_SET_DATA_VALUE_PTR(2U, &chartInstance->c16_sampleTime);
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
  return "3YRugGLz9U4vcTbxtI8cKG";
}

static void sf_opaque_initialize_c16_my_system(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc16_my_systemInstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c16_my_system((SFc16_my_systemInstanceStruct*)
    chartInstanceVar);
  initialize_c16_my_system((SFc16_my_systemInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c16_my_system(void *chartInstanceVar)
{
  enable_c16_my_system((SFc16_my_systemInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c16_my_system(void *chartInstanceVar)
{
  disable_c16_my_system((SFc16_my_systemInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c16_my_system(void *chartInstanceVar)
{
  sf_c16_my_system((SFc16_my_systemInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c16_my_system(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c16_my_system
    ((SFc16_my_systemInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c16_my_system();/* state var info */
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

extern void sf_internal_set_sim_state_c16_my_system(SimStruct* S, const mxArray *
  st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c16_my_system();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c16_my_system((SFc16_my_systemInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c16_my_system(SimStruct* S)
{
  return sf_internal_get_sim_state_c16_my_system(S);
}

static void sf_opaque_set_sim_state_c16_my_system(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c16_my_system(S, st);
}

static void sf_opaque_terminate_c16_my_system(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc16_my_systemInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_my_system_optimization_info();
    }

    finalize_c16_my_system((SFc16_my_systemInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc16_my_system((SFc16_my_systemInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c16_my_system(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c16_my_system((SFc16_my_systemInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c16_my_system(SimStruct *S)
{
  /* Actual parameters from chart:
     sampleTime
   */
  const char_T *rtParamNames[] = { "sampleTime" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for sampleTime*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_my_system_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      16);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,16,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,16,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,16);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,16,1);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,16,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 1; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,16);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2516734265U));
  ssSetChecksum1(S,(2867571396U));
  ssSetChecksum2(S,(1994147733U));
  ssSetChecksum3(S,(2072590030U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c16_my_system(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c16_my_system(SimStruct *S)
{
  SFc16_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc16_my_systemInstanceStruct *)utMalloc(sizeof
    (SFc16_my_systemInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc16_my_systemInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c16_my_system;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c16_my_system;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c16_my_system;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c16_my_system;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c16_my_system;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c16_my_system;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c16_my_system;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c16_my_system;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c16_my_system;
  chartInstance->chartInfo.mdlStart = mdlStart_c16_my_system;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c16_my_system;
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

void c16_my_system_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c16_my_system(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c16_my_system(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c16_my_system(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c16_my_system_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
