/* Include files */

#include <stddef.h>
#include "blas.h"
#include "my_system_sfun.h"
#include "c17_my_system.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "my_system_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c17_debug_family_names[13] = { "R1wrt0T", "R2wrt1T",
  "R2wrt0T", "nargin", "nargout", "th1", "th2", "ph1", "ph2", "l1", "l2", "p1",
  "p2" };

/* Function Declarations */
static void initialize_c17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance);
static void initialize_params_c17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance);
static void enable_c17_my_system(SFc17_my_systemInstanceStruct *chartInstance);
static void disable_c17_my_system(SFc17_my_systemInstanceStruct *chartInstance);
static void c17_update_debugger_state_c17_my_system
  (SFc17_my_systemInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c17_my_system(SFc17_my_systemInstanceStruct *
  chartInstance);
static void set_sim_state_c17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance, const mxArray *c17_st);
static void finalize_c17_my_system(SFc17_my_systemInstanceStruct *chartInstance);
static void sf_c17_my_system(SFc17_my_systemInstanceStruct *chartInstance);
static void c17_chartstep_c17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance);
static void initSimStructsc17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance);
static void registerMessagesc17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c17_machineNumber, uint32_T
  c17_chartNumber);
static const mxArray *c17_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static void c17_emlrt_marshallIn(SFc17_my_systemInstanceStruct *chartInstance,
  const mxArray *c17_p2, const char_T *c17_identifier, real_T c17_y[3]);
static void c17_b_emlrt_marshallIn(SFc17_my_systemInstanceStruct *chartInstance,
  const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId, real_T c17_y[3]);
static void c17_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static const mxArray *c17_b_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static real_T c17_c_emlrt_marshallIn(SFc17_my_systemInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId);
static void c17_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static const mxArray *c17_c_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static void c17_d_emlrt_marshallIn(SFc17_my_systemInstanceStruct *chartInstance,
  const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId, real_T c17_y[9]);
static void c17_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static void c17_info_helper(c17_ResolvedFunctionInfo c17_info[13]);
static void c17_eml_scalar_eg(SFc17_my_systemInstanceStruct *chartInstance);
static void c17_b_eml_scalar_eg(SFc17_my_systemInstanceStruct *chartInstance);
static const mxArray *c17_d_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static int32_T c17_e_emlrt_marshallIn(SFc17_my_systemInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId);
static void c17_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static uint8_T c17_f_emlrt_marshallIn(SFc17_my_systemInstanceStruct
  *chartInstance, const mxArray *c17_b_is_active_c17_my_system, const char_T
  *c17_identifier);
static uint8_T c17_g_emlrt_marshallIn(SFc17_my_systemInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId);
static void init_dsm_address_info(SFc17_my_systemInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance)
{
  chartInstance->c17_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c17_is_active_c17_my_system = 0U;
}

static void initialize_params_c17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance)
{
  real_T c17_d0;
  real_T c17_d1;
  sf_set_error_prefix_string(
    "Error evaluating data 'l1' in the parent workspace.\n");
  sf_mex_import_named("l1", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c17_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c17_l1 = c17_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'l2' in the parent workspace.\n");
  sf_mex_import_named("l2", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      &c17_d1, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c17_l2 = c17_d1;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c17_my_system(SFc17_my_systemInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c17_my_system(SFc17_my_systemInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c17_update_debugger_state_c17_my_system
  (SFc17_my_systemInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c17_my_system(SFc17_my_systemInstanceStruct *
  chartInstance)
{
  const mxArray *c17_st;
  const mxArray *c17_y = NULL;
  int32_T c17_i0;
  real_T c17_u[3];
  const mxArray *c17_b_y = NULL;
  int32_T c17_i1;
  real_T c17_b_u[3];
  const mxArray *c17_c_y = NULL;
  uint8_T c17_hoistedGlobal;
  uint8_T c17_c_u;
  const mxArray *c17_d_y = NULL;
  real_T (*c17_p2)[3];
  real_T (*c17_p1)[3];
  c17_p2 = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c17_p1 = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c17_st = NULL;
  c17_st = NULL;
  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_createcellarray(3), FALSE);
  for (c17_i0 = 0; c17_i0 < 3; c17_i0++) {
    c17_u[c17_i0] = (*c17_p1)[c17_i0];
  }

  c17_b_y = NULL;
  sf_mex_assign(&c17_b_y, sf_mex_create("y", c17_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_setcell(c17_y, 0, c17_b_y);
  for (c17_i1 = 0; c17_i1 < 3; c17_i1++) {
    c17_b_u[c17_i1] = (*c17_p2)[c17_i1];
  }

  c17_c_y = NULL;
  sf_mex_assign(&c17_c_y, sf_mex_create("y", c17_b_u, 0, 0U, 1U, 0U, 1, 3),
                FALSE);
  sf_mex_setcell(c17_y, 1, c17_c_y);
  c17_hoistedGlobal = chartInstance->c17_is_active_c17_my_system;
  c17_c_u = c17_hoistedGlobal;
  c17_d_y = NULL;
  sf_mex_assign(&c17_d_y, sf_mex_create("y", &c17_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c17_y, 2, c17_d_y);
  sf_mex_assign(&c17_st, c17_y, FALSE);
  return c17_st;
}

static void set_sim_state_c17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance, const mxArray *c17_st)
{
  const mxArray *c17_u;
  real_T c17_dv0[3];
  int32_T c17_i2;
  real_T c17_dv1[3];
  int32_T c17_i3;
  real_T (*c17_p1)[3];
  real_T (*c17_p2)[3];
  c17_p2 = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c17_p1 = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c17_doneDoubleBufferReInit = TRUE;
  c17_u = sf_mex_dup(c17_st);
  c17_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c17_u, 0)), "p1",
                       c17_dv0);
  for (c17_i2 = 0; c17_i2 < 3; c17_i2++) {
    (*c17_p1)[c17_i2] = c17_dv0[c17_i2];
  }

  c17_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c17_u, 1)), "p2",
                       c17_dv1);
  for (c17_i3 = 0; c17_i3 < 3; c17_i3++) {
    (*c17_p2)[c17_i3] = c17_dv1[c17_i3];
  }

  chartInstance->c17_is_active_c17_my_system = c17_f_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c17_u, 2)),
     "is_active_c17_my_system");
  sf_mex_destroy(&c17_u);
  c17_update_debugger_state_c17_my_system(chartInstance);
  sf_mex_destroy(&c17_st);
}

static void finalize_c17_my_system(SFc17_my_systemInstanceStruct *chartInstance)
{
}

static void sf_c17_my_system(SFc17_my_systemInstanceStruct *chartInstance)
{
  int32_T c17_i4;
  int32_T c17_i5;
  real_T *c17_th1;
  real_T *c17_th2;
  real_T *c17_ph1;
  real_T *c17_ph2;
  real_T (*c17_p2)[3];
  real_T (*c17_p1)[3];
  c17_p2 = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c17_p1 = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c17_ph2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c17_ph1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c17_th2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c17_th1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 4U, chartInstance->c17_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c17_th1, 0U);
  _SFD_DATA_RANGE_CHECK(*c17_th2, 1U);
  _SFD_DATA_RANGE_CHECK(*c17_ph1, 2U);
  _SFD_DATA_RANGE_CHECK(*c17_ph2, 3U);
  for (c17_i4 = 0; c17_i4 < 3; c17_i4++) {
    _SFD_DATA_RANGE_CHECK((*c17_p1)[c17_i4], 4U);
  }

  for (c17_i5 = 0; c17_i5 < 3; c17_i5++) {
    _SFD_DATA_RANGE_CHECK((*c17_p2)[c17_i5], 5U);
  }

  _SFD_DATA_RANGE_CHECK(chartInstance->c17_l1, 6U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c17_l2, 7U);
  chartInstance->c17_sfEvent = CALL_EVENT;
  c17_chartstep_c17_my_system(chartInstance);
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_my_systemMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c17_chartstep_c17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance)
{
  real_T c17_hoistedGlobal;
  real_T c17_b_hoistedGlobal;
  real_T c17_c_hoistedGlobal;
  real_T c17_d_hoistedGlobal;
  real_T c17_e_hoistedGlobal;
  real_T c17_f_hoistedGlobal;
  real_T c17_th1;
  real_T c17_th2;
  real_T c17_ph1;
  real_T c17_ph2;
  real_T c17_b_l1;
  real_T c17_b_l2;
  uint32_T c17_debug_family_var_map[13];
  real_T c17_R1wrt0T[9];
  real_T c17_R2wrt1T[9];
  real_T c17_R2wrt0T[9];
  real_T c17_nargin = 6.0;
  real_T c17_nargout = 2.0;
  real_T c17_p1[3];
  real_T c17_p2[3];
  real_T c17_x;
  real_T c17_b_x;
  real_T c17_c_x;
  real_T c17_d_x;
  real_T c17_a;
  real_T c17_b;
  real_T c17_y;
  real_T c17_e_x;
  real_T c17_f_x;
  real_T c17_g_x;
  real_T c17_h_x;
  real_T c17_b_a;
  real_T c17_b_b;
  real_T c17_b_y;
  real_T c17_i_x;
  real_T c17_j_x;
  real_T c17_k_x;
  real_T c17_l_x;
  real_T c17_m_x;
  real_T c17_n_x;
  real_T c17_o_x;
  real_T c17_p_x;
  real_T c17_q_x;
  real_T c17_r_x;
  real_T c17_c_a;
  real_T c17_c_b;
  real_T c17_c_y;
  real_T c17_s_x;
  real_T c17_t_x;
  real_T c17_u_x;
  real_T c17_v_x;
  real_T c17_d_a;
  real_T c17_d_b;
  real_T c17_d_y;
  real_T c17_w_x;
  real_T c17_x_x;
  real_T c17_y_x;
  real_T c17_ab_x;
  real_T c17_bb_x;
  real_T c17_cb_x;
  real_T c17_e_a;
  real_T c17_e_b;
  real_T c17_e_y;
  real_T c17_db_x;
  real_T c17_eb_x;
  real_T c17_fb_x;
  real_T c17_gb_x;
  real_T c17_f_a;
  real_T c17_f_b;
  real_T c17_f_y;
  real_T c17_hb_x;
  real_T c17_ib_x;
  real_T c17_jb_x;
  real_T c17_kb_x;
  real_T c17_lb_x;
  real_T c17_mb_x;
  real_T c17_nb_x;
  real_T c17_ob_x;
  real_T c17_pb_x;
  real_T c17_qb_x;
  real_T c17_g_a;
  real_T c17_g_b;
  real_T c17_g_y;
  real_T c17_rb_x;
  real_T c17_sb_x;
  real_T c17_tb_x;
  real_T c17_ub_x;
  real_T c17_h_a;
  real_T c17_h_b;
  real_T c17_h_y;
  real_T c17_vb_x;
  real_T c17_wb_x;
  int32_T c17_i6;
  real_T c17_i_a[9];
  int32_T c17_i7;
  real_T c17_i_b[9];
  int32_T c17_i8;
  int32_T c17_i9;
  int32_T c17_i10;
  real_T c17_C[9];
  int32_T c17_i11;
  int32_T c17_i12;
  int32_T c17_i13;
  int32_T c17_i14;
  int32_T c17_i15;
  int32_T c17_i16;
  int32_T c17_i17;
  int32_T c17_i18;
  int32_T c17_i19;
  real_T c17_j_b[3];
  int32_T c17_i20;
  int32_T c17_i21;
  int32_T c17_i22;
  real_T c17_b_C[3];
  int32_T c17_i23;
  int32_T c17_i24;
  int32_T c17_i25;
  int32_T c17_i26;
  int32_T c17_i27;
  int32_T c17_i28;
  int32_T c17_i29;
  int32_T c17_i30;
  int32_T c17_i31;
  int32_T c17_i32;
  int32_T c17_i33;
  int32_T c17_i34;
  real_T c17_i_y[3];
  int32_T c17_i35;
  int32_T c17_i36;
  int32_T c17_i37;
  int32_T c17_i38;
  int32_T c17_i39;
  real_T *c17_b_th1;
  real_T *c17_b_th2;
  real_T *c17_b_ph1;
  real_T *c17_b_ph2;
  real_T (*c17_b_p1)[3];
  real_T (*c17_b_p2)[3];
  c17_b_p2 = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c17_b_p1 = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c17_b_ph2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c17_b_ph1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c17_b_th2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c17_b_th1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 4U, chartInstance->c17_sfEvent);
  c17_hoistedGlobal = *c17_b_th1;
  c17_b_hoistedGlobal = *c17_b_th2;
  c17_c_hoistedGlobal = *c17_b_ph1;
  c17_d_hoistedGlobal = *c17_b_ph2;
  c17_e_hoistedGlobal = chartInstance->c17_l1;
  c17_f_hoistedGlobal = chartInstance->c17_l2;
  c17_th1 = c17_hoistedGlobal;
  c17_th2 = c17_b_hoistedGlobal;
  c17_ph1 = c17_c_hoistedGlobal;
  c17_ph2 = c17_d_hoistedGlobal;
  c17_b_l1 = c17_e_hoistedGlobal;
  c17_b_l2 = c17_f_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 13U, 13U, c17_debug_family_names,
    c17_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c17_R1wrt0T, 0U, c17_c_sf_marshallOut,
    c17_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c17_R2wrt1T, 1U, c17_c_sf_marshallOut,
    c17_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c17_R2wrt0T, 2U, c17_c_sf_marshallOut,
    c17_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_nargin, 3U, c17_b_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_nargout, 4U, c17_b_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c17_th1, 5U, c17_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c17_th2, 6U, c17_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c17_ph1, 7U, c17_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c17_ph2, 8U, c17_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_b_l1, 9U, c17_b_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_b_l2, 10U, c17_b_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c17_p1, 11U, c17_sf_marshallOut,
    c17_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c17_p2, 12U, c17_sf_marshallOut,
    c17_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 3);
  c17_x = c17_ph1;
  c17_b_x = c17_x;
  c17_b_x = muDoubleScalarCos(c17_b_x);
  c17_c_x = c17_th1;
  c17_d_x = c17_c_x;
  c17_d_x = muDoubleScalarCos(c17_d_x);
  c17_a = c17_b_x;
  c17_b = c17_d_x;
  c17_y = c17_a * c17_b;
  c17_e_x = c17_ph1;
  c17_f_x = c17_e_x;
  c17_f_x = muDoubleScalarSin(c17_f_x);
  c17_g_x = c17_th1;
  c17_h_x = c17_g_x;
  c17_h_x = muDoubleScalarCos(c17_h_x);
  c17_b_a = c17_f_x;
  c17_b_b = c17_h_x;
  c17_b_y = c17_b_a * c17_b_b;
  c17_i_x = c17_th1;
  c17_j_x = c17_i_x;
  c17_j_x = muDoubleScalarSin(c17_j_x);
  c17_k_x = c17_ph1;
  c17_l_x = c17_k_x;
  c17_l_x = muDoubleScalarSin(c17_l_x);
  c17_m_x = c17_ph1;
  c17_n_x = c17_m_x;
  c17_n_x = muDoubleScalarCos(c17_n_x);
  c17_o_x = c17_ph1;
  c17_p_x = c17_o_x;
  c17_p_x = muDoubleScalarCos(c17_p_x);
  c17_q_x = c17_th1;
  c17_r_x = c17_q_x;
  c17_r_x = muDoubleScalarSin(c17_r_x);
  c17_c_a = c17_p_x;
  c17_c_b = c17_r_x;
  c17_c_y = c17_c_a * c17_c_b;
  c17_s_x = c17_ph1;
  c17_t_x = c17_s_x;
  c17_t_x = muDoubleScalarSin(c17_t_x);
  c17_u_x = c17_th1;
  c17_v_x = c17_u_x;
  c17_v_x = muDoubleScalarSin(c17_v_x);
  c17_d_a = c17_t_x;
  c17_d_b = c17_v_x;
  c17_d_y = c17_d_a * c17_d_b;
  c17_w_x = c17_th1;
  c17_x_x = c17_w_x;
  c17_x_x = muDoubleScalarCos(c17_x_x);
  c17_R1wrt0T[0] = c17_y;
  c17_R1wrt0T[3] = c17_b_y;
  c17_R1wrt0T[6] = -c17_j_x;
  c17_R1wrt0T[1] = -c17_l_x;
  c17_R1wrt0T[4] = c17_n_x;
  c17_R1wrt0T[7] = 0.0;
  c17_R1wrt0T[2] = c17_c_y;
  c17_R1wrt0T[5] = c17_d_y;
  c17_R1wrt0T[8] = c17_x_x;
  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 6);
  c17_y_x = c17_ph2;
  c17_ab_x = c17_y_x;
  c17_ab_x = muDoubleScalarCos(c17_ab_x);
  c17_bb_x = c17_th2;
  c17_cb_x = c17_bb_x;
  c17_cb_x = muDoubleScalarCos(c17_cb_x);
  c17_e_a = c17_ab_x;
  c17_e_b = c17_cb_x;
  c17_e_y = c17_e_a * c17_e_b;
  c17_db_x = c17_ph2;
  c17_eb_x = c17_db_x;
  c17_eb_x = muDoubleScalarSin(c17_eb_x);
  c17_fb_x = c17_th2;
  c17_gb_x = c17_fb_x;
  c17_gb_x = muDoubleScalarCos(c17_gb_x);
  c17_f_a = c17_eb_x;
  c17_f_b = c17_gb_x;
  c17_f_y = c17_f_a * c17_f_b;
  c17_hb_x = c17_th2;
  c17_ib_x = c17_hb_x;
  c17_ib_x = muDoubleScalarSin(c17_ib_x);
  c17_jb_x = c17_ph2;
  c17_kb_x = c17_jb_x;
  c17_kb_x = muDoubleScalarSin(c17_kb_x);
  c17_lb_x = c17_ph2;
  c17_mb_x = c17_lb_x;
  c17_mb_x = muDoubleScalarCos(c17_mb_x);
  c17_nb_x = c17_ph2;
  c17_ob_x = c17_nb_x;
  c17_ob_x = muDoubleScalarCos(c17_ob_x);
  c17_pb_x = c17_th2;
  c17_qb_x = c17_pb_x;
  c17_qb_x = muDoubleScalarSin(c17_qb_x);
  c17_g_a = c17_ob_x;
  c17_g_b = c17_qb_x;
  c17_g_y = c17_g_a * c17_g_b;
  c17_rb_x = c17_ph2;
  c17_sb_x = c17_rb_x;
  c17_sb_x = muDoubleScalarSin(c17_sb_x);
  c17_tb_x = c17_th2;
  c17_ub_x = c17_tb_x;
  c17_ub_x = muDoubleScalarSin(c17_ub_x);
  c17_h_a = c17_sb_x;
  c17_h_b = c17_ub_x;
  c17_h_y = c17_h_a * c17_h_b;
  c17_vb_x = c17_th2;
  c17_wb_x = c17_vb_x;
  c17_wb_x = muDoubleScalarCos(c17_wb_x);
  c17_R2wrt1T[0] = c17_e_y;
  c17_R2wrt1T[3] = c17_f_y;
  c17_R2wrt1T[6] = -c17_ib_x;
  c17_R2wrt1T[1] = -c17_kb_x;
  c17_R2wrt1T[4] = c17_mb_x;
  c17_R2wrt1T[7] = 0.0;
  c17_R2wrt1T[2] = c17_g_y;
  c17_R2wrt1T[5] = c17_h_y;
  c17_R2wrt1T[8] = c17_wb_x;
  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 9);
  for (c17_i6 = 0; c17_i6 < 9; c17_i6++) {
    c17_i_a[c17_i6] = c17_R1wrt0T[c17_i6];
  }

  for (c17_i7 = 0; c17_i7 < 9; c17_i7++) {
    c17_i_b[c17_i7] = c17_R2wrt1T[c17_i7];
  }

  c17_eml_scalar_eg(chartInstance);
  c17_eml_scalar_eg(chartInstance);
  for (c17_i8 = 0; c17_i8 < 9; c17_i8++) {
    c17_R2wrt0T[c17_i8] = 0.0;
  }

  for (c17_i9 = 0; c17_i9 < 9; c17_i9++) {
    c17_R2wrt0T[c17_i9] = 0.0;
  }

  for (c17_i10 = 0; c17_i10 < 9; c17_i10++) {
    c17_C[c17_i10] = c17_R2wrt0T[c17_i10];
  }

  for (c17_i11 = 0; c17_i11 < 9; c17_i11++) {
    c17_R2wrt0T[c17_i11] = c17_C[c17_i11];
  }

  for (c17_i12 = 0; c17_i12 < 9; c17_i12++) {
    c17_C[c17_i12] = c17_R2wrt0T[c17_i12];
  }

  for (c17_i13 = 0; c17_i13 < 9; c17_i13++) {
    c17_R2wrt0T[c17_i13] = c17_C[c17_i13];
  }

  for (c17_i14 = 0; c17_i14 < 3; c17_i14++) {
    c17_i15 = 0;
    for (c17_i16 = 0; c17_i16 < 3; c17_i16++) {
      c17_R2wrt0T[c17_i15 + c17_i14] = 0.0;
      c17_i17 = 0;
      for (c17_i18 = 0; c17_i18 < 3; c17_i18++) {
        c17_R2wrt0T[c17_i15 + c17_i14] += c17_i_a[c17_i17 + c17_i14] *
          c17_i_b[c17_i18 + c17_i15];
        c17_i17 += 3;
      }

      c17_i15 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 13);
  for (c17_i19 = 0; c17_i19 < 9; c17_i19++) {
    c17_i_a[c17_i19] = c17_R1wrt0T[c17_i19];
  }

  c17_j_b[0] = c17_b_l1;
  c17_j_b[1] = 0.0;
  c17_j_b[2] = 0.0;
  c17_b_eml_scalar_eg(chartInstance);
  c17_b_eml_scalar_eg(chartInstance);
  for (c17_i20 = 0; c17_i20 < 3; c17_i20++) {
    c17_p1[c17_i20] = 0.0;
  }

  for (c17_i21 = 0; c17_i21 < 3; c17_i21++) {
    c17_p1[c17_i21] = 0.0;
  }

  for (c17_i22 = 0; c17_i22 < 3; c17_i22++) {
    c17_b_C[c17_i22] = c17_p1[c17_i22];
  }

  for (c17_i23 = 0; c17_i23 < 3; c17_i23++) {
    c17_p1[c17_i23] = c17_b_C[c17_i23];
  }

  for (c17_i24 = 0; c17_i24 < 3; c17_i24++) {
    c17_b_C[c17_i24] = c17_p1[c17_i24];
  }

  for (c17_i25 = 0; c17_i25 < 3; c17_i25++) {
    c17_p1[c17_i25] = c17_b_C[c17_i25];
  }

  for (c17_i26 = 0; c17_i26 < 3; c17_i26++) {
    c17_p1[c17_i26] = 0.0;
    c17_i27 = 0;
    for (c17_i28 = 0; c17_i28 < 3; c17_i28++) {
      c17_p1[c17_i26] += c17_i_a[c17_i27 + c17_i26] * c17_j_b[c17_i28];
      c17_i27 += 3;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 14);
  for (c17_i29 = 0; c17_i29 < 9; c17_i29++) {
    c17_i_a[c17_i29] = c17_R1wrt0T[c17_i29];
  }

  c17_j_b[0] = c17_b_l1;
  c17_j_b[1] = 0.0;
  c17_j_b[2] = 0.0;
  c17_b_eml_scalar_eg(chartInstance);
  c17_b_eml_scalar_eg(chartInstance);
  for (c17_i30 = 0; c17_i30 < 3; c17_i30++) {
    c17_b_C[c17_i30] = 0.0;
    c17_i31 = 0;
    for (c17_i32 = 0; c17_i32 < 3; c17_i32++) {
      c17_b_C[c17_i30] += c17_i_a[c17_i31 + c17_i30] * c17_j_b[c17_i32];
      c17_i31 += 3;
    }
  }

  for (c17_i33 = 0; c17_i33 < 9; c17_i33++) {
    c17_i_a[c17_i33] = c17_R2wrt0T[c17_i33];
  }

  c17_j_b[0] = c17_b_l2;
  c17_j_b[1] = 0.0;
  c17_j_b[2] = 0.0;
  c17_b_eml_scalar_eg(chartInstance);
  c17_b_eml_scalar_eg(chartInstance);
  for (c17_i34 = 0; c17_i34 < 3; c17_i34++) {
    c17_i_y[c17_i34] = 0.0;
    c17_i35 = 0;
    for (c17_i36 = 0; c17_i36 < 3; c17_i36++) {
      c17_i_y[c17_i34] += c17_i_a[c17_i35 + c17_i34] * c17_j_b[c17_i36];
      c17_i35 += 3;
    }
  }

  for (c17_i37 = 0; c17_i37 < 3; c17_i37++) {
    c17_p2[c17_i37] = c17_b_C[c17_i37] + c17_i_y[c17_i37];
  }

  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, -14);
  _SFD_SYMBOL_SCOPE_POP();
  for (c17_i38 = 0; c17_i38 < 3; c17_i38++) {
    (*c17_b_p1)[c17_i38] = c17_p1[c17_i38];
  }

  for (c17_i39 = 0; c17_i39 < 3; c17_i39++) {
    (*c17_b_p2)[c17_i39] = c17_p2[c17_i39];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 4U, chartInstance->c17_sfEvent);
}

static void initSimStructsc17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance)
{
}

static void registerMessagesc17_my_system(SFc17_my_systemInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c17_machineNumber, uint32_T
  c17_chartNumber)
{
}

static const mxArray *c17_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  int32_T c17_i40;
  real_T c17_b_inData[3];
  int32_T c17_i41;
  real_T c17_u[3];
  const mxArray *c17_y = NULL;
  SFc17_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc17_my_systemInstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  for (c17_i40 = 0; c17_i40 < 3; c17_i40++) {
    c17_b_inData[c17_i40] = (*(real_T (*)[3])c17_inData)[c17_i40];
  }

  for (c17_i41 = 0; c17_i41 < 3; c17_i41++) {
    c17_u[c17_i41] = c17_b_inData[c17_i41];
  }

  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", c17_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, FALSE);
  return c17_mxArrayOutData;
}

static void c17_emlrt_marshallIn(SFc17_my_systemInstanceStruct *chartInstance,
  const mxArray *c17_p2, const char_T *c17_identifier, real_T c17_y[3])
{
  emlrtMsgIdentifier c17_thisId;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_p2), &c17_thisId, c17_y);
  sf_mex_destroy(&c17_p2);
}

static void c17_b_emlrt_marshallIn(SFc17_my_systemInstanceStruct *chartInstance,
  const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId, real_T c17_y[3])
{
  real_T c17_dv2[3];
  int32_T c17_i42;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), c17_dv2, 1, 0, 0U, 1, 0U, 1, 3);
  for (c17_i42 = 0; c17_i42 < 3; c17_i42++) {
    c17_y[c17_i42] = c17_dv2[c17_i42];
  }

  sf_mex_destroy(&c17_u);
}

static void c17_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_p2;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  real_T c17_y[3];
  int32_T c17_i43;
  SFc17_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc17_my_systemInstanceStruct *)chartInstanceVoid;
  c17_p2 = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_p2), &c17_thisId, c17_y);
  sf_mex_destroy(&c17_p2);
  for (c17_i43 = 0; c17_i43 < 3; c17_i43++) {
    (*(real_T (*)[3])c17_outData)[c17_i43] = c17_y[c17_i43];
  }

  sf_mex_destroy(&c17_mxArrayInData);
}

static const mxArray *c17_b_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  real_T c17_u;
  const mxArray *c17_y = NULL;
  SFc17_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc17_my_systemInstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  c17_u = *(real_T *)c17_inData;
  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", &c17_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, FALSE);
  return c17_mxArrayOutData;
}

static real_T c17_c_emlrt_marshallIn(SFc17_my_systemInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId)
{
  real_T c17_y;
  real_T c17_d2;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), &c17_d2, 1, 0, 0U, 0, 0U, 0);
  c17_y = c17_d2;
  sf_mex_destroy(&c17_u);
  return c17_y;
}

static void c17_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_b_l2;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  real_T c17_y;
  SFc17_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc17_my_systemInstanceStruct *)chartInstanceVoid;
  c17_b_l2 = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_y = c17_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_b_l2),
    &c17_thisId);
  sf_mex_destroy(&c17_b_l2);
  *(real_T *)c17_outData = c17_y;
  sf_mex_destroy(&c17_mxArrayInData);
}

static const mxArray *c17_c_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  int32_T c17_i44;
  int32_T c17_i45;
  int32_T c17_i46;
  real_T c17_b_inData[9];
  int32_T c17_i47;
  int32_T c17_i48;
  int32_T c17_i49;
  real_T c17_u[9];
  const mxArray *c17_y = NULL;
  SFc17_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc17_my_systemInstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  c17_i44 = 0;
  for (c17_i45 = 0; c17_i45 < 3; c17_i45++) {
    for (c17_i46 = 0; c17_i46 < 3; c17_i46++) {
      c17_b_inData[c17_i46 + c17_i44] = (*(real_T (*)[9])c17_inData)[c17_i46 +
        c17_i44];
    }

    c17_i44 += 3;
  }

  c17_i47 = 0;
  for (c17_i48 = 0; c17_i48 < 3; c17_i48++) {
    for (c17_i49 = 0; c17_i49 < 3; c17_i49++) {
      c17_u[c17_i49 + c17_i47] = c17_b_inData[c17_i49 + c17_i47];
    }

    c17_i47 += 3;
  }

  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", c17_u, 0, 0U, 1U, 0U, 2, 3, 3), FALSE);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, FALSE);
  return c17_mxArrayOutData;
}

static void c17_d_emlrt_marshallIn(SFc17_my_systemInstanceStruct *chartInstance,
  const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId, real_T c17_y[9])
{
  real_T c17_dv3[9];
  int32_T c17_i50;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), c17_dv3, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c17_i50 = 0; c17_i50 < 9; c17_i50++) {
    c17_y[c17_i50] = c17_dv3[c17_i50];
  }

  sf_mex_destroy(&c17_u);
}

static void c17_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_R2wrt0T;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  real_T c17_y[9];
  int32_T c17_i51;
  int32_T c17_i52;
  int32_T c17_i53;
  SFc17_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc17_my_systemInstanceStruct *)chartInstanceVoid;
  c17_R2wrt0T = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_R2wrt0T), &c17_thisId,
    c17_y);
  sf_mex_destroy(&c17_R2wrt0T);
  c17_i51 = 0;
  for (c17_i52 = 0; c17_i52 < 3; c17_i52++) {
    for (c17_i53 = 0; c17_i53 < 3; c17_i53++) {
      (*(real_T (*)[9])c17_outData)[c17_i53 + c17_i51] = c17_y[c17_i53 + c17_i51];
    }

    c17_i51 += 3;
  }

  sf_mex_destroy(&c17_mxArrayInData);
}

const mxArray *sf_c17_my_system_get_eml_resolved_functions_info(void)
{
  const mxArray *c17_nameCaptureInfo;
  c17_ResolvedFunctionInfo c17_info[13];
  const mxArray *c17_m0 = NULL;
  int32_T c17_i54;
  c17_ResolvedFunctionInfo *c17_r0;
  c17_nameCaptureInfo = NULL;
  c17_nameCaptureInfo = NULL;
  c17_info_helper(c17_info);
  sf_mex_assign(&c17_m0, sf_mex_createstruct("nameCaptureInfo", 1, 13), FALSE);
  for (c17_i54 = 0; c17_i54 < 13; c17_i54++) {
    c17_r0 = &c17_info[c17_i54];
    sf_mex_addfield(c17_m0, sf_mex_create("nameCaptureInfo", c17_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c17_r0->context)), "context", "nameCaptureInfo",
                    c17_i54);
    sf_mex_addfield(c17_m0, sf_mex_create("nameCaptureInfo", c17_r0->name, 15,
      0U, 0U, 0U, 2, 1, strlen(c17_r0->name)), "name", "nameCaptureInfo",
                    c17_i54);
    sf_mex_addfield(c17_m0, sf_mex_create("nameCaptureInfo",
      c17_r0->dominantType, 15, 0U, 0U, 0U, 2, 1, strlen(c17_r0->dominantType)),
                    "dominantType", "nameCaptureInfo", c17_i54);
    sf_mex_addfield(c17_m0, sf_mex_create("nameCaptureInfo", c17_r0->resolved,
      15, 0U, 0U, 0U, 2, 1, strlen(c17_r0->resolved)), "resolved",
                    "nameCaptureInfo", c17_i54);
    sf_mex_addfield(c17_m0, sf_mex_create("nameCaptureInfo", &c17_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c17_i54);
    sf_mex_addfield(c17_m0, sf_mex_create("nameCaptureInfo", &c17_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c17_i54);
    sf_mex_addfield(c17_m0, sf_mex_create("nameCaptureInfo",
      &c17_r0->mFileTimeLo, 7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo",
                    c17_i54);
    sf_mex_addfield(c17_m0, sf_mex_create("nameCaptureInfo",
      &c17_r0->mFileTimeHi, 7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo",
                    c17_i54);
  }

  sf_mex_assign(&c17_nameCaptureInfo, c17_m0, FALSE);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c17_nameCaptureInfo);
  return c17_nameCaptureInfo;
}

static void c17_info_helper(c17_ResolvedFunctionInfo c17_info[13])
{
  c17_info[0].context = "";
  c17_info[0].name = "cos";
  c17_info[0].dominantType = "double";
  c17_info[0].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c17_info[0].fileTimeLo = 1343830372U;
  c17_info[0].fileTimeHi = 0U;
  c17_info[0].mFileTimeLo = 0U;
  c17_info[0].mFileTimeHi = 0U;
  c17_info[1].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  c17_info[1].name = "eml_scalar_cos";
  c17_info[1].dominantType = "double";
  c17_info[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m";
  c17_info[1].fileTimeLo = 1286818722U;
  c17_info[1].fileTimeHi = 0U;
  c17_info[1].mFileTimeLo = 0U;
  c17_info[1].mFileTimeHi = 0U;
  c17_info[2].context = "";
  c17_info[2].name = "mtimes";
  c17_info[2].dominantType = "double";
  c17_info[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c17_info[2].fileTimeLo = 1289516092U;
  c17_info[2].fileTimeHi = 0U;
  c17_info[2].mFileTimeLo = 0U;
  c17_info[2].mFileTimeHi = 0U;
  c17_info[3].context = "";
  c17_info[3].name = "sin";
  c17_info[3].dominantType = "double";
  c17_info[3].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c17_info[3].fileTimeLo = 1343830386U;
  c17_info[3].fileTimeHi = 0U;
  c17_info[3].mFileTimeLo = 0U;
  c17_info[3].mFileTimeHi = 0U;
  c17_info[4].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  c17_info[4].name = "eml_scalar_sin";
  c17_info[4].dominantType = "double";
  c17_info[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m";
  c17_info[4].fileTimeLo = 1286818736U;
  c17_info[4].fileTimeHi = 0U;
  c17_info[4].mFileTimeLo = 0U;
  c17_info[4].mFileTimeHi = 0U;
  c17_info[5].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c17_info[5].name = "eml_index_class";
  c17_info[5].dominantType = "";
  c17_info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c17_info[5].fileTimeLo = 1323166978U;
  c17_info[5].fileTimeHi = 0U;
  c17_info[5].mFileTimeLo = 0U;
  c17_info[5].mFileTimeHi = 0U;
  c17_info[6].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c17_info[6].name = "eml_scalar_eg";
  c17_info[6].dominantType = "double";
  c17_info[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c17_info[6].fileTimeLo = 1286818796U;
  c17_info[6].fileTimeHi = 0U;
  c17_info[6].mFileTimeLo = 0U;
  c17_info[6].mFileTimeHi = 0U;
  c17_info[7].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c17_info[7].name = "eml_xgemm";
  c17_info[7].dominantType = "char";
  c17_info[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c17_info[7].fileTimeLo = 1299073172U;
  c17_info[7].fileTimeHi = 0U;
  c17_info[7].mFileTimeLo = 0U;
  c17_info[7].mFileTimeHi = 0U;
  c17_info[8].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  c17_info[8].name = "eml_blas_inline";
  c17_info[8].dominantType = "";
  c17_info[8].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  c17_info[8].fileTimeLo = 1299073168U;
  c17_info[8].fileTimeHi = 0U;
  c17_info[8].mFileTimeLo = 0U;
  c17_info[8].mFileTimeHi = 0U;
  c17_info[9].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold";
  c17_info[9].name = "mtimes";
  c17_info[9].dominantType = "double";
  c17_info[9].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  c17_info[9].fileTimeLo = 1289516092U;
  c17_info[9].fileTimeHi = 0U;
  c17_info[9].mFileTimeLo = 0U;
  c17_info[9].mFileTimeHi = 0U;
  c17_info[10].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c17_info[10].name = "eml_index_class";
  c17_info[10].dominantType = "";
  c17_info[10].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  c17_info[10].fileTimeLo = 1323166978U;
  c17_info[10].fileTimeHi = 0U;
  c17_info[10].mFileTimeLo = 0U;
  c17_info[10].mFileTimeHi = 0U;
  c17_info[11].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c17_info[11].name = "eml_scalar_eg";
  c17_info[11].dominantType = "double";
  c17_info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  c17_info[11].fileTimeLo = 1286818796U;
  c17_info[11].fileTimeHi = 0U;
  c17_info[11].mFileTimeLo = 0U;
  c17_info[11].mFileTimeHi = 0U;
  c17_info[12].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  c17_info[12].name = "eml_refblas_xgemm";
  c17_info[12].dominantType = "char";
  c17_info[12].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m";
  c17_info[12].fileTimeLo = 1299073174U;
  c17_info[12].fileTimeHi = 0U;
  c17_info[12].mFileTimeLo = 0U;
  c17_info[12].mFileTimeHi = 0U;
}

static void c17_eml_scalar_eg(SFc17_my_systemInstanceStruct *chartInstance)
{
}

static void c17_b_eml_scalar_eg(SFc17_my_systemInstanceStruct *chartInstance)
{
}

static const mxArray *c17_d_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  int32_T c17_u;
  const mxArray *c17_y = NULL;
  SFc17_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc17_my_systemInstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  c17_u = *(int32_T *)c17_inData;
  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", &c17_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, FALSE);
  return c17_mxArrayOutData;
}

static int32_T c17_e_emlrt_marshallIn(SFc17_my_systemInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId)
{
  int32_T c17_y;
  int32_T c17_i55;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), &c17_i55, 1, 6, 0U, 0, 0U, 0);
  c17_y = c17_i55;
  sf_mex_destroy(&c17_u);
  return c17_y;
}

static void c17_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_b_sfEvent;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  int32_T c17_y;
  SFc17_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc17_my_systemInstanceStruct *)chartInstanceVoid;
  c17_b_sfEvent = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_y = c17_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_b_sfEvent),
    &c17_thisId);
  sf_mex_destroy(&c17_b_sfEvent);
  *(int32_T *)c17_outData = c17_y;
  sf_mex_destroy(&c17_mxArrayInData);
}

static uint8_T c17_f_emlrt_marshallIn(SFc17_my_systemInstanceStruct
  *chartInstance, const mxArray *c17_b_is_active_c17_my_system, const char_T
  *c17_identifier)
{
  uint8_T c17_y;
  emlrtMsgIdentifier c17_thisId;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_y = c17_g_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c17_b_is_active_c17_my_system), &c17_thisId);
  sf_mex_destroy(&c17_b_is_active_c17_my_system);
  return c17_y;
}

static uint8_T c17_g_emlrt_marshallIn(SFc17_my_systemInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId)
{
  uint8_T c17_y;
  uint8_T c17_u0;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), &c17_u0, 1, 3, 0U, 0, 0U, 0);
  c17_y = c17_u0;
  sf_mex_destroy(&c17_u);
  return c17_y;
}

static void init_dsm_address_info(SFc17_my_systemInstanceStruct *chartInstance)
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

void sf_c17_my_system_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2229334461U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3213330593U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1486039580U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2487768657U);
}

mxArray *sf_c17_my_system_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("mHwpEqptu722zW1gN1IX1B");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

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
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c17_my_system_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

static const mxArray *sf_get_sim_state_info_c17_my_system(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[10],T\"p1\",},{M[1],M[11],T\"p2\",},{M[8],M[0],T\"is_active_c17_my_system\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c17_my_system_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc17_my_systemInstanceStruct *chartInstance;
    chartInstance = (SFc17_my_systemInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _my_systemMachineNumber_,
           17,
           1,
           1,
           8,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"th1");
          _SFD_SET_DATA_PROPS(1,1,1,0,"th2");
          _SFD_SET_DATA_PROPS(2,1,1,0,"ph1");
          _SFD_SET_DATA_PROPS(3,1,1,0,"ph2");
          _SFD_SET_DATA_PROPS(4,2,0,1,"p1");
          _SFD_SET_DATA_PROPS(5,2,0,1,"p2");
          _SFD_SET_DATA_PROPS(6,10,0,0,"l1");
          _SFD_SET_DATA_PROPS(7,10,0,0,"l2");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,523);
        _SFD_TRANS_COV_WTS(0,0,0,1,0);
        if (chartAlreadyPresent==0) {
          _SFD_TRANS_COV_MAPS(0,
                              0,NULL,NULL,
                              0,NULL,NULL,
                              1,NULL,NULL,
                              0,NULL,NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c17_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c17_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c17_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c17_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c17_sf_marshallOut,(MexInFcnForType)
            c17_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c17_sf_marshallOut,(MexInFcnForType)
            c17_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c17_b_sf_marshallOut,(MexInFcnForType)
          c17_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c17_b_sf_marshallOut,(MexInFcnForType)
          c17_b_sf_marshallIn);

        {
          real_T *c17_th1;
          real_T *c17_th2;
          real_T *c17_ph1;
          real_T *c17_ph2;
          real_T (*c17_p1)[3];
          real_T (*c17_p2)[3];
          c17_p2 = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
          c17_p1 = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
          c17_ph2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c17_ph1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c17_th2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c17_th1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c17_th1);
          _SFD_SET_DATA_VALUE_PTR(1U, c17_th2);
          _SFD_SET_DATA_VALUE_PTR(2U, c17_ph1);
          _SFD_SET_DATA_VALUE_PTR(3U, c17_ph2);
          _SFD_SET_DATA_VALUE_PTR(4U, *c17_p1);
          _SFD_SET_DATA_VALUE_PTR(5U, *c17_p2);
          _SFD_SET_DATA_VALUE_PTR(6U, &chartInstance->c17_l1);
          _SFD_SET_DATA_VALUE_PTR(7U, &chartInstance->c17_l2);
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
  return "2ffaOmsPSm4ZfOuEwEkUBE";
}

static void sf_opaque_initialize_c17_my_system(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc17_my_systemInstanceStruct*) chartInstanceVar
    )->S,0);
  initialize_params_c17_my_system((SFc17_my_systemInstanceStruct*)
    chartInstanceVar);
  initialize_c17_my_system((SFc17_my_systemInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c17_my_system(void *chartInstanceVar)
{
  enable_c17_my_system((SFc17_my_systemInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c17_my_system(void *chartInstanceVar)
{
  disable_c17_my_system((SFc17_my_systemInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c17_my_system(void *chartInstanceVar)
{
  sf_c17_my_system((SFc17_my_systemInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c17_my_system(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c17_my_system
    ((SFc17_my_systemInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c17_my_system();/* state var info */
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

extern void sf_internal_set_sim_state_c17_my_system(SimStruct* S, const mxArray *
  st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c17_my_system();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c17_my_system((SFc17_my_systemInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c17_my_system(SimStruct* S)
{
  return sf_internal_get_sim_state_c17_my_system(S);
}

static void sf_opaque_set_sim_state_c17_my_system(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c17_my_system(S, st);
}

static void sf_opaque_terminate_c17_my_system(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc17_my_systemInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_my_system_optimization_info();
    }

    finalize_c17_my_system((SFc17_my_systemInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc17_my_system((SFc17_my_systemInstanceStruct*) chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c17_my_system(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c17_my_system((SFc17_my_systemInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c17_my_system(SimStruct *S)
{
  /* Actual parameters from chart:
     l1 l2
   */
  const char_T *rtParamNames[] = { "l1", "l2" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for l1*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);

  /* registration for l2*/
  ssRegDlgParamAsRunTimeParam(S, 1, 1, rtParamNames[1], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_my_system_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      17);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,17,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,17,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,17);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,17,4);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,17,2);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=2; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 4; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,17);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1160376076U));
  ssSetChecksum1(S,(556603698U));
  ssSetChecksum2(S,(1262913193U));
  ssSetChecksum3(S,(657210858U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c17_my_system(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c17_my_system(SimStruct *S)
{
  SFc17_my_systemInstanceStruct *chartInstance;
  chartInstance = (SFc17_my_systemInstanceStruct *)utMalloc(sizeof
    (SFc17_my_systemInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc17_my_systemInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c17_my_system;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c17_my_system;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c17_my_system;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c17_my_system;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c17_my_system;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c17_my_system;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c17_my_system;
  chartInstance->chartInfo.getSimStateInfo = sf_get_sim_state_info_c17_my_system;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c17_my_system;
  chartInstance->chartInfo.mdlStart = mdlStart_c17_my_system;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c17_my_system;
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

void c17_my_system_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c17_my_system(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c17_my_system(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c17_my_system(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c17_my_system_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
