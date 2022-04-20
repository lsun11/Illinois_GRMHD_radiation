#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL ICN_avg_theta;
  CCTK_REAL MoL_Tiny;
  CCTK_REAL RHS_error_weight;
  CCTK_REAL maximum_absolute_error;
  CCTK_REAL maximum_decrease;
  CCTK_REAL maximum_increase;
  CCTK_REAL maximum_relative_error;
  CCTK_REAL safety_factor;
  const char * Generic_Method_Descriptor;
  const char * Generic_Type;
  const char * ODE_Method;
  const char * verbose;
  CCTK_INT ICN_avg_swapped;
  CCTK_INT MoL_Intermediate_Steps;
  CCTK_INT MoL_Memory_Always_On;
  CCTK_INT MoL_NaN_Check;
  CCTK_INT adaptive_stepsize;
  CCTK_INT copy_ID_after_MoL_PostStep;
  CCTK_INT disable_prolongation;
  CCTK_INT initial_data_is_crap;
} PRIVATE_MOL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_MOL_STRUCT_PARAMS \
  CCTK_REAL const ICN_avg_theta = PRIVATE_MOL_STRUCT.ICN_avg_theta; \
  CCTK_REAL const MoL_Tiny = PRIVATE_MOL_STRUCT.MoL_Tiny; \
  CCTK_REAL const RHS_error_weight = PRIVATE_MOL_STRUCT.RHS_error_weight; \
  CCTK_REAL const maximum_absolute_error = PRIVATE_MOL_STRUCT.maximum_absolute_error; \
  CCTK_REAL const maximum_decrease = PRIVATE_MOL_STRUCT.maximum_decrease; \
  CCTK_REAL const maximum_increase = PRIVATE_MOL_STRUCT.maximum_increase; \
  CCTK_REAL const maximum_relative_error = PRIVATE_MOL_STRUCT.maximum_relative_error; \
  CCTK_REAL const safety_factor = PRIVATE_MOL_STRUCT.safety_factor; \
  const char * const Generic_Method_Descriptor = PRIVATE_MOL_STRUCT.Generic_Method_Descriptor; \
  const char * const Generic_Type = PRIVATE_MOL_STRUCT.Generic_Type; \
  const char * const ODE_Method = PRIVATE_MOL_STRUCT.ODE_Method; \
  const char * const verbose = PRIVATE_MOL_STRUCT.verbose; \
  CCTK_INT const ICN_avg_swapped = PRIVATE_MOL_STRUCT.ICN_avg_swapped; \
  CCTK_INT const MoL_Intermediate_Steps = PRIVATE_MOL_STRUCT.MoL_Intermediate_Steps; \
  CCTK_INT const MoL_Memory_Always_On = PRIVATE_MOL_STRUCT.MoL_Memory_Always_On; \
  CCTK_INT const MoL_NaN_Check = PRIVATE_MOL_STRUCT.MoL_NaN_Check; \
  CCTK_INT const adaptive_stepsize = PRIVATE_MOL_STRUCT.adaptive_stepsize; \
  CCTK_INT const copy_ID_after_MoL_PostStep = PRIVATE_MOL_STRUCT.copy_ID_after_MoL_PostStep; \
  CCTK_INT const disable_prolongation = PRIVATE_MOL_STRUCT.disable_prolongation; \
  CCTK_INT const initial_data_is_crap = PRIVATE_MOL_STRUCT.initial_data_is_crap; \
  enum { \
      dummy_PRIVATE_MOL_STRUCT_ICN_avg_theta = sizeof( ICN_avg_theta ) \
    , dummy_PRIVATE_MOL_STRUCT_MoL_Tiny = sizeof( MoL_Tiny ) \
    , dummy_PRIVATE_MOL_STRUCT_RHS_error_weight = sizeof( RHS_error_weight ) \
    , dummy_PRIVATE_MOL_STRUCT_maximum_absolute_error = sizeof( maximum_absolute_error ) \
    , dummy_PRIVATE_MOL_STRUCT_maximum_decrease = sizeof( maximum_decrease ) \
    , dummy_PRIVATE_MOL_STRUCT_maximum_increase = sizeof( maximum_increase ) \
    , dummy_PRIVATE_MOL_STRUCT_maximum_relative_error = sizeof( maximum_relative_error ) \
    , dummy_PRIVATE_MOL_STRUCT_safety_factor = sizeof( safety_factor ) \
    , dummy_PRIVATE_MOL_STRUCT_Generic_Method_Descriptor = sizeof( Generic_Method_Descriptor ) \
    , dummy_PRIVATE_MOL_STRUCT_Generic_Type = sizeof( Generic_Type ) \
    , dummy_PRIVATE_MOL_STRUCT_ODE_Method = sizeof( ODE_Method ) \
    , dummy_PRIVATE_MOL_STRUCT_verbose = sizeof( verbose ) \
    , dummy_PRIVATE_MOL_STRUCT_ICN_avg_swapped = sizeof( ICN_avg_swapped ) \
    , dummy_PRIVATE_MOL_STRUCT_MoL_Intermediate_Steps = sizeof( MoL_Intermediate_Steps ) \
    , dummy_PRIVATE_MOL_STRUCT_MoL_Memory_Always_On = sizeof( MoL_Memory_Always_On ) \
    , dummy_PRIVATE_MOL_STRUCT_MoL_NaN_Check = sizeof( MoL_NaN_Check ) \
    , dummy_PRIVATE_MOL_STRUCT_adaptive_stepsize = sizeof( adaptive_stepsize ) \
    , dummy_PRIVATE_MOL_STRUCT_copy_ID_after_MoL_PostStep = sizeof( copy_ID_after_MoL_PostStep ) \
    , dummy_PRIVATE_MOL_STRUCT_disable_prolongation = sizeof( disable_prolongation ) \
    , dummy_PRIVATE_MOL_STRUCT_initial_data_is_crap = sizeof( initial_data_is_crap ) \
  }; \

