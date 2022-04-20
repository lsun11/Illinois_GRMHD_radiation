#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT MoL_Max_Evolved_Array_Size;
  CCTK_INT MoL_Max_Evolved_ComplexArray_Size;
  CCTK_INT MoL_Num_ArrayConstrained_Vars;
  CCTK_INT MoL_Num_ArrayEvolved_Vars;
  CCTK_INT MoL_Num_ArraySaveAndRestore_Vars;
  CCTK_INT MoL_Num_ComplexArrayConstrained_Vars;
  CCTK_INT MoL_Num_ComplexArrayEvolved_Vars;
  CCTK_INT MoL_Num_ComplexArraySaveAndRestore_Vars;
  CCTK_INT MoL_Num_ComplexConstrained_Vars;
  CCTK_INT MoL_Num_ComplexEvolved_Vars;
  CCTK_INT MoL_Num_ComplexSaveAndRestore_Vars;
  CCTK_INT MoL_Num_Constrained_Vars;
  CCTK_INT MoL_Num_Evolved_Vars;
  CCTK_INT MoL_Num_SaveAndRestore_Vars;
  CCTK_INT MoL_Num_Scratch_Levels;
} RESTRICTED_METHODOFLINES_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_METHODOFLINES_STRUCT_PARAMS \
  CCTK_INT const MoL_Max_Evolved_Array_Size = RESTRICTED_METHODOFLINES_STRUCT.MoL_Max_Evolved_Array_Size; \
  CCTK_INT const MoL_Max_Evolved_ComplexArray_Size = RESTRICTED_METHODOFLINES_STRUCT.MoL_Max_Evolved_ComplexArray_Size; \
  CCTK_INT const MoL_Num_ArrayConstrained_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_ArrayConstrained_Vars; \
  CCTK_INT const MoL_Num_ArrayEvolved_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_ArrayEvolved_Vars; \
  CCTK_INT const MoL_Num_ArraySaveAndRestore_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_ArraySaveAndRestore_Vars; \
  CCTK_INT const MoL_Num_ComplexArrayConstrained_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_ComplexArrayConstrained_Vars; \
  CCTK_INT const MoL_Num_ComplexArrayEvolved_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_ComplexArrayEvolved_Vars; \
  CCTK_INT const MoL_Num_ComplexArraySaveAndRestore_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_ComplexArraySaveAndRestore_Vars; \
  CCTK_INT const MoL_Num_ComplexConstrained_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_ComplexConstrained_Vars; \
  CCTK_INT const MoL_Num_ComplexEvolved_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_ComplexEvolved_Vars; \
  CCTK_INT const MoL_Num_ComplexSaveAndRestore_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_ComplexSaveAndRestore_Vars; \
  CCTK_INT const MoL_Num_Constrained_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_Constrained_Vars; \
  CCTK_INT const MoL_Num_Evolved_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_Evolved_Vars; \
  CCTK_INT const MoL_Num_SaveAndRestore_Vars = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_SaveAndRestore_Vars; \
  CCTK_INT const MoL_Num_Scratch_Levels = RESTRICTED_METHODOFLINES_STRUCT.MoL_Num_Scratch_Levels; \
  enum { \
      dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Max_Evolved_Array_Size = sizeof( MoL_Max_Evolved_Array_Size ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Max_Evolved_ComplexArray_Size = sizeof( MoL_Max_Evolved_ComplexArray_Size ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_ArrayConstrained_Vars = sizeof( MoL_Num_ArrayConstrained_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_ArrayEvolved_Vars = sizeof( MoL_Num_ArrayEvolved_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_ArraySaveAndRestore_Vars = sizeof( MoL_Num_ArraySaveAndRestore_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_ComplexArrayConstrained_Vars = sizeof( MoL_Num_ComplexArrayConstrained_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_ComplexArrayEvolved_Vars = sizeof( MoL_Num_ComplexArrayEvolved_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_ComplexArraySaveAndRestore_Vars = sizeof( MoL_Num_ComplexArraySaveAndRestore_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_ComplexConstrained_Vars = sizeof( MoL_Num_ComplexConstrained_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_ComplexEvolved_Vars = sizeof( MoL_Num_ComplexEvolved_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_ComplexSaveAndRestore_Vars = sizeof( MoL_Num_ComplexSaveAndRestore_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_Constrained_Vars = sizeof( MoL_Num_Constrained_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_Evolved_Vars = sizeof( MoL_Num_Evolved_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_SaveAndRestore_Vars = sizeof( MoL_Num_SaveAndRestore_Vars ) \
    , dummy_RESTRICTED_METHODOFLINES_STRUCT_MoL_Num_Scratch_Levels = sizeof( MoL_Num_Scratch_Levels ) \
  }; \

