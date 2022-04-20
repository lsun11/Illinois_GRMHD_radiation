#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * new_mask_excised_name;
  const char * new_mask_field_name;
  const char * new_style_excision_var;
  const char * old_style_excision_var;
  CCTK_INT allow_grid_resize;
  CCTK_INT cartoon_active;
  CCTK_INT eno_order;
  CCTK_INT new_excision;
  CCTK_INT old_excision;
  CCTK_INT order;
  CCTK_INT stencil;
  CCTK_INT verbose;
} PRIVATE_CARTOON2D_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARTOON2D_STRUCT_PARAMS \
  const char * const new_mask_excised_name = PRIVATE_CARTOON2D_STRUCT.new_mask_excised_name; \
  const char * const new_mask_field_name = PRIVATE_CARTOON2D_STRUCT.new_mask_field_name; \
  const char * const new_style_excision_var = PRIVATE_CARTOON2D_STRUCT.new_style_excision_var; \
  const char * const old_style_excision_var = PRIVATE_CARTOON2D_STRUCT.old_style_excision_var; \
  CCTK_INT const allow_grid_resize = PRIVATE_CARTOON2D_STRUCT.allow_grid_resize; \
  CCTK_INT const cartoon_active = PRIVATE_CARTOON2D_STRUCT.cartoon_active; \
  CCTK_INT const eno_order = PRIVATE_CARTOON2D_STRUCT.eno_order; \
  CCTK_INT const new_excision = PRIVATE_CARTOON2D_STRUCT.new_excision; \
  CCTK_INT const old_excision = PRIVATE_CARTOON2D_STRUCT.old_excision; \
  CCTK_INT const order = PRIVATE_CARTOON2D_STRUCT.order; \
  CCTK_INT const stencil = PRIVATE_CARTOON2D_STRUCT.stencil; \
  CCTK_INT const verbose = PRIVATE_CARTOON2D_STRUCT.verbose; \
  enum { \
      dummy_PRIVATE_CARTOON2D_STRUCT_new_mask_excised_name = sizeof( new_mask_excised_name ) \
    , dummy_PRIVATE_CARTOON2D_STRUCT_new_mask_field_name = sizeof( new_mask_field_name ) \
    , dummy_PRIVATE_CARTOON2D_STRUCT_new_style_excision_var = sizeof( new_style_excision_var ) \
    , dummy_PRIVATE_CARTOON2D_STRUCT_old_style_excision_var = sizeof( old_style_excision_var ) \
    , dummy_PRIVATE_CARTOON2D_STRUCT_allow_grid_resize = sizeof( allow_grid_resize ) \
    , dummy_PRIVATE_CARTOON2D_STRUCT_cartoon_active = sizeof( cartoon_active ) \
    , dummy_PRIVATE_CARTOON2D_STRUCT_eno_order = sizeof( eno_order ) \
    , dummy_PRIVATE_CARTOON2D_STRUCT_new_excision = sizeof( new_excision ) \
    , dummy_PRIVATE_CARTOON2D_STRUCT_old_excision = sizeof( old_excision ) \
    , dummy_PRIVATE_CARTOON2D_STRUCT_order = sizeof( order ) \
    , dummy_PRIVATE_CARTOON2D_STRUCT_stencil = sizeof( stencil ) \
    , dummy_PRIVATE_CARTOON2D_STRUCT_verbose = sizeof( verbose ) \
  }; \

