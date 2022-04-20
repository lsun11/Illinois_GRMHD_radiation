#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL cH1;
  CCTK_REAL chi_exponent;
  CCTK_REAL phi_cap;
  CCTK_INT Gammai_enforcement_type;
  CCTK_INT bc_type;
  CCTK_INT chi_evolution_enable;
  CCTK_INT chi_floor_enable;
  CCTK_INT enable_second_order_inside_ah;
  CCTK_INT gambound;
  CCTK_INT phi_cap_enable;
} PRIVATE_BSSN_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_BSSN_STRUCT_PARAMS \
  CCTK_REAL const cH1 = PRIVATE_BSSN_STRUCT.cH1; \
  CCTK_REAL const chi_exponent = PRIVATE_BSSN_STRUCT.chi_exponent; \
  CCTK_REAL const phi_cap = PRIVATE_BSSN_STRUCT.phi_cap; \
  CCTK_INT const Gammai_enforcement_type = PRIVATE_BSSN_STRUCT.Gammai_enforcement_type; \
  CCTK_INT const bc_type = PRIVATE_BSSN_STRUCT.bc_type; \
  CCTK_INT const chi_evolution_enable = PRIVATE_BSSN_STRUCT.chi_evolution_enable; \
  CCTK_INT const chi_floor_enable = PRIVATE_BSSN_STRUCT.chi_floor_enable; \
  CCTK_INT const enable_second_order_inside_ah = PRIVATE_BSSN_STRUCT.enable_second_order_inside_ah; \
  CCTK_INT const gambound = PRIVATE_BSSN_STRUCT.gambound; \
  CCTK_INT const phi_cap_enable = PRIVATE_BSSN_STRUCT.phi_cap_enable; \
  enum { \
      dummy_PRIVATE_BSSN_STRUCT_cH1 = sizeof( cH1 ) \
    , dummy_PRIVATE_BSSN_STRUCT_chi_exponent = sizeof( chi_exponent ) \
    , dummy_PRIVATE_BSSN_STRUCT_phi_cap = sizeof( phi_cap ) \
    , dummy_PRIVATE_BSSN_STRUCT_Gammai_enforcement_type = sizeof( Gammai_enforcement_type ) \
    , dummy_PRIVATE_BSSN_STRUCT_bc_type = sizeof( bc_type ) \
    , dummy_PRIVATE_BSSN_STRUCT_chi_evolution_enable = sizeof( chi_evolution_enable ) \
    , dummy_PRIVATE_BSSN_STRUCT_chi_floor_enable = sizeof( chi_floor_enable ) \
    , dummy_PRIVATE_BSSN_STRUCT_enable_second_order_inside_ah = sizeof( enable_second_order_inside_ah ) \
    , dummy_PRIVATE_BSSN_STRUCT_gambound = sizeof( gambound ) \
    , dummy_PRIVATE_BSSN_STRUCT_phi_cap_enable = sizeof( phi_cap_enable ) \
  }; \

