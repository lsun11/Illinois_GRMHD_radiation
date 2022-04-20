#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL dR_GW;
  CCTK_REAL ddR_GW;
  CCTK_REAL phi_GW;
  CCTK_REAL radius_GW;
  CCTK_REAL radius_GW_phys;
  CCTK_REAL theta_GW;
} RESTRICTED_GW_EXTRACTION_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_GW_EXTRACTION_STRUCT_PARAMS \
  CCTK_REAL const dR_GW = RESTRICTED_GW_EXTRACTION_STRUCT.dR_GW; \
  CCTK_REAL const ddR_GW = RESTRICTED_GW_EXTRACTION_STRUCT.ddR_GW; \
  CCTK_REAL const phi_GW = RESTRICTED_GW_EXTRACTION_STRUCT.phi_GW; \
  CCTK_REAL const radius_GW = RESTRICTED_GW_EXTRACTION_STRUCT.radius_GW; \
  CCTK_REAL const radius_GW_phys = RESTRICTED_GW_EXTRACTION_STRUCT.radius_GW_phys; \
  CCTK_REAL const theta_GW = RESTRICTED_GW_EXTRACTION_STRUCT.theta_GW; \
  enum { \
      dummy_RESTRICTED_GW_EXTRACTION_STRUCT_dR_GW = sizeof( dR_GW ) \
    , dummy_RESTRICTED_GW_EXTRACTION_STRUCT_ddR_GW = sizeof( ddR_GW ) \
    , dummy_RESTRICTED_GW_EXTRACTION_STRUCT_phi_GW = sizeof( phi_GW ) \
    , dummy_RESTRICTED_GW_EXTRACTION_STRUCT_radius_GW = sizeof( radius_GW ) \
    , dummy_RESTRICTED_GW_EXTRACTION_STRUCT_radius_GW_phys = sizeof( radius_GW_phys ) \
    , dummy_RESTRICTED_GW_EXTRACTION_STRUCT_theta_GW = sizeof( theta_GW ) \
  }; \

