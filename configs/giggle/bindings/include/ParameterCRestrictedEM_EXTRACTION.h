#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL dR_EM;
  CCTK_REAL ddR_EM;
  CCTK_REAL phi_EM;
  CCTK_REAL radius_EM;
  CCTK_REAL radius_EM_phys;
  CCTK_REAL theta_EM;
} RESTRICTED_EM_EXTRACTION_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_EM_EXTRACTION_STRUCT_PARAMS \
  CCTK_REAL const dR_EM = RESTRICTED_EM_EXTRACTION_STRUCT.dR_EM; \
  CCTK_REAL const ddR_EM = RESTRICTED_EM_EXTRACTION_STRUCT.ddR_EM; \
  CCTK_REAL const phi_EM = RESTRICTED_EM_EXTRACTION_STRUCT.phi_EM; \
  CCTK_REAL const radius_EM = RESTRICTED_EM_EXTRACTION_STRUCT.radius_EM; \
  CCTK_REAL const radius_EM_phys = RESTRICTED_EM_EXTRACTION_STRUCT.radius_EM_phys; \
  CCTK_REAL const theta_EM = RESTRICTED_EM_EXTRACTION_STRUCT.theta_EM; \
  enum { \
      dummy_RESTRICTED_EM_EXTRACTION_STRUCT_dR_EM = sizeof( dR_EM ) \
    , dummy_RESTRICTED_EM_EXTRACTION_STRUCT_ddR_EM = sizeof( ddR_EM ) \
    , dummy_RESTRICTED_EM_EXTRACTION_STRUCT_phi_EM = sizeof( phi_EM ) \
    , dummy_RESTRICTED_EM_EXTRACTION_STRUCT_radius_EM = sizeof( radius_EM ) \
    , dummy_RESTRICTED_EM_EXTRACTION_STRUCT_radius_EM_phys = sizeof( radius_EM_phys ) \
    , dummy_RESTRICTED_EM_EXTRACTION_STRUCT_theta_EM = sizeof( theta_EM ) \
  }; \

