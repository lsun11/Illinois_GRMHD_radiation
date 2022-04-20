#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL r_out_flux1;
  CCTK_REAL r_out_flux2;
  CCTK_REAL r_out_flux3;
  CCTK_REAL sam;
} RESTRICTED_MAG_BONDI_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_MAG_BONDI_STRUCT_PARAMS \
  CCTK_REAL const r_out_flux1 = RESTRICTED_MAG_BONDI_STRUCT.r_out_flux1; \
  CCTK_REAL const r_out_flux2 = RESTRICTED_MAG_BONDI_STRUCT.r_out_flux2; \
  CCTK_REAL const r_out_flux3 = RESTRICTED_MAG_BONDI_STRUCT.r_out_flux3; \
  CCTK_REAL const sam = RESTRICTED_MAG_BONDI_STRUCT.sam; \
  enum { \
      dummy_RESTRICTED_MAG_BONDI_STRUCT_r_out_flux1 = sizeof( r_out_flux1 ) \
    , dummy_RESTRICTED_MAG_BONDI_STRUCT_r_out_flux2 = sizeof( r_out_flux2 ) \
    , dummy_RESTRICTED_MAG_BONDI_STRUCT_r_out_flux3 = sizeof( r_out_flux3 ) \
    , dummy_RESTRICTED_MAG_BONDI_STRUCT_sam = sizeof( sam ) \
  }; \

