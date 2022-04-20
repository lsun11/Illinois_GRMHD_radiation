#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL ipoison;
  CCTK_REAL poison;
  CCTK_INT barriers;
  CCTK_INT check_tree_search;
  CCTK_INT tree_search;
} PRIVATE_CARPETINTERP_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARPETINTERP_STRUCT_PARAMS \
  CCTK_REAL const ipoison = PRIVATE_CARPETINTERP_STRUCT.ipoison; \
  CCTK_REAL const poison = PRIVATE_CARPETINTERP_STRUCT.poison; \
  CCTK_INT const barriers = PRIVATE_CARPETINTERP_STRUCT.barriers; \
  CCTK_INT const check_tree_search = PRIVATE_CARPETINTERP_STRUCT.check_tree_search; \
  CCTK_INT const tree_search = PRIVATE_CARPETINTERP_STRUCT.tree_search; \
  enum { \
      dummy_PRIVATE_CARPETINTERP_STRUCT_ipoison = sizeof( ipoison ) \
    , dummy_PRIVATE_CARPETINTERP_STRUCT_poison = sizeof( poison ) \
    , dummy_PRIVATE_CARPETINTERP_STRUCT_barriers = sizeof( barriers ) \
    , dummy_PRIVATE_CARPETINTERP_STRUCT_check_tree_search = sizeof( check_tree_search ) \
    , dummy_PRIVATE_CARPETINTERP_STRUCT_tree_search = sizeof( tree_search ) \
  }; \

