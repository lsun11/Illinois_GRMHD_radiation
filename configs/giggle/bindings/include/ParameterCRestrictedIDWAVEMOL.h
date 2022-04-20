#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * initial_data;
} RESTRICTED_IDWAVEMOL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_IDWAVEMOL_STRUCT_PARAMS \
  const char * const initial_data = RESTRICTED_IDWAVEMOL_STRUCT.initial_data; \
  enum { \
      dummy_RESTRICTED_IDWAVEMOL_STRUCT_initial_data = sizeof( initial_data ) \
  }; \

