#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * initial_data_setup_method;
} RESTRICTED_INITBASE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_INITBASE_STRUCT_PARAMS \
  const char * const initial_data_setup_method = RESTRICTED_INITBASE_STRUCT.initial_data_setup_method; \
  enum { \
      dummy_RESTRICTED_INITBASE_STRUCT_initial_data_setup_method = sizeof( initial_data_setup_method ) \
  }; \

