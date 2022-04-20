#include "h5_check.h"
#include "h5_error.h"
#include "h5_pline.h"

ck_err_t
h5checker_obj(char *fname, ck_addr_t obj_addr, int format_num, ck_errmsg_t *errbuf)
{
    driver_t 		*thefile=NULL;
    global_shared_t     *shared;
    ck_err_t		ret_value=SUCCEED;
    ck_err_t		ret;
    ck_addr_t           ss;

    g_obj_api = TRUE;
    g_obj_api_err = 0;

    g_format_num = format_num;
    if(g_format_num != FORMAT_ONE_SIX && g_format_num != FORMAT_ONE_EIGHT) {
        printf("Invalid library version provided.  Default library version is assumed.\n");
        g_format_num = DEFAULT_FORMAT;
    }

    g_obj_addr = obj_addr;

    if (g_obj_addr == CK_ADDR_UNDEF)
	printf("VALIDATING %s ", fname);
    else
	printf("VALIDATING %s at object header address %llu ", fname, g_obj_addr);

    if (g_format_num == FORMAT_ONE_SIX)
	printf("according to library release version 1.6.6...\n");
    else if (g_format_num == FORMAT_ONE_EIGHT)
	printf("according to library release version 1.8.0...\n");
    else
	printf("...invalid library release version...shouldn't happen.\n");
	

    shared = calloc(1, sizeof(global_shared_t));	
    if((thefile = (driver_t *)FD_open(fname, shared, SEC2_DRIVER)) == NULL) {
	error_push(ERR_FILE, ERR_NONE_SEC,
	    "Failure in opening input file using the default driver. Validation discontinued.", -1, NULL);
	++g_obj_api_err;
	goto done;
   }

   /* superblock validation has to be all passed before proceeding further */
   if(check_superblock(thefile) != SUCCEED) {
      error_push(ERR_LEV_0, ERR_NONE_SEC,
	    "Errors found when checking superblock. Validation stopped.", -1, NULL);
      ++g_obj_api_err;
      goto done;
    }

    /* not using the default driver */
    if(thefile->shared->driverid != SEC2_DRIVER) {
	if(FD_close(thefile) != SUCCEED) {
	    error_push(ERR_FILE, ERR_NONE_SEC, "Errors in closing input file using the default driver", -1, NULL);
    	    ++g_obj_api_err;
       }

       printf("Switching to new file driver...\n");
       if((thefile = (driver_t *)FD_open(fname, shared, shared->driverid)) == NULL) {
	    error_push(ERR_FILE, ERR_NONE_SEC,
		"Errors in opening input file. Validation stopped.", -1, NULL);
	    ++g_obj_api_err;
	    goto done;
	}
    }

    ss = FD_get_eof(thefile);
    if ((ss == CK_ADDR_UNDEF) || (ss < shared->stored_eoa)) {
	error_push(ERR_FILE, ERR_NONE_SEC,
	    "Invalid file size or file size less than superblock eoa. Validation stopped.", 
	    -1, NULL);
	++g_obj_api_err;
	goto done;
    }

    if((g_obj_addr != CK_ADDR_UNDEF) && (g_obj_addr >= shared->stored_eoa)) {
	error_push(ERR_FILE, ERR_NONE_SEC,
	    "Invalid Object header address provided. Validation stopped.", -1, NULL);
	++g_obj_api_err;
	goto done;
    }

    if(table_init(&obj_table) != SUCCEED) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Errors in initializing hard link table", -1, NULL);
	++g_obj_api_err;
    }

    if(pline_init_interface() != SUCCEED) {
        error_push(ERR_LEV_0, ERR_NONE_SEC, "Problems in initializing filters...later validation may be affected",
            -1, NULL);
	++g_obj_api_err;
    }   


    /* check the whole file if g_obj_addr is undefined */
    if(g_obj_addr == CK_ADDR_UNDEF) {
	ret = check_obj_header(thefile, shared->root_grp->header, NULL);
    } else
    	ret = check_obj_header(thefile, g_obj_addr, NULL);

    if(ret != SUCCEED)
	++g_obj_api_err;

done:
    if(thefile && FD_close(thefile) != SUCCEED) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Errors in closing the input file", -1, NULL);
	++g_obj_api_err;
    }

    if(shared)
	free(shared);

    if((errbuf != NULL) && (g_obj_api_err))
	process_errors(errbuf);

    return(g_obj_api_err? -1: 0);
}
