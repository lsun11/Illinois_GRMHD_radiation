#include "h5_check.h"
#include "h5_error.h"


/* for handling hard links */
static int	table_search(ck_addr_t);
static void	table_insert(ck_addr_t, int);

/* for handling symbol table names */
static int name_list_init(name_list_t **name_list);
static int name_list_search(name_list_t *nl, char *symname);
static int name_list_insert(name_list_t *nl, char *symname);
static int name_list_dest(name_list_t *nl);

/* Validation routines */
static ck_err_t check_btree(driver_t *, ck_addr_t, unsigned, uint8_t *, name_list_t *, void *, void *);
static ck_err_t check_sym(driver_t *, ck_addr_t, uint8_t *, name_list_t *);
static ck_err_t check_lheap(driver_t *, ck_addr_t, uint8_t **);

#ifdef NOTUSED
static ck_err_t check_gheap(driver_t *, ck_addr_t, uint8_t **);
#endif

static 	ck_addr_t locate_super_signature(driver_t *);
static  ck_err_t  gp_ent_decode(global_shared_t *, const uint8_t **, GP_entry_t *);
static  ck_err_t  gp_ent_decode_vec(global_shared_t *, const uint8_t **, GP_entry_t *, unsigned);
static	type_t    *type_alloc(ck_addr_t);

static  ck_err_t 	OBJ_alloc_msgs(OBJ_t *, ck_size_t);
static  int 		find_in_ohdr(driver_t *, OBJ_t *, int);
static	void 		*OBJ_shared_read(driver_t *, OBJ_shared_t *, const obj_class_t *);
static  void 		*OBJ_shared_decode(driver_t *, const uint8_t *, const obj_class_t *, const uint8_t *, ck_addr_t);
static	ck_err_t 	decode_validate_messages(driver_t *, OBJ_t *);


/*
 *  Virtual file drivers
 */
static 	void 		set_driver_id(int */*out*/, char *);
static 	driver_class_t 	*get_driver_class(int);
static 	void 		*get_driver_info(int, global_shared_t *);
static 	ck_err_t 	decode_driver(global_shared_t *, const uint8_t *);

static 	driver_t 	*sec2_open(const char *, global_shared_t *, int);
static 	ck_err_t 	sec2_read(driver_t *, ck_addr_t, size_t, void */*out*/);
static 	ck_err_t 	sec2_close(driver_t *);
static 	ck_addr_t 	sec2_get_eof(driver_t *);
static 	char 		*sec2_get_fname(driver_t *, ck_addr_t);

static const driver_class_t sec2_g = {
    "sec2",                                /* name                  */
    NULL,                                  /* decode_driver         */
    sec2_open,                             /* open                  */
    sec2_close,                            /* close                 */
    sec2_read,                             /* read                  */
    sec2_get_eof,                          /* get_eof               */
    sec2_get_fname,			   /* get file name	    */
};


static 	ck_err_t 	multi_decode_driver(global_shared_t *, const unsigned char *);
static 	driver_t 	*multi_open(const char *, global_shared_t *, int);
static 	ck_err_t 	multi_close(driver_t *);
static 	ck_err_t 	multi_read(driver_t *, ck_addr_t, size_t, void */*out*/);
static 	ck_addr_t 	multi_get_eof(driver_t *file);
static 	char 		*multi_get_fname(driver_t *file, ck_addr_t logi_addr);


static void 	set_multi_driver_properties(driver_multi_fapl_t **, driver_mem_t *, const char **, ck_addr_t *);
static ck_err_t	compute_next(driver_multi_t *file);
static ck_err_t open_members(driver_multi_t *file, global_shared_t *);

static const driver_class_t multi_g = {
    "multi",                               /* name                 */
    multi_decode_driver,                   /* decode_driver        */
    multi_open,                            /* open                  */
    multi_close,                           /* close                 */
    multi_read,                            /* read                  */
    multi_get_eof,                         /* get_eof               */
    multi_get_fname,			   /* get file name 	    */
};

static 	ck_err_t 	family_decode_driver(global_shared_t *, const unsigned char *);
static 	driver_t 	*family_open(const char *, global_shared_t *, int);
static 	ck_err_t 	family_close(driver_t *);
static 	ck_err_t 	family_read(driver_t *, ck_addr_t, size_t, void */*out*/);
static 	ck_addr_t 	family_get_eof(driver_t *file);
static 	char 		*family_get_fname(driver_t *file, ck_addr_t logi_addr);

static void set_fami_driver_properties(driver_fami_fapl_t **, driver_mem_t *, const char **, ck_addr_t *);

static const driver_class_t family_g = {
    "family",                                   /* name                  */
    family_decode_driver,                       /* decode_driver         */
    family_open,                           	/* open                  */
    family_close,                          	/* close                 */
    family_read,                           	/* read                  */
    family_get_eof,                        	/* get_eof               */
    family_get_fname,				/* get file name */
};



/*
 * Header messages
 */

/* NIL: 0x0000 */
static const obj_class_t OBJ_NIL[1] = {{
    OBJ_NIL_ID,             /* message id number             */
    NULL,                   /* no decode method              */
}};

static void *OBJ_sds_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Simple Database: 0x0001 */
static const obj_class_t OBJ_SDS[1] = {{
    OBJ_SDS_ID,             	/* message id number         */
    OBJ_sds_decode,         	/* decode message            */
}};


static void *OBJ_linfo_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);
static ck_err_t G_dense_ck_fh_msg_cb(driver_t *file, const void *_record, void *_ck_msg_udata);

/* Link Info: 0x0002 */
static const obj_class_t OBJ_LINFO[1] = {{
    OBJ_LINFO_ID,             	/* message id number         */
    OBJ_linfo_decode,         	/* decode message            */
}};

static void *OBJ_dt_decode (driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);
static ck_err_t OBJ_dt_decode_helper(driver_t *, const uint8_t **, type_t *, const uint8_t *, ck_addr_t);

/* Datatype: 0x0003 */
static const obj_class_t OBJ_DT[1] = {{
    OBJ_DT_ID,              /* message id number          */
    OBJ_dt_decode,           /* decode message            */
}};

static void  *OBJ_fill_old_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Data Storage - Fill Value (old): 0x0004 */
static const obj_class_t OBJ_FILL_OLD[1] = {{
    OBJ_FILL_OLD_ID,            /* message id number               */
    OBJ_fill_old_decode,        /* decode message                  */
}};

static void  *OBJ_fill_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Data Storage - Fill Value: 0x0005 */
static const obj_class_t OBJ_FILL[1] = {{
    OBJ_FILL_ID,            /* message id number               */
    OBJ_fill_decode,        /* decode message                  */
}};

static void  *OBJ_link_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Link Message: 0x0006 */
static const obj_class_t OBJ_LINK[1] = {{
    OBJ_LINK_ID,             	/* message id number         */
    OBJ_link_decode,         	/* decode message            */
}};

static void *OBJ_edf_decode(driver_t *, const uint8_t *p, const uint8_t *, ck_addr_t);

/* Data Storage - External Data Files: 0x0007 */
static const obj_class_t OBJ_EDF[1] = {{
    OBJ_EDF_ID,                 /* message id number             */
    OBJ_edf_decode,             /* decode message                */
}};


static void *OBJ_ginfo_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);
/*  Group Information: 0x000A */
static const obj_class_t OBJ_GINFO[1] = {{
    OBJ_GINFO_ID,               /* message id number             */
    OBJ_ginfo_decode,	        /* decode message                */
}};

static void *OBJ_layout_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Data Storage - layout: 0x0008 */
static const obj_class_t OBJ_LAYOUT[1] = {{
    OBJ_LAYOUT_ID,              /* message id number             */
    OBJ_layout_decode,	        /* decode message                */
}};

static void *OBJ_bogus_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Bogus Message 0x0009 */
static const obj_class_t OBJ_BOGUS[1] = {{
    OBJ_BOGUS_ID,              /* message id number             */
    OBJ_bogus_decode,	       /* decode message                */
}};

static void *OBJ_filter_decode (driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Data Storage - Filter pipeline: 0x000B */
static const obj_class_t OBJ_FILTER[1] = {{
    OBJ_FILTER_ID,              /* message id number            */
    OBJ_filter_decode,          /* decode message               */
}};

static void *OBJ_attr_decode (driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Attribute: 0x000C */
static const obj_class_t OBJ_ATTR[1] = {{
    OBJ_ATTR_ID,                 /* message id number           */
    OBJ_attr_decode,             /* decode message              */
}};


static void *OBJ_comm_decode(driver_t *, const uint8_t *p, const uint8_t *, ck_addr_t);

/* Object Comment: 0x000D */
static const obj_class_t OBJ_COMM[1] = {{
    OBJ_COMM_ID,                /* message id number             */
    OBJ_comm_decode,            /* decode message                */
}};


static void *OBJ_mdt_old_decode(driver_t *, const uint8_t *p, const uint8_t *, ck_addr_t);

/* Track whether tzset routine was called */
static int      ntzset=0;

/* Object Modification Date & Time (old): 0x000E */
static const obj_class_t OBJ_MDT_OLD[1] = {{
    OBJ_MDT_OLD_ID,                 /* message id number             */
    OBJ_mdt_old_decode,             /* decode message                */
}};


static void *OBJ_shmesg_decode (driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Shared Object Message: 0x000F */
static const obj_class_t OBJ_SHMESG[1] = {{
    OBJ_SHMESG_ID,              /* message id number             */
    OBJ_shmesg_decode,          /* decode method                 */
}};



static void *OBJ_cont_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Object Header Continuation: 0x0010 */
static const obj_class_t OBJ_CONT[1] = {{
    OBJ_CONT_ID,                /* message id number             */
    OBJ_cont_decode,             /* decode message               */
}};


static void *OBJ_group_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Symbol Table Message: 0x0011 */
static const obj_class_t OBJ_GROUP[1] = {{
    OBJ_GROUP_ID,              /* message id number             */
    OBJ_group_decode,          /* decode message                */
}};


static void *OBJ_mdt_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Object Modification Date & Time: 0x0012 */
static const obj_class_t OBJ_MDT[1] = {{
    OBJ_MDT_ID,           	/* message id number             */
    OBJ_mdt_decode,       	/* decode message                */
}};

static void *OBJ_btreek_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Non-default v1 B-tree 'K' values: 0x0013 */
static const obj_class_t OBJ_BTREEK[1] = {{
    OBJ_BTREEK_ID,           	/* message id number             */
    OBJ_btreek_decode,       	/* decode message                */
}};


static void *OBJ_drvinfo_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Driver Info settings: 0x0014 */
static const obj_class_t OBJ_DRVINFO[1] = {{
    OBJ_DRVINFO_ID,           	/* message id number             */
    OBJ_drvinfo_decode,       	/* decode message                */
}};


static void *OBJ_ainfo_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);
static ck_err_t A_dense_ck_fh_msg_cb(driver_t *file, const void *_record, void *_ck_msg_udata);

/* Attribute Information : 0x0015 */
static const obj_class_t OBJ_AINFO[1] = {{
    OBJ_AINFO_ID,           	/* message id number             */
    OBJ_ainfo_decode,       	/* decode message                */
}};

static void *OBJ_refcount_decode(driver_t *, const uint8_t *, const uint8_t *, ck_addr_t);

/* Object's Reference Count: 0x0016 */
static const obj_class_t OBJ_REFCOUNT[1] = {{
    OBJ_REFCOUNT_ID,           	/* message id number             */
    OBJ_refcount_decode,       	/* decode message                */
}};

/*  Unknown Message : 0x0017 */
static const obj_class_t OBJ_UNKNOWN[1] = {{
    OBJ_UNKNOWN_ID,             /* message id number             */
    NULL,                       /* decode message                */
}};


const obj_class_t *const message_type_g[] = {
    OBJ_NIL,           	/* 0x0000 NIL                                   */
    OBJ_SDS,        	/* 0x0001 Simple Dataspace			*/
    OBJ_LINFO,          /* 0x0002 Link Info Message			*/
    OBJ_DT,          	/* 0x0003 Datatype                              */
    OBJ_FILL_OLD,       /* 0x0004 Data Storage - Fill Value (Old) 	*/
    OBJ_FILL,      	/* 0x0005 Data storage - Fill Value         	*/
    OBJ_LINK,           /* 0x0006 LINK 					*/
    OBJ_EDF,           	/* 0x0007 Data Storage - External Data Files    */
    OBJ_LAYOUT,        	/* 0x0008 Data Storage - Layout                 */
    OBJ_BOGUS,          /* 0x0009 Bogus				        */
    OBJ_GINFO,          /* 0x000A Group information Message   		*/
    OBJ_FILTER,       	/* 0x000B Data storage - Filter Pipeline        */
    OBJ_ATTR,          	/* 0x000C Attribute 				*/
    OBJ_COMM,       	/* 0x000D Object Comment                        */
    OBJ_MDT_OLD,      	/* 0x000E Object Modification Date & Time (Old) */
    OBJ_SHMESG,        	/* 0x000F File-wide shared message table 	*/
    OBJ_CONT,          	/* 0x0010 Object Header Continuation            */
    OBJ_GROUP,          /* 0x0011 Symbol Table Message				*/
    OBJ_MDT,		/* 0x0012 Object modification Date & Time  	*/
    OBJ_BTREEK,		/* 0x0013 Non-default v1 B-tree 'K' values	*/
    OBJ_DRVINFO,	/* 0x0014 Driver Info settings			*/
    OBJ_AINFO,		/* 0x0015 Attribute information 		*/
    OBJ_REFCOUNT,	/* 0x0016 Object's ref. count			*/
    OBJ_UNKNOWN,	/* 0x0017 Placeholder for unknown message	*/
};

/* end header messages */

/* group node (symbol table node) */
static ck_size_t gp_node_sizeof_rkey(global_shared_t *, unsigned);
static ck_err_t gp_node_decode_key(global_shared_t *, unsigned, const uint8_t **, void **_key);
static int gp_node_cmp_key(global_shared_t *shared, unsigned UNUSED ndims, uint8_t *, void *_lt_key, void *_rt_key);

BT_class_t BT_SNODE[1] = {
    BT_SNODE_ID,               	/* id                   */
    sizeof(GP_node_key_t),     	/* sizeof_nkey          */
    gp_node_sizeof_rkey,       	/* get_sizeof_rkey      */
    gp_node_decode_key,        	/* decode               */
    gp_node_cmp_key        	/* compare two keys 	*/
};

/* chunked raw data node */
static size_t raw_node_sizeof_rkey(global_shared_t *, unsigned);
static ck_err_t raw_node_decode_key(global_shared_t *shared, unsigned ndims, const uint8_t **p, void **_key);
static int raw_node_cmp_key(global_shared_t *shared, unsigned ndims, uint8_t *UNUSED heap_chunk, void *_lt_key, void *_rt_key);

static BT_class_t BT_ISTORE[1] = {
    BT_ISTORE_ID,             	/* id                   */
    sizeof(RAW_node_key_t),   	/* sizeof_nkey          */
    raw_node_sizeof_rkey,     	/* get_sizeof_rkey      */
    raw_node_decode_key,      	/* decode               */
    raw_node_cmp_key      	/* compare two keys 	*/
};

static const BT_class_t *const node_key_g[] = {
    BT_SNODE,  	/* group node: symbol table */
    BT_ISTORE	/* raw data chunk node */
};


ck_addr_t
get_logical_addr(const uint8_t *p, const uint8_t *start, ck_addr_t base)
{
    ck_addr_t	diff;

    diff = p - start;
    return(base + diff);
} /* get_logical_addr() */

static int
vector_cmp(unsigned n, const ck_hsize_t *v1, const ck_hsize_t *v2)
{
    int ret_value=0;    /* Return value */

    if (v1 == v2) CK_SET_RET_DONE(0)
    if (v1 == NULL) CK_SET_RET_DONE(-1)
    if (v2 == NULL) CK_SET_RET_DONE(1)
    while (n--) {
        if (*v1 < *v2) CK_SET_RET_DONE(-1)
        if (*v1 > *v2) CK_SET_RET_DONE(1)
        v1++;
        v2++;
    }

done:
    return(ret_value);
} /* vector_cmp() */

/*
 *  For handling hard links
 */

int
table_init(table_t **obj_table)
{
    int 	i;
    table_t	*tb;

    if((tb = malloc(sizeof(table_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Could not malloc() table_t", -1, NULL);
	return(FAIL);
    }
    tb->size = 20;
    tb->nobjs = 0;
    tb->objs = malloc(tb->size*sizeof(obj_t));
    for (i = 0; i < tb->size; i++) {
	tb->objs[i].objno = 0;
	tb->objs[i].nlink = 0;
    }
    *obj_table = tb;
    return(SUCCEED);
} /* table_init() */

static int
table_search(ck_addr_t obj_id)
{
    int	i;

    if (obj_table == NULL)
	return(FAIL);
    for (i = 0; i < obj_table->nobjs; i++) {
	if (obj_table->objs[i].objno == obj_id) /* FOUND */
		return(i);
    }
    return(FAIL);
} /* table_search() */

static void
table_insert(ck_addr_t objno, int nlink)
{
    int	i;

    if (obj_table == NULL)
	return;
    if (obj_table->nobjs == obj_table->size) {
	obj_table->size *= 2;
	obj_table->objs =  realloc(obj_table->objs, obj_table->size * sizeof(obj_t));
    }
    for (i = obj_table->nobjs; i < obj_table->size; i++) {
	obj_table->objs[i].objno = 0;
	obj_table->objs[i].nlink = 0;
    }
    i = obj_table->nobjs++;
    obj_table->objs[i].objno = objno;
    obj_table->objs[i].nlink = nlink;
} /* table_insert() */

/* 
 * Handle symbol table names 
 */
static int
name_list_init(name_list_t **name_list)
{
    name_list_t *nl;

    if((nl = malloc(sizeof(name_list_t))) == NULL) {
	return(-1);
    }
    nl->head = nl->tail = NULL;

    *name_list = nl;
    return(1);
} /* name_list_init() */


static int
name_list_search(name_list_t *nl, char *symname)
{
    name_t *ptr;
    int found = FALSE;

    if(nl->head == NULL && nl->tail == NULL) /* empty */
	return(found);

    ptr = nl->head;
    while(ptr != NULL) {
	if(!strcmp(ptr->name, symname)) {
	    found = TRUE;
	    break;
	}
	ptr = ptr->next;
    }
    return(found);
} /* name_list_search() */

static int
name_list_insert(name_list_t *nl, char *symname)
{
    name_t *name;

    if((name = malloc(sizeof(name_t))) == NULL) {
	return(-1);
    }

    name->name = strdup(symname);
    name->next = NULL;

    if(nl->head == NULL && nl->tail == NULL) /* empty */
	nl->head = nl->tail = name;
    else {
	nl->tail->next = name;
	nl->tail = name;
    }
    return(1);
} /* name_list_insert() */

static int
name_list_dest(name_list_t *nl)
{
    name_t *ptr;

    if(nl->head == NULL && nl->tail == NULL) /* empty */
	return(1);

    while(nl->head != NULL) {
	ptr = nl->head;
	free(nl->head->name);
	nl->head = nl->head->next;
	free(ptr);
    }
    return(1);
} /* name_list_dest() */

/*
 *  Virtual drivers
 */
/* 
 * Set the given driver id based on the given driver name 
 */
static void
set_driver_id(int *driverid/*out*/, char *driver_name) 
{
    if (strcmp(driver_name, "NCSAmult") == 0)
	*driverid = MULTI_DRIVER;
    else if (strcmp(driver_name, "NCSAfami") == 0)
	*driverid = FAMILY_DRIVER;
    else
	/* default driver */
	*driverid = SEC2_DRIVER;
} /* set_driver_id() */

/* 
 * get driver class specific methods based on the given driver id 
 */
static driver_class_t *
get_driver_class(int driver_id)
{
    driver_class_t	*driver=NULL;
    ck_size_t		size;

    size = sizeof(driver_class_t);

    switch (driver_id) {
    case SEC2_DRIVER:
	driver = calloc(1, size);
	memcpy(driver, &sec2_g, size);
	break;
    case MULTI_DRIVER:
	driver = calloc(1, size);
	memcpy(driver, &multi_g, size);
	break;
    case FAMILY_DRIVER:
	driver = calloc(1, size);
	memcpy(driver, &family_g, size);
	break;
    default:
	error_push(ERR_LEV_0, ERR_LEV_0B, "Unsupported file driver", -1, NULL);
	break;
    }  /* end switch */

    return(driver);
} /* get_driver_class() */


/* 
 * Get the driver specific info from shared.
 * It should only be called by multi or family driver.
 */
static void *
get_driver_info(int driver_id, global_shared_t *shared)
{
    void	*driverinfo=NULL;
    ck_size_t	size;

    if (driver_id == MULTI_DRIVER) {
	size = sizeof(driver_multi_fapl_t);
	driverinfo = calloc(1, size);
	memcpy(driverinfo, shared->fa, size);
    } else if (driver_id == FAMILY_DRIVER) {
	size = sizeof(driver_fami_fapl_t);
	driverinfo = calloc(1, size);
	memcpy(driverinfo, shared->fa, size);
    } else
	error_push(ERR_LEV_0, ERR_LEV_0B, "Unsupported file driver", -1, NULL);

    return(driverinfo);
} /* get_driver_info() */


/*
 * Decode the Driver Information field of the Driver Information Block or Driver Info Message
     *	Family driver is encoded in the super block's Driver Information block for
 *	library version 1.8 and after.
 */
static ck_err_t
decode_driver(global_shared_t *shared, const uint8_t *buf)
{
    ck_err_t	ret_value = SUCCEED;

    if (shared->driverid == MULTI_DRIVER)
	ret_value = multi_decode_driver(shared, buf);
    else if ((shared->driverid == FAMILY_DRIVER) && (g_format_num == FORMAT_ONE_EIGHT))
	ret_value = family_decode_driver(shared, buf);
    else {
	error_push(ERR_LEV_0, ERR_LEV_0B, "Unsupported file driver", -1, NULL);
	ret_value = FAIL;
    }

    return(ret_value);
} /* decode_driver() */

driver_t *
FD_open(const char *name, global_shared_t *shared, int driver_id)
{
    driver_t		*file= NULL;
    driver_class_t 	*driver;

    if((driver = get_driver_class(driver_id)) != NULL) {
    	if((file = driver->open(name, shared, driver_id)) != NULL) {
	   file->cls = driver;
	   file->driver_id = driver_id;
	   file->shared = shared;
	}
    }
    return(file);
} /* FD_open() */


ck_err_t
FD_close(driver_t *file)
{
    ck_err_t	ret_value=SUCCEED;
    
    assert(file);
    if(file->cls->close(file) == FAIL)
	ret_value = FAIL;

    return(ret_value);
} /* FD_close() */

ck_err_t
FD_read(driver_t *file, ck_addr_t addr, size_t size, void *buf/*out*/)
{
    int		ret = SUCCEED;
    ck_addr_t	new_addr;


    assert(file);
    assert(file->shared);
    assert(file->cls);

    /* adjust logical "addr" to be the physical address */
    new_addr = addr + file->shared->super_addr;
    if (file->cls->read(file, new_addr, size, buf) == FAIL)
	ret = FAIL;

    return(ret);
} /* FD_read() */

ck_addr_t
FD_get_eof(driver_t *file)
{
   ck_addr_t     ret_value;

   assert(file && file->cls);

   if (file->cls->get_eof)
	ret_value = file->cls->get_eof(file);
    return(ret_value);
} /* FD_get_eof() */

char *
FD_get_fname(driver_t *file, ck_addr_t logi_addr)
{
    ck_addr_t	new_logi_addr;

    assert(file);
    assert(file->cls);

    /* adjust logical "addr" to be the physical address */
    new_logi_addr = logi_addr + file->shared->super_addr;
    return(file->cls->get_fname(file, new_logi_addr));
} /* FD_get_fname() */


/*
 *  sec2 file driver
 */
static driver_t *
sec2_open(const char *name, global_shared_t UNUSED *shared, int UNUSED driver_id)
{
    driver_sec2_t	*file = NULL;
    driver_t	*ret_value;
    int		fd = -1;
    struct stat	sb;

    /* do unix open here */
    fd = open(name, O_RDONLY);
    if (fd < 0) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Unable to open the file", -1, NULL);
	goto done;
    }
    if (fstat(fd, &sb)<0) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Unable to fstat file", -1, NULL);
	goto done;
    }
    file = calloc(1, sizeof(driver_sec2_t));
    file->fd = fd;
    file->eof = sb.st_size;
    file->name = strdup(name);
done:
    ret_value = (driver_t *)file;
    return(ret_value);
} /* sec2_open() */

static ck_addr_t
sec2_get_eof(driver_t *file)
{
    driver_sec2_t *sec2_file = (driver_sec2_t*)file;
 
    return(sec2_file->eof);

} /* sec2_get_eof() */

static char *
sec2_get_fname(driver_t *file, ck_addr_t UNUSED logi_addr)
{
    return(((driver_sec2_t *)file)->name);
} /* sec2_get_fname() */

static ck_err_t
sec2_close(driver_t *_file)
{
    int	ret;
    driver_sec2_t   	*file=(driver_sec2_t*)_file;
 
    ret = close(file->fd);
    if (file->name) free(file->name);
    free(file);
    return(ret);
} /* sec2_close() */


/* NEED To find out about lseek64??? */
static ck_err_t
sec2_read(driver_t *file, ck_addr_t addr, size_t size, void *buf/*out*/)
{
    int	status;
    int	ret = SUCCEED;
    int	fd;
size_t bytes;

    
    if (addr == CK_ADDR_UNDEF)
	ret = FAIL;
    else if ((addr + size) > ((driver_sec2_t *)file)->eof)
	ret = FAIL;
    else {
	fd = ((driver_sec2_t *)file)->fd;
	if (lseek(fd, addr, SEEK_SET) < 0)
	    ret = FAIL;
	else if ((bytes=read(fd, buf, size)) < 0)
	    ret = FAIL;
    }

#ifdef TOBEFIXED
check for number of bytes returned from read??
#endif
    return(ret);
} /* sec2_read() */


/*
 *  multi file driver
 *  1. User needs to give the 1st file "multi-s.h5" so that h5check will get 
 *     superblock information for validation of all the multi files.
 *  2. (driver_multi_t *)file->name is the original file name supplied by the user
 *  3. (driver_multi_t *)file->fa.memb_name is the sprintf generated name 
 *     e.g. %s-s.h5, %s-b.h5, %s-r.h5 ....
 *  4. Problem: need to solve H5Tconvert(H5T_STD_U64LE, H5T_NATIVE_HADDR, nseen*2, x, NULL, H5P_DEFAULT)
 */

#define UNIQUE_MEMBERS(MAP,LOOPVAR) {                                         \
    driver_mem_t _unmapped, LOOPVAR;                                            \
    ck_bool_t _seen[FD_MEM_NTYPES];                                           \
                                                                              \
    memset(_seen, 0, sizeof _seen);                                           \
    for (_unmapped=FD_MEM_SUPER; _unmapped<FD_MEM_NTYPES; _unmapped=(driver_mem_t)(_unmapped+1)) {  \
        LOOPVAR = MAP[_unmapped];                                             \
        if (FD_MEM_DEFAULT==LOOPVAR) LOOPVAR=_unmapped;                     \
        assert(LOOPVAR>0 && LOOPVAR<FD_MEM_NTYPES);                         \
        if (_seen[LOOPVAR]++) continue;


#define ALL_MEMBERS(LOOPVAR) {                                                \
    driver_mem_t LOOPVAR;                                                       \
    for (LOOPVAR=FD_MEM_DEFAULT; LOOPVAR<FD_MEM_NTYPES; LOOPVAR=(driver_mem_t)(LOOPVAR+1)) {


#define END_MEMBERS     }}


static void
set_multi_driver_properties(driver_multi_fapl_t **fa, driver_mem_t map[], const char *memb_name[], ck_addr_t memb_addr[])
{
    *fa = calloc(1, sizeof(driver_multi_fapl_t));

    ALL_MEMBERS(mt) {
	(*fa)->memb_map[mt] = map[mt];
   	(*fa)->memb_addr[mt] = memb_addr[mt];
	if (memb_name[mt]) {
	    (*fa)->memb_name[mt] = malloc(strlen(memb_name[mt])+1);
	    strcpy((*fa)->memb_name[mt], memb_name[mt]);
	} else 
	    (*fa)->memb_name[mt] = NULL;
    } END_MEMBERS;
} /* set_multi_driver_properties() */

static ck_err_t
multi_decode_driver(global_shared_t *shared, const unsigned char *buf)
{
    char                x[2*FD_MEM_NTYPES*8];
    driver_mem_t        map[FD_MEM_NTYPES];
    int                 i;
    const char          *memb_name[FD_MEM_NTYPES];
    ck_addr_t           memb_addr[FD_MEM_NTYPES];
    ck_addr_t           memb_eoa[FD_MEM_NTYPES];
    ck_addr_t	    	xx;
    ck_err_t		ret_value=SUCCEED;

    /* Set default values */
    ALL_MEMBERS(mt) {
       	memb_addr[mt] = CK_ADDR_UNDEF;
       	memb_eoa[mt] = CK_ADDR_UNDEF;
       	memb_name[mt] = NULL;
    } END_MEMBERS;

    memset(map, 0, sizeof map);
    for (i=0; i<6; i++) {
       	map[i+1] = (driver_mem_t)buf[i];
    }

    buf += 8;

    /* Decode Address and EOA values */
    /* NOT SURE: the library does a H5Tconvert(...H5T_STD_U64LE, H5T_NATIVE_HADDR...)?? */
    assert(sizeof(ck_addr_t)<=8);
    UNIQUE_MEMBERS(map, mt) {
       	UINT64DECODE(buf, xx);
       	memb_addr[_unmapped] = xx;
       	UINT64DECODE(buf, xx);
       	memb_eoa[_unmapped] = xx;
    } END_MEMBERS;

    /* Decode name templates */
    UNIQUE_MEMBERS(map, mt) {
       	ck_size_t n = strlen((const char *)buf)+1;
       	memb_name[_unmapped] = (const char *)buf;
       	buf += (n+7) & ~((unsigned)0x0007);
    } END_MEMBERS;

    set_multi_driver_properties((driver_multi_fapl_t **)&(shared->fa), map, memb_name, memb_addr);
    return(ret_value);
} /* multi_decode_driver() */

static driver_t *
multi_open(const char *name, global_shared_t *shared, int driver_id)
{
    driver_multi_t      *file=NULL;
    driver_multi_fapl_t *fa;
    driver_mem_t        m;
    driver_t		*ret_ptr=NULL;
    int			ret_value=SUCCEED;
    int			len;

   if (!name || !*name) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Invalid file name", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    if ((file=calloc(1, sizeof(driver_multi_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    fa = get_driver_info(driver_id, shared);
    if (fa == NULL) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Unable to get driver information", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    ALL_MEMBERS(mt) {
	file->fa.memb_map[mt] = fa->memb_map[mt];
	file->fa.memb_addr[mt] = fa->memb_addr[mt];
	if (fa->memb_name[mt]) {
	    file->fa.memb_name[mt] = malloc(strlen(fa->memb_name[mt]));
	    strcpy(file->fa.memb_name[mt], fa->memb_name[mt]);
	} else
	    file->fa.memb_name[mt] = NULL;
    } END_MEMBERS;

    /* name of the file supplied by the user */
    file->name = malloc(strlen(name)+1);
    strcpy(file->name, name);

    if (compute_next(file) != SUCCEED) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Unable to compute member addresses", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    if (open_members(file, shared) != SUCCEED) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Unable to open member files", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    if(ret_value == SUCCEED)
	ret_ptr = (driver_t *)file;
done:
    if(ret_value == FAIL && file) {
	ALL_MEMBERS(mt) {
	    if(file->memb[mt]) 
		(void)FD_close(file->memb[mt]);
	    if(file->fa.memb_name[mt])
		free(file->fa.memb_name[mt]);
	} END_MEMBERS;
	if(file->name) free(file->name);
	free(file);
    }
    return(ret_ptr);
} /* multi_open() */

static ck_err_t 
multi_close(driver_t *_file)
{
    driver_multi_t   	*file=(driver_multi_t*)_file;
    int			errs=0;
    ck_err_t		ret_value=SUCCEED;

    /* Close as many members as possible */
    ALL_MEMBERS(mt) {
	if(file->memb[mt]) {
	    if (FD_close(file->memb[mt])<0)
		errs++;
	    else
		file->memb[mt] = NULL;
	}
    } END_MEMBERS;

    if(errs) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Error closing member file(s)", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

done:
    ALL_MEMBERS(mt) {
	if(file->fa.memb_name[mt]) 
	    free(file->fa.memb_name[mt]);
    } END_MEMBERS;

    if(file->name) free(file->name);
    free(file);

    return(ret_value);
} /* multi_close() */

static ck_err_t
compute_next(driver_multi_t *file)
{
    ALL_MEMBERS(mt) {
        file->memb_next[mt] = CK_ADDR_UNDEF;
    } END_MEMBERS;

    UNIQUE_MEMBERS(file->fa.memb_map, mt1) {
        UNIQUE_MEMBERS(file->fa.memb_map, mt2) {
            if (file->fa.memb_addr[mt1]<file->fa.memb_addr[mt2] &&
		(CK_ADDR_UNDEF==file->memb_next[mt1] ||
		file->memb_next[mt1]>file->fa.memb_addr[mt2])) {

                file->memb_next[mt1] = file->fa.memb_addr[mt2];
            }
        } END_MEMBERS;

        if (CK_ADDR_UNDEF==file->memb_next[mt1]) {
            file->memb_next[mt1] = CK_ADDR_MAX; /*last member*/
        }

    } END_MEMBERS;

    return(SUCCEED);
} /* compute_next() */


static ck_err_t 
open_members(driver_multi_t *file, global_shared_t *shared)
{
    char    	tmp[1024], newname[1024];
    char	*ptr;
    int		ret_value=SUCCEED;

    /* fix the name */
    strcpy(newname, file->name);
    ptr = strrchr(newname, '-');
    *ptr = '\0';

    UNIQUE_MEMBERS(file->fa.memb_map, mt) {
	assert(file->fa.memb_name[mt]);
	sprintf(tmp, file->fa.memb_name[mt], newname);
        file->memb[mt] = FD_open(tmp, shared, SEC2_DRIVER);
	if (file->memb[mt] == NULL)
	    ret_value = FAIL;
    } END_MEMBERS;

    return(ret_value);
} /* open_members() */


static ck_err_t
multi_read(driver_t *file, ck_addr_t addr, size_t size, void *buf/*out*/)
{
    driver_multi_t      *ff = (driver_multi_t*)file;
    driver_mem_t        mt, mmt, hi = FD_MEM_DEFAULT;
    ck_addr_t           start_addr = 0;
    ck_err_t		ret_value = SUCCEED;

    /* Find the file to which this address belongs */
    for(mt = FD_MEM_SUPER; mt < FD_MEM_NTYPES; mt = (driver_mem_t)(mt + 1)) {
	mmt = ff->fa.memb_map[mt];
	if(FD_MEM_DEFAULT == mmt) 
	    mmt = mt;
	if((mmt <= 0) || (mmt >= FD_MEM_NTYPES)) {
	    error_push(ERR_FILE, ERR_NONE_SEC, "Invalid member mapping type", -1, NULL);
	    CK_SET_RET_DONE(FAIL)
	}

        if(ff->fa.memb_addr[mmt]>addr) 
	    continue;
        if(ff->fa.memb_addr[mmt]>=start_addr) {
	    start_addr = ff->fa.memb_addr[mmt];
	    hi = mmt;
        }
    }
    assert(hi > 0);

    /* Read from that member */
    if(FD_read(ff->memb[hi], addr-start_addr, size, buf) == FAIL) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Error reading member file", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

done:
    return(ret_value);
} /* multi_read() */


static ck_addr_t
multi_get_eof(driver_t *file)
{
    driver_multi_t    *multi_file = (driver_multi_t *)file;
    ck_addr_t         eof=0, tmp;

    UNIQUE_MEMBERS(multi_file->fa.memb_map, mt) {
	if(multi_file->memb[mt]) {
	    tmp = FD_get_eof(multi_file->memb[mt]);
	    if(tmp == CK_ADDR_UNDEF) {
		error_push(ERR_FILE, ERR_NONE_SEC, "Member file has unknown eof", -1, NULL);
		return(CK_ADDR_UNDEF);
	    }
	    if(tmp > 0) 
		tmp += multi_file->fa.memb_addr[mt];
        } else {
	    error_push(ERR_FILE, ERR_NONE_SEC, "Bad eof", -1, NULL);
	    return(CK_ADDR_UNDEF);
        }

        if(tmp > eof) 
	    eof = tmp;
    } END_MEMBERS;

    return(eof);
} /* multi_get_eof() */

static char *
multi_get_fname(driver_t *_file, ck_addr_t logi_addr)
{
    char 		*tmp, *newname, *ptr;
    driver_multi_t      *file = (driver_multi_t*)_file;
    driver_mem_t        mt, mmt, hi=FD_MEM_DEFAULT;
    ck_addr_t           start_addr=0;

    /* Find the file to which this address belongs */
    for(mt = FD_MEM_SUPER; mt < FD_MEM_NTYPES; mt = (driver_mem_t)(mt+1)) {
	mmt = file->fa.memb_map[mt];
	if(FD_MEM_DEFAULT == mmt) 
	    mmt = mt;
        assert(mmt>0 && mmt<FD_MEM_NTYPES);

        if(file->fa.memb_addr[mmt] > logi_addr) 
	    continue;
        if(file->fa.memb_addr[mmt] >= start_addr) {
	    start_addr = file->fa.memb_addr[mmt];
	    hi = mmt;
	}
    }
    assert(hi>0);

    /* tmp is the 1st part of the name before "-" */
    tmp = strdup(file->name);
    ptr = strrchr(tmp, '-');
    *ptr = '\0';
    newname = malloc(strlen(tmp) + 5); 
    sprintf(newname, file->fa.memb_name[hi], tmp);
    return(newname);
} /* multi_get_fname() */

/*
 *  Family driver
 *
 *  1. User needs to give the 1st file "family00000.h5" so that h5check will get 
 *     superblock information for validation of all the family files.
 *  2. (driver_fami_t *)file->name is the original file name supplied by the user
 *  3. To generate family file names, the format follows the convention:
 *	family%05d.h5
 *  4. Problem: IA32 Architecture whose size_t's size is 4 bytes....for SIZET_MAX in family_read()
 * 	      : family_open() at the end of closing files when nerrors...goto done again??
 */
static void
set_family_driver_properties(driver_fami_fapl_t **fa, ck_hsize_t msize)
{
    *fa = calloc(1, sizeof(driver_fami_fapl_t));

    (*fa)->memb_size = msize;
} /* set_family_driver_properties() */

static ck_err_t
family_decode_driver(global_shared_t *shared, const unsigned char *buf)
{
    uint64_t            msize;
    ck_err_t		ret_value=SUCCEED;

    UINT64DECODE(buf, msize);
    set_family_driver_properties((driver_fami_fapl_t **)&(shared->fa), msize);	    
    return(ret_value);
} /* family_decode_driver() */

static driver_t *
family_open(const char *name, global_shared_t *shared, int driver_id)
{
    driver_fami_t      *file=NULL;
    driver_fami_fapl_t *fa;
    driver_mem_t        m;
    ck_hsize_t        	eof=CK_ADDR_UNDEF;
    driver_t		*ret_ptr=NULL;
    ck_err_t		ret_value=SUCCEED;
    char                memb_name[4096], temp[4096];
    int			len, ret;
    char		*ptr;

    if(!name || !*name) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Invalid file name", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    if((file=calloc(1, sizeof(driver_fami_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    fa = get_driver_info(driver_id, shared);
    if(fa == NULL) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Unable to get driver information", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    file->fa.memb_size = fa->memb_size;

    /* name of the file supplied by the user */
    file->name = malloc(strlen(name)+1);
    strcpy(file->name, name);

    /* fix up the name to be "name%05d.h5" */
    /* convention:family file's name contains 5 digits after "name" */
    strcpy(temp, name);
    ptr = strrchr(temp, '.');
    ptr = ptr - 5;
    *ptr = '\0';
    strcat(temp, "%05d.h5");

    while (1) {
	sprintf(memb_name, temp, file->nmembs);
        if(file->nmembs >= file->amembs) {
            unsigned n 	= MAX(64, 2*file->amembs);
            driver_t 	**x = realloc(file->memb, n*sizeof(driver_t *));

            if(!x) {
		error_push(ERR_INTERNAL, ERR_NONE_SEC, "Internal allocation error", -1, NULL);
		CK_SET_RET_DONE(FAIL)
	    }
            file->amembs = n;
            file->memb = x;
	}
	file->memb[file->nmembs] = FD_open(memb_name, shared, SEC2_DRIVER);
	if(!file->memb[file->nmembs]) {
	    /* if failure in opening the very first file, error */
	    /* otherwise, done with opening all member files */
	    if(file->nmembs == 0) {
		error_push(ERR_FILE, ERR_NONE_SEC, "Unable to open member file", -1, NULL);
		CK_SET_RET_DONE(FAIL)
	    }
	    /* clear the errors in opening the last file */
	    error_clear();
	    break;
	}
	file->nmembs++;
    }

    /* update member file size in case there is only one file */
    if(eof = FD_get_eof(file->memb[0])) 
	file->fa.memb_size = eof;

    if(ret_value == SUCCEED)
	ret_ptr = (driver_t *)file;
done:
    if(ret_value == FAIL && file) {
	unsigned nerrors = 0;
	unsigned u;     

	for(u = 0; u < file->nmembs; u++)
	    if(file->memb[u])
		if(FD_close(file->memb[u]) < 0)
		    nerrors++;

	/* NEED check on this */
	if(nerrors)
	    error_push(ERR_FILE, ERR_NONE_SEC, "Unable to close member file(s)", -1, NULL);
	if(file->memb)
	    free(file->memb);
	if(file->name)
	    free(file->name);
	free(file);
    }
    return(ret_ptr);
} /* family_open() */

static ck_err_t 
family_close(driver_t *_file)
{
    driver_fami_t   	*file=(driver_fami_t*)_file;
    unsigned            nerrors=0;   
    unsigned		u;
    ck_err_t		ret_value=SUCCEED;

    for(u = 0; u < file->nmembs; u++) {
        if (file->memb[u]) {
            if(FD_close(file->memb[u]) < 0)
                nerrors++;
            else
                file->memb[u] = NULL;
        }
    }
    /* need check on this */
    if(nerrors)
	CK_SET_RET_DONE(FAIL)
done:
    if (file->memb) free(file->memb);
    if (file->name) free(file->name);
    return(ret_value);
} /* family_close() */

static ck_err_t
family_read(driver_t *_file, ck_addr_t addr, size_t size, void *_buf/*out*/)
{
    driver_fami_t       *file=(driver_fami_t*)_file;
    ck_err_t            ret_value=SUCCEED; 
    ck_addr_t          	sub;
    ck_size_t   	req;
    ck_size_t          	tempreq;
    unsigned            u; 
    char *		buf=(char*)_buf;

    while(size > 0) {
	ASSIGN_OVERFLOW(u,addr/file->fa.memb_size,ck_hsize_t,unsigned);
        sub = addr % file->fa.memb_size;

#ifdef OUT
        /* This check is for mainly for IA32 architecture whose size_t's size
         * is 4 bytes, to prevent overflow when user application is trying to
         * write files bigger than 4GB. */
        tempreq = file->memb_size-sub;
        if(tempreq > SIZET_MAX)
            tempreq = SIZET_MAX;
        req = MIN(size, (size_t)tempreq);
#endif
        tempreq = file->fa.memb_size - sub;
        req = MIN(size, (ck_size_t)tempreq);

	if(u >= file->nmembs) {
	    error_push(ERR_FILE, ERR_NONE_SEC, "Error reading member file", -1, NULL);
	    CK_SET_RET_DONE(FAIL)
	}

        if(FD_read(file->memb[u], sub, req, buf) == FAIL) {
	    error_push(ERR_FILE, ERR_NONE_SEC, "Error reading member file", -1, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
        addr += req;
        buf += req;
        size -= req;
    }

done:
    return(ret_value);
} /* family_read() */

static ck_addr_t
family_get_eof(driver_t *_file)
{
    driver_fami_t    	*file = (driver_fami_t *)_file;
    ck_addr_t        	eof=0, tmp;
    ck_addr_t 		ret_value; 
    int			i;

    assert(file->nmembs>0);
    /*
     * Find the last member that has a non-zero EOF and break out of the loop
     * with `i' equal to that member. If all members have zero EOF then exit
     * loop with i==0.
     */
    for(i = (int)file->nmembs-1; i >= 0; --i) {
        if((eof = FD_get_eof(file->memb[i])) != 0)
            break;
        if (0 == i)
            break;
    }

    /*
     * The file size is the number of members before the i'th member plus the
     * size of the i'th member.
     */
    eof += ((unsigned)i)*file->fa.memb_size;
    ret_value = MAX(eof, file->eoa);

done:
    return(ret_value);
} /* family_get_eof() */

static char *
family_get_fname(driver_t *_file, ck_addr_t logi_addr)
{
    char *temp, *newname;
    char *ptr;
    unsigned u; 
    driver_fami_t *file = (driver_fami_t *)_file;

    /* u is the family file member # */
    ASSIGN_OVERFLOW(u, logi_addr/file->fa.memb_size, ck_hsize_t, unsigned);

    /* fix up the name to be "name%05d.h5" */
    /* convention:family file's name contains 5 digits after "name" */
    temp = strdup(file->name);
    ptr = strrchr(temp, '.');
    ptr = ptr - 5;
    *ptr = '\0';
    strcat(temp, "%05d.h5");
    newname = malloc(strlen(temp)+10);
    sprintf(newname, temp, u);

    return(newname);
} /* family_get_fname() */
/*
 * End virtual file drivers
 */


/* 
 * Decoder for various messages
 */

/* Dataspace */
static void *
OBJ_sds_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_sds_extent_t  	*sdim=NULL;
    unsigned        	i;              
    unsigned        	flags, version;
    ck_addr_t		logical;
    int			ret_value=SUCCEED;
    void          	*ret_ptr=NULL;
    int			badinfo;

    assert(file);
    assert(p);

    logical = get_logical_addr(p, start_buf, logi_base);

    sdim = calloc(1, sizeof(OBJ_sds_extent_t));
    if(sdim == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Dataspace Message:Internal allocation error", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    version = *p++;

    /* version adjusted if error for further continuation */
    if(g_format_num == FORMAT_ONE_SIX) {
	if(version != OBJ_SDS_VERSION_1) {
	    badinfo = version;
	    version = OBJ_SDS_VERSION_1;
	    error_push(ERR_LEV_2, ERR_LEV_2A2b, "Dataspace Message v.1:Wrong version number", logical, &badinfo);
	    CK_SET_RET(FAIL)
	}
    } else {
	if(version < OBJ_SDS_VERSION_1 || version > OBJ_SDS_VERSION_2) {
	   badinfo = version;
	   version = OBJ_SDS_VERSION_2;
	   error_push(ERR_LEV_2, ERR_LEV_2A2b, "Dataspace Message:Wrong version number", logical, &badinfo);
	   CK_SET_RET(FAIL)
	}
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    sdim->rank = *p++;
    if(sdim->rank > OBJ_SDS_MAX_RANK) {
       badinfo = sdim->rank;
       error_push(ERR_LEV_2, ERR_LEV_2A2b, "Dataspace Message:Dimensionality is too large", logical, &badinfo);
       CK_SET_RET(FAIL)
    }

    /* Get dataspace flags for later */
    logical = get_logical_addr(p, start_buf, logi_base);
    flags = *p++;

    if((version == OBJ_SDS_VERSION_1) && (flags > 0x3)) {  
	/* For version 1: Only bit 0 and bit 1 can be set */
	error_push(ERR_LEV_2, ERR_LEV_2A2b, "Dataspace Message v.1:Corrupt flags", logical, NULL);
	CK_SET_RET(FAIL)
    } else if((version == OBJ_SDS_VERSION_2) && (flags > 0x1)) { 
	/* For version 2: Only bit 0 can be set */
	error_push(ERR_LEV_2, ERR_LEV_2A2b, "Dataspace Message v.2:Corrupt flags", logical, NULL);
	CK_SET_RET(FAIL)
    }

    if(version >= OBJ_SDS_VERSION_2) {
	logical = get_logical_addr(p, start_buf, logi_base);
        sdim->type = (OBJ_sds_class_t)*p++;
	if((sdim->type != OBJ_SDS_SCALAR) && 
	   (sdim->type != OBJ_SDS_SIMPLE) && 
	   (sdim->type != OBJ_SDS_NULL)) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2b, "Dataspace Message v.2:Invalid type", logical, NULL);
	    CK_SET_RET(FAIL)
	}
    } else {
	if(sdim->rank > 0)
	    sdim->type = OBJ_SDS_SIMPLE;
	else
	    sdim->type = OBJ_SDS_SCALAR;
	p++;
    }

    if(version == OBJ_SDS_VERSION_1)
        p += 4; /*reserved*/

    logical = get_logical_addr(p, start_buf, logi_base);
    if(sdim->rank > 0) {
	if((sdim->size = (ck_hsize_t *)malloc(sizeof(ck_hsize_t) * (ck_size_t)sdim->rank))==NULL) {
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, "Dataspace Message:Internal allocation error", 
		logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
        for(i = 0; i < sdim->rank; i++)
	    DECODE_LENGTH(file->shared, p, sdim->size[i]);

	logical = get_logical_addr(p, start_buf, logi_base);
        if(flags & OBJ_SDS_VALID_MAX) {
             if((sdim->max = (ck_hsize_t *)malloc(sizeof(ck_hsize_t) * (ck_size_t)sdim->rank))==NULL) {
		error_push(ERR_INTERNAL, ERR_NONE_SEC, "Dataspace Message:Internal allocation error", 
		    logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }
	    for(i = 0; i < sdim->rank; i++)
            	DECODE_LENGTH(file->shared, p, sdim->max[i]);
        }
    }

    /* Compute the number of elements in the extent */
    if(sdim->type == OBJ_SDS_NULL)
        sdim->nelem = 0;
    else {
	for (i = 0, sdim->nelem = 1; i < sdim->rank; i++)
	    sdim->nelem *= sdim->size[i];
    }

    if(ret_value == SUCCEED)
	ret_ptr = (void*)sdim;

done:
    if(ret_value == FAIL && sdim) {
	if (sdim->size) free(sdim->size);
	if (sdim->max) free(sdim->max);
	free(sdim);
    }

    return(ret_ptr);
} /* OBJ_sds_decode() */

/* Link Info */
static void *
OBJ_linfo_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_linfo_t 	*mesg=NULL;  
    unsigned char 	index_flags;  		/* Flags for encoding link index info */
    int			badinfo;
    int			ret_value=SUCCEED;
    void        	*ret_ptr=NULL;
    unsigned        	version;
    ck_addr_t		logical;

    assert(file);
    assert(p);

    if(g_format_num == FORMAT_ONE_SIX) {
       error_push(ERR_LEV_2, ERR_LEV_2A2c, "Link Info Message:Unsupported message", -1, NULL);
       CK_SET_RET_DONE(FAIL)
    }

    /* Version of message */
    logical = get_logical_addr(p, start_buf, logi_base);
    version = *p++;
    if(version != OBJ_LINFO_VERSION) {
	badinfo = version;
	error_push(ERR_LEV_2, ERR_LEV_2A2c, 
	    "Link Info Message:Bad version number", logical, &badinfo);
	CK_SET_RET(FAIL)
    }

    /* Allocate space for message */
    if((mesg = calloc(1, sizeof(OBJ_linfo_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Link Info Message:Internal allocation error", logical, NULL);
        CK_SET_RET_DONE(FAIL)
    }

    /* Get the index flags for the group */
    index_flags = *p++;
    if(index_flags & ~OBJ_LINFO_ALL_FLAGS) {
	error_push(ERR_LEV_2, ERR_LEV_2A2c, "Link Info Message:Bad flag value", logical, NULL);
	CK_SET_RET(FAIL)
    }
    mesg->track_corder = (index_flags & OBJ_LINFO_TRACK_CORDER) ? TRUE : FALSE;
    mesg->index_corder = (index_flags & OBJ_LINFO_INDEX_CORDER) ? TRUE : FALSE;

    if(mesg->track_corder)
        UINT64DECODE(p, mesg->max_corder)
    else
        mesg->max_corder = 0;

    /* Address of fractal heap to store "dense" links */
    addr_decode(file->shared, &p, &(mesg->fheap_addr));

    /* Address of v2 B-tree to index names of links (names are always indexed) */
    addr_decode(file->shared, &p, &(mesg->name_bt2_addr));

    /* Address of v2 B-tree to index creation order of links, if there is one */
    if(mesg->index_corder)
        addr_decode(file->shared, &p, &(mesg->corder_bt2_addr));
    else
        mesg->corder_bt2_addr = CK_ADDR_UNDEF;

    if (ret_value == SUCCEED)
	ret_ptr = mesg;
done:
    if(ret_value == FAIL && mesg)
	free(mesg);

    return(ret_ptr);
} /* OBJ_linfo_decode() */

/* Datatype */
static ck_err_t
OBJ_dt_decode_helper(driver_t *file, const uint8_t **pp, type_t *dt, const uint8_t *start_buf, ck_addr_t logi_base)
{
    unsigned        	flags, version, tmp;
    unsigned        	i, j;
    size_t        	z;
    ck_err_t    	ret_value=SUCCEED;
    ck_addr_t		logical;
    int			badinfo;

    /* check args */
    assert(file);
    assert(pp && *pp);
    assert(dt && dt->shared);

    logical = get_logical_addr(*pp, start_buf, logi_base);
    UINT32DECODE(*pp, flags);
    version = (flags>>4) & 0x0f;

    /* version is adjusted if error for continuation */
    if(g_format_num == FORMAT_ONE_SIX) {
	if(version != DT_VERSION_1 && version != DT_VERSION_2) {
	   badinfo = version;
	   version = DT_VERSION_2;
	   error_push(ERR_LEV_2, ERR_LEV_2A2d, "Datatype Message:Bad version number", logical, &badinfo);
	   CK_SET_RET(FAIL)
	}
    } else {
	if(version < DT_VERSION_1 || version > DT_VERSION_3) {
	   badinfo = version;
	   version = DT_VERSION_3;
	   error_push(ERR_LEV_2, ERR_LEV_2A2d, "Datatype Message:Bad version number", logical, &badinfo);
	   CK_SET_RET(FAIL)
	}
    }

    dt->shared->type = (DT_class_t)(flags & 0x0f);
    if ((dt->shared->type < DT_INTEGER) ||(dt->shared->type >DT_ARRAY)) {
	error_push(ERR_LEV_2, ERR_LEV_2A2d, "Datatype Message:Invalid class value", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    flags >>= 8;
    UINT32DECODE(*pp, dt->shared->size);

    switch(dt->shared->type) {
	case DT_INTEGER:  /* Fixed-Point datatypes */
            dt->shared->u.atomic.order = (flags & 0x1) ? DT_ORDER_BE : DT_ORDER_LE;
            dt->shared->u.atomic.lsb_pad = (flags & 0x2) ? DT_PAD_ONE : DT_PAD_ZERO;
            dt->shared->u.atomic.msb_pad = (flags & 0x4) ? DT_PAD_ONE : DT_PAD_ZERO;
            dt->shared->u.atomic.u.i.sign = (flags & 0x8) ? DT_SGN_2 : DT_SGN_NONE;

	    if((flags>>4) != 0) {  /* should be reserved(zero) */
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Fixed-Point:Bits 4-23 should be 0 for class bit field", 
		  logical, NULL);
		CK_SET_RET(FAIL)
	    }

            UINT16DECODE(*pp, dt->shared->u.atomic.offset);
            UINT16DECODE(*pp, dt->shared->u.atomic.prec);
            break;

	case DT_FLOAT:  /* Floating-Point datatypes */
	    dt->shared->u.atomic.order = (flags & 0x1) ? DT_ORDER_BE : DT_ORDER_LE;
	    if(version == DT_VERSION_1 || version == DT_VERSION_2) {
                if (flags & 0x40) { /* bit 6 should be reserved */
		    error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		       	"Datatype Message:Floating-Point:Bit 6 should be reserved", logical, NULL);
		    CK_SET_RET(FAIL)
		}
	    } else {
		if(version == DT_VERSION_3) {
		    if((flags & 0x40) && (flags & 0x1))
			dt->shared->u.atomic.order = DT_ORDER_VAX;
		    else {
			error_push(ERR_LEV_2, ERR_LEV_2A2d, 
			    "Datatype Message:Floating-Point:Bad byte order for VAX-endian", logical, NULL);
			CK_SET_RET(FAIL)
		    }
		} 
	    }
            dt->shared->u.atomic.lsb_pad = (flags & 0x2) ? DT_PAD_ONE : DT_PAD_ZERO;
            dt->shared->u.atomic.msb_pad = (flags & 0x4) ? DT_PAD_ONE : DT_PAD_ZERO;
            dt->shared->u.atomic.u.f.pad = (flags & 0x8) ? DT_PAD_ONE : DT_PAD_ZERO;
            switch((flags >> 4) & 0x03) {
                case 0:
                    dt->shared->u.atomic.u.f.norm = DT_NORM_NONE;
                    break;

                case 1:
                    dt->shared->u.atomic.u.f.norm = DT_NORM_MSBSET;
                    break;

                case 2:
                    dt->shared->u.atomic.u.f.norm = DT_NORM_IMPLIED;
                    break;

                default:
		    error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		      "Datatype Message:Unknown Floating-Point normalization", logical, NULL);
		    CK_SET_RET(FAIL)
		    break;
            } /* end switch */

	    if(((flags >> 7) & 0x01) != 0) {
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Floating-Point:Bit 7 should be 0 for class bit field", 
		  logical, NULL);
		CK_SET_RET(FAIL)
	    }
	    if((flags >> 16) != 0) {
		error_push(ERR_LEV_2, ERR_LEV_2A2d,
		  "Datatype Message:Floating-Point:Bits 16-23 should be 0 for class bit field", 
		  logical, NULL);
		CK_SET_RET(FAIL)
	    }

            dt->shared->u.atomic.u.f.sign = (flags >> 8) & 0xff;
            UINT16DECODE(*pp, dt->shared->u.atomic.offset);
            UINT16DECODE(*pp, dt->shared->u.atomic.prec);
            dt->shared->u.atomic.u.f.epos = *(*pp)++;
            dt->shared->u.atomic.u.f.esize = *(*pp)++;

            if(dt->shared->u.atomic.u.f.esize <= 0) {
	       error_push(ERR_LEV_2, ERR_LEV_2A2d,
		  "Datatype Message:Floating-Point:size of exponent should be greater than 0", 
		  logical, NULL);
	       CK_SET_RET(FAIL)
	    }

            dt->shared->u.atomic.u.f.mpos = *(*pp)++;
            dt->shared->u.atomic.u.f.msize = *(*pp)++;

            if(dt->shared->u.atomic.u.f.msize <= 0) {
	       error_push(ERR_LEV_2, ERR_LEV_2A2d,
		  "Datatype Message:size of matissa should be greater than 0", logical, NULL);
	       CK_SET_RET(FAIL)
	    }

            UINT32DECODE(*pp, dt->shared->u.atomic.u.f.ebias);
            break;

	case DT_TIME:  /* Time datatypes */
            dt->shared->u.atomic.order = (flags & 0x1) ? DT_ORDER_BE : DT_ORDER_LE;

	    if((flags>>1) != 0) {  /* should be reserved(zero) */
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Time:Bits 1-23 should be 0 for class bit field", logical, NULL);
		CK_SET_RET(FAIL)
	    }

            UINT16DECODE(*pp, dt->shared->u.atomic.prec);
            break;

	case DT_STRING:  /* character string datatypes */
            dt->shared->u.atomic.order = DT_ORDER_NONE;
            dt->shared->u.atomic.prec = 8 * dt->shared->size;
            dt->shared->u.atomic.offset = 0;
            dt->shared->u.atomic.lsb_pad = DT_PAD_ZERO;
            dt->shared->u.atomic.msb_pad = DT_PAD_ZERO;

            dt->shared->u.atomic.u.s.pad = (DT_str_t)(flags & 0x0f);
            dt->shared->u.atomic.u.s.cset = (DT_cset_t)((flags>>4) & 0x0f);            

	    /* Padding Type supported: Null Terminate, Null Pad and Space Pad */
	    if((dt->shared->u.atomic.u.s.pad != DT_STR_NULLTERM) &&
	       (dt->shared->u.atomic.u.s.pad != DT_STR_NULLPAD) &&
	       (dt->shared->u.atomic.u.s.pad != DT_STR_SPACEPAD)) {
	        error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:String:Unsupported padding type for class bit field", logical, NULL);
	        CK_SET_RET(FAIL)
	    }

	    if((dt->shared->u.atomic.u.s.cset != DT_CSET_ASCII) && 
		(dt->shared->u.atomic.u.s.cset != DT_CSET_UTF8)) {
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:String:Unsupported character set for class bit field", logical, NULL);
		CK_SET_RET(FAIL)
	    }
	    if((flags>>8) != 0) {  /* should be reserved(zero) */
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:String:Bits 8-23 should be 0 for class bit field", logical, NULL);
		CK_SET_RET(FAIL)
	    }
	    break;
	
	case DT_BITFIELD:  /* Bit fields datatypes */
            dt->shared->u.atomic.order = (flags & 0x1) ? DT_ORDER_BE : DT_ORDER_LE;
            dt->shared->u.atomic.lsb_pad = (flags & 0x2) ? DT_PAD_ONE : DT_PAD_ZERO;
            dt->shared->u.atomic.msb_pad = (flags & 0x4) ? DT_PAD_ONE : DT_PAD_ZERO;
	    if((flags>>3) != 0) {  /* should be reserved(zero) */
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Bitfield:Bits 3-23 should be 0 for class bit field", 
		  logical, NULL);
		CK_SET_RET(FAIL)
	    }
            UINT16DECODE(*pp, dt->shared->u.atomic.offset);
            UINT16DECODE(*pp, dt->shared->u.atomic.prec);
            break;

	case DT_OPAQUE:  /* Opaque datatypes */
            z = flags & (DT_OPAQUE_TAG_MAX - 1);
	    if((z&0x7) != 0) { /* must be aligned */
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Opaque:Tag must be aligned", 
		  logical, NULL);
		CK_SET_RET(FAIL)
	    }

	    if((flags>>8) != 0) {  /* should be reserved(zero) */
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Opaque:Bits 8-23 should be 0 for class bit field", 
		  logical, NULL);
		CK_SET_RET(FAIL)
	    }

            if((dt->shared->u.opaque.tag=malloc(z+1))==NULL) {
		error_push(ERR_INTERNAL, ERR_NONE_SEC, "Datatype Message:Opaque:Internal allocation error", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }

            memcpy(dt->shared->u.opaque.tag, *pp, z);
            dt->shared->u.opaque.tag[z] = '\0';
            *pp += z;
            break;

	case DT_COMPOUND:  /* Compound datatypes */
	{
	    unsigned 	offset_nbytes;   
	    unsigned 	j;

            /* Compute the # of bytes required to store a member offset */
            offset_nbytes = (V_log2_gen((uint64_t)dt->shared->size) + 7) / 8;

	    dt->shared->u.compnd.nmembs = flags & 0xffff;
	    if(dt->shared->u.compnd.nmembs <= 0) {
	       error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Compound:Number of members should be greater than 0", logical, NULL);
	       CK_SET_RET_DONE(FAIL)
	    }

            dt->shared->u.compnd.packed = TRUE; /* Start off packed */
            dt->shared->u.compnd.nalloc = dt->shared->u.compnd.nmembs;
            dt->shared->u.compnd.memb = (DT_cmemb_t *)calloc(dt->shared->u.compnd.nalloc, sizeof(DT_cmemb_t));

            if(NULL == dt->shared->u.compnd.memb) {
		error_push(ERR_INTERNAL, ERR_NONE_SEC,
		    "Datatype Message:Compound:Internal allocation error", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }
            for(i = 0; i < dt->shared->u.compnd.nmembs; i++) {
		unsigned ndims = 0;     /* Number of dimensions of the array field */
                ck_hsize_t dim[OBJ_LAYOUT_NDIMS];  /* Dimensions of the array */
                type_t *array_dt;    /* Temporary pointer to the array datatype */
                type_t *temp_type;   /* Temporary pointer to the field's datatype */

                dt->shared->u.compnd.memb[i].name = strdup((const char *)*pp);

                /* Version 3 of the datatype message eliminated the padding to multiple of 8 bytes */
                if(version >= DT_VERSION_3)
		    *pp += strlen((const char *)*pp) + 1;
                else
		    *pp += ((strlen((const char *)*pp) + 8) / 8) * 8;

                if(version >= DT_VERSION_3) {
		   UINT32DECODE_VAR(*pp, dt->shared->u.compnd.memb[i].offset, offset_nbytes);
                } else {
                   UINT32DECODE(*pp, dt->shared->u.compnd.memb[i].offset);
		}

		if(version == DT_VERSION_1) {
		    /* Decode the number of dimensions */
		    ndims = *(*pp)++;
		    if(ndims > 4) {
		       error_push(ERR_LEV_2, ERR_LEV_2A2d, 
			    "Datatype Message:Compound:Number of dimensions should not exceed 4 for version 1", logical, NULL);
		       CK_SET_RET_DONE(FAIL)
		    }
		    *pp += 3;  /*reserved bytes */
		    *pp += 4;  /* skip dimension permutation */
		    *pp += 4;  /* reserved bytes */

                    /* Decode array dimension sizes */
                    for(j = 0; j < 4; j++)
			UINT32DECODE(*pp, dim[j]);
                } /* end if */

		/* Allocate space for the field's datatype */
		if((temp_type = type_alloc(logical)) == NULL) {
		    error_push(ERR_LEV_2, ERR_LEV_2A2d, 
			"Datatype Message:Compound:Internal allocation", logical, NULL);
		    CK_SET_RET_DONE(FAIL)
		}

                /* Decode the field's datatype information */
                if(OBJ_dt_decode_helper(file, pp, temp_type, start_buf, logi_base) != SUCCEED) {
                    for (j = 0; j <= i; j++)
                        free(dt->shared->u.compnd.memb[j].name);
                    free(dt->shared->u.compnd.memb);
                    error_push(ERR_LEV_2, ERR_LEV_2A2d,
		      "Datatype Message:Unable to decode Compound member type", logical, NULL);
		    CK_SET_RET_DONE(FAIL)
                }

		dt->shared->u.compnd.memb[i].size = temp_type->shared->size;
                dt->shared->u.compnd.memb[i].type = temp_type;
            } /* end for */
	    
	}
	    break;


	case DT_REFERENCE: /* Reference datatypes...  */
            dt->shared->u.atomic.order = DT_ORDER_NONE;
            dt->shared->u.atomic.prec = 8 * dt->shared->size;            
	    dt->shared->u.atomic.offset = 0;
            dt->shared->u.atomic.lsb_pad = DT_PAD_ZERO;
            dt->shared->u.atomic.msb_pad = DT_PAD_ZERO;

            /* Set reference type */
            dt->shared->u.atomic.u.r.rtype = (DTR_type_t)(flags & 0x0f);
	    /* Value of 2 is not implemented yet (see spec.) */
	    /* value of 3-15 is supposedly to be reserved */
	    if((flags&0x0f) >= 2) {
		error_push(ERR_LEV_2, ERR_LEV_2A2d, "Datatype Message:Reference:Invalid class bit field", logical, NULL);
		CK_SET_RET(FAIL)
	    }

	    if((flags>>4) != 0) {
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Reference:Bits 4-23 should be 0 for class bit field", logical, NULL);
		CK_SET_RET(FAIL)
	    }
            break;

	case DT_ENUM:  /* Enumeration datatypes */
	    dt->shared->u.enumer.nmembs = dt->shared->u.enumer.nalloc = flags & 0xffff;
	    if((flags>>16) != 0) {
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Enumeration:Bits 16-23 should be 0 for class bit field", 
		  logical, NULL);
		CK_SET_RET(FAIL)
	    }

	    if((dt->shared->parent=type_alloc(logical)) == NULL) {
              	error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Enumeration:Internal allocation error", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }

            if(OBJ_dt_decode_helper(file, pp, dt->shared->parent, start_buf, logi_base) != SUCCEED) {
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		    "Datatype Message:Unable to decode enumeration parent type", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }

	    if((dt->shared->u.enumer.name=calloc(dt->shared->u.enumer.nalloc, sizeof(char*)))==NULL) {
              	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
		  "Datatype Message:Enumeration:Internal allocation error", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }
            if((dt->shared->u.enumer.value=calloc(dt->shared->u.enumer.nalloc,
                 dt->shared->parent->shared->size))==NULL) {
		error_push(ERR_INTERNAL, ERR_NONE_SEC, 
		    "Datatype Message:Enumeration:Internal allocation error", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }


            /* Names, each a multiple of 8 with null termination */
            for(i=0; i<dt->shared->u.enumer.nmembs; i++) {
		size_t	tmp_len;

		 dt->shared->u.enumer.name[i] = strdup((const char*)*pp);

                /* Version 3 of the datatype message eliminated the padding to multiple of 8 bytes */
                if(version >= DT_VERSION_3)
                    *pp += strlen((const char *)*pp) + 1;
                else
                    /* Advance multiple of 8 w/ null terminator */
                    *pp += ((strlen((const char *)*pp) + 8) / 8) * 8;
            }

            /* Values */
            memcpy(dt->shared->u.enumer.value, *pp,
		dt->shared->u.enumer.nmembs * dt->shared->parent->shared->size);
            *pp += dt->shared->u.enumer.nmembs * dt->shared->parent->shared->size;
            break;

	case DT_VLEN:  /* Variable length datatypes */
	    /* Set the type of VL information, either sequence or string */
            dt->shared->u.vlen.type = (DT_vlen_type_t)(flags & 0x0f);
	    if((dt->shared->u.vlen.type != DT_VLEN_STRING) &&
	       (dt->shared->u.vlen.type != DT_VLEN_SEQUENCE)) {
		    error_push(ERR_LEV_2, ERR_LEV_2A2d, 
			"Datatype Message:Variable Length:Unsupported variable length datatype", logical, NULL);
		    CK_SET_RET(FAIL)
	    }


            if(dt->shared->u.vlen.type == DT_VLEN_STRING) {
                dt->shared->u.vlen.pad  = (DT_str_t)((flags>>4) & 0x0f);
		if((dt->shared->u.vlen.pad !=DT_STR_NULLTERM) &&
		   (dt->shared->u.vlen.pad !=DT_STR_NULLPAD) &&
		   (dt->shared->u.vlen.pad !=DT_STR_SPACEPAD)) {
		   error_push(ERR_LEV_2, ERR_LEV_2A2d, 
			"Datatype Message:Variable Length:Unsupported padding type", logical, NULL);
		   CK_SET_RET(FAIL)
		}

                dt->shared->u.vlen.cset = (DT_cset_t)((flags>>8) & 0x0f);
		if((dt->shared->u.vlen.cset != DT_CSET_ASCII) && 
		   (dt->shared->u.vlen.cset != DT_CSET_UTF8)) {
		   error_push(ERR_LEV_2, ERR_LEV_2A2d, 
			"Datatype Message:Variable Length:Unsupported character set", logical, NULL);
		   CK_SET_RET(FAIL)
		}
            } /* end if */


	    if((flags>>12) != 0) {
	       error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Variable-Length:Bits 12-23 should be 0 for class bit field", 
		  logical, NULL);
	       CK_SET_RET(FAIL)
	    }

            /* Decode base type of VL information */
            if((dt->shared->parent = type_alloc(logical))==NULL) {
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Variable Length:Internal allocation error", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }

            if(OBJ_dt_decode_helper(file, pp, dt->shared->parent, start_buf, logi_base) != SUCCEED) {
	       error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Unable to decode variable-length parent type", logical, NULL);
	       CK_SET_RET_DONE(FAIL)
	    }

            break;


	case DT_ARRAY:  /* Array datatypes */
	    /* Decode the number of dimensions */
            dt->shared->u.array.ndims = *(*pp)++;

            /* Double-check the number of dimensions */
            if(dt->shared->u.array.ndims > OBJ_SDS_MAX_RANK) {
	       error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Array:Dimension exceeds limit", logical, NULL);
	       CK_SET_RET_DONE(FAIL)
	    }

            /* Skip reserved bytes */
	    if(version < DT_VERSION_3)
		*pp += 3;

            /* Decode array dimension sizes & compute number of elements */
            for(j=0, dt->shared->u.array.nelem=1; j<(unsigned)dt->shared->u.array.ndims; j++) {
                UINT32DECODE(*pp, dt->shared->u.array.dim[j]);
                dt->shared->u.array.nelem *= dt->shared->u.array.dim[j];
            } /* end for */

	    /* Skip array dimension permutations, if version has them */
            if(version < DT_VERSION_3)
                *pp += dt->shared->u.array.ndims * 4;

            /* Decode base type of array */
            if((dt->shared->parent = type_alloc(logical))==NULL) {
		error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Array:Internal allocation error", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }
            if(OBJ_dt_decode_helper(file, pp, dt->shared->parent, start_buf, logi_base) != SUCCEED) {
	       error_push(ERR_LEV_2, ERR_LEV_2A2d, 
		  "Datatype Message:Unable to decode Array parent type", logical, NULL);
	       CK_SET_RET_DONE(FAIL)
	    }

            break;

	default:
	    /* shouldn't come here */
	    printf("OBJ_dt_decode_helper: datatype class not handled yet...type=%d\n", dt->shared->type);
    }

done:
    return(ret_value);
} /* OBJ_dt_decode_helper() */


/* Datatype */
static void *
OBJ_dt_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    type_t	*dt=NULL;
    void    	*ret_ptr=NULL;
    int		ret_value=SUCCEED;
    ck_addr_t	logical;

    assert(file);
    assert(p);

    logical = get_logical_addr(p, start_buf, logi_base);
    if((dt = calloc(1, sizeof(type_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Datatype Message:Internal allocation error", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* NOT SURE what to use with the group entry yet */
    memset(&(dt->ent), 0, sizeof(GP_entry_t));
    dt->ent.header = CK_ADDR_UNDEF;
    if((dt->shared = calloc(1, sizeof(DT_shared_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Datatype Message:Internal allocation error", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    if(OBJ_dt_decode_helper(file, &p, dt, start_buf, logi_base) != SUCCEED) {
	/* should have errors already pushed onto the error stack from decode_helper() */
	CK_SET_RET_DONE(FAIL)
    }

    if(ret_value == SUCCEED)
	ret_ptr = dt;
done:
    if(ret_value == FAIL && dt) {
	/* NEED: probably have more to be freed see decode_helper() */
	if (dt->shared) free(dt->shared);
	free(dt);
    }

    return(ret_ptr);
}  /* OBJ_dt_decode() */

/* Data Storage - Fill Value (old) */
static void *
OBJ_fill_old_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_fill_t 	*mesg=NULL;            /* Decoded fill value message */
    void 	*ret_ptr=NULL;                  
    int         ret_value=SUCCEED;
    ck_addr_t   logical;
    int         version, badinfo;

    assert(file);
    assert(p);

    logical = get_logical_addr(p, start_buf, logi_base);
    mesg = calloc(1, sizeof(OBJ_fill_t));
    if(mesg == NULL) {
        error_push(ERR_INTERNAL, ERR_NONE_SEC,
          "Fill Value (old) Message:Internal allocation error", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Set non-zero default fields */
    mesg->version = OBJ_FILL_VERSION_2;
    mesg->alloc_time = FILL_ALLOC_TIME_LATE;
    mesg->fill_time = FILL_TIME_IFSET;

    /* Fill value size */
    UINT32DECODE(p, mesg->size);

    /* Only decode the fill value itself if there is one */
    if(mesg->size > 0) {
        if(NULL == (mesg->buf = malloc((ck_size_t)mesg->size))) {
	    error_push(ERR_INTERNAL, ERR_NONE_SEC,
		"Fill Value (old) Message:Internal allocation error", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
        memcpy(mesg->buf, p, (ck_size_t)mesg->size);
        mesg->fill_defined = TRUE;
    } /* end if */
    else
        mesg->size = (-1);

    if(ret_value == SUCCEED)
        ret_ptr = (void*)mesg;

done:
    if(ret_value == FAIL && mesg) {
	if (mesg->buf) free(mesg->buf);
	free(mesg);
    }

    return(ret_ptr);
} /* OBJ_fill_old_decode() */


/* Data Storage - Fill Value */
static void *
OBJ_fill_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_fill_t  *mesg=NULL;
    void        *ret_ptr=NULL;
    int		ret_value=SUCCEED;
    ck_addr_t	logical;
    int		badinfo;
    unsigned	version;

    assert(file);
    assert(p);

    logical = get_logical_addr(p, start_buf, logi_base);
    mesg = calloc(1, sizeof(OBJ_fill_t));
    if(mesg == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	  "Fill Value Message:Internal allocation error", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    mesg->version = version = *p++;

    /* version is adjusted if error for further continuation */
    if(g_format_num == FORMAT_ONE_SIX) {
	if(version != OBJ_FILL_VERSION && version != OBJ_FILL_VERSION_2) {
	    badinfo = version;
	    version = OBJ_FILL_VERSION_2;
	    error_push(ERR_LEV_2, ERR_LEV_2A2f, "FIll Value Message:Bad version number", logical, &badinfo);
		CK_SET_RET(FAIL)
	}
    } else {
	if(version < OBJ_FILL_VERSION || version > OBJ_FILL_VERSION_LATEST) {
	    badinfo = version;
	    version = OBJ_FILL_VERSION_LATEST;
	    error_push(ERR_LEV_2, ERR_LEV_2A2f, "FIll Value Message:Bad version number", logical, &badinfo);
		CK_SET_RET(FAIL)
	}
    }

    if(version < OBJ_FILL_VERSION_3) {
       /* Space allocation time */
       logical = get_logical_addr(p, start_buf, logi_base);
       mesg->alloc_time = (fill_alloc_time_t)*p++;
       if((mesg->alloc_time != FILL_ALLOC_TIME_EARLY) &&
	  (mesg->alloc_time != FILL_ALLOC_TIME_LATE) &&
	  (mesg->alloc_time != FILL_ALLOC_TIME_INCR)) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2f,
		"Fill Value Message:Invalid Space Allocation Time", logical, NULL);
	    CK_SET_RET(FAIL)
	}

        /* Fill value write time */
	logical = get_logical_addr(p, start_buf, logi_base);
        mesg->fill_time = (fill_time_t)*p++;
	if((mesg->fill_time != FILL_TIME_ALLOC) &&
	   (mesg->fill_time != FILL_TIME_NEVER) &&
	   (mesg->fill_time != FILL_TIME_IFSET)) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2f, "Fill Value Message:Invalid Fill Value Write Time", logical, NULL);
	    CK_SET_RET(FAIL)
	}

        /* Whether fill value is defined */
	logical = get_logical_addr(p, start_buf, logi_base);
        mesg->fill_defined = *p++;
	if((mesg->fill_defined != 0) && (mesg->fill_defined != 1)) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2f, "Fill Value Message:Invalid Fill Value Defined", logical, NULL);
	    CK_SET_RET(FAIL)
	}

        /* Only decode fill value information if one is defined */
        if(mesg->fill_defined) {
	    logical = get_logical_addr(p, start_buf, logi_base);
            INT32DECODE(p, mesg->size);
            if(mesg->size > 0) {
	       CHECK_OVERFLOW(mesg->size, ssize_t, ck_size_t);
	       if(NULL == (mesg->buf=malloc((ck_size_t)mesg->size))) {
		    error_push(ERR_INTERNAL, ERR_NONE_SEC, "Fill Value Message:Internal allocation error", logical, NULL);
		    CK_SET_RET_DONE(FAIL)
		}
		memcpy(mesg->buf, p, (ck_size_t)mesg->size);
            } 
        } else
            mesg->size = (-1);
    } else {
	unsigned flags;

        /* Flags */
	logical = get_logical_addr(p, start_buf, logi_base);
        flags = *p++;

        /* Check for unknown flags */
        if(flags & (unsigned)~OBJ_FILL_FLAGS_ALL) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2f, "Fill Value Message:Unknown flag", logical, NULL);
	    CK_SET_RET(FAIL)
	}

        /* Space allocation time */
        mesg->alloc_time = (flags >> OBJ_FILL_SHIFT_ALLOC_TIME) & OBJ_FILL_MASK_ALLOC_TIME;

        /* Fill value write time */
        mesg->fill_time = (flags >> OBJ_FILL_SHIFT_FILL_TIME) & OBJ_FILL_MASK_FILL_TIME;

        /* Check for undefined fill value */
        if(flags & OBJ_FILL_FLAG_UNDEFINED_VALUE) {
            if (flags & OBJ_FILL_FLAG_HAVE_VALUE) {
		error_push(ERR_LEV_2, ERR_LEV_2A2f, "Fill Value Message:Invalid Fill Value Defined", logical, NULL);
		CK_SET_RET(FAIL)
	    }

            /* Set value for "undefined" fill value */
            mesg->size = (-1);
        } else if(flags & OBJ_FILL_FLAG_HAVE_VALUE) {
            /* Fill value size */
	    logical = get_logical_addr(p, start_buf, logi_base);
            UINT32DECODE(p, mesg->size);

            /* Fill value */
            CHECK_OVERFLOW(mesg->size, ssize_t, ck_size_t);
            if(NULL == (mesg->buf = malloc((ck_size_t)mesg->size))) {
		error_push(ERR_INTERNAL, ERR_NONE_SEC, "Fill Value Message:Internal allocation error", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }
            memcpy(mesg->buf, p, (ck_size_t)mesg->size);

            /* Set the "defined" flag */
            mesg->fill_defined = TRUE;
        } else {
            /* Set the "defined" flag */
            mesg->fill_defined = TRUE;
        } 
    } /* end else */

    if(ret_value == SUCCEED)
	ret_ptr = (void*)mesg;
done:
    if(ret_value == FAIL && mesg) {
	if(mesg->buf) free(mesg->buf);
	free(mesg);
    }

    return(ret_ptr);
} /* OBJ_fill_decode() */

/* Link Message */
static void *
OBJ_link_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_link_t          *lnk=NULL;    /* Pointer to link message */
    ck_size_t           len;          /* Length of a string in the message */
    unsigned char       link_flags;   /* Flags for encoding link info */
    int			badinfo;
    unsigned        	version;
    void                *ret_ptr=NULL;
    int			ret_value=SUCCEED;
    ck_addr_t		logical;

    assert(file);
    assert(p);

    logical = get_logical_addr(p, start_buf, logi_base);
    if(g_format_num == FORMAT_ONE_SIX) {
	error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Unsupported message", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Allocate space for message */
    if((lnk = calloc(1, sizeof (OBJ_link_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Link Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* decode */
    logical = get_logical_addr(p, start_buf, logi_base);
    version = *p++;
    if(version != OBJ_LINK_VERSION) {
	badinfo = version;
	error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Bad version number", logical, &badinfo);
	CK_SET_RET(FAIL)
    }

    /* Get the encoding flags for the link */
    logical = get_logical_addr(p, start_buf, logi_base);
    link_flags = *p++;
    if(link_flags & ~OBJ_LINK_ALL_FLAGS) {
	error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Bad Flag Value", logical, NULL);
	CK_SET_RET(FAIL)
    }

    /* Check for non-default link type */
    logical = get_logical_addr(p, start_buf, logi_base);
    if(link_flags & OBJ_LINK_STORE_LINK_TYPE) {
        /* Get the type of the link */
        lnk->type = *p++;
        if(lnk->type < L_TYPE_HARD || lnk->type > L_TYPE_MAX) {
	   badinfo = lnk->type;
	   error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Bad Link Type", logical, &badinfo);
	   CK_SET_RET_DONE(FAIL)
	}
    } else
        lnk->type = L_TYPE_HARD;

    /* Get the link creation time from the file */
    logical = get_logical_addr(p, start_buf, logi_base);
    if(link_flags & OBJ_LINK_STORE_CORDER) {
        INT64DECODE(p, lnk->corder);
        lnk->corder_valid = TRUE;
    } /* end if */
    else {
        lnk->corder = 0;
        lnk->corder_valid = FALSE;
    } /* end else */

    /* Check for non-default name character set */ 
    logical = get_logical_addr(p, start_buf, logi_base);
    if(link_flags & OBJ_LINK_STORE_NAME_CSET) {
        /* Get the link name's character set */
        lnk->cset = (DT_cset_t)*p++;   
        if(lnk->cset < DT_CSET_ASCII || lnk->cset > DT_CSET_UTF8) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Invalid character set for link name", logical, &badinfo);
	    CK_SET_RET(FAIL)
	}
    } /* end if */
    else
        lnk->cset = DT_CSET_ASCII;
    
    /* Get the length of the link's name */
    logical = get_logical_addr(p, start_buf, logi_base);
    switch(link_flags & OBJ_LINK_NAME_SIZE) {
        case 0:     /* 1 byte size */
            len = *p++;
            break;
    
        case 1:     /* 2 byte size */
            UINT16DECODE(p, len); 
            break;
    
        case 2:     /* 4 byte size */
            UINT32DECODE(p, len);
            break;

        case 3:     /* 8 byte size */
            UINT64DECODE(p, len);
            break; 
        
        default:
            assert(0 && "bad size for name"); 
    } /* end switch */

    if(len == 0) {
	error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Invalid name length for link", logical, &badinfo);
	CK_SET_RET_DONE(FAIL)
    }

    /* Get the link's name */
    if(NULL == (lnk->name = malloc(len + 1))) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Link Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    memcpy(lnk->name, p, len);
    lnk->name[len] = '\0';
    p += len;

    /* Get the appropriate information for each type of link */
    logical = get_logical_addr(p, start_buf, logi_base);
    switch(lnk->type) {
        case L_TYPE_HARD:
            /* Get the address of the object the link points to */
            addr_decode(file->shared, &p, &(lnk->u.hard.addr));
            break;

        case L_TYPE_SOFT:
            /* Get the link value */
            UINT16DECODE(p, len)
            if(len == 0) {
		error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Invalid name length for link", logical, &badinfo);
		CK_SET_RET_DONE(FAIL)
	    }
            if(NULL == (lnk->u.soft.name = malloc((ck_size_t)len + 1))) {
		error_push(ERR_INTERNAL, ERR_NONE_SEC, "Link Message:Internal allocation error", -1, NULL);
		CK_SET_RET_DONE(FAIL)
	    }
            memcpy(lnk->u.soft.name, p, len);
            lnk->u.soft.name[len] = '\0';
            p += len;
            break;

        /* External & User-defined links */
        default: /* user-defined link */
            if(lnk->type < L_TYPE_UD_MIN || lnk->type > L_TYPE_MAX) {
		error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Invalid user-defined link type", logical, &badinfo);
		CK_SET_RET(FAIL)
	    }

            /* A UD link.  Get the user-supplied data */
	    logical = get_logical_addr(p, start_buf, logi_base);
            UINT16DECODE(p, len)
            lnk->u.ud.size = len;
            if(len > 0) {
		if(NULL == (lnk->u.ud.udata = malloc((ck_size_t)len))) {
		    error_push(ERR_INTERNAL, ERR_NONE_SEC, "Link Message:Internal allocation error", -1, NULL);
		    CK_SET_RET_DONE(FAIL)
		}

		logical = get_logical_addr(p, start_buf, logi_base);
                memcpy(lnk->u.ud.udata, p, len);
		if(lnk->type == L_TYPE_EXTERNAL) { /* external link */
		    uint8_t    	*s; 
		    char 	*file_name, *obj_name;
		    ck_size_t	fname_len, obj_len;

		    s = lnk->u.ud.udata;
		    /* Check external link version & flags */
		    if(((*s >> 4) & 0x0F) > L_EXT_VERSION) {
			error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Bad version # for external link type", 
			    logical, &badinfo);
			CK_SET_RET(FAIL)
		    }
		    if((*s & 0x0F) & ~L_EXT_FLAGS_ALL) {
			error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Bad flags for external link type", 
			    logical, &badinfo);
			CK_SET_RET(FAIL)
		    }
		    s++;
		    file_name = (char *)s;
		    fname_len = strlen((const char *)file_name) + 1;
		    if(1 + fname_len > len) {
			error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Invalid file length for external link type", 
			    logical, &badinfo);
			CK_SET_RET(FAIL)
		    }
		    obj_name = (char *)s + fname_len;
		    obj_len = strlen((const char *)obj_name) + 1;

		    if ((1 + fname_len + obj_len) > len) {
			error_push(ERR_LEV_2, ERR_LEV_2A2g, "Link Message:Invalid object length for external link type", 
			    logical, &badinfo);
			CK_SET_RET(FAIL)
		    }
		} else
		    p += len;
            } else
                lnk->u.ud.udata = NULL;
    } /* end switch */

    if (ret_value == SUCCEED)
	ret_ptr = lnk;

done:
    if(ret_value == FAIL && lnk) {
	if(lnk->name != NULL)
	    free(lnk->name);
	if(lnk->type == L_TYPE_SOFT && lnk->u.soft.name != NULL)
	    free(lnk->u.soft.name);
	if(lnk->type >= L_TYPE_UD_MIN && lnk->u.ud.size > 0 && lnk->u.ud.udata != NULL)
	    free(lnk->u.ud.udata);
	free(lnk);
    }

    return(ret_ptr);
} /* OBJ_link_decode() */


/* Data Storage - External Data Files */
static void *
OBJ_edf_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_edf_t   *mesg=NULL;
    int         version, badinfo;
    const char  *s=NULL;
    size_t      u;      	
    void 	*ret_ptr=NULL;    
    int		ret_value=SUCCEED;
    ck_addr_t	logical;

    assert(file);
    assert(p);

    logical = get_logical_addr(p, start_buf, logi_base);
    if((mesg = calloc(1, sizeof(OBJ_edf_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "External Data Files Message:Internal allocation error", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Version */
    logical = get_logical_addr(p, start_buf, logi_base);
    version = *p++;
    if(version != OBJ_EDF_VERSION) {
       badinfo = version;
       error_push(ERR_LEV_2, ERR_LEV_2A2h, "External Data Files Message:Bad version number", logical, &badinfo);
       CK_SET_RET(FAIL)
    }

    /* Reserved */
    p += 3;

    /* Number of slots */
    logical = get_logical_addr(p, start_buf, logi_base);
    UINT16DECODE(p, mesg->nalloc);
    assert(mesg->nalloc>0);
    UINT16DECODE(p, mesg->nused);
    assert(mesg->nused <= mesg->nalloc);

    if(!(mesg->nalloc >= mesg->nused)) {
	error_push(ERR_LEV_2, ERR_LEV_2A2h, 
	    "External Data Files Message:Inconsistent number of Allocated Slots", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    addr_decode(file->shared, &p, &(mesg->heap_addr));
    if(!(addr_defined(mesg->heap_addr))) {
       error_push(ERR_LEV_2, ERR_LEV_2A2h, 
	    "External Data Files Message:Undefined heap address", logical, NULL);
       CK_SET_RET(FAIL)
    }

    /* Decode the file list */
    logical = get_logical_addr(p, start_buf, logi_base);
    mesg->slot = calloc(mesg->nalloc, sizeof(OBJ_edf_entry_t));
    if(NULL == mesg->slot) {
       error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "External Data Files Message:Internal allocation error", logical, NULL);
       CK_SET_RET_DONE(FAIL)
    }

    for(u = 0; u < mesg->nused; u++) {
	/* Name offset: will be validated later in decode_validate_messages() */
	DECODE_LENGTH(file->shared, p, mesg->slot[u].name_offset);
        /* File offset */
        DECODE_LENGTH(file->shared, p, mesg->slot[u].offset);
	if(!(mesg->slot[u].offset >=0)) {
	   error_push(ERR_LEV_2, ERR_LEV_2A2h, 
		"External Data Files Messag:Invalid file offset", logical, NULL);
	   CK_SET_RET(FAIL)	
	}
        /* size */
        DECODE_LENGTH(file->shared, p, mesg->slot[u].size);
	if(mesg->slot[u].size <= 0) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2h, "External Data Files Message:Invalid size", logical, NULL);
	    CK_SET_RET(FAIL)	
	}
	if(!(mesg->slot[u].offset <= mesg->slot[u].size)) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2h, 
		"External data Files Message:Inconsistent file offset/size", logical, NULL);
	    CK_SET_RET(FAIL)	
	}
    }

    if(ret_value == SUCCEED)
	ret_ptr = mesg;
done:
    if(ret_value == FAIL && mesg) {
	if(mesg->slot)	 
	    free(mesg->slot);
    }
    return(ret_ptr);
} /* OBJ_edf_decode() */


/* Data Storage - Layout */
static void *
OBJ_layout_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_layout_t    	*mesg=NULL;
    unsigned        	u;
    void            	*ret_ptr=NULL;
    int			ret_value=SUCCEED;
    ck_addr_t		logical;
    int			badinfo;
    unsigned		version;

    assert(file);
    assert(p);

    logical = get_logical_addr(p, start_buf, logi_base);

    if((mesg=calloc(1, sizeof(OBJ_layout_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Layout Message:Internal allocation error", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Version. 1 when space allocated; 2 when space allocation is delayed */
    logical = get_logical_addr(p, start_buf, logi_base);
    mesg->version = version = *p++;

    /* version is adjusted if error for further continuation */
    if(version < OBJ_LAYOUT_VERSION_1 || version>OBJ_LAYOUT_VERSION_3) {
	badinfo = version;
	version = OBJ_LAYOUT_VERSION_3;
	error_push(ERR_LEV_2, ERR_LEV_2A2i, "Layout Message:Bad version number", logical, &badinfo);
	CK_SET_RET(FAIL)
    }

    if(version < OBJ_LAYOUT_VERSION_3) {  /* version 1 & 2 */
        unsigned        ndims;    /* Num dimensions in chunk */

        /* Dimensionality */
	logical = get_logical_addr(p, start_buf, logi_base);
        ndims = *p++;
	/* this is not specified in the format specification ??? */
        if(ndims > OBJ_LAYOUT_NDIMS) {
	    badinfo = ndims;
	    error_push(ERR_LEV_2, ERR_LEV_2A2i, 
		"Layout Message:Dimensionality is too large", logical, &badinfo);
	    CK_SET_RET(FAIL)
	}

        /* Layout class */
	logical = get_logical_addr(p, start_buf, logi_base);
        mesg->type = (DATA_layout_t)*p++;

        if((DATA_CONTIGUOUS != mesg->type) && (DATA_CHUNKED != mesg->type) && (DATA_COMPACT != mesg->type)) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2i, "Layout Message:invalid layout class", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}

        /* Reserved bytes */
        p += 5;

       	/* Address */
       	if(mesg->type == DATA_CONTIGUOUS) {
	    addr_decode(file->shared, &p, &(mesg->u.contig.addr));
       	} else if(mesg->type==DATA_CHUNKED) {
	    addr_decode(file->shared, &p, &(mesg->u.chunk.addr));
	}

        /* Read the size */
        if(mesg->type!=DATA_CHUNKED) {
	    mesg->unused.ndims = ndims;
            for (u = 0; u < ndims; u++)
		UINT32DECODE(p, mesg->unused.dim[u]);
        } else {
	    mesg->u.chunk.ndims=ndims;
	    for (u = 0; u < ndims; u++) {
               	UINT32DECODE(p, mesg->u.chunk.dim[u]);
	    }

            /* Compute chunk size */
            for(u=1, mesg->u.chunk.size=mesg->u.chunk.dim[0]; u<ndims; u++)
		mesg->u.chunk.size *= mesg->u.chunk.dim[u];
        }

       	if (mesg->type == DATA_COMPACT) {
	    logical = get_logical_addr(p, start_buf, logi_base);
            UINT32DECODE(p, mesg->u.compact.size);
	    if(mesg->u.compact.size > 0) {
		if(NULL==(mesg->u.compact.buf=malloc(mesg->u.compact.size))) {
		    error_push(ERR_INTERNAL, ERR_NONE_SEC, 
			"Layout Message:Internal allocation error", -1, NULL);
		    CK_SET_RET_DONE(FAIL)
		}
                memcpy(mesg->u.compact.buf, p, mesg->u.compact.size);
                p += mesg->u.compact.size;
	    }
        }
    } else {  /* version 3 */
        /* Layout class */
	logical = get_logical_addr(p, start_buf, logi_base);
        mesg->type = (DATA_layout_t)*p++;

        /* Interpret the rest of the message according to the layout class */
        switch(mesg->type) {
	    case DATA_CONTIGUOUS:
                addr_decode(file->shared, &p, &(mesg->u.contig.addr));
		DECODE_LENGTH(file->shared, p, mesg->u.contig.size);
                break;

            case DATA_CHUNKED:
                /* Dimensionality */
		logical = get_logical_addr(p, start_buf, logi_base);
                mesg->u.chunk.ndims = *p++;
                if(mesg->u.chunk.ndims>OBJ_LAYOUT_NDIMS) {
		   badinfo = mesg->u.chunk.ndims;
		   error_push(ERR_LEV_2, ERR_LEV_2A2i, 
			"Layout Message:Chunked layout:Dimensionality is too large", 
			logical, &badinfo);
		   CK_SET_RET(FAIL)
		}

                /* B-tree address */
                addr_decode(file->shared, &p, &(mesg->u.chunk.addr));

                /* Chunk dimensions */
                for(u = 0; u < mesg->u.chunk.ndims; u++)
		    UINT32DECODE(p, mesg->u.chunk.dim[u]);

                /* Compute chunk size */
                for (u=1, mesg->u.chunk.size=mesg->u.chunk.dim[0]; u<mesg->u.chunk.ndims; u++)
		    mesg->u.chunk.size *= mesg->u.chunk.dim[u];
                break;

	    case DATA_COMPACT:
                UINT16DECODE(p, mesg->u.compact.size);
                if(mesg->u.compact.size > 0) {
		    if(NULL==(mesg->u.compact.buf=malloc(mesg->u.compact.size))) {
			error_push(ERR_INTERNAL, ERR_NONE_SEC, 
			    "Layout Message:Compact layout:Internal allocation error", -1, NULL);
			CK_SET_RET_DONE(FAIL)
		    }
		    memcpy(mesg->u.compact.buf, p, mesg->u.compact.size);
		    p += mesg->u.compact.size;
		} /* end if */
                break;

	    default:
		error_push(ERR_LEV_2, ERR_LEV_2A2i, 
		  "Layout Message:Invalid Layout Class", logical, NULL);
		CK_SET_RET(FAIL)
	} /* end switch */
    } /* version 3 */

    if(ret_value == SUCCEED)
	ret_ptr = mesg;
done:
    if(ret_value == FAIL && mesg) {
	if(mesg->u.compact.buf) free(mesg->u.compact.buf);
	free(mesg);
    }

    return(ret_ptr);
} /* OBJ_layout_decode() */

/* Bogus Message */
static void *
OBJ_bogus_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_bogus_t 	*mesg = NULL;
    int                 ret_value=SUCCEED;
    void                *ret_ptr=NULL;
    ck_addr_t		logical;

    assert(file);
    assert(p);

    if((mesg = calloc(1, sizeof(OBJ_bogus_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Bogus Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* decode */
    logical = get_logical_addr(p, start_buf, logi_base);
    UINT32DECODE(p, mesg->u);

    /* Validate the bogus info */
    if(mesg->u != OBJ_BOGUS_VALUE) {
	error_push(ERR_LEV_2, ERR_LEV_2A2k, "Bogus Message:Invalid bogus value", logical, NULL);
	CK_SET_RET(FAIL)
    }

    if(ret_value == SUCCEED)
	ret_ptr = mesg;
done:
    if(ret_value == FAIL && mesg)
	free(mesg);

    return(ret_ptr);
}  /* OBJ_bogus_decode() */

/* Group Info Message */
static void *
OBJ_ginfo_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_ginfo_t         *mesg=NULL;  		/* Pointer to group information message */
    unsigned char       flags;          	/* Flags for encoding group info */
    ck_addr_t		logical;
    int			version, badinfo;
    int			ret_value=SUCCEED;
    void                *ret_ptr=NULL;

    assert(file);
    assert(p);

    if(g_format_num == FORMAT_ONE_SIX) {
	error_push(ERR_LEV_2, ERR_LEV_2A2k, "Group Info Message:Unsupported message", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    if((mesg = calloc(1, sizeof(OBJ_ginfo_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Group Info Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    version = *p++;
    if(version != OBJ_GINFO_VERSION) {
	badinfo = version;
	error_push(ERR_LEV_2, ERR_LEV_2A2k, 
	    "Group Info Message:Bad version number", logical, &badinfo);
	CK_SET_RET(FAIL)
    }


    /* Get the flags for the group */
    logical = get_logical_addr(p, start_buf, logi_base);
    flags = *p++;
    if(flags & ~OBJ_GINFO_ALL_FLAGS) {
	error_push(ERR_LEV_2, ERR_LEV_2A2k, "Group Info Message:Bad flag value", logical, NULL);
	CK_SET_RET(FAIL)
    }

    mesg->store_link_phase_change = (flags & OBJ_GINFO_STORE_PHASE_CHANGE) ? TRUE : FALSE;
    mesg->store_est_entry_info = (flags & OBJ_GINFO_STORE_EST_ENTRY_INFO) ? TRUE : FALSE;

    if(mesg->store_link_phase_change) {
        UINT16DECODE(p, mesg->max_compact)
        UINT16DECODE(p, mesg->min_dense)
    } else {
        mesg->max_compact = OBJ_CRT_GINFO_MAX_COMPACT;
        mesg->min_dense = OBJ_CRT_GINFO_MIN_DENSE;
    }

    /* Get the estimated # of entries & name lengths */
    if(mesg->store_est_entry_info) {
        UINT16DECODE(p, mesg->est_num_entries)
        UINT16DECODE(p, mesg->est_name_len)
    } else {
        mesg->est_num_entries = OBJ_CRT_GINFO_EST_NUM_ENTRIES;
        mesg->est_name_len = OBJ_CRT_GINFO_EST_NAME_LEN;
    } 

    if(ret_value == SUCCEED)
	ret_ptr = mesg;
done:
    if(ret_value == FAIL && mesg)
	free(mesg);

    return(ret_ptr);
} /* OBJ_ginfo_decode() */

/* Data Storage - Filter Pipeline */
static void *
OBJ_filter_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_filter_t     	*pline=NULL;
    void            	*ret_ptr=NULL;
    int			ret_value=SUCCEED;
    unsigned        	version;
    size_t          	i, j, n, name_length;
    ck_addr_t		logical;
    int			badinfo;

    assert(file);
    assert(p);

    logical = get_logical_addr(p, start_buf, logi_base);

   /* Decode */
   if((pline = calloc(1, sizeof(OBJ_filter_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	  "Filter Pipeline Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)	
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    version = *p++;
    
    /* version is adjusted if error for further continuation */
    if(g_format_num == FORMAT_ONE_SIX) {
	if(version != OBJ_FILTER_VERSION_1) {
	   badinfo = version;
	   version = OBJ_FILTER_VERSION_1;
	   error_push(ERR_LEV_2, ERR_LEV_2A2l, 
		"Filter Pipeline Message:Bad version number", logical, &badinfo);
	   CK_SET_RET(FAIL)
	}
    } else {
	if(version < OBJ_FILTER_VERSION_1 || version > OBJ_FILTER_VERSION_LATEST) {
	   badinfo = version;
	   version = OBJ_FILTER_VERSION_LATEST;
	   error_push(ERR_LEV_2, ERR_LEV_2A2l, "Filter Pipeline Message:Bad version number", logical, &badinfo);
	   CK_SET_RET(FAIL)
	}
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    pline->nused = *p++;
    if((pline->nused > OBJ_MAX_NFILTERS) || (pline->nused <= 0)) {
	badinfo = pline->nused;
	error_push(ERR_LEV_2, ERR_LEV_2A2l, 
	    "Filter Pipeline Message:Invalid # of filters", logical, &badinfo);
	CK_SET_RET(FAIL)
    }

    if(version == OBJ_FILTER_VERSION_1)
        p += 6;  /* reserved */

    pline->nalloc = pline->nused;
    pline->filter = calloc(pline->nalloc, sizeof(pline->filter[0]));
    if(pline->filter == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Filter Pipeline Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)	
    } 

    for(i = 0; i < pline->nused; i++) {
        UINT16DECODE(p, pline->filter[i].id);
	logical = get_logical_addr(p, start_buf, logi_base);

	/* Length of filter name */
        if (version > OBJ_FILTER_VERSION_1 && pline->filter[i].id < OBJ_FILTER_RESERVED)
            name_length = 0;
        else {
            UINT16DECODE(p, name_length);
            if(version == OBJ_FILTER_VERSION_1 && name_length % 8) {
		error_push(ERR_LEV_2, ERR_LEV_2A2l, 
		    "Filter Pipeline Message:Filter name length is not a multiple of eight", logical, NULL);
		CK_SET_RET(FAIL)
	    }
        } /* end if */

        UINT16DECODE(p, pline->filter[i].flags);
        UINT16DECODE(p, pline->filter[i].cd_nelmts);

	logical = get_logical_addr(p, start_buf, logi_base);
        if(name_length) {
	    ck_size_t 	actual_name_length;          /* Actual length of name */

            /* Determine actual name length (without padding, but with null terminator) */
            actual_name_length = strlen((const char *)p) + 1;
	    if (actual_name_length > name_length) {
		error_push(ERR_LEV_2, ERR_LEV_2A2l, "Filter Pipeline Message:Inconsistent name length", logical, NULL);
		CK_SET_RET(FAIL)
	    }

            /* Allocate space for the filter name, or use the internal buffer */
            if(actual_name_length > Z_COMMON_NAME_LEN) {
                pline->filter[i].name = (char *)malloc(actual_name_length);
                if(NULL == pline->filter[i].name) {
		    error_push(ERR_INTERNAL, ERR_NONE_SEC, 
			"Filter Pipeline Message:Internal allocation error", logical, NULL);
		    CK_SET_RET_DONE(FAIL)
		}
            } else
                pline->filter[i].name = pline->filter[i]._name;

            strcpy(pline->filter[i].name, (const char *)p);
            p += name_length;
        }

        if((n=pline->filter[i].cd_nelmts)) {
	    ck_size_t      j;

            if(n > Z_COMMON_CD_VALUES) {
                pline->filter[i].cd_values = (unsigned *)malloc(n*sizeof(unsigned));
                if (NULL == pline->filter[i].cd_values) {
		    error_push(ERR_INTERNAL, ERR_NONE_SEC, 
			"Filter Pipeline Message:Internal allocation error", -1, NULL);
		    CK_SET_RET_DONE(FAIL)
		}
            } else
                pline->filter[i].cd_values = pline->filter[i]._cd_values;

            /*
             * Read the client data values and the padding
             */
            for(j = 0; j < pline->filter[i].cd_nelmts; j++)
                UINT32DECODE(p, pline->filter[i].cd_values[j]);

            if(version == OBJ_FILTER_VERSION_1)
                if(pline->filter[i].cd_nelmts % 2)
                    p += 4; /*padding*/
        }
    }  /* end for */

    if(ret_value == SUCCEED)
	ret_ptr = pline;

done:
    if(ret_value == FAIL && pline) {
	if(pline->filter) {
	    for (i = 0; i < pline->nused; i++) {
		free(pline->filter[i].name);
		free(pline->filter[i].cd_values);
	    }
	    free(pline->filter);
	}
	free(pline);
    }
    return(ret_ptr);
} /* OBJ_filter_decode() */

/* Attribute */
static void *
OBJ_attr_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_attr_t      	*attr=NULL;
    OBJ_sds_extent_t   	*extent;           /* extent dimensionality information  */
    DT_cset_t		encoding;
    size_t          	name_len;          /* attribute name length */
    unsigned        	flags=0;	   /* Attribute flags */
    OBJ_attr_t     	*ret_ptr=NULL; 
    int			ret_value=SUCCEED;
    ck_addr_t		logical;
    int			badinfo, version;

    assert(file);
    assert(p);

    logical = get_logical_addr(p, start_buf, logi_base);
    if((attr = calloc(1, sizeof(OBJ_attr_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Attribute Message:Internal allocation error", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    logical = get_logical_addr(p, start_buf, logi_base);

    version = *p++;

    /* version is adjusted if error for further continuation */
    if(g_format_num == FORMAT_ONE_SIX) {
	if (version != OBJ_ATTR_VERSION_1 && version != OBJ_ATTR_VERSION_2) {
	    badinfo = version;
	    version = OBJ_ATTR_VERSION_2;
	    error_push(ERR_LEV_2, ERR_LEV_2A2m, "Attribute Message:Bad version number", logical, &badinfo);
	    CK_SET_RET(FAIL)
	}
    } else {
	if (version < OBJ_ATTR_VERSION_1 || version > OBJ_ATTR_VERSION_LATEST) {
	    badinfo = version;
	    version = OBJ_ATTR_VERSION_LATEST;
	    error_push(ERR_LEV_2, ERR_LEV_2A2m, "Attribute Message:Bad version number", logical, &badinfo);
	    CK_SET_RET(FAIL)
	}
    }

    logical = get_logical_addr(p, start_buf, logi_base);

    /* version 1 does not support flags; it is reserved and set to zero */
    /* version 2 & 3 used this flag to indicate data type and/or space are shared or not shared */
    if(version >= OBJ_ATTR_VERSION_2) {
	flags = *p++;
	if (flags & (unsigned)~OBJ_ATTR_FLAG_ALL) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2m, "Attribute Message:Unknown flag", logical, NULL);
	    CK_SET_RET(FAIL)
	}
    } else /* byte is unused when version < 2 */
	p++;

    /*
     * Decode the sizes of the parts of the attribute.  The sizes stored in
     * the file are exact but the parts are aligned on 8-byte boundaries.
     */
    logical = get_logical_addr(p, start_buf, logi_base);
    UINT16DECODE(p, name_len); /*including null*/
    UINT16DECODE(p, attr->dt_size);
    UINT16DECODE(p, attr->ds_size);
/* NEED CHECK:to check FOR THOSE 2 fields above to be greater than 0 or >=??? */
/* NEED CHECK:to validate this field ? */
    if(version >= OBJ_ATTR_VERSION_3)
        encoding = *p++;

    /* Decode and store the name */
    logical = get_logical_addr(p, start_buf, logi_base);
    attr->name = NULL;
    if((attr->name = strdup((const char *)p)) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Attribute Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }
    
    /* should be null-terminated */
    if(attr->name[name_len-1] != '\0') {
	error_push(ERR_LEV_2, ERR_LEV_2A2m, 
	    "Attribute Message:Name should be null-terminated", logical, NULL);
	CK_SET_RET(FAIL)
    }

    if(version < OBJ_ATTR_VERSION_2)
        p += OBJ_ALIGN_OLD(name_len);
    else
        p += name_len;

    logical = get_logical_addr(p, start_buf, logi_base);
	
    /* decode the attribute datatype */
    if(flags & OBJ_ATTR_FLAG_TYPE_SHARED) {
	OBJ_shared_t *sh_shared; 
	
        /* Get the shared information */
        if(NULL == (sh_shared = OBJ_shared_decode(file, p, OBJ_DT, start_buf, logi_base))) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2m, 
		"Attribute Message:Errors found when decoding shared datatype", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
	/* Get the actual datatype information */
	if((attr->dt = OBJ_shared_read(file, sh_shared, OBJ_DT))==NULL) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2m, 
		"Attribute Message:Errors found when reading shared datatype", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
    } else {
	if((attr->dt=(OBJ_DT->decode)(file, p, start_buf, logi_base))==NULL) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2m, 
		"Attribute Message:Errors found when decoding datatype description", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
    }

    if(version < OBJ_ATTR_VERSION_2)
        p += OBJ_ALIGN_OLD(attr->dt_size);
    else
        p += attr->dt_size;

    if ((attr->ds = calloc(1, sizeof(OBJ_space_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	  "Attribute Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* decode the attribute datatype */
    logical = get_logical_addr(p, start_buf, logi_base);
    if(flags & OBJ_ATTR_FLAG_SPACE_SHARED) {
	OBJ_shared_t *sh_shared; 
	
        /* Get the shared information */
        if(NULL == (sh_shared = OBJ_shared_decode(file, p, OBJ_SDS, start_buf, logi_base))) {
	   error_push(ERR_LEV_2, ERR_LEV_2A2m, 
		"Attribute Message:Errors found when decoding shared dataspace", 
		logical, NULL);
	   CK_SET_RET_DONE(FAIL)
	}
	/* Get the actual dataspace information */
	if((extent = OBJ_shared_read(file, sh_shared, OBJ_SDS))==NULL) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2m, 
		"Attribute Message:Errors found when reading shared dataspace", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
    } else {
	if((extent = (OBJ_SDS->decode)(file, p, start_buf, logi_base))==NULL) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2m, 
		"Attribute Message:Errors found when decoding dataspace description", 
		logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
    }

    /* Copy the extent information */
    memcpy(&(attr->ds->extent),extent, sizeof(OBJ_sds_extent_t));

    /* Release temporary extent information */ 
    if(extent)
   	free(extent);

    if(version < OBJ_ATTR_VERSION_2)
        p += OBJ_ALIGN_OLD(attr->ds_size);
    else
        p += attr->ds_size;


    /* Compute the size of the data */
    attr->data_size = attr->ds->extent.nelem * attr->dt->shared->size;
    ASSIGN_OVERFLOW(attr->data_size,attr->ds->extent.nelem*attr->dt->shared->size,ck_size_t,size_t);

    /* Go get the data */
    logical = get_logical_addr(p, start_buf, logi_base);
    if(attr->data_size) {
	if((attr->data = malloc(attr->data_size))==NULL) {
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, 
		"Attribute Message:Internal allocation error", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
	/* How the data is interpreted is not stated in the format specification. */
	memcpy(attr->data, p, attr->data_size);
    }

    if(ret_value == SUCCEED)
	ret_ptr = attr;
done:
    if (ret_value == FAIL && attr) {
	if(attr->name) free(attr->name);
	if(attr->ds) free(attr->ds);
	if(attr->data) free(attr->data);
	free(attr);
    }

    return(ret_ptr);
} /* OBJ_attr_decode() */


/* Object Comment */
static void *
OBJ_comm_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_comm_t  *mesg=NULL;
    void        *ret_ptr=NULL;
    int		ret_value=SUCCEED;
    int		len;
    ck_addr_t	logical;


    /* check args */
    assert(file);
    assert(p);

    len = strlen((const char *)p);

    if((mesg=calloc(1, sizeof(OBJ_comm_t))) == NULL || (mesg->s = malloc(len+1))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	  "Could not calloc() OBJ_comm_t for Object Comment Message", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    strcpy(mesg->s, (const char*)p);
    if(mesg->s[len] != '\0') {
	error_push(ERR_LEV_2, ERR_LEV_2A2n, 
	  "Object Comment Message:The comment string should be null-terminated", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    if(ret_value == SUCCEED)
	ret_ptr = mesg;

done:
    if(ret_value == FAIL && mesg) {
	if (mesg->s) free(mesg->s);
	free(mesg);
    }
    return(ret_ptr);
} /* OBJ_comm_decode() */

/* Object Modification Date & Time (OLD) */
static void *
OBJ_mdt_old_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    time_t      *mesg=NULL;
    time_t	the_time;
    int 	i;
    struct tm   tm;
    int		ret_value=SUCCEED;
    void        *ret_ptr=NULL;
    ck_addr_t	logical;

    assert(file);
    assert(p);

    /* Initialize time zone information */
    if(!ntzset) {
        tzset();
        ntzset = 1;
    } 

    /* decode */
    logical = get_logical_addr(p, start_buf, logi_base);
    for(i = 0; i < 14; i++) {
	if(!isdigit(p[i])) {
	    error_push(ERR_LEV_2, ERR_LEV_2A2o, 
		"Object Modification Time (old) Message:Badly formatted time", logical, NULL);
	    CK_SET_RET(FAIL)
	}
    }

    memset(&tm, 0, sizeof tm);
    tm.tm_year = (p[0]-'0')*1000 + (p[1]-'0')*100 +
                 (p[2]-'0')*10 + (p[3]-'0') - 1900;
    tm.tm_mon = (p[4]-'0')*10 + (p[5]-'0') - 1;
    tm.tm_mday = (p[6]-'0')*10 + (p[7]-'0');
    tm.tm_hour = (p[8]-'0')*10 + (p[9]-'0');
    tm.tm_min = (p[10]-'0')*10 + (p[11]-'0');
    tm.tm_sec = (p[12]-'0')*10 + (p[13]-'0');
    tm.tm_isdst = -1; /*figure it out*/
    if ((time_t)-1 == (the_time=mktime(&tm))) {
	error_push(ERR_LEV_2, ERR_LEV_2A2o, 
	    "Object Modification Time (old) Message:Badly formatted time", logical, NULL);
	CK_SET_RET(FAIL)
    }

    if (NULL == (mesg = calloc(1, sizeof(time_t)))) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	  "Object Modification Time (old) Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }
    *mesg = the_time;

    if(ret_value == SUCCEED)
	ret_ptr = mesg;

done:
    if(ret_value == FAIL && mesg)
	free(mesg);

    return(ret_ptr);
} /* OBJ_mdt_old_decode() */


/* Shared message table Message */
static void *
OBJ_shmesg_decode(driver_t *file, const uint8_t *buf, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_shmesg_table_t  *mesg=NULL; 
    void                *ret_ptr=NULL;
    int			ret_value=SUCCEED;
    ck_addr_t		logical;
    int			badinfo;

    /* check arguments */
    assert(file);
    assert(buf);

    if(g_format_num == FORMAT_ONE_SIX) {
	error_push(ERR_LEV_2, ERR_LEV_2A2p, "Shared Message Table Message:Unsupported message", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    if ((mesg = calloc(1, sizeof(OBJ_shmesg_table_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Shared Message Table Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Retrieve version, table address, and number of indexes */
    logical = get_logical_addr(buf, start_buf, logi_base);
    mesg->version = *buf++;
    if (mesg->version != SHAREDHEADER_VERSION) {
	badinfo = mesg->version;
	error_push(ERR_LEV_2, ERR_LEV_2A2p, "Shared Message Table Message:Bad version number", 
	    logical, &badinfo);
	CK_SET_RET(FAIL)
    }

    logical = get_logical_addr(buf, start_buf, logi_base);
    addr_decode(file->shared, &buf, &(mesg->addr));
    if(mesg->addr == CK_ADDR_UNDEF) {
	error_push(ERR_LEV_2, ERR_LEV_2A2p, "Shared Message Table Message:Undefined address",
	    logical, NULL);
	CK_SET_RET(FAIL)
    }
    
    logical = get_logical_addr(buf, start_buf, logi_base);
    mesg->nindexes = *buf++;
    /* nindexes < 256--1 byte to hold nindexes */
    if((mesg->nindexes <= 0) && (mesg->nindexes > OBJ_SHMESG_MAX_NINDEXES)) {
	badinfo = mesg->nindexes;
	error_push(ERR_LEV_2, ERR_LEV_2A2p, "Shared Message Table Message:Invalid value for number of indices", 
	    logical, NULL);
	CK_SET_RET(FAIL)
    }

    if(ret_value == SUCCEED)
	ret_ptr = (void *)mesg;

done:
    if(ret_value == FAIL && mesg) 
	free(mesg);
    return(ret_ptr);
} /* OBJ_shmesg_decode() */

/* Object Header Continutaion */
static void *
OBJ_cont_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_cont_t 	*cont = NULL;
    void       	*ret_ptr=NULL;
    int		ret_value=SUCCEED;
    ck_addr_t	logical;

    assert(file);
    assert(p);


    cont = calloc(1, sizeof(OBJ_cont_t));
    if(cont == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	  "Object Header Continuation Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }
	
    logical = get_logical_addr(p, start_buf, logi_base);
    addr_decode(file->shared, &p, &(cont->addr));

    if(cont->addr == CK_ADDR_UNDEF) {
	error_push(ERR_LEV_2, ERR_LEV_2A2p, 
	  "Object Header Continuation Message:Undefined offset", logical, NULL);
	CK_SET_RET(FAIL)
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    DECODE_LENGTH(file->shared, p, cont->size);

    if(cont->size < 0) {
	error_push( ERR_LEV_2, ERR_LEV_2A2p, 
	  "Object Header continuation Message:Invalid length", logical, NULL);
	CK_SET_RET(FAIL)
    }

    cont->chunkno = 0;

    if(ret_value == SUCCEED)
	ret_ptr = cont;

done:
    if(ret_value == FAIL && cont)
	free(cont);
    return(ret_ptr);
} /* OBJ_cont_decode() */


/* Symbol Table Message */
static void *
OBJ_group_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_stab_t  	*stab=NULL;
    void    		*ret_ptr=NULL;
    int			ret_value=SUCCEED;
    ck_addr_t		logical;

    assert(file);
    assert(p);


    stab = calloc(1, sizeof(OBJ_stab_t));
    if(stab == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	  "Symbol Table Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    addr_decode(file->shared, &p, &(stab->btree_addr));
    if(stab->btree_addr == CK_ADDR_UNDEF) {
	error_push(ERR_LEV_2, ERR_LEV_2A2r, "Symbol Table Message:Undefined version 1 btree address", logical, NULL);
	CK_SET_RET(FAIL)
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    addr_decode(file->shared, &p, &(stab->heap_addr));
    if(stab->heap_addr == CK_ADDR_UNDEF) {
	error_push(ERR_LEV_2, ERR_LEV_2A2r, "Symbol Table Message:Undefined local heap address", logical, NULL);
	CK_SET_RET(FAIL)
    }

    if(ret_value == SUCCEED)
	ret_ptr = stab;

done:
    if(ret_value == FAIL && stab)
	free(stab);

    return(ret_ptr);
} /* OBJ_group_decode() */


/* Object Modification Time */
static void *
OBJ_mdt_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    time_t      *mesg=NULL;
    uint32_t    tmp_time;    
    void        *ret_ptr=NULL;
    int		ret_value=SUCCEED;
    ck_addr_t	logical;
    int		version, badinfo;

    assert(file);
    assert(p);

    logical = get_logical_addr(p, start_buf, logi_base);
    version = *p++;
    if(version != OBJ_MTIME_VERSION) {
	badinfo = version;
	error_push(ERR_LEV_2, ERR_LEV_2A2s, 
	  "Object Modification Time Message:Bad version number", logical, &badinfo);
	CK_SET_RET(FAIL)
    }

    /* Skip reserved bytes */
    p+=3;

    /* Get the time_t from the file */
    logical = get_logical_addr(p, start_buf, logi_base);
    UINT32DECODE(p, tmp_time);

    if((mesg=calloc(1, sizeof(time_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	  "Object Modification Time Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }
    *mesg = (time_t)tmp_time;

    if(ret_value == SUCCEED)
	ret_ptr = mesg;

done:
    if(ret_value == FAIL && mesg)
	free(mesg);

    return(ret_ptr);
} /* OBJ_mdt_decode() */

/* Non-default v1 B-tree 'K' values message */
static void *
OBJ_btreek_decode(driver_t *file, const uint8_t *buf, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_btreek_t  	*mesg=NULL;   
    void          	*ret_ptr=NULL;   
    int			ret_value=SUCCEED;
    ck_addr_t		logical;
    unsigned		version;
    int			badinfo;

    /* check arguments */
    assert(file);
    assert(buf);
	
    if(g_format_num == FORMAT_ONE_SIX) {
	error_push(ERR_LEV_2, ERR_LEV_2A2t, "B-tree 'K' Values Message:Unsupported message", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    if ((mesg = calloc(1, sizeof(OBJ_btreek_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "B-tree 'K' Values Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    logical = get_logical_addr(buf, start_buf, logi_base);
    version = *buf++;
    if(version != OBJ_BTREEK_VERSION) {
	badinfo = version;
	error_push(ERR_LEV_2, ERR_LEV_2A2t, "B-tree 'K' Values Message:Bad Version number", logical, &badinfo);
	CK_SET_RET(FAIL)
    }

    logical = get_logical_addr(buf, start_buf, logi_base);
    UINT16DECODE(buf, mesg->btree_k[BT_ISTORE_ID]);
    if(mesg->btree_k[BT_ISTORE_ID] <= 0) {
	error_push(ERR_LEV_2, ERR_LEV_2A2t, 
	    "B-tree 'K' Values Message:Invalid value for Indexed Storage Internal Node K", logical, NULL);
	CK_SET_RET(FAIL)
    }

    logical = get_logical_addr(buf, start_buf, logi_base);
    UINT16DECODE(buf, mesg->btree_k[BT_SNODE_ID]);
    if(mesg->btree_k[BT_SNODE_ID] <= 0) {
	error_push(ERR_LEV_2, ERR_LEV_2A2t, "B-tree 'K' Values Message:Invalid value for Group Internal Node K", logical, NULL);
	    CK_SET_RET(FAIL)
    }

    logical = get_logical_addr(buf, start_buf, logi_base);
    UINT16DECODE(buf, mesg->sym_leaf_k);
    if(mesg->sym_leaf_k <= 0) {
	error_push(ERR_LEV_2, ERR_LEV_2A2t, "B-tree 'K' Values Message:Invalid value for Group Leaf Node K", logical, NULL);
	    CK_SET_RET(FAIL)
    }

    if(ret_value == SUCCEED)
	ret_ptr = (void *)mesg;

done:
    if(ret_value == FAIL && mesg)
	free(mesg);
    return(ret_ptr);
} /* OBJ_btreek_decode() */

/* Driver Info Message */
static void *
OBJ_drvinfo_decode(driver_t *file, const uint8_t *buf, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_drvinfo_t  	*mesg=NULL;
    void          	*ret_ptr=NULL;
    int			ret_value=SUCCEED;
    ck_addr_t		logical;
    int			badinfo;
    unsigned		version;

    /* check arguments */
    assert(file);
    assert(buf);

    if(g_format_num == FORMAT_ONE_SIX) {
	error_push(ERR_LEV_2, ERR_LEV_2A2u, "Driver Info Message:Unsupported message", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    if((mesg = calloc(1, sizeof(OBJ_drvinfo_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Driver Info Message: Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

     /* Version of message */
    logical = get_logical_addr(buf, start_buf, logi_base);
    version = *buf++;
    if(version != OBJ_DRVINFO_VERSION) {
	badinfo = version;
	error_push(ERR_LEV_2, ERR_LEV_2A2u, "Driver Info Message: Bad version number", logical, &badinfo);
	CK_SET_RET(FAIL)
    }

    /* Retrieve driver name */
    memcpy(mesg->name, buf, 8);
    mesg->name[8] = '\0';
    buf += 8;

    /* Decode Driver Information Size */
    logical = get_logical_addr(buf, start_buf, logi_base);
    UINT16DECODE(buf, mesg->len);
    if(mesg->len <=0) {
	error_push(ERR_LEV_2, ERR_LEV_2A2u, "Driver Info Message:Invalid driver information size", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Allocate space for buffer */
    if(NULL == (mesg->buf = malloc(mesg->len))) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Driver Info Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Copy encoded driver info into buffer */
    memcpy(mesg->buf, buf, mesg->len);

    if(ret_value == SUCCEED)
	ret_ptr = (void *)mesg;

done:
    if(ret_value == FAIL && mesg)
	free(mesg);

    return(ret_ptr);
} /* OBJ_drvinfo_decode() */

/* Attribute Info Message */
static void *
OBJ_ainfo_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_ainfo_t 	*ainfo=NULL;  /* Attribute info */
    unsigned char 	flags;        /* Flags for encoding attribute info */
    unsigned		version;
    int			badinfo;
    void        	*ret_ptr=NULL;
    int			ret_value=SUCCEED;
    ck_addr_t		logical;

    /* check args */
    assert(file);
    assert(p);

    if(g_format_num == FORMAT_ONE_SIX) {
	error_push(ERR_LEV_2, ERR_LEV_2A2v, "Attribute Info Message:Unsupported message", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Allocate space for message */
    if(NULL == (ainfo = calloc(1, sizeof(OBJ_ainfo_t)))) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Attribute Info Message: Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    logical = get_logical_addr(p, start_buf, logi_base);
    version = *p++;
    if(version != OBJ_AINFO_VERSION) {
	badinfo = version;
	error_push(ERR_LEV_2, ERR_LEV_2A2v, "Attribute Info Message: Bad version number", logical, &badinfo);
	CK_SET_RET(FAIL)
    }


    /* Get the flags for the message */
    logical = get_logical_addr(p, start_buf, logi_base);
    flags = *p++;
    if(flags & ~OBJ_AINFO_ALL_FLAGS) {
	error_push(ERR_LEV_2, ERR_LEV_2A2v, "Attribute Info Message: Bad flag value", logical, &badinfo);
	CK_SET_RET(FAIL)
    }

    ainfo->track_corder = (flags & OBJ_AINFO_TRACK_CORDER) ? TRUE : FALSE;
    ainfo->index_corder = (flags & OBJ_AINFO_INDEX_CORDER) ? TRUE : FALSE;

    /* Max. creation order value for the object */
    if(ainfo->track_corder)
        UINT16DECODE(p, ainfo->max_crt_idx)
    else
        ainfo->max_crt_idx = OBJ_MAX_CRT_ORDER_IDX;

    /* Address of fractal heap to store "dense" attributes */
    addr_decode(file->shared, &p, &(ainfo->fheap_addr));

    /* Address of v2 B-tree to index names of attributes (names are always indexed) */
    addr_decode(file->shared, &p, &(ainfo->name_bt2_addr));

    /* Address of v2 B-tree to index creation order of links, if there is one */
    if(ainfo->index_corder)
        addr_decode(file->shared, &p, &(ainfo->corder_bt2_addr));
    else
        ainfo->corder_bt2_addr = CK_ADDR_UNDEF;

    if (ret_value == SUCCEED)
	ret_ptr = ainfo;

done:
    if(ret_value == FAIL && ainfo)
        free(ainfo);

    return(ret_ptr);
} /* OBJ_ainfo_decode() */

/* Object Reference Count Message */
static void *
OBJ_refcount_decode(driver_t *file, const uint8_t *p, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_refcount_t 	*refcount=NULL;  /* Reference count */
    void        	*ret_ptr =NULL;
    int			badinfo;
    unsigned		version;
    int			ret_value=SUCCEED;
    ck_addr_t		logical;

    /* check args */
    assert(file);
    assert(p);

    if(g_format_num == FORMAT_ONE_SIX) {
	error_push(ERR_LEV_2, ERR_LEV_2A2w, "Object Reference Count Message:Unsupported message", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Allocate space for message */
    if(NULL == (refcount = malloc(sizeof(OBJ_refcount_t)))) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Object Reference Count Message:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Version of message */
    logical = get_logical_addr(p, start_buf, logi_base);
    version = *p++;
    if(version != OBJ_REFCOUNT_VERSION) {
	badinfo = version;
	error_push(ERR_LEV_2, ERR_LEV_2A2w, "Object Reference Count Message: Bad version number", logical, &badinfo);
	CK_SET_RET(FAIL)
    }

    /* Get ref. count for object */
    UINT32DECODE(p, *refcount)

    if(ret_value == SUCCEED)
	ret_ptr = refcount;

done:
    if(ret_value == FAIL && refcount)
        free(refcount);

    return(ret_ptr);
} /* OBJ_refcount_decode() */

/*
 *  end decoder for header messages
 */


static type_t *
type_alloc(ck_addr_t logical)
{
    type_t *dt=NULL;                  /* Pointer to datatype allocated */
    type_t *ret_value=NULL;           /* Return value */


    /* Allocate & initialize new datatype info */
    if((dt=calloc(1, sizeof(type_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Could not calloc() type_t in type_alloc()", logical, NULL);
	CK_SET_RET_DONE(NULL)
    }

    assert(&(dt->ent));
    memset(&(dt->ent), 0, sizeof(GP_entry_t));
    dt->ent.header = CK_ADDR_UNDEF;

    if((dt->shared=calloc(1, sizeof(DT_shared_t)))==NULL) {
	if(dt != NULL) free(dt);
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, 
		"Could not calloc() H5T_shared_t in type_alloc()", logical, NULL);
	    CK_SET_RET_DONE(NULL)
    }
   
    /* Assign return value */
    ret_value = dt;

done:
    if(ret_value == NULL && dt) {
	if (dt->shared) free(dt->shared);
	free(dt);
    }
    return(ret_value);
} /* type_alloc() */


/* size of a symbol table node */
static size_t
gp_node_size(global_shared_t *shared)
{
    return (SNODE_SIZEOF_HDR(shared) +
	    (2 * SYM_LEAF_K(shared)) * GP_SIZEOF_ENTRY(shared));
} /* gp_node_size */

/*
 * Decode a symbol table group entry
 */
static ck_err_t
gp_ent_decode(global_shared_t *shared, const uint8_t **pp, GP_entry_t *ent)
{
    const uint8_t *p_ret = *pp;
    uint32_t      tmp;
    ck_err_t    ret_value=SUCCEED;       /* Return value */

    assert(pp);
    assert(ent);

    /* decode header */
    DECODE_LENGTH(shared, *pp, ent->name_off);
    addr_decode(shared, pp, &(ent->header));
    UINT32DECODE(*pp, tmp);
    *pp += 4; /*reserved*/
    ent->type=(GP_type_t)tmp;

    /* decode scratch-pad */
    switch (ent->type) {
	case GP_NOTHING_CACHED:
            break;

        case GP_CACHED_STAB:
            assert(2*SIZEOF_ADDR(shared) <= GP_SIZEOF_SCRATCH);
            addr_decode(shared, pp, &(ent->cache.stab.btree_addr));
            addr_decode(shared, pp, &(ent->cache.stab.heap_addr));
            break;

        case GP_CACHED_SLINK:
            UINT32DECODE (*pp, ent->cache.slink.lval_offset);
            break;

        default:  /* allows other cache values; will be validated later */
	    CK_SET_RET_DONE(FAIL);
    }

    *pp = p_ret + GP_SIZEOF_ENTRY(shared);

done:
    return(ret_value);

} /* gp_ent_decode() */


/*
 * Decode a vector of symbol table group entries 
 */
static ck_err_t
gp_ent_decode_vec(global_shared_t *shared, const uint8_t **pp, GP_entry_t *ent, unsigned n)
{
    unsigned    u;
    ck_err_t    ret_value=SUCCEED;       /* Return value */

    /* check arguments */
    assert(pp);
    assert(ent);

    for (u = 0; u < n; u++)
	if ((ret_value = gp_ent_decode(shared, pp, ent + u)) < SUCCEED) {
	    error_push(ERR_LEV_1, ERR_LEV_1C, 
		"Symbol table node:Unable to decode node entries", -1, NULL);
            CK_SET_RET_DONE(FAIL);
	}
done:
    return(ret_value);
} /* gp_ent_decode_vec() */

/* support routine for decoding an address */
void
addr_decode(global_shared_t *shared, const uint8_t **pp/*in,out*/, ck_addr_t *addr_p/*out*/)
{
    unsigned    	i;
    ck_addr_t          	tmp;
    uint8_t            	c;
    ck_bool_t          	all_zero = TRUE;
    ck_addr_t		ret;

    assert(pp && *pp);
    assert(addr_p);

    *addr_p = 0;

    for (i=0; i<SIZEOF_ADDR(shared); i++) {
	c = *(*pp)++;
	if (c != 0xff)
	    all_zero = FALSE;

       	if (i<sizeof(*addr_p)) {
	    tmp = c;
	    tmp <<= (i * 8);    /*use tmp to get casting right */
	    *addr_p |= tmp;
        } else if (!all_zero) {
	    if (**pp != 0) {
		*addr_p = CK_ADDR_UNDEF;
		break;
	    }
	    assert(0 == **pp);  /*overflow */
        }
    }
    if (all_zero)
	*addr_p = CK_ADDR_UNDEF;
} /* addr_decode() */


/* size of the key for a group node (symbol table node) */
static size_t
gp_node_sizeof_rkey(global_shared_t *shared, unsigned UNUSED ndims)
{

    return (SIZEOF_SIZE(shared));       /*the name offset */
} /* gp_node_sizeof_rkey() */

/* decode the key for group node (symbol table node) */
static ck_err_t
gp_node_decode_key(global_shared_t *shared, unsigned UNUSED ndims, const uint8_t **p/*in,out*/, void **_key)
{
    GP_node_key_t *key;

    assert(p);
    assert(*p);

    /* decode */
    if((key = calloc(1, sizeof(GP_node_key_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Key for group node in version 1 B-tree:Internal allocation error", -1, NULL);
	return(FAIL);
    }

    DECODE_LENGTH(shared, *p, key->offset);

    *_key = key;
    return(SUCCEED);
} /* gp_node_decode_key() */

static int
gp_node_cmp_key(global_shared_t *shared, unsigned UNUSED ndims, uint8_t *heap_chunk, void *_lt_key, void *_rt_key)
{
    GP_node_key_t *lt_key = (GP_node_key_t *)_lt_key;
    GP_node_key_t *rt_key = (GP_node_key_t *)_rt_key;
    const char *s1 = NULL, *s2 = NULL;
    int ret_value;

    assert(lt_key);
    assert(rt_key);
    assert(heap_chunk);

    if(heap_chunk) {
	s1 = (char *)(heap_chunk + HL_SIZEOF_HDR(shared) + lt_key->offset);
	s2 = (char *)(heap_chunk + HL_SIZEOF_HDR(shared) + rt_key->offset);

	/* Set return value */
	ret_value = strcmp(s1, s2);
    }

    return(ret_value);

} /* gp_node_cmp_key() */



/* size of the key for a chunked raw data node (indexed storage node) */
static size_t
raw_node_sizeof_rkey(global_shared_t *shared, unsigned ndims)
{
	
    size_t 	nbytes;

    assert(shared);
    assert(ndims>0 && ndims<=OBJ_LAYOUT_NDIMS);

    nbytes = 4 +           /*storage size          */
	     4 +           /*filter mask           */
             ndims*8;      /*dimension indices     */

    return(nbytes);
} /* end raw_node_sizeof_rkey() */

/* decode the key for chunked raw data node */
static ck_err_t
raw_node_decode_key(global_shared_t *shared, unsigned ndims, const uint8_t **p/*in, out*/, void **_key)
{
    RAW_node_key_t *key;
    unsigned  u;

    /* check args */
    assert(shared);
    assert(p);
    assert(*p);
    assert(ndims>0);
    assert(ndims<=OBJ_LAYOUT_NDIMS);

    /* decode */
    if((key = calloc(1, sizeof(RAW_node_key_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Key for chunked node in version 1 B-tree:Internal allocation error", -1, NULL);
	return(FAIL);
    }
    
    UINT32DECODE(*p, key->nbytes);
    UINT32DECODE(*p, key->filter_mask);
    for (u = 0; u < ndims; u++)
	UINT64DECODE(*p, key->offset[u]);

    *_key = key;
    return(SUCCEED);
} /* end of raw_node_decode_key() */

static int
raw_node_cmp_key(global_shared_t *shared, unsigned ndims, uint8_t *UNUSED heap_chunk, void *_lt_key, void *_rt_key)
{
    RAW_node_key_t *lt_key = (RAW_node_key_t *)_lt_key;
    RAW_node_key_t *rt_key = (RAW_node_key_t *)_rt_key;
    int ret_value;

    assert(lt_key);
    assert(rt_key);
    assert(ndims > 0);
    assert(ndims <= OBJ_LAYOUT_NDIMS);

    /* Compare the offsets but ignore the other fields */
    ret_value = vector_cmp(ndims, lt_key->offset, rt_key->offset);

    return(ret_value);
} /* raw_node_cmp_key() */

/*
 * Validate superblock for both version 0, 1 & 2
 */
static ck_addr_t
locate_super_signature(driver_t *file)
{
    ck_addr_t       addr = 0;
    uint8_t         buf[HDF_SIGNATURE_LEN];
    unsigned        n, maxpow;

    /* find file size */
    addr = FD_get_eof(file);
	
    /* Find the smallest N such that 2^N is larger than the file size */
    for (maxpow=0; addr; maxpow++)
	addr>>=1;
    maxpow = MAX(maxpow, 9);

    /*
     * Search for the file signature at format address zero followed by
     * powers of two larger than 9.
     */
    for (n=8; n<maxpow; n++) {
       	addr = (8==n) ? 0 : (ck_addr_t)1 << n;
	if (FD_read(file, addr, (size_t)HDF_SIGNATURE_LEN, buf) == FAIL) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, 
		"Superblock:Errors when reading superblock signature", LOGI_SUPER_BASE, NULL);
	    addr = CK_ADDR_UNDEF;
	    break;
	}
	if (memcmp(buf, HDF_SIGNATURE, (size_t)HDF_SIGNATURE_LEN) == 0) {
	    if (debug_verbose())
		printf("FOUND super block signature\n");
	    break;
	}
    }
    if (n >= maxpow) {
	error_push(ERR_LEV_0, ERR_LEV_0A, 
	    "Superblock:Unable to find super block signature", LOGI_SUPER_BASE, NULL);
	addr = CK_ADDR_UNDEF;
    }

    /* Set return value */
    return(addr);
} /* end locate_super_signature() */

ck_err_t
check_superblock(driver_t *file)
{
    uint8_t		*p, *start_buf;
    uint8_t		buf[MAX_SUPERBLOCK_SIZE];
    size_t		fixed_size = SUPERBLOCK_FIXED_SIZE;
    size_t		variable_size, remain_size;

    unsigned		super_vers; 
    unsigned		freespace_vers;     /* version # of global free-space storage */
    unsigned		root_sym_vers;      /* version # of root group symbol table entry */
    unsigned		shared_head_vers;   /* version # of shared header message format */

    GP_entry_t 		*root_ent;
    global_shared_t	*lshared;
    ck_addr_t		logical;	/* logical address */
    ck_err_t		ret_value = SUCCEED;

    char 		driver_name[9];
    int			badinfo;

    /* 
     * The super block may begin at 0, 512, 1024, 2048 etc. 
     * or may not be there at all.
     */
    lshared = file->shared;
    lshared->super_addr = 0;
    lshared->super_addr = locate_super_signature(file);
    if(lshared->super_addr == CK_ADDR_UNDEF) {
	if (!object_api()) {
	    error_print(stderr, file);
	    error_clear();
	}
	printf("ASSUMING super block at physical address 0.\n");
	lshared->super_addr = 0;
    }
	
    if(debug_verbose())
	printf("VALIDATING the super block at physical address %llu...\n", 
	    lshared->super_addr);

    if(FD_read(file, LOGI_SUPER_BASE, fixed_size, buf) == FAIL) {
	error_push(ERR_FILE, ERR_NONE_SEC, 
	    "Superblock:Unable to read in the fixed size portion of the superblock", 
	    LOGI_SUPER_BASE, NULL);
	CK_SET_RET_DONE(FAIL)
    } 
    /* superblock signature already checked */
    p = buf + HDF_SIGNATURE_LEN;
    start_buf = buf;

    logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
    super_vers = *p++;

    /* if super_vers is wrong, it would be adjusted accordingly for validation to continue */
    if(g_format_num == FORMAT_ONE_SIX) {
	if((super_vers != SUPERBLOCK_VERSION_0) && (super_vers != SUPERBLOCK_VERSION_1)) {
	    badinfo = super_vers;
	    super_vers = SUPERBLOCK_VERSION_1;
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock:Version number should be 0 or 1", logical, &badinfo);
	    CK_SET_RET(FAIL)
	}
    } else if(g_format_num == DEFAULT_FORMAT) {
	if (super_vers > SUPERBLOCK_VERSION_LATEST) {
	    badinfo = super_vers;
	    super_vers = SUPERBLOCK_VERSION_LATEST;
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock:Version number should be 0, 1 or 2", logical, &badinfo);
	    CK_SET_RET(FAIL)
	}
    } else
	printf("Invalid library version...shouldn't happen\n");

    logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);

    /* Determine the size of the variable-length part of the superblock */
    variable_size = SUPERBLOCK_VARLEN_SIZE(super_vers);

    if((fixed_size + variable_size) > sizeof(buf)) {
	error_push(ERR_LEV_0, ERR_LEV_0A, 
	    "Superblock:Total size of super block is incorrect", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }
    if(FD_read(file, LOGI_SUPER_BASE+fixed_size, variable_size, p) == FAIL) {
	error_push(ERR_FILE, ERR_NONE_SEC, 
	    "Superblock:Unable to read in the variable size portion of the superblock", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* set default driver to SEC2 unless changed later */
    driver_name[0] = '\0';
    set_driver_id(&(lshared->driverid)/*out*/, driver_name);

    /* Check for older version of superblock format */
    if(super_vers == SUPERBLOCK_VERSION_0 || super_vers == SUPERBLOCK_VERSION_1) {
	if(debug_verbose())
	    printf("Validating version 0/1 superblock...\n");
     	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	freespace_vers = *p++;

	if(freespace_vers != FREESPACE_VERSION) {
	    badinfo = freespace_vers;
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.0/1:Version number of Global Free-space Storage should be 0", 
		logical, &badinfo);
	    CK_SET_RET(FAIL)
	}

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	root_sym_vers = *p++;
	if(root_sym_vers != OBJECTDIR_VERSION) {
	    badinfo = root_sym_vers;
	    error_push(ERR_LEV_0, ERR_LEV_0A, 
		"Superblock 0/1:Version number of the Root Group Symbol Table Entry should be 0", 
		logical, &badinfo);
	    CK_SET_RET(FAIL)
	}

	/* skip over reserved byte */
	p++;

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	shared_head_vers = *p++;
	if(shared_head_vers != SHAREDHEADER_VERSION) {
	    badinfo = shared_head_vers;
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.0/1:Version number of Shared Header Message Format should be 0", 
		logical, &badinfo);
	    CK_SET_RET(FAIL)
	}

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	lshared->size_offsets = *p++;
	if(lshared->size_offsets != 2 && lshared->size_offsets != 4 &&
            lshared->size_offsets != 8 && lshared->size_offsets != 16 && 
	    lshared->size_offsets != 32) {
		error_push(ERR_LEV_0, ERR_LEV_0A, 
		    "Superblock v.0/1:Invalid Size of Offsets", logical, NULL);
		CK_SET_RET(FAIL)
	}

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	lshared->size_lengths = *p++;
	if(lshared->size_lengths != 2 && lshared->size_lengths != 4 &&
            lshared->size_lengths != 8 && lshared->size_lengths != 16 && 
	    lshared->size_lengths != 32) {
		error_push(ERR_LEV_0, ERR_LEV_0A, 
		    "Superblock v.0/1:Invalid Size of Lengths", logical, NULL);
		CK_SET_RET(FAIL)
	}
	/* skip over reserved byte */
	p++;

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	UINT16DECODE(p, lshared->gr_leaf_node_k);
	if (lshared->gr_leaf_node_k <= 0) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, 
		"Superblock v.0/1:Invalid value for Group Leaf Node K", logical, NULL);
	    CK_SET_RET(FAIL)
	}

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	UINT16DECODE(p, lshared->btree_k[BT_SNODE_ID]);
	if(lshared->btree_k[BT_SNODE_ID] <= 0) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, 
		"Superblock v.0/1:Invalid value for Group Internal Node K", logical, NULL);
	    CK_SET_RET(FAIL)
	}

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
    	UINT32DECODE(p, lshared->file_consist_flg);
	if(lshared->file_consist_flg > 255) { /* 2 to the power 8 */
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.0/1:Invalid value for file consistency flags.", logical, NULL);
	    CK_SET_RET(FAIL)
	}
	if(lshared->file_consist_flg & ~SUPER_ALL_FLAGS) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.0/1:Invalid file consistency flags.", logical, NULL);
	    CK_SET_RET(FAIL)
	}

    	/*
     	 * If the superblock version # is greater than 0, read in the indexed
     	 * storage B-tree internal 'K' value
     	 */
    	if(super_vers > SUPERBLOCK_VERSION_0) {
	    logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	    UINT16DECODE(p, lshared->btree_k[BT_ISTORE_ID]);
	    p += 2;   /* reserved */
	    if (lshared->btree_k[BT_ISTORE_ID] == 0) {
		error_push(ERR_LEV_0, ERR_LEV_0A, 
		    "Superblock v.1:Invalid value for Indexed Storage Internal Node K", logical, NULL);
		CK_SET_RET(FAIL)
	    }
	} else
	    lshared->btree_k[BT_ISTORE_ID] = BT_ISTORE_K;


	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);

	/* Determine the remaining size of the superblock for V0 and V1 */
	remain_size = SUPERBLOCK_REMAIN_SIZE(super_vers, lshared);
	if((fixed_size+variable_size+remain_size) > sizeof(buf)) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, 
		"Superblock v0/1:Total size of super block is incorrect", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	if(FD_read(file, LOGI_SUPER_BASE+fixed_size+variable_size, remain_size, p) == FAIL) {
	    error_push(ERR_FILE, ERR_NONE_SEC, 
		"Superblock v.0/1:Unable to read in the remaining size portion of the superblock", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}

	addr_decode(lshared, (const uint8_t **)&p, &lshared->base_addr);
	if(lshared->base_addr != lshared->super_addr) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.0/1:Invalid base address", 
		logical, NULL);
	    CK_SET_RET(FAIL)
	} 

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
    	addr_decode(lshared, (const uint8_t **)&p, &lshared->extension_addr);
	if(lshared->extension_addr != CK_ADDR_UNDEF) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.0/1:Address of global Free-space Index should be undefined", 
		logical, NULL);
	    CK_SET_RET(FAIL)
	}

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
    	addr_decode(lshared, (const uint8_t **)&p, &lshared->stored_eoa);
	if((lshared->stored_eoa == CK_ADDR_UNDEF) || (lshared->base_addr >= lshared->stored_eoa)) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.0/1:Invalid End of File Address", logical, NULL);
	    CK_SET_RET(FAIL)
	}

    	addr_decode(lshared, (const uint8_t **)&p, &lshared->driver_addr);

	root_ent = calloc(1, sizeof(GP_entry_t));
	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	if(gp_ent_decode(lshared, (const uint8_t **)&p, root_ent/*out*/) < SUCCEED) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.0/1:Unable to read root symbol table entry", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
	lshared->root_grp = root_ent;

#ifdef DEBUG
	printf("header_addr=%llu, stored_eoa=%llu, base_addr=%llu, super_addr=%llu\n",
		lshared->root_grp->header,lshared->stored_eoa, lshared->base_addr, lshared->super_addr);
#endif

	if(lshared->root_grp->header == CK_ADDR_UNDEF) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, 
		"Superblock v.0/1:Undefined object header address in root group symbol table entry", logical, NULL);
	    CK_SET_RET(FAIL)
	}

	driver_name[0] = '\0';

	/* Read in driver information block and validate */
	/* Defined driver information block address or not */
    	if(lshared->driver_addr != CK_ADDR_UNDEF) {
	    /* since base_addr may not be valid, use super_addr to be sure */
	    ck_addr_t drv_addr = lshared->super_addr + lshared->driver_addr;
	    uint8_t 	dbuf[DRVINFOBLOCK_SIZE];     /* Local buffer */
	    uint8_t	*drv_start;
	    size_t  	driver_size;   		     /* Size of driver info block, in bytes */
	    int		drv_version;

	    /* read in the first 16 bytes of the driver information block */
	    if(FD_read(file, lshared->driver_addr, 16, dbuf) == FAIL) {
		error_push(ERR_FILE, ERR_NONE_SEC, 
		    "Superblock v.0/1:Unable to read in the first 16 bytes of Driver Information Block.", 
		    LOGI_SUPER_BASE+lshared->driver_addr, NULL);
		CK_SET_RET_DONE(FAIL)
	    }

	    drv_start = dbuf;
	    p = dbuf;
	    logical = get_logical_addr(p, drv_start, lshared->driver_addr);
	    drv_version = *p++;
	    if(drv_version != DRIVERINFO_VERSION) {
		badinfo = drv_version;
		error_push(ERR_LEV_0, ERR_LEV_0B, "Superblock v.0/1:Driver Information Block version number should be 0", 
		    logical, &badinfo);
		CK_SET_RET(FAIL)
	    }
	    p += 3; /* reserved */

	    /* Driver information size */
	    UINT32DECODE(p, driver_size);

	    /* Driver name and/or version */
	    strncpy(driver_name, (const char *)p, (size_t)8);
	    driver_name[8] = '\0';
	    set_driver_id(&(lshared->driverid)/*out*/, driver_name);
	    p += 8; /* advance past driver identification */

	    logical = get_logical_addr(p, drv_start, lshared->driver_addr);

	    if((driver_size+DRVINFOBLOCK_HDR_SIZE) > sizeof(dbuf)) {
		error_push(ERR_LEV_0, ERR_LEV_0B, "Superblock v.0/1:Invalid size for Driver Information Block", 
		    logical, &badinfo);
		CK_SET_RET(FAIL)
	    }

	    if(FD_read(file, lshared->driver_addr+16, driver_size, p) == FAIL) {
		error_push(ERR_FILE, ERR_NONE_SEC, "Superblock v.0/1:Unable to read Driver Information", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }
	    if(decode_driver(lshared, p) != SUCCEED) {
		error_push(ERR_FILE, ERR_NONE_SEC, "Superblock v.0/1:Unable to decode Driver Information", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }
	}  /* DONE with driver information block */
    } else if(super_vers == SUPERBLOCK_VERSION_2) { /* version 2 format */
	uint32_t	read_chksum;
	uint32_t	computed_chksum;

	if(debug_verbose())
	    printf("Validating version 2 superblock...\n");
	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	lshared->size_offsets = *p++;
	if(lshared->size_offsets != 2 && lshared->size_offsets != 4 &&
            lshared->size_offsets != 8 && lshared->size_offsets != 16 && 
	    lshared->size_offsets != 32) {
		error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.2:Invalid Size of Offsets", logical, NULL);
		CK_SET_RET(FAIL)
	}

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	lshared->size_lengths = *p++;
	if(lshared->size_lengths != 2 && lshared->size_lengths != 4 &&
            lshared->size_lengths != 8 && lshared->size_lengths != 16 && 
	    lshared->size_lengths != 32) {
		error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.2:Invalid Size of Lengths", logical, NULL);
		CK_SET_RET(FAIL)
	}

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	lshared->file_consist_flg = *p++;
	if(lshared->file_consist_flg & ~SUPER_ALL_FLAGS) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.2:Invalid file consistency flags.", logical, NULL);
	    CK_SET_RET(FAIL)
	}
	/* Determine the remaining size of the superblock for V0 and V1 */
	remain_size = SUPERBLOCK_REMAIN_SIZE(super_vers, lshared);
	if((fixed_size+variable_size+remain_size) > sizeof(buf)) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, 
		"Superblock v.2:Total size of super block is incorrect", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}

	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	if(FD_read(file, LOGI_SUPER_BASE+fixed_size+variable_size, remain_size, p) == FAIL) {
	    error_push(ERR_FILE, ERR_NONE_SEC, 
		"Superblock v.2:Unable to read in the remaining size portion of the superblock", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
	/* Base address */
	addr_decode(lshared, (const uint8_t **)&p, &lshared->base_addr/*out*/);
	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);

	/* superblock extension address */
    	addr_decode(lshared, (const uint8_t **)&p, &lshared->extension_addr/*out*/);
	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);

	/* end of file address */
    	addr_decode(lshared, (const uint8_t **)&p, &lshared->stored_eoa);
	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);
	if((lshared->stored_eoa == CK_ADDR_UNDEF) || (lshared->base_addr >= lshared->stored_eoa)) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.2:Invalid End of File Address", logical, NULL);
	    CK_SET_RET(FAIL)
	}

	/* root group object header address */
	root_ent = calloc(1, sizeof(GP_entry_t));
	lshared->root_grp = root_ent;
    	addr_decode(lshared, (const uint8_t **)&p, &lshared->root_grp->header);
	logical = get_logical_addr(p, start_buf, LOGI_SUPER_BASE);

	computed_chksum = checksum_metadata(buf, (ck_size_t)(p - buf), 0);
        UINT32DECODE(p, read_chksum);
	if(computed_chksum != read_chksum) {
	    error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock v.2:Bad checksum", logical, NULL);
	    CK_SET_RET(FAIL)
	}

#ifdef DEBUG
	printf("base_addr=%llu, super_addr=%llu, stored_eoa=%llu, extension_addr=%llu\n",
		lshared->base_addr, lshared->super_addr, lshared->stored_eoa, 
		lshared->extension_addr);
	printf("size_offsets=%d, size_lengths=%d\n",	
		lshared->size_offsets, lshared->size_lengths);
	printf("header_address=%llu\n", lshared->root_grp->header);
#endif
    } else  /* error should have pushed already onto the error stack */
	CK_SET_RET(FAIL)

    /* check for errors before proceeding further */
    if (ret_value != SUCCEED)
	goto done;

/* should I return mesg ptr instead of idx from find_in_ohdr()? */
    lshared->sohm_tbl = NULL;
    
    if(addr_defined(lshared->extension_addr) && g_format_num != FORMAT_ONE_EIGHT) {
	error_push(ERR_LEV_0, ERR_LEV_0A, "Superblock:extension should not exist for this library version", logical, &badinfo);
	CK_SET_RET(FAIL)
    }

    if(addr_defined(lshared->extension_addr)) {
	OBJ_t	*oh = NULL;
	int	idx;

	if(debug_verbose()) {
	    printf("VALIDATING Superblock extension...\n");
#ifdef DEBUG
	    printf("GOING to validate the object header for the superblock extension at %llu...\n", lshared->extension_addr);
#endif
	}

	if(check_obj_header(file, lshared->extension_addr, &oh) < 0 || !oh)
	    CK_SET_RET_DONE(FAIL)

	if((idx = find_in_ohdr(file, oh, OBJ_SHMESG_ID)) < SUCCEED)
	    ;
#ifdef DEBUG
	    printf("There is no OBJ_SHMESG_ID in superblock extension.\n");
#endif
	else {
	    assert(oh->mesg[idx].native); 
#ifdef DEBUG
	    printf("Found OBJ_SHMESG_ID in superblock extension.\n");
#endif
	}
	if ((idx = find_in_ohdr(file, oh, OBJ_BTREEK_ID)) < SUCCEED) {
#ifdef DEBUG
	    printf("Unable to find OBJ_BTREEK_ID in superblock extension...set default SNODE/ISTORE values\n");
#endif
	    lshared->btree_k[BT_SNODE_ID] = BT_SNODE_K;
	    lshared->btree_k[BT_ISTORE_ID] = BT_ISTORE_K;
	    lshared->gr_leaf_node_k = CRT_SYM_LEAF_DEF;
	} else {
#ifdef DEBUG
	    printf("Found OBJ_BTREEK_ID in superblock extension...set SNODE/ISTORE values\n");
#endif
	    assert(oh->mesg[idx].native); 
	    lshared->btree_k[BT_SNODE_ID] = ((OBJ_btreek_t *)(oh->mesg[idx].native))->btree_k[BT_SNODE_ID];
	    lshared->btree_k[BT_ISTORE_ID] = ((OBJ_btreek_t *)(oh->mesg[idx].native))->btree_k[BT_ISTORE_ID];
	    lshared->gr_leaf_node_k = ((OBJ_btreek_t *)(oh->mesg[idx].native))->sym_leaf_k;
	}

	if ((idx = find_in_ohdr(file, oh, OBJ_DRVINFO_ID)) < SUCCEED) {
#ifdef DEBUG
	    printf("Unable to find OBJ_DRVINFO_ID in superblock extension.\n");
#endif
	} else {
#ifdef DEBUG
	    printf("Found OBJ_DRVINFO_ID in superblock extension...set_driver()\n");
#endif
	    assert(oh->mesg[idx].native); 
	    set_driver_id(&(lshared->driverid)/*out*/, ((OBJ_drvinfo_t *)(oh->mesg[idx].native))->name);

	    p = ((OBJ_drvinfo_t *)(oh->mesg[idx].native))->buf;

	    if(decode_driver(lshared, p) != SUCCEED) {
		error_push(ERR_FILE, ERR_NONE_SEC, "Superblock v.0/1:Unable to decode Driver Information", logical, NULL);
		CK_SET_RET_DONE(FAIL)
	    }
	}
    }

done:
    return(ret_value);
} /* check_superblock */

/*
 * Validate  SNOD
 */
static ck_err_t
check_sym(driver_t *file, ck_addr_t sym_addr, uint8_t *heap_chunk, name_list_t *name_list) 
{
    size_t 	size = 0;
    uint8_t	*buf = NULL;
    uint8_t	*p;
    unsigned	nsyms, u;
    GP_node_t	*sym = NULL;
    GP_entry_t	*ent;
    GP_entry_t	*prev_ent;
    int		sym_version, badinfo, ret;

    ck_err_t 	ret_value = SUCCEED;
    ck_err_t	ret_err = 0;
    ck_err_t	ret_other_err = 0;

    assert(file);
    assert(addr_defined(sym_addr));

    if(debug_verbose())
	printf("VALIDATING the Symbol table node at logical address %llu...\n", sym_addr);
    size = gp_node_size(file->shared);

    if((buf = malloc(size)) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Symbol table node:Internal allocation error", sym_addr, NULL);
	CK_INC_ERR_DONE
    }

    if((sym = calloc(1, sizeof(GP_node_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Symbol table node:Internal allocation error", 
	    sym_addr, NULL);
	CK_INC_ERR_DONE
    }

    sym->entry = malloc(2*SYM_LEAF_K(file->shared)*sizeof(GP_entry_t));
    if(sym->entry == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Symbol table node:Internal allocation error", 
	    sym_addr, NULL);
	CK_INC_ERR_DONE
    }

    if(FD_read(file, sym_addr, size, buf) == FAIL) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Symbol table node:Unable to read in the node", 
	    sym_addr, NULL);
	CK_INC_ERR_DONE
    }

    p = buf;
    if(memcmp(p, SNODE_MAGIC, SNODE_SIZEOF_MAGIC)) {
	error_push(ERR_LEV_1, ERR_LEV_1B, "Symbol table node:Could not find signature.", sym_addr, NULL);
	CK_INC_ERR
    } else if(debug_verbose())
	printf("FOUND Symbol table node signature.\n");

    p += 4;
    sym_version = *p++;
    if(SNODE_VERS != sym_version) {
	badinfo = sym_version;
	error_push(ERR_LEV_1, ERR_LEV_1B, "Symbol table node:Version should be 1", sym_addr, &badinfo);
	CK_INC_ERR
    }

    /* reserved */
    p++;

    /* number of symbols */
    UINT16DECODE(p, sym->nsyms);
    if(sym->nsyms > (2 * SYM_LEAF_K(file->shared))) {
	error_push(ERR_LEV_1, ERR_LEV_1B, "Symbol table node:Number of symbols exceeds (2*Group Leaf Node K)", 
	    sym_addr, NULL);
	CK_INC_ERR
    }

    /* reading the vector of symbol table group entries */
    if (gp_ent_decode_vec(file->shared, (const uint8_t **)&p, sym->entry, sym->nsyms) != SUCCEED)
	CK_INC_ERR_DONE

    if(!heap_chunk)
	printf("Warning: Symbol table node: invalid heap address--name not validated\n");

    /* validate symbol table group entries here  */
    prev_ent = sym->entry;
    for(u = 0, ent = sym->entry; u < sym->nsyms; u++, ent++) {
	char *s1, *s2, *sym_name;

	sym_name = (char *)(heap_chunk + HL_SIZEOF_HDR(file->shared) + ent->name_off);
	if(name_list_search(name_list, sym_name)) {
	    error_push(ERR_LEV_1, ERR_LEV_1C, "Symbol table node entry:Duplicate name", sym_addr, NULL);
		CK_INC_ERR
	} else
	    name_list_insert(name_list, sym_name);

	if(u && heap_chunk) {
	    s1 = (char *)(heap_chunk + HL_SIZEOF_HDR(file->shared) + prev_ent->name_off);
	    s2 = (char *)(heap_chunk + HL_SIZEOF_HDR(file->shared) + ent->name_off);

	    if(strcmp(s1, s2) >= 0) {
		error_push(ERR_LEV_1, ERR_LEV_1C, "Symbol table node entry:Name out of order", sym_addr, NULL);
		CK_INC_ERR
	    }
	}

	prev_ent = ent;

	/* for symbolic link: the object header address is undefined */
	if(ent->type != GP_CACHED_SLINK) {
	    if(ent->header==CK_ADDR_UNDEF) {
		error_push(ERR_LEV_1, ERR_LEV_1C, "Symbol table node entry:Undefined object header address.", 
		    sym_addr, NULL);
		CK_INC_ERR
	    }
	
	    if (ent->type < 0) {
		error_push(ERR_LEV_1, ERR_LEV_1C, "Symbol table node entry:Invalid cache type", sym_addr, NULL);
		CK_INC_ERR
	    }
	}
    } /* end for */

    for(u = 0, ent = sym->entry; u < sym->nsyms; u++, ent++) {
	if(ent->type != GP_CACHED_SLINK && ent->header != CK_ADDR_UNDEF)
	    if(check_obj_header(file, ent->header, NULL) < 0)
		++ret_other_err;
    } /* end for */

done:
    if(buf) free(buf);
    if(sym) {
	if(sym->entry) free(sym->entry);
	free(sym);
    }

    /* Flush out errors if any */
    if(ret_err && !object_api()) {
	error_print(stderr, file);
	error_clear();
    }

    if(ret_err || ret_other_err)
	ret_value = FAIL;

    return(ret_value);
} /* check_sym() */


/*
 * Validate version 1 B-tree
 */
static ck_err_t
check_btree(driver_t *file, ck_addr_t btree_addr, unsigned ndims, uint8_t *heap_chunk, name_list_t *name_list, void *lt_key, void *rt_key)
{
    ck_err_t ret_value = SUCCEED;
    ck_err_t ret_err = 0;	/* error from this routine */
    ck_err_t ret_other_err = 0;	/* track errors returned from other routines */

    uint8_t	*buf=NULL, *buffer=NULL;
    uint8_t	*p, nodetype;
    unsigned	nodelev, entries, u;
    ck_addr_t   left_sib;   /*address of left sibling */
    ck_addr_t   right_sib;  /*address of left sibling */
    size_t	hdr_size, key_size, key_ptr_size;
    ck_addr_t	child;
    uint8_t	*start_buf;
    ck_addr_t	logical;
    int		badinfo;

    assert(file);
    assert(addr_defined(btree_addr));

    hdr_size = BT_SIZEOF_HDR(file->shared);

    if(debug_verbose())
	printf("VALIDATING version 1 btree at logical address %llu...\n", btree_addr);

    if((buf = malloc(hdr_size)) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "version 1 B-tree:Internal allocation error", btree_addr, NULL);
	CK_INC_ERR
    }

    if(FD_read(file, btree_addr, hdr_size, buf) == FAIL) {
	error_push(ERR_LEV_1, ERR_LEV_1A1, 
	    "version 1 B-tree:Unable to read B-tree header", logical, NULL);
	CK_INC_ERR_DONE
    }

    start_buf = buf;
    logical = get_logical_addr(p, start_buf, btree_addr);
    p = buf;

    /* magic number */
    if(memcmp(p, BT_MAGIC, (size_t)BT_SIZEOF_MAGIC) != 0) {
	error_push(ERR_LEV_1, ERR_LEV_1A1, "version 1 B-tree:Could not find B-tree signature", btree_addr, NULL);
	CK_INC_ERR
    } else if (debug_verbose())
	printf("FOUND version 1 btree signature.\n");

    p += 4;
    logical = get_logical_addr(p, start_buf, btree_addr);
    nodetype = *p++;
    if((nodetype != 0) && (nodetype != 1)) {
	badinfo = nodetype;
	error_push(ERR_LEV_1, ERR_LEV_1A1, "Version 1 B-tree:Node Type should be 0 or 1", logical, &badinfo);
	CK_INC_ERR
    }

    logical = get_logical_addr(p, start_buf, btree_addr);
    nodelev = *p++;
    if(nodelev < 0) {
	badinfo = nodelev;
	error_push(ERR_LEV_1, ERR_LEV_1A1, 
	    "Version 1 B-tree:Node Level should not be negative", logical, &badinfo);
	CK_INC_ERR
    }

    logical = get_logical_addr(p, start_buf, btree_addr);
    UINT16DECODE(p, entries);

    if(entries > ((2*file->shared->btree_k[nodetype])+1)) {
	error_push(ERR_LEV_1, ERR_LEV_1A1, "Version 1 B-tree: Entries should not exceed 2K+1", logical, &badinfo);
	CK_INC_ERR
    }

    logical = get_logical_addr(p, start_buf, btree_addr);
    addr_decode(file->shared, (const uint8_t **)&p, &left_sib/*out*/);

    logical = get_logical_addr(p, start_buf, btree_addr);
    addr_decode(file->shared, (const uint8_t **)&p, &right_sib/*out*/);

    /* the remaining node size: key + child pointer */
    key_size = node_key_g[nodetype]->get_sizeof_rkey(file->shared, ndims);
    key_ptr_size = (entries)*SIZEOF_ADDR(file->shared) + (entries+1)*key_size;

#ifdef BUG_CHECK
    key_ptr_size = (2*file->shared->btree_k[nodetype])*SIZEOF_ADDR(file->shared) +
	(2*(file->shared->btree_k[nodetype]+1))*key_size;
#endif

#ifdef DEBUG
	printf("btree_k[SNODE]=%u, btree_k[ISTORE]=%u\n", file->shared->btree_k[BT_SNODE_ID], 
		file->shared->btree_k[BT_ISTORE_ID]);
	printf("key_size=%d, key_ptr_size=%d\n", key_size, key_ptr_size);
#endif

    if((buffer = malloc(key_ptr_size)) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Version 1 B-tree:Internal allocation error for key+child", btree_addr, NULL);
	CK_INC_ERR_DONE
    }

    if(FD_read(file, btree_addr+hdr_size, key_ptr_size, buffer) == FAIL) {
	error_push(ERR_LEV_1, ERR_LEV_1A1, "Version 1 B-tree:Unable to read key+child", btree_addr, NULL);
	CK_INC_ERR_DONE
    }

    /* SNOD */
    if(!nodetype && !heap_chunk)
	printf("Warning: Version 1 B-tree: invalid heap address--name not validated\n");

    p = start_buf = buffer;
		
    for(u = 0; u < entries; u++) {
	/* decode left key */
	if(lt_key == NULL &&
	   node_key_g[nodetype]->decode(file->shared, ndims, (const uint8_t **)&p, &lt_key) == FAIL) {
		error_push(ERR_LEV_1, ERR_LEV_1A1, "Version 1 B-tree:Errors when decoding left key", 
		       logical, NULL);
		CK_INC_ERR
	}

	logical = get_logical_addr(p, start_buf, btree_addr+hdr_size);

       	/* Decode child address */
  	addr_decode(file->shared, (const uint8_t **)&p, &child/*out*/);

	/* decode right key */
	if(rt_key == NULL &&
	   node_key_g[nodetype]->decode(file->shared, ndims, (const uint8_t **)&p, &rt_key) == FAIL) {

		error_push(ERR_LEV_1, ERR_LEV_1A1, "Version 1 B-tree:Errors when decoding right key", 
			logical, NULL);
		CK_INC_ERR
	}

	if(node_key_g[nodetype]->cmp(file->shared, ndims, heap_chunk, lt_key, rt_key) >= 0) {
	    error_push(ERR_LEV_1, ERR_LEV_1A1, "Version 1 B-tree:left & right keys are out of order", logical, NULL);
	    CK_INC_ERR
	}
    
	if(nodelev > 0) {
	    if(check_btree(file, child, ndims, heap_chunk, name_list, NULL, NULL) < 0)
		++ret_other_err;
	} else if(!nodetype) /* SNOD */
	    if(check_sym(file, child, heap_chunk, name_list) < 0)
		++ret_other_err;

	lt_key = rt_key;
	rt_key = NULL;
    }  /* end for */

done:
    if(buf) free(buf);
    if(buffer) free(buffer);

    /* Flush out errors if any */
    if(ret_err && !object_api()) {
	error_print(stderr, file);
	error_clear();
    }

    if(ret_err || ret_other_err)
	ret_value = FALSE;

    return(ret_value);
} /* check_btree() */


/*
 * Validate local heap
 */
static ck_err_t
check_lheap(driver_t *file, ck_addr_t lheap_addr, uint8_t **ret_heap_chunk)
{
    uint8_t	hdr[52];
    size_t	hdr_size, data_seg_size;
    size_t  	next_free_off = HL_FREE_NULL;
    size_t	size_free_block, saved_offset;
    uint8_t	*heap_chunk=NULL;
    ck_addr_t	addr_data_seg;
    uint8_t	*start_buf = NULL;
    ck_addr_t	logical;
    int		lheap_version, badinfo;
    const uint8_t *p=NULL;

    ck_err_t ret_err = 0;
    ck_err_t ret_value = SUCCEED;

    assert(file);
    assert(addr_defined(lheap_addr));

    hdr_size = HL_SIZEOF_HDR(file->shared);

    if(debug_verbose())
	printf("VALIDATING the local heap at logical address %llu...\n", lheap_addr);

    if(hdr_size > sizeof(hdr)) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Local Heap:Invalid header size", lheap_addr, NULL);
	CK_INC_ERR_DONE
    }
    if(FD_read(file, lheap_addr, hdr_size, hdr) == FAIL) {
	error_push(ERR_FILE, ERR_NONE_SEC, 
	    "Local Heap:Unable to read local heap header", lheap_addr, NULL);
	CK_INC_ERR_DONE
    }

    start_buf = hdr;
    p = hdr;
    logical = get_logical_addr(p, start_buf, lheap_addr);

    /* magic number */
    if(memcmp(p, HL_MAGIC, HL_SIZEOF_MAGIC) != 0) {
	error_push(ERR_LEV_1, ERR_LEV_1D, 
	    "Local Heap:Could not find local heap signature", logical, NULL);
	CK_INC_ERR
    } else if (debug_verbose())
	printf("FOUND local heap signature.\n");

    p += HL_SIZEOF_MAGIC;

    logical = get_logical_addr(p, start_buf, lheap_addr);
    /* Version */
    lheap_version = *p++;
    if(lheap_version != HL_VERSION) {
	badinfo = lheap_version;
	error_push(ERR_LEV_1, ERR_LEV_1D, "Local Heap:version number should be 0", logical, &badinfo);
	CK_INC_ERR
    }

    /* Reserved */
    p += 3;

    logical = get_logical_addr(p, start_buf, lheap_addr);
    /* data segment size */
    DECODE_LENGTH(file->shared, p, data_seg_size);
    if(data_seg_size <= 0) {
	error_push(ERR_LEV_1, ERR_LEV_1D, 
	    "Local Heap:Invalid data segment size", logical, NULL);
	CK_INC_ERR_DONE
    }

    /* offset to head of free-list */
    DECODE_LENGTH(file->shared, p, next_free_off);

    /* address of data segment */
    logical = get_logical_addr(p, start_buf, lheap_addr);
    addr_decode(file->shared, &p, &addr_data_seg);

#if 0
    addr_data_seg = addr_data_seg + file->shared->base_addr;
#endif
    if(addr_data_seg == CK_ADDR_UNDEF) {
	error_push(ERR_LEV_1, ERR_LEV_1D, 
	    "Local Heap:Address of data segment is undefined", logical, NULL);
	CK_INC_ERR_DONE
    }

    if((heap_chunk = malloc(hdr_size+data_seg_size)) == NULL) {
	error_push(ERR_LEV_1, ERR_LEV_1D, 
	    "Local Heap:Internal allocation error for data segment", logical, NULL);
	CK_INC_ERR_DONE
    }
	
#ifdef DEBUG
    printf("data_seg_size=%d, next_free_off =%d, addr_data_seg=%d\n",
	    data_seg_size, next_free_off, addr_data_seg);
#endif

    if(data_seg_size) {
	if(FD_read(file, addr_data_seg, data_seg_size, heap_chunk+hdr_size) == FAIL) {
	    error_push(ERR_FILE, ERR_NONE_SEC, 
		"Local Heap:Unable to read data segment", logical, NULL);
	    CK_INC_ERR_DONE
	}
    } 

    start_buf = heap_chunk+hdr_size;

    /* traverse the free list */
    while(next_free_off != HL_FREE_NULL) {
	if(next_free_off >= data_seg_size) {
	    error_push(ERR_LEV_1, ERR_LEV_1D, 
		"Local Heap:Offset of the next free block is invalid", logical, NULL);
	    CK_INC_ERR_DONE
	}
	saved_offset = next_free_off;
	p = heap_chunk + hdr_size + next_free_off;
	DECODE_LENGTH(file->shared, p, next_free_off);
	get_logical_addr(p, start_buf, addr_data_seg);
	DECODE_LENGTH(file->shared, p, size_free_block);
#if DEBUG
	printf("next_free_off=%d, size_free_block=%d\n",
	    next_free_off, size_free_block);
#endif
	if(!(size_free_block >= (2*SIZEOF_SIZE(file->shared)))) {
	    error_push(ERR_LEV_1, ERR_LEV_1D, 
		"Local Heap:Offset of the next free block is invalid", logical, NULL);
	    CK_INC_ERR_DONE
	}
        if(saved_offset + size_free_block > data_seg_size) {
	    error_push(ERR_LEV_1, ERR_LEV_1D, 
		"Local Heap:Bad heap free list", logical, NULL);
	    CK_INC_ERR_DONE
	}
    } /* end while */

done:
    if(!ret_err && ret_heap_chunk)
	*ret_heap_chunk = (uint8_t *)heap_chunk;
    else {
	if(heap_chunk)
	    free(heap_chunk);
	if(ret_err) {
	    ret_value = FAIL;
	    if(ret_heap_chunk)
		*ret_heap_chunk = NULL;
	    /* Flush out the errors if any */
	    if(!object_api()) {
		error_print(stderr, file);
		error_clear();
	    }
	}
    }
    return(ret_value);
} /* check_lheap() */


#ifdef NOTUSED
/* NEED CHECK: NOT CALLED anywhere yet */
static ck_err_t
check_gheap(driver_t *file, ck_addr_t gheap_addr, uint8_t **ret_heap_chunk)
{
    H5HG_heap_t *heap = NULL;
    uint8_t     *p = NULL;
    size_t      nalloc, need;
    size_t      max_idx=0;
    int		i, ret_value=SUCCEED; 

    uint8_t	*start_buf;
    ck_addr_t	logical;
    int		gheap_version, badinfo;

    /* check arguments */
    assert(addr_defined(gheap_addr));

    /* Read the initial 4k page */
    if((heap = calloc(1, sizeof(H5HG_heap_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Global Heap:Internal allocation error", gheap_addr, NULL);
	    CK_SET_RET_DONE(FAIL)
    }
    heap->addr = gheap_addr;
    if((heap->chunk = malloc(H5HG_MINSIZE))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Gloabl Heap:Internal allocation error", gheap_addr, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    if(debug_verbose())
	printf("VALIDATING the global heap at logical address %llu...\n", gheap_addr);
        
    if(FD_read(file, gheap_addr, H5HG_MINSIZE, heap->chunk) == FAIL) {
	error_push(ERR_FILE, ERR_NONE_SEC, 
	    "Global Heap:Unable to read collection", gheap_addr, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Magic number */
    if(memcmp(heap->chunk, H5HG_MAGIC, H5HG_SIZEOF_MAGIC)) {
	error_push(ERR_LEV_1, ERR_LEV_1E, 
	    "Global Heap:Could not find GCOL signature", gheap_addr, NULL);
	CK_SET_RET_DONE(FAIL)
    } else if(debug_verbose())
	printf("FOUND GLOBAL HEAP SIGNATURE\n");

    start_buf = heap->chunk;
    p = heap->chunk + H5HG_SIZEOF_MAGIC;
    logical = get_logical_addr(p, start_buf, gheap_addr);

    /* Version */
    gheap_version = *p++;
    if(H5HG_VERSION != gheap_version) {
	badinfo = gheap_version;
	error_push(ERR_LEV_1, ERR_LEV_1E, 
	    "Global Heap:version number should be 1", logical, &badinfo);
	CK_SET_RET(FAIL)
    } else if (debug_verbose())
	printf("Version 1 of global heap is detected\n");

    /* Reserved */
    p += 3;

    logical = get_logical_addr(p, start_buf, gheap_addr);
    /* Size */
    DECODE_LENGTH(file->shared, p, heap->size);
    assert(heap->size>=H5HG_MINSIZE);

    if(heap->size > H5HG_MINSIZE) {
       	ck_addr_t next_addr = gheap_addr + (ck_size_t)H5HG_MINSIZE;

    	if((heap->chunk = realloc(heap->chunk, heap->size))==NULL) {
	    error_push(ERR_LEV_1, ERR_LEV_1E, 
		"Global Heap:Internal allocation error", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
	if (FD_read(file, next_addr, heap->size-H5HG_MINSIZE, heap->chunk+H5HG_MINSIZE) == FAIL) {
	    error_push(ERR_FILE, ERR_NONE_SEC, 
		"Global Heap:Unable to read global heap collection", logical, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
    }  /* end if */

    /* Decode each object */
    p = heap->chunk + H5HG_SIZEOF_HDR(file->shared);
    logical = get_logical_addr(p, start_buf, gheap_addr);
    nalloc = H5HG_NOBJS(file->shared, heap->size);
    if((heap->obj = malloc(nalloc*sizeof(H5HG_obj_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Global Heap:Internal allocation error", logical, NULL);
	CK_SET_RET_DONE(FAIL)
    }
    heap->obj[0].size = heap->obj[0].nrefs=0;
    heap->obj[0].begin = NULL;

    heap->nalloc = nalloc;
    while(p < heap->chunk+heap->size) {
        if(p+H5HG_SIZEOF_OBJHDR(file->shared) > heap->chunk+heap->size) {
           /*
            * The last bit of space is too tiny for an object header, so we
            * assume that it's free space.
            */
           assert(NULL==heap->obj[0].begin);
           heap->obj[0].size = (heap->chunk+heap->size) - p;
           heap->obj[0].begin = p;
           p += heap->obj[0].size;
	} else {
	    unsigned idx;
	    uint8_t *begin = p;

	    UINT16DECODE (p, idx);

            /* Check if we need more room to store heap objects */
            if(idx >= heap->nalloc) {
		size_t new_alloc;       /* New allocation number */
		H5HG_obj_t *new_obj;    /* New array of object descriptions */

		/* Determine the new number of objects to index */
		new_alloc = MAX(heap->nalloc*2,(idx+1));

		logical = get_logical_addr(p, start_buf, gheap_addr);
		/* Reallocate array of objects */
		if((new_obj=realloc(heap->obj, new_alloc*sizeof(H5HG_obj_t)))==NULL) {
		    error_push(ERR_LEV_1, ERR_LEV_1E, 
			"Global Heap:Internal allocation error", logical, NULL);
		    CK_SET_RET_DONE(FAIL)
		}

                /* Update heap information */
                heap->nalloc = new_alloc;
                heap->obj = new_obj;
            } /* end if */

	    UINT16DECODE (p, heap->obj[idx].nrefs);
	    p += 4; /*reserved*/
	    DECODE_LENGTH (file->shared, p, heap->obj[idx].size);
	    heap->obj[idx].begin = begin;
#ifdef DEBUG
	    if(idx > 0)
		printf("global:%s\n", p);
#endif
            /*
             * The total storage size includes the size of the object header
             * and is zero padded so the next object header is properly
             * aligned. The last bit of space is the free space object whose
             * size is never padded and already includes the object header.
             */
            if (idx > 0) {
		need = H5HG_SIZEOF_OBJHDR(file->shared) + H5HG_ALIGN(heap->obj[idx].size);

		/* Check for "gap" in index numbers (caused by deletions) 
	  	   and fill in heap object values */
          	if(idx > (max_idx+1))
		    memset(&heap->obj[max_idx+1],0,sizeof(H5HG_obj_t)*(idx-(max_idx+1)));
		    max_idx = idx;
            } else {
		need = heap->obj[idx].size;
            }
	    p = begin + need;
	}  /* end else */
    }  /* end while */

    assert(p==heap->chunk+heap->size);
    assert(H5HG_ISALIGNED(heap->obj[0].size));

    /* Set the next index value to use */
    if(max_idx > 0)
       	heap->nused = max_idx + 1;
    else
       	heap->nused = 1;

done:
    if((ret_value == SUCCEED) && ret_heap_chunk) {
	*ret_heap_chunk = (uint8_t *)heap;
    } else {
	if (ret_value == FAIL) {
	    *ret_heap_chunk = NULL;
	    if(heap) {
		if (heap->chunk) free(heap->chunk);
		if (heap->obj) free(heap->obj);
		free(heap);
	    }
	}
    }
    return(ret_value);
} /* end check_gheap() */
#endif

/*
 * Callback function from check_btree2() for indexed group:
 * 	Get info about the message from the btree record _RECORD
 *	Read the message from fractal heap 
 *	Decode the message
 */
static ck_err_t
G_dense_ck_fh_msg_cb(driver_t *file, const void *_record, void *_ck_udata)
{
    G_dense_bt2_name_rec_t *rec = (G_dense_bt2_name_rec_t *)_record;
    HF_hdr_t *fhdr = (HF_hdr_t *)_ck_udata;
    ck_err_t ret_value = SUCCEED;      /* Return value */
    obj_info_t objinfo;
    uint8_t	*mesg_ptr;

    assert(fhdr);
    assert(fhdr->man_dtable.table_addr);

    if(HF_get_obj_info(file, fhdr, &rec->id, &objinfo) < SUCCEED) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Dense msg cb:cannot get fractal heap ID info", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    mesg_ptr = malloc(objinfo.size);
    if(HF_read(file, fhdr, &rec->id, mesg_ptr, &objinfo) < SUCCEED) {
	error_push(ERR_FILE, ERR_NONE_SEC, "Dense msg cb:Unable to read message from fractal heap", 
		fhdr->heap_addr, NULL);
	CK_SET_RET_DONE(FAIL)
    }
    if(message_type_g[OBJ_LINK_ID]->decode(file, mesg_ptr, NULL, -1) == NULL) {
	error_push(ERR_LEV_2, ERR_LEV_2A, "Dense msg cb:Errors found when decoding message from fractal heap", fhdr->heap_addr, NULL);
	CK_SET_RET_DONE(FAIL)
    }

done:
    return(ret_value);
} /* G_dense_ck_fh_msg_cb() */

/*
 * Callback function from check_btree2() for indexed attribute:
 *	1. Shared indexed attribute: not implemented yet 
 *		(need to get it from the shared message table)
 * 	2. Get info about the message from the btree record _RECORD
 *	   Read the message from fractal heap 
 *	   Decode the message
 */
static ck_err_t
A_dense_ck_fh_msg_cb(driver_t *file, const void *_record, void *_ck_udata)
{
    A_dense_bt2_name_rec_t *rec = (A_dense_bt2_name_rec_t *)_record;
    HF_hdr_t *fhdr = (HF_hdr_t *)_ck_udata;
    ck_err_t ret_value = SUCCEED;      /* Return value */
    obj_info_t objinfo;
    uint8_t	*mesg_ptr;

    if(rec->flags & OBJ_MSG_FLAG_SHARED)
	printf("Warning: Callback for shared indexed attributes not implemented yet...\n");
    else {
	assert(fhdr);
	assert(fhdr->man_dtable.table_addr);
	if(HF_get_obj_info(file, fhdr, &rec->id, &objinfo) < SUCCEED) {
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, 
		"Indexed attribute cb:cannot get fractal heap ID info", -1, NULL);
	    CK_SET_RET_DONE(FAIL)
	}

	mesg_ptr = malloc(objinfo.size);
	if(HF_read(file, fhdr, &rec->id, mesg_ptr, &objinfo) < SUCCEED) {
	    error_push(ERR_FILE, ERR_NONE_SEC, "Indexed attribute cb:Unable to read message from fractal heap", 
		fhdr->heap_addr, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
    }

done:
    return(ret_value);
} /* A_dense_ck_fh_msg_cb() */


/*
 * Call routines to decode annd validate messages
 * Do further validation needed for some messages
 * Both version 1 & 2
 */
static ck_err_t
decode_validate_messages(driver_t *file, OBJ_t *oh)
{
    ck_err_t ret_value = SUCCEED;
    ck_err_t ret_err = 0;
    ck_err_t ret_other_err = 0;

    int		i, k;
    unsigned	id;
    void	*mesg;
    uint8_t	*start_buf, *p;
    ck_addr_t	logical, logi_base;

    uint8_t	*heap_chunk;
    unsigned	ndims = 0, local_ndims = 0;
    ck_addr_t	btree_addr;
    name_list_t *name_list;

    assert(file);
    assert(oh);

    for(i = 0; i < oh->nmesgs; i++) {
	start_buf = oh->chunk[oh->mesg[i].chunkno].image;
	p = oh->mesg[i].raw;
	logi_base = oh->chunk[oh->mesg[i].chunkno].addr;
	logical = get_logical_addr(p, start_buf, logi_base);


	id = oh->mesg[i].type->id;
	if(id == OBJ_CONT_ID || id == OBJ_NIL_ID)
	    continue;
	else if(id == OBJ_UNKNOWN_ID) {
	    error_push(ERR_LEV_2, ERR_LEV_2A, "Unsupported message encountered", logical, NULL);
	    CK_INC_ERR_CONTINUE
	} else if(oh->mesg[i].flags & OBJ_FLAG_SHARED) {
	    if((mesg = OBJ_shared_decode(file, oh->mesg[i].raw, oh->mesg[i].type, start_buf, logi_base)) != NULL)
		mesg = OBJ_shared_read(file, mesg, oh->mesg[i].type);
	} else
	    mesg = message_type_g[id]->decode(file, oh->mesg[i].raw, start_buf, logi_base);

	if((oh->mesg[i].native = mesg) == NULL) {
	    error_push(ERR_LEV_2, ERR_LEV_2A, "Errors found when decoding message", logical, NULL);
	    CK_INC_ERR_CONTINUE
	}

	switch(id) {
	    case OBJ_SDS_ID:
	    case OBJ_DT_ID:
	    case OBJ_FILL_OLD_ID:
	    case OBJ_FILL_ID:
	    case OBJ_FILTER_ID:
	    case OBJ_BOGUS_ID:
	    case OBJ_MDT_ID:
	    case OBJ_MDT_OLD_ID:
	    case OBJ_COMM_ID:
	    case OBJ_ATTR_ID:
		break;

	    case OBJ_EDF_ID:
	    {
		OBJ_edf_t *edf = (OBJ_edf_t *)mesg;

		ret_other_err = check_lheap(file, edf->heap_addr, &heap_chunk);
		if(heap_chunk) {
		    char  *s = NULL;

		    for(k = 0; k < edf->nused; k++) {
			s = (char*)(heap_chunk + HL_SIZEOF_HDR(file->shared) + edf->slot[k].name_offset);
			if(s) {
			    if(!(*s)) {
				error_push(ERR_LEV_2, ERR_LEV_2A2h, 
				    "Invalid external file name found in local heap", logical, NULL);
				CK_INC_ERR
			    }
			} else {
			    error_push(ERR_LEV_2, ERR_LEV_2A2h, "Invalid name offset into local heap", logical, NULL);
			    CK_INC_ERR
			}
		    }  
		    free(heap_chunk);
		}
	    }
	    break;

	    case OBJ_LAYOUT_ID:
	    {
		OBJ_layout_t *layout = (OBJ_layout_t *)mesg;

		if(layout->type == DATA_CHUNKED) {
		    local_ndims = layout->u.chunk.ndims;
		    btree_addr = layout->u.chunk.addr;

		    /* may be undefined if storage is not allocated yet */
		    if(addr_defined(btree_addr) && check_btree(file, btree_addr, local_ndims, NULL, NULL, NULL, NULL) < 0)
			++ret_other_err;
		}
	    }
	    break;

	    case OBJ_GROUP_ID:
	    {
		OBJ_stab_t *stab = (OBJ_stab_t *)mesg;

		name_list_init(&name_list);
		if(check_lheap(file, stab->heap_addr, &heap_chunk) < 0)
		    ++ret_other_err;
		if(check_btree(file, stab->btree_addr, ndims, heap_chunk, name_list, NULL, NULL) < 0)
		    ++ret_other_err;

		name_list_dest(name_list);
		if(heap_chunk) free(heap_chunk);
	    }
	    break;

	    case OBJ_LINFO_ID:
	    {
		OBJ_linfo_t *linfo = (OBJ_linfo_t *)mesg;
		ck_op_t ck_fh_msg_op = NULL;
		HF_hdr_t *fhdr = NULL;

		if(g_format_num == FORMAT_ONE_SIX)
		    break;

		if(addr_defined(linfo->fheap_addr)) {
		    if(check_fheap(file, linfo->fheap_addr) != SUCCEED)
			++ret_other_err;
		    else {
			if((fhdr = HF_open(file, linfo->fheap_addr)) == NULL) {
			    error_push(ERR_INTERNAL, ERR_NONE_SEC, "Internal: Unable to open fractal heap", -1, NULL);
			    CK_INC_ERR
			} else 
			    ck_fh_msg_op = G_dense_ck_fh_msg_cb;
		    }
		} 

		if(addr_defined(linfo->corder_bt2_addr) &&
		   check_btree2(file, linfo->corder_bt2_addr, G_BT2_CORDER, ck_fh_msg_op, fhdr) != SUCCEED)
			++ret_other_err;

		if(addr_defined(linfo->name_bt2_addr) &&
		   check_btree2(file, linfo->name_bt2_addr, G_BT2_NAME, ck_fh_msg_op, fhdr) != SUCCEED)
			++ret_other_err;
		if(fhdr)
		    (void) HF_close(fhdr);
	    }

	    break;

	    case OBJ_GINFO_ID:
		break;

	    case OBJ_SHMESG_ID:
	    {
		OBJ_shmesg_table_t *shm = (OBJ_shmesg_table_t *)mesg;

		if(g_format_num == FORMAT_ONE_SIX)
		    break;

		if(addr_defined(shm->addr)) {
#ifdef DEBUG
		    printf("The SOHM table address from message is %llu; num of indices=%u\n",
			shm->addr, shm->nindexes);
#endif
		    if((check_SOHM(file, shm->addr, shm->nindexes)) != SUCCEED)
			++ret_other_err;
		}
	    }
	    break;

	    case OBJ_BTREEK_ID:
	    case OBJ_DRVINFO_ID:
		break;

	    case OBJ_AINFO_ID:
	    {
		OBJ_ainfo_t *ainfo = (OBJ_ainfo_t *)mesg;
		ck_op_t ck_fh_msg_op = NULL;
		HF_hdr_t *fhdr = NULL;

		if(g_format_num == FORMAT_ONE_SIX)
		    break;

		if(addr_defined(ainfo->fheap_addr)) {
		    if(check_fheap(file, ainfo->fheap_addr) != SUCCEED)
			++ret_other_err;
		    else {
			if((fhdr = HF_open(file, ainfo->fheap_addr)) == NULL) {
			    error_push(ERR_INTERNAL, ERR_NONE_SEC, "Internal: Unable to open fractal heap", -1, NULL);
			    CK_INC_ERR
			} else 
			    ck_fh_msg_op = A_dense_ck_fh_msg_cb;
		    }
		} 

		if(addr_defined(ainfo->corder_bt2_addr)) {
		    if(check_btree2(file, ainfo->corder_bt2_addr, A_BT2_CORDER, ck_fh_msg_op, fhdr) != SUCCEED)
			++ret_other_err;
		}

		if(addr_defined(ainfo->name_bt2_addr)) {
		    if(check_btree2(file, ainfo->name_bt2_addr, A_BT2_NAME, ck_fh_msg_op, fhdr) != SUCCEED)
			++ret_other_err;
		}
		if(fhdr)
		    (void) HF_close(fhdr);

	    }
	    break;

	    case OBJ_REFCOUNT_ID:
	    case OBJ_LINK_ID:
	    case OBJ_UNKNOWN_ID:
	    default:
		break;

	}  /* end switch on id */
    }  /* end for nmesgs */

done:
    /* flush out errors if any */
    if(ret_err && !object_api()) {
	error_print(stderr, file);
	error_clear();
    }

    if(ret_err || ret_other_err)
	ret_value = FAIL;

    return(ret_value);
} /* decode_validate_messages() */


/*
 * Find a message in the given object header with the desired type_id
 * failure: return FAIL
 * success: return the index in oh->mesg[] array
 */
static int
find_in_ohdr(driver_t *file, OBJ_t *oh, int type_id)
{
    const obj_class_t   *type = NULL;
    int			u, ret_value=SUCCEED;
    const uint8_t	*start_buf;
    ck_addr_t		logi_base;

    assert(file);
    assert(oh);

    /* Scan through the messages looking for the right one */
    for(u = 0; u < oh->nmesgs; u++) {
	if(oh->mesg[u].type->id == type_id)
	    break;
    }

    if(u == oh->nmesgs) {
	CK_SET_RET_DONE(FAIL)
    } else
	ret_value = u;
	
    start_buf = oh->chunk[oh->mesg[u].chunkno].image;
    logi_base = oh->chunk[oh->mesg[u].chunkno].addr;

    if(oh->mesg[u].native == NULL) {
	if(oh->mesg[u].flags & OBJ_FLAG_SHARED)
	    oh->mesg[u].native = OBJ_shared_decode(file, oh->mesg[u].raw, oh->mesg[u].type, start_buf, logi_base);
	else
	    oh->mesg[u].native = ((oh->mesg[u].type)->decode)(file, oh->mesg[u].raw, start_buf, logi_base);

        if(oh->mesg[u].native==NULL) {
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, 
		"find_in_ohdr:Unable to decode message", -1, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
    }

done:
    return(ret_value);
} /* end find_in_ohdr() */

/*
 * Read a message that has the OBJ_FLAG_SHARED bit on
 *   indicating that it is shared.
 * Handles version 1 & 2
 */
static void *
OBJ_shared_decode(driver_t *file, const uint8_t *buf, const obj_class_t *type, const uint8_t *start_buf, ck_addr_t logi_base)
{
    OBJ_shared_t    	*mesg=NULL;
    unsigned    	version;
    void        	*ret_ptr=NULL;
    int			ret_value=SUCCEED;
    int			badinfo;

    assert(file);
    assert(buf);
    assert(type);

    if((mesg = calloc(1, sizeof(OBJ_shared_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Internal Shared Message: Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Version */
    version = *buf++;
    /* version is adjusted if error for further continuation */
    if(g_format_num == FORMAT_ONE_SIX) {
	if(version != OBJ_SHARED_VERSION_1 && version != OBJ_SHARED_VERSION_2) {
	    badinfo = version;
	    version = OBJ_SHARED_VERSION_2;
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, "Shared Message:Bad version number", -1, &badinfo);
	    CK_SET_RET(FAIL)
	}
    } else {
	if(version < OBJ_SHARED_VERSION_1 || version > OBJ_SHARED_VERSION_LATEST) {
	    badinfo = version;
	    version = OBJ_SHARED_VERSION_LATEST;
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, "Shared Message:Bad version number", -1, &badinfo);
	    CK_SET_RET(FAIL)
	}
    }

    /* 
     * Get the shared information type
     * Flags are unused before version 3.
     */
    if(version >= OBJ_SHARED_VERSION_2)
        mesg->type = *buf++;
    else {
        mesg->type = OBJ_SHARE_TYPE_COMMITTED;
        buf++;
    } 

    /* Skip reserved bytes (for version 1) */
    if(version==OBJ_SHARED_VERSION_1)
       	buf += 6;

    if(version==OBJ_SHARED_VERSION_1) {
        mesg->u.loc.index = 0;
        buf += SIZEOF_SIZE(file->shared);  
        addr_decode(file->shared, &buf, &(mesg->u.loc.oh_addr));
    } else if(version >= OBJ_SHARED_VERSION_2) {
        if(mesg->type == OBJ_SHARE_TYPE_SOHM) {
	    if (version < OBJ_SHARED_VERSION_3) {
		error_push(ERR_INTERNAL, ERR_NONE_SEC, 
		    "Shared Message:Inconsistent message type and version", -1, &badinfo);
		CK_SET_RET(FAIL)
	    }
            memcpy(&mesg->u.heap_id, buf, sizeof(mesg->u.heap_id));
        } else {
            if(version < OBJ_SHARED_VERSION_3)
                mesg->type = OBJ_SHARE_TYPE_COMMITTED;

            mesg->u.loc.index = 0;
            addr_decode(file->shared, &buf, &mesg->u.loc.oh_addr);
        }
    }
    mesg->msg_type_id = type->id;

    if(mesg->type != OBJ_SHARE_TYPE_SOHM) {
	if (mesg->u.loc.oh_addr == CK_ADDR_UNDEF) {
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, "Shared Message:Invalid object header address", 
		-1, NULL);
	    CK_SET_RET_DONE(FAIL)
	}
    }
    
    if(ret_value == SUCCEED)
	ret_ptr = mesg;

done:
    if(ret_value == FAIL && mesg)
	free(mesg);
    return(ret_ptr);
} /* end OBJ_shared_decode() */

/*
 * Read the message pointed to by the message that has OBJ_FLAG_SHARED turned on
 * both version 1 & 2
 */
static void *
OBJ_shared_read(driver_t *file, OBJ_shared_t *obj_shared, const obj_class_t *type)
{
    OBJ_t	*oh = NULL;
    int		idx, status;
    ck_addr_t	fheap_addr;
    void	*ret_value=NULL;
    void	*mesg;

    HF_hdr_t 	*fhdr=NULL;
    ck_size_t	mesg_size;
    uint8_t	*mesg_ptr;
    uint8_t	*start_buf;
    ck_addr_t	logi_base;

    obj_info_t	objinfo;

    assert(file);
    assert(obj_shared);
    assert(type);

/* NEEd to initialize logi_base */
    if(obj_shared->type == OBJ_SHARE_TYPE_SOHM) {
	if(SM_get_fheap_addr(file, type->id, &fheap_addr) < SUCCEED) {
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, 
		"Internal Shared Read:Cannot get fractal heap address for shared message", -1, NULL);
	    CK_SET_RET_DONE(NULL)
	}
	if((fhdr = HF_open(file, fheap_addr)) == NULL) {
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, 
		"Internal Shared Read:Cannot open fractal heap header", -1, NULL);
	    CK_SET_RET_DONE(NULL)
	}

	if(HF_get_obj_info(file, fhdr, &(obj_shared->u.heap_id), &objinfo) < SUCCEED) {
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, 
		"Internal Shared Read:Cannot get info from fractal heap ID", -1, NULL);
	    CK_SET_RET_DONE(NULL)
	}

	mesg_ptr = malloc(objinfo.size);
	if(HF_read(file, fhdr, &(obj_shared->u.heap_id), mesg_ptr, &objinfo) < SUCCEED) {
	    error_push(ERR_FILE, ERR_NONE_SEC, "Internal Shared Read:Unable to read object from fractal heap", 
		-1, NULL);
	    CK_SET_RET_DONE(NULL)
	}

	ret_value = type->decode(file, mesg_ptr, start_buf, logi_base);

    } else if(obj_shared->type == OBJ_SHARE_TYPE_COMMITTED) {
	if(check_obj_header(file, obj_shared->u.loc.oh_addr, &oh) < 0 || !oh)
	   CK_SET_RET_DONE(NULL)

	idx = find_in_ohdr(file, oh, type->id);
	if(idx < SUCCEED) {
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, 
		"Internal Shared Read:Cannot find message type in object header", -1, NULL);
	    CK_SET_RET_DONE(NULL)
	}

	if(oh->mesg[idx].flags & OBJ_FLAG_SHARED) {
	    OBJ_shared_t *sh_shared;

	    sh_shared = (OBJ_shared_t *)(oh->mesg[idx].native);
	    ret_value = OBJ_shared_read(file, sh_shared, type);
	} else 
	    ret_value = oh->mesg[idx].native;
    } else {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Internal Shared Read:Invalid type", -1, NULL);
	CK_SET_RET_DONE(NULL)
    }
	
done:
    /* NEED CHECK:for return value? */
    if(fhdr)
	(void) HF_close(fhdr);
    return(ret_value);
} /* OBJ_shared_read() */


/* support routine for check_obj_header() */
static ck_err_t
OBJ_alloc_msgs(OBJ_t *oh, ck_size_t min_alloc)
{
    ck_size_t 	old_alloc;              /* Old number of messages allocated */
    ck_size_t 	na;                     /* New number of messages allocated */
    OBJ_mesg_t 	*new_mesg;              /* Pointer to new message array */
    ck_err_t 	ret_value = SUCCEED;    /* Return value */

    /* check args */
    assert(oh);

    old_alloc = oh->alloc_nmesgs;
    na = oh->alloc_nmesgs + MAX(oh->alloc_nmesgs, min_alloc);   /* At least double */

    /* Attempt to allocate more memory */
    if(NULL == (new_mesg = realloc(oh->mesg, na*sizeof(OBJ_mesg_t)))) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Object Header:Internal allocation error", -1, NULL);
	CK_SET_RET_DONE(FAIL)
    }

    /* Update ohdr information */
    oh->alloc_nmesgs = na;
    oh->mesg = new_mesg;

    /* Set new object header info to zeros */
    memset(&oh->mesg[old_alloc], 0, (oh->alloc_nmesgs - old_alloc) * sizeof(OBJ_mesg_t));
done:
    return(ret_value);
} /* OBJ_alloc_msgs() */

/*
 *  Checking object header for both version 1 & 2
 */
ck_err_t
check_obj_header(driver_t *file, ck_addr_t obj_head_addr, OBJ_t **ret_oh)
{
    size_t chunk_size;
    size_t spec_read_size;
    size_t prefix_size;
    uint8_t *p, flags, *start_buf;
    uint8_t buf[OBJ_SPEC_READ_SIZE];
    unsigned nmesgs;
    unsigned curmesg = 0;

    ck_addr_t chunk_addr;
    ck_addr_t logical, logi_base;
    ck_addr_t rel_eoa, abs_eoa;

    int version, nlink, idx, badinfo;

    ck_err_t ret_value = SUCCEED;  /* return value */
    ck_err_t ret_err = 0;	/* errors from the current routine */
    ck_err_t ret_other_err = 0;	/* track errors from other routines */

    OBJ_t *oh = NULL;
    OBJ_cont_t *cont;
    global_shared_t *lshared;


    /* FORMAT_ONE_EIGHT && object header version 2 SIGNATURE */
    int	format_objvers_two=FALSE;

    assert(file);
    assert(addr_defined(obj_head_addr));

    if(debug_verbose())
	printf("VALIDATING the object header at logical address %llu...\n", obj_head_addr);

    assert(addr_defined(obj_head_addr));
    idx = table_search(obj_head_addr);
    if(idx >= 0) { /* FOUND the object */
	if(obj_table->objs[idx].nlink > 0) {
	   obj_table->objs[idx].nlink--;
	    if (!ret_oh) goto done;
	} else {
	    error_push(ERR_LEV_2, ERR_LEV_2A, 
		"Object Header:Inconsistent reference count", obj_head_addr, NULL);
	    CK_INC_ERR_DONE
	}
    }

    lshared = file->shared;
    start_buf = buf;
	
    p = buf;
    if(CK_ADDR_UNDEF == (abs_eoa = FD_get_eof(file))) {
	error_push(ERR_FILE, ERR_NONE_SEC, 
	    "Object Header:Unable to determine file size", obj_head_addr, NULL);
	CK_INC_ERR_DONE
    }
    rel_eoa = abs_eoa - lshared->base_addr;
    spec_read_size = MIN(rel_eoa-obj_head_addr, OBJ_SPEC_READ_SIZE);

    if(FD_read(file, obj_head_addr, spec_read_size, buf) == FAIL) {
	error_push(ERR_FILE, ERR_NONE_SEC, 
	    "Object Header:Unable to read object header", obj_head_addr, NULL);
	CK_INC_ERR_DONE
    }

    if((oh = calloc(1, sizeof(OBJ_t))) == NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, 
	    "Object Header:Internal allocation error", obj_head_addr, NULL);
	CK_INC_ERR_DONE
    }

    if((g_format_num == FORMAT_ONE_EIGHT) && (!memcmp(p, OBJ_HDR_MAGIC, (size_t)OBJ_SIZEOF_MAGIC)))
	format_objvers_two = TRUE;

    if(format_objvers_two) {
	if(debug_verbose()) {
	   printf("VALIDATING version 2 object header ...\n");	
	   printf("FOUND Version 2 object header signature\n");	
	}

	start_buf = buf;
	p = buf;
	logical = get_logical_addr(p, start_buf, obj_head_addr);

	p += OBJ_SIZEOF_MAGIC;
        oh->version = *p++;
        if(OBJ_VERSION_2 != oh->version) {
    	    badinfo = oh->version;
	    error_push(ERR_LEV_2, ERR_LEV_2A1b, "version 2 Object Header:Bad version number", logical, &badinfo);
	    CK_INC_ERR
	}

	logical = get_logical_addr(p, start_buf, obj_head_addr);
	oh->flags = *p++;
	if(oh->flags & ~OBJ_HDR_ALL_FLAGS) {
	    error_push(ERR_LEV_2, ERR_LEV_2A1b, 
		"version 2 Object Header:Unknown object header status flags", logical, NULL);
	    CK_INC_ERR
	}

        /* Number of messages (to allocate initially) */
        nmesgs = 1;
	oh->nlink = 1;

        /* Time fields */
        if(oh->flags & OBJ_HDR_STORE_TIMES) {
            UINT32DECODE(p, oh->atime);
            UINT32DECODE(p, oh->mtime);
            UINT32DECODE(p, oh->ctime);
            UINT32DECODE(p, oh->btime);
        } /* end if */
        else
            oh->atime = oh->mtime = oh->ctime = oh->btime = 0;

        /* Attribute fields */
        if(oh->flags & OBJ_HDR_ATTR_STORE_PHASE_CHANGE) {
	    logical = get_logical_addr(p, start_buf, obj_head_addr);
            UINT16DECODE(p, oh->max_compact);
            UINT16DECODE(p, oh->min_dense);
            if(oh->max_compact < oh->min_dense) {
		error_push(ERR_LEV_2, ERR_LEV_2A1b, 
		    "version 2 Object Header:Invalid attribute phase changed values", logical, NULL);
		CK_INC_ERR
	    }
        } else {
            oh->max_compact = OBJ_CRT_ATTR_MAX_COMPACT_DEF;
            oh->min_dense = OBJ_CRT_ATTR_MIN_DENSE_DEF;
        } /* end else */

	/* First chunk size */
	logical = get_logical_addr(p, start_buf, obj_head_addr);
        switch(oh->flags & OBJ_HDR_CHUNK0_SIZE) {
            case 0:     /* 1 byte size */
                chunk_size = *p++;
                break;

            case 1:     /* 2 byte size */
                UINT16DECODE(p, chunk_size);
                break;

            case 2:     /* 4 byte size */
                UINT32DECODE(p, chunk_size);
                break;

            case 3:     /* 8 byte size */
                UINT64DECODE(p, chunk_size);
                break;

            default:
                printf("bad size for chunk 0\n");
        } /* end switch */
        if(chunk_size && chunk_size < OBJ_SIZEOF_MSGHDR_VERS(OBJ_VERSION_2, oh->flags & OBJ_HDR_ATTR_CRT_ORDER_TRACKED)) {
	    error_push(ERR_LEV_2, ERR_LEV_2A1b, 
		"version 2 Object Header:Bad object header size", logical, &badinfo);
	    CK_INC_ERR_DONE
	}
    } else { /* FORMAT_ONE_SIX or version 1 object header encountered for FORMAT_ONE_EIGHT */
	if(debug_verbose())
	   printf("VALIDATING version 1 object header...\n");	

	start_buf = buf;
	p = buf;
	logical = get_logical_addr(p, start_buf, obj_head_addr);

	oh->version = *p++;
	if(oh->version != OBJ_VERSION_1) {
    	    badinfo = oh->version;
	    error_push(ERR_LEV_2, ERR_LEV_2A1a, "Version 1 Object Header:Bad version number", logical, &badinfo);
	    CK_INC_ERR
	} else if(debug_verbose())
	    printf("Version 1 object header encountered\n");	

	/* Flags */ 
        oh->flags = OBJ_CRT_OHDR_FLAGS_DEF;

	p++;  /* reserved */

	logical = get_logical_addr(p, start_buf, obj_head_addr);
	UINT16DECODE(p, nmesgs);
	if((int)nmesgs < 0) {
	    badinfo = nmesgs;
	    error_push(ERR_LEV_2, ERR_LEV_2A1a, 
		"version 1 Object Header:number of header messages should not be negative", logical, &badinfo);
	    CK_INC_ERR
	}

	UINT32DECODE(p, oh->nlink);

	/* NEED CHECK: why nlink-1? */
	table_insert(obj_head_addr, oh->nlink-1);

	UINT32DECODE(p, chunk_size);
	p += 4;
    }

    prefix_size = (size_t)(p - buf);
    chunk_addr = obj_head_addr + prefix_size;

    logical = get_logical_addr(p, start_buf, chunk_addr);

    oh->alloc_nmesgs = (nmesgs > 0) ? nmesgs: 1;
    if((oh->mesg = calloc(oh->alloc_nmesgs, sizeof(OBJ_mesg_t)))==NULL) {
	error_push(ERR_INTERNAL, ERR_NONE_SEC, "Object Header:Internal allocation error", logical, NULL);
	CK_INC_ERR_DONE
    }

#ifdef DEBUG
    printf("oh->version =%d, nmesgs=%d, oh->nlink=%d\n", oh->version, nmesgs, oh->nlink);
    printf("chunk_addr=%llu, chunk_size=%u\n", chunk_addr, chunk_size);
#endif

    /* read each chunk from disk */
    while(addr_defined(chunk_addr)) {
	unsigned 	chunkno;	/* current chunk's index */
	uint8_t		*eom_ptr;	/* pointer to end of messages for a chunk */

	logical = get_logical_addr(p, start_buf, chunk_addr);

	/* increase chunk array size */
	if(oh->nchunks >= oh->alloc_nchunks) {

	    unsigned na = MAX(OBJ_NCHUNKS, oh->alloc_nchunks * 2);
	    OBJ_chunk_t *x = realloc(oh->chunk, (size_t)(na*sizeof(OBJ_chunk_t)));
	    if(!x) {
		error_push(ERR_INTERNAL, ERR_NONE_SEC, "Object Header:Internal allocation error", logical, NULL);
		CK_INC_ERR_DONE
	    }

	    oh->alloc_nchunks = na;
	    oh->chunk = x;
        }  /* end if */


        /* read the chunk raw data */
        chunkno = oh->nchunks++;
	if(chunkno == 0) {
	    oh->chunk[chunkno].addr = obj_head_addr;
	    if (format_objvers_two)
		oh->chunk[chunkno].size = chunk_size + OBJ_SIZEOF_HDR_VERS(OBJ_VERSION_2, oh);
	    else
		oh->chunk[chunkno].size = chunk_size + OBJ_SIZEOF_HDR_VERS(OBJ_VERSION_1, oh);
	} else {
	    oh->chunk[chunkno].addr = chunk_addr;
	    oh->chunk[chunkno].size = chunk_size;
	}

	if((oh->chunk[chunkno].image = calloc(1, oh->chunk[chunkno].size)) == NULL) {
	    error_push(ERR_INTERNAL, ERR_NONE_SEC, "Object Header:Internal allocation error", logical, NULL);
	    CK_INC_ERR_DONE
	}

	/* Handle chunk 0 as special case */
        if(chunkno == 0) {
            /* Check for speculative read of first chunk containing all the data needed */
            if(spec_read_size >= oh->chunk[0].size)
                memcpy(oh->chunk[0].image, buf, oh->chunk[0].size);
            else {
                /* Copy the object header prefix into chunk 0's image */
                memcpy(oh->chunk[0].image, buf, prefix_size);
		/* read the message data */
		if(FD_read(file, chunk_addr, (oh->chunk[0].size-prefix_size), (oh->chunk[chunkno].image+prefix_size)) == FAIL) {
		    error_push(ERR_FILE, ERR_NONE_SEC, "Object Header:Unable to read object header data", logical, NULL);
		    CK_INC_ERR_DONE
		}
            } /* end else */

            /* Point into chunk image to decode */
	    start_buf = oh->chunk[chunkno].image;
            p = oh->chunk[0].image + prefix_size;
        } /* end if */
        else {
	    if(FD_read(file, chunk_addr, chunk_size, oh->chunk[chunkno].image) == FAIL) {
		error_push(ERR_FILE, ERR_NONE_SEC, 
		    "Object Header:Unable to read object header data", logical, NULL);
		CK_INC_ERR_DONE
	    }
            /* Point into chunk image to decode */
	    start_buf = oh->chunk[chunkno].image;
            p = oh->chunk[chunkno].image;
        } /* end else */

	logical = get_logical_addr(p, start_buf, chunk_addr);

	/* Check for "CONT" magic # on chunks > 0 in later versions of the format */
	if(chunkno > 0 && format_objvers_two) {
            /* Magic number */
            if(memcmp(p, OBJ_CHK_MAGIC, (size_t)OBJ_SIZEOF_MAGIC)) {
		error_push(ERR_LEV_2, ERR_LEV_2A1b, 
		    "version 2 Object Header:Couldn't find CONT signature", logical, NULL);
		CK_INC_ERR_DONE
	    }
            p += OBJ_SIZEOF_MAGIC;
        } /* end if */

	/* Decode messages from this chunk */
	if(format_objvers_two)
	    eom_ptr = oh->chunk[chunkno].image + (oh->chunk[chunkno].size - OBJ_SIZEOF_CHKSUM_VERS(OBJ_VERSION_2));
	else
	    eom_ptr = oh->chunk[chunkno].image + (oh->chunk[chunkno].size - OBJ_SIZEOF_CHKSUM_VERS(OBJ_VERSION_1));

	while(p < eom_ptr) {
	    unsigned    mesgno;         /* Current message to operate on */
            size_t      mesg_size;      /* Size of message read in */
            unsigned    id;             /* ID (type) of current message */
            uint8_t     flags;          /* Flags for current message */
            OBJ_msg_crt_idx_t crt_idx = 0;  /* Creation index for current message */

	    /* Decode message prefix info */

	    logical = get_logical_addr(p, start_buf, chunk_addr);
            /* Version # */
            if(oh->version == OBJ_VERSION_1)
                UINT16DECODE(p, id)
            else
                id = *p++;

            if(id == OBJ_UNKNOWN_ID) {
		error_push(ERR_LEV_2, ERR_LEV_2A, 
		    "Object Header:unknown message ID encoded in file", logical, NULL);
		CK_INC_ERR
	    }

	    UINT16DECODE(p, mesg_size);
	    assert(mesg_size==OBJ_ALIGN_OH(oh, mesg_size));
	    flags = *p++;

	    if(flags & ~OBJ_MSG_FLAG_BITS) {
		error_push(ERR_LEV_2, ERR_LEV_2A, "Object Header:invalid message flag", logical, NULL);
		CK_INC_ERR
	    }

	    /* Reserved bytes/creation index */
	    if(!format_objvers_two)
                p += 3; /*reserved*/
            else {
                /* Only encode creation index if they are being tracked */
                if(oh->flags & OBJ_HDR_ATTR_CRT_ORDER_TRACKED)
                    UINT16DECODE(p, crt_idx);
            } /* end else */

	    /* Try to detect invalidly formatted object header message that
             *  extends past end of chunk.
             */
	    logical = get_logical_addr(p, start_buf, chunk_addr);
            if(p + mesg_size > eom_ptr) {
		error_push(ERR_LEV_2, ERR_LEV_2A, "Object Header:corrupt object header", logical, NULL);
		CK_INC_ERR_DONE
	    }
	    
	    if((!format_objvers_two) && (oh->nmesgs >= nmesgs)) {
		error_push(ERR_LEV_2, ERR_LEV_2A, "Object Header:corrupt object header", logical, NULL);
		CK_INC_ERR
	    }

            if(oh->nmesgs >= oh->alloc_nmesgs)
		if(OBJ_alloc_msgs(oh, (size_t)1) < 0)
		    CK_INC_ERR_DONE

            /* Get index for message */
	    mesgno = oh->nmesgs++;
            oh->mesg[mesgno].flags = flags;
	    oh->mesg[mesgno].native = NULL;
	    oh->mesg[mesgno].raw = (uint8_t *)p;
	    oh->mesg[mesgno].raw_size = mesg_size;
	    oh->mesg[mesgno].chunkno = chunkno;

	    if(id >= NELMTS(message_type_g) || NULL == message_type_g[id]) {
		oh->mesg[mesgno].type = message_type_g[OBJ_UNKNOWN_ID];
	    } else
		oh->mesg[mesgno].type = message_type_g[id];

	    p+= mesg_size;

	    logical = get_logical_addr(p, start_buf, chunk_addr);
	    /* Check for 'gap' at end of chunk */
	    if((eom_ptr - p) > 0 && 
	       (eom_ptr - p) < OBJ_SIZEOF_MSGHDR_VERS(OBJ_VERSION_2,oh->flags&OBJ_HDR_ATTR_CRT_ORDER_TRACKED)) {

                /* Gaps can only occur in later versions of the format */
		if(format_objvers_two)
		    p += eom_ptr - p;
		else {
		    error_push(ERR_LEV_2, ERR_LEV_2A, "Object Header:corrupt object header", logical, NULL);
		    CK_INC_ERR_DONE
		}
            } /* end if */
	} /* end while */

	if(format_objvers_two) {
            uint32_t stored_chksum;     /* Checksum from file */
            uint32_t computed_chksum;   /* Checksum computed in memory */

	    computed_chksum = checksum_metadata(oh->chunk[chunkno].image, (oh->chunk[chunkno].size-OBJ_SIZEOF_CHKSUM), 0);
	    logical = get_logical_addr(p, start_buf, chunk_addr);
            UINT32DECODE(p, stored_chksum);

	    if(computed_chksum != stored_chksum) {
		error_push(ERR_LEV_2, ERR_LEV_2A1b, "version 2 Object Header:Bad checksum", logical, NULL);
		CK_INC_ERR
	    }
        } /* end if */

	assert(p == oh->chunk[chunkno].image + oh->chunk[chunkno].size);

       	/* decode next object header continuation message */
        for(chunk_addr = CK_ADDR_UNDEF; !addr_defined(chunk_addr) && curmesg < oh->nmesgs; ++curmesg) {
	    if(oh->mesg[curmesg].type->id == OBJ_CONT_ID) {
		start_buf = oh->chunk[oh->mesg[curmesg].chunkno].image;
		logi_base = oh->chunk[oh->mesg[curmesg].chunkno].addr;
		cont = (OBJ_CONT->decode) (file, oh->mesg[curmesg].raw, start_buf, logi_base);
		logical = get_logical_addr(oh->mesg[curmesg].raw, start_buf, logi_base);
		if (cont == NULL) {
		    error_push(ERR_LEV_2, ERR_LEV_2A, 
			"Object Header:Corrupt continuation message...skipped", logical, NULL);
		    CK_INC_ERR_CONTINUE
		}
                oh->mesg[curmesg].native = cont;
                chunk_addr = cont->addr;
                chunk_size = cont->size;
                cont->chunkno = oh->nchunks;    /*the next chunk to allocate */
	    }  /* end if */
        }  /* end for */
    }  /* end while(addr_defined(chunk_addr)) */

done:
    if(ret_err && !object_api()) {
	error_print(stderr, file);
	error_clear();
    }

    if(oh && ret_oh)
	*ret_oh = oh;

    if(oh && oh->nmesgs)
	if(decode_validate_messages(file, oh) < 0)
	    ++ret_other_err;

    if(ret_err || ret_other_err)
	ret_value = FAIL;

    return(ret_value);
} /* check_obj_header() */


/* 
 * Checksum: copy from the library, see description there
 */
#define lookup3_rot(x,k) (((x)<<(k)) ^ ((x)>>(32-(k))))

#define lookup3_mix(a,b,c) \
{ \
  a -= c;  a ^= lookup3_rot(c, 4);  c += b; \
  b -= a;  b ^= lookup3_rot(a, 6);  a += c; \
  c -= b;  c ^= lookup3_rot(b, 8);  b += a; \
  a -= c;  a ^= lookup3_rot(c,16);  c += b; \
  b -= a;  b ^= lookup3_rot(a,19);  a += c; \
  c -= b;  c ^= lookup3_rot(b, 4);  b += a; \
}

#define lookup3_final(a,b,c) \
{ \
  c ^= b; c -= lookup3_rot(b,14); \
  a ^= c; a -= lookup3_rot(c,11); \
  b ^= a; b -= lookup3_rot(a,25); \
  c ^= b; c -= lookup3_rot(b,16); \
  a ^= c; a -= lookup3_rot(c,4);  \
  b ^= a; b -= lookup3_rot(a,14); \
  c ^= b; c -= lookup3_rot(b,24); \
}

uint32_t
checksum_lookup3(const void *key, ck_size_t length, uint32_t initval)
{
    const uint8_t *k = (const uint8_t *)key;
    uint32_t a, b, c;           /* internal state */

    assert(key);
    assert(length > 0);

    /* Set up the internal state */
    a = b = c = 0xdeadbeef + ((uint32_t)length) + initval;

    /*--------------- all but the last block: affect some 32 bits of (a,b,c) */
    while (length > 12)
    {
      a += k[0];
      a += ((uint32_t)k[1])<<8;
      a += ((uint32_t)k[2])<<16;
      a += ((uint32_t)k[3])<<24;
      b += k[4];
      b += ((uint32_t)k[5])<<8;
      b += ((uint32_t)k[6])<<16;
      b += ((uint32_t)k[7])<<24;
      c += k[8];
      c += ((uint32_t)k[9])<<8;
      c += ((uint32_t)k[10])<<16;
      c += ((uint32_t)k[11])<<24;
      lookup3_mix(a, b, c);
      length -= 12;
      k += 12;
    }

    /*-------------------------------- last block: affect all 32 bits of (c) */
    switch(length)                   /* all the case statements fall through */
    {
        case 12: c+=((uint32_t)k[11])<<24;
        case 11: c+=((uint32_t)k[10])<<16;
        case 10: c+=((uint32_t)k[9])<<8;
        case 9 : c+=k[8];
        case 8 : b+=((uint32_t)k[7])<<24;
        case 7 : b+=((uint32_t)k[6])<<16;
        case 6 : b+=((uint32_t)k[5])<<8;
        case 5 : b+=k[4];
        case 4 : a+=((uint32_t)k[3])<<24;
        case 3 : a+=((uint32_t)k[2])<<16;
        case 2 : a+=((uint32_t)k[1])<<8;
        case 1 : a+=k[0];
                 break;
        case 0 : goto done;
    }

    lookup3_final(a, b, c);

done:
    return(c);
} /* checksum_lookup3() */

uint32_t
checksum_metadata(const void *data, ck_size_t len, uint32_t initval)
{
    assert(data);
    assert(len > 0);

    return(checksum_lookup3(data, len, initval));
}

/* end checksum routines from the library */


void
print_version(const char *prog_name)
{
    fprintf(stdout, "%s: Version %s\n", prog_name, H5Check_VERSION);
} /* print_version() */

void
usage(char *prog_name)
{
    fflush(stdout);
    fprintf(stdout, "usage: %s [OPTIONS] file\n", prog_name);
    fprintf(stdout, "  OPTIONS\n");
    fprintf(stdout, "     -h,  --help   	Print a usage message and exit\n");
    fprintf(stdout, "     -V,  --version	Print version number and exit\n");
    fprintf(stdout, "     -vn, --verbose=n	Verboseness mode\n");
    fprintf(stdout, "     		n=0	Terse--only indicates if the file is compliant or not\n");
    fprintf(stdout, "     		n=1	Default--prints its progress and all errors found\n");
    fprintf(stdout, "     		n=2	Verbose--prints everything it knows, usually for debugging\n");
    fprintf(stdout, "     -fn, --format=n	File Format\n");
    fprintf(stdout, "     		n=16	Validate according to library release version 1.6.6\n");
    fprintf(stdout, "     		n=18	Default--Validate according to library release version 1.8.0\n");
    fprintf(stdout, "     -oa, --object=a	Check object header\n");
    fprintf(stdout, "     		a	Address of the object header to be validated\n");
    fprintf(stdout, "\n");

} /* usage() */

void
leave(int ret)
{
    exit(ret);
} /* leave() */

int
debug_verbose(void)
{
   return(g_verbose_num==DEBUG_VERBOSE);
} /* debug_verbose() */

int 
object_api(void)
{
    if(g_obj_api)
	g_obj_api_err++;
    return(g_obj_api);
} /* object_api() */
