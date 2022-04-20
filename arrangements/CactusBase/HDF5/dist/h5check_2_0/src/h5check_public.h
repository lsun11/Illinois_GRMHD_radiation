typedef unsigned long long              ck_addr_t;
typedef unsigned int                    ck_bool_t;
typedef int                             ck_err_t;

#define CK_ADDR_UNDEF           	((ck_addr_t)(-1))
#define	NSLOTS 	32

typedef struct errmsg_t {
	const	char	*desc;
	ck_addr_t	addr;
}errmsg_t;

typedef struct ck_errmsg_t {
	int	nused;
	errmsg_t	slot[NSLOTS];
} ck_errmsg_t;

ck_err_t	h5checker_obj(char *, ck_addr_t, int, ck_errmsg_t *);
