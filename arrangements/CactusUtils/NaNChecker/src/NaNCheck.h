 /*@@
   @header    NaNCheck.h
   @date      Mon 10 Dec 2001
   @author    Thomas Radke
   @desc
              Prototypes for NaNChecker API routines
   @enddesc
   @version   $Header: /cactusdevcvs/CactusUtils/NaNChecker/src/NaNCheck.h,v 1.3 2003/04/07 20:20:12 tradke Exp $
 @@*/


#ifndef NANCHECKER_NANCHECK_H
#define NANCHECKER_NANCHECK_H 1

#ifdef  __cplusplus
extern "C"
{
#endif

/* prototypes of exported functions */
int NaNChecker_CheckVarsForNaN (const cGH *GH,
                                int report_max,
                                const char *vars,
                                const char *check_for,
                                const char *action_if_found);
int NaNChecker_SetVarsToNaN (const cGH *GH,
                             const char *vars);

#ifdef  __cplusplus
}
#endif

#endif /* NANCHECKER_NANCHECK_H */
