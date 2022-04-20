/*@@
  @file      ApplyCartoon.c
  @date      6 Jan 2003
  @author    David Rideout
  @desc 
  Applies the Cartoon boundary condition to all variables 
  selected for a bc.
  @enddesc 
  @version   $Header: /cactusdevcvs/BetaThorns/Cartoon2D/src/ApplyCartoon.c,v 1.12 2008/02/27 05:20:47 schnetter Exp $
  @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"
#include "Cartoon2D.h"
#include "Cartoon2D_tensors.h"
#include "cctk_Groups.h"
#include <stdlib.h>
#include <stdio.h>


static const char *rcsid = "$Header: /cactusdevcvs/BetaThorns/Cartoon2D/src/ApplyCartoon.c,v 1.12 2008/02/27 05:20:47 schnetter Exp $";

CCTK_FILEVERSION(BetaThorns_Cartoon2D_ApplyCartoon_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/
/* #define DEBUG 1 */
#define TENSORTYPE_BUFF_SIZE 7
#define PROLONG_BUFF_SIZE 1000

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Aliased Routine Prototypes ***********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void Cartoon_ApplyBoundaries(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     Aliased Routines   **********************
 ********************************************************************/
 
/********************************************************************
 *********************     Scheduled Routines   **********************
 ********************************************************************/

/*@@
  @routine    Cartoon_ApplyBoundaries
  @date       6 Jan 2003
  @author     David Rideout
  @desc  
  This will apply the Cartoon boundary condition to all
  variables selected for any (physical) boundary condition.
  @enddesc 
  @calls      
  @history 
  @endhistory 
  @var        CCTK_ARGUMENTS
  @vdesc      Cactus argument list
  @vtype      CCTK_*
  @vio        in
  @endvar
  @returntype void
  @returndesc
  @endreturndesc
  @@*/

void Cartoon_ApplyBoundaries(CCTK_ARGUMENTS) 
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int num_vars, err, i, gi, group_tags_table;
  CCTK_INT * vars;
  char tensortype[TENSORTYPE_BUFF_SIZE];
  char prolongtype[PROLONG_BUFF_SIZE];
  int prolongmethod;

  int whiskycartoon, len;

  


  /* Check grid size */
  if(cctk_gsh[1] != 2*cctk_nghostzones[1]+1)
    {
      CCTK_WARN(0, "Grid size in y direction inappropriate.");
    }

  /* Allocate memory to hold selected bcs */
  num_vars = Boundary_SelectedGVs(cctkGH, 0, NULL, NULL, NULL, NULL, NULL);
#ifdef DEBUG
  printf("Cartoon_ApplyBoundaries: num_vars is %d\n",num_vars);
#endif
  vars = malloc(num_vars*sizeof *vars);

  /* get all selected vars */
  err = Boundary_SelectedGVs(cctkGH, num_vars, vars, NULL, NULL, NULL, NULL);
  if (err != num_vars)
    {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		 "Boundary_SelectedGVs returned %d selected variables, but %d "
		 "expected", err, num_vars);
    }

  /* Apply CartoonBC to each of them. */
  /* One should probably check to see that the entire group has been
   * selected, since Cartoon operates on an entire group at a time.
   * For now I'll just skip a variable if it is not a 'group leader'.
   */
  for (i=0; i<num_vars; ++i) {
#ifdef DEBUG
    printf("Cartoon_ApplyBoundaries: i=%d applying cartoon to vi %d\n",i,
	   vars[i]);
#endif
    gi = CCTK_GroupIndexFromVarI(vars[i]);
    if (gi<0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		 "Invalid variable index %d selected for a boundary "
		 "condition", gi);
    }
    if (vars[i] != CCTK_FirstVarIndexI(gi)) 
      {
	/* not a group leader -- skip to next var index */
	continue;
      }

    /* Here one should check that the entire group is registered,
     * using CCTK_NumVarsInGroupI.
     */

    /* Get table handle for group tags table */
    group_tags_table = CCTK_GroupTagsTableI(gi);
    if (group_tags_table<0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		 "Tags table for variable group %s is %d", CCTK_GroupName(gi),
		 group_tags_table);
    } 
   

    /* This lines are writen for Whisky2D to avoid to treat the hydrovariables with Cartoon */
 
    len = Util_TableGetString(group_tags_table,0,NULL,"whiskycartoon");

    whiskycartoon = 1;          /* default is yes */
    if (len >= 0)
      {
        char* value = malloc (len + 1);
        Util_TableGetString (group_tags_table, len + 1, value, "whiskycartoon");
        if (CCTK_Equals(value, "yes"))
          {
            whiskycartoon = 1;
          }
        else if (CCTK_Equals(value, "no"))
          {
            whiskycartoon = 0;
          }
        else
          {
            CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Tags table entry \"whiskycartoon\" for variable group %s should be either \"yes\" or \"no\", but is \"%s\"",
                       CCTK_GroupName(gi), value);
          }
        free (value);
      }

    if (whiskycartoon)
      {
  
	/*###################################################################*/


	/* Get tensor type from group tags table */
	err = Util_TableGetString(group_tags_table, TENSORTYPE_BUFF_SIZE,
				  tensortype, "tensortypealias");
	if (err<0) 
	  {
	    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		       "Error in TAGS table for variable group %s",
		       CCTK_GroupName(gi));
	  }
#ifdef DEBUG
	printf("Cartoon_ApplyBoundaries: tensor type is %s\n",tensortype);
#endif

	/* Get prolongation type from group tags table */

	err = Util_TableGetString(group_tags_table, PROLONG_BUFF_SIZE,
				  prolongtype, "Prolongation");

	if (err == UTIL_ERROR_TABLE_NO_SUCH_KEY) 
	  {
	    /* Use the default */

	    prolongmethod = PROLONG_LAGRANGE;
	  }
	else if (err < 0)
	  {
	    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
		       "Error (%d) in TAGS table for variable group %s", err,
		       CCTK_GroupName(gi));
	  }
	else
	  {
	    if (CCTK_Equals(prolongtype, "None")) 
	      {
		prolongmethod = PROLONG_NONE; /* But why? */
	      } 
	    else if (CCTK_Equals(prolongtype, "Lagrange")) 
	      {
		prolongmethod = PROLONG_LAGRANGE;
	      } 
	    else if (CCTK_Equals(prolongtype, "TVD")) 
	      {
		prolongmethod = PROLONG_ENO;
	      } 
	    else if (CCTK_Equals(prolongtype, "ENO")) 
	      {
		prolongmethod = PROLONG_ENO;
	      } 
	    else 
	      {
		CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
			   "Error in TAGS table for variable group %s",
			   CCTK_GroupName(gi));
	      }
	  }    
    
	/* Here one should check that the group's size is correct for the
	 * specified tensortype.
	 */

	/* Call BndCartoon2DVI, passing the appropriate tensor type integer 
	   macro */
	if (CCTK_Equals(tensortype, "scalar"))
	  {
	    BndCartoon2DVI(cctkGH, TENSORTYPE_SCALAR, prolongmethod, vars[i]);
	  } else if ((CCTK_Equals(tensortype, "u")) ||
		     (CCTK_Equals(tensortype, "d")))
	  {
	    BndCartoon2DVI(cctkGH, TENSORTYPE_U, prolongmethod, vars[i]);
	  } else if (CCTK_Equals(tensortype, "dd_sym"))
	  {
	    BndCartoon2DVI(cctkGH, TENSORTYPE_DDSYM, prolongmethod, vars[i]);
	  } else 
	  {
	    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
		       "invalid tensor type for group %s", CCTK_GroupName(gi));
	  }  
      }
  }
  /* Free data */
  free(vars);
}
