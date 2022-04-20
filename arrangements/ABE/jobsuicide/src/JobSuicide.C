/*
  This subroutine is useful for queuing systems.
*/

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <fstream>

#include <time.h>
#include <math.h>

#include <unistd.h>

#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

static char *rcsid="$Header: /peter/piper/picked/a/peck/of/pickled/whatever $";

CCTK_FILEVERSION(JobSuicideCheck)

  extern "C" void JobSuicideCheck(CCTK_ARGUMENTS)
{
  using namespace std;
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT handle;

  int kill_run_flag = 0;
  if (cctk_iteration%20 == 0) {
    if (CCTK_MyProc(cctkGH) == 0) {
      ifstream time_file("start_time");
      double new_start_time = 0,new_end_time = 0;
      double end_time = *start_time + run_time * 3600.0;
      time_file >> new_start_time;
      time_file >> new_end_time;
      time_file.close();

      printf("Hello! rt = %e new_start_time = %.16f start_time = %.16f \n Hi. new_end_time = %.16f end_time = %.16f \n",run_time,new_start_time,*start_time,new_end_time,end_time);

      //Following lines handle the case when another job designated to run for shorter time writes to file after we do.
      //This case is possible when two jobs start within one second of each other, but shorter job writes last!
      if(end_time > new_end_time) {
	ofstream new_time_file("start_time");
	new_time_file << setprecision(15);
	new_time_file << *start_time << endl;
	new_time_file << end_time << endl;
	new_time_file.close();
      }

      if (new_start_time > *start_time || new_end_time > end_time) {
	cout << "\n\n\nO happy dagger! \nThis is thy sheath; there rust, and let me die.\n\n"
	     << "\t-Romeo and Juliet, Act V\n\n\n"
	     << "2Looks like another job (that will finish after me) has started. I'd better quit." << endl
	     << "To disable the suicide feature, please set the \njobsuicide::run_time parameter in your .par file to some NEGATIVE integer"
	     << endl;
	kill_run_flag = 1;
      } // end if (new_start_time > *start_time)
      if (kill_run_flag > 0) {
	CCTK_Abort(cctkGH, EXIT_FAILURE);  
	exit(1);
      }
    } // end if (CCTK_MyProc(cctkGH) == 0)

    /*    
	  int sum_kill_run_flag = 0;
	  handle = CCTK_ReductionHandle("sum");
	  CCTK_ReduceLocalScalar(cctkGH,-1,handle,&kill_run_flag,&sum_kill_run_flag,CCTK_VARIABLE_INT);
	  
	  if (sum_kill_run_flag > 0) {
	  CCTK_Abort(cctkGH, EXIT_FAILURE);  
	  exit(1);
	  }
    */
  }
  return;
}

  extern "C" void Setup_JobSuicide(CCTK_ARGUMENTS)
{
  using namespace std;
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  int kill_run_flag = 0;
  CCTK_INT handle;

  double local_start_time = 0.0;
  *start_time = 0.0;

  if (CCTK_MyProc(cctkGH) == 0) {
    *start_time = (double)time(NULL);
    double end_time = *start_time + run_time * 3600.0;
    ifstream time_file("start_time");
    double old_start_time = 0, old_end_time = 0;
    if (!time_file) {
      // do nothing
    } else {
      time_file >> old_start_time;
      time_file >> old_end_time; 
    }
    time_file.close();

    printf("Hello! rt = %e end_time = %.16f old_end_time = %.16f start_time = %.16f \n",run_time,end_time,old_end_time,*start_time);

    if ((end_time > old_end_time + 900) || *start_time > end_time) {	
      ofstream new_time_file("start_time");
      new_time_file << setprecision(15);
      new_time_file << *start_time << endl;
      new_time_file << end_time << endl;
      new_time_file.close();
    } else {
	cout << "\n\n\nO happy dagger! \nThis is thy sheath; there rust, and let me die.\n\n"
	     << "\t-Romeo and Juliet, Act V\n\n\n"
	     << "2Looks like another job (that will finish after me) has started. I'd better quit." << endl
	     << "To disable the suicide feature, please set the \njobsuicide::run_time parameter in your .par file to some NEGATIVE integer"
	   << endl;
      kill_run_flag = 1;
    }
    if (kill_run_flag > 0) {
      CCTK_Abort(cctkGH, EXIT_FAILURE);  
      exit(1);
    }
    local_start_time = *start_time;
  }

  CCTK_ReduceLocalScalar(cctkGH,-1,CCTK_ReductionHandle("sum"),&local_start_time,start_time,CCTK_VARIABLE_REAL);
  
  /*
    int sum_kill_run_flag = 0;
    handle = CCTK_ReductionHandle("sum");
    CCTK_ReduceLocalScalar(cctkGH,-1,handle,&kill_run_flag,&sum_kill_run_flag,CCTK_VARIABLE_INT);
    if (sum_kill_run_flag > 0) {
    CCTK_Abort(cctkGH, EXIT_FAILURE);  
    exit(1);
    }
  */
  return;
}
