/*@@
   @file    CarpetRegrid_Parameters.c
   @author  Automatically generated by CreateParameterBindings.pl
   @desc
            Creates/extends parameters for this thorn
   @enddesc
 @@*/


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "cctk_Config.h"
#include "cctk_Constants.h"
#include "ParameterBindings.h"
#include "CParameterStructNames.h"

/* structure containing all private parameters of thorn CarpetRegrid */
struct
{
  CCTK_REAL l1xmax;
  CCTK_REAL l1xmin;
  CCTK_REAL l1ymax;
  CCTK_REAL l1ymin;
  CCTK_REAL l1zmax;
  CCTK_REAL l1zmin;
  CCTK_REAL l2xmax;
  CCTK_REAL l2xmin;
  CCTK_REAL l2ymax;
  CCTK_REAL l2ymin;
  CCTK_REAL l2zmax;
  CCTK_REAL l2zmin;
  CCTK_REAL l3xmax;
  CCTK_REAL l3xmin;
  CCTK_REAL l3ymax;
  CCTK_REAL l3ymin;
  CCTK_REAL l3zmax;
  CCTK_REAL l3zmin;
  CCTK_REAL maxerror;
  CCTK_REAL minfraction;
  CCTK_REAL moving_centre_x;
  CCTK_REAL moving_centre_y;
  CCTK_REAL moving_centre_z;
  CCTK_REAL moving_circle_frequency;
  CCTK_REAL moving_circle_radius;
  CCTK_REAL moving_region_radius;
  const char * activate_levels_on_regrid;
  const char * coordinates;
  const char * errorvar;
  const char * gridpoints;
  const char * moving_trajectory;
  const char * outerbounds;
  const char * refined_regions;
  CCTK_INT activate_next;
  CCTK_INT keep_same_grid_structure;
  CCTK_INT l1ixmax;
  CCTK_INT l1ixmin;
  CCTK_INT l1iymax;
  CCTK_INT l1iymin;
  CCTK_INT l1izmax;
  CCTK_INT l1izmin;
  CCTK_INT l2ixmax;
  CCTK_INT l2ixmin;
  CCTK_INT l2iymax;
  CCTK_INT l2iymin;
  CCTK_INT l2izmax;
  CCTK_INT l2izmin;
  CCTK_INT l3ixmax;
  CCTK_INT l3ixmin;
  CCTK_INT l3iymax;
  CCTK_INT l3iymin;
  CCTK_INT l3izmax;
  CCTK_INT l3izmin;
  CCTK_INT merge_overlapping_components;
  CCTK_INT minwidth;
  CCTK_INT num_new_levels;
  CCTK_INT refinement_levels;
  CCTK_INT regrid_every;
  CCTK_INT smart_outer_boundaries;
  CCTK_INT symmetry_x;
  CCTK_INT symmetry_y;
  CCTK_INT symmetry_z;
  CCTK_INT tracking;
  CCTK_INT verbose;
  CCTK_INT veryverbose;
} PRIVATE_CARPETREGRID_STRUCT;


/* structure containing all restricted parameters of thorn CarpetRegrid */
struct
{
  CCTK_REAL offsetx[10];
  CCTK_REAL offsety[10];
  CCTK_REAL offsetz[10];
  CCTK_INT num_offsets;
  CCTK_INT offset_firstlevel;
} RESTRICTED_CARPETREGRID_STRUCT;


int CCTKi_BindingsCreateCarpetRegridParameters(void);
int CCTKi_BindingsCreateCarpetRegridParameters(void)
{
  CCTKi_ParameterCreate("num_offsets",
                        "CarpetRegrid",
                        "INT",
                        "RESTRICTED",
                        CCTK_STEERABLE_ALWAYS,
                        "Number of given offsets",
                        "0",
                        &(RESTRICTED_CARPETREGRID_STRUCT.num_offsets),
                        0,
                        NULL,
                        1,
                        "0:10", "");

  CCTKi_ParameterCreate("offset_firstlevel",
                        "CarpetRegrid",
                        "INT",
                        "RESTRICTED",
                        CCTK_STEERABLE_ALWAYS,
                        "First (lowest) refinement level that should have an offset applied",
                        "1",
                        &(RESTRICTED_CARPETREGRID_STRUCT.offset_firstlevel),
                        0,
                        NULL,
                        1,
                        "1:*", "");

  CCTKi_ParameterCreate("offsetx",
                        "CarpetRegrid",
                        "REAL",
                        "RESTRICTED",
                        CCTK_STEERABLE_ALWAYS,
                        "x-coordinate of offset",
                        "0.0",
                        (RESTRICTED_CARPETREGRID_STRUCT.offsetx),
                        10,
                        NULL,
                        1,
                        "(*:*)", "");

  CCTKi_ParameterCreate("offsety",
                        "CarpetRegrid",
                        "REAL",
                        "RESTRICTED",
                        CCTK_STEERABLE_ALWAYS,
                        "y-coordinate of offset",
                        "0.0",
                        (RESTRICTED_CARPETREGRID_STRUCT.offsety),
                        10,
                        NULL,
                        1,
                        "(*:*)", "");

  CCTKi_ParameterCreate("offsetz",
                        "CarpetRegrid",
                        "REAL",
                        "RESTRICTED",
                        CCTK_STEERABLE_ALWAYS,
                        "z-coordinate of offset",
                        "0.0",
                        (RESTRICTED_CARPETREGRID_STRUCT.offsetz),
                        10,
                        NULL,
                        1,
                        "(*:*)", "");

  CCTKi_ParameterCreate("activate_levels_on_regrid",
                        "CarpetRegrid",
                        "KEYWORD",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Whether to activate or deactivate new levels on regridding",
                        "none",
                        &(PRIVATE_CARPETREGRID_STRUCT.activate_levels_on_regrid),
                        0,
                        NULL,
                        3,
                        "none", "Do not activate or deactivate any levels",
                        "fixed", "Activate or deactivate a fixed number of levels",
                        "function", "Activate or deactivate a variable number of levels, determined by a user-specified function.  When this option is used, the parameters num_new_levels and activate_next have no effect and should not be set.");

  CCTKi_ParameterCreate("activate_next",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "The next iteration at which new levels should be activated",
                        "1",
                        &(PRIVATE_CARPETREGRID_STRUCT.activate_next),
                        0,
                        NULL,
                        1,
                        "0:", "Note that this parameter is steered when new levels are activated");

  CCTKi_ParameterCreate("coordinates",
                        "CarpetRegrid",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "List of bounding box coordinates",
                        "",
                        &(PRIVATE_CARPETREGRID_STRUCT.coordinates),
                        0,
                        NULL,
                        2,
                        "^$", "leave empty for no refinement",
                        ".*", "[ [ ([<xmin>,<ymin>,<zmin>]:[<xmax>,<ymax>,<zmax>]:[<xstride>,<ystride>,<zstride>]), ... ], ... ]");

  CCTKi_ParameterCreate("errorvar",
                        "CarpetRegrid",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Name of grid function that contains the error",
                        "",
                        &(PRIVATE_CARPETREGRID_STRUCT.errorvar),
                        0,
                        NULL,
                        1,
                        ".*", "must be the name of a grid function");

  CCTKi_ParameterCreate("gridpoints",
                        "CarpetRegrid",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "List of bounding box gridpoints",
                        "",
                        &(PRIVATE_CARPETREGRID_STRUCT.gridpoints),
                        0,
                        NULL,
                        2,
                        "^$", "leave empty for no refinement",
                        ".*", "[ [ ([<imin>,<jmin>,<kmin>]:[<imax>,<jmax>,<kmax>]:[<istride>,<jstride>,<kstride>]), ... ], ... ]");

  CCTKi_ParameterCreate("keep_same_grid_structure",
                        "CarpetRegrid",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Do not allow the grid structure to change; only allow levels to be switched on or off",
                        "no",
                        &(PRIVATE_CARPETREGRID_STRUCT.keep_same_grid_structure),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("l1ixmax",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 1 box in x-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1ixmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l1ixmin",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 1 box in x-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1ixmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l1iymax",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 1 box in y-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1iymax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l1iymin",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 1 box in y-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1iymin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l1izmax",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 1 box in z-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1izmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l1izmin",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 1 box in z-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1izmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l1xmax",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 1 box in x-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1xmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l1xmin",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 1 box in x-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1xmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l1ymax",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 1 box in y-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1ymax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l1ymin",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 1 box in y-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1ymin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l1zmax",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 1 box in z-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1zmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l1zmin",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 1 box in z-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l1zmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2ixmax",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 2 box in x-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2ixmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2ixmin",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 2 box in x-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2ixmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2iymax",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 2 box in y-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2iymax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2iymin",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 2 box in y-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2iymin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2izmax",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 2 box in z-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2izmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2izmin",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 2 box in z-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2izmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2xmax",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 2 box in x-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2xmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2xmin",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 2 box in x-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2xmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2ymax",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 2 box in y-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2ymax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2ymin",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 2 box in y-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2ymin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2zmax",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 2 box in z-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2zmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l2zmin",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 2 box in z-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l2zmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3ixmax",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 3 box in x-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3ixmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3ixmin",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 3 box in x-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3ixmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3iymax",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 3 box in y-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3iymax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3iymin",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 3 box in y-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3iymin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3izmax",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 3 box in z-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3izmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3izmin",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 3 box in z-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3izmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3xmax",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 3 box in x-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3xmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3xmin",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 3 box in x-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3xmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3ymax",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 3 box in y-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3ymax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3ymin",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 3 box in y-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3ymin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3zmax",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Upper boundary of level 3 box in z-direction",
                        "-1",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3zmax),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("l3zmin",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Lower boundary of level 3 box in z-direction",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.l3zmin),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("maxerror",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Maximum allowed error for non-refined grid points",
                        "1.0",
                        &(PRIVATE_CARPETREGRID_STRUCT.maxerror),
                        0,
                        NULL,
                        1,
                        "*:*", "everything goes");

  CCTKi_ParameterCreate("merge_overlapping_components",
                        "CarpetRegrid",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Merge overlapping components",
                        "no",
                        &(PRIVATE_CARPETREGRID_STRUCT.merge_overlapping_components),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("minfraction",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Minimum fraction of points in need of refinement in a refined region",
                        "0.75",
                        &(PRIVATE_CARPETREGRID_STRUCT.minfraction),
                        0,
                        NULL,
                        1,
                        "0:1", "must be positive and less than one");

  CCTKi_ParameterCreate("minwidth",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Minimum width of refined region",
                        "8",
                        &(PRIVATE_CARPETREGRID_STRUCT.minwidth),
                        0,
                        NULL,
                        1,
                        "1:*", "must be positive");

  CCTKi_ParameterCreate("moving_centre_x",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "x-coordinate of the centre",
                        "0.0",
                        &(PRIVATE_CARPETREGRID_STRUCT.moving_centre_x),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("moving_centre_y",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "y-coordinate of the centre",
                        "0.0",
                        &(PRIVATE_CARPETREGRID_STRUCT.moving_centre_y),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("moving_centre_z",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "z-coordinate of the centre",
                        "0.0",
                        &(PRIVATE_CARPETREGRID_STRUCT.moving_centre_z),
                        0,
                        NULL,
                        1,
                        ":", "");

  CCTKi_ParameterCreate("moving_circle_frequency",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Angular frequency on the circle",
                        "1.0",
                        &(PRIVATE_CARPETREGRID_STRUCT.moving_circle_frequency),
                        0,
                        NULL,
                        1,
                        "0:", "");

  CCTKi_ParameterCreate("moving_circle_radius",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Radius of the circle",
                        "1.0",
                        &(PRIVATE_CARPETREGRID_STRUCT.moving_circle_radius),
                        0,
                        NULL,
                        1,
                        "0:", "");

  CCTKi_ParameterCreate("moving_region_radius",
                        "CarpetRegrid",
                        "REAL",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Radius of the moving region (on the first refined level)",
                        "1.0",
                        &(PRIVATE_CARPETREGRID_STRUCT.moving_region_radius),
                        0,
                        NULL,
                        1,
                        "(0:", "");

  CCTKi_ParameterCreate("moving_trajectory",
                        "CarpetRegrid",
                        "KEYWORD",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Type of trajectory",
                        "point",
                        &(PRIVATE_CARPETREGRID_STRUCT.moving_trajectory),
                        0,
                        NULL,
                        2,
                        "point", "Do not move",
                        "circle", "Move in a circle");

  CCTKi_ParameterCreate("num_new_levels",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "When regridding, activate this many new levels (if possible).  Note that this will steer the parameter refinement_levels.",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.num_new_levels),
                        0,
                        NULL,
                        1,
                        ":", "Number of new levels to activate (negative numbers deactivate)");

  CCTKi_ParameterCreate("outerbounds",
                        "CarpetRegrid",
                        "STRING",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Outer boundaries",
                        "",
                        &(PRIVATE_CARPETREGRID_STRUCT.outerbounds),
                        0,
                        NULL,
                        2,
                        "^$", "leave empty for no outer boundaries",
                        ".*", "[ [ [[?,?],[?,?],[?,?]], ... ], ...]");

  CCTKi_ParameterCreate("refined_regions",
                        "CarpetRegrid",
                        "KEYWORD",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Regions where the grid is refined",
                        "centre",
                        &(PRIVATE_CARPETREGRID_STRUCT.refined_regions),
                        0,
                        NULL,
                        8,
                        "none", "Don't refine",
                        "centre", "Refine around the centre of the grid only",
                        "manual-gridpoints", "Refine the regions specified by integer grid points l[123]i[xyz]{min,max}",
                        "manual-coordinates", "Refine the regions specified by coordinates l[123][xyz]{min,max}",
                        "manual-gridpoint-list", "Refine the regions specified by integer grid points in the parameter 'gridpoints'",
                        "manual-coordinate-list", "Refine the regions specified by coordinates in the parameter 'coordinates'",
                        "moving", "Refine a moving region",
                        "automatic", "Refine automatically");

  CCTKi_ParameterCreate("refinement_levels",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Number of refinement levels (including the base level)",
                        "1",
                        &(PRIVATE_CARPETREGRID_STRUCT.refinement_levels),
                        0,
                        NULL,
                        1,
                        "1:*", "must be positive, and must not be larger than Carpet::max_refinement_levels");

  CCTKi_ParameterCreate("regrid_every",
                        "CarpetRegrid",
                        "INT",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Regrid every n time steps",
                        "0",
                        &(PRIVATE_CARPETREGRID_STRUCT.regrid_every),
                        0,
                        NULL,
                        3,
                        "-1", "regrid never",
                        "0", "regrid during initial data calculation only",
                        "1:*", "regrid every n time steps");

  CCTKi_ParameterCreate("smart_outer_boundaries",
                        "CarpetRegrid",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Use the CoordBase interface for outer boundaries",
                        "no",
                        &(PRIVATE_CARPETREGRID_STRUCT.smart_outer_boundaries),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("symmetry_x",
                        "CarpetRegrid",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Refine the lower half in x-direction",
                        "no",
                        &(PRIVATE_CARPETREGRID_STRUCT.symmetry_x),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("symmetry_y",
                        "CarpetRegrid",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Refine the lower half in y-direction",
                        "no",
                        &(PRIVATE_CARPETREGRID_STRUCT.symmetry_y),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("symmetry_z",
                        "CarpetRegrid",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Refine the lower half in z-direction",
                        "no",
                        &(PRIVATE_CARPETREGRID_STRUCT.symmetry_z),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("tracking",
                        "CarpetRegrid",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_ALWAYS,
                        "Enable tracking",
                        "no",
                        &(PRIVATE_CARPETREGRID_STRUCT.tracking),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("verbose",
                        "CarpetRegrid",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Print screen output while running",
                        "no",
                        &(PRIVATE_CARPETREGRID_STRUCT.verbose),
                        0,
                        NULL,
                        0);

  CCTKi_ParameterCreate("veryverbose",
                        "CarpetRegrid",
                        "BOOLEAN",
                        "PRIVATE",
                        CCTK_STEERABLE_NEVER,
                        "Print much screen output while running",
                        "no",
                        &(PRIVATE_CARPETREGRID_STRUCT.veryverbose),
                        0,
                        NULL,
                        0);

  return 0;
}

int CCTKi_BindingsCarpetRegridParameterExtensions(void);
int CCTKi_BindingsCarpetRegridParameterExtensions(void)
{
  return 0;
}

