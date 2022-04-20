/* $Header: /cactusdevcvs/CactusBase/CoordBase/src/Domain.c,v 1.5 2008/02/20 00:19:32 schnetter Exp $ */



#include "cctk.h"
#include "cctk_Parameters.h"



CCTK_INT CoordBase_GetBoundarySpecification
  (CCTK_INT const size,
   CCTK_INT * const nboundaryzones,
   CCTK_INT * const is_internal,
   CCTK_INT * const is_staggered,
   CCTK_INT * const shiftout)
{
  DECLARE_CCTK_PARAMETERS;
  
  if(! (size>=0))
    CCTK_WARN (0, "size is less than zero");
  if(! (nboundaryzones))
    CCTK_WARN (0, "nboundaryzones is out of bounds");
  if(! (is_internal))
    CCTK_WARN (0, "is_internal is out of bounds");
  if(! (is_staggered))
    CCTK_WARN (0, "is_staggered is out of bounds");
  if(! (shiftout))
    CCTK_WARN (0, "shiftout is out of bounds");
  
  if(! (size == 6))
    CCTK_WARN (0, "size is out of bounds");

  nboundaryzones[0] = boundary_size_x_lower;
  nboundaryzones[1] = boundary_size_x_upper;
  nboundaryzones[2] = boundary_size_y_lower;
  nboundaryzones[3] = boundary_size_y_upper;
  nboundaryzones[4] = boundary_size_z_lower;
  nboundaryzones[5] = boundary_size_z_upper;
  
  is_internal[0] = boundary_internal_x_lower;
  is_internal[1] = boundary_internal_x_upper;
  is_internal[2] = boundary_internal_y_lower;
  is_internal[3] = boundary_internal_y_upper;
  is_internal[4] = boundary_internal_z_lower;
  is_internal[5] = boundary_internal_z_upper;
  
  is_staggered[0] = boundary_staggered_x_lower;
  is_staggered[1] = boundary_staggered_x_upper;
  is_staggered[2] = boundary_staggered_y_lower;
  is_staggered[3] = boundary_staggered_y_upper;
  is_staggered[4] = boundary_staggered_z_lower;
  is_staggered[5] = boundary_staggered_z_upper;
  
  shiftout[0] = boundary_shiftout_x_lower;
  shiftout[1] = boundary_shiftout_x_upper;
  shiftout[2] = boundary_shiftout_y_lower;
  shiftout[3] = boundary_shiftout_y_upper;
  shiftout[4] = boundary_shiftout_z_lower;
  shiftout[5] = boundary_shiftout_z_upper;
  
  return 0;
}



CCTK_INT CoordBase_GetDomainSpecification
  (CCTK_INT const size,
   CCTK_REAL * const physical_min,
   CCTK_REAL * const physical_max,
   CCTK_REAL * const interior_min,
   CCTK_REAL * const interior_max,
   CCTK_REAL * const exterior_min,
   CCTK_REAL * const exterior_max,
   CCTK_REAL * const thespacing)
{
  DECLARE_CCTK_PARAMETERS;
  
  int ierr;
  
  if(! (size>=0))
    CCTK_WARN (0, "size is out of bounds");
  if(! (physical_min))
    CCTK_WARN (0, "physical_min is out of bounds");
  if(! (physical_max))
    CCTK_WARN (0, "physical_max is out of bounds");
  if(! (interior_min))
    CCTK_WARN (0, "interior_min is out of bounds");
  if(! (interior_max))
    CCTK_WARN (0, "interior_max is out of bounds");
  if(! (exterior_min))
    CCTK_WARN (0, "exterior_min is out of bounds");
  if(! (exterior_max))
    CCTK_WARN (0, "exterior_max is out of bounds");
  if(! (thespacing))
    CCTK_WARN (0, "thespacing is out of bounds");

  if(! (size == 3))
    CCTK_WARN (0, "size is out of bounds");  
  
  
  if (CCTK_EQUALS (domainsize, "minmax")) {
    
    physical_min[0] = xmin;
    physical_min[1] = ymin;
    physical_min[2] = zmin;
    physical_max[0] = xmax;
    physical_max[1] = ymax;
    physical_max[2] = zmax;
    if (CCTK_EQUALS (spacing, "gridspacing")) {
      thespacing[0] = dx;
      thespacing[1] = dy;
      thespacing[2] = dz;
    } else if (CCTK_EQUALS (spacing, "numcells")) {
      thespacing[0] = (physical_max[0] - physical_min[0]) / ncells_x;
      thespacing[1] = (physical_max[1] - physical_min[1]) / ncells_y;
      thespacing[2] = (physical_max[2] - physical_min[2]) / ncells_z;
    }
    
  } else if (CCTK_EQUALS (domainsize, "extent")) {
    
    if (zero_origin_x) {
      physical_min[0] = xmin;
      physical_max[0] = xmin + xextent;
    } else {
      physical_min[0] = - xextent / 2;
      physical_max[0] = + xextent / 2;
    }
    if (zero_origin_y) {
      physical_min[1] = ymin;
      physical_max[1] = ymin + yextent;
    } else {
      physical_min[1] = - yextent / 2;
      physical_max[1] = + yextent / 2;
    }
    if (zero_origin_z) {
      physical_min[2] = zmin;
      physical_max[2] = zmin + zextent;
    } else {
      physical_min[2] = - zextent / 2;
      physical_max[2] = + zextent / 2;
    }
    if (CCTK_EQUALS (spacing, "gridspacing")) {
      thespacing[0] = dx;
      thespacing[1] = dy;
      thespacing[2] = dz;
    } else if (CCTK_EQUALS (spacing, "numcells")) {
      thespacing[0] = (physical_max[0] - physical_min[0]) / ncells_x;
      thespacing[1] = (physical_max[1] - physical_min[1]) / ncells_y;
      thespacing[2] = (physical_max[2] - physical_min[2]) / ncells_z;
    }
    
  } else if (CCTK_EQUALS (domainsize, "spacing")) {
    
    if (zero_origin_x) {
      physical_min[0] = xmin;
      physical_max[0] = xmin + dx * ncells_x;
    } else {
      physical_min[0] = - dx * ncells_x / 2;
      physical_max[0] = + dx * ncells_x / 2;
    }
    if (zero_origin_y) {
      physical_min[1] = ymin;
      physical_max[1] = ymin + dy * ncells_y;
    } else {
      physical_min[1] = - dy * ncells_y / 2;
      physical_max[1] = + dy * ncells_y / 2;
    }
    if (zero_origin_z) {
      physical_min[2] = zmin;
      physical_max[2] = zmin + dz * ncells_z;
    } else {
      physical_min[2] = - dz * ncells_z / 2;
      physical_max[2] = + dz * ncells_z / 2;
    }
    thespacing[0] = dx;
    thespacing[1] = dy;
    thespacing[2] = dz;

  } else {
    CCTK_WARN (0, "domainsize is out of bounds");
  }
  
  ierr = ConvertFromPhysicalBoundary
    (size, physical_min, physical_max, interior_min, interior_max,
     exterior_min, exterior_max, thespacing);
  if(ierr)
    CCTK_WARN (0, "Error returned from ConvertFromPhysicalBoundary");  
  
  return 0;
}



CCTK_INT CoordBase_ConvertFromPhysicalBoundary
  (CCTK_INT const size,
   CCTK_REAL const * const physical_min,
   CCTK_REAL const * const physical_max,
   CCTK_REAL * const interior_min,
   CCTK_REAL * const interior_max,
   CCTK_REAL * const exterior_min,
   CCTK_REAL * const exterior_max,
   CCTK_REAL const * const thespacing)
{
  CCTK_INT nboundaryzones[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  
  int d;
  int ierr;
  
  if(! (size>=0))
    CCTK_WARN (0, "size is out of bounds");
  if(! (physical_min))
    CCTK_WARN (0, "physical_min is out of bounds");
  if(! (physical_max))
    CCTK_WARN (0, "physical_max is out of bounds");
  if(! (interior_min))
    CCTK_WARN (0, "interior_min is out of bounds");
  if(! (interior_max))
    CCTK_WARN (0, "interior_max is out of bounds");
  if(! (exterior_min))
    CCTK_WARN (0, "exterior_min is out of bounds");
  if(! (exterior_max))
    CCTK_WARN (0, "exterior_max is out of bounds");
  if(! (thespacing))
    CCTK_WARN (0, "thespacing is out of bounds");

  if(! (size == 3))
    CCTK_WARN (0, "size is out of bounds");  
  
  ierr = CoordBase_GetBoundarySpecification
    (6, nboundaryzones, is_internal, is_staggered, shiftout);
  if (ierr)
    CCTK_WARN (0, "error returned from function CoordBase_GetBoundarySpecification");
  
  for (d=0; d<3; ++d) {
    
    exterior_min[d] = physical_min[d] - thespacing[d] *
      (+ (is_internal[2*d] ? 0 : nboundaryzones[2*d] - 1)
       + (is_staggered[2*d] ? 0.5 : 0.0)
       + shiftout[2*d]);
    exterior_max[d] = physical_max[d] + thespacing[d] *
      (+ (is_internal[2*d+1] ? 0 : nboundaryzones[2*d+1] - 1)
       + (is_staggered[2*d+1] ? 0.5 : 0.0)
       + shiftout[2*d+1]);
  
    interior_min[d] = exterior_min[d] + thespacing[d] * nboundaryzones[2*d];
    interior_max[d] = exterior_max[d] - thespacing[d] * nboundaryzones[2*d+1];

  }
  
  return 0;
}



CCTK_INT CoordBase_ConvertFromInteriorBoundary
  (CCTK_INT const size,
   CCTK_REAL * const physical_min,
   CCTK_REAL * const physical_max,
   CCTK_REAL const * const interior_min,
   CCTK_REAL const * const interior_max,
   CCTK_REAL * const exterior_min,
   CCTK_REAL * const exterior_max,
   CCTK_REAL const * const thespacing)
{
  CCTK_INT nboundaryzones[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  
  int d;
  int ierr;
  
  if(! (size>=0))
    CCTK_WARN (0, "size is out of bounds");
  if(! (physical_min))
    CCTK_WARN (0, "physical_min is out of bounds");
  if(! (physical_max))
    CCTK_WARN (0, "physical_max is out of bounds");
  if(! (interior_min))
    CCTK_WARN (0, "interior_min is out of bounds");
  if(! (interior_max))
    CCTK_WARN (0, "interior_max is out of bounds");
  if(! (exterior_min))
    CCTK_WARN (0, "exterior_min is out of bounds");
  if(! (exterior_max))
    CCTK_WARN (0, "exterior_max is out of bounds");
  if(! (thespacing))
    CCTK_WARN (0, "thespacing is out of bounds");

  if(! (size == 3))
    CCTK_WARN (0, "size is out of bounds");  
  
  ierr = CoordBase_GetBoundarySpecification
    (6, nboundaryzones, is_internal, is_staggered, shiftout);
  if (ierr)
    CCTK_VWarn (0, __LINE__, __FILE__, "CactusBase", "error returned from function CoordBase_GetBoundarySpecification");
  
  for (d=0; d<3; ++d) {
    
    exterior_min[d] = interior_min[d] - thespacing[d] * nboundaryzones[2*d];
    exterior_max[d] = interior_max[d] + thespacing[d] * nboundaryzones[2*d+1];

    physical_min[d] = exterior_min[d] + thespacing[d] *
      (+ (is_internal[2*d] ? 0 : nboundaryzones[2*d] - 1)
       + (is_staggered[2*d] ? 0.5 : 0.0)
       + shiftout[2*d]);
    physical_max[d] = exterior_max[d] - thespacing[d] *
      (+ (is_internal[2*d+1] ? 0 : nboundaryzones[2*d+1] - 1)
       + (is_staggered[2*d+1] ? 0.5 : 0.0)
       + shiftout[2*d+1]);
  
  }
  
  return 0;
}



CCTK_INT CoordBase_ConvertFromExteriorBoundary
  (CCTK_INT const size,
   CCTK_REAL * const physical_min,
   CCTK_REAL * const physical_max,
   CCTK_REAL * const interior_min,
   CCTK_REAL * const interior_max,
   CCTK_REAL const * const exterior_min,
   CCTK_REAL const * const exterior_max,
   CCTK_REAL const * const thespacing)
{
  CCTK_INT nboundaryzones[6];
  CCTK_INT is_internal[6];
  CCTK_INT is_staggered[6];
  CCTK_INT shiftout[6];
  
  int d;
  int ierr;
  
  if(! (size>=0))
    CCTK_WARN (0, "size is out of bounds");
  if(! (physical_min))
    CCTK_WARN (0, "physical_min is out of bounds");
  if(! (physical_max))
    CCTK_WARN (0, "physical_max is out of bounds");
  if(! (interior_min))
    CCTK_WARN (0, "interior_min is out of bounds");
  if(! (interior_max))
    CCTK_WARN (0, "interior_max is out of bounds");
  if(! (exterior_min))
    CCTK_WARN (0, "exterior_min is out of bounds");
  if(! (exterior_max))
    CCTK_WARN (0, "exterior_max is out of bounds");
  if(! (thespacing))
    CCTK_WARN (0, "thespacing is out of bounds");

  if(! (size == 3))
    CCTK_WARN (0, "size is out of bounds");  
  
  ierr = CoordBase_GetBoundarySpecification
    (6, nboundaryzones, is_internal, is_staggered, shiftout);
  if (ierr)
    CCTK_WARN (0, "Error returned from CoordBase_GetBoundarySpecification");  
  
  for (d=0; d<3; ++d) {
    
    physical_min[d] = exterior_min[d] + thespacing[d] *
      (+ (is_internal[2*d] ? 0 : nboundaryzones[2*d] - 1)
       + (is_staggered[2*d] ? 0.5 : 0.0)
       + shiftout[2*d]);
    physical_max[d] = exterior_max[d] - thespacing[d] *
      (+ (is_internal[2*d+1] ? 0 : nboundaryzones[2*d+1] - 1)
       + (is_staggered[2*d+1] ? 0.5 : 0.0)
       + shiftout[2*d+1]);
  
    interior_min[d] = exterior_min[d] + thespacing[d] * nboundaryzones[2*d];
    interior_max[d] = exterior_max[d] - thespacing[d] * nboundaryzones[2*d+1];

  }
  
  return 0;
}
