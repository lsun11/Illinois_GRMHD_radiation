#include <stdio.h>

#include "cctki_ActiveThorns.h"

void CCTKi_BindingsThorn_bssn(void);
void CCTKi_BindingsThorn_bssn(void)
{
  const char *name[] = {"bssn", 0};
  const char *implementation[] = {"BSSN", 0};
  const char *ancestors[] =
  {
    "COORDBASE",
    "GRID",
    0,
  };

  const char *friends[] =
  {
    "diagnostics_vacuum",
    "fisheye",
    "lapse",
    "shift",
    "excision",
    0,
  };

  /*
   * Should be able to do below with a constant initializer
   * but sr8000 compiler doesn't like it.
   * So have to laboriously assign values to each member of array.
   */
  struct iAttributeList attributes[5];

  attributes[0].attribute                = "name";
  attributes[0].AttributeData.StringList = name;
  attributes[1].attribute                = "implementation";
  attributes[1].AttributeData.StringList = implementation;
  attributes[2].attribute                = "ancestors";
  attributes[2].AttributeData.StringList = ancestors;
  attributes[3].attribute                = "friends";
  attributes[3].AttributeData.StringList = friends;
  attributes[4].attribute                = 0;
  attributes[4].AttributeData.StringList = 0;

  CCTKi_RegisterThorn(attributes);
}

