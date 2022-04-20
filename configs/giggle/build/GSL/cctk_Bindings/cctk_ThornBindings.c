#include <stdio.h>

#include "cctki_ActiveThorns.h"

void CCTKi_BindingsThorn_GSL(void);
void CCTKi_BindingsThorn_GSL(void)
{
  const char *name[] = {"GSL", 0};
  const char *implementation[] = {"GSL", 0};
  /*
   * Should be able to do below with a constant initializer
   * but sr8000 compiler doesn't like it.
   * So have to laboriously assign values to each member of array.
   */
  struct iAttributeList attributes[3];

  attributes[0].attribute                = "name";
  attributes[0].AttributeData.StringList = name;
  attributes[1].attribute                = "implementation";
  attributes[1].AttributeData.StringList = implementation;
  attributes[2].attribute                = 0;
  attributes[2].AttributeData.StringList = 0;

  CCTKi_RegisterThorn(attributes);
}
