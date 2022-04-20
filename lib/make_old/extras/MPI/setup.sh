#! /bin/sh
# /*@@
#   @file      setup.sh
#   @date      Wed Jul 21 11:18:40 1999
#   @author    Tom Goodale
#   @desc 
#   Setup MPI
#   @enddesc 
#   @version $Header: /cactusdevcvs/Cactus/lib/make/extras/MPI/setup.sh,v 1.9 2004/04/10 18:12:14 schnetter Exp $
# @@*/

if test -n "$MPI" -a "$MPI" != 'none' ; then

echo 'Configuring with MPI'

# Work out which variation of MPI is installed

if test -r $srcdir/extras/MPI/$MPI ; then
  . $srcdir/extras/MPI/$MPI
else
  echo "  MPI selected, but no known MPI method - what is \"$MPI\" ?"
  exit 2
fi

# Write the data out to the header and make files.

if test -z "$MPI_VERSION" ; then
  MPI_VERSION="$MPI"
fi
CCTK_WriteLine cctk_Extradefs.h "#define CCTK_MPI_$MPI_VERSION 1"
CCTK_WriteLine cctk_Extradefs.h "#define CCTK_MPI CCTK_MPI_$MPI_VERSION"

CCTK_WriteLine make.extra.defn "HAVE_MPI     = 1"
CCTK_WriteLine make.extra.defn "MPI_LIBS     = $MPI_LIBS"
CCTK_WriteLine make.extra.defn "MPI_LIB_DIRS = $MPI_LIB_DIRS"
CCTK_WriteLine make.extra.defn "MPI_INC_DIRS = $MPI_INC_DIRS"
CCTK_WriteLine make.extra.defn ''
CCTK_WriteLine make.extra.defn 'LIBS         += $(MPI_LIBS)'
CCTK_WriteLine make.extra.defn 'LIBDIRS      += $(MPI_LIB_DIRS)'
CCTK_WriteLine make.extra.defn 'SYS_INC_DIRS += $(MPI_INC_DIRS)'

fi
