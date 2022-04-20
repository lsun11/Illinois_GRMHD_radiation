#! /bin/sh
# /*@@
#  @file   hdf5.sh
#  @date   Fri Jul 30 1999
#  @author Thomas Radke, Yaakoub El-Khamra
#  @desc
#          Setup HDF5 as a thorn
#  @enddesc
# @@*/



# /*@@
#   @routine    CCTK_Search
#   @date       Wed Jul 21 11:16:35 1999
#   @author     Tom Goodale
#   @desc
#   Used to search for something in various directories
#   @enddesc
#@@*/

CCTK_Search()
{
  eval  $1=""
  if test $# -lt 4 ; then
    cctk_basedir=""
  else
    cctk_basedir="$4/"
  fi
  for cctk_place in $2
    do
#      echo $ac_n "  Looking in $cctk_place""...$ac_c" #1>&6
      if test -r "$cctk_basedir$cctk_place/$3" ; then
#        echo "$ac_t""... Found" #1>&6
        eval $1="$cctk_place"
        break
      fi
      if test -d "$cctk_basedir$cctk_place/$3" ; then
#        echo "$ac_t""... Found" #1>&6
        eval $1="$cctk_place"
        break
      fi
#      echo "$ac_t"" No" #1>&6
    done

  return
}

# if HDF5_DIR isn't set take it from a previous configure (if it exists)
if [ -z "$HDF5_DIR" ]; then
  hdf5_configfile="${TOP}/bindings/Configuration/make.HDF5.defn"
  if [ -r ${hdf5_configfile} ]; then
    HDF5_DIR=`grep HDF5_DIR ${hdf5_configfile} | cut -f2 -d'='`
  fi
fi

# Work out which variation of HDF5 is installed
# and set the HDF5 libs, libdirs and includedirs
if [ -z "$HDF5_DIR" ]; then
  echo "BEGIN MESSAGE"
  echo 'HDF5 selected but no HDF5_DIR set. Checking some places...'
  echo "END MESSAGE"

  CCTK_Search HDF5_DIR "/usr /usr/local /usr/local/hdf5 /usr/local/packages/hdf5 /usr/local/apps/hdf5 $HOME c:/packages/hdf5" include/hdf5.h
  if [ -z "$HDF5_DIR" ]; then
     echo "BEGIN ERROR" 
     echo 'Thorn HDF5 requires an external installation of the HDF5 library.  '\
          'Please set HDF5_DIR to the directory of this installation or ' \
	  'remove HDF5 from your configuration ThornList.'
     echo "END ERROR" 
     exit 2
  fi
  echo "BEGIN MESSAGE"
  echo "Found an HDF5 package in $HDF5_DIR"
  echo "END MESSAGE"

  # don't explicitely add standard include and library search paths
  if [ "$HDF5_DIR" != '/usr' -a "$HDF5_DIR" != '/usr/local' ]; then
    HDF5_LIB_DIRS="$HDF5_DIR/lib"
    HDF5_INC_DIRS="$HDF5_DIR/include"
  fi
else
  echo "BEGIN MESSAGE" 
  echo "Using HDF5 package in $HDF5_DIR"
  echo "END MESSAGE"

  HDF5_LIB_DIRS="$HDF5_DIR/lib"
  HDF5_INC_DIRS="$HDF5_DIR/include"
fi
HDF5_LIBS=hdf5


# Check if HDF5 was built with Linux Large File Support (LFS)
grep -qe _LARGEFILE_SOURCE ${HDF5_LIB_DIRS}/libhdf5.settings 2> /dev/null
test_LFS=$?

# Check if the Stream VFD was compiled in
grep -qe '#define H5_HAVE_STREAM 1' ${HDF5_DIR}/include/H5pubconf.h 2> /dev/null
test_stream_vfd=$?

# Check if the MPI I/O VFD was compiled in
grep -qe '#define H5_HAVE_PARALLEL 1' ${HDF5_DIR}/include/H5pubconf.h 2> /dev/null
test_phdf5=$?

if [ -n "$MPI" ]; then
  if [ $test_phdf5 -eq 0 ]; then
    echo "BEGIN MESSAGE" 
    echo 'Found parallel HDF5 library, so Cactus will potentially make use of parallel HDF5 support.'
    echo "END MESSAGE" 
#  else
#    echo "BEGIN MESSAGE" 
#    echo 'Found serial HDF5 library, so Cactus cannot make use of parallel HDF5 support.'
#    echo "END MESSAGE" 
  fi
else
  if [ $test_phdf5 -eq 0 ]; then
    echo "BEGIN ERROR" 
    echo "Found parallel HDF5 library, but Cactus wasn't configured with MPI."
    echo 'Please set HDF5_DIR to point to a serial HDF5 package, or configure Cactus with MPI.'
    echo "END ERROR" 
    exit 2
  fi
fi


# check that we have the right version of HDF5 under 32/64 bit IRIX
# This should better be checked by some autoconf script.
if [ -n "$IRIX_BITS" ]; then
  if [ -r "$HDF5_LIB_DIRS/libhdf5.a" ]; then
    hdf5_lib="$HDF5_LIB_DIRS/libhdf5.a"
  elif [ -r "$HDF5_LIB_DIRS/libhdf5.so" ]; then
    hdf5_lib="$HDF5_LIB_DIRS/libhdf5.so"
  else
    hdf5_lib=
  fi

  if [ -n "$hdf5_lib" ]; then
    file $hdf5_lib | grep -qe $IRIX_BITS 2> /dev/null
    if [ $? -ne 0 ]; then
      echo "BEGIN ERROR" 
      echo "The HDF5 library found in \"$HDF5_LIB_DIRS\" was not compiled as $IRIX_BITS bits !"
      echo 'Please reconfigure Cactus with the correct setting for HDF5_DIR !'
      echo "END ERROR" 
      exit 1
    fi
  fi
fi


# check whether we run Windows or not
perl -we 'exit (`uname` =~ /^CYGWIN/)'
is_windows=$?


# check whether we run MacOS or not
perl -we 'exit (`uname` =~ /^Darwin/)' 
is_macos=$?



# Check whether we have to link with libsz.a
grep -qe '#define H5_HAVE_LIBSZ 1' ${HDF5_DIR}/include/H5pubconf.h 2> /dev/null
test_szlib=$?
if [ $test_szlib -eq 0 ]; then
  if [ $is_windows -ne 0 ]; then
    libsz='szlib.lib'
  elif [ $is_macos -ne 0 ]; then
    libsz='libsz.dylib'
  else
    libsz='libsz.a'
  fi

  if [ -z "$LIBSZ_DIR" -a ! -r /usr/lib/$libsz ]; then
    echo "BEGIN MESSAGE"
    echo "HDF5 library was built with external szlib I/O filter, searching for library $libsz ..."
    echo "END MESSAGE"

    CCTK_Search LIBSZ_DIR '/usr/local/lib c:/packages/libsz/lib c:/packages/hdf5/lib' $libsz
    if [ -z "$LIBSZ_DIR" ]; then
      echo "BEGIN ERROR" 
      echo "Unable to locate the library $libsz - please set LIBSZ_DIR"
      echo "END ERROR" 
      exit 2
    fi
#    echo "  Found library $libsz in $LIBSZ_DIR"
  fi
  if [ $is_windows -eq 0 ]; then
    HDF5_LIBS="$HDF5_LIBS sz"
  else
    HDF5_LIBS="$HDF5_LIBS szlib"
  fi
  HDF5_LIB_DIRS="$HDF5_LIB_DIRS $LIBSZ_DIR"
fi


# Check whether we have to link with libz.a
# this is for current versions of HDF5 (starting from 1.4.x)
grep -qe '#define H5_HAVE_LIBZ 1' ${HDF5_DIR}/include/H5pubconf.h 2> /dev/null
test_zlib=$?

# this is for old versions of HDF5 (before 1.4.x)
if [ $test_zlib -ne 0 ]; then
  grep -qe '#define HAVE_LIBZ 1' ${HDF5_DIR}/include/H5config.h 2> /dev/null
  test_zlib=$?
fi

if [ $test_zlib -eq 0 ]; then
  if [ $is_windows -ne 0 ]; then
    libz='zlib.lib'
  elif [ $is_macos -ne 0 ]; then
    libz='libz.dylib'
  else
    libz='libz.a'
  fi
  if [ -z "$LIBZ_DIR" -a ! -r /usr/lib/$libz -a ! -r /usr/lib64/$libz ]; then
    echo "BEGIN MESSAGE"
    echo "HDF5 library was built with external deflate I/O filter, searching for library $libz ..."
    echo "END MESSAGE"
    CCTK_Search LIBZ_DIR '/usr/local/lib c:/packages/libz/lib c:/packages/hdf5/lib' $libz
    if [ -z "$LIBZ_DIR" ]; then
       echo "BEGIN ERROR" 
       echo "Unable to locate the library $libz - please set LIBZ_DIR"
       echo "END ERROR" 
       exit 2
    fi
    echo "BEGIN MESSAGE"
    echo "Found library $libz in $LIBZ_DIR"
    echo "END MESSAGE"
  fi
  if [ $is_windows -eq 0 ]; then
    HDF5_LIBS="$HDF5_LIBS z"
  else
    HDF5_LIBS="$HDF5_LIBS zlib"
  fi
  HDF5_LIB_DIRS="$HDF5_LIB_DIRS $LIBZ_DIR"
fi

# Add MPI libraries for parallel HDF5
if [ $test_phdf5 -eq 0 ]; then
  HDF5_LIBS="$HDF5_LIBS \$(MPI_LIBS)" 
  HDF5_INC_DIRS="$HDF5_INC_DIRS \$(MPI_INC_DIRS)"
  HDF5_LIB_DIRS="$HDF5_LIB_DIRS \$(MPI_LIB_DIRS)"
fi

# Finally, add the math lib which might not be linked against by default
if [ $is_windows -eq 0 ]; then
  HDF5_LIBS="$HDF5_LIBS m"
fi

echo "INCLUDE_DIRECTORY $HDF5_INC_DIRS"
echo "LIBRARY $HDF5_LIBS"
echo "LIBRARY_DIRECTORY $HDF5_LIB_DIRS"

# remember HDF5_DIR setting so that it can be used in later reconfigs
echo 'BEGIN MAKE_DEFINITION'
echo 'HAVE_HDF5 = 1'
echo "HDF5_DIR = $HDF5_DIR"
echo 'END MAKE_DEFINITION'

if [ $test_stream_vfd -eq 0 ]; then
  echo 'BEGIN MAKE_DEFINITION'
  echo 'HAVE_HDF5_STREAM_VFD = 1'
  echo 'END MAKE_DEFINITION'
fi
if [ $test_LFS -eq 0 ]; then
  echo 'BEGIN MAKE_DEFINITION'
  echo 'HDF5_LFS_FLAGS       = -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64'
  echo 'END MAKE_DEFINITION'
fi  
