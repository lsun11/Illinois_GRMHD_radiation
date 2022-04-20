#! /bin/sh
#
# Copyright by The HDF Group.
# All rights reserved.
#
# Tests for the h5checker tool

TOOL=h5check               # The tool name
TOOL_BIN=`pwd`/../tool/$TOOL    # The path of the tool binary

CMP='cmp -s'
DIFF='diff -c'
NLINES=20			# Max. lines of output to display if test fails

nerrors=0
verbose=yes

# The build (current) directory might be different than the source directory.
if test -z "$srcdir"; then
    srcdir=.
fi
test -d ../testfiles || mkdir ../testfiles

# Print a line-line message left justified in a field of 70 characters
# beginning with the word "Testing".
TESTING() {
    SPACES="                                                               "
    echo "Testing $* $SPACES" |cut -c1-70 |tr -d '\012'
}

# Run a test that the tool should pass (exit with 0).  It will print " PASS "
# and "*FAIL*" acoordingly.  When it fails, also display up to $NLINES lines
# of the actual output from the tool test.  The actual output file is usually
# removed unless $HDF5_NOCLEANUP is defined to any non-null value.  
# $* arguments to the TOOL.
TOOLPASS() {
    tmpout=$TOOL.$$.out
    tmperr=$TOOL.$$.err

    # Run test.
    # Stderr is included in stdout so that the diff can detect
    # any unexpected output from that stream too.
    TESTING $TOOL $@
    (
        $TOOL_BIN "$@"
    ) 2> $tmperr > $tmpout
    exitcode=$?
    if [ $exitcode -eq 0 ]; then
        echo " PASSED"
    else
	echo "*FAILED*"
	nerrors="`expr $nerrors + 1`"
	if [ yes = "$verbose" ]; then
	    echo "test returned with exit code $exitcode"
	    echo "test error output:"
	    cat $tmperr
	    echo "***end of test error output***"
	    echo "test output: (up to $NLINES lines)"
	    head -$NLINES $tmpout
	    echo "***end of test output***"
	    echo ""
	fi
    fi

    # Clean up output file
    if test -z "$HDF5_NOCLEANUP"; then
	rm -f $tmpout $tmperr
    fi
}


# Run a test that the tool should fail (exit with non-zero).  It will print
# " PASS " # and "*FAIL*" acoordingly.  When it fails, also display up to
# $NLINES lines of the actual output from the tool test.  The actual output
# file is usually removed unless $HDF5_NOCLEANUP is defined to any non-null
# value.  
# $* arguments to the TOOL.
TOOLFAIL() {
    tmpout=$TOOL.$$.out
    tmperr=$TOOL.$$.err

    # Run test.
    # Stderr is included in stdout so that the diff can detect
    # any unexpected output from that stream too.
    TESTING $TOOL $@
    (
        $TOOL_BIN "$@"
    ) 2> $tmperr > $tmpout
    exitcode=$?
    if [ $exitcode -ne 0 ]; then
        echo " PASSED"
    else
	echo "*FAILED*"
	nerrors="`expr $nerrors + 1`"
	if [ yes = "$verbose" ]; then
	    echo "test returned with exit code $exitcode"
	    echo "test error output:"
	    cat $tmperr
	    echo "***end of test error output***"
	    echo "test output: (up to $NLINES lines)"
	    head -$NLINES $tmpout
	    echo "***end of test output***"
	    echo ""
	fi
    fi

    # Clean up output file
    if test -z "$HDF5_NOCLEANUP"; then
	rm -f $tmpout $tmperr
    fi
}


# Run a test to see if the output matches with expected output and print " PASS " or
# "*FAIL*" accordingly.
# When it fails, also display up to $NLINES lines of the actual output from the tool
# test.  The actual output file is usually removed unless $HDF5_NOCLEANUP is defined
# to any non-null value.
TOOLMATCH() {
    expect="$srcdir/../testfiles/$1"
    actual="../testfiles/`basename $1 .ls`.out"
    shift

    # Run test.
    # Stderr is included in stdout so that the diff can detect
    # any unexpected output from that stream too.
    TESTING $TOOL $@
    (
	echo "#############################"
	echo " output for '$TOOL $@'" 
	echo "#############################"
	cd $srcdir/../testfiles
        $RUNSERIAL $TOOL_BIN "$@"
    ) 2>&1 |sed 's/Modified:.*/Modified:  XXXX-XX-XX XX:XX:XX XXX/' >$actual
    
    exitcode=$?
    if [ $exitcode -ne 0 ]; then
	echo "*FAILED*"
	nerrors="`expr $nerrors + 1`"
	if [ yes = "$verbose" ]; then
	    echo "test returned with exit code $exitcode"
	    echo "test output: (up to $NLINES lines)"
	    head -$NLINES $actual
	    echo "***end of test output***"
	    echo ""
	fi
    elif [ ! -f $expect ]; then
	# Create the expect file if it doesn't yet exist.
        echo " CREATED"
        cp $actual $expect
    elif $CMP $expect $actual; then
        echo " PASSED"
    else
        echo "*FAILED*"
	echo "    Expected result differs from actual result"
	nerrors="`expr $nerrors + 1`"
	test yes = "$verbose" && $DIFF $expect $actual |sed 's/^/    /'
    fi

    # Clean up output file
    if test -z "$HDF5_NOCLEANUP"; then
	rm -f $actual
    fi
}

##############################################################################
##############################################################################
###			  T H E   T E S T S                                ###
##############################################################################
##############################################################################

# Toss in a bunch of tests.  Not sure if they are the right kinds.
# test the help syntax

echo ========================================
echo The following tests are expected to pass.
echo ========================================
TOOLPASS basic_types.h5
TOOLPASS alternate_sb.h5
TOOLPASS array.h5
TOOLPASS attr.h5
TOOLPASS basic_types.h5
TOOLPASS compound.h5
TOOLPASS cyclical.h5
TOOLPASS enum.h5
TOOLPASS external_empty.h5
TOOLPASS external_full.h5
TOOLPASS filters.h5
TOOLPASS group_dsets.h5
TOOLPASS hierarchical.h5
TOOLPASS linear.h5
TOOLPASS log.h5
TOOLPASS multipath.h5
TOOLPASS rank_dsets_empty.h5
TOOLPASS rank_dsets_full.h5
TOOLPASS refer.h5
TOOLPASS root.h5
TOOLPASS stdio.h5
TOOLPASS time.h5
TOOLPASS vl.h5

# these 2 files are generated only when 1.8 library is used
if test -r new_grat.h5; then
TOOLPASS new_grat.h5
fi
if test -r sohm.h5; then
TOOLPASS sohm.h5
fi

# future tests for non-default VFD files
#TOOLPASS family00000.h5
#TOOLPASS family00001.h5
#TOOLPASS family00002.h5
#TOOLPASS split-m.h5
#TOOLPASS split-r.h5
#TOOLPASS multi-b.h5
#TOOLPASS multi-g.h5
#TOOLPASS multi-l.h5
#TOOLPASS multi-o.h5
#TOOLPASS multi-r.h5
#TOOLPASS multi-s.h5

echo ========================================
echo The following tests are expected to fail.
echo ========================================
TOOLFAIL invalidfiles/base_addr.h5
# Temporary block out since this file is not really invalid.
#TOOLFAIL invalidfiles/leaf_internal_k.h5
TOOLFAIL invalidfiles/offsets_lengths.h5
TOOLFAIL invalidfiles/sb_version.h5
TOOLFAIL invalidfiles/signature.h5
# this is a valid 1.8 file
# this should fail when checked against 1.6 format
TOOLFAIL --format=16 invalidfiles/vms_data.h5




if test $nerrors -eq 0 ; then
	echo "All $TOOL tests passed."
fi

exit $nerrors
