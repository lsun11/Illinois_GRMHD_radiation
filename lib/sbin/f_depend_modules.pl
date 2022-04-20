#!/usr/bin/perl -w

#/*@@
#  @file    f_depend_modules.pl
#  @author  Erik Schnetter
#  @date    19 January 2004
#  @desc
#           Create dependencies for Fortran 90 "use" and "include" statements
#  @enddesc
#  @version $Header: /cactusdevcvs/Cactus/lib/sbin/f_depend_modules.pl,v 1.11 2007/10/11 22:51:24 rideout Exp $
# @@*/



use strict;

my $srcfile = $ARGV[0];
my $dest = $ARGV[1];
my $srcdir = $ARGV[2];
my @otherdirs = @ARGV[3..$#ARGV];

my @suffixes = (".f77", ".f", ".f90", ".F77", ".F", ".F90");

print "$dest:";

my %modules;

my $line = 0;
while (<STDIN>)
{
  ++ $line;
  if (/^\s*#\s*(\d+)/)
  {
    # line number directive from C preprocessor
    $line = $1 - 1;
  }
  elsif (/^\s*include\s*['"]([^'"]+)['"]/i)
  {
    # include statement
    my $name = $1;
    my $found = 0;
    if (! $found)
    {
      # reference to an include file in this thorn?
      my $dirhdl;
      if( opendir( DIRHDL, "$srcdir" ) )
      {
        file: while( defined( my $filename = readdir( DIRHDL ) ) )
        {
          if( $filename eq "$name" )
          {
            $found = 1;
            print " \\\n  $filename";
            last file;
          }
        }
        closedir DIRHDL;
      }
    }
    if (! $found)
    {
      # reference to an include file in another thorn?
      dir: foreach my $dir (@otherdirs)
      {
        # note: we could also use the SUBDIRS from the make.code.defn here
        foreach my $subdir (".", "include")
        {
          if( opendir( DIRHDL, "$dir/$subdir" ) )
          {
            while( defined( my $filename = readdir( DIRHDL ) ) )
            {
              if( $filename eq "$name" )
              {
                $found = 1;
                print " \\\n  $dir/$subdir/$filename";
                last dir;
              }
            }
            closedir DIRHDL;
          }
        }
      }
    }
    if (! $found)
    {
      print STDERR "$srcfile:$line: Warning: While tracing include dependencies: Include file \"$name\" not found\n";
      if (@otherdirs)
      {
        print STDERR "   Searched in thorn directory and in [" . join(', ', @otherdirs) . "]\n";
      }
      else
      {
        print STDERR "   Searched in thorn directory only.\n";
      }
    }
  }
  elsif (/^\s*module\s+(\w+)/i)
  {
    # begin of a module
    my $name = $1;
    $modules{$name} = 1;
  }
  elsif (/^\s*use\s+(\w+)/i)
  {
    # use statement
    my $name = $1;
    my $found = 0;
    if (! $found)
    {
      # reference to a module in this file?
      if ($modules{$name})
      {
        $found = 1;
      }
    }
    if (! $found)
    {
      # reference to a module in this thorn?
      my $dirhdl;
      if( opendir( DIRHDL, "$srcdir" ) )
      {
        file: while( defined( my $filename = readdir( DIRHDL ) ) )
        {
          foreach my $suffix (@suffixes)
          {
            if( $filename eq "$name$suffix" )
            {
              $found = 1;
              print " \\\n  $filename.o";
              last file;
            }
          }
        }
        closedir DIRHDL;
      }
    }
    if (! $found)
    {
      # reference to a module in another thorn?
      dir: foreach my $dir (@otherdirs)
      {
        # note: we could also use the SUBDIRS from the make.code.defn here
        foreach my $subdir (".", "include")
        {
          if( opendir( DIRHDL, "$dir/$subdir" ) )
          {
            while( defined( my $filename = readdir( DIRHDL ) ) )
            {
              foreach my $suffix (@suffixes)
              {
                if( $filename eq "$name$suffix.o" )
                {
                  $found = 1;
                  print " \\\n  $dir/$subdir/$filename";
                  last dir;
                }
              }
            }
            closedir DIRHDL;
          }
        }
      }
    }
    if (! $found)
    {
      print STDERR "$srcfile:$line: Warning: While tracing module dependencies: Source file for module \"$name\" not found\n";
      if (@otherdirs)
      {
        print STDERR "   Searched in thorn directory and in [" . join(', ', @otherdirs) . "]\n";
      }
      else
      {
        print STDERR "   Searched in thorn directory only.\n";
      }
    }
  }
}

print "\n";
