#! /usr/bin/perl -w

#/*@@
#  @file     ConfigurationParser.pl
#  @date     Tue Feb  8 17:36:48 2000
#  @author   Tom Goodale
#  @desc
#            Parser for configuration.ccl files
#  @enddesc
#  @version  $Header: /cactusdevcvs/Cactus/lib/sbin/ConfigurationParser.pl,v 1.18 2006/05/15 14:38:25 tradke Exp $
#@@*/


#/*@@
#  @routine    CreateConfigurationDatabase
#  @date       Tue Feb  8 17:47:26 2000
#  @author     Tom Goodale
#  @desc
#              Parses the information in all the thorns' configuration.ccl files
#              and creates a database from it
#  @enddesc
#@@*/
sub CreateConfigurationDatabase
{
  my($config_dir, %thorns) = @_;
  my(%cfg) = ();
  my(%thorn_dependencies);

  # Loop through each thorn's configuration file.
  foreach my $thorn (sort keys %thorns)
  {
    my $filename = "$thorns{$thorn}/configuration.ccl";
    next if (! -r $filename);

    # Get the configuration data from it
    &ParseConfigurationCCL($config_dir, $thorn, \%cfg, \%thorns, $filename);

#    if($debug)
#    {
#      print "   $thorn\n";
#      print "           Provides: ", $cfg{"\U$thorn\E PROVIDES"}, "\n"
#        if ($cfg{"\U$thorn\E PROVIDES"});
#      print "           Requires: ", $cfg{"\U$thorn\E REQUIRES"}, "\n"
#        if ($cfg{"\U$thorn\E REQUIRES"});
#    }

    $cfg{"\U$thorn\E USES THORNS"} = '';

    # verify that all required thorns are there in the ThornList
    next if (! $cfg{"\U$thorn\E REQUIRES THORNS"});

    my @missing = ();
    foreach my $required (split (' ', $cfg{"\U$thorn\E REQUIRES THORNS"}))
    {
      push (@missing, $required)
        if ((! $thorns{"$required"}) && (! $thorns{"\U$required\E"}));
    }
    if (@missing == 1)
    {
      &CST_error (0, "Thorn '$thorn' requires thorn '@missing'. " .
                     'Please add this thorn to your ThornList or remove ' .
                     "'$thorn' from it !");
    }
    elsif (@missing > 1)
    {
      &CST_error (0, "Thorn '$thorn' requires thorns '@missing'. " .
                     'Please add these thorns to your ThornList or ' .
                     "remove '$thorn' from it !");
    }

    $cfg{"\U$thorn\E USES THORNS"} .= $cfg{"\U$thorn\E REQUIRES THORNS"} . ' ';
  }

  foreach my $thorn (sort keys %thorns)
  {
    # verify that all required capabilities are there in the ThornList
    next if (! $cfg{"\U$thorn\E REQUIRES"});

    foreach my $requiredcap (split (' ', $cfg{"\U$thorn\E REQUIRES"}))
    {
      my @found = ();
      foreach my $thorncap (sort keys %thorns)
      {
        foreach my $cap (split (' ', $cfg{"\U$thorncap\E PROVIDES"}))
        {
          push (@found, $thorncap)
            if ("\U$cap\E" eq "\U$requiredcap\E");
        }
      }

      # there must be exactly one thorn providing a required capability
      if (@found == 0)
      {
        &CST_error (0, "Thorn '$thorn' requires the capability " .
                       "'$requiredcap'.\n" .
                       "     Please add a thorn that provides '$requiredcap' " .
                       "to your ThornList or remove '$thorn' from it !")
      }
      elsif (@found > 1)
      {
        &CST_error (0, "More than one thorn provides the capability " .
                       "'$requiredcap'.\n" .
                       "     These thorns are: '@found'.\n" .
                       "     Please use only one !\n");
      }
      else
      {
        $cfg{"\U$thorn\E USES THORNS"} .= $found[0] . ' ';
      }
    }
  }

# Check for cyclic dependencies
# create a hash with thorn-> used thorns (no prefix)
  foreach my $thorn (sort keys %thorns)
  {
    $thorn_dependencies{uc($thorn)}=$cfg{"\U$thorn\E USES THORNS"};
    $thorn_dependencies{uc($thorn)} =~ s/\b$thorn\b//i;
  }

  my $message = &find_dep_cycles(%thorn_dependencies);
  if ("" ne $message)
  {
    $message  =~ s/^\s*//g;
    $message  =~ s/\s*$//g;
    $message =~ s/\s+/->/g;
    $message = "Found a cyclic dependency in configuration requirements:".$message."\n";
    &CST_error(0, $message);
  }

  return \%cfg;
}


#/*@@
#  @routine    ParseConfigurationCCL
#  @date       Tue Feb  8 19:23:18 2000
#  @author     Tom Goodale
#  @desc
#  Parses a configuration.ccl file and generates a database of the values
#  @enddesc
#@@*/
sub ParseConfigurationCCL
{
  my($config_dir, $thorn, $cfg, $thorns, $filename) = @_;
  my(@data);
  my($line_number, $line);
  my($provides, $script, $lang, $options);
  my($optional, $define);

  # Initialise some stuff to prevent perl -w from complaining.

  $cfg->{"\U$thorn\E PROVIDES"} = '';
  $cfg->{"\U$thorn\E REQUIRES"} = '';
  $cfg->{"\U$thorn\E REQUIRES THORNS"} = '';
  $cfg->{"\U$thorn\E OPTIONAL"} = '';
  $cfg->{"\U$thorn\E OPTIONS"}  = '';

  # Read the data
  @data = &read_file($filename);

  for($line_number = 0; $line_number < @data; $line_number++)
  {
    $line = $data[$line_number];
    # Parse the line
    if($line =~ m/^\s*PROVIDES\s*/i)
    {
      $lang = $script = '';
      ($provides, $script, $lang, $options, $line_number) = &ParseProvidesBlock($line_number, \@data);
      $cfg->{"\U$thorn\E PROVIDES"} .= "$provides ";
      if($script)
      {
        $cfg->{"\U$thorn\E PROVIDES \U$provides\E SCRIPT"} = "$thorns->{$thorn}/$script";
      }
      $cfg->{"\U$thorn\E PROVIDES \U$provides\E LANG"} = $lang;
      $cfg->{"\U$thorn\E PROVIDES \U$provides\E OPTIONS"} = $options;

#      if ($script)
#      {
#        print "Running configuration script '$script'\n";
#
#        &ParseConfigScript($config_dir, $provides, $lang, $script,
#                           $thorn, $cfg);
#        print "\n";
#      }

      next;
    }
    elsif($line =~ m/^\s*REQUIRES\s+THORNS\s*:\s*(.*)/i)
    {
      my $temp = $1;
      $temp =~ s/\b$thorn\b//i;
      $cfg->{"\U$thorn\E REQUIRES THORNS"} .= $temp;
#      if ($cfg->{"\U$thorn\E REQUIRES THORNS"})
#      {
#        &CST_error (3, '\'Requires Thorns\' will not be supported in release beta-14' .
#        "\n Please adjust thorn \U$thorn\E to use \'Requires\' instead.");
#      }
    }
    elsif($line =~ m/^\s*REQUIRES\s*(.*)/i)
    {
      $cfg->{"\U$thorn\E REQUIRES"} .= "$1 ";
    }
    elsif($line =~ m/^\s*OPTIONAL\s*/i)
    {
      ($optional, $define, $line_number) = &ParseOptionalBlock($line_number, \@data);
      $cfg->{"\U$thorn\E OPTIONAL"} .= "$optional ";
      $cfg->{"\U$thorn\E OPTIONAL \U$optional\E DEFINE"} = $define;
    }
    elsif($line =~ m/^\s*NO_SOURCE\s*/i)
    {
      $cfg->{"\U$thorn\E OPTIONS"} .= "NO_SOURCE";
    }
    else
    {
      chomp($line);
      &CST_error (0, "Unrecognised line '$line' in configure.ccl of thorn '$thorn'");
    }
  }
}


#/*@@
#  @routine    ParseProvidesBlock
#  @date       Mon May  8 15:52:40 2000
#  @author     Tom Goodale
#  @desc
#  Parses the PROVIDES block in a configuration.ccl file.
#  @enddesc
#@@*/
sub ParseProvidesBlock
{
  my ($line_number, $data) = @_;
  my ($provides, $script, $lang, $options);

  $provides = "";
  $script   = "";
  $lang     = "";
  $options  = [];

  $data->[$line_number] =~ m/^\s*PROVIDES\s*(.*)/i;

  $provides = $1;

  $line_number++;
  if($data->[$line_number] !~ m/^\s*\{\s*$/)
  {
    &CST_error (0, "Error parsing provides block line '$data->[$line_number]'.".
                   'Missing { at start of block');
    $line_number++ while($data->[$line_number] !~ m:\s*\}\s*:);
  }
  else
  {
    while($data->[$line_number] !~ m:\s*\}\s*:)
    {
      $line_number++;
      if($data->[$line_number] =~ m/^\s*SCRIPT\s*(.*)$/i)
      {
        $script = $1;
        next;
      }
      elsif($data->[$line_number] =~ m/^\s*LANG[^\s]*\s*(.*)$/i)
      {
        $lang = $1;
        next;
      }
      elsif($data->[$line_number] =~ m/^\s*OPTIONS[^\s]*\s*(.*)$/i)
      {
        push(@$options, split(' ',$1));
        next;
      }
      elsif($data->[$line_number] =~ m:\s*\}\s*:)
      {
        # do nothing.
      }
      else
      {
        print STDERR "Error parsing provides block line '$data->[$line_number]'\n";
        &CST_error (0, 'Unrecognised statement');
      }
    }
  }

  return ($provides, $script, $lang, $options, $line_number);
}


#/*@@
#  @routine    ParseOptionalBlock
#  @date       Mon May  8 15:52:40 2000
#  @author     Tom Goodale
#  @desc
#  Parses the OPTIONAL block in a configuration.ccl file.
#  @enddesc
#@@*/
sub ParseOptionalBlock
{
  my ($line_number, $data) = @_;
  my ($optional, $define);

  $data->[$line_number] =~ m/^\s*OPTIONAL\s*(.*)/i;

  $optional = $1;

  $define = "";

  $line_number++;

  if($data->[$line_number] !~ m/^\s*\{\s*$/)
  {
    &CST_error (0, "Error parsing optional block line '$data->[$line_number]'".
                ' Missing { at start of block.');
    $line_number++ while($data->[$line_number] !~ m:\s*\}\s*:);
  }
  else
  {
    while($data->[$line_number] !~ m:\s*\}\s*:)
    {
      $line_number++;
      if($data->[$line_number] =~ m/^\s*DEFINE\s*(.*)$/i)
      {
        if($define eq "")
        {
          $define = $1;
          next;
        }
        else
        {
          &CST_error (0, "Error parsing optional block line '$data->[$line_number]' " . 'Only one define allowed.');
        }
      }
      elsif($data->[$line_number] =~ m:\s*\}\s*:)
      {
        # do nothing.
      }
      else
      {
        &CST_error (0, "Error parsing provides block line '$data->[$line_number]' " . 'Unrecognised statement.');
      }
    }
  }

  return ($optional, $define, $line_number);
}

sub print_configuration_database
{
  my($database) = @_;
  my($field);

  foreach $field ( sort keys %$database )
  {
    print "$field has value $database->{$field}\n";
  }
}

1;
