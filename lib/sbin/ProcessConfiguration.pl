#/*@@
#  @file      ProcessConfiguration.pl
#  @date      Mon May  8 15:52:08 2000
#  @author    Tom Goodale
#  @desc 
#  
#  @enddesc 
#@@*/

#/*@@
#  @routine    SplitThorns
#  @date       Mon May  8 16:04:59 2000
#  @author     Tom Goodale
#  @desc 
#  Splits the thorns hash into those with source and those without source
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#
#@@*/
sub SplitThorns
{
  my ($configuration_database, $thorns, $source_thorns, $nosource_thorns) = @_;

  foreach $thorn (sort keys %$thorns)
  {
    if($configuration_database->{"\U$thorn OPTIONS\E"} =~ m/NO_SOURCE/i)
    {
      $nosource_thorns->{"$thorn"} = $thorns->{"$thorn"};
    }
    else
    {
      $source_thorns->{"$thorn"} = $thorns->{"$thorn"};
    }
  }
}

#/*@@
#  @routine    ProcessConfiguration
#  @date       Thu Aug 26 22:09:26 2004
#  @author     Tom Goodale
#  @desc 
#  Runs all the configuration scripts belonging to thorns.
#  Has to setup the environment correctly, and as a result
#  modifies the config-info file for future reference.
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#
#@@*/
sub ProcessConfiguration
{
  my($config_dir,$config_database, $thorns, $config_file) = @_;

  my $thorn;
  my $provides;
  my @allowed_opts;

  # Find the master list of allowed options
  foreach $thorn (sort keys %thorns)
  {
#    print "DEBUG: Thorn $thorn\n";

    foreach $provides (split(' ',$config_database->{"\U$thorn\E PROVIDES"}))
    {
#      print "DEBUG: Provides $provides\n";

      if(@{$config_database->{"\U$thorn\E PROVIDES \U$provides\E OPTIONS"}} != 0)
      {
#        print "DEBUG: $thorn provides $provides with options\n";
        push(@allowed_opts, @{$config_database->{"\U$thorn\E PROVIDES \U$provides\E OPTIONS"}});
      }
    }
  }

#  print "DEBUG: allowed options are @allowed_opts\n";

  # Now get all the configuration options.
  my $env;
  my $optfile;
  my $configinfo;
  my $headers;

  ($configinfo,$headers) = ParseConfigInfo($config_file);

  if($ENV{"options"})
  {
    $optfile = ParseOptionsFile($ENV{"options"})
  }
  else
  {
    $optfile = {};
  }

  $env = GetOptionsFromEnv(\%ENV, \@allowed_opts);

  my $modified = AmalgamateOptions($env,$optfile,$configinfo,\@allowed_opts);

  # Write a new config-info file if anything has changed
  if($modified)
  {
    WriteNewConfigInfo($config_file,$headers,$configinfo);
  }

  # Now setup the environment
  # FIXME: Would like to restrict this to just the @allowed_opts, but then
  # flesh configuration options like MPI or arch specific ones like
  # IRIX_BITS would not be propogated 8-(

  foreach $option (keys %$configinfo)
  {
    $ENV{$option} = $configinfo->{$option};
  }

  # Ok, can now run the configuration scripts.

  foreach $thorn (sort keys %thorns)
  {
    foreach $provides (split(' ',$config_database->{"\U$thorn\E PROVIDES"}))
    {
      my $script = $config_database->{"\U$thorn\E PROVIDES \U$provides\E SCRIPT"};
      my $lang   = $config_database->{"\U$thorn\E PROVIDES \U$provides\E LANG"};

      if ($script)
      {
        print "Running configuration script '$script'\n";

        &ParseConfigScript($config_dir, $provides, $lang, $script,
                           $thorn, $config_database);
        print "\n";
      }
    }
  }
}

1;
