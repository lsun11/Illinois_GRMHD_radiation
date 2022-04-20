#/*@@
#  @file      CreateConfigurationBindings.pl
#  @date      Thu Mar 25 14:25:13 2004
#  @author    Yaakoub Y El-Khamra
#  @desc
#             New Configuration.ccl script processing
#  @enddesc
#  @version   $Header: /cactusdevcvs/Cactus/lib/sbin/CreateConfigurationBindings.pl,v 1.10 2005/11/21 18:41:03 schnetter Exp $
#@@*/

require "$sbin_dir/CSTUtils.pl";

#/*@@
#  @routine    CreateConfigurationBindings
#  @date       Thu Mar 25 14:25:13 2004
#  @author     Yaakoub Y El-Khamra
#  @desc
#              Creates the configuration bindings.
#  @enddesc
#@@*/
sub CreateConfigurationBindings
{
  my($bindings_dir, $cfg, $thorns)=@_;
  my($field, $providedcap, $providedcaplist, $thorn, $temp,$defs,$incs,$deps,$linkerflagdirs,$linkerflaglibs);
  my(%linker_thorns, %linker_cfg, $linker_list, $linkerdirs, $linkerlibs);

  if(! $build_dir)
  {
    $build_dir = "$bindings_dir/build";
  }

  if(! -d $bindings_dir)
  {
    mkdir("$bindings_dir", 0755) || die "Unable to create $bindings_dir";
  }

  $start_dir = `pwd`;

  chdir $bindings_dir;

  if(! -d 'Configuration')
  {
    mkdir('Configuration', 0755) || die 'Unable to create Configuration directory';
  }

  if(! -d 'include')
  {
    mkdir('include', 0755) || die 'Unable to create include directory';
  }

  chdir 'Configuration';

  if(! -d "Capabilities")
  {
    mkdir("Capabilities", 0755) || die "Unable to create Capabilities directory";
  }
  if(! -d 'Thorns')
  {
    mkdir('Thorns', 0755) || die "Unable to create Thorns directory";
  }


  # this string goes into the cactus executable directly 
  $linkerflagdirs = '';
  $linkerflaglibs = '';

  $providedcaplist = '';
  # here we put all the provided capabilities where they belong
  foreach $thorn (sort keys %thorns)
  {
    # we know that all the requirements have been satisfied
    # so all we need to do is put the provides where they belong
    # and make references to them from the requirements and optional
    # since we can have multiple provides, we make each capability
    # separate
    if ($cfg->{"\U$thorn\E PROVIDES"})
    {
      foreach $providedcap (split (' ', $cfg->{"\U$thorn\E PROVIDES"}))
      {
        $providedcaplist .= "$providedcap ";
        $temp = '';
        # put include_dirs  and make.definition in one file: make.capability.defn
        if ($cfg->{"\U$thorn $providedcap\E INCLUDE_DIRECTORY"})
        {
          $temp .="INC_DIRS +=" . $cfg->{"\U$thorn $providedcap\E INCLUDE_DIRECTORY"} . "\n";
        }

        if ($cfg->{"\U$thorn $providedcap\E MAKE_DEFINITION"})
        {
          $temp .= $cfg->{"\U$thorn $providedcap\E MAKE_DEFINITION"};
        }

        &WriteFile("Capabilities/make.\U$providedcap\E.defn",\$temp);

        $temp = '';
        # put include and DEFINE in one file: capability.h
        if ($cfg->{"\U$thorn $providedcap\E INCLUDE"})
        {
          $temp .= $cfg->{"\U$thorn $providedcap\E INCLUDE"};
        }

        if ($cfg->{"\U$thorn $providedcap\E DEFINE"})
        {
          $temp .=  "#define " . $cfg->{"\U$thorn $providedcap\E DEFINE"};
        }

        &WriteFile("Capabilities/cctki_\U$providedcap\E.h",\$temp);
        $temp = '';

        # put make.capability.deps in one file: make.capabiltiy.deps
        if ($cfg->{"\U$thorn $providedcap\E MAKE_DEPENDENCY"})
        {
          $temp .= $cfg->{"\U$thorn $providedcap\E MAKE_DEPENDENCY"};
        }

        &WriteFile("Capabilities/make.\U$providedcap\E.deps",\$temp);

        if ( $cfg->{"\U$thorn $providedcap\E LIBRARY"} )
        {
          $linker_thorns{"$thorn"} = $thorn;
          $linker_cfg{"\U$thorn\E USES"} = $cfg->{"\U$thorn\E USES THORNS"};
        }

        if ( $cfg->{"\U$thorn $providedcap\E LIBRARY_DIRECTORY"} )
        {
          $linker_thorns{"$thorn"} = $thorn;
          $linker_cfg{"\U$thorn\E USES"} = $cfg->{"\U$thorn\E USES THORNS"};
        }
      } 
    }
  }

  # here we add the files to the thorns that require capabilities
  foreach $thorn (sort keys %thorns)
  {
    # we know that all the requirements have been satisfied
    # so all we need to do is make references to the capabilities
    # from the requirements since we can have multiple provides,
    # we make each capability separate

    $defs = '';
    $incs = '';
    $deps = '';

    if ($cfg->{"\U$thorn\E REQUIRES"})
    {
      foreach $providedcap (split (' ', $cfg->{"\U$thorn\E REQUIRES"}))
      {
        # put reference to provided capability
        $defs .= "include $bindings_dir/Configuration/Capabilities/make.\U$providedcap\E.defn\n";
        $incs .= "#include \"../Capabilities/cctki_\U$providedcap\E.h\"\n";
        $deps .= "include $bindings_dir/Configuration/Capabilities/make.\U$providedcap\E.deps\n";
      }
    }

   if ($cfg->{"\U$thorn\E OPTIONAL"})
    {
      foreach $providedcap (split (' ', $cfg->{"\U$thorn\E OPTIONAL"}))
      {
        if ($providedcaplist =~ m/$providedcap/i)
        {
          $defs .= $providedcap . " = 1\n";
          $defs .= "include $bindings_dir/Configuration/Capabilities/make.\U$providedcap\E.defn\n";
          $incs .= "#define " . $cfg->{"\U$thorn\E OPTIONAL \U$providedcap\E DEFINE"} . " 1\n";
          $incs .= "#include \"../Capabilities/cctki_\U$providedcap\E.h\"\n";
          $deps .= "include $bindings_dir/Configuration/Capabilities/make.\U$providedcap\E.deps\n";
        }
      }
    }

    if ($cfg->{"\U$thorn\E REQUIRES"} || $cfg->{"\U$thorn\E OPTIONAL"})
    {
      # write everything to file
      # (write the files even if they are empty)
      &WriteFile("./Thorns/make.$thorn.defn",\$defs);
      &WriteFile("./Thorns/cctki_$thorn.h",\$incs);
      &WriteFile("./Thorns/make.$thorn.deps",\$deps);
    }
    else
    {
      # remove the files
      # (we cannot have old files staying around)
      unlink "./Thorns/make.$thorn.defn";
      unlink "./Thorns/cctki_$thorn.h";
      unlink "./Thorns/make.$thorn.deps";
    }
    
  }

  # Sort the linker thorns
  $linkerdirs = 'LIBDIRS +=';
  $linkerlibs = 'LIBS +=';

  $linker_list = &TopoSort(\%linker_thorns, \%linker_cfg);
  foreach $thorn (split (' ', $linker_list))
  {
    foreach $providedcap (split (' ', $cfg->{"\U$thorn\E PROVIDES"}))
    {
      $linkerdirs .= ' ' . $cfg->{"\U$thorn $providedcap\E LIBRARY_DIRECTORY"};
      $linkerlibs .= ' ' . $cfg->{"\U$thorn $providedcap\E LIBRARY"};
    }
  }
  $temp = $linkerdirs . "\n" . $linkerlibs . "\n";
  &WriteFile("make.link",\$temp);

  # write cctki_Capabilities.h file to bindings/include
  # this file adds the if_i_am_thorn stuff
  $temp = '';
  foreach $thorn (sort keys %thorns)
  {
    if ($cfg->{"\U$thorn\E REQUIRES"} || $cfg->{"\U$thorn\E OPTIONAL"})
    {
      $temp .= "#ifdef THORN_IS_$thorn\n";
      $temp .= "#include \"../Configuration/Thorns/cctki_$thorn.h\"\n";
      $temp .= "#endif\n";
      $temp .= "\n";
    }
  }
  &WriteFile("../include/cctki_Capabilities.h",\$temp);
}

return 1;
