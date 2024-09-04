#!/usr/bin/perl
#
# sqp_opts.pm
# Eric Nawrocki
# EPN, Tue Oct 28 14:16:44 2014 [rnavore (rvr-options.pm)]
# EPN, Tue Feb 23 13:34:31 2016 [incorporation into dnaorg]
# EPN, Tue Jul  2 11:40:49 2019 [migrated from epn-options]
# version: 0.10
#
# Many aspects of this module are modelled after Easel's esl_getopts module
# by Sean Eddy (bioeasel.org).
# 
use strict;
use warnings;

# NOTE: do not add any 'require' statements here, e.g. 'require
# sqp_utils.pm' because the program that uses sequip must handle that
# so each program can specify sequip from a specific directory defined
# by a specific environment variable. This is how, for example,
# ribovore can require a specific version of sequip on the same file
# system that has vadr installed with a potentially different version
# of sequip.

#####################################################################
# Data structures used in this module:
#
# opt_HH: 2D hash:
#         1D key: option name (e.g. "-h")
#         2D key: string denoting type of information 
#                 (one of "type", "default", "group", "requires", "incompatible", "preamble", "help")
#         value:  string explaining 2D key:
#                 "type":          "boolean", "string", "integer" or "real"
#                 "default":       default value for option, in double or single quotes
#                 "group":         positive integer denoting group number this option belongs to
#                 "requires":      string of 0 or more other options this option requires to work, each separated by a ','
#                                  special case: "*" for 'required' option ('incompatible' must be "*" as well)
#                 "incompatible":  string of 0 or more other options this option is incompatible with, each separated by a ','
#                                  special case: "*" for 'required' option ('required' must be "*" as well)
#                 "preamble":      string describing option for preamble section (beginning of output from script)
#                 "help":          string describing option for help section (printed if -h used)
#                 "setby":         '1' if option set by user, else 'undef'
#                 "value":         value for option, can be undef if default is undef
#
# opt_order_A: array of options in the order they should be processed
# 
# opt_group_desc_H: key: group number (integer), value: description of group for help output
#
#####################################################################

#####################################################################
# Subroutine: opt_Add()
# Incept:     EPN, Tue Oct 28 14:16:54 2014
# 
# Purpose:    Add an option to a growing 2D hash (<opt_HHR>) and 
#             'options order' array (<opt_order_AR>).
#
# Arguments:
# $optname:      string for option (e.g. "-h")
# $takes_arg:    '1' if option takes an argument, else '0'
# $default:      default value for option
# $group:        group this option belongs to (integer)
# $requires:     string of >= 0 other options required in combo with this option
#                each separated by a ','
# $incompatible: string of >= 0 other options incompatible with this option
#                each separated by a ','
# $preamble:     string describing this option for preamble section
# $help:         string describing this option for help section
# $opt_HHR:      ref to %opt_HH (see description at top of this file)
# $opt_order_AR: ref to @opt_order_A (see description at top of this file)
#
# Returns:    void
# Dies:       if somethings wrong with option being added 
#             (This must be a coding error, not a user error)
#
####################################################################
sub opt_Add {
  my $narg_expected = 10;
  my $sub_name = "opt_Add()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($optname, $type, $default, $group, $requires, $incompatible, $preamble, $help, $opt_HHR, $opt_order_AR) = @_;

  my $naming_convention_desc = "\tOptions must be of the form of:\n";
  $naming_convention_desc   .= "\t\t-a   (single dash followed by single char) or\n";
  $naming_convention_desc   .= "\t\t--aa (double dash followed by > 1 char)\n";
  $naming_convention_desc   .= "\tAnd no internal dashes or commas allowed.\n"; 
  my $tmp_option;
  my $is_required_flag = 0; # set to '1' if this option is 'required'

  # contract checks
  if(exists $opt_HHR->{$optname})                                    { print STDERR ("ERROR, $sub_name trying to add $optname, but it already exists.\n"); exit(1); } 
  if($optname !~ m/^\-[^\-,]$/ && $optname !~ m/^\-\-[^\-,][^\-,]+/) { print STDERR ("ERROR, $sub_name option $optname violates naming convention:\n$naming_convention_desc."); exit(1); }
  if($group !~ m/^\d+$/)                                             { print STDERR ("ERROR, $sub_name trying to add $optname, group is not an integer ($group).\n"); exit(1); } 


  # check requires string
  if(defined $requires) { 
    if($requires eq "*") { 
      if(! defined $incompatible) { 
        printf STDERR ("ERROR, $sub_name, option $tmp_option has special value $requires in requires string, but incompatible string is undefined (it should be $requires also).\n"); 
        exit(1);
      }
      if($incompatible ne "*") { 
        printf STDERR ("ERROR, $sub_name, option $tmp_option has special value $requires in requires string, but incompatible string is $incompatible (it should be $requires also).\n"); 
        exit(1);
      }
      if(defined $default) { 
        printf STDERR ("ERROR, $sub_name, with required options, default value must be 'undef', but for $optname default value is $default\n");
        exit(1);
      }
      # if we get here, both $requires and $incompatible are "*", make this option required
      $is_required_flag = 1;
    }
    else { 
      foreach $tmp_option (split(",", $requires)) { 
        if($tmp_option !~ m/^\-[^\-]$/ && 
           $tmp_option !~ m/^\-\-[^\-][^\-]+/) { 
          printf STDERR ("ERROR, $sub_name, option $tmp_option in requires string ($requires) violates naming convention:\n$naming_convention_desc.\n"); 
          exit(1); 
        }
      }
    }
  }

  # check incompatible string
  if(defined $incompatible) { 
    if($incompatible eq "*") { 
      if(! defined $incompatible) { 
        printf STDERR ("ERROR, $sub_name, option $tmp_option has special value $incompatible in incompatible string, but requires string is undefined (it should be $incompatible also).\n"); 
        exit(1);
      }
      if($requires ne "*") { 
        printf STDERR ("ERROR, $sub_name, option $tmp_option has special value $incompatible in incompatible string, but requires string is $requires (it should be $incompatible also).\n"); 
        exit(1);
      }
      # if we get here, both $requires and $incompatible are "*", make this option required
      if(! $is_required_flag) { 
        printf STDERR ("ERROR, $sub_name, internal error making $tmp_option required.\n");
        exit(1);
      }
    }
    else { 
      foreach $tmp_option (split(",", $incompatible)) { 
        if($tmp_option !~ m/^\-[^\-]$/ && 
           $tmp_option !~ m/^\-\-[^\-][^\-]+/) { 
          printf STDERR ("ERROR, $sub_name, option $tmp_option in incompatible string ($incompatible) violates naming convention:\n$naming_convention_desc.\n"); 
          exit(1); 
        }
      }
    }
  }
  
  # make sure $type makes sense in conjuction with $default and $help
  # do this after we set $is_required_flag, so we can ignore type/default inconsistencies for boolean if $is_required_flag
  if($type eq "boolean") { 
    if(! $is_required_flag) { 
      if(! defined $default) { 
        print STDERR ("ERROR, $sub_name trying to add $optname of 'boolean' type but default is undefined.\n"); exit(1); 
      }
      if($default != 0 && $default != 1) { 
        print STDERR ("ERROR, $sub_name trying to add $optname of 'boolean' type but default is not '0' or '1'.\n"); exit(1); 
      } 
    }
    if(defined $help && $help =~ m/\<\w\>/) { 
      print STDERR ("ERROR, $sub_name trying to add $optname of 'boolean' type but help string has <> in it ($help).\n"); exit(1); 
    }
  }
  elsif($type eq "integer") { 
    if(defined $default && (! verify_integer($default))) { 
      print STDERR ("ERROR, $sub_name trying to add $optname of 'integer' type but default is $default.\n"); exit(1); 
    }
    if(defined $help && $help !~ m/\<n\>/) { 
      print STDERR ("ERROR, $sub_name trying to add $optname of 'integer' type but help string does not have <n> in it ($help).\n"); exit(1); 
    }
  }
  elsif($type eq "real") { 
    if(defined $default && (! verify_real($default))) { 
      print STDERR ("ERROR, $sub_name trying to add $optname of 'real' type but default is $default.\n"); exit(1); 
    }
    if(defined $help && $help !~ m/\<x\>/) { 
      print STDERR ("ERROR, $sub_name trying to add $optname of 'real' type but help string does not have <x> in it ($help).\n"); exit(1); 
    }
  }
  elsif($type eq "string") { 
    if(defined $help && $help !~ m/\<s\>/) { 
      print STDERR ("ERROR, $sub_name trying to add $optname of 'string' type but help string does not have <s> in it ($help).\n"); exit(1); 
    }
  }
  else { 
    print STDERR ("ERROR, $sub_name trying to add $optname of unknown type: $type.\n"); exit(1); 
  }

  # initialize
  %{$opt_HHR->{$optname}} = ();
  # undefine $incompatible and $requires if $is_required_flag raised
  if($is_required_flag) { 
    $incompatible = undef;
    $requires     = undef; 
  }
  # fill
  $opt_HHR->{$optname}{"type"}         = $type;
  $opt_HHR->{$optname}{"default"}      = $default;
  $opt_HHR->{$optname}{"group"}        = $group;
  $opt_HHR->{$optname}{"requires"}     = $requires;
  $opt_HHR->{$optname}{"incompatible"} = $incompatible;
  $opt_HHR->{$optname}{"preamble"}     = $preamble;
  $opt_HHR->{$optname}{"help"}         = $help;

  $opt_HHR->{$optname}{"setby"}        = "default";
  $opt_HHR->{$optname}{"value"}        = $default;
  if($is_required_flag) { 
    $opt_HHR->{$optname}{"is_required"} = 1;
  }
  else { 
    $opt_HHR->{$optname}{"is_required"} = 0;
  }

  push(@{$opt_order_AR}, $optname);

  return;
}

#####################################################################
# Subroutine: opt_IsDefault()
# Incept:     EPN, Wed Oct 29 06:21:43 2014
# 
# Purpose:    Returns '1' if an option remained at default setting,
#             i.e. was not set at command line, else returns '0'.
#             Returns '0' if option was set at command line to its
#             default value.
#    
#             Modelled after easel's esl_getopts.c::esl_opt_IsDefault().
#
# Arguments:
# $optname:       string for option (e.g. "-h")
# $opt_HHR:      ref to %opt_HH (see description at top of this file)
#
# Returns:    void
# Dies:       if $optname does not exist in $opt_HHR.
#             (This must be a coding error, not a user error)
#
####################################################################
sub opt_IsDefault {
  my $narg_expected = 2;
  my $sub_name = "opt_IsDefault()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($optname, $opt_HHR) = @_;

  if(! exists $opt_HHR->{$optname}) { printf STDERR ("ERROR, $sub_name option $optname does not exist.\n"); exit(1); } 

  return ($opt_HHR->{$optname}{"setby"} eq "default") ? 1 : 0;
}

#####################################################################
# Subroutine: opt_IsOn()
# Incept:     EPN, Wed Oct 29 06:25:02 2014
# 
# Purpose:    Returns '1' if an option is set to a defined value
#             that is not '0'.
#
#             Most useful when using options that take an argument
#             also as boolean switches, where they can either be
#             OFF, or they can be turned ON by having a value. 
#             Caller can test <opt_IsOn()> to see if the option's
#             active at all, then use <opt_Get()> to extract
#             the option value.
#
#             Modelled after easel's esl_getopts.c::esl_opt_IsOn().
#
# Arguments:
# $optname:   string for option (e.g. "-h")
# $opt_HHR:   ref to %opt_HH (see description at top of this file)
#
# Returns:    void
# Dies:       if $optname does not exist in $opt_HHR.
#             (This must be a coding error, not a user error)
#
####################################################################
sub opt_IsOn {
  my $narg_expected = 2;
  my $sub_name = "opt_IsOn()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($optname, $opt_HHR) = @_;

  if(! exists $opt_HHR->{$optname}) { printf STDERR ("ERROR, $sub_name option $optname does not exist.\n"); exit(1); } 

  return (defined $opt_HHR->{$optname}{"value"} && $opt_HHR->{$optname}{"value"} ne 0) ? 1 : 0;
}

#####################################################################
# Subroutine: opt_IsUsed()
# Incept:     EPN, Wed Oct 29 08:58:37 2014
# 
# Purpose:    Returns '1' if option <$opt_name> is in use: it has
#             been set by the user, and that value corresponds to
#             the option being "on" (a defined value).
# 
#             This is useful in at least the following context:
#             when we print a preamble to a program that includes
#             all the options that are in effect that weren't already
#             on by default, i.e. in options_OutputPreamble().
#
# Arguments:
# $optname:   string for option (e.g. "-h")
# $opt_HHR:   ref to %opt_HH (see description at top of this file)
#
# Returns:    void
# Dies:       if $optname does not exist in $opt_HHR.
#             (This must be a coding error, not a user error)
#
####################################################################
sub opt_IsUsed {
  my $narg_expected = 2;
  my $sub_name = "opt_IsUsed()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($optname, $opt_HHR) = @_;

  if(! exists $opt_HHR->{$optname}) { printf STDERR ("ERROR, $sub_name option $optname does not exist.\n"); exit(1); } 

  return ($opt_HHR->{$optname}{"setby"} eq "user") ? 1 : 0;
}


#####################################################################
# Subroutine: opt_Get()
# Incept:     EPN, Wed Oct 29 06:12:46 2014
# 
# Purpose:    Returns value of option set by user.
#             Returns default value if not set by user,
#             which may be 'undef'. 
#             Dies if option does not exist.
#
# Arguments:
# $optname:   string for option (e.g. "-n")
# $opt_HHR:   ref to %opt_HH (see description at top of this file)
#
# Returns:    see Purpose.
# Dies:       if $optname does not exist in $opt_HHR.
#             (This must be a coding error, not a user error)
#
####################################################################
sub opt_Get { 
  my $narg_expected = 2;
  my $sub_name = "opt_Get()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($optname, $opt_HHR) = @_;

  if(! exists $opt_HHR->{$optname}) { printf STDERR ("ERROR, $sub_name option $optname does not exist.\n"); exit(1); } 

  return $opt_HHR->{$optname}{"value"}; # may be 'undef';
}

#####################################################################
# Subroutine: opt_SetFromUserHash()
# Incept:     EPN, Wed Oct 29 05:59:12 2014
# 
# Purpose:    Given a hash ($optset_HR) indicating which options have been set
#             by the user, and to what value, set the same options to
#             the same values in $opt_HHR.
#             $optset_HR is usually a reference to a hash returned from
#             GetOptions(). 
#             If any options that exist in $opt_HHR do not exist in 
#             $opt_set_HR, then we exit (options that are not set
#             will exist in $opt_set_HR, but will have value 'undef').
#
# Arguments:
# $optset_HR: ref to a hash: key: an option string, value: value option is set to
# $opt_HHR:   ref to %opt_HH (see description at top of this file)
#
# Returns:    void
# Dies:       if an option in $optset_HR does not exist in $opt_HHR
#             or if an option in $opt_HHR does not exist in $optset_HR.
#             (This must be a coding error, not a user error)
#
####################################################################
sub opt_SetFromUserHash {
  my $narg_expected = 2;
  my $sub_name = "opt_SetFromUserHash()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($optset_HR, $opt_HHR) = @_;

  # check that each option in $opt_HHR also exists (defined or not) in $optset_HR
  foreach my $optname ( keys %{$opt_HHR}) { 
    if(! exists $optset_HR->{$optname}) { 
      printf STDERR ("ERROR, $sub_name option $optname does not exist in GetOptions() user hash.\n", scalar(@_), $narg_expected); 
      exit(1); 
    }       
  }

  # for each option defined in $optset_HR, set it in $opt_HHR
  foreach my $optname ( sort keys %{$optset_HR} ) { 
    if(defined $optset_HR->{$optname}) { 
      opt_SetByUser($optname, $optset_HR->{$optname}, $opt_HHR);
    }
  }

  return;
}

#####################################################################
# Subroutine: opt_SetByUser()
# Incept:     EPN, Wed Oct 29 06:03:07 2014
# 
# Purpose:    Given an option ($optname) and a value ($value), 
#             update $opt_HHR to reflect that $optname is set to $value
#             by the user.
#
# Arguments:
# $optname:    option to update
# $value:     value to set $optname to 
# $opt_HHR:   ref to %opt_HH (see description at top of this file)
#
# Returns:    void
# Dies:       if $optname does not exist, or $value is incompatible.
#             (This must be a coding error, not a user error)
#
####################################################################
sub opt_SetByUser {
  my $narg_expected = 3;
  my $sub_name = "opt_SetByUser()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($optname, $value, $opt_HHR) = @_;

  if(! exists $opt_HHR->{$optname})             { printf STDERR ("ERROR, $sub_name option $optname does not exist.\n"); exit(1); } 
  if($opt_HHR->{$optname}{"setby"} eq "user")   { printf STDERR ("ERROR, $sub_name option $optname already set by user.\n"); exit(1); } 

  if($opt_HHR->{$optname}{"type"} eq "boolean") { 
    if($value != 0 && $value != 1) { 
      printf STDERR ("ERROR, $sub_name option $optname trying to set to $value, but it's a boolean.\n"); exit(1); } 
  }
  if($opt_HHR->{$optname}{"type"} eq "integer") { 
    if(! verify_integer($value)) { 
      printf STDERR ("ERROR, $sub_name option $optname trying to set to $value, but it's an integer.\n"); exit(1); } 
  }
  if($opt_HHR->{$optname}{"type"} eq "real") { 
    if(! verify_real($value)) { 
      printf STDERR ("ERROR, $sub_name option $optname trying to set to $value, but it's a real.\n"); exit(1); } 
  }

  $opt_HHR->{$optname}{"value"} = $value;
  $opt_HHR->{$optname}{"setby"} = "user";

  return;
}

#####################################################################
# Subroutine: opt_Exists()
# Incept:     EPN, Wed Mar  2 14:37:02 2016
# 
# Purpose:    Given an option ($optname) check if it exists
#             in $opt_HHR. Return '1' if it does, else return '0'.
#
# Arguments:
# $optname:   option to check for
# $opt_HHR:   ref to %opt_HH (see description at top of this file)
#
# Returns:    void
#
# Dies:       neve
#
####################################################################
sub opt_Exists {
  my $narg_expected = 2;
  my $sub_name = "opt_Exists()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($optname, $opt_HHR) = @_;

  return ((defined $opt_HHR) && (exists $opt_HHR->{$optname})) ? 1 : 0;
}

#####################################################################
# Subroutine: opt_ValidateSet()
# Incept:     EPN, Thu Oct 30 08:53:54 2014
# 
# Purpose:    Validate the options that have been set by the user
#             by enforcing that none are incompatible and that
#             all required combinations are set.
#
# Arguments:
# $opt_HHR:      ref to %opt_HH (see description at top of this file)
# $opt_order_AR: ref to @opt_order_A (see description at top of this file)
#
# Returns:    void, if we return then all is good
# Dies:       if an incompatible combination is set, or one of 
#             a required combination is not set.
#
####################################################################
sub opt_ValidateSet {
  my $narg_expected = 2;
  my $sub_name = "opt_ValidateSet()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($opt_HHR, $opt_order_AR) = @_;

  my $optname;
  my $optname2;
  my @optname2_A = ();
  my $optname3;
  my @optname3_A = ();

  # before checking the options that the user set, check that all the 
  # options are internally consistent:
  # 1. no required option should be incompatible or required by another option.
  #    (it is required by all other options, so having it required by any one option is silly.)
  # 2. no option should be incompatible with a required option
  # 3. no option should require and be incompatible with the same option
  foreach $optname (@{$opt_order_AR}) { 
    # 1. no required option should be incompatible or required by another option.
    #    (it is required by all other options, so having it required by any one option is silly.)
    if(defined $opt_HHR->{$optname}{"incompatible"}) { 
      @optname2_A = split(",", $opt_HHR->{$optname}{"incompatible"});
      foreach $optname2 (@optname2_A) { 
        if(! exists $opt_HHR->{$optname2}) { 
            printf STDERR ("ERROR, in $sub_name, option $optname2 listed as incompatible with $optname does not exist.\n"); 
            exit(1);
        }
        if($opt_HHR->{$optname2}{"is_required"} == 1) { 
          printf STDERR ("ERROR, in $sub_name, option $optname is incompatible with the required option $optname2, this is not allowed.\n"); 
          exit(1);
        }
      }
    }
    # 2. no option should be incompatible with a required option
    if(defined $opt_HHR->{$optname}{"requires"}) { 
      @optname2_A = split(",", $opt_HHR->{$optname}{"requires"});
      foreach $optname2 (@optname2_A) { 
        if(! exists $opt_HHR->{$optname2}) { 
            printf STDERR ("ERROR, in $sub_name, option $optname2 listed as required with $optname does not exist.\n"); 
            exit(1);
        }
        if($opt_HHR->{$optname2}{"is_required"} == 1) { 
          printf STDERR ("ERROR, in $sub_name, option $optname is listed as required in combination with the required option $optname2, this is not allowed (all options are implicitly required in combination with required options).\n"); 
          exit(1);
        }
      }
    }
    # 3. no option should require and be incompatible with the same option
    if((defined $opt_HHR->{$optname}{"incompatible"}) && 
       (defined $opt_HHR->{$optname}{"requires"})) { 
      @optname2_A = split(",", $opt_HHR->{$optname}{"incompatible"});
      @optname3_A = split(",", $opt_HHR->{$optname}{"requires"});
      foreach $optname2 (@optname2_A) { 
        foreach $optname3 (@optname3_A) { 
          if($optname2 eq $optname3) { 
            printf STDERR ("ERROR, in $sub_name, option $optname2 is listed as required and incompatible with $optname, this is not allowed.\n"); 
            exit(1);
          }
        }
      }
    }
  }
       
  # check that all required options were set
  foreach $optname (@{$opt_order_AR}) { 
    if($opt_HHR->{$optname}{"is_required"} == 1) { 
      if($opt_HHR->{$optname}{"setby"} ne "user") { 
        printf STDERR ("ERROR, REQUIRED option $optname was not used.\n"); 
        exit(1);
      }
    }
  }

  # check for incompatible options
  foreach $optname (@{$opt_order_AR}) { 
    if($opt_HHR->{$optname}{"setby"} eq "user") { 
      if(defined $opt_HHR->{$optname}{"incompatible"}) { 
        @optname2_A = split(",", $opt_HHR->{$optname}{"incompatible"});
        foreach $optname2 (@optname2_A) { 
          if(! exists $opt_HHR->{$optname2}) { 
            printf STDERR ("ERROR, $sub_name option $optname2 listed as incompatible with $optname does not exist.\n"); 
            exit(1);
          }
          if($opt_HHR->{$optname2}{"setby"} eq "user") { 
            printf STDERR ("ERROR, $sub_name option $optname and $optname2 are incompatible, choose one or the other.\n"); 
            exit(1);
          }
        }
      }
    }
  }

  # check required options
  foreach $optname (@{$opt_order_AR}) { 
    if($opt_HHR->{$optname}{"setby"} eq "user") { 
      if(defined $opt_HHR->{$optname}{"requires"}) { 
        @optname2_A = split(",", $opt_HHR->{$optname}{"requires"});
        foreach $optname2 (@optname2_A) { 
          if(! exists $opt_HHR->{$optname2}) { 
            printf STDERR ("ERROR, $sub_name option $optname2 listed as required with $optname does not exist.\n"); 
            exit(1);
          }
          if($opt_HHR->{$optname2}{"setby"} ne "user") { 
            printf STDERR ("ERROR, $sub_name option $optname2 is required in combination with $optname, but is not set.\n"); 
            exit(1);
          }
        }
      }
    }
  }

  return;
}

#####################################################################
# Subroutine: opt_EnabledString()
# Incept:     EPN, Thu Oct 30 09:03:09 2014
# 
# Purpose:    Return a string of all options enabled by the user 
#             with their arguments. The order of the options may
#             not be the same as the order from the cmdline.
#
# Arguments:
# $opt_HHR:      ref to %opt_HH (see description at top of this file)
# $opt_order_AR: ref to @opt_order_A (see description at top of this file)
#
# Returns:    String of all enabled options with their arguments.
#
####################################################################
sub opt_EnabledString {
  my $narg_expected = 2;
  my $sub_name = "opt_EnabledString()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($opt_HHR, $opt_order_AR) = @_;

  my $ret_string = "";
  foreach my $optname (@{$opt_order_AR}) { 
    if($opt_HHR->{$optname}{"setby"} eq "user") { 
      $ret_string .= "$optname ";
      if($opt_HHR->{$optname}{"type"} ne "boolean") { 
        $ret_string .= $opt_HHR->{$optname}{"value"} . " ";
      }
    }
  }
  $ret_string =~ s/\s+$//; # remove trailing spaces

  return $ret_string;
}

#####################################################################
# Subroutine: opt_OutputPreamble()
# Incept:     EPN, Thu Oct 30 09:46:09 2014
# 
# Purpose:    Output preamble section including a list of all
#             enabled options.
#
# Arguments:
# $FH:           file handle to print to (may be *STDOUT)
# $arg_AR:       ref to array of cmdline arguments
# $arg_desc_AR:  ref to array of cmdline argument descriptions
# $opt_HHR:      ref to %opt_HH (see description at top of this file)
# $opt_order_AR: ref to @opt_order_A (see description at top of this file)
#
# Returns:    void
#
####################################################################
sub opt_OutputPreamble { 
  my $narg_expected = 5;
  my $sub_name = "opt_OutputPreamble()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($FH, $arg_AR, $arg_desc_AR, $opt_HHR, $opt_order_AR) = @_;

  my $narg = scalar(@{$arg_AR});
  if(scalar(@{$arg_desc_AR}) != $narg) { printf STDERR ("ERROR, $sub_name entered with different sizes for arg_AR and arg_desc_AR.\n", scalar(@_), $narg_expected); exit(1); } 

  my @lhs_A = ();
  my @rhs_A = ();
  my $nlines = 0;
  my $i;
  for($i = 0; $i < $narg; $i++) { 
    push(@lhs_A,  "# " . $arg_AR->[$i] . ":");
    push(@rhs_A,  $arg_desc_AR->[$i]);
    $nlines++;
  }
  foreach my $optname (@{$opt_order_AR}) {
    if(opt_IsUsed($optname, $opt_HHR)) { 
      if(defined $opt_HHR->{$optname}{"preamble"}) { 
        push(@lhs_A, "# " . $opt_HHR->{$optname}{"preamble"} . ":");
        push(@rhs_A, opt_PreambleRHS($optname, $opt_HHR)); 
        $nlines++;
      }
    }
  }

  my $width_lhs = maxLengthScalarInArray(\@lhs_A);
  #my $width_rhs = maxLengthScalarInArray(\@rhs_A);
  for($i = 0; $i < $nlines; $i++) { 
    printf $FH ("%-*s  %-s\n", $width_lhs, $lhs_A[$i], $rhs_A[$i]); 
  }
  print $FH ("\# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  return;
}

#####################################################################
# Subroutine: opt_OutputHelp()
# Incept:     EPN, Wed Oct 29 06:03:07 2014
# 
# Purpose:    Output help, with list of options with descriptions.
#
# Arguments:
# $FH:                file handle to print to (may be *STDOUT)
# $usage:             name of executable
# $opt_HHR:           ref to %opt_HH (see description at top of this file)
# $opt_order_AR:      ref to @opt_order_A (see description at top of this file)
# $opt_group_desc_HR: ref to hash of group descriptions
#
# Returns:    void
#
####################################################################
sub opt_OutputHelp { 
  my $narg_expected = 5;
  my $sub_name = "opt_OutputHelp()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($FH, $usage, $opt_HHR, $opt_order_AR, $opt_group_desc_HR) = @_;

  my @lhs_A = ();
  my @rhs_A = ();
  my $optname;
  my $i;

  printf $FH $usage;

  # for each group
  foreach my $g (sort {$a <=> $b} keys %{$opt_group_desc_HR}) { 
    print $FH "\n";
    print $FH $opt_group_desc_HR->{$g} . ":\n";
    @lhs_A = ();
    @rhs_A = ();
    # find all options in the group
    my $nopt = 0;
    foreach $optname (@{$opt_order_AR}) { 
      if($opt_HHR->{$optname}{"group"} == $g) { 
        # create string for left hand side, e.g.: "-Z <s>"
        my $lhs_str = $optname;
        if($opt_HHR->{$optname}{"type"} eq "integer") { $lhs_str .= " <n>"; }
        if($opt_HHR->{$optname}{"type"} eq "real")    { $lhs_str .= " <x>"; }
        if($opt_HHR->{$optname}{"type"} eq "string")  { $lhs_str .= " <s>"; }
        push(@lhs_A, $lhs_str);

        # create string for right hand side, e.g.: "set database size as <x> Mb";
        my $rhs_str = $opt_HHR->{$optname}{"help"};
        # add default value unless it's a boolean 
        if($opt_HHR->{$optname}{"type"} ne "boolean" && defined $opt_HHR->{$optname}{"default"}) { 
          $rhs_str .= " [" . $opt_HHR->{$optname}{"default"} . "]";
        }
        push(@rhs_A, $rhs_str);

        $nopt++;
      }
    }
    my $width_lhs = maxLengthScalarInArray(\@lhs_A);
    # output formatted string for each option
    for($i = 0; $i < $nopt; $i++) { 
      printf $FH ("  %-*s : %s\n", $width_lhs, $lhs_A[$i], $rhs_A[$i]);
    }
  }

  return;
}

#####################################################################
# Subroutine: opt_PreambleRHS()
# Incept:     EPN, Thu Oct 30 09:52:35 2014
# 
# Purpose:    Create the 'right hand side' of a preamble line for 
#             for a <$optname>.
#
# Arguments:
# $optname:     option name
# $opt_HHR:     ref to %opt_HH (see description at top of this file)
#
# Returns:    void
# Dies:       if $optname does not exist
#             (This must be a coding error, not a user error)
#
####################################################################
sub opt_PreambleRHS { 
  my $narg_expected = 2;
  my $sub_name = "opt_PreambleRHS()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($optname, $opt_HHR) = @_;

  if(! exists $opt_HHR->{$optname})    { print STDERR ("ERROR, $sub_name for $optname, but it does not exist.\n");     exit(1); } 
  if(! opt_IsUsed($optname, $opt_HHR)) { print STDERR ("ERROR, $sub_name for $optname, but it is not set by user.\n"); exit(1); } 

  my $ret_str = "";
  if($opt_HHR->{$optname}{"type"} eq "boolean") { 
    if($opt_HHR->{$optname}{"value"}    == 0) { $ret_str .= "yes "; }
    elsif($opt_HHR->{$optname}{"value"} == 1) { $ret_str .= "yes "; }
    else { print STDERR ("ERROR, $sub_name for $optname, boolean type but value is %s.\n", $opt_HHR->{$optname}{"value"}); exit(1); } 
  }
  else { # not a boolean
    $ret_str .= $opt_HHR->{$optname}{"value"} . " ";
  }
  $ret_str .= "[$optname]";

  return $ret_str;
}

#####################################################################
# Subroutine: verify_integer()
# Incept:     EPN, Wed Oct 29 10:14:06 2014
# 
# Purpose:  Verify that $x is a (positive or negative) integer. Return
#           '1' if it is, else return '0'.
#
# Arguments:
# $x:         value to check
#
# Returns:    see purpose
####################################################################
sub verify_integer { 
  my $narg_expected = 1;
  my $sub_name = "verify_integer()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($x) = @_;

  return ($x =~ m/^\-?\d+$/) ? 1 : 0;
}

#####################################################################
# Subroutine: verify_real()
# Incept:     EPN, Wed Oct 29 10:15:31 2014
# 
# Purpose:    Verify that $x is a real number. Return '1' if it is,
#             else return '0'.
#
# Arguments:
# $x:         value to check
#
# Returns:    see purpose
####################################################################
sub verify_real { 
  my $narg_expected = 1;
  my $sub_name = "verify_real()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($x) = @_;

  if(! defined $x)                { return 0; }
  if($x eq "")                    { return 0; }
  if($x =~ m/^\-?\d*\.?\d*$/)     { return 1; } # matches: 1, 3.4, 53.43, .3, .333, 0., 155., 150, -1, -3.4, -53.43, -.3, -.333, -0., -155., -150
  if($x =~ m/^\-?\d+[Ee]\-?\d*$/) { return 1; } # matches: 3E5, -3E5, -3E-5, 3e5, -3e5, -3e-5
  return 0;
}

#################################################################
# Subroutine : maxLengthScalarInArray()
# Incept:      EPN, Tue Nov  4 15:19:44 2008 [ssu-align]
# 
# Purpose:     Return the maximum length of a scalar in an array
#
# Arguments: 
#   $AR: reference to the array
# 
# Returns:     The length of the maximum length scalar.
#
################################################################# 
sub maxLengthScalarInArray {
  my $nargs_expected = 1;
  my $sub_name = "maxLengthScalarInArray()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR) = $_[0];

  my $max = 0;
  my $len = 0;
  foreach my $el (@{$AR}) { 
    $len = length($el);
    if($len > $max) { $max = $len; }
  }
  return $max;
}

####################################################################
# the next line is critical, a perl module must return a true value
return 1;
####################################################################
