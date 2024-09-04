#!/usr/bin/perl
#
# sqp_utils.pm
# Eric Nawrocki
# EPN, Tue Mar 19 13:35:06 2019 [incept, in vadr]
# EPN, Tue Jul  2 11:53:49 2019 [migrated from vadr's epn-utils.pm (as of commit 69b003d)]]
# version: 0.10
#
use strict;
use warnings;
use Time::HiRes qw(gettimeofday);

# NOTE: do not add any 'require' statements here, e.g. 'require
# sqp_utils.pm' because the program that uses sequip must handle that
# so each program can specify sequip from a specific directory defined
# by a specific environment variable. This is how, for example,
# ribovore can require a specific version of sequip on the same file
# system that has vadr installed with a potentially different version
# of sequip.

#################################################################
# Subroutine: utl_RunCommand()
# Incept:     EPN, Thu Feb 11 13:32:34 2016 [dnaorg.pm]
#
# Purpose:     Runs a command using system() and exits in error 
#              if the command fails. If $be_verbose, outputs
#              the command to stdout. If $FH_HR->{"cmd"} is
#              defined, outputs command to that file handle.
#
# Arguments:
#   $cmd:         command to run, with a "system" command;
#   $be_verbose:  '1' to output command to stdout before we run it, '0' not to
#   $do_failok:   '1' to NOT exit if command fails, '0' to exit if command fails
#   $FH_HR:       REF to hash of file handles, including "cmd"
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if $cmd fails and $do_failok is '0'
#################################################################
sub utl_RunCommand {
  my $sub_name = "utl_RunCommand()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $be_verbose, $do_failok, $FH_HR) = (@_);
  
  my $cmd_FH = undef;
  if(defined $FH_HR && defined $FH_HR->{"cmd"}) { 
    $cmd_FH = $FH_HR->{"cmd"};
  }

  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }

  if(defined $cmd_FH) { 
    print $cmd_FH ("$cmd\n");
  }

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));

  system($cmd);

  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));

  if(($? != 0) && (! $do_failok)) { 
    ofile_FAIL("ERROR in $sub_name, the following command failed:\n$cmd\n", $?, $FH_HR); 
  }

  return ($stop_time - $start_time);
}

#################################################################
# Subroutine: utl_ArrayOfHashesToArray()
# Incept:     EPN, Wed Mar 20 09:07:06 2019
#
# Purpose:    Fill @{$AR} with all values in $AHR->[]{$key}.
# Arguments:
#   $AHR:      REF to array of hashes
#   $AR:       REF to array to add to
#   $key:      key of interest
# 
# Returns: number of elements added to @{$AR}
#
#################################################################
sub utl_ArrayOfHashesToArray {
  my $sub_name = "utl_ArrayOfHashesToArray()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($AHR, $AR, $key) = (@_);

  my $ret_n = 0;
  my $n = scalar(@{$AHR});
  for(my $i = 0; $i < $n; $i++) { 
    if(defined $AHR->[$i]{$key}) { 
      push(@{$AR}, $AHR->[$i]{$key}); 
      $ret_n++;
    }
  }
  
  return $ret_n;
}

#################################################################
# Subroutine:  utl_ConcatenateListOfFiles()
# Incept:      EPN, Sun Apr 24 08:08:15 2016
#
# Purpose:     Concatenate a list of files into one file.  If the list is
#              more than 20K characters, break it down into multiple 
#              cat calls of at most 20K characters each, and call this 
#              subroutine recursively to concanenate them.
#
#              To allow recursive calls without clobbering files
#              created in previous calls, we use a unique string to
#              append to the output file names. This is derived from
#              the $caller_sub_name array. If that is undef, we use
#              'rec0' as the unique string, if it ends with 'rec<d>'
#              we know we have a recursive call and use 'rec<d+1>' in
#              the recursive call.
# 
#              We remove all files that we concatenate unless
#              --keep option is on in %{$opt_HHR}.
#
# Arguments: 
#   $file_AR:          REF to array of all files to concatenate
#   $outfile:          name of output file to create by concatenating
#                      all files in @{$file_AR}.
#   $caller_sub_name:  name of calling subroutine (can be undef)
#   $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:            ref to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If one of the cat commands fails.
#              If $outfile is in @{$file_AR}
#
################################################################# 
sub utl_ConcatenateListOfFiles { 
  my $nargs_expected = 5;
  my $sub_name = "utl_ConcatenateListOfFiles()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($file_AR, $outfile, $caller_sub_name, $opt_HHR, $FH_HR) = (@_);

  if((! defined $file_AR) || (scalar(@{$file_AR}) == 0)) { 
    ofile_FAIL(sprintf("ERROR in $sub_name%s, array of file names to concatenate is undefined or empty", 
                        (defined $caller_sub_name) ? " called by $caller_sub_name" : ""), 1, $FH_HR);
  }

  if(utl_AFindNonNumericValue($file_AR, $outfile, $FH_HR) != -1) { 
    ofile_FAIL(sprintf("ERROR in $sub_name%s, output file name $outfile exists in list of files to concatenate", 
                        (defined $caller_sub_name) ? " called by $caller_sub_name" : ""), 1, $FH_HR);
  }

  # determine unique string for tmp naming files that allows
  # proper handling of recursive calls without namespace clashes
  my $rec_idx = 0;
  if((defined $caller_sub_name) && ($caller_sub_name =~ /rec(\d+)$/)) { 
    # caller was this subroutine, using $1 as the $rec_idx
    $rec_idx = $1 + 1;
    # e.g. if $caller_sub_name = 'utl_ConcatenateListOfFiles.rec2', then $rec_key set to 3
  }

  my $nchar_limit = 20000; # 20K characters, hard-coded
  my $tmp_outfile = undef; # name of temporary outfile, only used if needed
  my $tot_nfiles = scalar(@{$file_AR}); # total number of files in @file_A
  my $cur_nfiles = 0; # current number of files in cat call
  my $i = 0;          # index of file in @{$file_AR}

  my $outfile_len = length($outfile); # length of output file name
  my $cat_cmd = "cat ";
  my $tmp_idx = 1;
  my @tmp_outfile_A = (); # array of temporary files we created to avoid too long of a command

  for($i = 0; $i < $tot_nfiles; $i++) { 
    $cat_cmd .= $file_AR->[$i] . " ";
    $cur_nfiles++;
    if(($cur_nfiles >= 2) && # we're cat'ing at least 2 files
       (length($cat_cmd) + 2 + $outfile_len + 4 + length ($rec_idx) + 1 + length($tmp_idx)) > $nchar_limit) { # length of command exceeds $nchar_limit (2 = length("> "), 4 = length(".rec")), 1 = length(".")
      # finish the cat command and execute it to create a tmp file we'll cat in a recursive call
      $tmp_outfile = $outfile . ".rec" . $rec_idx . "." . $tmp_idx; 
      $cat_cmd .= "> $tmp_outfile";
      
      # execute the command
      utl_RunCommand($cat_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
      push(@tmp_outfile_A, $tmp_outfile);

      $tmp_idx++;
      $cat_cmd = "cat ";
      $cur_nfiles = 0;
    }
  }
  if(scalar(@tmp_outfile_A) == 0) { 
    # we did not create any temporary files
    # finish the command and execute it
    $cat_cmd .= "> $outfile";
    utl_RunCommand($cat_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
  }    
  else { 
    # we did create at least one temporary file, we need to call this subroutine recursively
    # first finish the final tmp file if nec
    if($cur_nfiles > 0) { 
      $tmp_outfile = $outfile . ".rec" . $rec_idx . "." . $tmp_idx; 
      $cat_cmd .= "> $tmp_outfile";
      
      # execute the command
      utl_RunCommand($cat_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
      push(@tmp_outfile_A, $tmp_outfile);
    }
    my $tmp_caller_str = $sub_name . ".rec" . $rec_idx;
    utl_ConcatenateListOfFiles(\@tmp_outfile_A, $outfile, (defined $caller_sub_name) ? $caller_sub_name . ":" . $tmp_caller_str : $tmp_caller_str, $opt_HHR, $FH_HR);
  }

  # remove all original files (recursive call(s) will remove any tmp files we created)
  if(! opt_Get("--keep", $opt_HHR)) { 
    utl_FileRemoveList($file_AR, $sub_name, $opt_HHR, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine:  utl_AFindNonNumericValue()
# Incept:      EPN, Tue Feb 16 10:40:57 2016 [dnaorg.pm]
#
# Purpose:     Returns (first) index in @{$AR} that has the 
#              nonnumeric value $value. Returns -1 
#              if it does not exist.
#
# Arguments: 
#   $AR:       REF to array 
#   $value:    the value we're checking exists in @{$AR}
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
# 
# Returns:     index ($i) '1' if $value exists in @{$AR}, '-1' if not
#
# Dies:        if $value is numeric, or @{$AR} is not defined.
################################################################# 
sub utl_AFindNonNumericValue { 
  my $nargs_expected = 3;
  my $sub_name = "utl_AFindNonNumericValue()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($AR, $value, $FH_HR) = (@_);

  if(utl_IsReal($value)) { 
    ofile_FAIL("ERROR in $sub_name, value $value seems to be numeric, we can't compare it for equality", 1, $FH_HR);
  }

  if(! defined $AR) { 
    ofile_FAIL("ERROR in $sub_name, array reference is not defined", 1, $FH_HR);
  }

  for(my $i = 0; $i < scalar(@{$AR}); $i++) {
    if($AR->[$i] eq $value) { 
      return $i; 
    }
  }

  return -1; # did not find it
}

#################################################################
# Subroutine: utl_AFindValue()
# Incept:     EPN, Tue Mar  8 11:26:03 2016
# Synopsis:   Look for a value in an array and return the index
#             of it, if found, else return -1. If it exists more than
#             once, return the minimum index.
#
# Arguments:
#  $value:   value to look for
#  $AR:      array to look in
#
# Returns:     index ($i) '1' if $value exists in @{$AR}, '-1' if not
#
# Dies:        if $value is numeric, or @{$AR} is not defined.
# 
#################################################################
sub utl_AFindValue { 
  my $sub_name = "utl_AFindValue()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($value, $AR, $FH_HR) = @_;

  if(! defined $AR) { 
    ofile_FAIL("ERROR in $sub_name, array reference is not defined", 1, $FH_HR);
  }

  if(utl_IsReal($value)) { # compare with ==
    for(my $i = 0; $i < scalar(@{$AR}); $i++) { 
      my $el = $AR->[$i];
      if(utl_IsReal($el) && ($value == $el)) { 
        return $i;
      }
    }
  }
  else { # compare with 'eq'
    for(my $i = 0; $i < scalar(@{$AR}); $i++) { 
      my $el = $AR->[$i];
      if((! utl_IsReal($el)) && ($value eq $el)) { 
        return $i;
      }
    }
  }
  return -1;
}  

#################################################################
# Subroutine:  utl_ACountNonNumericValue()
# Incept:      EPN, Fri Mar 11 06:34:51 2016
#
# Purpose:     Returns number of times nonnumeric value 
#              $value exists in @{$AR}. Returns 0 if
#              it doesn't exist.
#
# Arguments: 
#   $AR:       REF to array 
#   $value:    the value we're looking for in @{$AR}
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
# 
# Returns:     Number of occurrences of $value in @{$AR}.
#
# Dies:        if $value is numeric, or @{$AR} is not defined.
#
################################################################# 
sub utl_ACountNonNumericValue { 
  my $nargs_expected = 3;
  my $sub_name = "utl_ACountNonNumericValue()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($AR, $value, $FH_HR) = @_;

  if(utl_IsReal($value)) { 
    ofile_FAIL("ERROR in $sub_name, value $value seems to be numeric, we can't compare it for equality", 1, $FH_HR);
  }

  if(! defined $AR) { 
    ofile_FAIL("ERROR in $sub_name, array reference is not defined", 1, $FH_HR);
  }

  my $ct = 0;
  for(my $i = 0; $i < scalar(@{$AR}); $i++) {
    if($AR->[$i] eq $value) { 
      $ct++;
    }
  }

  return $ct;
}

#################################################################
# Subroutine: utl_FileRemoveUsingSystemRm
# Incept:     EPN, Fri Mar  4 15:57:25 2016 [dnaorg.pm]
#
# Purpose:    Remove a file from the filesystem by using
#             the system rm command.
# Arguments:
#   $file:            file to remove
#   $caller_sub_name: name of caller, can be undef
#   $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    - if the file does not exist
#
#################################################################
sub utl_FileRemoveUsingSystemRm {
  my $sub_name = "utl_FileRemoveUsingSystemRm";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($file, $caller_sub_name, $opt_HHR, $FH_HR) = (@_);
  
  if(! -e $file) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, %s trying to remove file $file but it does not exist", 
                (defined $caller_sub_name) ? "called by $caller_sub_name," : 0), 1, $FH_HR); 
  }

  utl_RunCommand("rm $file", opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}

#################################################################
# Subroutine:  utl_RemoveDirPath()
# Incept:      EPN, Mon Nov  9 14:30:59 2009 [ssu-align]
#
# Purpose:     Given a full path of a file remove the directory path.
#              For example: "foodir/foodir2/foo.stk" becomes "foo.stk".
#
# Arguments: 
#   $fullpath: name of original file
# 
# Returns:     The string $fullpath with dir path removed.
#
################################################################# 
sub utl_RemoveDirPath {
  my $sub_name = "utl_RemoveDirPath()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($fullpath) = (@_);

  $fullpath =~ s/^.+\///;

  return $fullpath;
}

#################################################################
# Subroutine:  utl_HHFromAH()
# Incept:      EPN, Fri Mar 22 09:43:34 2019
#
# Purpose:     Create a 2D hash %{$HHR} from an array of hashes
#              @{$AHR}.
#              First dim keys in %{$HHR} will be values from 
#              $AHR->[]{$AH_key_for_HH_key}.
#              Second dim keys in %{$HHR} will be all other keys from
#              $AHR->[]{}.
#              
# Arguments: 
#   $HHR:               ref to 2D hash to create
#   $AHR:               ref to array of hashes to copy from
#   $AH_key_for_HH_key: key in @{$AHR} to get to use value from to use
#                       as key in 1st dim of %{$HHR}
#   $call_str:          string describing caller, to output if we die
#   $FH_HR:             ref to hash of file handles, including "log" and "cmd"
# 
# Returns:     void
# 
# Dies:        If not all elements of @{$AHR} have 
#              $AHR->[]{$AH_key_for_HH_key} defined. 
#              If more than one elements of @{$AHR} have same
#              value for @{$AHR->[]{$key}}.
#
################################################################# 
sub utl_HHFromAH {
  my $sub_name = "utl_HHFromAH()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($HHR, $AHR, $AH_key_for_HH_key, $call_str, $FH_HR) = (@_);

  %{$HHR} = ();
  my $n = scalar(@{$AHR});
  my $AH_value_for_HH_key = undef;
  for(my $i = 0; $i < $n; $i++) { 
    if(! defined $AHR->[$i]{$AH_key_for_HH_key}) {
      ofile_FAIL("ERROR in $sub_name,%selement $i does not have key $AH_key_for_HH_key", 
                 (defined $call_str) ? "$call_str" : "", 1, $FH_HR); 
    }
    $AH_value_for_HH_key = $AHR->[$i]{$AH_key_for_HH_key};
    if(defined $HHR->{$AH_value_for_HH_key}) { 
      ofile_FAIL("ERROR in $sub_name,%stwo elements have same value for $AH_key_for_HH_key key ($AH_value_for_HH_key)", 
                 (defined $call_str) ? "$call_str" : "", 1, $FH_HR); 
    }
    %{$HHR->{$AH_value_for_HH_key}} = ();
    foreach my $AH_key2_for_HH_key (keys (%{$AHR->[$i]})) { 
      if($AH_key2_for_HH_key ne $AH_key_for_HH_key) { 
        $HHR->{$AH_value_for_HH_key}{$AH_key2_for_HH_key} = $AHR->[$i]{$AH_key2_for_HH_key};
      }
    }
  }

  return;
}

#################################################################
# Subroutine:  utl_HHFromAHAddIdx()
# Incept:      EPN, Thu Mar 21 06:35:28 2019
#
# Purpose:     Create a 2D hash %{$HHR} from an array of hashes
#              @{$AHR} by calling utL_HHFromAH() 
#              (see that sub's Purpose for more details)
#              
#              And then add $HHR->{$key}{"idx"}, that gives index <i> of 
#              @{$AHR} for which $AHR->[<i>]{$AH_key_for_HH_key} == $key.
#              
# Arguments: 
#   $HHR:               ref to 2D hash to create
#   $AHR:               ref to array of hashes to copy from
#   $AH_key_for_HH_key: key in @{$AHR} to get to use value from to use
#                       as key in 1st dim of %{$HHR}
#   $call_str:          string describing caller, to output if we die
#   $FH_HR:             ref to hash of file handles, including "log" and "cmd"
# 
# Returns:     void
# 
# Dies:        If not all elements of @{$AHR} have 
#              $AHR->[]{$AH_key_for_HH_key} defined. 
#              If more than one elements of @{$AHR} have same
#              value for @{$AHR->[]{$key}}.
#              If $AHR->[]{"idx"} is defined for any element.
#              If $AH_key_for_HH_key is "idx";
#
################################################################# 
sub utl_HHFromAHAddIdx {
  my $sub_name = "utl_HHFromAHAddIdx()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($HHR, $AHR, $AH_key_for_HH_key, $call_str, $FH_HR) = (@_);

  # make sure "idx" is not the $AH_key_for_HH_key
  if($AH_key_for_HH_key eq "idx") { 
    ofile_FAIL("ERROR in $sub_name,%skey to choose is \"idx\"", 
               (defined $call_str) ? "$call_str" : "", 1, $FH_HR); 
  }
  # make sure "idx" 2D key does not exist for any element
  my $n = scalar(@{$AHR});
  for(my $i = 0; $i < $n; $i++) { 
    if(defined $AHR->[$i]{"idx"}) {
      ofile_FAIL("ERROR in $sub_name,%selement $i already has key \"idx\" upon entry", 
                 (defined $call_str) ? "$call_str" : "", 1, $FH_HR); 
    }
  }

  # do most of the work with utl_HHFromAH
  utl_HHFromAH($HHR, $AHR, $AH_key_for_HH_key, $call_str, $FH_HR);

  # add idx 
  my $AH_value_for_HH_key = undef;
  for(my $i = 0; $i < $n; $i++) { 
    $AH_value_for_HH_key = $AHR->[$i]{$AH_key_for_HH_key};
    $HHR->{$AH_value_for_HH_key}{"idx"} = $i;
  }

  return;
}

#################################################################
# Subroutine:  utl_HFromAH()
# Incept:      EPN, Fri Mar 22 06:20:18 2019
#
# Purpose:     Create a 1D hash %{$HR} using key/value pairs
#              from @{$AHR}. 
#              Keys in %{$HR} will be values from 
#              $AHR->[]{$AH_key_for_H_key}.
#              Values in %{$HR} will be values from 
#              $AHR->[]{$AH_key_for_H_value}.
#              
# Arguments: 
#   $HR:                  ref to 1D hash to create
#   $AHR:                 ref to array of hashes to copy from
#   $AH_key_for_H_key:    key in @{$AHR} to get value from to use as key in %{$HR}
#   $AH_key_for_H_value:  key in @{$AHR} to get value from to use as value in %{$HR}
#   $call_str:            string describing caller, to output if we die
#   $FH_HR:               ref to hash of file handles, including "log" and "cmd"
# 
# Returns:     void
# 
# Dies:        If not all elements of @{$AHR} have 
#              $AHR->[]{$AH_key_for_H_key}}
#              If not all elements of @{$AHR} have 
#              $AHR->[]{$AH_key_for_H_value}}
#              If more than one elements of @{$AHR} have same
#              value for @{$AHR->[]{$AH_key_for_H_key}}.
#
################################################################# 
sub utl_HFromAH {
  my $sub_name = "utl_HFromAH()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($HR, $AHR, $AH_key_for_H_key, $AH_key_for_H_value, $call_str, $FH_HR) = (@_);

  printf("in $sub_name, call_str: $call_str\n");

  %{$HR} = ();
  my $n = scalar(@{$AHR});
  my $AH_value_for_H_key = undef;
  for(my $i = 0; $i < $n; $i++) { 
    if(! defined $AHR->[$i]{$AH_key_for_H_key}) { 
      ofile_FAIL("ERROR in $sub_name,%selement $i does not have key $AH_key_for_H_key", 
                 (defined $call_str) ? "$call_str" : "", 1, $FH_HR); 
    }
    if(! defined $AHR->[$i]{$AH_key_for_H_value}) { 
      ofile_FAIL("ERROR in $sub_name,%selement $i does not have key $AH_key_for_H_value", 
                 (defined $call_str) ? "$call_str" : "", 1, $FH_HR); 
    }
    my $H_key   = $AHR->[$i]{$AH_key_for_H_key};
    my $H_value = $AHR->[$i]{$AH_key_for_H_value};
    if(defined $HR->{$H_key}) { 
      ofile_FAIL("ERROR in $sub_name,%stwo elements have same value for key $AH_key_for_H_key ($H_key)", 
                 (defined $call_str) ? "$call_str" : "", 1, $FH_HR); 
    }
    $HR->{$H_key} = $H_value;
  }

  return;
}

#################################################################
# Subroutine:  utl_IdxHFromA()
# Incept:      EPN, Fri Mar 22 10:01:38 2019
#
# Purpose:     Create a 1D 'index' hash %{$idx_HR} such that
#              $idx_HR->{$key} == $i if
#              $AR->[$i] == $key
#              from @{$AHR}. 
#              
# Arguments: 
#   $idx_HR:              ref to 1D hash to create
#   $AR:                  ref to array
#   $call_str:            string describing caller, to output if we die
#   $FH_HR:               ref to hash of file handles, including "log" and "cmd"
# 
# Returns:     void
# 
# Dies:        Two elements of @{$AR} are identical:
#              (if AR->[$i] == AR->[$j] and $i != $j for any $i, $j)
#
################################################################# 
sub utl_IdxHFromA {
  my $sub_name = "utl_IdxHFromA()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($HR, $AR, $call_str, $FH_HR) = (@_);

  %{$HR} = ();
  my $n = scalar(@{$AR});
  for(my $i = 0; $i < $n; $i++) { 
    my $key = $AR->[$i];
    if(defined $HR->{$key}) { 
      # should I check here that $key is not a number? 
      ofile_FAIL("ERROR in $sub_name,%stwo elements in array have same value ($key)",
                 (defined $call_str) ? "$call_str" : "", 1, $FH_HR); 
    }
    $HR->{$key} = $i;
  }

  return;
}

#################################################################
# Subroutine:  utl_IdxHFromAH()
# Incept:      EPN, Fri Mar 22 10:01:38 2019
#
# Purpose:     Create a 1D 'index' hash %{$idx_HR} from @{$AHR} 
#              such that
#              $idx_HR->{$key} == $i if
#              $AHR->[$i]{$AH_key} == $key
#              
# Arguments: 
#   $idx_HR:              ref to 1D hash to create
#   $AHR:                 ref to hash of arrays
#   $AH_key:              ref to key in %{$AHR->[$i]}
#   $call_str:            string describing caller, to output if we die
#   $FH_HR:               ref to hash of file handles, including "log" and "cmd"
# 
# Returns:     void
# 
# Dies:        Two values we are trying to add as keys to %{$idx_HR} are identical
#              (if AHR->[$i]{$AH_key} == AHR->[$j]{$AH_key} and $i != $j for any $i, $j)
#              If AHR->[$i]{$AH_key} is undefined for any $i
################################################################# 
sub utl_IdxHFromAH {
  my $sub_name = "utl_IdxHFromAH()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($HR, $AHR, $AH_key, $call_str, $FH_HR) = (@_);

  %{$HR} = ();
  my $n = scalar(@{$AHR});
  for(my $i = 0; $i < $n; $i++) { 
    if(! defined $AHR->[$i]{$AH_key}) { 
      ofile_FAIL("ERROR in $sub_name,%shash that is array element $i does not have key $AH_key",
                 (defined $call_str) ? "$call_str" : "", 1, $FH_HR); 
    }
    my $H_key = $AHR->[$i]{$AH_key};
    if(defined $HR->{$H_key}) { 
      ofile_FAIL("ERROR in $sub_name,%stwo elements of source array of hashes we are trying to use as keys in destination hash have same value ($H_key)",
                 (defined $call_str) ? "$call_str" : "", 1, $FH_HR); 
    }
    $HR->{$H_key} = $i;
  }

  return;
}

#################################################################
# Subroutine: utl_AHCountKeyValue()
# Incept:     EPN, Tue Mar 19 11:37:32 2019
#
# Synopsis: Return the number of elements in @{AHR} that 
#           have a key $key in %{$AHR->[]} with value $value.
#
# Arguments:
#  $AHR:      ref to array of hashes
#  $key:      hash key
#  $value:    hash value
#
# Returns:    Number of array elements for which $AHR->[]{$key} eq $value
#
#################################################################
sub utl_AHCountKeyValue {
  my $sub_name = "utl_AHCountKeyValue";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($AHR, $key, $value) = (@_);

  my $ret_n = 0;
  for(my $i = 0; $i < scalar(@{$AHR}); $i++) { 
    # printf("in $sub_name: AHR->[$i]{$key} $AHR->[$i]{$key}\n");
    if((defined $AHR->[$i]{$key}) && 
       ($AHR->[$i]{$key} eq $value)) { 
      $ret_n++;
    }
  }

  # printf("in $sub_name: returning $ret_n\n");
  return $ret_n;
}

#################################################################
# Subroutine: utl_AHValidate()
# Incept:     EPN, Wed Mar 13 13:24:38 2019
#
# Purpose:    Validate an array of hashes, by making sure it
#             includes a key/value for all keys in @{$keys_AR}.
# Arguments:
#   $AHR:      REF to array of hashes to validate
#   $keys_AR:  REF to array of keys that may be excluded from the hash
#   $fail_str: extra string to output if we die
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
# 
# Returns: scalar(@{$AHR});
#
# Dies:    - if one of the keys in @{$keys_AR} does not exist in all hashes
#            of the array
#
#################################################################
sub utl_AHValidate {
  my $sub_name = "utl_AHValidate()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($AHR, $keys_AR, $fail_str, $FH_HR) = (@_);

  my $n = scalar(@{$AHR});

  for(my $i = 0; $i < $n; $i++) { 
    utl_HValidate($AHR->[$i], $keys_AR, $fail_str, $FH_HR);
  }

  return $n;
}

#################################################################
# Subroutine: utl_HValidate()
# Incept:     EPN, Fri Mar 15 09:37:19 2019
#
# Purpose:    Validate a hash, by making sure a defined value
#             exists for each key in @{$keys_AR}.
# Arguments:
#   $HR:       REF to the hash
#   $keys_AR:  REF to array of keys that may be excluded from the hash
#   $fail_str: extra string to output if we die
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    - if one of the keys in @{$keys_AR} does not exist in the hash
#            of the array
#
#################################################################
sub utl_HValidate {
  my $sub_name = "utl_HValidate()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($HR, $keys_AR, $fail_str, $FH_HR) = (@_);

  foreach my $key (@{$keys_AR}) { 
    if(! exists $HR->{$key}) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, required hash key $key does not exist\n%s", (defined $fail_str) ? $fail_str : ""), 1, $FH_HR); 
    }
    if(! defined $HR->{$key}) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, required hash key $key exists but its value is undefined\n%s", (defined $fail_str) ? $fail_str : ""), 1, $FH_HR); 
    }
  }

  return;
}

#################################################################
# Subroutine:  utl_HMaxLengthKey()
# Incept:      EPN, Thu Dec 13 15:52:09 2018
# 
# Purpose:     Return the maximum length of a scalar key
#              in a hash.
#
# Arguments: 
#   $HR: reference to the hash
# 
# Returns:     The length of the maximum length scalar key.
#
################################################################# 
sub utl_HMaxLengthKey { 
  my $nargs_expected = 1;
  my $sub_name = "utl_HMaxLengthKey()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HR) = $_[0];

  my $max = 0;
  my $len = 0;
  foreach my $key (keys (%{$HR})) { 
    $len = length($key);
    if($len > $max) { $max = $len; }
  }
  return $max;
}

#################################################################
# Subroutine:  utl_HMaxLengthValue()
# Incept:      EPN, Mon Nov  3 09:09:59 2014 [rnavore]
# 
# Purpose:     Return the maximum length of a scalar value
#              in a hash.
#
# Arguments: 
#   $HR: reference to the hash
# 
# Returns:     The length of the maximum length scalar.
#
################################################################# 
sub utl_HMaxLengthValue { 
  my $nargs_expected = 1;
  my $sub_name = "utl_HMaxLengthValue()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HR) = $_[0];

  my $max = 0;
  my $len = 0;
  foreach my $key (keys (%{$HR})) { 
    $len = length($HR->{$key});
    if($len > $max) { $max = $len; }
  }
  return $max;
}


#################################################################
# Subroutine:  utl_HHMaxLengthValueGiven2DKey()
# Incept:      EPN, Fri Mar 29 12:01:18 2019
# 
# Purpose:     Return the maximum length of a scalar value
#              in a 2D hash for a given 1D key hash.
#              max(length($HHR->{}{$key}))
#
# Arguments: 
#   $HHR:  reference to the hash of hashes
#   $key2: 2D key
# 
# Returns: The length of the maximum length scalar.
#
################################################################# 
sub utl_HHMaxLengthValueGiven2DKey { 
  my $nargs_expected = 2;
  my $sub_name = "utl_HHMaxLengthValueGiven2DKey()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HHR, $key2) = (@_);

  my $max = 0;
  foreach my $key1 (keys (%{$HHR})) { 
    if(defined $HHR->{$key1}{$key2}) { 
      $max = utl_Max($max, length($HHR->{$key1}{$key2}));
    }
  }
  return $max;
}

#################################################################
# Subroutine:  utl_AMaxLengthValue()
# Incept:      EPN, Thu Mar 17 12:38:53 2016
# 
# Purpose:     Return the maximum length of a scalar value
#              in an array.
#
# Arguments: 
#   $AR: reference to the array
# 
# Returns:     The length of the maximum length scalar.
#
################################################################# 
sub utl_AMaxLengthValue { 
  my $nargs_expected = 1;
  my $sub_name = "utl_AMaxLengthValue()";
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

#################################################################
# Subroutine:  utl_AArgMax()
# Incept:      EPN, Fri Jan 24 15:09:46 2020
# 
# Purpose:     Return the index of the maximum value numeric 
#              element in an array.
#
# Arguments: 
#   $AR: reference to the array
# 
# Returns:     The index of the maximum value.
#
################################################################# 
sub utl_AArgMax { 
  my $nargs_expected = 1;
  my $sub_name = "utl_AArgMax()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR) = $_[0];

  my $argmax = undef;
  my $max    = undef;
  my $nel = scalar(@{$AR});
  if($nel >= 1) { 
    $argmax = 0;
    $max    = $AR->[0];
  }
  for(my $i = 1; $i < $nel; $i++) { 
    if($AR->[$i] > $max) { 
      $argmax = $i; 
      $max    = $AR->[$i];
    }
  }

  return $argmax;
}

#################################################################
# Subroutine:  utl_NumberOfDigits()
# Incept:      EPN, Tue May  9 11:33:50 2017 [ribovore]
#              EPN, Fri Nov 13 06:17:25 2009 [ssu-align:ssu.pm:NumberOfDigits()]
# 
# Purpose:     Return the number of digits in a number before
#              the decimal point. (ex: 1234.56 would return 4).
# Arguments:
# $num:        the number
# 
# Returns:     the number of digits before the decimal point
#
################################################################# 
sub utl_NumberOfDigits { 
  my $nargs_expected = 1;
  my $sub_name = "utl_NumberOfDigits()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($num) = (@_);
  
  my $ndig = 1; 
  while($num >= 10) { $ndig++; $num /= 10.; }
  
  return $ndig;
}

#################################################################
# Subroutine:  utl_Max()
# Incept:      EPN, Tue Mar 26 11:52:54 2019
# 
# Purpose:     Returns the maximum of $x and $y.
# Arguments:
# $x:          first number
# $y:          second number
# 
# Returns:     maximum of $x and $y
#
################################################################# 
sub utl_Max { 
  my $nargs_expected = 2;
  my $sub_name = "utl_Max()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($x, $y) = (@_);
  
  return ($x > $y) ? $x : $y;
}

#################################################################
# Subroutine:  utl_Min()
# Incept:      EPN, Tue Mar 26 11:54:09 2019
# 
# Purpose:     Returns the minimum of $x and $y.
# Arguments:
# $x:          first number
# $y:          second number
# 
# Returns:     minimum of $x and $y
#
################################################################# 
sub utl_Min { 
  my $nargs_expected = 2;
  my $sub_name = "utl_Min()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($x, $y) = (@_);
  
  return ($x < $y) ? $x : $y;
}

#################################################################
# Subroutine:  utl_Swap()
# Incept:      EPN, Wed Apr  3 06:37:38 2019
# 
# Purpose:     Swaps $$xR and $$yR in place.
# Arguments:
# $xR:         ref to first scalar
# $yR:         ref to second scalar
# 
# Returns:     void
#
################################################################# 
sub utl_Swap { 
  my $nargs_expected = 2;
  my $sub_name = "utl_Swap()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($xR, $yR) = (@_);
  
  my $tmp = $$xR;
  $$xR = $$yR;
  $$yR = $tmp;
  
  return;
}

#################################################################
# Subroutine: utl_ADump()
# Incept:     EPN, Fri Apr  5 12:12:16 2019
#
# Purpose:    Dump the contents of an array,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $AR:          ref of the array
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub utl_ADump { 
  my $sub_name = "utl_ADump()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $AR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");
  
  my $nel = scalar(@{$AR});
  my $undef2print = "!SEQUIP:undef!";
  for(my $i = 0; $i < $nel; $i++) { 
    printf $FH ("\tA[$i]: %s\n", (defined $AR->[$i]) ? $AR->[$i] : $undef2print);
  }

  return;
}

#################################################################
# Subroutine: utl_HDump()
# Incept:     EPN, Thu Apr  4 06:12:39 2019
#
# Purpose:    Dump the contents of a hash,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $HR:          ref of the hash
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub utl_HDump { 
  my $sub_name = "utl_HDump()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $HR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");
  
  my $undef2print = "!SEQUIP:undef!";
  foreach my $key (sort keys %{$HR}) { 
    printf $FH ("\tH{$key}: %s\n", (defined $HR->{$key}) ? $HR->{$key} : $undef2print);
  }

  return;
}

#################################################################
# Subroutine: utl_AADump()
# Incept:     EPN, Fri Apr  5 10:09:57 2019
#
# Purpose:    Dump the contents of a 2D array
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $AAR:         ref of the array of arrays
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub utl_AADump { 
  my $sub_name = "utl_AADump()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $AAR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $undef2print = "!SEQUIP:undef!";
  my $nel1 = scalar(@{$AAR});
  for(my $i1 = 0; $i1 < $nel1; $i1++) { 
    if(! defined $AAR->[$i1]) { 
      print $FH ("*A*A[$i1]: $undef2print\n");
    }
    else { 
      printf $FH ("*A*A[$i1]\n");
      my $nel2 = scalar(@{$AAR->[$i1]});
      for(my $i2 = 0; $i2 < $nel2; $i2++) { 
        printf $FH ("\tA*A*[$i1][$i2]: %s\n", (defined $AAR->[$i1][$i2]) ? $AAR->[$i1][$i2] : $undef2print);
      }
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: utl_HHDump()
# Incept:     EPN, Thu Dec 20 13:36:00 2018
#
# Purpose:    Dump the contents of  hashes of hashes,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $HHR:         ref of the hash of hashes
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub utl_HHDump { 
  my $sub_name = "utl_HHDump()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $HHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");
  
  my $undef2print = "!SEQUIP:undef!";
  foreach my $key1 (sort keys %{$HHR}) { 
    if(! defined $HHR->{$key1}) { 
      print $FH ("*H*H{$key1}: $undef2print\n");
    }
    else { 
      printf $FH ("*H*H{$key1}\n");
      my $nel = scalar(keys %{$HHR->{$key1}});
      foreach my $key2 (sort keys %{$HHR->{$key1}}) { 
        printf $FH ("\tH*H*{$key1}{$key2}: %s\n", (defined $HHR->{$key1}{$key2}) ? $HHR->{$key1}{$key2} : $undef2print);
      }
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: utl_AHHDump()
# Incept:     EPN, Fri Mar  4 16:02:28 2016
#
# Purpose:    Dump the contents of an array of hashes of hashes,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $AHHR:        ref of the array of hashes of hashes
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub utl_AHHDump { 
  my $sub_name = "utl_AHHDump()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $AHHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $undef2print = "!SEQUIP:undef!";
  my $nel1 = scalar(@{$AHHR});
  for(my $i1 = 0; $i1 < $nel1; $i1++) { 
    if(! defined $AHHR->[$i1]) { 
      print $FH ("*A*HH[$i1]: $undef2print\n");
    }
    else { 
      printf $FH ("*A*HH[$i1]\n");
      foreach my $key2 (sort keys %{$AHHR->[$i1]}) { 
        if(! defined $AHHR->[$i1]{$key2}) { 
          print $FH ("A*H*H[$i1]{$key2}: $undef2print\n");
        }
        else { 
          printf $FH ("\tA*H*H[$i1]{$key2}\n");
          foreach my $key3 (sort keys %{$AHHR->[$i1]{$key2}}) { 
            printf $FH ("\tAH*H*[$i1]{$key2}{$key3}: %s\n", (defined $AHHR->[$i1]{$key2}{$key3}) ? $AHHR->[$i1]{$key2}{$key3} : $undef2print);
          }
        }
        printf $FH ("\n");
      }
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: utl_AHDump()
# Incept:     EPN, Thu Feb  8 11:01:29 2018
#
# Purpose:    Dump the contents of an array of hashes,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $AHR:         ref of the array of hashes
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub utl_AHDump { 
  my $sub_name = "utl_AHDump()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $AHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $undef2print = "!SEQUIP:undef!";
  my $nel1 = scalar(@{$AHR});
  for(my $i1 = 0; $i1 < $nel1; $i1++) { 
    if(! defined $AHR->[$i1]) { 
      print $FH ("*A*H[$i1]: $undef2print\n");
    }
    else { 
      printf $FH ("*A*H[$i1]\n");
      my $nel2 = scalar(keys %{$AHR->[$i1]}); 
      foreach my $key2 (sort keys %{$AHR->[$i1]}) { 
        printf $FH ("\tA*H*[$i1]{$key2}: %s\n", (defined $AHR->[$i1]{$key2}) ? $AHR->[$i1]{$key2} : $undef2print);
      }
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: utl_HADump()
# Incept:     EPN, Tue Apr 21 06:04:07 2020
#
# Purpose:    Dump the contents of a hash of arrays
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $HAR:         ref of the hash of arrays
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub utl_HADump { 
  my $sub_name = "utl_HADump()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $HAR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $undef2print = "!SEQUIP:undef!";
  foreach my $key1 (sort keys %{$HAR}) { 
    if(defined $HAR->{$key1}) { 
      print $FH ("*H*A{$key1}: $undef2print\n");
    }
    else { 
      printf $FH ("*H*A{$key1}:\n");
      my $nel2 = scalar(@{$HAR->{$key1}});
      for(my $i2 = 0; $i2 < $nel2; $i2++) { 
        printf $FH ("\tH*A*{$key1}[$i2]: %s\n", (defined $HAR->{$key1}[$i2]) ? $HAR->{$key1}[$i2] : $undef2print);
      }
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: utl_HAHDump()
# Incept:     EPN, Tue Mar 19 12:30:24 2019
#
# Purpose:    Dump the contents of a hash of arrays of hashes,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $HAHR:        ref of the hash of array of hashes
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub utl_HAHDump { 
  my $sub_name = "utl_HAHDump()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $HAHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $undef2print = "!SEQUIP:undef!";
  foreach my $key1 (sort keys %{$HAHR}) { 
    if(! defined $HAHR->{$key1}) { 
      print $FH ("*H*AH{$key1}: $undef2print\n");
    }
    else { 
      printf $FH ("*H*AH{$key1}\n");
      my $nel2 = scalar(@{$HAHR->{$key1}});
      for (my $i2 = 0; $i2 < $nel2; $i2++) { 
        if(! defined $HAHR->{$key1}[$i2]) { 
          print $FH ("H*A*H{$key1}[$i2]: $undef2print\n");
        }
        else { 
          printf $FH ("\tH*A*H{$key1}[$i2]:\n", $i2);
          foreach my $key3 (sort keys %{$HAHR->{$key1}[$i2]}) { 
            printf $FH ("\t\tHA*H*{$key1}[$i2]{$key3}: %s\n", (defined $HAHR->{$key1}[$i2]{$key3}) ? $HAHR->{$key1}[$i2]{$key3} : $undef2print);
          }
        }
        printf $FH ("\n");
      }
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: utl_HHHDump()
# Incept:     EPN, Wed Mar 20 14:44:54 2019
#
# Purpose:    Dump the contents of a hash of hashes of hashes,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of hash of hashes of hashes
#             $HHHR:        ref of the hash of hashes of hashes
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub utl_HHHDump { 
  my $sub_name = "utl_HHHDump()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $HHHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $undef2print = "!SEQUIP:undef!";
  foreach my $key1 (sort keys %{$HHHR}) { 
    if(! defined $HHHR->{$key1}) { 
      print $FH ("*H*HH{$key1}: $undef2print\n");
    }
    else { 
      printf $FH ("*H*HH{$key1}\n");
      foreach my $key2 (sort keys %{$HHHR->{$key1}}) { 
        if(! defined $HHHR->{$key1}{$key2}) { 
          printf $FH ("\tH*H*H{$key1}{$key2}: $undef2print\n");
        }
        else{ 
          printf $FH ("\tH*H*H{$key1}{$key2}\n");
          foreach my $key3 (sort keys %{$HHHR->{$key1}{$key2}}) { 
            printf $FH ("\t\tHH*H*{$key1}{$key2}{$key3}: %s\n", (defined $HHHR->{$key1}{$key2}{$key3}) ? $HHHR->{$key1}{$key2}{$key3} : $undef2print);
          }
        }
        printf $FH ("\n");
      }
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: utl_HHADump()
# Incept:     EPN, Tue Apr 21 05:59:05 2020
#
# Purpose:    Dump the contents of a hash of hashes of arrays,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of hash of hashes of hashes
#             $HHAR:        ref of the hash of hashes of arrays
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub utl_HHADump { 
  my $sub_name = "utl_HHADump()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $HHAR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $undef2print = "!SEQUIP:undef!";
  foreach my $key1 (sort keys %{$HHAR}) { 
    if(! defined $HHAR->{$key1}) { 
      print $FH ("*H*HA{$key1}: $undef2print\n");
    }
    else { 
      printf $FH ("*H*HA{$key1}\n");
      foreach my $key2 (sort keys %{$HHAR->{$key1}}) { 
        if(! defined $HHAR->{$key1}{$key2}) { 
          print $FH ("H*H*A{$key1}{$key2}: $undef2print\n");
        }        
        else { 
          printf $FH ("\tH*H*A{$key1}{$key2}\n");
          my $nel3 = scalar(@{$HHAR->{$key1}{$key2}});
          for(my $i3 = 0; $i3 < $nel3; $i3++) { 
            printf $FH ("\t\tHH*A*{$key1}{$key2}[$i3]: %s\n", (defined $HHAR->{$key1}{$key2}[$i3]) ? $HHAR->{$key1}{$key2}[$i3] : $undef2print);
          }
        }
        printf $FH ("\n");
      }
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: utl_HHAHDump()
# Incept:     EPN, Mon May  4 21:52:00 2020
#
# Purpose:    Dump the contents of a hash of hashes of arrays of hashes,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of hash of hashes of hashes
#             $HHAHR:       ref of the hash of hashes of arrays of hashes
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub utl_HHAHDump { 
  my $sub_name = "utl_HHAHDump()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $HHAHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $undef2print = "!SEQUIP:undef!";
  foreach my $key1 (sort keys %{$HHAHR}) { 
    if(! defined $HHAHR->{$key1}) { 
      print $FH ("*H*HAH{$key1}: $undef2print\n");
    }
    else { 
      printf $FH ("*H*HAH{$key1}\n");
      foreach my $key2 (sort keys %{$HHAHR->{$key1}}) { 
        if(! defined $HHAHR->{$key1}{$key2}) { 
          print $FH ("H*H*AH{$key1}{$key2}: $undef2print\n");
        }
        else { 
          printf $FH ("\tH*H*AH{$key1}{$key2}\n");
          my $nel3 = scalar(@{$HHAHR->{$key1}{$key2}});
          for(my $i3 = 0; $i3 < $nel3; $i3++) { 
            printf $FH ("\t\tHH*A*H{$key1}{$key2}[$i3]\n");
            if(! defined $HHAHR->{$key1}{$key2}[$i3]) { 
              print $FH ("HH*A*H{$key1}{$key2}[$i3]: $undef2print\n");
            }
            else { 
              foreach my $key4 (sort keys %{$HHAHR->{$key1}{$key2}[$i3]}) { 
                printf $FH ("\t\t\tHHA*H*{$key1}{$key2}[{$i3]{$key4}: %s\n", (defined $HHAHR->{$key1}{$key2}[$i3]{$key4}) ? $HHAHR->{$key1}{$key2}[$i3]{$key4} : $undef2print);
              }
            }
            printf $FH ("\n");
          }
        }
        printf $FH ("\n");
      }
    }
    printf $FH ("\n");
  }
  
  return;
}

#####################################################################
# Subroutine: utl_IsInteger()
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
sub utl_IsInteger { 
  my $narg_expected = 1;
  my $sub_name = "utl_IsInteger()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($x) = @_;

  return ($x =~ m/^\-?\d+$/) ? 1 : 0;
}

#####################################################################
# Subroutine: utl_IsReal()
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
sub utl_IsReal { 
  my $narg_expected = 1;
  my $sub_name = "utl_IsReal()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($x) = @_;

  if(! defined $x)                { return 0; }
  if($x eq "")                    { return 0; }
  if($x =~ m/^\-?\d*\.?\d*$/)     { return 1; } # matches: 1, 3.4, 53.43, .3, .333, 0., 155., 150, -1, -3.4, -53.43, -.3, -.333, -0., -155., -150
  if($x =~ m/^\-?\d+[Ee]\-?\d*$/) { return 1; } # matches: 3E5, -3E5, -3E-5, 3e5, -3e5, -3e-5
  return 0;
}

#################################################################
# Subroutine:  utl_ASum()
# Incept:      EPN, Wed Feb  7 13:53:07 2018
# 
# Purpose:     Sum the scalar values in an array
#
# Arguments: 
#   $AR: reference to the array
# 
# Returns:     Sum of all values in the array
#
################################################################# 
sub utl_ASum {
  my $nargs_expected = 1;
  my $sub_name = "utl_ASum()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR) = $_[0];

  my $sum = 0;
  foreach my $el (@{$AR}) { 
    $sum += $el; 
  }
  return $sum;
}

#################################################################
# Subroutine:  utl_HSumValues()
# Incept:      EPN, Wed Feb  7 13:58:25 2018
# 
# Purpose:     Sum the values for all keys in a hash
#
# Arguments: 
#   $HR: reference to the hash
# 
# Returns:     Sum of all values in the hash
#
################################################################# 
sub utl_HSumValues {
  my $nargs_expected = 1;
  my $sub_name = "utl_HSumValues()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HR) = $_[0];

  my $sum = 0;
  foreach my $key (keys (%{$HR})) { 
    $sum += $HR->{$key};
  }
  return $sum;
}

#################################################################
# Subroutine:  utl_HSumValuesSubset()
# Incept:      EPN, Thu Apr 23 06:52:57 2020
# 
# Purpose:     Sum the values for a subset of keys in a hash.
#              Does not check that all values in @{$AR} are 
#              actually keys in %{$HR}.
#
# Arguments: 
#   $HR: reference to the hash
#   $AR: reference to the array with the subset of keys to sum
# 
# Returns:     Sum of all values in the hash
#
################################################################# 
sub utl_HSumValuesSubset {
  my $nargs_expected = 2;
  my $sub_name = "utl_HSumValuesSubset()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HR, $AR) = @_;

  my $sum = 0;
  foreach my $key (@{$AR}) { 
    if(defined $HR->{$key}) { 
      $sum += $HR->{$key};
    }
  }
  return $sum;
}

#################################################################
# Subroutine: utl_StringMonoChar()
# Incept:     EPN, Thu Mar 10 21:02:35 2016
#
# Purpose:    Return a string of length $len of repeated instances
#             of the character $char.
#
# Arguments:
#   $len:   desired length of the string to return
#   $char:  desired character
#   $FH_HR: REF to hash of file handles, including "log" and "cmd"
#
# Returns:  A string of $char repeated $len times.
# 
# Dies:     if $len is not a positive integer
#
#################################################################
sub utl_StringMonoChar {
  my $sub_name = "utl_StringMonoChar";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($len, $char, $FH_HR) = @_;

  if(! utl_IsInteger($len)) { 
    ofile_FAIL("ERROR in $sub_name, passed in length ($len) is not a non-negative integer", 1, $FH_HR);
  }
  if($len < 0) { 
    ofile_FAIL("ERROR in $sub_name, passed in length ($len) is a negative integer", 1, $FH_HR);
  }
    
  my $ret_str = "";
  for(my $i = 0; $i < $len; $i++) { 
    $ret_str .= $char;
  }

  return $ret_str;
}

#################################################################
# Subroutine:  utl_FileCountLines()
# Incept:      EPN, Tue Mar  1 09:36:56 2016
#
# Purpose:     Count the number of lines in a file
#              by opening it and reading in all lines.
#
# Arguments: 
#   $filename:         file that we are checking on
#   $FH_HR:            ref to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If $filename does not exist or cannot be opened for reading.
#
################################################################# 
sub utl_FileCountLines { 
  my $nargs_expected = 2;
  my $sub_name = "utl_FileCountLines()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($filename, $FH_HR) = @_;

  my $nlines = 0;
  open(IN, $filename) || ofile_FileOpenFailure($filename, $sub_name, $!, "reading", $FH_HR);
  while(<IN>) { 
    $nlines++;
  }
  close(IN);

  return $nlines;
}

#################################################################
# Subroutine:  utl_FileLinesToArray()
# Incept:      EPN, Tue Nov 21 10:26:58 2017
#
# Purpose:     Store each non-blank line in a file as an element
#              in an array, after removing newline.
#
# Arguments: 
#   $filename:                   file that we are parsing
#   $remove_trailing_whitespace: '1' to remove trailing whitespace in each line, '0' not to
#   $AR:                         ref to array to add to
#   $FH_HR:                      ref to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If $filename does not exist or cannot be opened for reading.
#
################################################################# 
sub utl_FileLinesToArray { 
  my $nargs_expected = 4;
  my $sub_name = "utl_FileLinesToArray()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($filename, $remove_trailing_whitespace, $AR, $FH_HR) = @_;

  open(IN, $filename) || ofile_FileOpenFailure($filename, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <IN>) { 
    if($line =~ /\S/) { 
      chomp $line;
      if($remove_trailing_whitespace) { $line =~ s/\s*$//; }
      push(@{$AR}, $line);
    }
  }
  close(IN);

  return;
}

#################################################################
# Subroutine:  utl_FileLinesToHash()
# Incept:      EPN, Tue Jul  2 09:20:13 2019
#
# Purpose:     Store each non-blank line in a file as an key
#              in an array, after removing newline. Values for
#              all keys will be '1'.
#
# Arguments: 
#   $filename:                   file that we are parsing
#   $remove_trailing_whitespace: '1' to remove trailing whitespace in each line, '0' not to
#   $HR:                         ref to array to add to
#   $FH_HR:                      ref to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If $filename does not exist or cannot be opened for reading.
#
################################################################# 
sub utl_FileLinesToHash { 
  my $nargs_expected = 4;
  my $sub_name = "utl_FileLinesToHash()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($filename, $remove_trailing_whitespace, $HR, $FH_HR) = @_;

  open(IN, $filename) || ofile_FileOpenFailure($filename, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <IN>) { 
    if($line =~ /\S/) { 
      chomp $line;
      if($remove_trailing_whitespace) { $line =~ s/\s*$//; }
      $HR->{$line} = 1;
    }
  }
  close(IN);

  return;
}

#################################################################
# Subroutine:  utl_FileRemoveList()
# Incept:      EPN, Fri Oct 19 12:44:05 2018 [ribovore]
#
# Purpose:     Remove each file in an array of file
#              names. If there are more than 100 files, then
#              remove 100 at a time.
# 
# Arguments: 
#   $files2remove_AR:  REF to array with list of files to remove
#   $caller_sub_name:  name of calling subroutine (can be undef)
#   $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:            ref to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If one of the rm -rf commands fails.
#
################################################################# 
sub utl_FileRemoveList { 
  my $nargs_expected = 4;
  my $sub_name = "utl_FileRemoveList()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($files2remove_AR, $caller_sub_name, $opt_HHR, $FH_HR) = @_;

  my $i = 0; 
  my $nfiles = scalar(@{$files2remove_AR});

  while($i < $nfiles) { 
    my $file_list = "";
    my $up = $i+100;
    if($up > $nfiles) { $up = $nfiles; }
    for(my $j = $i; $j < $up; $j++) { 
      $file_list .= " " . $files2remove_AR->[$j];
    }
    my $rm_cmd = "rm $file_list"; 
    utl_RunCommand($rm_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
    $i = $up;
  }
  
  return;
}

#################################################################
# Subroutine:  utl_AToNewLineDelimString()
# Incept:      EPN, Fri Dec 14 09:21:41 2018
#
# Purpose:     Return a newline delimited string with all values of an array.
#
# Arguments: 
#   $AR:     ref to array
# 
# Returns:     string
# 
# Dies:        Never.
#
################################################################# 
sub utl_AToNewLineDelimString {
  my $nargs_expected = 1;
  my $sub_name = "utl_AToNewLineDelimString";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR) = @_;

  my $retstr = "";
  foreach my $el (@{$AR}) { 
    $retstr .= $el . "\n";
  }
  if($retstr eq "") { 
    $retstr = "\n";
  }
  return $retstr;
}

#################################################################
# Subroutine:  utl_AToString()
# Incept:      EPN, Tue Mar 12 13:07:03 2019
#
# Purpose:     Return a string with all values of an array delimited
#              by $delim_char;
#
# Arguments: 
#   $AR:         ref to array
#   $delim_char: character to delimit with
# 
# Returns:     string
# 
# Dies:        Never.
#
################################################################# 
sub utl_AToString { 
  my $nargs_expected = 2;
  my $sub_name = "utl_AToString";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR, $delim_char) = (@_);

  my $retstr = "";
  if(@{$AR}) { 
    $retstr .= $AR->[0];
    for(my $i = 1; $i < scalar(@{$AR}); $i++) { 
      $retstr .= $delim_char . $AR->[$i];
    }
  }
  return $retstr;
}

#################################################################
# Subroutine:  utl_HKeysToNewLineDelimString()
# Incept:      EPN, Fri Dec 14 09:25:24 2018
#
# Purpose:     Return a newline delimited string with all (sorted) keys in a hash.
#
# Arguments: 
#   $HR:     ref to hash
# 
# Returns:     string
# 
# Dies:        Never.
#
################################################################# 
sub utl_HKeysToNewLineDelimString {
  my $nargs_expected = 1;
  my $sub_name = "utl_HKeysToNewLineDelimString";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HR) = @_;

  my $retstr = "";
  foreach my $el (sort keys %{$HR}) { 
    $retstr .= $el . "\n";
  }
  if($retstr eq "") { 
    $retstr = "\n";
  }
  return $retstr;
}

#################################################################
# Subroutine:  utl_HValuesToNewLineDelimString()
# Incept:      EPN, Fri Dec 14 09:26:16 2018
#
# Purpose:     Return a newline delimited string with all values in a hash.
#
# Arguments: 
#   $HR:     ref to hash
# 
# Returns:     string
# 
# Dies:        Never.
#
################################################################# 
sub utl_HValuesToNewLineDelimString {
  my $nargs_expected = 1;
  my $sub_name = "utl_HValuesToNewLineDelimString";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HR) = @_;

  my $retstr = "";
  foreach my $el (sort keys %{$HR}) { 
    $retstr .= $HR->{$el} . "\n";
  }
  if($retstr eq "") { 
    $retstr = "\n";
  }
  return $retstr;
}


#################################################################
# Subroutine : utl_FileValidateExistsAndNonEmpty()
# Incept:      EPN, Thu May  4 09:30:32 2017
#
# Purpose:     Check if a file exists and is non-empty. 
#
# Arguments: 
#   $filename:         file that we are checking on
#   $filedesc:         description of file
#   $calling_sub_name: name of calling subroutine (can be undef)
#   $do_die:           '1' if we should die if it does not exist.  
#   $FH_HR:            ref to hash of file handles, can be undef
# 
# Returns:     Return '1' if it does and is non empty
#              Return '0' if it does not exist (and ! $do_die)
#              Return '-1' if it exists but is empty (and ! $do_die)
#              Return '-2' if it exists as a directory (and ! $do_die)
#
# Dies:        If file does not exist or is empty and $do_die is 1.
#              if $filename is undefined (regardless of value of $do_die).
# 
################################################################# 
sub utl_FileValidateExistsAndNonEmpty { 
  my $nargs_expected = 5;
  my $sub_name = "utl_FileValidateExistsAndNonEmpty()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($filename, $filedesc, $calling_sub_name, $do_die, $FH_HR) = @_;

  if(! defined $filename) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, %sfilename%s is undef", 
                         (defined $calling_sub_name ? "called by $calling_sub_name," : ""),
                         (defined $filedesc         ? " ($filedesc)" : "")),
                 1, $FH_HR); 
  }

  if(-d $filename) {
    if($do_die) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, %sfile $filename%s exists but is a directory.", 
                         (defined $calling_sub_name ? "called by $calling_sub_name," : ""),
                         (defined $filedesc         ? " ($filedesc)" : "")),
                 1, $FH_HR); 
    }
    return -2;
  }
  elsif(! -e $filename) { 
    if($do_die) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, %sfile $filename%s does not exist.", 
                         (defined $calling_sub_name ? "called by $calling_sub_name," : ""),
                         (defined $filedesc         ? " ($filedesc)" : "")),
                 1, $FH_HR); 
    }
    return 0;
  }
  elsif(! -s $filename) { 
    if($do_die) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, %sfile $filename%s exists but is empty.", 
                         (defined $calling_sub_name ? "called by $calling_sub_name," : ""),
                         (defined $filedesc         ? " ($filedesc)" : "")),
                 1, $FH_HR); 
    }
    return -1;
  }
  
  return 1;
}


#################################################################
# Subroutine:  utl_FileMd5()
# Incept:      EPN, Fri May 27 14:02:30 2016
#
# Purpose:     Use md5sum to get a checksum of a file, return
#              the checksum. Not efficient. Creates a temporary
#              file and then deletes it.
# 
# Arguments: 
#   $file:             REF to array of all files to concatenate
#   $caller_sub_name:  name of calling subroutine (can be undef)
#   $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:            ref to hash of file handles
# 
# Returns:     md5sum of the file.
# 
# Dies:        If the file doesn't exist or the command fails.
#
################################################################# 
sub utl_FileMd5 { 
  my $nargs_expected = 4;
  my $sub_name = "utl_FileMd5()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($file, $caller_sub_name, $opt_HHR, $FH_HR) = @_;

  if(! -s $file) { 
    ofile_FAIL(sprintf("ERROR in $sub_name%s, file to get md5 checksum of ($file) does no exist or is empty", 
                        (defined $caller_sub_name) ? " called by $caller_sub_name" : ""), 1, $FH_HR);
  }

  my $out_file = utl_RemoveDirPath($file . ".md5sum");
  utl_RunCommand("md5sum $file > $out_file", opt_Get("-v", $opt_HHR), 0, $FH_HR);

  open(MD5, $out_file) || ofile_FileOpenFailure($out_file, $sub_name, $!, "reading", $FH_HR);
  #194625f7c3e2a5129f9880c7e29f63de  wnv.lin2.matpept.in
  my $md5sum = <MD5>;
  chomp $md5sum;
  if($md5sum =~ /^(\S+)\s+(\S+)$/) { 
    $md5sum = $1;
  }
  else { 
    ofile_FAIL(sprintf("ERROR in $sub_name%s, unable to parse md5sum output: $md5sum", 
                        (defined $caller_sub_name) ? " called by $caller_sub_name" : ""), 1, $FH_HR);
  }
  close(MD5);

  utl_FileRemoveUsingSystemRm($out_file, $caller_sub_name, $opt_HHR, $FH_HR);

  return $md5sum;
}

#################################################################
# Subroutine:  utl_DirEnvVarValid()
# Incept:      EPN, Wed Oct 25 10:09:28 2017 [ribo.pm]
#
# Purpose:     Verify that the environment variable $envvar exists 
#              and that it is a valid directory. Return directory path.
#              
# Arguments: 
#   $envvar:  environment variable
#
# Returns:    directory path $ENV{'$envvar'}
#
################################################################# 
sub utl_DirEnvVarValid { 
  my $nargs_expected = 1;
  my $sub_name = "utl_DirEnvVarValid()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($envvar) = $_[0];

  if(! exists($ENV{"$envvar"})) { 
    die "ERROR, the environment variable $envvar is not set";
    # it's okay this isn't ofile_FAIL because this is called before ofile_info_HH is set-up
  }
  my $envdir = $ENV{"$envvar"};
  if(! (-d $envdir)) { 
    die "ERROR, the directory specified by your environment variable $envvar does not exist.\n"; 
    # it's okay this isn't ofile_FAIL because this is called before ofile_info_HH is set-up
  }    

  return $envdir;
}

#################################################################
# Subroutine: utl_ExistsHFromCommaSepString
# Incept:     EPN, Mon Mar 18 06:52:21 2019
#
# Synopsis: Given a hash reference and a comma separated string
#           fill the hash with keys for each token in the string,
#           with all values set as 1.
#
# Arguments:
#  $HR:      hash reference
#  $string:  comma separated string
#
# Returns:    void
#
# Dies:       never
#################################################################
sub utl_ExistsHFromCommaSepString {
  my $sub_name = "utl_ExistsHFromCommaSepString";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($HR, $string) = @_;

  %{$HR} = ();
  my @key_A = split(",", $string);
  foreach my $key (@key_A) { 
    $HR->{$key} = 1; 
  }

  return;
}

#################################################################
# Subroutine:  utl_ExecHValidate()
# Incept:      EPN, Sat Feb 13 06:27:51 2016
#
# Purpose:     Given a reference to a hash in which the 
#              values are paths to executables, validate
#              those files are executable.
#
# Arguments: 
#   $execs_HR: REF to hash, keys are short names to executable
#              e.g. "cmbuild", values are full paths to that
#              executable, e.g. "/usr/local/infernal/1.1.1/bin/cmbuild"
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
# 
# Returns:     void
#
# Dies:        if one or more executables does not exist#
#
################################################################# 
sub utl_ExecHValidate { 
  my $nargs_expected = 2;
  my $sub_name = "utl_ExecHValidate()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($execs_HR, $FH_HR) = @_;

  my $fail_str = undef;
  foreach my $key (sort keys %{$execs_HR}) { 
    if(! -e $execs_HR->{$key}) { 
      $fail_str .= "\t$execs_HR->{$key} does not exist.\n"; 
    }
    elsif(! -x $execs_HR->{$key}) { 
      $fail_str .= "\t$execs_HR->{$key} exists but is not an executable file.\n"; 
    }
  }
  
  if(defined $fail_str) { 
    ofile_FAIL("ERROR in $sub_name(),\n$fail_str", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine:  utl_HCompare()
# Incept:      EPN, Mon Jun 24 15:35:19 2019
#
# Purpose:     Compare two single dimension hashes for equality.
#              
# Arguments: 
#   $H1R:      ref to hash 1
#   $H2R:      ref to hash 2
#   $name1:    name for hash 1, for return string
#   $name2:    name for hash 2, for return string
# 
# Returns:     string explaining differences between the two hashes
#              empty string if hashes are equal (same set of 
#              keys and values)
# 
################################################################# 
sub utl_HCompare { 
  my $sub_name = "utl_HCompare()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($H1R, $H2R, $name1, $name2) = (@_);

  my @A1 = keys (%{$H1R});
  my @A2 = keys (%{$H1R});

  my $ret_str = "";
  foreach my $key1 (@A1) { 
    if(! exists $H2R->{$key1}) { 
      $ret_str .= "$key1 exists as a key in $name1 hash but not in $name2 hash\n";
    }
    else { 
      if($H1R->{$key1} ne $H2R->{$key1}) { 
        $ret_str .= "values differ for key $key1 ($name1 hash value: " . $H1R->{$key1} . "; $name2 hash value: " . $H2R->{$key1} . ")\n";
      }
    }
  }
  foreach my $key2 (@A2) { 
    if(! exists $H1R->{$key2}) { 
      $ret_str .= "$key2 exists as a key in $name2 hash but not in $name1 hash\n";
    }
    # already checked all keys that exist in both above
  } 

  return $ret_str; 
}

#################################################################
# Subroutine:  utl_AToFile()
# Incept:      EPN, Wed Mar 18 14:02:09 2020
#
# Purpose:     Print lines of an array to a file.
#              
# Arguments: 
#   $AR:          ref to array
#   $filename:    output file to write array to
#   $add_newline: '1' to add newline after each element, 0 not to
#   $FH_HR:       REF to hash of file handles, including "log" and "cmd"
# 
# Returns:     void
# 
################################################################# 
sub utl_AToFile { 
  my $sub_name = "utl_AToFile()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($AR, $filename, $add_newline, $FH_HR) = (@_);

  open(OUT, ">", $filename) || ofile_FileOpenFailure($filename, $sub_name, $!, "writing", $FH_HR);

  if($add_newline) { 
    foreach my $el (@{$AR}) { 
      print OUT $el . "\n";
    }
  }
  else { 
    foreach my $el (@{$AR}) { 
      print OUT $el;
    }
  }

  close(OUT);

  return;
}

####################################################################
# the next line is critical, a perl module must return a true value
return 1;
####################################################################
