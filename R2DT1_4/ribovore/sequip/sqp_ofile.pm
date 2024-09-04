#!/usr/bin/perl
#
# sqp_ofile.pm
# Eric Nawrocki
# EPN, Tue Jul  2 11:16:31 2019 [epn-ofile.pm started]
# EPN, Tue Jul  2 11:43:02 2019 [migrated from epn-ofile]
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

#####################################################################
# Data structures used in this module:
#
#
# - $ofile_info_HHR: Reference to an hash of hashes in which all 2nd
#                    dim hashes are the same size and have the same
#                    set of keys, and contain information on an output
#                    file.
#                    
#                    The 1st dim hash keys describe the type of
#                    information, e.g. "fullpath" is the full path of
#                    the file, and 2nd dim keys pertain to which file,
#                    e.g. "log" for the log file which contains all
#                    the output printed to stdout during the course of
#                    the execution of the script. A special 1d key is
#                    'order' which is used for keeping track of the
#                    order that the files are added in, mainly so we
#                    can output information on them in the same
#                    order. The HH values for "order" are
#                    1..$num_ofiles, where $num_ofiles is the number
#                    of total output files (scalar(keys
#                    %{$ofile_info_HHR{"order"}})). Another special
#                    1d key 
#
#                    The set of 1D keys: 
#                    "order":     integer, the order in which this element was added
#                    "fullpath":  full path to the file 
#                    "nodirpath": file name, "fullpath" minus directories
#                    "mainout":   '1' if this file should be listed in the main output,
#                                 '0' if it should only be listed in the .list file
#                    "desc":      short description of the file
#                    "FH":        open file handle for this file, or undef             
#                                 See ofile_ValidateOutputFileInfoHashOfHashes()
#                                 for a list and explanation of all of the keys.
#
#####################################################################
#
# List of subroutines:
# 
#   ofile_OpenAndAddFileToOutputInfo()
#   ofile_AddClosedFileToOutputInfo()
#   ofile_HelperAddFileToOutputInfo()
#   ofile_ValidateOutputFileInfoHashOfHashes()
#   ofile_OutputProgressPrior()
#   ofile_OutputProgressComplete()
#   ofile_OutputConclusionAndCloseFilesOk()
#   ofile_OutputConclusionAndCloseFilesFail()
#   ofile_OutputTiming()
#   ofile_OutputString()
#   ofile_OutputBanner()
#   ofile_OutputDividingLine()
#   ofile_RemoveDirPath()
#   ofile_GetDirPath()
#   ofile_SecondsSinceEpoch()
#   ofile_FormatTimeString()
#   ofile_MaxLengthScalarValueInHash()
#   ofile_MaxLengthScalarValueInArray()
#   ofile_FileOpenFailure()
#   ofile_FAIL()
#
#################################################################
# Subroutine: ofile_OpenAndAddFileToOutputInfo()
# Incept:     EPN, Thu May 24 15:59:43 2018
#             EPN, Fri Feb 26 11:11:09 2016 [dnaorg_scripts:dnaorg.pm:openAndAddFileToOutputInfo()]
# 
# Purpose:    Add information about an output file to the
#             %{$ofile_info_HHR} and open that output file. Eventually
#             we'll output information about this file with
#             ofile_OutputConclusionAndCloseFiles().
#
#             Most of the work is done by helperAddFileToOutputInfo().
#
# Arguments:
#   $ofile_info_HHR:        REF to the 2D hash of output file information, ADDED TO HERE 
#   $key2d:                 2D key for the file we're adding and opening, e.g. "log"
#   $fullpath:              full path to the file we're adding and opening
#   $mainout:               '1' to output description of this file to 'main' when script ends, '0' not to
#   $listout:               '1' to output description of this file to {"FH"}{"list"} list file, '0' not to
#   $desc:                  description of the file we're adding and opening
#
# Returns:    void
# 
# Dies:       If $ofile_info_HHR{*}{$key} already exists.
#
#################################################################
sub ofile_OpenAndAddFileToOutputInfo { 
  my $sub_name = "ofile_OpenAndAddFileToOutputInfo";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($ofile_info_HHR, $key2d, $fullpath, $mainout, $listout, $desc) = @_;

  # this helper function does everything but open the file handle
  ofile_HelperAddFileToOutputInfo($ofile_info_HHR, $key2d, $fullpath, $mainout, $listout, $desc);

  # and open the file handle
  # we can only pass $FH_HR to ofile_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  if(! open($ofile_info_HHR->{"FH"}{$key2d}, ">", $fullpath)) { 
    ofile_FAIL("ERROR in $sub_name, unable to open $fullpath for writing.", 1, $FH_HR); 
  }

  return;
}

#################################################################
# Subroutine: ofile_AddClosedFileToOutputInfo()
# Incept:     EPN, Thu May 24 15:52:11 2018
#             EPN, Tue Feb 16 14:22:36 2016 [dnaorg_scripts:dnaorg.pm:addClosedFileToOutputInfo()]
# 
# Purpose:    Add information about a created output file (not open) to
#             the %{$ofile_info_HHR data structure, for eventual
#             output in ofile_OutputConclusionAndCloseFiles().
#
#             Most of the work is done by ofile_HelperAddFileToOutputInfo().
#
# Arguments:
#   $ofile_info_HHR:        REF to the 2D hash of output file information, ADDED TO HERE 
#                           for 1D key $key
#   $key2d:                 2D key for the file we're adding and opening, e.g. "fasta"
#   $fullpath:              full path to the closed file we're adding
#   $mainout:               '1' to output description of this file to 'main' when script ends, '0' not to
#   $listout:               '1' to output description of this file to {"FH"}{"list"} list file, '0' not to
#   $desc:                  description of the closed file we're adding 
#
# Returns:    void
# 
# Dies:       If $ofile_desc_HR->{$key} already exists.
#
#################################################################
sub ofile_AddClosedFileToOutputInfo { 
  my $sub_name = "ofile_AddClosedFileToOutputInfo()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($ofile_info_HHR, $key2d, $fullpath, $mainout, $listout, $desc) = @_;

  # this helper function does everything but set the file handle ("FH") value
  ofile_HelperAddFileToOutputInfo($ofile_info_HHR, $key2d, $fullpath, $mainout, $listout, $desc);

  # set FH value to undef
  $ofile_info_HHR->{"FH"}{$key2d} = undef;

  return;
}

#################################################################
# Subroutine: ofile_HelperAddFileToOutputInfo()
# Incept:     EPN, Thu May 24 15:54:03 2018
#             EPN, Fri Feb 26 14:35:36 2016 [dnaorg_scripts:dnaorg.pm:helperAddClosedFileToOutputInfo()]
# 
# Purpose:    Add information about an output file to the %{$ofile_info_HHR}
#             data structure. Helper function that's called by both 
#             openAndAddFileToOutputInfo() and addClosedFileToOutputInfo().
#             Also, if $ofile_info_HHR->{"FH"}{"list"} is defined, 
#             output the description of this file to the list file.
#
# Arguments:
#   $ofile_info_HHR:        REF to the 2D hash of output file information, ADDED TO HERE 
#                           for 1D key $key
#   $key2d:                 2D key for the file we're adding and opening, e.g. "log"
#   $fullpath:              full path to the file we're adding and opening
#   $mainout:               '1' to output description of this file to 'main' when script ends, '0' not to
#   $listout:               '1' to output description of this file to {"FH"}{"list"} list file, '0' not to
#   $desc:                  description of the file we're adding and opening
#
# Returns:    void
# 
# Dies:       If $ofile_info_HHR{*}{$key} already exists.
#
#################################################################
sub ofile_HelperAddFileToOutputInfo { 
  my $sub_name = "ofile_HelperAddFileToOutputInfo";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($ofile_info_HHR, $key2d, $fullpath, $mainout, $listout, $desc) = @_;

  # we can only pass $FH_HR to ofile_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  # make sure $mainout value is 0 or 1
  if(($mainout ne "0") && ($mainout ne "1")) { 
    ofile_FAIL("ERROR in $sub_name, entered with invalid 'mainout' value of $mainout (should be 0 or 1)", 1, $FH_HR);
  }
  # make sure $listout value is 0 or 1
  if(($listout ne "0") && ($listout ne "1")) { 
    ofile_FAIL("ERROR in $sub_name, entered with invalid 'listout' value of $listout (should be 0 or 1)", 1, $FH_HR);
  }
  # if $mainout is '1', $listout must be '1' also
  if(($mainout eq "1") && ($listout ne "1")) { 
    ofile_FAIL("ERROR in $sub_name, entered with invalid 'mainout'/'listout' value combination mainout:$mainout listout:$listout, if 'mainout' is '1', 'listout' must be as well", 1, $FH_HR);
  }

  # make sure we don't already have any information for this 2nd dim key $key2d:
  foreach my $key1d (keys (%{$ofile_info_HHR})) { 
    if(exists $ofile_info_HHR->{$key1d}{$key2d}) { 
      ofile_FAIL("ERROR in $sub_name, trying to add file $fullpath with key $key2d, but that key already exists for first dim key $key1d", 1, $FH_HR);
    }
  }

  # set the values of the 2D hash
  my $nodirpath = ofile_RemoveDirPath($fullpath);
  my $nidx = (defined $ofile_info_HHR->{"order"}) ? (scalar(keys %{$ofile_info_HHR->{"order"}})) : 0;
  $ofile_info_HHR->{"order"}{$key2d}     = $nidx+1; # first 2d key added will be '1', 2nd will be '2', etc.
  $ofile_info_HHR->{"fullpath"}{$key2d}  = $fullpath;
  $ofile_info_HHR->{"nodirpath"}{$key2d} = $nodirpath;
  $ofile_info_HHR->{"desc"}{$key2d}      = $desc;
  $ofile_info_HHR->{"mainout"}{$key2d}   = $mainout;
  $ofile_info_HHR->{"listout"}{$key2d}   = $listout;

  # output the description of this file to the list file
  my $list_FH = ((defined $ofile_info_HHR) && (defined $ofile_info_HHR->{"FH"}) && (exists $ofile_info_HHR->{"FH"}{"list"})) ? 
      $ofile_info_HHR->{"FH"}{"list"} : undef;

  if((defined $list_FH) && ($ofile_info_HHR->{"listout"}{$key2d})) { 
    my $width_desc = length("# ") + ofile_MaxLengthScalarValueInHash($ofile_info_HHR->{"desc"}) + length(" saved in:");
    if($width_desc < 80) { $width_desc = 80; }
    ofile_OutputString($list_FH, 0, sprintf("# %-*s %s\n", $width_desc, $ofile_info_HHR->{"desc"}{$key2d} . " saved in:", $ofile_info_HHR->{"nodirpath"}{$key2d}));
  }

  # validate that we've correctly updated the output info 2D hash
  ofile_ValidateOutputFileInfoHashOfHashes($ofile_info_HHR);

  return;
}

#################################################################
# Subroutine: ofile_ValidateOutputFileInfoHashOfHashes()
# Incept:     EPN, Thu May 24 15:43:05 2018 
#             EPN, Mon Feb 29 09:21:51 2016 [dnaorg_scripts:dnaorg.pm:validateOutputFileInfoHashOfHashes()]
#
# Purpose:    Validate an 'output file' info hash of hashes.
#             A valid info hash of hashes has the same set of 2d
#             keys for each 1d key except for "FH". The set of 1d keys is 
#             "order": integer, the order in which this element was added
#             "fullpath":  full path to the file 
#             "nodirpath": file name, "fullpath" minus directories
#             "mainout":   '1' to output description of this file to 'main' when script ends, '0' not to
#             "listout":   '1' to output description of this file to {"FH"}{"list"} list file, '0' not to
#             "desc":      short description of the file
#             "FH":        open file handle for this file, or undef             
#
#             For "FH", the value for each $ofile_info_HH{"FH"}{$key2d} can
#             be either defined or not defined.
#
# Arguments:
#   $ofile_info_HHR:  REF to hash of hashes of output file information
# 
# Returns: Number of elements in each and every 2d hash (except possibly %{$ofile_info_HHR->{"FH"}})
#
# Dies:    - if one of the expected keys (listed above) does not exist in $ofile_info_HHR
#          - if two 2d hashes in $ofile_info_HHR (besides %{$ofile_info_HHR->{"FH"}}) are of different sizes
#          - if two 2d hashes in $ofile_info_HHR (besides %{$ofile_info_HHR->{"FH"}}) have different set of keys
#################################################################
sub ofile_ValidateOutputFileInfoHashOfHashes { 
  my $sub_name = "ofile_ValidateOutputFileInfoHashOfHashes()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($ofile_info_HHR) = (@_);
  
  # we can only pass $FH_HR to ofile_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my @same_keys1d_A = ("order", "fullpath", "nodirpath", "mainout", "listout", "desc"); # all of these 2nd dim hashes should have same set of keys
  my @all_keys1d_A   = (@same_keys1d_A, "FH");             # all 1d keys
  my $i;     # a counter
  my $key1d; # a 1st dim key
  my $key2d; # a 2nd dim key

  # make sure we don't have any extra 1d keys we don't expect
  foreach $key1d (keys %{$ofile_info_HHR}) { 
    my $found_it = 0;
    foreach my $tmp_key1d (@all_keys1d_A) { 
      if($key1d eq $tmp_key1d) { 
        $found_it = 1;
      }
    }
    if($found_it == 0) { 
      ofile_FAIL("ERROR in $sub_name, unexpected 1d key $key1d exists.", 1, $FH_HR);
    }     
  } 
  
  # make sure all 2nd dim keys for all 1st dim keys are the same as the 2nd dim keys for 1st dim key "order"
  if(! defined $ofile_info_HHR->{"order"}) { 
    ofile_FAIL("ERROR in $sub_name, expected 1d key order does not exist.", 1, $FH_HR);
  }
  foreach my $key1d (@same_keys1d_A) { 
    if($key1d ne "order") { # skip "order"
      if(! defined $ofile_info_HHR->{$key1d}) { 
        ofile_FAIL("ERROR in $sub_name, expected 1d key $key1d does not exist.", 1, $FH_HR);
      }
      # we make sure the set of 2d keys in $ofile_info_HHR->{"order"} and $ofile_info_HHR->{$key1d} are 
      # identical in 2 steps:
      # 1) make sure all 2d keys from $ofile_info_HHR->{"order"} are also in $ofile_info_HHR->{"order"}
      foreach $key2d (keys %{$ofile_info_HHR->{"order"}}) { 
        if(! defined $ofile_info_HHR->{$key1d}{$key2d}) { 
          ofile_FAIL("ERROR in $sub_name, 2nd dim key $key2d exists for ofile_info_HHR->{order} but not for ofile_info_HHR->{$key1d}", 1, $FH_HR); 
        }
      }
      # 2) make sure all the 2d keys in $ofile_info_HHR->{$key1d} are also in $ofile_info_HHR->{"order"}
      foreach $key2d (keys %{$ofile_info_HHR->{$key1d}}) { 
        if(! defined $ofile_info_HHR->{"order"}{$key2d}) { 
          ofile_FAIL("ERROR in $sub_name, 2nd dim key $key2d exists for ofile_info_HHR->{order} but not for ofile_info_HHR->{$key1d}", 1, $FH_HR); 
        }
      }
    }
  }

  # make sure that $ofile_info_HHR->{"order"} has all values 1..$nkey2d
  my $nkey2d = scalar(keys %{$ofile_info_HHR->{"order"}});
  my @check_A = (); 
  for ($i = 1; $i <= $nkey2d; $i++) { 
    $check_A[$i] = 0; # changed to 1 when we observe it below
  }
  foreach $key2d (keys %{$ofile_info_HHR->{"order"}}) { 
    $check_A[$ofile_info_HHR->{"order"}{$key2d}] = 1;
  }
  for ($i = 1; $i <= $nkey2d; $i++) { 
    if($check_A[$i] != 1) { 
      ofile_FAIL("ERROR in $sub_name, invalid values for ofile_info_HH{order}, $nkey2d 2nd dim keys, but value $i does not exist", 1, $FH_HR);
    }
  }

  return $nkey2d;
}

#################################################################
# Subroutine : ofile_OutputProgressPrior()
# Incept:      EPN, Fri May 25 09:38:50 2018 
#              EPN, Fri Feb 12 17:22:24 2016 [dnaorg_scripts:dnaorg.pm:outputProgressPrior()]
#
# Purpose:      Output to $FH1 (and possibly $FH2) a message indicating
#               that we're about to do 'something' as explained in
#               $outstr.  
#
#               Caller should call *this* function, then do
#               the 'something', then call outputProgressComplete().
#
#               We return the number of seconds since the epoch, which
#               should be passed into the downstream
#               outputProgressComplete() call if caller wants to
#               output running time.
#
# Arguments: 
#   $outstr:     string to print to $FH
#   $progress_w: width of progress messages
#   $FH1:        file handle to print to
#   $FH2:        another file handle to print to, can be undef
# 
# Returns:     Number of seconds and microseconds since the epoch.
#
################################################################# 
sub ofile_OutputProgressPrior { 
  my $nargs_expected = 4;
  my $sub_name = "ofile_OutputProgressPrior()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($outstr, $progress_w, $FH1, $FH2) = @_;

  if(defined $FH1) { printf $FH1 ("# %-*s ... ", $progress_w, $outstr); }
  if(defined $FH2) { printf $FH2 ("# %-*s ... ", $progress_w, $outstr); }

  return ofile_SecondsSinceEpoch();
}

#################################################################
# Subroutine : ofile_OutputProgressComplete()
# Incept:      EPN, Fri May 25 09:39:40 2018
#              EPN, Fri Feb 12 17:28:19 2016 [dnaorg_scripts:dnaorg.pm:outputProgressComplete()]
#
# Purpose:     Output to $FH1 (and possibly $FH2) a 
#              message indicating that we've completed 
#              'something'.
#
#              Caller should call *this* function,
#              after both a call to outputProgressPrior()
#              and doing the 'something'.
#
#              If $start_secs is defined, we determine the number
#              of seconds the step took, output it, and 
#              return it.
#
# Arguments: 
#   $start_secs:    number of seconds either the step took
#                   (if $secs_is_total) or since the epoch
#                   (if !$secs_is_total)
#   $extra_desc:    extra description text to put after timing
#   $FH1:           file handle to print to
#   $FH2:           another file handle to print to, can be undef
# 
# Returns:     Number of seconds the step took (if $secs is defined,
#              else 0)
#
################################################################# 
sub ofile_OutputProgressComplete { 
  my $nargs_expected = 4;
  my $sub_name = "ofile_OutputProgressComplete()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($start_secs, $extra_desc, $FH1, $FH2) = @_;

  my $total_secs = undef;
  if(defined $start_secs) { 
    $total_secs = ofile_SecondsSinceEpoch() - $start_secs;
  }

  if(defined $FH1) { printf $FH1 ("done."); }
  if(defined $FH2) { printf $FH2 ("done."); }

  if(defined $total_secs || defined $extra_desc) { 
    if(defined $FH1) { printf $FH1 (" ["); }
    if(defined $FH2) { printf $FH2 (" ["); }
  }
  if(defined $total_secs) { 
    if(defined $FH1) { printf $FH1 (sprintf("%7.1f seconds%s", $total_secs, (defined $extra_desc) ? ", " : "")); }
    if(defined $FH2) { printf $FH2 (sprintf("%7.1f seconds%s", $total_secs, (defined $extra_desc) ? ", " : "")); }
  }
  if(defined $extra_desc) { 
    if(defined $FH1) { printf $FH1 $extra_desc };
    if(defined $FH2) { printf $FH2 $extra_desc };
  }
  if(defined $total_secs || defined $extra_desc) { 
    if(defined $FH1) { printf $FH1 ("]"); }
    if(defined $FH2) { printf $FH2 ("]"); }
  }

  if(defined $FH1) { printf $FH1 ("\n"); }
  if(defined $FH2) { printf $FH2 ("\n"); }
  
  return (defined $total_secs) ? $total_secs : 0.;
}

#######################################################################
# Subroutine: ofile_OutputConclusionAndCloseFilesOk()
# Incept:     EPN, Fri May 25 09:40:26 2018
#             EPN, Thu Nov  5 18:25:31 2009 [ssu-align] 
# 
# Purpose:    Output a list of the main output files created 
#             and the final few lines of output and optionally the 
#             run time timing to the summary file. Print date and
#             system information to the log file. 
# 
#             Close all open file handles.
#
#             Final line will be '[ok]', this is only difference with
#             ofile_OutputConclusionsAndCloseFilesFail()
# Arguments: 
#  $total_secs:            total number of seconds, "" to not print timing
#  $odir:                  output directory, if "", files were put in cwd
#  $ofile_info_HHR:        REF to the 2D hash of output file information
#
# Returns:   Nothing.
# 
# Dies:      Never.
#
####################################################################
sub ofile_OutputConclusionAndCloseFilesOk { 
  my $nargs_expected = 3;
  my $sub_name = "ofile_OutputConclusionAndCloseFilesOk";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($total_secs, $odir, $ofile_info_HHR) = @_;

  return _ofile_helper_output_conclusion_and_close_files($total_secs, $odir, "[ok]", $ofile_info_HHR);
}

#######################################################################
# Subroutine: ofile_OutputConclusionAndCloseFilesFail()
# Incept:     EPN, Fri May 25 09:40:26 2018
#             EPN, Thu Nov  5 18:25:31 2009 [ssu-align] 
# 
# Purpose:    Output a list of the main output files created 
#             and the final few lines of output and optionally the 
#             run time timing to the summary file. Print date and
#             system information to the log file. 
# 
#             Close all open file handles.
#
#             Final line will be '[FAIL]', this is only difference with
#             ofile_OutputConclusionsAndCloseFilesOk()
# Arguments: 
#  $total_secs:            total number of seconds, "" to not print timing
#  $odir:                  output directory, if "", files were put in cwd
#  $ofile_info_HHR:        REF to the 2D hash of output file information
#
# Returns:   Nothing.
# 
# Dies:      Never.
#
####################################################################
sub ofile_OutputConclusionAndCloseFilesFail { 
  my $nargs_expected = 3;
  my $sub_name = "ofile_OutputConclusionsndCloseFilesFail";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($total_secs, $odir, $ofile_info_HHR) = @_;

  return _ofile_helper_output_conclusion_and_close_files($total_secs, $odir, "[FAIL]", $ofile_info_HHR);
}

#######################################################################
# Subroutine: _ofile_helper_output_conclusion_and_close_files()
# Incept:     EPN, Fri May 25 09:40:26 2018
#             EPN, Thu Nov  5 18:25:31 2009 [ssu-align] 
# 
# Purpose:    Output a list of the main output files created 
#             and the final few lines of output and optionally the 
#             run time timing to the summary file. Print date and
#             system information to the log file. 
#
#             Close all open file handles.
#
# Arguments: 
#  $total_secs:            total number of seconds, "" to not print timing
#  $odir:                  output directory, if "", files were put in cwd
#  $final_msg:             string to put on final line: "[ok]" or "[FAIL]"
#  $ofile_info_HHR:        REF to the 2D hash of output file information
#
# Returns:   Nothing.
# 
# Dies:      Never.
#
####################################################################
sub _ofile_helper_output_conclusion_and_close_files { 
  my $nargs_expected = 4;
  my $sub_name = "_ofile_helper_output_conclusion_and_close_files";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($total_secs, $odir, $final_msg, $ofile_info_HHR) = @_;

  ofile_ValidateOutputFileInfoHashOfHashes($ofile_info_HHR);

  my $log_FH  = $ofile_info_HHR->{"FH"}{"log"};
  my $cmd_FH  = $ofile_info_HHR->{"FH"}{"cmd"};
  my $list_FH = $ofile_info_HHR->{"FH"}{"list"};

  my $key2d; # a key in the 2nd dimension of $ofile_info_HHR

  # output the list of files that we created for which the 'mainout' variable 
  # is 1 to the $log file and stdout, we already printed the descriptions
  # to the list file in helperAddFileToOutputInfo().
  ofile_OutputString($log_FH, 1, sprintf("#\n"));
  # create a temporary array with description of files with 'mainout' set to 1 (we'll only print these)
  # so we get pretty formatting
  my @tmp_A = ();
  foreach $key2d (keys (%{$ofile_info_HHR->{"desc"}})) { 
    if($ofile_info_HHR->{"mainout"}{$key2d}) { 
      push(@tmp_A, $ofile_info_HHR->{"desc"}{$key2d});
    }
  }
  my $width_desc = length("# ") + ofile_MaxLengthScalarValueInArray(\@tmp_A) + length(" saved in:");
  my $cur_idx = 1;
  my $num_ofile = ofile_ValidateOutputFileInfoHashOfHashes($ofile_info_HHR); # this function validates we have exactly 1 of each "order" value of 1..$num_ofile
  for(my $i = 1; $i <= $num_ofile; $i++) { 
    foreach $key2d (keys (%{$ofile_info_HHR->{"order"}})) { 
      if(($ofile_info_HHR->{"order"}{$key2d} == $i) && 
         ($ofile_info_HHR->{"mainout"}{$key2d})) { # only print out files for which the special "mainout" value is '1'
        ofile_OutputString($log_FH, 1, sprintf("# %-*s %s\n", $width_desc, $ofile_info_HHR->{"desc"}{$key2d} . " saved in:", $ofile_info_HHR->{"nodirpath"}{$key2d}));
      }
    }
  }
  ofile_OutputString($log_FH, 1, sprintf("#\n"));
  ofile_OutputString($log_FH, 1, sprintf("# All output files created in %s\n", ($odir eq "") ? "the current working directory" : "directory \.\/$odir\/"));
  ofile_OutputString($log_FH, 1, sprintf("#\n"));
  if($total_secs ne "") { # don't print this if rvr-align is caller
    ofile_OutputTiming("# Elapsed time: ", $total_secs, 1, $log_FH); 
    ofile_OutputString($log_FH, 1, "#                hh:mm:ss\n");
    ofile_OutputString($log_FH, 1, "# \n");
    ofile_OutputString($log_FH, 1, $final_msg . "\n");
  }

  if(defined $cmd_FH) { 
    ofile_OutputString($cmd_FH, 0, "# " . `date`);      # prints date,        e.g.: 'Mon Feb 22 16:37:09 EST 2016'
    ofile_OutputString($cmd_FH, 0, "# " . `uname -a`);  # prints system info, e.g.: 'Linux cbbdev13 2.6.32-573.7.1.el6.x86_64 #1 SMP Tue Sep 22 22:00:00 UTC 2015 x86_64 x86_64 x86_64 GNU/Linux'
    if($total_secs ne "") { 
      ofile_OutputString($cmd_FH, 0, $final_msg . "\n");
    }
  }
  
  # close any open file handles
  foreach $key2d (keys (%{$ofile_info_HHR->{"FH"}})) { 
    if(defined $ofile_info_HHR->{"FH"}{$key2d}) { 
      close $ofile_info_HHR->{"FH"}{$key2d};
    }
  }

  return;
}

#####################################################################
# Subroutine: ofile_OutputTiming()
# Incept:     EPN, Fri May 25 09:42:44 2018
#             EPN, Tue Jun 16 08:52:08 2009 [ssu-align]
#
# Purpose:    Output elapsed time in hhhh:mm:ss format.
# 
# Arguments:
#   $prefix:               string to print before the hhhh:mm:ss time info.
#   $inseconds:            number of seconds
#   $print_to_stdout:      '1' to print to stdout, '0' not to
#   $FH:                   file handle to print to
#
# Returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub ofile_OutputTiming { 
  my $nargs_expected = 4;
  my $sub_name = "ofile_OutputTiming()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($prefix, $inseconds, $print_to_stdout, $FH) = @_;

  my $time_str = ofile_FormatTimeString($inseconds);
  ofile_OutputString($FH, $print_to_stdout, $prefix . " " . $time_str . "\n"); # this will always end with a newline
  
  return;
}

###########################################################
# Subroutine: ofile_OutputString()
# Incept: EPN, Fri May 25 09:35:58 2018
#         EPN, Wed Oct 29 20:42:16 2014 [rnavore]
#
# Purpose: Given a string and an open file handle <$FH>, 
#          output the string to the file handle 
#          and potentially to stdout as well. 
#          If <$FH> is not defined then do not 
#          print to a file. 
#
# Arguments:
#   $FH:              file handle to output to, can be undef
#   $print_to_stdout: if '1' also output string to stdout
#   $string:          the string to output
#
# Returns: Nothing. 
#
# Dies:    Never.
#
###########################################################
sub ofile_OutputString {
  my $nargs_expected = 3;
  my $sub_name = "ofile_OutputString()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($FH, $print_to_stdout, $string) = @_;

  if(defined $FH)      { print $FH $string; }
  if($print_to_stdout) { print     $string; }

  return;
}

#####################################################################
# Subroutine: ofile_OutputBanner()
# Incept:     EPN, Thu Oct 30 09:43:56 2014 (rnavore)
# 
# Purpose:    Output the banner.
#
# Arguments: 
#    $FH:                file handle to print to
#    $pkg:               name of package
#    $version:           version
#    $releasedate:       month/year of version (e.g. "Feb 2016")
#    $synopsis:          string reporting the date
#    $date:              date information to print
#    $extra_HR:          ref to hash of additional information to output
#                        key is left hand column, value is right hand column
#                        can be undefined
#
# Returns:    Nothing
# 
# Dies: never
####################################################################
sub ofile_OutputBanner {
  my $nargs_expected = 7;
  my $sub_name = "ofile_OutputBanner()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($FH, $pkg, $version, $releasedate, $synopsis, $date, $extra_HR) = @_;

  print $FH ("\# $synopsis\n");
  print $FH ("\# $pkg $version ($releasedate)\n");
#  print $FH ("\# Copyright (C) 2014 HHMI Janelia Research Campus\n");
#  print $FH ("\# Freely distributed under the GNU General Public License (GPLv3)\n");
  print $FH ("\# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  my $width = 9;
  if(defined $extra_HR) { 
    foreach my $key (sort keys (%{$extra_HR})) { 
      if((length($key)+3) > $width) { 
        $width = length($key)+3;
      }
    }
  }
  if(defined $date) {
    printf $FH ("%-*s  $date\n", $width, "# date:"); 
  }
  if(defined $extra_HR) { 
    foreach my $key (sort keys (%{$extra_HR})) { 
      printf $FH ("%-*s  $extra_HR->{$key}\n", $width, "# $key:"); 
    }
  }
  printf $FH ("#\n");

  return;
}

#################################################################
# Subroutine:  ofile_OutputDividingLine()
# Incept:      EPN, Tue Apr 12 14:40:13 2016
#
# Purpose:     Print a line of dashes followed by single spaces
#              with $ndash dashes to file handle $FH.
#              if $ndash is undefined, set it to 66.
#
# Arguments: 
#   $ndashes:  number of dashes in output dividing line
#   $FH:       file handle to print to
# 
# Returns:     Nothing.
# 
# Dies:        Never.
#
################################################################# 
sub ofile_OutputDividingLine { 
  my $nargs_expected = 2;
  my $sub_name = "ofile_OutputDividingLine()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ndashes, $FH) = @_;

  if(! defined $ndashes) { 
    $ndashes = 66;
  }

  my $div_line = "#";
  for(my $i = 0; $i < $ndashes; $i++) { 
    $div_line .= " -";
  }
  $div_line .= "\n";

  print $FH $div_line;
  
  return;
}

#################################################################
# Subroutine : ofile_RemoveDirPath()
# Incept:      EPN, Fri May 25 10:14:16 2018
#              EPN, Mon Nov  9 14:30:59 2009 [ssu-align]
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
sub ofile_RemoveDirPath {
  my $sub_name = "ofile_RemoveDirPath()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my $fullpath = $_[0];

  $fullpath =~ s/^.+\///;

  return $fullpath;
}

#################################################################
# Subroutine : ofile_GetDirPath()
# Incept:      EPN, Thu May  4 09:39:06 2017
#              EPN, Mon Mar 15 10:17:11 2010 [ssu.pm:ReturnDirPath()]
#
# Purpose:     Given a file name return the directory path, with the final '/'
#              For example: "foodir/foodir2/foo.stk" becomes "foodir/foodir2/".
#
# Arguments: 
#   $orig_file: name of original file
# 
# Returns:     The string $orig_file with actual file name removed 
#              or "./" if $orig_file is "".
#
################################################################# 
sub ofile_GetDirPath {
  my $narg_expected = 1;
  my $sub_name = "ofile_GetDirPath()";
  if(scalar(@_) != $narg_expected) { printf STDERR ("ERROR, in $sub_name, entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
  my $orig_file = $_[0];
  
  $orig_file =~ s/[^\/]+$//; # remove final string of non '/' characters
  
  if($orig_file eq "") { return "./";       }
  else                 { return $orig_file; }
}

#################################################################
# Subroutine : ofile_SecondsSinceEpoch()
# Incept:      EPN, Sat Feb 13 06:17:03 2016
#              EPN, Sat Feb 13 06:17:03 2016 [dnaorg_scripts:dnaorg.pm:ofile_SecondsSinceEpoch()]
#
# Purpose:     Return the seconds and microseconds since the 
#              Unix epoch (Jan 1, 1970) using 
#              Time::HiRes::gettimeofday().
#
# Arguments:   NONE
# 
# Returns:     Number of seconds and microseconds
#              since the epoch.
#
################################################################# 
sub ofile_SecondsSinceEpoch { 
  my $nargs_expected = 0;
  my $sub_name = "ofile_SecondsSinceEpoch()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($seconds, $microseconds) = gettimeofday();
  return ($seconds + ($microseconds / 1000000.));
}

#####################################################################
# Subroutine: ofile_FormatTimeString()
# Incept:     EPN, Fri May 25 10:16:35 2018
#             EPN, Fri Oct 24 13:18:23 2014 [rnavore]
#
# Purpose:    Get a timing in hhhh:mm:ss format.
# 
# Arguments:
# $inseconds: number of seconds
#
# Returns:    string that describes time in hhhh:mm:ss format
# 
# Dies:       Never.
#
####################################################################
sub ofile_FormatTimeString { 
  my $nargs_expected = 1;
  my $sub_name = "ofile_FormatTimeString()";
   if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($inseconds) = @_;
  my ($i, $hours, $minutes, $seconds, $thours, $tminutes, $tseconds, $ndig_hours);

  $hours = int($inseconds / 3600);
  $inseconds -= ($hours * 3600);
  $minutes = int($inseconds / 60);
  $inseconds -= ($minutes * 60);
  $seconds = $inseconds;
  $thours   = sprintf("%02d", $hours);
  $tminutes = sprintf("%02d", $minutes);
  $ndig_hours = length($hours);
  if($ndig_hours < 2) { $ndig_hours = 2; }

  $tseconds = sprintf("%05.2f", $seconds);
  my $ret_str = sprintf("%*s:%2s:%5s", $ndig_hours, $thours, $tminutes, $tseconds);
  # %*s covers two of the arguments: $ndig_hours specifies width of string, $thours is the string
  
  return $ret_str;
}

#################################################################
# Subroutine : ofile_MaxLengthScalarValueInHash()
# Incept:      EPN, Fri May 25 10:30:21 2018
#              EPN, Mon Nov  3 09:09:59 2014 [rnavore]
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
sub ofile_MaxLengthScalarValueInHash { 
  my $nargs_expected = 1;
  my $sub_name = "ofile_MaxLengthScalarValueInHash()";
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
# Subroutine : ofile_MaxLengthScalarValueInArray()
# Incept:      EPN, Fri May 25 10:29:59 2018
#              EPN, Thu Mar 17 12:38:53 2016 [dnaorg_scripts:dnaorg.pm:maxLengthScalarValueInArray()]
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
sub ofile_MaxLengthScalarValueInArray { 
  my $nargs_expected = 1;
  my $sub_name = "ofile_MaxLengthScalarValueInArray()";
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
# Subroutine : ofile_FileOpenFailure()
# Incept:      EPN, Fri May 25 11:49:18 2018
#              EPN, Wed Nov 11 05:39:56 2009 (rnavore)
#
# Purpose:     Called if an open() call fails on a file.
#              Print an informative error message
#              to $FH_HR->{"cmd"} and $FH_HR->{"log"}
#              and to STDERR, then exit with <$status>.
#
# Arguments: 
#   $filename:   file that we couldn't open
#   $c_sub_name: name of calling subroutine name
#   $status:     error status
#   $action:     "reading", "writing", "appending"
#   $FH_HR:      ref to hash of open file handles to close
# 
# Returns:     Nothing, this function will exit the program.
#
################################################################# 
sub ofile_FileOpenFailure { 
  my $nargs_expected = 5;
  my $sub_name = "ofile_FileOpenFailure()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($filename, $c_sub_name, $status, $action, $FH_HR) = @_;

  if(($action eq "reading") && (! (-e $filename))) { 
    ofile_FAIL(sprintf("ERROR, could not open %s%s for reading. It does not exist.", $filename, (defined $c_sub_name) ? " in subroutine $c_sub_name" : ""), $status, $FH_HR);
  }
  else { 
    ofile_FAIL(sprintf("ERROR, could not open %s%s for %s", $filename, (defined $c_sub_name) ? " in subroutine $c_sub_name" : "", $action), $status, $FH_HR);
  }

  return; # never reached
}

#################################################################
# Subroutine:  ofile_FAIL()
# Incept:      EPN, Wed Nov 11 06:22:59 2009 (rnavore) 
#
# Purpose:     Print an error message to STDERR and sum and 
#              log files in $FH_HR (ref to hash of file handles)
#              then close those file handles and exit.
#
# Arguments: 
#   $errmsg:  the error message to write
#   $status:  error status 
#   $FH_HR:   ref to hash of file handles, including "log" and "cmd"
# 
# Returns:     Nothing, this function will exit the program.
#
################################################################# 
sub ofile_FAIL { 
  my $nargs_expected = 3;
  my $sub_name = "ofile_FAIL()";
  if(scalar(@_) != $nargs_expected) { 
    if(scalar(@_) > 0) { 
      printf STDERR ("ERROR, ofile_FAIL() entered with %d != %d input arguments.\n(errmsg: $_[0])\n\n", scalar(@_), $nargs_expected); 
    }
    else { 
      printf STDERR ("ERROR, ofile_FAIL() entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); 
    }
    exit(1); 
  }
  my ($errmsg, $status, $FH_HR) = @_;

  if($errmsg !~ m/\n$/) { $errmsg .= "\n\n"; }
  else                  { $errmsg .= "\n"; }
  if($errmsg !~ m/^\n/) { $errmsg = "\n" . $errmsg; }
  
  if(defined $FH_HR) { 
    my $cmd_FH = $FH_HR->{"cmd"};
    my $log_FH = $FH_HR->{"log"};
    if(defined $cmd_FH) { 
      print $cmd_FH $errmsg;
      print $cmd_FH ("[fail]\n");
    }
    if(defined $log_FH) {
      print $log_FH $errmsg;
      print $log_FH ("[fail]\n");
    }
    # close each file handle
    foreach my $key (keys %{$FH_HR}) { 
      if(defined $FH_HR->{$key}) { 
        close $FH_HR->{$key};
      }
    }
  }
  
  print STDERR $errmsg; 
  exit 1;
}

#################################################################
# Subroutine:  ofile_TableHumanOutput()
# Incept:      EPN, Fri Apr  5 05:34:54 2019
#
# Purpose:     Output a human readable table with appropriate
#              width columns.
#
# Arguments: 
#   $data_AAR:   ref to 2D array of table data [0..$nrow-1][0..$ncol-1], 
#                must be defined, if $data_AAR->[$x] is undef, we will 
#                print a blank line (prefixed with $pfx_com)
#   $head_AAR:   ref to 2D array of header info [0..$nrow-1][0..$ncol-1], 
#                can be undef, will print no header if undefined
#   $cljust_AR:  ref to '1'/'0' array of indicating if a column is left justified
#                or not, can be undef, all columns will be left just if undefined
#                $cljust_A[0] will be set to '1' regardless of input
#   $bcom_AR:    ref to array of comment lines to add before table, 
#                can be undef, no lines if undefined
#   $acom_AR:    ref to array of comment lines to add after table, 
#                can be undef, no lines if undefined
#   $csep:       column separation string, can be undef
#                set to "  " if undef
#   $hsep:       header separation line character, "" for no line, 
#                can be undef, set to "-" if undef
#   $pfx_head:   header row prefix string, "" for none, 
#                can be undef set to "#" if undef
#   $pfx_com:    comment row prefix string, "" for none, 
#                can be undef, set to "# " if undef
#   $pfx_data:   data row prefix string, "" for none, can be undef
#                set to "" if undef
#   $empty_flag: '1' to output empty lines for empty arrays of data, '0'
#                to output empty lines as header separation lines
#   $out_FH1:    ref to output file handle 1, must be defined
#   $out_FH2:    ref to output file handle 1, can be undef
#   $FH_HR:      ref to hash of other file handles, including "log" and "cmd"
#  
# Returns:     Nothing.
#
# Dies:        
################################################################# 
sub ofile_TableHumanOutput { 
  my $nargs_expected = 14;
  my $sub_name = "ofile_TableHumanOutput";
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($data_AAR, $head_AAR, $cljust_AR, $bcom_AR, $acom_AR, $csep, $hsep, 
      $pfx_head, $pfx_com, $pfx_data, $empty_flag, $out_FH1, $out_FH2, $FH_HR) = (@_);

  # contract checks to make sure we have our required input
  # $data_AAR must be defined and have at least 1 row
  if(! defined $data_AAR) { 
    ofile_FAIL("ERROR in $sub_name, data 2D array is not defined", 1, $FH_HR);
  }
  # $out_FH1 must be defined
  if(! defined $out_FH1) { 
    ofile_FAIL("ERROR in $sub_name, output file handle 1 is not defined", 1, $FH_HR);
  }

  # make sure all per-column data that are not empty have same number of columns
  my $r; # row counter
  my $c; # column counter
  my $ncol_data = 0;
  my $nrow_data = scalar(@{$data_AAR});
  if($nrow_data > 0) { 
    $r = 0;
    while(($r < $nrow_data) && (scalar(@{$data_AAR->[$r]}) == 0)) { 
      $r++; 
    }
    if($r < $nrow_data) { 
      $ncol_data = scalar(@{$data_AAR->[$r]});
      my $first_r = $r;
      for($r = ($first_r+1); $r < $nrow_data; $r++) { 
        if((scalar(@{$data_AAR->[$r]}) > 0) && (scalar(@{$data_AAR->[$r]}) != $ncol_data)) { # allow empty arrays --> blank lines
          ofile_FAIL("ERROR in $sub_name, data row %d has $ncol_data columns, but row " . ($first_r+1) . " has " . scalar(@{$data_AAR->[$r]}) . " columns", 1, $FH_HR);
        }
      }
    }
  }

  my $nrow_head = (defined $head_AAR) ? scalar(@{$head_AAR}) : 0;
  my $ncol_head = 0;
  if($nrow_head > 0) { 
    $ncol_head = scalar(@{$head_AAR->[0]});
    if((defined $head_AAR) && (scalar(@{$head_AAR}) == 0)) { 
      ofile_FAIL("ERROR in $sub_name, header defined but is empty", 1, $FH_HR);
    }
    for($r = 1; $r < $nrow_head; $r++) { 
      if(scalar(@{$head_AAR->[$r]}) != $ncol_head) { 
        ofile_FAIL("ERROR in $sub_name, header row 1 has $ncol_head columns, but header row " . ($r+1) . " has " . scalar(@{$head_AAR->[$r]}) . " columns", 1, $FH_HR);
      }
    }
  }
  my $ncol;
  if(($ncol_head == 0) && ($ncol_data == 0)) { 
    ofile_FAIL("ERROR in $sub_name, data 2D array and header 2D arrays are both empty", 1, $FH_HR);
  }
  elsif(($ncol_head > 0) && ($ncol_data == 0)) { 
    $ncol = $ncol_head;
  }
  elsif(($ncol_head == 0) && ($ncol_data > 0)) { 
    $ncol = $ncol_data;
  }
  elsif($ncol_head == $ncol_data) { 
    $ncol = $ncol_data;
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, data num columns ($ncol_data) differs from header num columns ($ncol_head)", 1, $FH_HR);
  }

  # cljust_AR must be undef, or have $ncol values all of which must be '0' or '1'
  if(defined $cljust_AR) { 
    if(scalar(@{$cljust_AR}) != $ncol) { 
      ofile_FAIL("ERROR in $sub_name, data row 1 has $ncol columns, but column left justification array has data for " . scalar(@{$cljust_AR}) . " columns", 1, $FH_HR);
    }
    for($c = 0; $c < $ncol; $c++) { 
      if(($cljust_AR->[$c] != "0") && ($cljust_AR->[$c] != "1")) { 
        ofile_FAIL("ERROR in $sub_name, cljust_AR column " . ($c+1) . " value is " . $cljust_AR->[$c] . ", but it must be 0 or 1", 1, $FH_HR);
      }
    }
  }
  # if we get here, input has passed all contract checks

  # enforce defaults
  if(! defined $csep)      { $csep     = "  "; }
  if(! defined $hsep)      { $hsep     = "-"; }
  if(! defined $pfx_head)  { $pfx_head = "#";  }
  if(! defined $pfx_com)   { $pfx_com  = "# "; }
  if(! defined $pfx_data)  { $pfx_data = "";   }

  # create the left justification array, either by copying input or making default
  my @cljust_A = ();
  if(! defined $cljust_AR) { 
    for($c = 0; $c < $ncol; $c++) { $cljust_A[$c] = 1; } # default to left justification
  }
  else { 
    @cljust_A = @{$cljust_AR};
  }
  $cljust_A[0] = 1; # first column must be left justified
    
  # determine column widths
  my @w_A = ();
  for($c = 0; $c < $ncol; $c++) { $w_A[$c] = 0; } # initialize
  # consider data in each data row
  for($r = 0; $r < $nrow_data; $r++) { 
    if(scalar(@{$data_AAR->[$r]}) > 0) { 
      $w_A[0] = utl_Max($w_A[0], length($data_AAR->[$r][0]) + length($pfx_data)); 
      for($c = 1; $c < ($ncol-1); $c++) { # don't update w_A for final column
        $w_A[$c] = utl_Max($w_A[$c], length($data_AAR->[$r][$c])); 
      }
    }
  }
  # consider headers in each header row
  for($r = 0; $r < $nrow_head; $r++) { 
    $w_A[0] = utl_Max($w_A[0], (length($head_AAR->[$r][0]) + length($pfx_head))); 
    for($c = 1; $c < $ncol; $c++) { 
      $w_A[$c] = utl_Max($w_A[$c], length($head_AAR->[$r][$c])); 
    }
  }

  # output before table comment lines
  my $line; 
  my $out_line;
  if(defined $bcom_AR) { 
    foreach $line (@{$bcom_AR}) { 
      $out_line = $pfx_com . $line; 
      if($out_line !~ m/\n$/) { $out_line .= "\n"; } # add newline if nec
      print $out_FH1 $out_line; 
      if(defined $out_FH2) { print $out_FH2 $out_line; }
    }
  }

  # create the header separation line
  my $head_sep_line = undef;
  if($hsep ne "") { 
    $head_sep_line = $pfx_head . utl_StringMonoChar($w_A[0]-1, $hsep, undef);
    if($ncol > 1) { $head_sep_line .= $csep; }
    for($c = 1; $c < ($ncol-1); $c++) { 
      $head_sep_line .= utl_StringMonoChar($w_A[$c], $hsep, undef);
      $head_sep_line .= $csep;
    }
    # determine width of final head separation line token
    my $w_final = 0; # max width of final header column over all rows 
    for($r = 0; $r < $nrow_head; $r++) { 
      $w_final = utl_Max($w_final, length($head_AAR->[$r][($ncol-1)]));
    }
    $head_sep_line .= utl_StringMonoChar($w_final, $hsep, undef) . "\n";
  }
  else { # $hsep eq ""
    $head_sep_line = $pfx_head . "\n" 
  }
  my $empty_sep_line = ($empty_flag) ? ($pfx_com . "\n") : $head_sep_line;
  
  # output header lines
  if(defined $head_AAR) { 
    _ofile_helper_output_table_data($head_AAR, \@w_A, \@cljust_A, $pfx_head, $csep, $empty_sep_line, $out_FH1, $out_FH2);
  }

  # output header sep line
  if($nrow_head > 0) { 
    print $out_FH1 $head_sep_line;
    if(defined $out_FH2) { print $out_FH2 $head_sep_line; }
  }

  # output data lines
  _ofile_helper_output_table_data($data_AAR, \@w_A, \@cljust_A, $pfx_data, $csep, $empty_sep_line, $out_FH1, $out_FH2);

  # output after table comment lines
  if(defined $acom_AR) { 
    foreach $line (@{$acom_AR}) { 
      $out_line = $pfx_com . $line; 
      if($out_line !~ m/\n$/) { $out_line .= "\n"; } # add newline if nec
      print $out_FH1 $out_line; 
      if(defined $out_FH2) { print $out_FH2 $line; }
    }
  }

  return;
}
    
#################################################################
# Subroutine:  _ofile_helper_output_table_data()
# Incept:      EPN, Fri Apr  5 06:59:23 2019
#
# Purpose:     Helper subroutine for ofile_TableHumanOutput()
#              that prints 2D array text. (Called separately for
#              header text and data text.)
#              This subroutine does not validate anything (e.g.
#              that all rows have same number columns) -- the 
#              caller should have already done that.
#
# Arguments: 
#   $AAR:        ref to 2D array of output text [0..$nrow-1][0..$ncol-1], 
#                must be defined
#   $w_AR:       ref to array of widths per column
#   $lj_AR:      ref to '1'/'0' array of indicating if a column is left justified
#   $pfx:        prefix for all lines, "" for none
#   $csep:       column separation string, typically "  "
#   $empty_line: line to output for empty arrays
#   $FH1:        output file handle 1
#   $FH2:        output file handle 2
# 
# Returns:       Nothing.
#
# Dies: Never. 
################################################################# 
sub _ofile_helper_output_table_data { 
  my $nargs_expected = 8;
  my $sub_name = "_ofile_helper_output_table_data";
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($AAR, $w_AR, $lj_AR, $pfx, $csep, $empty_line, $FH1, $FH2) = (@_);

  my $nrow = scalar(@{$AAR});
  my $r; # counter over rows
  my $c; # counter over columns
  my $out_line = undef; # output line
  for($r = 0; $r < $nrow; $r++) { 
    my $ncol = scalar(@{$AAR->[$r]});
    if($ncol > 0) { 
      $out_line = sprintf("%-*s", $w_AR->[0], $pfx . $AAR->[$r][0]);
      if($ncol > 1) { 
        $out_line .= $csep; 
        for($c = 1; $c < $ncol; $c++) { 
          if($lj_AR->[$c]) { $out_line .= sprintf("%-*s", $w_AR->[$c], $AAR->[$r][$c]); }
          else             { $out_line .= sprintf("%*s",  $w_AR->[$c], $AAR->[$r][$c]); }
          if($c < ($ncol-1)) { $out_line .= $csep; }
        }
      }
      $out_line .= "\n";
      print $FH1 $out_line;
      if(defined $FH2) { print $FH2 $out_line; }
    }
    else { # empty array, print empty line
      print $FH1 $empty_line;
      if(defined $FH2) { print $FH2 $empty_line; }
    }
  }

  return;
}


####################################################################
# the next line is critical, a perl module must return a true value
return 1;
####################################################################
