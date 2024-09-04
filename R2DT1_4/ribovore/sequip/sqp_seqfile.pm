#!/usr/bin/perl
#
# sqp_seqfile.pm
# Eric Nawrocki
# EPN, Wed Apr  3 06:13:49 2019 [incept, in vadr]
# EPN, Tue Jul  2 11:45:46 2019 [migrated from vadr's epn-seqfile.pm (as of commit 69b003d)]]
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
#####################################################################
#
# List of subroutines:
# 
#################################################################
# Subroutine : sqf_EslSeqstatOptAParse()
# Incept:      EPN, Wed Dec 14 16:16:22 2016 [ribo.pm]
#
# Purpose:     Parse an esl-seqstat -a output file.
#              
# Arguments: 
#   $seqstat_file:  file to parse
#   $seq_name_AR:   REF to array of sequences in order to fill here
#   $seq_len_HR     REF to hash of sequence names to fill here
#   $FH_HR:         REF to hash of file handles, including "cmd"
#
# Returns:     Total number of nucleotides read (summed length of all sequences). 
# 
# Dies:        If the sequence file has two sequences with identical names.
#              Error message will list all duplicates.
#              If no sequences were read.
#
################################################################# 
sub sqf_EslSeqstatOptAParse { 
  my $nargs_expected = 4;
  my $sub_name = "sqf_EslSeqstatOptAParse";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($seqstat_file, $seq_name_AR, $seq_len_HR, $FH_HR) = @_;

  open(IN, $seqstat_file) || ofile_FileOpenFailure($seqstat_file, $sub_name, $!, "reading", $FH_HR);

  my $nread = 0;            # number of sequences read
  my $tot_length = 0;       # summed length of all sequences
  my $seq_name;             # a sequence name
  my $length;               # length of a seq
  my %seq_exists_H = ();    # key is name of a read sequence, value is always '1'
  my %seq_dups_H = ();      # key is a sequence name that exists more than once in seq file, value is number of occurences
  my $at_least_one_dup = 0; # set to 1 if we find any duplicate sequence names
  @{$seq_name_AR} = ();
  %{$seq_len_HR} = ();

  # parse the seqstat -a output 
  # sequences must have non-empty names (else esl-seqstat call would have failed)
  # lengths must be >= 0 (lengths of 0 are okay)
  while(my $line = <IN>) { 
    # = lcl|dna_BP331_0.3k:467     1232 
    # = lcl|dna_BP331_0.3k:10     1397 
    # = lcl|dna_BP331_0.3k:1052     1414 
    chomp $line;
    #print $line . "\n";
    if($line =~ /^\=\s+(\S+)\s+(\d+)/) { 
      ($seq_name, $length) = ($1, $2);

      if(exists($seq_exists_H{$seq_name})) { 
        if(exists($seq_dups_H{$seq_name})) { 
          $seq_dups_H{$seq_name}++; 
        }
        else { 
          $seq_dups_H{$seq_name} = 2;
        }
        $at_least_one_dup = 1;
      }
      else { 
        $seq_exists_H{$seq_name} = 1; 
      }

      push(@{$seq_name_AR}, $seq_name);
      $seq_len_HR->{$seq_name} = $length;
      $tot_length += $length;
      $nread++;
    }
  }
  close(IN);
  if($nread == 0) { 
    ofile_FAIL("ERROR in $sub_name, did not read any sequence lengths in esl-seqstat file $seqstat_file, did you use -a option with esl-seqstat", 1, $FH_HR);
  }
  if($at_least_one_dup) { 
    my $i = 1;
    my $die_string = "\nERROR, not all sequences in input sequence file have a unique name. They must.\nList of sequences that occur more than once, with number of occurrences:\n";
    foreach $seq_name (sort keys %seq_dups_H) { 
      $die_string .= "\t($i) $seq_name $seq_dups_H{$seq_name}\n";
      $i++;
    }
    $die_string .= "\n";
    ofile_FAIL($die_string, 1, $FH_HR);
  }

  return $tot_length;
}

#################################################################
# Subroutine : sqf_FastaFileRemoveDescriptions()
# Incept:      EPN, Mon Mar 25 11:30:16 2019
#
# Purpose:     Given an FASTA file, create a new one that is
#              identical but with sequence descriptions removed.
#              DOES NOT VALIDATE INPUT FILE IS IN PROPER FASTA FORMAT.
#              
# Arguments: 
#   $in_file:        input FASTA file
#   $out_file:       output FASTA file, with sequence descriptions removed
#   $ofile_info_HHR: ref to the ofile info 2D hash, can be undef
# 
# Returns:     void
#
# Dies:        If unable to open $in_file for reading or $out_file
#              for writing.
#
################################################################# 
sub sqf_FastaFileRemoveDescriptions {
  my $nargs_expected = 3;
  my $sub_name = "sqf_FastaFileRemoveDescriptions()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($in_file, $out_file, $ofile_info_HHR) = (@_);

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef; # for convenience

  open(IN,       $in_file)  || ofile_FileOpenFailure($in_file,  $sub_name, $!, "reading", $FH_HR);
  open(OUT, ">", $out_file) || ofile_FileOpenFailure($out_file, $sub_name, $!, "writing", $FH_HR);

  while(my $line = <IN>) { 
    chomp $line;
    if($line =~ /^>(\S+)/) { 
      print OUT ">" . $1 . "\n";
    }
    else { 
      print OUT $line . "\n";
    }
  }
  close(IN);
  close(OUT);

  return;
}

#################################################################
# Subroutine: sqf_FeatureTableParse()
# Incept:     EPN, Mon May 20 09:36:09 2019
#
# Synopsis: Parse a INSDC feature table format file.
#
# Arguments:
#  $infile:         feature table file to parse
#  $ftr_info_HAHR:  feature information, filled here
#                   1D key: accession
#                   2D:     feature index
#                   3D key: qualifer, value: qualifier value
#  $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if we have trouble parsing the file
#             if $allow_incomplete is '1' and we read an incomplete feature
#
# Reference: https://www.ncbi.nlm.nih.gov/Sequin/table.html
#################################################################
sub sqf_FeatureTableParse { 
  my $sub_name = "sqf_FeatureTableParse";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $ftr_info_HAHR, $FH_HR) = @_;

  open(IN, $infile) || ofile_FileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

  my $long_accver = undef;   # full accession, e.g. ref|NC_000883.2|
  my $acc         = undef;   # accession, e.g. 'NC_000883.2'
  my $ver         = undef;   # version, e.g. '2'
  my $qname       = undef;   # a qualifier name,  e.g. 'organism'
  my $qval        = undef;   # a qualifier value, e.g. 'Paramecia bursaria Chlorella virus 1'
  my $feature     = undef;   # a feature name, e.g. "CDS", or "gene"
  my $ftr_idx     = -1;      # number of features read for current sequence
  my $coords      = undef;   # coordinates
  my $strand      = undef;   # strand of last segment read 
  my $trunc5      = undef;   # '1' if current feature is 5' truncated (start carrot, e.g. NC_031327:"<3281..4207")
  my $trunc3      = undef;   # '1' if current feature is 5' truncated (start carrot, e.g. "3281..>4207")
  my $line_idx    = 0;       # count of number of lines read in ftable
  my $prv_was_accn           = 0; # set to '1' if previous line was an accession line
  my $prv_was_coords_feature = 0; # set to '1' if previous line was a coordinates line with a feature name
  my $prv_was_coords_only    = 0; # set to '1' if previous line was a coordinates line without a feature name
  my $prv_was_quals          = 0; # set to '1' if previous line was a qualifier_name qualifier value line

  while(my $line = <IN>) { 
    $line_idx++;
    chomp $line;
    if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
    if($line =~ m/\w/) { 
      # parse each of the 4 line types differently
      # -------------------------------------------------------
      if($line =~ /^\>Feature\s+(\S+)$/) { 
        # ACCESSION LINE
        # example:
        #>Feature ref|NC_001359.1|    
        # or
        #>Feature anyseqname
        $long_accver = $1;
        # accession line can occur after any other line type, so we don't have to check if line order makes sense for this case

        # if our previous line was coords_feature or coords_only, we need to store the feature from that previous line
        if(($prv_was_coords_feature) || ($prv_was_coords_only)) { 
          $ftr_idx++;
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",   $feature, $FH_HR);
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "coords", $coords,  $FH_HR);
          if($trunc5) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc5", 1, $FH_HR); }
          if($trunc3) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc3", 1, $FH_HR); }
        }

        # determine accession and version, e.g. NC_001359.1 in above example
        if($long_accver =~ /[^\|]*\|([^\|]+)\.(\d+)\|/) { 
          $acc = $1;
          $ver = $2;
        }
        else { 
          $acc = $long_accver;
          #ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, unable to parse header line:\n$line\n", 1, $FH_HR);
        }
        @{$ftr_info_HAHR->{$acc}} = (); # initialize array
        $feature = undef; 
        $coords  = undef;
        $trunc5  = undef;
        $trunc3  = undef;
        $strand  = undef;
        $ftr_idx = -1;

        # update '$prv_*' values that we use to make sure line order makes sense
        $prv_was_accn           = 1;
        $prv_was_coords_feature = 0;
        $prv_was_coords_only    = 0;
        $prv_was_quals          = 0;
        #printf("set prv_was_accn\n");
      }
      # -------------------------------------------------------
      elsif($line =~ /^(\<?)(\d+\^?)\t(\>?)(\d+)\t(\S+)$/) { 
        # COORDINATES LINE WITH A FEATURE NAME (coords_feature)
        # example:
        # 230   985     gene
        my ($start_carrot, $start_coord, $stop_carrot, $stop_coord, $tmp_feature) = ($1, $2, $3, $4, $5);
        # coords_feature line can occur after any other line type, so we don't have to check if line order makes sense for this case

        # if our previous line was coords_feature or coords_only, we need to store the feature from that previous line
        if(($prv_was_coords_feature) || ($prv_was_coords_only)) { 
          $ftr_idx++;
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",   $feature, $FH_HR);
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "coords", $coords,  $FH_HR);
          if($trunc5) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc5", 1, $FH_HR); }
          if($trunc3) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc3", 1, $FH_HR); }
        }
        $feature = $tmp_feature;
        $trunc5 = ($start_carrot eq "<") ? 1 : 0;
        $trunc3 = ($stop_carrot  eq ">") ? 1 : 0;

        $strand = ($start_coord <= $stop_coord) ? "+" : "-"; # sets + strand for single nt
        $coords = $start_coord . ".." . $stop_coord . ":" . $strand;

        # update '$prv_*' values that we use to make sure line order makes sense
        $prv_was_accn           = 0;
        $prv_was_coords_feature = 1;
        $prv_was_coords_only    = 0;
        $prv_was_quals          = 0;
        #printf("set prv_was_coords_feature\n");
      }
      # -------------------------------------------------------
      elsif($line =~ /^(\<?)(\d+)\t(\>?)(\d+)$/) {  
        # COORDINATES LINE WITHOUT A FEATURE NAME (coords_only) 
        # example:
        # 154   183
        my ($start_carrot, $start_coord, $stop_carrot, $stop_coord) = ($1, $2, $3, $4);

        # a coords_only line can only occur after a coords_feature line or coords_only line, 
        # check to make sure that's the case
        if($prv_was_coords_feature || # previous line was a coords line with a feature (common)
           $prv_was_coords_only) {    # previous line was a coords line without a feature (common)
          # line order makes sense, keep going...

          if($start_carrot ne "") { 
            ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read start_carrot indicating 5' truncated feature from coordinate line that wasn't the first coordinate line for the feature, line:\n$line\n", 1, $FH_HR);
          }
          $trunc3 = ($stop_carrot eq ">") ? 1 : 0;

          if($start_coord == $stop_coord) { 
            ; # single nt, use strand of previous segment
          }
          else { 
            $strand = ($start_coord < $stop_coord) ? "+" : "-";
          }
          $coords .= "," . $start_coord . ".." . $stop_coord . ":" . $strand;

          # update '$prv_*' values that we use to make sure line order makes sense
          $prv_was_accn           = 0;
          $prv_was_coords_feature = 0;
          $prv_was_coords_only    = 1;
          $prv_was_quals          = 0;
          #printf("set prv_was_coords_only\n");
        }
        else { # line order is unexpected
          ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, unexpected line order (coords_only line), line:\n$line\n", 1, $FH_HR);
        }
      }
      # -------------------------------------------------------
      elsif(($line =~ /^\t\t\t[^\t]+\t[^\t]+$/) || 
            ($line =~ /^\t\t\t[^\t]+$/)) { 
        # QUALIFIER LINE
        # examples:
        #gene       AR1
        #locus_tag  PhyvvsAgp1

        # before parsing it, do two sanity checks
        if(! defined $coords)  { 
          ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read qualifier line but coords is not yet defined, line:\n$line\n", 1, $FH_HR);
        }
        if(! defined $feature)  { 
          ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read qualifier line but feature is not yet defined, line:\n$line\n", 1, $FH_HR);
        }
        # does line order make sense?
        if($prv_was_coords_feature || 
           $prv_was_coords_only    ||
           $prv_was_quals) { 
          # line order makes sense, keep going...
          if(! $prv_was_quals) { 
            # first quals line for this feature, store the feature
            $ftr_idx++;
            sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",   $feature, $FH_HR);
            sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "coords", $coords,  $FH_HR);
            if($trunc5) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc5", 1, $FH_HR); }
            if($trunc3) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc3", 1, $FH_HR); }
          }
          # parse the line
          if($line =~ /^\t\t\t([^\t]+)\t([^\t]+)$/) { 
            ($qname, $qval) = ($1, $2);
          }
          elsif($line =~ /^\t\t\t([^\t]+)$/) { 
            ($qname, $qval) = ($1, "");
          }
          else { 
            ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, unable to parse qualifier line, line:\n$line\n", 1, $FH_HR);
          }
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qname, $qval, $FH_HR);
        } # end of 'if() that checks line order makes sense
        else { # unexpected line order
          ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, unexpected line order (quals line), line:\n$line\n", 1, $FH_HR);
        }          

        # update '$prv_*' values that we use to make sure line order makes sense
        $prv_was_accn           = 0;
        $prv_was_coords_feature = 0;
        $prv_was_coords_only    = 0;
        $prv_was_quals          = 1;
        #printf("set prv_was_quals\n");
      }
      # -------------------------------------------------------
      else { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, unable to parse line, line:\n$line\n", 1, $FH_HR);
      }
      # -------------------------------------------------------
    }
  }
  if(! defined $acc) { 
    ofile_FAIL("ERROR in $sub_name, problem parsing $infile, did not read any accession lines\n", 1, $FH_HR);
  }

  # add the final feature, if we haven't already, we can tell based on previous line type
  if(($prv_was_coords_feature) || ($prv_was_coords_only)) { 
    $ftr_idx++;
    sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",   $feature, $FH_HR);
    sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "coords", $coords,  $FH_HR);
    if($trunc5) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc5", 1, $FH_HR); }
    if($trunc3) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc3", 1, $FH_HR); }
  }

  return;
}

#################################################################
# Subroutine: sqf_GenbankParse()
# Incept:     EPN, Tue Mar 12 14:04:14 2019
#
# Synopsis: Parse a GenBank format file.
#
# Arguments:
#  $infile:        GenBank file to parse
#  $seq_info_HHR:  sequence info, filled here, can be undef
#                  1D key: accession
#                  2D key: "len", "ver", "def", "seq"
#  $ftr_info_HAHR: feature information, filled here
#                  1D key: accession
#                  2D:     feature index
#                  3D key: qualifer, value: qualifier value
#  $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if we have trouble parsing the file
#
# Reference: https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
#            http://www.insdc.org/documents/feature-table
#################################################################
sub sqf_GenbankParse { 
  my $sub_name = "sqf_GenbankParse";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $seq_info_HHR, $ftr_info_HAHR, $FH_HR) = @_;

  my $line_idx  = 0;     # line index of input file
  my $acc       = undef; # accession, read from LOCUS line
  my $tmp_acc   = undef; # accession, read from ACCESSION or VERSION line
  my $len       = undef; # length, read from LOCUS line
  my $def       = undef; # seq definition, read from DEFINITION line
  my $ver       = undef; # sequence version, read from VERSION line
  my $feature   = undef; # a feature   read from a feature/location line in the FEATURES section
  my $location  = undef; # a location  read from a feature/location line in the FEATURES section
  my $qualifier = undef; # a qualifier read from a qualifier/value  line in the FEATURES section
  my $value     = undef; # a value     read from a qualifier/value  line in the FEATURES section
  my $seq       = undef; # sequence, read from the ORIGIN section
  my $seqline   = undef; # single line of sequence
  my $seq_idx   = 0;     # number of sequences read
  my $ftr_idx   = -1;    # number of features read for current sequence
  my $line      = undef; # a line

  open(IN, $infile) || ofile_FileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

  $line = <IN>; 
  while(defined $line) { 
    chomp $line; $line_idx++;
    if($line =~ /^LOCUS\s+(\S+)\s+(\d+)/) { 
      #LOCUS       NC_039477               7567 bp    RNA     linear   VRL 22-FEB-2019
      if((defined $acc) || (defined $len)) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple LOCUS lines for single record ($acc), line:\n$line\n", 1, $FH_HR);
      }
      ($acc, $len) = ($1, $2);
      # initialize the array of hashes for this accession's features
      if(defined $ftr_info_HAHR->{$acc}) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, trying to add feature info for accession $acc, but it already exists, line:\n$line\n", 1, $FH_HR);
      }
      @{$ftr_info_HAHR->{$acc}} = ();
      $line = <IN>; 
    }
    elsif($line =~ /^DEFINITION\s+(.*)$/) { 
      #DEFINITION  Norovirus GII isolate strain Hu/GBR/2016/GII.P16-GII.4_Sydney/226,
      #            complete genome.
      if(defined $def) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple DEFINITION lines for single record ($acc), line:\n$line\n", 1, $FH_HR);
      }
      $def = $1;
      # read remainder of the definition (>= 0 lines)
      $line = <IN>; 
      while((defined $line) && ($line =~ /^\s+(.+)$/)) {
        chomp $line; $line_idx++;
        $def .= $1;
        $line = <IN>; 
      }
    }
    elsif($line =~ /^ACCESSION\s+(\S+)$/) { 
      # ACCESSION   NC_039477
      # verify this matches what we read in the LOCUS line
      $tmp_acc = $1;
      if((! defined $acc) || ($tmp_acc ne $acc)) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, accession mismatch for $tmp_acc, line:\n$line\n", 1, $FH_HR);
      }
      $line = <IN>;
    }
    elsif($line =~ /^VERSION\s+(\S+)$/) { 
      #VERSION     NC_039477.1
      if(defined $ver) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple VERSION lines for single record ($acc), line:\n$line\n", 1, $FH_HR);
      }
      # verify this matches what we read in the LOCUS line
      $ver = $1;
      $tmp_acc = $ver;
      seq_StripVersion(\$tmp_acc);
      if((! defined $acc) || ($tmp_acc ne $acc)) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, version/accession mismatch for $tmp_acc, line:\n$line\n", 1, $FH_HR);
      }
      $line = <IN>;
    }
    elsif($line =~ /^FEATURES\s+Location\/Qualifiers$/) { 
      # parse the features and then the sequence
      # FEATURES section
      # 4 types of line:
      # feature/location line
      #        example1:      gene            5..5104
      #        example2:      misc_feature    join(2682..2689,1..2)
      #        example3:      misc_feature    join(161990..162784,complement(88222..88806),complement(86666..87448))
      # qualifier, no value, single line
      #        example: /ribosomal_slippage
      # qualifier/value type A, single line of controlled vocab or enumerated value
      #        example: /codon_start=1
      # qualifier/value type B, single line of free text
      #        example1: /gene="ORF1"
      # qualifier/value type C, multiple lines of free text
      #       example1: /translation="MMMASKDVVPTAASSENANNNSSIKSRLLARLKGSGGATSPPNS
      #       example1: IKITNQDMALGLIGQVPAPKATSVDVPKQQRDRPPRTVAEVQQNLRWTERPQDQNVKT
      #       example1: QNVIDPWIRNNFVQAPGGEFTVSPRNAPGEILWSAPLGPDLNPYLSHLARMYNGYAGG
      #       example1: IPPNGYFRFDSWVNQFYTLAPMGNGTGRRRVV"
      #
      #       example2: /note="The mature peptides were added by the NCBI staff
      #       example2: following publications Liu et all (1996, 1999), that
      #       example2: (tentatively) determined processing map by site-directed
      #       example2: mutagenesis or by direct sequencing of the cleavage
      #       example2: products for closely related Southampton virus (a
      #       example2: Norwalk-like virus); ORF1; sequence homologies to 2C
      #       example2: helicase, 3C protease, and 3D RNA-dependent RNA polymerase
      #       example2: of picornavirus"
      #
      if($ftr_idx != -1) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple FEATURES lines for single record ($acc), line:\n$line\n", 1, $FH_HR);
      }
      $line = <IN>;
      while((defined $line) && ($line !~ /^ORIGIN/)) { 
        chomp $line; $line_idx++;
        if($line =~ /^\s+\/(\S+)\=(.+)$/) { # first token must start with '/'
          my ($save_qualifier, $save_value) = ($1, $2);
          # determine type of qualifier value
          # A: single line of controlled vocabulary or enumerated value, does not start with '"'
          # B: single line of free text: starts and ends with single '"'
          # C: first line of multiple lines of free text: starts with single '"', does not end with single '"'
          if($save_value  =~ /^[^\"].*/) { 
            # A: single line of controlled vocabulary or enumerated value, does not start with '"'
            # example:
            # /codon_start=1
            ; # do nothing, we just need this if to catch a valid line (so we don't get to the else below that leads to a failure)
          }
          elsif(($save_value =~ /^\"[^\"]\"$/) || # single non-'"' character in between two '"' characters
                ($save_value =~ /^\"[^\"].*[^\"]\"$/)) { # '"' then a non-'"', then >= 0 characters, then a non-'"' character then ends with a '"' character
            # B: single line of free text: starts and ends with single '"'
            # example:
            # /product="nonstructural polyprotein"
            ; # do nothing, we just need this elsif to catch a valid line (so we don't get to the else below that leads to a failure)
          }
          elsif(($save_value =~ /^\"[^\"]$/) || # single non-'"' character after '"' character
                ($save_value =~ /^\"[^\"].*[^\"]$/)) { # '"' then a non-'"', then >= 0 characters, then ends with a non-'"' character
            # C: first line of multiple lines of free text: starts with single '"', does not end with single '"'
            # example:
            # /note="The mature peptides were added by the NCBI staff
            # 
            # continue to read lines until we read one that ends with a single '"'
            #
            # example of remaining lines in value:
            #  following publications Liu et all (1996, 1999), that
            #  (tentatively) determined processing map by site-directed
            #  mutagenesis or by direct sequencing of the cleavage
            #  products for closely related Southampton virus (a
            #  Norwalk-like virus); ORF1; sequence homologies to 2C
            #  helicase, 3C protease, and 3D RNA-dependent RNA polymerase
            #  of picornavirus"
            $line = <IN>; chomp $line; $line_idx++;
            while((defined $line) && ($line !~ /[^\"]\"$/)) { # ends with a non-'"' character followed by '"' character
              $line =~ s/^\s+//; # remove leading whitespace
              # add a single white space if there is at least one white space character already exising in $save_value
              if($save_value =~ m/\s/) { $save_value .= " "; }
              $save_value .= $line;
              $line = <IN>; chomp $line; $line_idx++;
            }
            # add final line we read
            $line =~ s/^\s+//; # remove leading whitespace
            # add a single white space if there is at least one white space character already exising in $save_value
            if($save_value =~ m/\s/) { $save_value .= " "; }
            $save_value .= $line;
          }
          else { 
            ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, in FEATURES section read unparseable qualifier value line:\n$line\n", 1, $FH_HR);
          }
          if(defined $value) { # we are finished with previous value
            sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qualifier, $value, $FH_HR);
          }
          ($qualifier, $value) = ($save_qualifier, $save_value);
        }
        elsif($line =~ /^\s+\/([^\s\=]+)$/) { # first token must start with '/'
          my ($save_qualifier, $save_value) = ($1, "");
          if(defined $value) { # we are finished with previous value
            sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qualifier, $value, $FH_HR);
            ($qualifier, $value) = (undef, undef);
          }
          # qualifier, no value, single line
          #        example: /ribosomal_slippage
          ($qualifier, $value) = ($save_qualifier, $save_value);
        }
        elsif($line =~ /^\s+(\S+)\s+(\S+)$/) { 
          if(defined $value) { # we are finished with previous value
            sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qualifier, $value, $FH_HR);
            ($qualifier, $value) = (undef, undef);
          }
          # feature/location line, examples:
          #   gene            5..5104
          ($feature, $location) = ($1, $2);
          $ftr_idx++;
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",     $feature,  $FH_HR);
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "location", $location, $FH_HR);
        }
        $line = <IN>;
      }
      if(! defined $line) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, expected to read ORIGIN line after FEATURES but did not\n", 1, $FH_HR);
      }
      # if we get here we just read the ORIGIN line
      # first store final qualifier/value
      if(defined $value) { 
        sqf_StoreQualifierValue($ftr_info_HAHR->{$acc}, $ftr_idx, $qualifier, $value, $FH_HR);
      }
      # parse the ORIGIN sequence
      $line = <IN>;
      # sanity check
      if(defined $seq) { 
        ofile_FAIL("ERROR in $sub_name, read multiple ORIGIN lines for single record ($acc), line:\n$line\n", 1, $FH_HR);
      }
      $seq = "";
      while((defined $line) && ($line !~ /^\/\/$/)) { 
        chomp $line; $line_idx++;
        # sequence lines
        # examples:
        # 7501 gtcacgggcg taatgtgaaa agacaaaact gattatcttt ctttttcttt agtgtctttt
        # 7561 aaaaaaa
        if($line =~ /^\s+\d+\s+(.+)$/) { 
          $seqline = $1;
          $seqline =~ s/\s+//g; # remove spaces
          $seq .= $seqline;
        }
        $line = <IN>;
      }
      if(! defined $line) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, expected to find a // line after ORIGIN but did not, line $line_idx\n", 1, $FH_HR);
      }
      # if we get here we just read the // line
      # we are finished with this sequence, store the information
      if(! defined $acc) { ofile_FAIL(        "ERROR in $sub_name, failed to read accession, line: $line_idx\n", 1, $FH_HR); }
      if(! defined $len) { ofile_FAIL(sprintf("ERROR in $sub_name, failed to read length (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), 1, $FH_HR); }
      if(! defined $ver) { ofile_FAIL(sprintf("ERROR in $sub_name, failed to read version (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), 1, $FH_HR); }
      if(! defined $def) { ofile_FAIL(sprintf("ERROR in $sub_name, failed to read definition (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), 1, $FH_HR); }
      if(! defined $seq) { ofile_FAIL(sprintf("ERROR in $sub_name, failed to read sequence (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), 1, $FH_HR); }

      # store sequence info
      if(defined $seq_info_HHR) { 
        %{$seq_info_HHR->{$acc}} = ();
        $seq_info_HHR->{$acc}{"len"} = $len;
        $seq_info_HHR->{$acc}{"ver"} = $ver;
        $seq_info_HHR->{$acc}{"def"} = $def;
        $seq_info_HHR->{$acc}{"seq"} = $seq;
      }

      # reset variables
      $seq = undef;
      $len = undef;
      $acc = undef;
      $ver = undef;
      $def = undef;
      $feature   = undef;
      $location  = undef;
      $qualifier = undef;
      $value     = undef;

      $seq_idx++;
      $ftr_idx = -1;
      $line = <IN>;
    } # end of 'elsif($line =~ /^FEATURES\s+Location\/Qualifiers$/) {' 
    else { 
      # not a line we will parse, read the next line
      $line = <IN>;
    }
  }

  if($seq_idx == 0) { 
    ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, failed to read any sequence data\n", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: sqf_StoreQualifierValue()
# Incept:     EPN, Wed Mar 13 09:42:22 2019
#
# Synopsis: Store a genbank qualifier and value.
#
# Arguments:
#  $ftr_info_AHR: REF to the array of hashes to store data in
#  $ftr_idx:      feature index
#  $qualifier:    qualifier
#  $value:        qualifier value
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#
# Returns:    '1' if $ftr_info_AHR->[$ftr_idx]{$qualifier} created
#             '0' if $ftr_info_AHR->[$ftr_idx]{$qualifier} exists upon entering function
#
# Dies:       If $value includes the string ":GPSEP:, which we use 
#             to separate multiple qualifier values for the same qualifier.
#             
#################################################################
sub sqf_StoreQualifierValue { 
  my $sub_name = "sqf_StoreQualifierValue";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($ftr_info_AHR, $ftr_idx, $qualifier, $value, $FH_HR) = @_;

  if($value =~ /\:GPSEP\:/) { 
    ofile_FAIL("ERROR in $sub_name, qualifier value $value includes the special string :GPSEP:, this is not allowed", 1, $FH_HR);
  }
  if($value =~ /^\"GBNULL\"$/) { 
    ofile_FAIL("ERROR in $sub_name, qualifier value $value is GBNULL, this is not allowed", 1, $FH_HR);
  }

  # remove leading and trailing " in the value, if they exist
  # GenBank format uses "" as a substitute for " in these strings
  $value =~ s/^\"//;
  $value =~ s/\"$//;

  if($value eq "") { $value = "GBNULL"; }

  # printf("in $sub_name q: $qualifier v: $value\n");
  if(! defined ($ftr_info_AHR->[$ftr_idx])) { 
    %{$ftr_info_AHR->[$ftr_idx]} = (); 
  }
  if(! defined $ftr_info_AHR->[$ftr_idx]{$qualifier}) { 
    $ftr_info_AHR->[$ftr_idx]{$qualifier} = $value;
  }
  else { 
    $ftr_info_AHR->[$ftr_idx]{$qualifier} .= ":GBSEP:" . $value;
  }

  return;
}

#################################################################
# Subroutine: sqf_BlastDbCreate
# Incept:     EPN, Mon Mar 18 09:40:28 2019
# 
# Purpose:    Create a blast database from a fasta file of either
#             type 'prot' or 'nucl'
#
# Arguments:
#   $makeblastdb:    path to 'makeblastdb' executable
#   $dbtype:         'prot' for protein or 'nucl' for nucleotide
#   $fa_file:        FASTA file of protein sequences to make blast db from
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#                    
# Dies: if $dbtype is not 'nucl' or 'prot' or if makeblastdb fails
# 
# Returns:    void
#
#################################################################
sub sqf_BlastDbCreate {
  my $sub_name = "sqf_BlastDbCreate";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($makeblastdb, $dbtype, $fa_file, $opt_HHR, $FH_HR) = @_;

  if(($dbtype ne "prot") && ($dbtype ne "nucl")) { 
    ofile_FAIL("ERROR in $sub_name, dbtype is $dbtype, expected to be 'prot' or 'nucl'", 1, $FH_HR);
  }

  utl_RunCommand($makeblastdb . " -in $fa_file -dbtype $dbtype > /dev/null", opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}

#################################################################
# Subroutine: sqf_EslTranslateCdsToFastaFile()
# Incept:     EPN, Thu Mar 14 12:30:28 2019
# 
# Purpose:    Use esl-translate to translate a fasta file with
#             CDS sequences pertaining to the CDS features in 
#             @{$ftr_info_AHR} into fasta protein files.
#
# Arguments:
#   $out_FH:         output file handle to print to 
#   $esl_translate:  path to esl-translate executable
#   $cds_fa_file:    fasta file with CDS sequences
#   $out_root:       string that is the 'root' for naming output files
#   $ftr_info_AHR:   REF to the feature info, pre-filled
#   $opt_HHR:        command line options
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#                    
# Returns: void
#
# Dies:    if we have trouble fetching a sequence
#
#################################################################
sub sqf_EslTranslateCdsToFastaFile { 
  my $sub_name = "sqf_EslTranslateCdsToFastaFile";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_FH, $esl_translate, $cds_fa_file, $out_root, $ftr_info_AHR, $opt_HHR, $FH_HR) = @_;

  my $tmp1_translate_fa_file  = $out_root . ".cds.esl-translate.1.fa";
  my $tmp2_translate_fa_file  = $out_root . ".cds.esl-translate.2.fa";
  my $tmp1_translate_ssi_file = $out_root . ".cds.esl-translate.1.fa.ssi";
  my $tmp2_translate_ssi_file = $out_root . ".cds.esl-translate.2.fa.ssi";
  if(-e $tmp1_translate_ssi_file) { unlink $tmp1_translate_ssi_file; }
  if(-e $tmp2_translate_ssi_file) { unlink $tmp2_translate_ssi_file; }

  my $c_opt = "";
  if((opt_IsUsed("--ttbl", $opt_HHR)) && (opt_Get("--ttbl", $opt_HHR) != 1)) { 
    $c_opt = "-c " . opt_Get("--ttbl", $opt_HHR);
  }
  my $translate_cmd = "$esl_translate $c_opt -M -l 3 --watson $cds_fa_file > $tmp1_translate_fa_file";
  utl_RunCommand($translate_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  # go through output fasta file and rewrite names, so we can fetch 
  open(IN,       $tmp1_translate_fa_file) || ofile_FileOpenFailure($tmp1_translate_fa_file, $sub_name, $!, "reading", $FH_HR);
  open(OUT, ">", $tmp2_translate_fa_file) || ofile_FileOpenFailure($tmp2_translate_fa_file, $sub_name, $!, "writing", $FH_HR);
  while(my $line = <IN>) { 
    if($line =~ m/^\>/) { 
      #>orf58 source=NC_039477.1/5..5104:+ coords=1..5097 length=1699 frame=1  
      chomp $line;
      if($line =~ /^\>orf\d+\s+(source\=\S+)\s+(coords\=\S+)\s+length\=\d+\s+frame\=\S+/) { 
        # rename as 'source=NC_039477.1/5..5104:+,coords=1..5097'
        print OUT (">" . $1 . "," . $2 . "\n");
      }
      else { 
        ofile_FAIL("ERROR in $sub_name, problem parsing esl-translate output file $tmp1_translate_fa_file, line:\n$line\n", 1, $FH_HR);
      }
    }
    else { 
      print OUT $line; 
    }
  }
  close(IN);
  close(OUT);

  # $tmp2_translate_fa_file now includes renamed translated sequences from esl-translate
  # fetch expected translated seqs and print to $out_FH
  my $cds_sqfile     = Bio::Easel::SqFile->new({ fileLocation => $cds_fa_file });
  my $protein_sqfile = Bio::Easel::SqFile->new({ fileLocation => $tmp2_translate_fa_file });

  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $seq_idx = 0; $seq_idx < $cds_sqfile->nseq_ssi; $seq_idx++) { 
    my ($seq_name, $seq_length) = $cds_sqfile->fetch_seq_name_and_length_given_ssi_number($seq_idx);
    my $fetch_name = "source=" . $seq_name . ",coords=1.." . ($seq_length - 3); # subtract length of stop codon
    if(! $protein_sqfile->check_seq_exists($fetch_name)) { 
      ofile_FAIL("ERROR in $sub_name, problem translating CDS feature, unable to find expected translated sequence in $tmp2_translate_fa_file:\n\tseq: $seq_name\n\texpected sequence:$fetch_name\n", 1, $FH_HR);
    }
    print $out_FH ">" . $seq_name . "\n";
    print $out_FH seq_SqstringAddNewlines($protein_sqfile->fetch_seq_to_sqstring($fetch_name), 60);
  }
  # remove temporary files unless --keep
  if(! opt_Get("--keep", $opt_HHR)) { 
    utl_FileRemoveUsingSystemRm($tmp1_translate_fa_file, $sub_name, $opt_HHR, $FH_HR);
    utl_FileRemoveUsingSystemRm($tmp2_translate_fa_file, $sub_name, $opt_HHR, $FH_HR);
    utl_FileRemoveUsingSystemRm($tmp2_translate_fa_file . ".ssi", $sub_name, $opt_HHR, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: sqf_FastaWriteSequence()
# Incept:     EPN, Thu Mar 14 06:06:59 2019
#
# Synopsis: Print a sequence to a fasta file.
#
# Arguments:
#  $out_FH:    output file handle
#  $name:      sequence name
#  $def:       sequence definition, can be undef
#  $seq:       sequence string
#  $FH_HR:     REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $name or $seq is undef
#################################################################
sub sqf_FastaWriteSequence {
  my $sub_name = "sqf_FastaWriteSequence";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_FH, $name, $def, $seq, $FH_HR) = @_;

  if(! defined $name) { ofile_FAIL("ERROR in $sub_name, name is undefined", 1, $FH_HR); }
  if(! defined $seq)  { ofile_FAIL("ERROR in $sub_name, name is undefined", 1, $FH_HR); }

  # capitalize and DNAize $seq
  seq_SqstringCapitalize(\$seq);
  seq_SqstringDnaize(\$seq);
  printf $out_FH (">%s%s\n%s", 
                  $name, 
                  (defined $def) ? " " . $def : "",
                  seq_SqstringAddNewlines($seq, 60));
  
  return;
}

#################################################################
# Subroutine: sqf_FastaFileSplitRandomly
# Incept:     EPN, Fri Jul  6 09:56:37 2018
#
# Purpose:    Given a fasta file and a hash with sequence lengths
#             for all sequences in the file, split the file into
#             <n> files randomly, such that each sequence is randomly
#             placed in one of the <n> files with the exception that
#             the first i=1 to <n> sequences are placed in files 
#             1 to <n> (so that each file gets at least one sequence).
#
# Arguments:
#   $fa_file:     the fasta file
#   $seqlen_HR:   ref to hash, key is sequence name, value is sequence length
#   $out_dir:     output directory for placing sequence files
#   $tot_nseq:    total number of sequences
#   $tot_nres:    total number of nucleotides in all sequences
#   $targ_nres:   target number of residues per file
#   $rng_seed:    seed for srand(), to seed RNG, undef to not seed it
#   $FH_HR:       ref to hash of file handles, including "cmd"
# 
# Returns:  Number of files created.
# 
# Dies:     If we trouble parsing/splitting the fasta file
#
#################################################################
sub sqf_FastaFileSplitRandomly { 
  my $sub_name = "sqf_FastaFileSplitRandomly";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  #my $random_number = int(rand(100));
  my ($fa_file, $seqlen_HR, $out_dir, $tot_nseq, $tot_nres, $targ_nres, $rng_seed, $FH_HR) = @_;

  my $do_debug = 0;

  my $in_FH = undef;
  open($in_FH, $fa_file) || ofile_FileOpenFailure($fa_file, $sub_name, $!, "reading", $FH_HR);
  
  if(defined $rng_seed) { srand($rng_seed); }

  # determine number of files to create
  my $nfiles = int($tot_nres / $targ_nres);
  if($nfiles < 1) { $nfiles = 1; }
  # nfiles must be less than or equal to: $nseq and 300
  if($nfiles > $tot_nseq) { $nfiles = $tot_nseq; }
  if($nfiles > 300)       { $nfiles = 300; }
  if($nfiles <= 0) { 
    ofile_FAIL("ERROR in $sub_name, trying to make $nfiles files", 1, $FH_HR);
  }
  $targ_nres = int($tot_nres / $nfiles);

  # We need to open up all output file handles at once, we'll randomly
  # choose which one to print each sequence to. We need to keep track
  # of total length of all sequences output to each file so we know
  # when to close them. Once a file is closed, we won't choose to
  # write to it anymore, using the @map_A array as follows:
  #
  # We define an array @r2f_map_A with an element for each of the $nfiles
  # output files. For each sequence, we randomly choose a number
  # between 0 and $nfiles-1 to pick which output file to write the
  # sequence to. Initially $r2f_map_A[$i] == $i, but when if we close file
  # $i we set $r2f_map_A[$i] to $r2f_map_A[$nremaining-1], then choose a
  # random int between 0 and $nremaining-1. This gets us a random
  # sample without replacement. 
  #
  # @f2r_map_A is the inverse of @r2f_mapA, which we need only so that
  # we can guarantee that each file gets at least 1 sequence.
  # 
  my @r2f_map_A = ();  # map of random index to file number, $r2f_map_A[$ridx] = file number that random choice $ridx pertains to
  my @f2r_map_A = ();  # map of file number to random index, $f2r_map_A[$fidx] = random choice $ridx that file number $fidx pertains to
  my @nres_per_out_A = ();
  my $nres_tot_out = 0;  # total number of sequences output thus far
  my @nseq_per_out_A = ();
  my @out_filename_A = (); # array of file names
  my @out_FH_A = (); # [0..$nfiles-1], the actual open file handles
  my @isopen_A = (); # [0..$i..$nfiles-1], '1' if file $i is open, '0' if it has been closed
  my $nopen = 0; # number of files that are still open
  my $checkpoint_fraction_step = 0.05; # if($do_randomize) we will output update each time this fraction of total sequence has been output
  my $checkpoint_fraction = $checkpoint_fraction_step;
  my $checkpoint_nres = $checkpoint_fraction * $tot_nres;
  my $fidx; # file index of current file in @out_filename_A and file handle in @out_FH_A
  my $nres_this_seq = 0; # number of residues in current file
  
  # variables only used if $do_randomize
  my $ridx; # randomly selected index in @map_A for current sequence
  my $FH; # pointer to current file handle to print to
  my $nseq_remaining = $tot_nseq;
  my $nseq_output    = 0;
  my $fa_file_tail = ofile_RemoveDirPath($fa_file);

  for($fidx = 0; $fidx < $nfiles; $fidx++) { $r2f_map_A[$fidx] = $fidx; }
  for($fidx = 0; $fidx < $nfiles; $fidx++) { $f2r_map_A[$fidx] = $fidx; }
  for($fidx = 0; $fidx < $nfiles; $fidx++) { $nres_per_out_A[$fidx] = 0; }
  for($fidx = 0; $fidx < $nfiles; $fidx++) { $nseq_per_out_A[$fidx] = 0; }
  for($fidx = 0; $fidx < $nfiles; $fidx++) { $out_filename_A[$fidx] = $out_dir . "/" . $fa_file_tail . "." . ($fidx+1); } 

  # open up all output file handles, else open only the first
  for($fidx = 0; $fidx < $nfiles; $fidx++) { 
    open($out_FH_A[$fidx], ">", $out_filename_A[$fidx]) || ofile_FileOpenFailure($out_filename_A[$fidx], $sub_name, $!, "writing", $FH_HR);
    $isopen_A[$fidx] = 1;
  }
  $nopen = $nfiles; # will be decremented as we close files

  # read file until we see the first header line
  my ($next_header_line, $next_seqname) = sqf_FastaFileReadAndOutputNextSeq($in_FH, undef, $FH_HR); 
  # this will die if any non-whitespace characters exist before first header line

  while($nseq_remaining > 0) { 
    if(! defined $next_header_line) { 
      ofile_FAIL("ERROR in $sub_name, read too few sequences in $fa_file, read expected $tot_nseq", 1, $FH_HR); 
    }
    if(! exists $seqlen_HR->{$next_seqname}) { 
      ofile_FAIL("ERROR in $sub_name, no sequence length information exists for $next_seqname", 1, $FH_HR);
    }
    $nres_this_seq = $seqlen_HR->{$next_seqname};

    # first $nfiles sequences go to file $nseq_output so we guarantee we have >= 1 seq per file, 
    # remaining seqs are randomly placed
    if($nseq_output < $nfiles) { 
      $fidx = $nseq_output;
      $ridx = $f2r_map_A[$fidx];
    }
    else { 
      $ridx = int(rand($nopen)); 
      $fidx = $r2f_map_A[$ridx];
    }
    $FH = $out_FH_A[$fidx];

    # output seq
    print $FH $next_header_line; 
    ($next_header_line, $next_seqname) = sqf_FastaFileReadAndOutputNextSeq($in_FH, $FH, $FH_HR);
    
    $nseq_remaining--;
    $nseq_output++;
    
    # update counts of sequences and residues for the file we just printed to
    $nres_per_out_A[$fidx] += $nres_this_seq;
    $nseq_per_out_A[$fidx]++;
    $nres_tot_out += $nres_this_seq;

    # if we've reached our checkpoint output update
    if($do_debug && ($nres_tot_out > $checkpoint_nres)) { 
      my $nfiles_above_fract = 0;
      for(my $tmp_fidx = 0; $tmp_fidx < $nfiles; $tmp_fidx++) { 
        if($nres_per_out_A[$tmp_fidx] > ($checkpoint_fraction * $targ_nres)) { $nfiles_above_fract++; }
      }
      $checkpoint_fraction += $checkpoint_fraction_step;
      $checkpoint_nres = $checkpoint_fraction * $tot_nres;
    }

    # check if we need to close this file now, if so close it and open a new one (if nec)
    if(($nres_per_out_A[$fidx] >= $targ_nres) || ($nseq_remaining == 0)) { 
      if(($nopen > 1) || ($nseq_remaining == 0)) { 
        # don't close the final file unless we have zero sequences left
        close($out_FH_A[$fidx]);
        $isopen_A[$fidx] = 0;
        if($do_debug) { printf("$out_filename_A[$fidx] finished (%d seqs, %d residues)\n", $nseq_per_out_A[$fidx], $nres_per_out_A[$fidx]); }
        # update r2f_map_A so we can no longer choose the file handle we just closed 
        if($ridx != ($nopen-1)) { # edge case
          $r2f_map_A[$ridx] = $r2f_map_A[($nopen-1)];
        }
        $f2r_map_A[$r2f_map_A[($nopen-1)]] = $ridx; 
        $r2f_map_A[($nopen-1)] = -1; # this random index is now invalid
        $f2r_map_A[$fidx]      = -1; # this file is now closed
        $nopen--;
      }
    }
  }

  # go through and close any files that are still open
  for($fidx = 0; $fidx < $nfiles; $fidx++) { 
    if($isopen_A[$fidx] == 1) { 
      # file still open, close it
      close($out_FH_A[$fidx]);
      if($do_debug) { printf("$out_filename_A[$fidx] finished (%d seqs, %d residues)\n", $nseq_per_out_A[$fidx], $nres_per_out_A[$fidx]); }
    }
  }

  return $nfiles;
}

#################################################################
# Subroutine: sqf_FastaFileReadAndOutputNextSeq
# Incept:     EPN, Fri Jul  6 11:04:52 2018
#
# Purpose:    Given an open input file handle for a fasta sequence
#             and an open output file handle, read the next sequence and
#             and output it to the output file, by outputting all lines
#             we read until the next header line. Then return the
#             header line we stopped reading on.
#             If $out_FH is undef, do not output, which allows this
#             function to be called to return the first header line
#             of the file, but if $out_FH is undef, then require
#             that all lines read before the header line are empty,
#             else die.
#
# Arguments:
#   $in_FH:       input file handle
#   $out_FH:      output file handle, can be undef
#   $FH_HR:       ref to hash of file handles, for printing errors if we die"
# 
# Returns:  2 values:
#           $next_header_line: next header line in the file, undef 
#                              if we do not read one before end of the file
#           $next_seqname:     sequence name on next header line
#
# Dies:     If $out_FH is undef and there are non-whitespace characters
#           prior to the first header line read.
#           If $next_header_line is defined and we can't parse it to get the
#           sequence name.
#################################################################
sub sqf_FastaFileReadAndOutputNextSeq { 
  my $sub_name = "sqf_FastaFileReadAndOutputNextSeq";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($in_FH, $out_FH, $FH_HR) = @_;

  my $line = undef;
  $line = <$in_FH>;
  while((defined $line) && ($line !~ m/^\>/)) { 
    # does this line have any nonwhitespace characters? 
    if(! defined $out_FH) { 
      chomp $line;
      if($line =~ m/\S/) { 
        ofile_FAIL("ERROR in $sub_name, read line with non-whitespace character when none were expected:\n$line", 1, $FH_HR);
      }
    }
    else { 
      print $out_FH $line;
    }
    $line = <$in_FH>; 
  }
  my $seqname = undef;
  if(defined $line) { 
    if($line =~ m/^\>(\S+)/) { 
      $seqname = $1;
    }
    else { 
      ofile_FAIL("ERROR in $sub_name, unable to parse sequence name from header line: $line", 1, $FH_HR);
    }
  }
  
  return ($line, $seqname);
}

#################################################################
# Subroutine: sqf_EslReformatRun()
# Incept:     EPN, Fri Mar 15 12:56:04 2019
#
# Synopsis: Use esl-reformat to convert a file from one format to another
#
# Arguments:
#  $esl_reformat: esl-reformat executable file
#  $opt_string:   string to options to supply to esl-reformat
#                 (may not contain --informat)
#  $in_file:      input file for esl-reformat
#  $out_file:     output file from esl-reformat
#  $in_format:    input format
#  $out_format:   output format
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $opt_string contains '--informat'
#             if esl-reformat fails
#################################################################
sub sqf_EslReformatRun { 
  my $sub_name = "sqf_EslReformatRun";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_reformat, $opt_string, $in_file, $out_file, $in_format, $out_format, $opt_HHR, $FH_HR) = @_;

  if(defined $opt_string) { 
    if($opt_string =~ /\s*\-\-informat\s*/) { 
      ofile_FAIL("ERROR in $sub_name, opt_string $opt_string contains '--informat'", 1, $FH_HR);
    }
  }
  
  if(! defined $opt_string) { 
    $opt_string = "";
  }

  my $cmd = $esl_reformat . " --informat $in_format $opt_string $out_format $in_file > $out_file";
  utl_RunCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  # remove a .ssi file for newly created file if it exists
  my $ssi_file = $out_file . ".ssi";
  if(-e $ssi_file) { unlink $ssi_file; }

  return;
}

#################################################################
# Subroutine: sqf_EslAlimergeListRun()
# Incept:     EPN, Sun Apr 19 07:07:53 2020
#
# Synopsis: Use esl-alimerge to merge alignment files listed in
#           a list file
#
# Arguments:
#  $esl_alimerge: esl-reformat executable file
#  $list_file:    input file for esl-alimerge --list
#  $opt_string:   string of options to supply to esl-alimerge
#                 (may not contain, '-o', '--list' or '--outformat')
#  $out_file:     output alignment file from esl-alimerge
#  $out_format:   output format
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $opt_string contains '--list' or '--outformat' or -o
#             if esl-alimerge fails
#################################################################
sub sqf_EslAlimergeListRun { 
  my $sub_name = "sqf_EslAlimergeListRun";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_alimerge, $list_file, $opt_string, $out_file, $out_format, $opt_HHR, $FH_HR) = @_;

  if(defined $opt_string) { 
    if($opt_string =~ /\s*\-\-list\s*/) { 
      ofile_FAIL("ERROR in $sub_name, opt_string $opt_string contains '--list'", 1, $FH_HR);
    }
    if($opt_string =~ /\s*\-\-outformat\s*/) { 
      ofile_FAIL("ERROR in $sub_name, opt_string $opt_string contains '--outformat'", 1, $FH_HR);
    }
    if($opt_string =~ /\s*\-o\s*/) { 
      ofile_FAIL("ERROR in $sub_name, opt_string $opt_string contains '-o'", 1, $FH_HR);
    }
  }
  
  if(! defined $opt_string) { 
    $opt_string = "";
  }
  my $cmd = $esl_alimerge . " --list --outformat $out_format $opt_string $list_file > $out_file";
  utl_RunCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}

####################################################################
# the next line is critical, a perl module must return a true value
return 1;
####################################################################
