use strict;
use warnings FATAL => 'all';
use Test::More tests => 15;

# max length of a file name is 255 bytes
# max length of a file name plus path is 4096 bytes
# This length matches the PATH_MAX that is supported by the operating system.

BEGIN {
    use_ok( 'sqp_opts'  ) || print "Bail out!\n";
    use_ok( 'sqp_utils' ) || print "Bail out!\n";
    use_ok( 'sqp_ofile' ) || print "Bail out!\n";
}
my %opt_HH = ();
my @opt_order_A = ();
# need to define options that are used by any subroutines we test
opt_Add("-v",     "boolean", 0, 1, undef, undef, "be verbose",      "be verbose",      \%opt_HH, \@opt_order_A);
opt_Add("--keep", "boolean", 0, 1, undef, undef, "leave tmp files", "leave tmp files", \%opt_HH, \@opt_order_A);
opt_SetByUser("-v",     0, \%opt_HH);
opt_SetByUser("--keep", 0, \%opt_HH);

# arrays that store info for each test [..$t..$ntests-1]
my @nfiles_A        = (); # [0..$t..$ntests-1] number of files for each test
my @nlines_A        = (); # [0..$t..$ntests-1] number of lines per file for each test
my @outfile_A       = (); # [0..$t..$ntests-1] final output file name for each test
my @desc_A          = (); # [0..$t..$ntests-1] description of each test

# variables for a single, current test
my @exp_filelines_A = (); # expected lines from concatenated file for current test
my @cur_filelines_A = (); # actual lines from concatenated files
my @cur_file_A      = (); # list of files to concatenate for current test
my $exp_nlines      = undef; # expected number of lines in output file
my $cur_nlines      = undef; # actual number of lines in output file
my $exp_str         = undef; # expected string created when all line from output file are concatenated
my $cur_str         = undef; # actual string created when all line from output file are concatenated

my @to_remove_A = (); # array of files to remove before exiting
my ($t, $i, $j); # counters
#######################################
# utl_ConcatenateListOfFiles() tests
# utl_FileLinesToArray() tests
# utl_ANewLineDelimString() tests
@nfiles_A  = ();
@nlines_A  = ();
@outfile_A = ();
@desc_A    = ();
push(@desc_A,     "1 files with 7+ char name of 100 lines");
push(@nfiles_A,   "1");
push(@nlines_A,   "100");
push(@outfile_A,  "f1l100");

push(@desc_A,     "10 files with 7+ char name of 100 lines");
push(@nfiles_A,   "10");
push(@nlines_A,   "100");
push(@outfile_A,  "f10l100");

push(@desc_A,     "100 files with 9+ char name of 1000 lines");
push(@nfiles_A,   "100");
push(@nlines_A,   "1000");
push(@outfile_A,  "f100l1000");

push(@desc_A,     "1000 files with 200+ char name of 3 lines");
push(@nfiles_A,   "1000");
push(@nlines_A,   "3");
push(@outfile_A,  "thisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200cha");

push(@desc_A,     "10000 files with 200+ char name of 3 lines");
push(@nfiles_A,   "10000");
push(@nlines_A,   "3");
push(@outfile_A,  "thisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200chathisfilenameis200cha");

push(@desc_A,     "5552 files with 7+ char name of 2 lines (tests skipping rarely skipped loop)");
push(@nfiles_A,   "5552");
push(@nlines_A,   "2");
push(@outfile_A,  "f5kl2");

my $ntests = scalar(@desc_A);
for($t = 0; $t < $ntests; $t++) { 
  @exp_filelines_A = ();
  @cur_filelines_A = ();
  @cur_file_A = ();
  for($i = 1; $i <= $nfiles_A[$t]; $i++) { 
    my $cur_outfile = $outfile_A[$t] . "." . $i;
    open(OUT, ">", $cur_outfile) || die "ERROR unable to open $cur_outfile for writing"; 
    for($j = 1; $j <= $nlines_A[$t]; $j++) { 
      print OUT $j . "\n";
      push(@exp_filelines_A, $j);
    }
    close(OUT);
    push(@cur_file_A, $cur_outfile);
  }
  # concatenate files into one file
  utl_ConcatenateListOfFiles(\@cur_file_A, $outfile_A[$t], "01-iss2-catlist.t", \%opt_HH, undef);
  # read concatenate file into an array
  utl_FileLinesToArray($outfile_A[$t], 0, \@cur_filelines_A, undef);
  
  # check number of lines matches expected
  $exp_nlines = scalar(@exp_filelines_A);
  $cur_nlines = scalar(@cur_filelines_A);
  is($cur_nlines, $exp_nlines, "utl_ConcatenateListOfFiles() correct line numbers: $desc_A[$t]");
  
  # check actual lines match expected
  $exp_str = utl_AToNewLineDelimString(\@exp_filelines_A);
  $cur_str = utl_AToNewLineDelimString(\@cur_filelines_A);
  is($exp_str, $cur_str, "utl_ConcatenateListOfFiles() correct output: $desc_A[$t]");
  
  push(@to_remove_A, @cur_file_A); # these will be removed by utl_ConcatenateListOfFiles() unless something goes wrong, in which case we still want to remove them, because it's probably a lot of files
  push(@to_remove_A, $outfile_A[$t]); 
}
foreach my $file (@to_remove_A) { 
  if(-e $file) { 
    unlink $file;
  }
}
