# sequip 0.x release notes 

### sequip 0.10 release (Sep 2023): Minor bug-fix update
  * Fixes bug in sqp_seqfile.pm:sqf_FeatureTableParse() related to
    single nucleotide coordinate spans.

### sequip 0.09 release (Dec 2021): Minor bug-fix update
  * Replaces printf with print in ofile_FAIL because printf
    would format strings if they had '%' characters in them. 

### sequip 0.08 release (Feb 2021): Minor update
  * Adds ofile_GetDirPath()
  * Removes 'require' statements and adds explanatory comment

### sequip 0.07 release (Nov 2020): Minor bug-fix update
  * Update of test file t/do-prove-all-tests-from-vadr.sh
  * No changes to any perl module files relative to 0.06.

### sequip 0.06 release (July 2020): Minor bug-fix update
  * sqp_utils.pm:utl_ConcatenateListOfFiles() is now more robust. 
    Previous version could, in extreme cases, create 'cat' calls that
    exceeded number of allowed characters in a command. Documented
    as github issue#2. Tests added to the testsuite to cover this.

---

For more information, see the [git log for the develop
branch](https://github.com/nawrockie/vadr/commits/develop).

