## Here, we define several environment variables required by our unit tests
##

# include current directory into our Python path variable
PYTHONPATH='.:${PYTHONPATH}'

export PYTHONPATH

# include path to the built executables to check their functionality later on
PATH=../src/bin:${PATH}

export PATH

# the path to various data files required to run our tests
export DATADIR=./data

# current ViennaRNA Package version
export CURRENT_VERSION=2.4.18

# set DIFF variable
export DIFF=/usr/bin/diff

# set results directories
export RNAFOLD_RESULTSDIR=./RNAfold/results
export RNAALIFOLD_RESULTSDIR=./RNAalifold/results
export RNACOFOLD_RESULTSDIR=./RNAcofold/results

# misc/ directory
export MISC_DIR=../misc
