#!/bin/bash

RETVAL=0;
CURRETVAL=0;

# we require 0 or 1 args, if 1 arg, it must be 'teamcity'
do_teamcity=0
if [ "$#" -ne 0 ]; then
    if [ "$#" -gt 1 ] || [ "$1" != "teamcity" ]; then 
        echo "Usage:"
        echo "$0"
        echo "OR"
        echo "$0 teamcity"
        exit 1;
    fi
    # if we get here, there's 1 arg and it's 'teamcity'
    do_teamcity=1;
fi

if [ -z "${VADRSEQUIPDIR}" ] ; then
    echo "VADRSEQUIPDIR environment variable is not set, set it to top-level vadr/ dir and rerun"
    exit 1
fi

for test in \
    01-iss2-catlist.t \
; do
    if [ "$do_teamcity" -eq 1 ]; then
        echo "##teamcity[testStarted name=\"$test\" captureStandardOutput='true']"
    fi

    prove -v $VADRSEQUIPDIR/t/$test;
    CURRETVAL=$?

    if [ "$do_teamcity" -eq 1 ]; then 
        if [ "$CURRETVAL" -ne 0 ]; then
            echo "##teamcity[testFailed name=\"$test\" message=\"v-test.pl failure\"]"
        fi
        echo "##teamcity[testFinished name=\"$test\"]"
    fi

    if [ "$CURRETVAL" -ne 0 ]; then
        RETVAL=1
    fi
done

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed"
   exit 0
else 
   echo "FAIL: at least one test failed"
   exit 1
fi
