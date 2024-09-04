#!/bin/bash

RETVAL=0;

for test in \
    01-iss2-catlist.t \
; do
    prove -v ./$test;
    if [ $? != 0 ]; then
        RETVAL=1;
    fi   
done
if [ $RETVAL = 0 ]; then
   echo "Success: all tests passed"
   exit 0
else 
   echo "FAIL: at least one test failed"
   exit 1
fi
