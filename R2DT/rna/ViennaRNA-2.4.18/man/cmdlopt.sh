#/bin/sh
# this is a wrapper for help2man
# it enables it to obtain help, detailed-help and version information
# directly from a given .ggo file.
#
# the first argument passed to this script must be the .ggo file which
# should be parsed, while the second argument may be one of the following
# options


case $2 in
  --help)           no --show-help --set-version=2.4.18 < $1.ggo | grep -v "in doubt" | grep -v "Comments should be sent";;
  --version)        echo 2.4.18;;
  --detailed-help)  no --show-detailed-help --set-version=2.4.18 < $1.ggo | grep -v "in doubt" | grep -v "Comments should be sent";;
esac
