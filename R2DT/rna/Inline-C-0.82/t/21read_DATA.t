# This file checks that a bug in Inline::read_DATA() has been fixed.
# The bug existed up to and including Inline-0.52.

use strict; use warnings; use diagnostics;
use FindBin '$Bin';
use lib $Bin;
use TestInlineSetup;
use Inline Config => DIRECTORY => $TestInlineSetup::DIR;

use Inline C => 'DATA';

print "1..1\n";

my $foo = foo1() + foo2();

if($foo == 15) {print "ok 1\n"}
else {
  warn "\$foo: $foo\n";
  print "not ok 1\n";
}

  __DATA__

  __C__

#define __SYMBOL1__

#define __SYMBOL2__ 8

int foo1() {
  int ret = __SYMBOL2__
            - 1;
   return ret;
}

int foo2() {
  return __SYMBOL2__;
}

