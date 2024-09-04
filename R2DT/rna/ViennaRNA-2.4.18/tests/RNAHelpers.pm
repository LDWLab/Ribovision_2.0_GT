package RNAHelpers;

use strict;
use Exporter;
use warnings;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = "2.4.18";
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(MsgChecking getDataDirPath);
%EXPORT_TAGS = ( VERSION => $VERSION,
                 Messages => [qw(&MsgChecking)],
                 Paths => [qw(&getDataDirPath)] );

sub MsgChecking
{
  my $msg = shift;

  printf("Checking %s...\n", $msg);
}

sub getDataDirPath
{
    return "../tests/data/";
}

1;
