#!/usr/bin/env perl

use FindBin qw($Bin);
use lib "/home/kl23/code/global";
use setenv;
use framework;
# also make environment changes
setenv::setenv();

use Getopt::Std;
use Cwd;

# for sanity's sake
use strict;

# for hi res timing - from Li-San
use Time::HiRes qw(gettimeofday);

# for temp files
use File::Temp qw/ tempfile tempdir /;

my $file = shift;
my $start = shift; # 24185
my $length = shift; # 200
my $outfile = shift;

my $alignmentRef = framework::readAlignment($file);

foreach my $k (sort {$a cmp $b} keys(%$alignmentRef)) {
    $$alignmentRef{$k} = substr($$alignmentRef{$k}, $start, $length);
}

framework::filePrintAlignment($alignmentRef, $outfile);

