#!/usr/bin/env perl

use constant SUBSTR_LENGTH => 6000;

my $file = shift;
my $start = shift; # 24185
my $length = shift; # 200

my $ih;
open ($ih, $file);
# print <$ih>;
while (<$ih>) {
    my $line = $_;
    chomp($line);
    my @fields = split(/\s+/, $line);
    if (scalar(@fields) == 2) {
	# length($fields[1]) - SUBSTR_LENGTH
	# SUBSTR_LENGTH
	print ($fields[0] . " " . substr($fields[1], $start, $length) . "\n");
    } else {
	print "$line\n";
    }
}
close ($ih);
