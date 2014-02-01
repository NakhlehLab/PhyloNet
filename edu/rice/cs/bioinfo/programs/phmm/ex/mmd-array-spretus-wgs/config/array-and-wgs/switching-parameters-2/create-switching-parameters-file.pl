#!/usr/bin/env perl

#my @choices;

#     "wt" => "2 4 2 4", # self-transitions unnecessary

#     "sft" => "1 1,2,3 2 1,2,3,5,6,7,8,9,10,11,12,13,14,15",
#     "stf" => "2 4 1 4,5,6,7,8,9,10,11,12,13,14,15",
#     "uft" => "1 4,5,6,7,8,9,10,11,12,13,14,15 2 4",
#     "utf" => "2 1,2,3,5,6,7,8,9,10,11,12,13,14,15 1 1,2,3",
#     "vft" => "1 4,5,6,7,8,9,10,11,12,13,14,15 2 1,2,3,5,6,7,8,9,10,11,12,13,14,15",
#     "vtf" => "2 1,2,3,5,6,7,8,9,10,11,12,13,14,15 1 4,5,6,7,8,9,10,11,12,13,14,15",
#,
#     "wft" => "1 1,2,3 2 4",
#     "wtf" => "2 4 1 1,2,3"

sub printUsage {
    die "Configuration line format: <parameter name> <initial weight> <optimize flag> <source parental-tree p1a> <source gene-genealogy g1a> <sink parental-tree p2a> <sink gene genealogy g2a> ...\n" .
	"where we take all possible choices from 4 subsets corresponding to a 4-tuple suffix.\n";
}

my @choices = 
    (
     "w 1.0 false 1 1,2,3 1 1,2,3          2 4 2 4",
     "v 1.0 false 1 4,5,6,7,8,9,10,11,12,13,14,15 1 4,5,6,7,8,9,10,11,12,13,14,15           2 1,2,3,5,6,7,8,9,10,11,12,13,14,15 2 1,2,3,5,6,7,8,9,10,11,12,13,14,15",
     "s 1.0 true 1 1,2,3 1 4,5,6,7,8,9,10,11,12,13,14,15           2 4 2 1,2,3,5,6,7,8,9,10,11,12,13,14,15",
     "u 1.0 true 1 4,5,6,7,8,9,10,11,12,13,14,15 1 1,2,3           2 1,2,3,5,6,7,8,9,10,11,12,13,14,15 2 4"
    );
foreach my $entry (@choices) {
    my @args = split(/\s+/, $entry);
    if (scalar(@args) < 7) {
	print STDERR "ERROR: config line doesn't have at least seven fields: $entry\n";
	printUsage;
    }
    my $parameterName = shift(@args);
    my $initialWeight = shift(@args);
    my $optimizeFlag = shift(@args);
    print "$parameterName $initialWeight $optimizeFlag";
    if (scalar(@args) % 4 != 0) {
	print STDERR "ERROR: config line suffix incorrectly formatted: $entry\n";
	printUsage;
    }
    while (scalar(@args) > 0) {
	my @paChoices = split(/,/, shift(@args));
	my @gaChoices = split(/,/, shift(@args));
	my @pbChoices = split(/,/, shift(@args));
	my @gbChoices = split(/,/, shift(@args));
	foreach my $pa (@paChoices) {
	    foreach my $ga (@gaChoices) {
		foreach my $pb (@pbChoices) {
		    foreach my $gb (@gbChoices) {
			my $sga = $ga;
			if ($ga < 10) {
			    $sga = "0" . $sga;
			}

			my $sgb = $gb;
			if ($gb < 10) {
			    $sgb = "0" . $sgb;
			}
		
			my $ha = "p$pa,g$sga";
			my $hb = "p$pb,g$sgb";
				
			print " $ha|$hb";
		    }
		}
	    }
	}
    }
    print "\n";
}

