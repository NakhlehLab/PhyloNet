#!/usr/bin/env perl


my %choices = 
    (
     "sf" => "p1,p1",
     "st" => "p2,p2",
     "sft" => "p1,p2",
     "stf" => "p2,p1"
    );
foreach my $key (keys(%choices)) {
    print "$key 0.01";
    my @pa = split(",", $choices{$key});
    for (my $i = 4; $i <= 15; $i++) {
	my $si = $i;
	if ($i < 10) {
	    $si = "0" . $si;
	}
	
	print " $pa[0],g01|$pa[1],g$si";
	print " $pa[0],g02|$pa[1],g$si";
	print " $pa[0],g03|$pa[1],g$si";
    }
    print "\n";
}

my %choices2 = 
    (
     "uf" => "p1,p1",
     "ut" => "p2,p2",
     "uft" => "p1,p2",
     "utf" => "p2,p1"
    );
foreach my $key (keys(%choices2)) {
    print "$key 0.01";
    my @pa = split(",", $choices2{$key});
    for (my $i = 4; $i <= 15; $i++) {
	my $si = $i;
	if ($i < 10) {
	    $si = "0" . $si;
	}
	
	print " $pa[0],g$si|$pa[1],g01";
	print " $pa[0],g$si|$pa[1],g02";
	print " $pa[0],g$si|$pa[1],g03";
    }
    print "\n";
}

my %choices3 = 
    (
     "vf" => "p1,p1",
     "vt" => "p2,p2",
     "vft" => "p1,p2",
     "vtf" => "p2,p1"
    );
foreach my $key (keys(%choices3)) {
    print "$key 0.01";
    my @pa = split(",", $choices3{$key});
    for (my $i = 4; $i <= 15; $i++) {
	for (my $j = 4; $j <= 15; $j++) {
	    if ($i != $j) {
		my $si = $i;
		my $sj = $j;
		if ($i < 10) {
		    $si = "0" . $si;
		}
		if ($sj < 10) {
		    $sj = "0" . $sj;
		}
		print " $pa[0],g$si|$pa[1],g$sj";
	    }
	}
    }
    print "\n";
}

my %choices4 = 
    (
     "wf" => "p1,p1",
     "wt" => "p2,p2",
     "wft" => "p1,p2",
     "wtf" => "p2,p1"
    );
foreach my $key (keys(%choices4)) {
    print "$key 0.01";
    my @pa = split(",", $choices4{$key});
    for (my $i = 1; $i <= 3; $i++) {
	for (my $j = 1; $j <= 3; $j++) {
	    if ($i != $j) {
		my $si = $i;
		my $sj = $j;
		if ($i < 10) {
		    $si = "0" . $si;
		}
		if ($sj < 10) {
		    $sj = "0" . $sj;
		}
		print " $pa[0],g$si|$pa[1],g$sj";
	    }
	}
    }
    print "\n";
}

my %choices5 = 
    (
     "nft" => "p1,p2",
     "ntf" => "p2,p1"
    );
foreach my $key (keys(%choices5)) {
    print "$key 0.01";
    my @pa = split(",", $choices5{$key});
    for (my $i = 1; $i <= 15; $i++) {
	my $si = $i;
	if ($i < 10) {
	    $si = "0" . $si;
	}
	
	print " $pa[0],g$si|$pa[1],g$si";
    }
    print "\n";
}

