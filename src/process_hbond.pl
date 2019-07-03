#!/usr/bin/env perl
## res_ignore is to specify ending phosphate groups

use warnings;
use strict;

if ($#ARGV < 3) {
	print "Need to provide pdb, dssr and whatif output files, as well as hydbond output file\n";
	print "Usage: $0	<pdb>    <dssr>   <whatif>    <hydbond_dat>    <res_ignore>\n";
	exit;
}

#my $pdb     = $ARGV[0];
#my $dssr    = $ARGV[1];
#my $whatif  = $ARGV[2];
#my $hydbond = $ARGV[3];

my $pdb     = shift (@ARGV);
my $dssr    = shift (@ARGV);
my $whatif  = shift (@ARGV);
my $hydbond = shift (@ARGV);
my @res_ig  = @ARGV;

my (@A, @U, @G, @C);
open PDB, $pdb || die "Unable to open file $pdb   $!\n";
my $res_old = 0;
while (<PDB>) {
    @_ = split;
    next if $_[5] == $res_old || $_[0] eq "TER";
    if    ($_[3] eq 'U') { push @U, $_[5]; }
    elsif ($_[3] eq 'A') { push @A, $_[5]; }
    elsif ($_[3] eq 'G') { push @G, $_[5]; }
    elsif ($_[3] eq 'C') { push @C, $_[5]; }
    $res_old = $_[5];
}

my ($N_basepair, $N_canonical, $count) = (0, 0, 0);
my (@canonical_type, @bp_index1, @bp_index2);

open DSSR, $dssr || die "Unable to open file $dssr    $!\n";
while (<DSSR>) {
    # Only read (N_basepair + 1) lines after the pattern
    if (/^List of/ && /base pairs$/) {
	    @_ = split;
        $N_basepair = $_[2];
        next;
    }

    if ($N_basepair > 0) { $count++; }
    else { next; }
    next if $count == 1;  # Skip the second line
    last if $count > $N_basepair + 1;
    
    @_ = split;
    next if ($_[4] ne "WC" && $_[4] ne "Wobble");

    $N_canonical++;
    push @canonical_type, "$_[3]";

    $_[1] = substr($_[1], 1);
    $_[2] = substr($_[2], 1);
    push @bp_index1, "$_[1]";
    push @bp_index2, "$_[2]";
}
close DSSR;

my (@Hbond, @atm_index1, @atm_index2);

open WHATIF, $whatif || die "Unable to open file $whatif   $!\n";
while (<WHATIF>) {
    @_ = split;

    ####### non-canonical hydbond
    next if checkWC ($_[2], $_[8], \@bp_index1, \@bp_index2);
    my $new_bond = 1;
    for (my $i = 0; $i <= $#Hbond; $i++) {
        if (($atm_index1[$i] == $_[2]) && ($atm_index2[$i] == $_[8])) {
            $Hbond[$i] .= "$_[4] $_[10] ";
            $new_bond = 0;
            last;
        } elsif (($atm_index1[$i] == $_[8]) && ($atm_index2[$i] == $_[2])) {
            $Hbond[$i] .= "$_[10] $_[4] ";
            $new_bond = 0;
            last;
        }
    }

    next if not $new_bond;
    push @atm_index1, "$_[2]";
    push @atm_index2, "$_[8]";
    push @Hbond, "$_[4] $_[10] ";
}
close WHATIF;

open OUTFILE, ">", $hydbond || die "Unable to write to file $hydbond   $!\n";

for (my $i = 0; $i < $N_canonical; $i++) {
    print OUTFILE "NAT", $canonical_type[$i], " ", $bp_index1[$i], " ", $bp_index2[$i];
    if ($canonical_type[$i] eq "A-U") {
        print OUTFILE " N6 O4 N1 N3 TER\n";
    } elsif ($canonical_type[$i] eq "U-A") {
        print OUTFILE " O4 N6 N3 N1 TER\n";
    } elsif ($canonical_type[$i] eq "G-C") {
        print OUTFILE " N1 N3 N2 O2 O6 N4 TER\n";
    } elsif ($canonical_type[$i] eq "C-G") {
        print OUTFILE " N3 N1 O2 N2 N4 O6 TER\n";
    } elsif ($canonical_type[$i] eq "G-U") {
        print OUTFILE " N1 O2 O6 N3 TER\n";
    } elsif ($canonical_type[$i] eq "U-G") {
        print OUTFILE " O2 N1 N3 O6 TER\n";
    }
}

for (my $i = 0; $i <= $#Hbond; $i++) {
    # Do not consider interactions between neighboring residues
    next if $atm_index1[$i] == $atm_index2[$i] || $atm_index1[$i] == $atm_index2[$i] + 1 || $atm_index1[$i] + 1 == $atm_index2[$i];
    print OUTFILE "NATIVE ", $atm_index1[$i], " ", $atm_index2[$i], " ", $Hbond[$i], "TER\n";
}

for (my $i = 0; $i <= $#A; $i++) {
    next if (grep {$_ eq $A[$i]} @res_ig);
    for (my $j = 0; $j <= $#U; $j++) {
        next if (grep {$_ eq $U[$j]} @res_ig);
        next if checkWC ($A[$i], $U[$j], \@bp_index1, \@bp_index2);
        next if abs ($A[$i] - $U[$j]) < 5;
        #next if $A[$i] + 1 == $U[$j] || $A[$i] == $U[$j] + 1;
        print OUTFILE "A-U ", $A[$i], " ", $U[$j], " N6 O4 N1 N3 TER\n";
    }
}

for (my $i = 0; $i <= $#G; $i++) {
    next if (grep {$_ eq $G[$i]} @res_ig);
    for (my $j = 0; $j <= $#C; $j++) {
        next if (grep {$_ eq $C[$j]} @res_ig);
        next if checkWC ($G[$i], $C[$j], \@bp_index1, \@bp_index2);
        next if abs ($G[$i] - $C[$j]) < 5;
        #next if $G[$i] + 1 == $C[$j] || $G[$i] == $C[$j] + 1;
        print OUTFILE "G-C ", $G[$i], " ", $C[$j], " N1 N3 N2 O2 O6 N4 TER\n";
    }
}

for (my $i = 0; $i <= $#G; $i++) {
    next if (grep {$_ eq $G[$i]} @res_ig);
    for (my $j = 0; $j <= $#U; $j++) {
        next if (grep {$_ eq $U[$j]} @res_ig);
        next if checkWC ($G[$i], $U[$j], \@bp_index1, \@bp_index2);
        next if abs ($G[$i] - $U[$j]) < 5;
        #next if $G[$i] + 1 == $U[$j] || $G[$i] == $U[$j] + 1;
        print OUTFILE "G-U ", $G[$i], " ", $U[$j], " N1 O2 O6 N3 TER\n";
    }
}

close OUTFILE;

################################################
sub checkWC {
    my $res1 = shift;
    my $res2 = shift;
    my @index1 = @{$_[0]};
    my @index2 = @{$_[1]};

    my $WC = 0;
    for (my $i = 0; $i <= $#index1; $i++) {
        if ((($index1[$i] == $res1) && ($index2[$i] == $res2)) || (($index1[$i] == $res2) && ($index2[$i] == $res1))) {
            $WC = 1;
            last;
        }
    }
    return $WC;
}
