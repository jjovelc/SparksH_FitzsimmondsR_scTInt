#!/usr/bin/perl

use warnings;
use strict;

my ($file1, $file2) = @ARGV;

open my $Bt_in, "<", $file1;
open my $Bt_out, ">", $file1 . ".uniq";
open my $Rt_in, "<", $file2;
my $outfile = $file1 . "_and_" . $file2;
open my $Rt_out, ">", $file2 . ".uniq";
open my $both_out, ">", $outfile;

my %Bt;
while (my $id = <$Bt_in>) {
    chomp($id);
    $Bt{$id} = 1;
}

my %Rt;
while (my $id = <$Rt_in>) {
    chomp($id);
    $Rt{$id} = 1;
}

foreach my $key (sort keys %Bt) {
    if (exists $Rt{$key}) {
        print $both_out $key, "\n";
    }
    else {
        print $Bt_out $key, "\n";
    }
}

foreach my $key (sort keys %Rt) {
    if (!exists $Bt{$key}) {
        print $Rt_out $key, "\n";
    }
}

close($Bt_in);
close($Rt_in);
close($Bt_out);
close($Rt_out);
close($both_out);
