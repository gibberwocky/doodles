#!/usr/bin/perl -w

use strict;
use warnings;

my $file1 =$ARGV[0];
my $file2 =$ARGV[1];

chomp($file2);
open (FILE, $file2);
my @arr_intervals = <FILE>;
my %hsh_intervals = map { $_ => 1 } @arr_intervals;
close FILE;

open BAM, "samtools view $file1 |";
while(<BAM>){
	next if(/^(\@)/);  
        s/\n//;  s/\r//;  ## removing new line
        my @sam = split(/\t+/);


	if (exists($hsh_intervals{$sam[0]})){
		print "$sam[0]\n";
	}	
}



exit;



