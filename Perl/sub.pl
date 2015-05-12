#!/usr/bin/perl

use warnings;
use List::MoreUtils qw/ uniq /;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $find = $ARGV[2];
my $replace = $ARGV[3];

# Input file
open(INPUT_FILE, "<$infile") or die "Couldn't open $infile";

# Output file
open (OUTPUT_FILE, ">$outfile");

# Write header line to output file
my $junk=<INPUT_FILE>;
$junk =~ s/ /\t/g;
print OUTPUT_FILE "$junk";

# Read rest of file
while (<INPUT_FILE>) {

	my $currentLine = $_;
	@fields = split(/\t/, $currentLine);
	# Extract genotypes
	@alleles = "@fields[4..$#fields]";
	$alleles = "@alleles";
	# Substitute characters
	$alleles =~ s/$find/$replace/g;
	# Replaces spaces with tabs
	$alleles =~ s/ /\t/g;
	# Output to file
	print OUTPUT_FILE "@fields[0..3] \t $alleles";

}	

close(INPUT_FILE);
close(OUTPUT_FILE); 


