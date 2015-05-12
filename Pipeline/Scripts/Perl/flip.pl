#!/usr/bin/perl

use warnings;
use List::MoreUtils qw/ uniq /;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

# Input file
open(INPUT_FILE, "<$infile") or die "Couldn't open $infile";

# Output file
open (OUTPUT_FILE, ">$outfile");

# Write header line to output file
my $junk=<INPUT_FILE>;
$junk =~ s/ /\t/g;
print OUTPUT_FILE "$junk";

# Set row counter
$count = 0;

# Read rest of file
while (<INPUT_FILE>) {

	$count ++;

	my $currentLine = $_;
	@fields = split(/ /, $currentLine);
	# fields[0] == ProbeID
	# fields[1] == alleleA
	# fields[2] == alleleB
	# fields[3] == genotype

	# Extract alleles
	@alleles = split(//, "@fields[4..$#fields]");
	# Convert from array to string
	$alleles = "@alleles";
	# Strip white space
	$alleles =~ tr/ //ds;
	# Strip hyphen
	$alleles =~ tr/-//ds;
	# Strip zero if there are any (shouldn't be!)
	$alleles =~ tr/0//ds;
	# Strip linebreaks
	$alleles =~ s/\R//g;
	# Create an array of $alleles
	@alleles = split(//, $alleles);

	# Count unique alleles
	@unique = uniq @alleles;
	# If there are >2 then the alleles mismatch the annotation	
	# Because of the way the array is indexed (starting from 0)
	# Then we look for more than 1 alleles (first will be 0)...
	if($#unique>1) {
		print "Row $count: Annotation mismatch: @unique \n";

		if("$alleles[0]" eq "A") {$currentLine =~ s/T/A/g;}
		if("$alleles[0]" eq "T") {$currentLine =~ s/A/T/g;}
		if("$alleles[0]" eq "C") {$currentLine =~ s/G/C/g;}
		if("$alleles[0]" eq "G") {$currentLine =~ s/C/G/g;}

		if("$alleles[1]" eq "A") {$currentLine =~ s/T/A/g;}
		if("$alleles[1]" eq "T") {$currentLine =~ s/A/T/g;}
		if("$alleles[1]" eq "C") {$currentLine =~ s/G/C/g;}
		if("$alleles[1]" eq "G") {$currentLine =~ s/C/G/g;}
	}
	$currentLine =~ s/ /\t/g;
	print OUTPUT_FILE "$currentLine";

	# Line to break loop after nominated row
	#if($count==10) {last;}
}	

close(INPUT_FILE);
close(OUTPUT_FILE); 
