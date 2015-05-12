#!/usr/bin/perl

use warnings;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $chr = $ARGV[2];


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

	# Get current line
	my $currentLine = $_;
	# Split on space delimiter / / for HD, on \t for MD
#	@fields = split(/ /, $currentLine);	
	@fields = split(/\t/, $currentLine);
	# Print to file if chromosome matches
	if($fields[1] eq $chr) {print OUTPUT_FILE "$currentLine";}

}	

close(INPUT_FILE);
close(OUTPUT_FILE); 


