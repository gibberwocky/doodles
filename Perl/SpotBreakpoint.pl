#!/usr/bin/perl -w

use strict;
use warnings;

my $file =$ARGV[0];
chomp($file);
open (FILE, $file);

# Chop input file string to identify filename
my @path = split(/\//, $file);
# Chop filename up
my @tmp = split(/_/, $path[-1]);
# Identify chromosome
my $chr = $tmp[0];
# Identify interval
my @pos = split(/-/, $tmp[1]);
# Identify left position
my $left = $pos[0];
# Identify right position
@tmp = split(/.pslx/, $pos[1]);
my $right = $tmp[0];

chomp $chr;
chomp $left;
chomp $right;

print "Contig\tChr\tStart\tEnd\n";
#print "Interval\t$chr\t$left\t$right\n";

while (my $data = <FILE>) {
	next if $. < 6;

	chomp $data;
	my @datapoints = split(/\t/, $data);

#print "$chr:$datapoints[13]\t$left:$datapoints[15]\t$right:$datapoints[16]\n";

	if ($datapoints[13] eq $chr){
		if ( ($datapoints[15] >= $left && $datapoints[15] <= $right) ||
			($datapoints[16] >= $left && $datapoints[16] <= $right) ) {

			print "$datapoints[9]\t$datapoints[13]\t$datapoints[15]\t$datapoints[16]\n";

		}
	}
}


close FILE;

exit;








#perl /home/plxdw1/Dropbox/Perl/SpotBreakpoint.pl ${OUT}/${SAMPLE[${sample_id}]}/Intervals/${tmp}/${tmp}.pslx

