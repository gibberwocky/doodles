#!/usr/bin/perl -w
my $file =$ARGV[0];
chomp($file);
open (FILE, $file);

my $position=0;

while ($data = <FILE>) {
#	next if $. < 3;

	$position++;
	chomp $data;

	if (substr($data, 0, 1) eq "@") {
	} else {
	@datapoints = split(/\t/, $data);

	# Get strand, requires processing the bit FLAG
	# Bit 5 (16) = Strand of query (0=forward, 1=reverse)
#	print "$datapoints[3]\t";
#	$bin = dec2bin($datapoints[3]);
#	print "$bin\t";
#	if (substr($bin, length($bin)-5,1) == 0){
#		print "F\t";
#	} else {
#		print "R\t";
#	}

	# Get CGIAR string
	$cgiar = $datapoints[5];
#	print "$cgiar\n";

	# Get position of S and M
	$cgiarS = index($cgiar, 'S');
	$cgiarM = index($cgiar, 'M');

	if ($cgiarS == -1) {
#		print "\n";
	} elsif ($cgiarS > $cgiarM) {
		$tmp = substr($cgiar, $cgiarM+1, $cgiarS-$cgiarM-1);
#		print "$tmp\t";
		# Get Sequence
		if ($tmp >= 20) {
			$sequence = substr($datapoints[9],length($datapoints[9])-$tmp,$tmp);
			print ">$position-$cgiar\n";
			print "$sequence\n";
		}
	} else {
		$tmp = substr($cgiar, 0, $cgiarS);
#		print "$tmp\t";
		# Get Sequence
		if ($tmp >= 20) {
			$sequence = substr($datapoints[9],0,$tmp);
			print ">$position-$cgiar\n";
			print "$sequence\n";
		}
	} 
}

}

sub dec2bin { 
    my $str = unpack('B32', pack('N', shift)); 
    $str =~ s/^0+(?=\d)//; 
    return $str;
}


close FILE;

exit;








#perl /home/plxdw1/Dropbox/Perl/SeqSoftClip2.pl Z_79945647-79946302_map.eav.sam 


