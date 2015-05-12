#!/usr/bin/perl -w

use strict;  
use warnings;  
use Getopt::Long;  

# Declare variables 
my @nucs;  
my @data;  
my $data="";
my $seq;
my $GC_counter;  
my $sequence;  

# Define defaults for variables  
my $chr="";
my $fasta="";  
my $winlen="";  
my $increment="";  

parseArgs();   # subroutine to call GetOptions  

# Open the fasta file, read it in, pass only sequences to @data  
open( FASTA, $fasta ) || die "Doh! Can't open $fasta: $!\n";  
while( my $seq = <FASTA> ){  
chomp $seq;  
   if(substr($seq, 0, 1) eq ">") {
	if ($data ne "") {
		sliding_window();
		$data="";
	}
	$chr = substr($seq, 1, 3);
    }else{  
	$data .= $seq; 
    }  
}           
close FASTA;  

sliding_window();

exit; 


###########################SUBROUTINES##################################

sub parseArgs{  
#Message to print if mandatory variables not declared  
my $usage ="\nUsage: $0 --fasta /path/to/*.fasta --winlen <integer> --increment <integer> [options]  

Mandatory options:  
--fasta       -  path to the first input file (fasta) in .fasta or .fsa format  
--winlen      -  specify size of the sliding window  
--increment   -  specify increment length for the sliding window  
\n";  

my $options = GetOptions  
    ( 
    'fasta=s{1,1}'      =>  \$fasta,  
    'winlen=i{1,1}'     =>  \$winlen,  
    'increment=i{1,1}'  =>  \$increment,  
    );  
    if ( $fasta eq "" ){ die "\n\nDoh!: fasta input file must be specified!\n\n$usage\n"};  
    if ( $winlen eq "" ){die "\n\nDoh!: length of sliding window must be specified\n\n$usage\n"};  
    if ( $increment eq "" ){die "\n\nDoh!: increment length must be specified\n\n$usage\n"};  
};  

#########################################################################   

sub sliding_window{  
    my $start = 0;  
    while( $start < length( $data ) ){;  
        my $sequence = substr( $data, $start, $winlen );  
        print $chr, "\t", $start + 1, "\t", $start + $increment, "\t";  
        getGC( $sequence );  
        $start = $start + $increment;  
    }  
}  

#########################################################################  

sub getGC{  
    my $sequence = $_[0];  
    my $GC_counter = 0;  
    my @nucs = split( //, $sequence );  
    for ( my $i = 0; $i < @nucs; $i++ ){  
        if( $nucs[$i] =~ /G|C/g ){  
            $GC_counter++;  
        }  
    }  
    printf ( "%.2f\n", $GC_counter / $winlen );
}  

#########################################################################  
