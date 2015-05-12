#!/usr/bin/perl -w

# Reads in a VCF file and replaces missing genotypes with most frequent allele

use strict;
use warnings;
use List::Util qw(first);
use List::MoreUtils qw(uniq);

my $file =$ARGV[0];
my $chr_info =$ARGV[1];
my $win_size =$ARGV[2];
my $win_step =$ARGV[3];


# =========================================================================================
# 1. Configure bins
# =========================================================================================

# Read in chromosome info
open(my $fh, '<', $chr_info) or die "Could not open file '$chr_info' $!";
my @chrs = <$fh>;
close $fh;

# Initialize arrays
my @tmp_chr;
my @tmp_start;
my @tmp_stop;
my @win_chr;
my @win_start;
my @win_stop;

# Chromosome iteration
for my $chr_n (0 .. $#chrs){

        # Identify chromosome
        my @var_chr = split('\t', $chrs[$chr_n]);
        chomp(@var_chr);

        # Reset temporary arrays
        my @tmp_chr = ();
        my @tmp_start = ();
        my @tmp_stop = ();

        # Mark start of each window
        for (my $i = 1; $i <= $var_chr[1]; $i += $win_step) { push(@tmp_start, $i); push(@tmp_chr, $var_chr[0]); }

        # Mark end of each window
        foreach my $tmp_pos (@tmp_start) { push(@tmp_stop, $tmp_pos + $win_size); }
        if ($tmp_stop[-1] > $var_chr[1]) { $tmp_stop[-1] = $var_chr[1]; }

        # Merge results into @win_* arrays
        push(@win_chr, @tmp_chr);
        push(@win_start, @tmp_start);
        push(@win_stop, @tmp_stop);

}


# Concatenate win data into single matrix
my @win_matrix;
for my $i (0 .. $#win_chr)
{
  $win_matrix[$i] = [ $win_chr[$i], $win_start[$i], $win_stop[$i] ];
}





# =========================================================================================
# 2. Read in VCF file, record counts and calculate Hp in bins
# =========================================================================================

# Open a file to output allele counts to
my $file_counts ="$ARGV[0].counts";
open(my $fh_counts, '>', $file_counts) or die "Could not open file '$file_counts' $!";

# Open a file to output bin results to
my $file_bins ="$ARGV[0].hp";
open(my $fh_bins, '>', $file_bins) or die "Could not open file '$file_bins' $!";

# Clear arrays to store data in
@win_chr =();
my @win_pos;
my @win_ref;
my @win_alt;
my $Nsnps = 0;
my $index = 0;
my $sum_Maj =0;
my $sum_Min =0;
my $Hp = 0;

# Open VCF file to read allele counts from
chomp($file);
open (FILE, $file); 

while (my $data = <FILE>) {
	chomp $data;
	if (substr($data, 0, 1) eq "#") {
	} else {
		# Grab AC and AN values to calculate #REF an #ALT counts
		$_ = $data;
		my ($AA) = /AC\=(.+?)[^0-9]/; 
		my ($AN) = /AN\=(.+?)[^0-9]/;
		my $AR = ($AN-$AA);
                my @vcf = split(/\t/, $data);
		print $fh_counts "$vcf[0] \t $vcf[1] \t $AR \t $AA\n";

		# Check values before writing to bin
		if ( $vcf[0] == $win_matrix[$index][0] &&
			$vcf[1] >= $win_matrix[$index][1] &&
			$vcf[1] <= $win_matrix[$index][2] )
		{
	                # Push data to array
        	        push(@win_chr, $vcf[0]);
	                push(@win_pos, $vcf[1]);
	                push(@win_ref, $AR);
	                push(@win_alt, $AA);

		} else {

			# Calculate sum of major anbd minor alleles
			$sum_Maj = 0;
			$sum_Min = 0;
			for my $i ( reverse 0 .. $#win_pos)
			{
				if($win_ref[$i] >= $win_alt[$i]) { 
					$sum_Maj = $sum_Maj + $win_ref[$i];
					$sum_Min = $sum_Min + $win_alt[$i];
				} else {
                                      $sum_Maj = $sum_Maj + $win_alt[$i];
                                      $sum_Min = $sum_Min + $win_ref[$i];
				}
			}

			# Calculate and record pooled heterozgosity = (2*aMj*aMn) / (aMj+aMn)^2
			if ($sum_Maj >0) { $Hp = (2 * $sum_Maj * $sum_Min) / ($sum_Maj + $sum_Min)**2; }
	
			# Output results to file			
			print $fh_bins "$win_matrix[$index][0] \t $win_matrix[$index][1] \t $win_matrix[$index][2] \t $win_chr[0] \t $win_pos[0] \t $win_pos[-1] \t $#win_chr \t $sum_Maj \t $sum_Min \t $Hp \n";
			$index += 1;

			# Remove elements from array that are out of bin dimensions
			for my $i (reverse 0 .. $#win_pos){
				if ( $win_pos[$i] <= $win_matrix[$index][1]) {
                                        shift @win_chr;
                                        shift @win_pos;
                                        shift @win_ref;
                                        shift @win_alt;
				}
			}

			# If moving to different chromosome then clear the array
			if(  $vcf[0] != $win_matrix[$index][0] ) {
				@win_chr =();
				@win_pos =();
				@win_ref =();
				@win_alt =();
			}

                        # Push new data to array
                       	push(@win_chr, $vcf[0]);
                       	push(@win_pos, $vcf[1]);
                       	push(@win_ref, $AR);
                       	push(@win_alt, $AA);
		}




	}
}
close FILE;
close $fh_counts;
close $fh_bins;








exit;









