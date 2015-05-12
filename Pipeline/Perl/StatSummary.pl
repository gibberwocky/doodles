#!/usr/bin/perl -w

use strict;
use warnings;
use List::Util qw(first);
use List::MoreUtils qw(uniq);

my $path = $ARGV[0];
my $samples = $ARGV[1];
my @SAMPLE_LIST = ();
my $sample = ();
my $file = ();
my $fh = ();
my @tmp = ();
my @data = ();
my $data = ();


# ==============================================================================
# Sample List
# ==============================================================================
chomp($samples);
open($fh, $samples);
	@SAMPLE_LIST = <$fh>;
close $fh;



# Header
print "ID\tTOTAL_READS\tMAPPED_IN_PAIR\tSINGLETONS\tDUPLICATION\tHQ_ALIGNED\tINS_MEAN\tINS_SD\tDEPTH_COV\tCALLABLE\tNO_COV\tLOW_COV\tEXCESS_COV\tMQ_POOR\tSNPS_GATK\tSNPS_PILEUP\tSNPS_PLATYPUS\tBAYSIC_GATK\tBAYSIC_PILEUP\tBAYSIC_PLATYPUS\tBAYSIC_FINAL\n";


foreach (@SAMPLE_LIST) {


	$sample = $_;
	$sample =~ s/[\n\r\s]+//g;

	print "$sample\t";

	# ==============================================================================
	# Flagstat
	# ==============================================================================
	$file = join "", $path, "/",  $sample, "/metrics/", $sample, ".flagstat";
	my $READS = 0;
	my $PAIRS = 0;
	my $SINGLES = 0;
	chomp($file);
	open ($fh, $file);
	        @data = <$fh>;
	        @tmp = split(/\s/, $data[0]);
	        $READS = $tmp[0];
		shift @data;
		shift @data;
		shift @data;
		shift @data;
		shift @data;
		shift @data;
		shift @data;
		shift @data;
		shift @data;
		@tmp = split(/\s/, $data[0]);
		$PAIRS = $tmp[0];
		shift @data;
	        @tmp = split(/\s/, $data[0]);
	        $SINGLES = $tmp[0];
	        printf "$READS\t %.3f \t %.3f \t", $PAIRS/$READS, $SINGLES/$READS;
	close $fh;






	# ==============================================================================
	# Duplication metrics
	# ==============================================================================
	$file = join "", $path, "/",  $sample, "/metrics/", $sample, "_dup.metrics";
	my $DUP = 0;
	chomp($file);
	open ($fh, $file);
	        @data = <$fh>;
	        shift @data;
	        shift @data;
	        shift @data;
	        shift @data;
	        shift @data;
	        shift @data;
	        shift @data;
	        @tmp = split('\t', $data[0]);
	        $DUP = $tmp[7];
	        printf "%.3f \t", $DUP;
	close $fh;
	
	
	
	
	
	
	# ==============================================================================
	# Alignment metrics
	# ==============================================================================
	$file = join "", $path, "/",  $sample, "/metrics/", $sample, "_aln.metrics";
	my $HQ_ALN = 0;
	chomp($file);
	open ($fh, $file);
	        @data = <$fh>;
	        shift @data;
	        shift @data;
	        shift @data;
	        shift @data;
	        shift @data;
	        shift @data;
	        shift @data;
		shift @data;
		shift @data;
	        @tmp = split('\t', $data[0]);
	        $HQ_ALN = $tmp[8] / $tmp[5];
	        printf "%.3f \t", $HQ_ALN;
	close $fh;
	
	
	
	
	# ==============================================================================
	# Insert metrics
	# ==============================================================================
	$file = join "", $path, "/",  $sample, "/metrics/", $sample, "_ins.metrics";
	my $INS_MEAN = 0;
	my $INS_SD = 0;
	chomp($file);
	open ($fh, $file);
	        @data = <$fh>;
	        shift @data;
	        shift @data;
		shift @data;
		shift @data;
		shift @data;
		shift @data;
		shift @data;
	        @tmp = split('\t', $data[0]);
	        $INS_MEAN = $tmp[4];
	        $INS_MEAN =~ s/[\n\r\s]+//g;
		$INS_SD = $tmp[5];
		$INS_SD =~ s/[\n\r\s]+//g;
	        printf "%.3f \t %.3f \t", $INS_MEAN, $INS_SD;
	close $fh;
	
	
	
	
	
	
	# ==============================================================================
	# GATK.sample_summary
	# ==============================================================================
	$file = join "", $path, "/",  $sample, "/metrics/", $sample, "_GATKcov.sample_summary";
	my $DCOV = 0;
	chomp($file);
	open ($fh, $file);
	        @data = <$fh>;
	        shift @data;
		shift @data;
		@tmp = split('\t', $data[0]);
		$DCOV = $tmp[2];
		$DCOV =~ s/[\n\r\s]+//g;
		printf "%.3f \t", $DCOV;
	close $fh;
	
	
	
	
	
	# ==============================================================================
	# CALLABLE
	# ==============================================================================
	my $REF_N = 0;
	my $CALLABLE = 0;
	my $NOCOV = 0;
	my $LOCOV = 0;
	my $EXCOV = 0;
	my $POORMQ = 0;
	$file = join "", $path, "/",  $sample, "/metrics/", $sample, "_callable.summary";
	chomp($file);
	open ($fh, $file);
		@data = <$fh>;
	
		shift @data;
		@tmp = split ('REF_N', $data[0]);
		$REF_N = $tmp[1];
		$REF_N =~ s/[\n\r\s]+//g;
	
	        shift @data;
	        @tmp = split ('CALLABLE', $data[0]);
	        $CALLABLE = $tmp[1];
	        $CALLABLE =~ s/[\n\r\s]+//g;
	
	        shift @data;
	        @tmp = split ('NO_COVERAGE', $data[0]);
	        $NOCOV = $tmp[1];
	        $NOCOV =~ s/[\n\r\s]+//g;
	
	        shift @data;
	        @tmp = split ('LOW_COVERAGE', $data[0]);
	        $LOCOV = $tmp[1];
	        $LOCOV =~ s/[\n\r\s]+//g;
	
	        shift @data;
	        @tmp = split ('EXCESSIVE_COVERAGE', $data[0]);
	        $EXCOV = $tmp[1];
	        $EXCOV =~ s/[\n\r\s]+//g;
	
	        shift @data;
	        @tmp = split ('POOR_MAPPING_QUALITY', $data[0]);
	        $POORMQ = $tmp[1];
	        $POORMQ =~ s/[\n\r\s]+//g;
	
	close $fh;
	
	my $BASES = $REF_N + $CALLABLE + $LOCOV + $NOCOV + $EXCOV + $POORMQ;
	$CALLABLE = $CALLABLE/$BASES;
	$NOCOV = $NOCOV/$BASES;
	$LOCOV = $LOCOV/$BASES;
	$EXCOV = $EXCOV/$BASES;
	$POORMQ = $POORMQ/$BASES;
	printf "%.3f \t", $CALLABLE;
	printf "%.3f \t", $NOCOV;
	printf "%.3f \t", $LOCOV;
	printf "%.3f \t", $EXCOV;
	printf "%.3f \t", $POORMQ;
	
	
	
	
	
	
	
	# ==============================================================================
	# GATK
	# ==============================================================================
	my $GATK = 0;
	$file = join "", $path, "/", $sample, "/vcfs/", $sample, "_GATK_UG_pass.stats";
	chomp($file);
	open (FILE, $file);
	while ($data = <FILE>) {
	        chomp $data;
	        my $tmp = grep (/'snp_count'/, $data);
	        if ($tmp == 1) {
	                my @snps = split(/ => /, $data);
	                $GATK = $snps[1];
	                $GATK =~ s/,//g;
			last;
	        }
	}
	print "$GATK\t";
	close FILE;
	
	# ==============================================================================
	# PILEUP
	# ==============================================================================
	my $PILEUP = 0;
	$file = join "", $path, "/", $sample, "/vcfs/", $sample, "_pileup_pass.stats";
	chomp($file);
	open (FILE, $file);
	while (my $data = <FILE>) {
	        chomp $data;
	        my $tmp = grep (/'snp_count'/, $data);
	        if ($tmp == 1) {
	                my @snps = split(/ => /, $data);
	                $PILEUP = $snps[1];
	                $PILEUP =~ s/,//g;
	                last;
	        }
	}
	print "$PILEUP\t";
	
	# ==============================================================================
	# PLATYPUS
	# ==============================================================================
	my $PLATYPUS = 0;
	$file = join "", $path, "/", $sample, "/vcfs/", $sample, "_platypus_pass.stats";
	chomp($file);
	open (FILE, $file);
	while (my $data = <FILE>) {
	        chomp $data;
	        my $tmp = grep (/'snp_count'/, $data);
	        if ($tmp == 1) {
	                my @snps = split(/ => /, $data);
	                $PLATYPUS = $snps[1];
	                $PLATYPUS =~ s/,//g;
	                last;
	        }
	}
	print "$PLATYPUS\t";
	




        # ==============================================================================
        # BAYSIC
        # ==============================================================================
        $file = join "", $path, "/", $sample, "/vcfs/", $sample, "_BAYSIC_vcf.stats";
	my $BAYSIC = 0;
	my $BAYSIC_GATK = 0;
	my $BAYSIC_PILEUP = 0;
	my $BAYSIC_PLATYPUS = 0;
	my $tmp = ();
        chomp($file);
        open (FILE, $file);
        while (my $data = <FILE>) {
                chomp $data;
		# GATK
                $tmp = grep (/'$sample'/, $data);
                if ($tmp == 1) {
			$data = <FILE>;
			chomp $data;
                       	my @snps = split(/ => /, $data);
                        $BAYSIC_GATK = $snps[1];
                        $BAYSIC_GATK =~ s/,//g;
                }
                # Pileup
                $tmp = grep (/'.*pileup.*'/, $data);
                if ($tmp == 1) {
                        $data = <FILE>;
                        chomp $data;
                        my @snps = split(/ => /, $data);
                        $BAYSIC_PILEUP = $snps[1];
                        $BAYSIC_PILEUP =~ s/,//g;
                }
                # Platypus
                $tmp = grep (/'.*platypus.*'/, $data);
                if ($tmp == 1) {
                        $data = <FILE>;
                        chomp $data;
                        my @snps = split(/ => /, $data);
                        $BAYSIC_PLATYPUS = $snps[1];
                        $BAYSIC_PLATYPUS =~ s/,//g;
                }
                # BAYSIC
                $tmp = grep (/'all'/, $data);
                if ($tmp == 1) {
			my $FLAG = 0;
			do {
	                        $data = <FILE>;
				chomp $data;
				$tmp = grep (/'snp_count'/, $data);
				if ($tmp==1) {$FLAG=1};
			} while ($FLAG==0);
                	my @snps = split(/ => /, $data);
	                $BAYSIC = $snps[1];
                        $BAYSIC =~ s/,//g;
                }
        }
#	print "$BAYSIC_GATK\t$BAYSIC_PILEUP\t$BAYSIC_PLATYPUS\t$BAYSIC\n";
	printf "%.3f \t %.3f \t %.3f \t $BAYSIC\n", $BAYSIC_GATK/$GATK, $BAYSIC_PILEUP/$PILEUP, $BAYSIC_PLATYPUS/$PLATYPUS;
        close FILE;




}



exit;




