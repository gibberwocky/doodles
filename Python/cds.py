#!/usr/bin/env python
# David Wragg, GenPhySE INRA Toulouse

import sqlite3, sys, getopt, csv, numpy, os
import scipy.stats as st
import matplotlib.pyplot as plt
from math import *
from numpy import genfromtxt
from scipy import mean
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
Rstats = importr('stats')

#import warnings
#warnings.simplefilter('ignore')

# Continuous STDOUT flush
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

#following from Python cookbook, #475186
def has_colours(stream):
    if not hasattr(stream, "isatty"):
        return False
    if not stream.isatty():
        return False # auto color only on TTYs
    try:
        import curses
        curses.setupterm()
        return curses.tigetnum("colors") > 2
    except:
        # guess false in case of error
        return False
has_colours = has_colours(sys.stdout)


def printout(text, colour=WHITE):
        if has_colours:
                seq = "\x1b[1;%dm" % (30+colour) + text + "\x1b[0m"
                sys.stdout.write(seq)
        else:
                sys.stdout.write(text)






try:
    myopts, args = getopt.getopt(sys.argv[1:],"e:g:v:o:")
except getopt.GetoptError as e:
    print (str(e))
    print("Usage: %s -e Ensembl VEP output -g Species GTF file -v VCF file -o OUTPUT name" % sys.argv[0])
    sys.exit(2)

for opt, val in myopts:
    if opt == '-e':
        vepfile=val
    elif opt == '-g':
        gtffile=val
    elif opt =='-v' :
        vcffile=val
    elif opt =='-o' :
	outfile=val



printout ("\n================================================\n", BLUE)
printout ("Frequency-weighted analysis of CDS polymorphisms\n", GREEN)
printout ("================================================\n\n", BLUE)

# Connect to database
f1 = open((outfile+'.db'), 'w')
printout ("Creating database : %s\n" % (outfile+'.db'), RED)
conn = sqlite3.connect(outfile+'.db')
curs = conn.cursor()

try:
        # Create OUTPUT table
        curs.executescript('DROP TABLE IF EXISTS output;')
        curs.execute('CREATE TABLE output (id INTEGER PRIMARY KEY AUTOINCREMENT, gene TEX, chr TEXT, start INT, end INT, cds INT, missense INT, synonymous INT, meannon REAL, meansyn REAL, maxnon REAL, maxsyn REAL, sumnon REAL, sumsyn REAL, snpmis TEXT, snpsyn TEXT, meanA REAL, maxA REAL, sumA REAL) ')
        conn.commit()

	# Populate VEP
	printout ("Importing VEP file : %s\n" % (vepfile), RED )
	curs.executescript('DROP TABLE IF EXISTS vep;')
	curs.execute('CREATE TABLE vep (id INTEGER PRIMARY KEY AUTOINCREMENT, snp TEXT, pos INTEGER, allele TEXT, gene TEXT, transcript TEXT, feature TEXT, consequence TEXT, U1 TEXT, U2 TEXT, U3 TEXT, U4 TEXT, U5 TEXT, U6 TEXT, strand TEXT) ')
	conn.commit()
	vep = genfromtxt(vepfile, delimiter='\t', comments="#", dtype='str')
	for row in vep:
		tmp = row
		curs.execute('INSERT INTO vep (SNP, Pos, Allele, Gene, Transcript, Feature, Consequence, U1, U2, U3, U4, U5, U6, Strand) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', row )
	curs.execute('DELETE FROM vep WHERE Gene = "-"')
	curs.execute('CREATE INDEX vep_snp ON vep(snp)')
	conn.commit()

	# Count missense variants
	curs.execute('SELECT COUNT (DISTINCT snp) FROM vep WHERE consequence LIKE "%issens%"')
	tmp = curs.fetchall()
	vep_mis = [int(x[0]) for x in tmp][0]

	# Count synonymous variants
        curs.execute('SELECT COUNT (DISTINCT snp) FROM vep WHERE consequence LIKE "%ynonymou%"')
        tmp = curs.fetchall()
        vep_syn	= [int(x[0]) for x in tmp][0]

	printout ("...%s missense variants\n" % (vep_mis), BLUE)
	printout ("...%s synonymous variants\n" % (vep_syn), BLUE)

	# Populate GTF
	printout ("Importing GTF file : %s\n" % (gtffile), RED )
	curs.executescript('DROP TABLE IF EXISTS gtf;')
	curs.execute('CREATE TABLE gtf (id INTEGER PRIMARY KEY AUTOINCREMENT, chr TEXT, type TEXT, feature TEXT, start INTEGER, end INTEGER, strand TEXT, gene TEXT) ')
	conn.commit()
	curs.executescript('CREATE INDEX gtf_id ON gtf(id);')
	gtf = genfromtxt(gtffile, delimiter='\t', dtype='str')
	for row in gtf:
		tmp = row
		curs.execute('INSERT INTO gtf (chr, type, feature, start, end, strand, gene) VALUES (?,?,?,?,?,?,?)', row )
	conn.commit()


	# Populate VCF
	printout ("Importing VCF file : %s\n" % (vcffile), RED )
	curs.executescript('DROP TABLE IF EXISTS vcf;')
	curs.execute('CREATE TABLE vcf (id INTEGER PRIMARY KEY AUTOINCREMENT, chr TEXT, pos INTEGER, snp TEXT, ref TEXT, alt TEXT, freq REAL) ')
	conn.commit()
	vcf = genfromtxt(vcffile, delimiter="\t", dtype='str')
	for row in vcf:
		tmp = row
		curs.execute('INSERT INTO vcf (chr, pos, snp, ref, alt, freq) VALUES (?,?,?,?,?,?)', row )
	curs.execute('CREATE INDEX vcf_snp ON vcf(snp)')
	conn.commit()

except sqlite3.Error, e:
        if conn:
                conn.rollback()
                print "Error %s:" % e.args[0]
                print tmp
                sys.exit(1)



# Get list of distinct chromosomes
printout("Identifying distinct chromosomes in VCF file", GREEN)
curs.execute('SELECT DISTINCT chr FROM vcf')
chr_list = curs.fetchall()

# Loop through each chromosome and analyse data
for chr in chr_list:

	# Recover CDS SNPs from VEP table
	printout ("\nExtracting data and performing analysis for chromosome : %s\n" % (chr), GREEN)
	curs.executescript('DROP TABLE IF EXISTS tmp;')
	curs.execute('CREATE TABLE tmp AS SELECT * FROM vcf WHERE vcf.chr LIKE (?)', chr)
	conn.commit()
	curs.executescript('CREATE INDEX tmp_snp ON tmp(snp);')
	curs.execute('SELECT * FROM vep INNER JOIN tmp ON vep.snp=tmp.snp WHERE vep.consequence LIKE "%issens%" OR vep.consequence LIKE "%ynonymou%" ' )
	cds_snps = curs.fetchall()
	printout ("...%s records identified\n" % (len(cds_snps)), BLUE)
	gene_unique = [str(x[4]) for x in cds_snps]
	printout ("...%s distinct genes\n" % (len(set(gene_unique))), BLUE)

	# Calculte CDS lengths
	idx = 0
	for row in set(gene_unique):
		idx = idx + 1
		printout ("\r...%0.0f%%" % (100 * (float(idx) / float(len(set(gene_unique))))), BLUE)
		curs.execute('SELECT gene, start, end FROM gtf WHERE gene LIKE (?) AND feature="CDS"', [row])
		tmp = curs.fetchall()
		curs.execute('SELECT * FROM gtf WHERE gene LIKE (?)', [row])
		tmp_all = curs.fetchall()

		tmp_gene = row
		tmp_chr = min([int(x[1]) for x in tmp_all])
		tmp_start = min([int(x[4]) for x in tmp_all])
		tmp_end = max([int(x[5]) for x in tmp_all]) 
		tmp_cds = int(sum([int(x[2])-int(x[1]) for x in tmp]) + 1)

		curs.execute('SELECT * FROM vep INNER JOIN tmp ON vep.snp=tmp.snp WHERE vep.consequence LIKE "%issens%" AND vep.gene LIKE (?)', [row] )
		tmp = curs.fetchall()
		tmp_missense = len(tmp)
		tmp_meanmis = float(0)
		tmp_maxmis = float(0)
		tmp_summis = float(0)
		tmp_snpmis = ",".join([str(x[1]) for x in tmp])
		if tmp_missense > 0 :
			tmp_meanmis = float( sum([float(x[21]) for x in tmp]) / float(tmp_missense) )
			tmp_maxmis = float( max([float(x[21]) for x in tmp]) )
			tmp_summis = float( sum([float(x[21]) for x in tmp]) )

		curs.execute('SELECT * FROM vep INNER JOIN tmp ON vep.snp=tmp.snp WHERE vep.consequence LIKE "%ynonymou%" AND vep.gene LIKE (?)', [row])
		tmp = curs.fetchall()
		tmp_synonymous = len(tmp)
		tmp_meansyn = float(0)
		tmp_maxsyn = float(0)
		tmp_sumsyn = float(0)
		tmp_snpsyn = ",".join([str(x[1]) for x in tmp])
		if tmp_synonymous >0 :
			tmp_meansyn = float(sum([float(x[21]) for x in tmp]) / float(tmp_synonymous))
                        tmp_maxsyn = float( max([float(x[21]) for x in tmp]) )
                        tmp_sumsyn = float( sum([float(x[21]) for x in tmp]) )

        

		# Calcualte alpha
		alpha_mean = 0
		alpha_max = 0
		alpha_sum = 0
		alpha_mean = (float(tmp_missense) / (float(tmp_missense)+float(tmp_synonymous)) * float(tmp_meanmis)) - (float(tmp_synonymous) / (float(tmp_missense)+float(tmp_synonymous)) * float(tmp_meansyn))
		alpha_max = (float(tmp_missense) / (float(tmp_missense)+float(tmp_synonymous)) * float(tmp_maxmis)) - (float(tmp_synonymous) / (float(tmp_missense)+float(tmp_synonymous)) * float(tmp_maxsyn)) 
		alpha_sum = (float(tmp_missense) / (float(tmp_missense)+float(tmp_synonymous)) * float(tmp_summis)) - (float(tmp_synonymous) / (float(tmp_missense)+float(tmp_synonymous)) * float(tmp_sumsyn))

		# Write to table
		curs.execute('INSERT INTO output (gene, chr, start, end, cds, missense, synonymous, meannon, meansyn, maxnon, maxsyn, sumnon, sumsyn, snpmis, snpsyn, meana, maxa, suma) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', (tmp_gene, tmp_chr, tmp_start, tmp_end, tmp_cds, tmp_missense, tmp_synonymous, tmp_meanmis, tmp_meansyn, tmp_maxmis, tmp_maxsyn, tmp_summis, tmp_sumsyn, tmp_snpmis, tmp_snpsyn, alpha_mean, alpha_max, alpha_sum))
		conn.commit()



# Fit to exponential distribution and calculate P values
#printout("\nFinding best fit to distribution", YELLOW)

#curs.execute('SELECT gene, chr, start, end, cds, missense, synonymous, freqmis, freqsyn, alpha, pval, pfdr, zalpha, snpmis, snpsyn FROM output')
#tmp = curs.fetchall()

# Fit exponential distribution and pull p-values
#alpha = [float(x[9]) for x in tmp]


# Distributions to assess
#cdfs = [
#    "norm",            #Normal (Gaussian)
#    "alpha",           #Alpha
#    "anglit",          #Anglit
#    "arcsine",         #Arcsine
#    "beta",            #Beta
#    "betaprime",       #Beta Prime
#    "bradford",        #Bradford
#    "burr",            #Burr
#    "cauchy",          #Cauchy
#    "chi",             #Chi
#    "chi2",            #Chi-squared
#    "cosine",          #Cosine
#    "dgamma",          #Double Gamma
#    "dweibull",        #Double Weibull
#    "erlang",          #Erlang
#    "expon",           #Exponential
#    "exponweib",       #Exponentiated Weibull
#    "exponpow",        #Exponential Power
#    "fatiguelife",     #Fatigue Life (Birnbaum-Sanders)
#    "foldcauchy",      #Folded Cauchy
#    "f",               #F (Snecdor F)
#    "fisk",            #Fisk
#    "foldnorm",        #Folded Normal
#    "frechet_r",       #Frechet Right Sided, Extreme Value Type II
#    "frechet_l",       #Frechet Left Sided, Weibull_max
#    "gamma",           #Gamma
#    "gausshyper",      #Gauss Hypergeometric
#    "genexpon",        #Generalized Exponential
#    "genextreme",      #Generalized Extreme Value
#    "gengamma",        #Generalized gamma
#    "genlogistic",     #Generalized Logistic
#    "genpareto",       #Generalized Pareto
#    "genhalflogistic", #Generalized Half Logistic
#    "gilbrat",         #Gilbrat
#    "gompertz",        #Gompertz (Truncated Gumbel)
#    "gumbel_l",        #Left Sided Gumbel, etc.
#    "gumbel_r",        #Right Sided Gumbel
#    "halfcauchy",      #Half Cauchy
#    "halflogistic",    #Half Logistic
#    "halfnorm",        #Half Normal
#    "hypsecant",       #Hyperbolic Secant
#    "invgamma",        #Inverse Gamma
#    "invweibull",      #Inverse Weibull
#    "johnsonsb",       #Johnson SB
#    "johnsonsu",       #Johnson SU
#    "laplace",         #Laplace
#    "logistic",        #Logistic
#    "loggamma",        #Log-Gamma
#    "loglaplace",      #Log-Laplace (Log Double Exponential)
#    "lognorm",         #Log-Normal
#    "lomax",           #Lomax (Pareto of the second kind)
#    "maxwell",         #Maxwell
#    "mielke",          #Mielke's Beta-Kappa
#    "nakagami",        #Nakagami
#    "ncx2",            #Non-central chi-squared
#    "ncf",             #Non-central F
#    "nct",             #Non-central Student's T
#    "pareto",          #Pareto
#    "powerlaw",        #Power-function
#    "powerlognorm",    #Power log normal
#    "powernorm",       #Power normal
#    "rdist",           #R distribution
#    "reciprocal",      #Reciprocal
#    "rayleigh",        #Rayleigh
#    "rice",            #Rice
#    "recipinvgauss",   #Reciprocal Inverse Gaussian
#    "semicircular",    #Semicircular
#    "t"]               #Student's T
#    "triang",          #Triangular
#    "truncexpon",      #Truncated Exponential
#    "truncnorm",       #Truncated Normal
#    "tukeylambda",     #Tukey-Lambda
#    "uniform",         #Uniform
#    "vonmises",        #Von-Mises (Circular)
#    "wald",            #Wald
#    "weibull_min",     #Minimum Weibull (see Frechet)
#    "weibull_max",     #Maximum Weibull (see Frechet)
#    "wrapcauchy",      #Wrapped Cauchy
#    "ksone",           #Kolmogorov-Smirnov one-sided (no stats)
#    "kstwobign"]       #Kolmogorov-Smirnov two-sided test for Large N


#result = []
#for cdf in cdfs:
#    # Attempt to fit data to probability distributions
#    parameters = eval("st."+cdf+".fit(alpha)");
#    # Record maximum likelihood estimate
#    mle = eval("st."+cdf+".nnlf")(parameters, alpha)
#    tmp_nom = cdf.ljust(16)
#    printout("\n...%s\tMLE: %s" % (tmp_nom, mle), BLUE)
#    result.append([cdf, mle])

#best = sorted(result, key=lambda elem : mean(elem[1]), reverse=False)
#printout("\nBest MLE fit is %s, MLE: %s" % (best[0][0], best[0][1]), YELLOW)

# plot data (requires X server, cannot run on cluster)
#plt.hist(alpha, normed=True, bins=10)
#params = eval("st."+best[0][0]+".fit(alpha)")
#fct = eval("st."+best[0][0]+".freeze"+str(params))
#x = numpy.linspace(fct.ppf(0.0), max(alpha), 500)
#plt.plot(x, fct.pdf(x), lw=3, label=best[0][0])
#plt.axvline(fct.ppf(1-0.05), color='b', linestyle='dashed', lw=2)
#plt.legend(loc='best', frameon=False)
#plt.title("(alpha) Fit to distribution")
#plt.savefig(outfile+'-fit.pdf', format='pdf')

# Calculate P value, FDR adjusted P value, and Z(alpha)
#p_alpha = [1-fct.cdf(x) for x in alpha]
#p_fdr = list(Rstats.p_adjust(FloatVector(p_alpha), method='BH'))
#z_alpha = st.mstats.zscore(alpha)


# Update SQL table
#x = 0
#for row in tmp:
#	curs.execute('UPDATE output SET pval = (?), pfdr = (?), zalpha = (?) WHERE gene = (?)', (float(p_alpha[x]), float(p_fdr[x]), float(z_alpha[x]), row[0]))
#	x = x + 1
#	conn.commit()

# Write output to file
printout("\nWriting ALPHA OUTPUT to : %s\n\n" % (outfile+'-alpha.csv'), RED)
curs.execute('SELECT gene, chr, start, end, cds, missense, synonymous, meannon, meansyn, maxnon, maxsyn, sumnon, sumsyn, meana, maxa, suma FROM output')
tmp = curs.fetchall()
with open((outfile+'-alpha.csv'), 'wb') as fh:
	writer = csv.writer(fh, delimiter='\t')
	writer.writerow(['GENE', 'CHR', 'START', 'END', 'CDS', 'NON-SYN', 'SYN', 'Q_mean_N', 'Q_mean_S', 'Q_max_N', 'Q_max_S', 'Q_sum_N', 'Q_sum_S', 'A_mean', 'A_max', 'A_sum'])
	writer.writerows(tmp)

printout("\nWriting SNPs OUTPUT to : %s\n\n" % (outfile+'-snps.csv'), RED)
curs.execute('SELECT gene, snpmis, snpsyn FROM output')
tmp = curs.fetchall()
with open((outfile+'-snps.csv'), 'wb') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['GENE', 'SNPs_NON-SYN', 'SNPs_SYNON'])
        writer.writerows(tmp)



