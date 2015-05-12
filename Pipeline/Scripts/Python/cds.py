#!/usr/bin/python
# David Wragg, GenPhySE INRA Toulouse

import sqlite3, sys, getopt, csv, numpy, os
from numpy import genfromtxt

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
        curs.execute('CREATE TABLE output (id INTEGER PRIMARY KEY AUTOINCREMENT, gene TEX, chr TEXT, start INT, end INT, cds INT, missense INT, synonymous INT, freqmis REAL, freqsyn REAL, alpha REAL, snpmis TEXT, snpsyn TEXT) ')
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
		tmp_freqmis = float(0)
		tmp_snpmis = ",".join([str(x[1]) for x in tmp])
		if tmp_missense > 0 :
			tmp_freqmis = float(sum([float(x[21]) for x in tmp]) / float(tmp_missense))

		curs.execute('SELECT * FROM vep INNER JOIN tmp ON vep.snp=tmp.snp WHERE vep.consequence LIKE "%ynonymou%" AND vep.gene LIKE (?)', [row])
		tmp = curs.fetchall()
		tmp_synonymous = len(tmp)
		tmp_freqsyn = float(0)
		tmp_snpsyn = ",".join([str(x[1]) for x in tmp])
		if tmp_synonymous >0 :
			tmp_freqsyn = float(sum([float(x[21]) for x in tmp]) / float(tmp_synonymous))
        
		tmp_d = 0
		tmp_d = float( (tmp_freqmis - tmp_freqsyn) / tmp_cds )


		# Write to table
		curs.execute('INSERT INTO output (gene, chr, start, end, cds, missense, synonymous, freqmis, freqsyn, alpha, snpmis, snpsyn) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)', (tmp_gene, tmp_chr, tmp_start, tmp_end, tmp_cds, tmp_missense, tmp_synonymous, tmp_freqmis, tmp_freqsyn, tmp_d, tmp_snpmis, tmp_snpsyn))
		conn.commit()


# Output to file
printout("\nWriting OUTPUT to : %s\n\n" % (outfile+'.csv'), RED)
curs.execute('SELECT gene, chr, start, end, cds, missense, synonymous, freqmis, freqsyn, alpha, snpmis, snpsyn FROM output')
tmp = curs.fetchall()
with open((outfile+'.csv'), 'wb') as fh:
	writer = csv.writer(fh, delimiter='\t')
	writer.writerow(['GENE', 'CHR', 'START', 'END', 'CDS', 'N_MISSENSE', 'N_SYNON', 'FQ_MISSENSE', 'FQ_SYNON', 'ALPHA', 'SNPs_MISSENSE', 'SNPs_SYNON'])
	writer.writerows(tmp)


