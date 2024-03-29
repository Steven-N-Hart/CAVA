#!env/bin/python3

import datetime
import gzip
import os
import sys
from optparse import OptionParser

import pysam


# Process data file
def processData(options, datafile, IDs):
    outfile = open(options.output, 'w')

    counter = 0
    found = 0
    chrom_proc = ''
    for line in datafile:

        if options.input is not None and found == len(IDs): break

        line = line.strip()
        if line.startswith('#'): continue

        cols = line.split('\t')

        chrom = cols[0]
        if chrom != chrom_proc:
            chrom_proc = chrom
            sys.stdout.write('\rProcessing SNPs on chr' + chrom_proc + ' ...  ')
            sys.stdout.flush()

        pos = cols[1]
        ID = cols[2]
        ref = cols[3]
        alts = cols[4].split(',')

        if ID in IDs: found += 1

        if not len(ref) == 1: continue
        if not all(len(x) == 1 for x in alts): continue

        if options.input is not None and ID not in IDs: continue

        info = cols[7]
        info = info[info.index('dbSNPBuildID=') + 13:]
        if ';' in info:
            build = info[:info.index(';')]
        else:
            build = info

        if int(build) > int(options.release): continue

        outfile.write('\t'.join([ID, chrom, pos, cols[4]]) + '\n')
        counter += 1

    sys.stdout.write('\rProcessing SNPs on all chrs ... OK')
    sys.stdout.flush()
    print('')
    datafile.close()
    outfile.close()
    return counter


# Use Tabix to index output file
def indexFile(options):
    sys.stdout.write('Compressing output file ... ')
    sys.stdout.flush()
    pysam.tabix_compress(options.output, options.output + '.gz', force=True)
    sys.stdout.write('OK\n')
    sys.stdout.write('Indexing output file ... ')
    sys.stdout.flush()
    pysam.tabix_index(options.output + '.gz', seq_col=1, start_col=2, end_col=2, meta_char='#', force=True)
    sys.stdout.write('OK\n')


# Read records from file as a list
def readIDs(inputfn):
    ret = set()
    for line in open(inputfn): ret.add(line.strip())
    return list(ret)


# Print out meta data
def printMetaData(datafile):
    ret = None
    for line in datafile:
        line = line.strip()
        if line.startswith('#'):
            if line.startswith('##dbSNP_BUILD_ID='):
                print('Input file build: ' + line[17:])
                ret = int(line[17:])
            if line.startswith('##reference='): print('Reference: ' + line[12:] + '\n')
        else:
            return ret


#######################################################################################################################

# Version
with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as version_file:
    version = version_file.read().strip()

descr = 'CAVA dbsnp_db v' + version
epilog = '\nExample usage: dbsnp_db.py -d 00-All.vcf.gz -s 138 -o out\n\n'.format(version)
OptionParser.format_epilog = lambda self, formatter: self.epilog
parser = OptionParser(usage='CAVA-{}/dbsnp_db <options>'.format(version), version=version, description=descr,
                      epilog=epilog)
parser.add_option('-i', "--in", default=None, dest='input', action='store', help="Input filename (list of dbSNP IDs)")
parser.add_option('-o', "--out", default='output', dest='output', action='store', help="Output filename prefix")
parser.add_option('-s', "--snp", default=None, dest='release', action='store', help="dbSNP release version.py")
parser.add_option('-d', "--data", default=None, dest='data', action='store', help="Data file (00-All.vcf.gz)")
(options, args) = parser.parse_args()

# Check options
if options.release is None:
    print('\nError: no dbSNP release specified')
    print('Please use option -s to specify dbSNP release version.py\n')
    quit()
if options.data is None:
    print('\nError: no 00-All.vcf.gz data file specified')
    print('Please use option -d to specify path to data file\n')
    quit()
if not os.path.isfile(options.data):
    print('\nError: 00-All.vcf.gz file (' + options.data + ') cannot be found.\n')
    quit()

# Print out version.py information
print("\n---------------------------------------------------------------------------------------")
print('CAVA ' + version + ' dbSNP database preparation tool (dbsnp_db) is now running.')
print('Started: ', datetime.datetime.now(), '\n')

# Open data fie for reading
datafile = gzip.open(options.data, 'rt')

# Print out meta data
build = printMetaData(datafile)
datafile.seek(0)
if int(options.release) > build:
    print('Error: requested release must be <=' + str(build) + '\n')
    quit()

# Print out info
print('Requested dbSNP release: ' + str(options.release))
IDs = []
if options.input is not None:
    IDs = readIDs(options.input)
    print('\nInput file contains ' + str(len(IDs)) + ' dbSNP IDs to be included in the database')
print('')

# Create compressed output file
N = processData(options, datafile, IDs)

# Index output file with Tabix
indexFile(options)

# Remove uncompressed output file
os.remove(options.output)

# Print out summary information
print('\nA total of ' + str(N) + ' SNPs have been retrieved\n')
print('---------------------')
print('Output files created:')
print('---------------------')
print(options.output + '.gz (SNP database)')
print(options.output + '.gz.tbi (index file)')
print('')
print('dbsnp_db.py successfully finished: ', datetime.datetime.now())
print("---------------------------------------------------------------------------------------\n")
