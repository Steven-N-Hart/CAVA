#!/usr/bin/env python3


# Collection of basic classes and functions
#######################################################################################################################


import gzip
import logging
import os
import time
import sys


#######################################################################################################################

# Class representing a single variant call
class Variant(object):
    # Constructor
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = Sequence(ref)
        self.alt = Sequence(alt)
        shiftpos, a, b = self.trimCommonStart(self.ref, self.alt)
        self.pos = self.pos + shiftpos
        shiftpos, x, y = self.trimCommonEnd(a, b)
        self.ref = Sequence(x)
        self.alt = Sequence(y)
        self.flags = []
        self.flagvalues = []

    # Getting basic information about variant
    def info(self):
        return self.getFullPosition() + '_' + self.ref + '>' + self.alt

    # Getting chromosome name and position combined
    def getFullPosition(self):
        return self.chrom + ':' + str(self.pos)

    # Checking if the variant is a base substitution
    def isSubstitution(self):
        return len(self.ref) == 1 and len(self.alt) == 1

    # Checking if the variant is an insertion
    def isInsertion(self):
        return len(self.ref) == 0 and len(self.alt) > 0

    # Checking if the variant is a deletion
    def isDeletion(self):
        return len(self.ref) > 0 and len(self.alt) == 0

    # Checking if the variant is a complex indel
    def isComplex(self):
        return len(self.ref) > 0 and len(self.alt) > 0 and not self.isSubstitution()

    # Checking if the variant is in-frame
    def isInFrame(self):
        return (len(self.alt) - len(self.ref)) % 3 == 0

    # Checking if the variant overlaps with a genomic region
    def overlap(self, start, end):
        if self.isInsertion(): return (self.pos - 1 >= start) and (self.pos <= end)
        int1_start = self.pos
        int1_end = self.pos + len(self.ref) - 1
        int2_start = start
        int2_end = end
        if int1_start == int2_start: return True
        if int1_start > int2_start:
            return int2_end >= int1_start
        else:
            return int1_end >= int2_start

    # Getting the value of a given annotation flag
    def getFlag(self, flag):
        return self.flagvalues[self.flags.index(flag)]

    # Adding annotation flag to variant
    def addFlag(self, flag, value):
        self.flags.append(flag)
        self.flagvalues.append(value)

    # Annotating variant
    def annotate(self, ensembl, dbsnp, reference, impactdir):
        self.annotateWithType()
        if not ensembl is None: self = ensembl.annotate(self, reference, impactdir)
        if not dbsnp is None: self = dbsnp.annotate(self)

    # Annotating variant with its type
    def annotateWithType(self):
        if self.isSubstitution(): self.addFlag('TYPE', 'Substitution')
        if self.isInsertion(): self.addFlag('TYPE', 'Insertion')
        if self.isDeletion(): self.addFlag('TYPE', 'Deletion')
        if self.isComplex(): self.addFlag('TYPE', 'Complex')

    # Aligning variant on the plus strand
    def alignOnPlusStrand(self, reference):
        seq1 = reference.getReference(self.chrom, self.pos, self.pos + len(self.ref) - 1 + 100)
        seq2 = self.alt + reference.getReference(self.chrom, self.pos + len(self.ref),
                                                 self.pos + len(self.ref) + 100 - 1)
        left, seq1, seq2 = self.rightAlign(seq1, seq2)
        if len(seq1) == 0 or len(seq2) == 0:
            left -= 1
            base = reference.getReference(self.chrom, self.pos + left, self.pos + left)
            seq1, seq2 = base + seq1, base + seq2
        ret = Variant(self.chrom, self.pos + left, seq1, seq2)
        ret.flags = self.flags
        ret.flagvalues = self.flagvalues
        return ret

    # Aligning variant on the minus strand
    def alignOnMinusStrand(self, reference):
        seq1 = reference.getReference(self.chrom, self.pos - 100, self.pos + len(self.ref) - 1)
        s = reference.getReference(self.chrom, self.pos - 100, self.pos - 1)
        seq2 = s + self.alt
        N = len(s)
        left, seq1, seq2 = self.leftAlign(seq1, seq2)
        if len(seq1) == 0 or len(seq2) == 0:
            left = left - 1
            base = reference.getReference(self.chrom, self.pos + left - N, self.pos + left - N)
            seq1, seq2 = base + seq1, base + seq2
        ret = Variant(self.chrom, self.pos + left - N, seq1, seq2)
        ret.flags = self.flags
        ret.flagvalues = self.flagvalues
        return ret

    # Right-aligning two sequences
    def rightAlign(self, seq1, seq2):
        left, seq1, seq2 = self.trimCommonStart(seq1, seq2)
        right, seq1, seq2 = self.trimCommonEnd(seq1, seq2)
        return left, seq1, seq2

    # Left-aligning two sequences
    def leftAlign(self, seq1, seq2):
        right, seq1, seq2 = self.trimCommonEnd(seq1, seq2)
        left, seq1, seq2 = self.trimCommonStart(seq1, seq2)
        return left, seq1, seq2

    # Trimming common starting subsequence of two sequences
    def trimCommonStart(self, s1, s2):
        counter = 0
        while True:
            if len(s1) == 0 or len(s2) == 0: return counter, s1, s2
            if s1[0] != s2[0]: return counter, s1, s2
            s1, s2 = s1[1:], s2[1:]
            counter += 1

    # Trimming common ending subsequence of two sequences
    def trimCommonEnd(self, s1, s2):
        counter = 0
        while True:
            if len(s1) == 0 or len(s2) == 0: return counter, s1, s2
            if s1[-1] != s2[-1]: return counter, s1, s2
            s1, s2 = s1[:-1], s2[:-1]
            counter += 1


#######################################################################################################################

# Class representing a single VCF record     
class Record(object):
    # Constructor
    def __init__(self, line, options, targetBED):

        self.targetBED = targetBED

        # Parsing VCF format
        if options.args['inputformat'].upper() == 'VCF':
            cols = line.strip().split("\t")
            self.chrom = cols[0]

            if self.chrom.startswith('chr'):
                self.chrom_chr_prefix = True
                self.chrom = self.chrom[3:]
            else:
                self.chrom_chr_prefix = False

            self.pos = int(cols[1])
            self.id = cols[2]
            self.ref = cols[3]
            alts = cols[4].split(",")
            self.alts = []
            self.qual = cols[5]
            self.filter = cols[6]
            if len(cols) > 7:
                self.info = cols[7]
            else:
                self.info = ''
            if len(cols) > 8:
                self.rest = cols[8:]
            else:
                self.rest = []

        # Parsing TXT format
        if options.args['inputformat'].upper() == 'TXT':
            cols = line.strip().split("\t")
            self.id = cols[0]
            self.chrom = cols[1]
            self.pos = int(cols[2])
            self.ref = cols[3]
            self.alts = []
            alts = cols[4].split(",")
            self.qual = ''
            self.filter = 'PASS'
            self.info = ''
            self.rest = []

        # Creating a Variant object for each variant call in the record
        self.variants = []
        for alt in alts:
            # Initializing each Variant object with different alt allele
            var = Variant(self.chrom, self.pos, self.ref, alt)
            # DO  NOT FILTER
            """
            if 'N' in self.ref or 'N' in alt:
                logging.info('Variant ignored as allele contains unknown base (\'N\'): ' + self.chrom + ':' + str(
                    self.pos) + ' ' + self.ref + '>' + alt)
                continue
           
            if alt == '.':
                logging.info("Variant ignored because it is monomorphic reference: " + self.chrom + ':' + str(
                    self.pos) + ' ' + self.ref + '>' + alt)
                continue

            if not alt.strip('ACGT') is '' or not self.ref.strip('ACGT') is '':
                logging.info("Variant ignored as format of alt allele is not supported: " + self.chrom + ':' + str(
                    self.pos) + ' ' + self.ref + '>' + alt)
                continue
            """
            # Filtering by variant type (i.e. substitution, indel, insertion, deletion, complex indel), if required
            if options.args['type'].upper() == 'SUBSTITUTION' and not var.isSubstitution():
                continue
            if options.args[
                'type'].upper() == 'INDEL' and not var.isInsertion() and not var.isDeletion() and not var.isComplex():
                continue
            if options.args['type'].upper() == 'INSERTION' and not var.isInsertion():
                continue
            if options.args['type'].upper() == 'DELETION' and not var.isDeletion():
                continue
            if options.args['type'].upper() == 'COMPLEX' and not var.isComplex():
                continue

            # Filtering by BED file, if required
            if not targetBED is None:
                goodchrom = var.chrom
                if not goodchrom in targetBED.contigs:
                    goodchrom = 'chr' + goodchrom
                    if not goodchrom in targetBED.contigs: continue
                if not var.isInsertion():
                    start = var.pos
                    end = var.pos + len(var.ref) - 1
                else:
                    start = var.pos - 1
                    end = var.pos
                foundstart = False
                for _ in self.targetBED.fetch(region=goodchrom + ':' + str(start) + '-' + str(start)): foundstart = True
                foundend = False
                for _ in self.targetBED.fetch(region=goodchrom + ':' + str(end) + '-' + str(end)): foundend = True
                if not (foundstart and foundend): continue

            # Adding Variant object to this record
            self.variants.append(var)
            self.alts.append(alt)

    # Annotating record
    def annotate(self, ensembl, dbsnp, reference, impactdir):
        # Annotating each variant in the record
        for variant in self.variants:
            variant.annotate(ensembl, dbsnp, reference, impactdir)

    # Writing record (a ref, multiple alts) to output file
    def output(self, outformat, outfile, options, genelist, transcriptlist, snplist, stdout):
        outvariants = []
        outalts = []
        # Values of Entries for Separate Alleles are separated by ','
        # Entries for Multiple transcripts (for 1 allele) are separated by ':'

        for i in range(len(self.variants)):
            variant = self.variants[i]
            alt = self.alts[i]

            isTRANSCRIPT = ('TRANSCRIPT' in variant.flags)
            if isTRANSCRIPT:
                annTRANSCRIPT = variant.flagvalues[variant.flags.index('TRANSCRIPT')]
            else:
                annTRANSCRIPT = ''
            isDBSNP = ('DBSNP' in variant.flags)
            if isDBSNP:
                annDBSNP = variant.flagvalues[variant.flags.index('DBSNP')]
            else:
                annDBSNP = ''

            if len(genelist) > 0 or len(transcriptlist) > 0:
                if annTRANSCRIPT == '': continue

            # Removing non-annotated variants, if required
            if not options.args['nonannot'] and not (
                    (isTRANSCRIPT and not annTRANSCRIPT == '') or (isDBSNP and not annDBSNP == '')):
                continue

            # Filtering by gene, transcript or snp list, if required
            if isDBSNP:
                if len(snplist) > 0 and annDBSNP not in snplist: continue

            outvariants.append(variant)
            outalts.append(alt)

        # Skipping record if all variants have been removed
        if len(outvariants) == 0: return

        # Writing output in VCF format
        if outformat.upper() == 'VCF':

            # Creating first part of the VCF record (up to FILTER field)
            chromstr = 'chr' + self.chrom if self.chrom_chr_prefix else self.chrom
            record = [chromstr, str(self.pos), self.id, self.ref, ",".join(outalts), self.qual, self.filter]

            # Preparing components of the String to be added to the INFO field
            flags = []
            flagvalues = []
            for variant in outvariants:
                for i in range(len(variant.flags)):
                    key = variant.flags[i]
                    value = variant.flagvalues[i]
                    if value == '': value = '.'
                    if key in flags:
                        flagvalues[flags.index(key)].append(value)
                    else:
                        flags.append(key)
                        flagvalues.append([value])

            # Creating String added to the INFO field, separating multiple alleles by ',' (within alt alleles, multiple
            # transcripts separated by ':')
            added = ''
            for i in range(len(flags)):
                key = flags[i]
                value = ','.join(flagvalues[i])
                if len(added) > 0: added += ';'
                if options.args['prefix']:
                    added += 'CAVA_' + key + '=' + value
                else:
                    added += key + '=' + value

            ########################################################################
            # Add HGVS
            # HS: 10/21/2020 Added HGVSP and multi-transcript support to HGVSc
            #                added prefix support

            if options.args['prefix']:
                HGVSC_key = 'CAVA_HGVSc='
                HGVSP_key = 'CAVA_HGVSp='
            else:
                HGVSC_key = 'HGVSc='
                HGVSP_key = 'HGVSp='

            HGVSC = ''
            HGVSP = ''
            if ('TRANSCRIPT' in flags) and ('GENE' in flags) and ('CSN' in flags):
                # List of 1 entry per variant (alt-allele) .. and whithin each multiple transcripts are separated by ','
                transcripts_list = flagvalues[flags.index('TRANSCRIPT')]
                gene_list = flagvalues[flags.index('GENE')]
                csn_list = flagvalues[flags.index('CSN')]
                csn_hgvs_list = flagvalues[flags.index('CSNHGVS')]
                if not (len(transcripts_list) == 0 or len(gene_list) != len(transcripts_list) or len(gene_list) != len(
                        csn_list)):
                    for ihg in range(0, len(transcripts_list)):  # Loop over alt-alleles
                        hgtranscripts = transcripts_list[ihg].split(":")
                        hggenes = gene_list[ihg].split(":")
                        hgcsns = csn_list[ihg].split(":")
                        hgcsns_hgvs = csn_hgvs_list[ihg].split(":")
                        if not (len(hgtranscripts) == len(hggenes) and len(hgtranscripts) == len(hgcsns)):
                            print("ERROR transcripts GENE and CSN not same length\n")
                            print("--------------------------------------------------------------------\n")
                            if options.args['logfile']:
                                logging.error("Bug: TRANSCRIPT GENE and CSN not same length\n")
                            sys.exit(-1)
                        else:
                            trHGVSC = ''
                            trHGVSP = ''
                            for itr in range(0, len(hgtranscripts)):  # Loop over transcripts within each alt-alleles
                                hgtranscript = hgtranscripts[itr]
                                hggene = hggenes[itr]
                                hgcsn = hgcsns[itr]
                                hgcsn_hgvs = hgcsns_hgvs[itr]
                                tHGVSC = hgtranscript
                                tHGVSC += '(' + hggene + '):'
                                try:
#                                    cdna, prot = hgcsn.split('_p.')
#                                    prot = prot.replace('extX', "extTer")
                                    cdna, prot = hgcsn_hgvs.split('_p.')
                                except ValueError:  # Example c.802-51_802-14del38, splice
                                    cdna = hgcsn_hgvs
                                    prot = '.'
                                tHGVSC += cdna
                                if tHGVSC == '.(.):.':
                                    tHGVSC = '.'
                                if prot == '.' or prot == '':
                                    tHGVSP = '.'
                                else:
                                    # Get Protein matching Transcript from options
                                    if hgtranscript in options.transcript2protein:
                                        protid = options.transcript2protein[hgtranscript]
                                        tHGVSP = protid + ':p.(' + prot + ')'
                                    else:
                                        # tHGVSP=hgtranscript+':p.('+prot+')'
                                        # If user wants HGVSP, they can always look in the CSN... but don't want to
                                        # provide incorrect HGVSP
                                        tHGVSP = '.'
                                        # Only report warning if transcript2protein mapping file was provided.
                                        if len(options.transcript2protein) > 0:
                                            print(
                                                "WARNING: transcript " + hgtranscript + " not in transcript2protein "
                                                                                        "file, HGVSp will be invalid\n")
                                            print(
                                                "--------------------------------------------------------------------"
                                                "\n")
                                            if options.args['logfile']:
                                                logging.info(
                                                    "WARNING: transcript " + hgtranscript + " not in transcript2protein"
                                                                                            " file\n")
                                if trHGVSC == '':
                                    trHGVSC = tHGVSC
                                    trHGVSP = tHGVSP
                                else:
                                    trHGVSC = trHGVSC + ':' + tHGVSC
                                    trHGVSP = trHGVSP + ':' + tHGVSP
                            if HGVSC == '':
                                HGVSC = trHGVSC
                                HGVSP = trHGVSP
                            else:
                                HGVSC = HGVSC + ',' + trHGVSC
                                HGVSP = HGVSP + ',' + trHGVSP

            if HGVSC == '':
                HGVSC = HGVSC_key + 'None'
                HGVSP = HGVSP_key + 'None'
            else:
                HGVSC = HGVSC_key + HGVSC
                HGVSP = HGVSP_key + HGVSP

            # Add multi-transcripts/multi-allele HGVS to output record
            if added == '':
                added = HGVSC + ';' + HGVSP
            else:
                added += ';' + HGVSC + ';' + HGVSP

            # Adding second part of the VCF record (starting from the INFO field)
            if self.info == '.' or self.info == '':
                record += [added]
                record += self.rest
            else:
                record += [self.info + ';' + added]
                record += self.rest

            # Writing record to the output file
            if stdout:
                print('\t'.join(record))
            else:
                outfile.write('\t'.join(record) + '\n')

        # Writing output in TSV format
        # 10/22/2020 added HGVS and HGVSp between core record and extras.

        if outformat.upper() == 'TSV':

            # Iterating through variants 
            #   Print one line for each alt-allele*Variant Combination
            c = 0
            flags = []
            flagvalues = []
            for variant in outvariants:
                for i in range(len(variant.flags)):
                    key = variant.flags[i]
                    value = variant.flagvalues[i]
                    if value == '': value = '.'
                    if key in flags:
                        flagvalues[flags.index(key)].append(value)
                    else:
                        flags.append(key)
                        flagvalues.append([value])

                # Standardize chromosome M notation if options specify
                if options.args['normalized_mitochondrial_chrom'] == 'MT':
                    if self.chrom == 'chrM':
                        self.chrom = 'chrMT'
                    if self.chrom == 'M':
                        self.chrom = 'MT'
                if options.args['normalized_mitochondrial_chrom'] == 'M':
                    if self.chrom == 'chrMT':
                        self.chrom = 'chrM'
                    if self.chrom == 'MT':
                        self.chrom = 'M'

                # Creating first part of the TSV record (up to FILTER field)
                # Common to all transcripts for that variant
                record = self.id + '\t' + self.chrom + '\t' + str(self.pos) + '\t' + self.ref + '\t' + outalts[
                    c] + '\t' + self.qual + '\t' + self.filter

                # Number of transcripts overlapping with the variant
                if 'TRANSCRIPT' in variant.flags:
                    transcripts_list = variant.flagvalues[variant.flags.index('TRANSCRIPT')].split(":")
                    N = len(transcripts_list)
                    if N > 0 and ('GENE' in variant.flags) and ('CSN' in variant.flags):
                        gene_list = variant.flagvalues[flags.index('GENE')].split(":")
                        csn_list = variant.flagvalues[flags.index('CSN')].split(":")
                    else:
                        if len(transcripts_list) > 0:
                            print("ERROR: transcript " + transcripts_list[0] + " has no matching GENE and/or CSN\n")
                            print("--------------------------------------------------------------------\n")
                            if options.args['logfile']:
                                logging.info(
                                    "ERROR: transcript " + transcripts_list[0] + " has no matching GENE and/or CSN\n")
                            quit()
                        gene_list = []
                        csn_list = []
                        transcripts_list = []
                else:  # Allow having No matching transcripts (for variant outside genes/transcripts)
                    N = 1
                    transcripts_list = []
                    gene_list = []
                    csn_list = []

                # Iterating through the transcripts
                for i in range(N):
                    # Creating second part of the TSV record - Transcript-specific per line.
                    rest = ''
                    for j in range(len(variant.flags)):
                        # First, do not split annotation that is not transcript specific
                        if not variant.flags[j] in ['TRANSCRIPT', 'GENE', 'GENEID', 'TRINFO', 'LOC', 'CSN', 'CLASS',
                                                    'SO',
                                                    'IMPACT', 'ALTANN', 'ALTCLASS', 'ALTSO', 'ALTFLAG', 'PROTPOS',
                                                    'PROTREF', 'PROTALT']:
                            value = variant.flagvalues[j]
                            if value == '': value = '.'
                            rest += '\t' + value
                            continue

                        values = variant.flagvalues[j].split(':')
                        value = values[i]
                        if value == '': value = '.'
                        rest += '\t' + value
                    if len(transcripts_list) == 0 or len(gene_list) != len(transcripts_list) or len(gene_list) != len(
                            csn_list) or len(transcripts_list) != N:
                        HGVSC = 'None'
                        HGVSP = 'None'
                    else:
                        hgtranscript = transcripts_list[i]
                        HGVSC = hgtranscript
                        HGVSC += '(' + gene_list[i] + '):'
                        try:
                            cdna, prot = csn_list[i].split('_p.')
                            prot = prot.replace('X', "Ter")
                        except ValueError:  # Example c.802-51_802-14del38, splice
                            cdna = csn_list[i]
                            prot = '.'
                        HGVSC += cdna
                        if HGVSC == '.(.):.' or HGVSC == '():':
                            HGVSC = '.'
                        if prot == '.' or prot == '':
                            HGVSP = '.'
                        else:
                            # Get Protein matching Transcript from options
                            if hgtranscript in options.transcript2protein:
                                protid = options.transcript2protein[hgtranscript]
                                HGVSP = protid + ':p.(' + prot + ')'
                            else:
                                #                                HGVSP=hgtranscript+':p.('+prot+')'
                                HGVSP = '.'  # If user wants HGVSP, they can always look in the CSN... but don't want to provide incorrect HGVSP
                                if len(
                                        options.transcript2protein) > 0:  # Only report warning if transcript2protein mapping file was provided.
                                    print(
                                        "WARNING: transcript " + hgtranscript + " not in transcript2protein file, HGVSp will be invalid\n")
                                    print("--------------------------------------------------------------------\n")
                                    if options.args['logfile']:
                                        logging.info(
                                            "WARNING: transcript " + hgtranscript + " not in transcript2protein file\n")
                    # Writing record to the output file
                    if stdout:
                        print(record + rest + "\t" + HGVSC + "\t" + HGVSP)
                    else:
                        outfile.write(record + rest + "\t" + HGVSC + "\t" + HGVSP + '\n')

                c += 1


#######################################################################################################################

# Class representing a single Ensembl transcript
class Transcript(object):
    # Constructor
    def __init__(self, line, reference):
        self.exons = []
        cols = line.split('\t')
        self.TRANSCRIPT = cols[0]
        self.geneSymbol = cols[1]
        self.geneID = cols[2]
        self.TRINFO = cols[3]
        self.chrom = cols[4]
        self.strand = int(cols[5])
        self.transcriptStart = int(cols[6])
        self.transcriptEnd = int(cols[7])
        self.codingStart = int(cols[8])
        self.codingStartGenomic = int(cols[9])
        self.codingEndGenomic = int(cols[10])
        # Initializing and adding exons
        for i in range(1, len(cols) - 11, 2):
            self.exons.append(Exon(int((i + 1) / 2), int(cols[10 + i]), int(cols[11 + i])))

    # self.reference=reference.getReference(self.chrom,self.transcriptStart,self.transcriptEnd)

    # Getting the full coding sequence of the transcript
    def getCodingSequence(self, reference, variant, exonseqs):
        ret = ''
        ret_exonseqs = []

        for i in range(len(self.exons)):
            exon = self.exons[i]

            if not variant is None:
                if exon.start < variant.pos <= exon.end:
                    if self.strand == 1:
                        ret += reference.getReference(self.chrom, exon.start + 1,
                                                      variant.pos - 1) + variant.alt + reference.getReference(
                            self.chrom, variant.pos + len(variant.ref), exon.end)
                    else:
                        temp = Sequence(reference.getReference(self.chrom, exon.start + 1,
                                                               variant.pos - 1) + variant.alt + reference.getReference(
                            self.chrom, variant.pos + len(variant.ref), exon.end))
                        ret += temp.reverseComplement()
                    continue

            if self.strand == 1:
                if exonseqs is None:
                    exonseq = reference.getReference(self.chrom, exon.start + 1, exon.end)
                    ret += exonseq
                    ret_exonseqs.append(exonseq)
                else:
                    ret += exonseqs[i]
            else:
                if exonseqs is None:
                    exonseq = reference.getReference(self.chrom, exon.start + 1, exon.end).reverseComplement()
                    ret += exonseq
                    ret_exonseqs.append(exonseq)
                else:
                    ret += exonseqs[i]

        ret = ret[self.codingStart - 1:]

        if exonseqs is None:
            return ret, ret_exonseqs
        else:
            return ret, exonseqs

    # Getting the translated protein sequence of the transcript with the variant included as well as a list of the (reference) exon sequence.
    def getProteinSequence(self, reference, variant, exonseqs, codon_usage):
        # Translating coding sequence
        codon_usage = codon_usage
        codingsequence, exonseqa = self.getCodingSequence(reference, variant, exonseqs)
        if (self.chrom == 'chrM' or
                self.chrom == 'chrMT' or
                self.chrom == 'M' or
                self.chrom == 'MT'):
            ret = Sequence(codingsequence).translate('4')
        else:
            ret = Sequence(codingsequence).translate(codon_usage)
        return ret, exonseqa

    # Checking if a given position is outside the region between the start and stop codon
    def isPositionOutsideCDS(self, pos):
        if self.strand == 1:
            return (pos < self.codingStartGenomic) or (pos > self.codingEndGenomic)
        else:
            return (pos > self.codingStartGenomic) or (pos < self.codingEndGenomic)

    # Checking if a given position is upstream the region between the start and stop codon
    def isPositionOutsideCDS_5prime(self, pos):
        return (self.strand == 1 and pos < self.codingStartGenomic) or (
                self.strand == -1 and pos > self.codingStartGenomic)

    # Checking if a given position is downstream the region between the start and stop codon
    def isPositionOutsideCDS_3prime(self, pos):
        return (self.strand == 1 and pos > self.codingEndGenomic) or (self.strand == -1 and pos < self.codingEndGenomic)

    # Checking if the given variant is outside of the translated region of the transcript
    def isOutsideTranslatedRegion(self, variant):
        if self.strand == 1:
            if variant.isInsertion():
                if variant.pos <= self.codingStartGenomic: return True
                if variant.pos - 1 >= self.codingEndGenomic: return True
                return False
            else:
                if variant.pos + len(variant.ref) - 1 < self.codingStartGenomic: return True
                if variant.pos > self.codingEndGenomic: return True
                return False
        else:
            if variant.isInsertion():
                if variant.pos <= self.codingEndGenomic: return True
                if variant.pos - 1 >= self.codingStartGenomic: return True
                return False
            else:
                if variant.pos + len(variant.ref) - 1 < self.codingEndGenomic: return True
                if variant.pos > self.codingStartGenomic: return True
                return False

    # Checking if the given variant is outside of the translated region of the transcript, +/- the first and last 3 bases of the coding sequence
    def isOutsideTranslatedRegionPlus3(self, variant):
        if self.strand == 1:
            if variant.isInsertion():
                if variant.pos <= self.codingStartGenomic + 3: return True
                if variant.pos - 1 >= self.codingEndGenomic - 3: return True
                return False
            else:
                if variant.pos + len(variant.ref) - 1 < self.codingStartGenomic + 3: return True
                if variant.pos > self.codingEndGenomic - 3: return True
                return False
        else:
            if variant.isInsertion():
                if variant.pos <= self.codingEndGenomic + 3: return True
                if variant.pos - 1 >= self.codingStartGenomic - 3: return True
                return False
            else:
                if variant.pos + len(variant.ref) - 1 < self.codingEndGenomic + 3: return True
                if variant.pos > self.codingStartGenomic - 3: return True
                return False

    # Checking if the given variant overlaps with splicing region
    def isInSplicingRegion(self, variant, ssrange):
        if not self.isOutsideTranslatedRegion(variant):
            for exon in self.exons:
                if variant.overlap(exon.end + 1, exon.end + ssrange): return True
                if variant.overlap(exon.start - (ssrange - 1), exon.start): return True
            return False
        else:
            return False

    # Checking if the given variant affects an essential splice site
    def isInEssentialSpliceSite(self, variant):
        if not self.isOutsideTranslatedRegion(variant):
            for exon in self.exons:
                if variant.overlap(exon.end + 1, exon.end + 2): return True
                if variant.overlap(exon.start - 1, exon.start): return True
            return False
        else:
            return False

    # Checking if the given variant affects a +5 essential splice site
    def isIn_SS5_Site(self, variant):
        if not self.isOutsideTranslatedRegion(variant):
            if self.strand == 1:
                for exon in self.exons:
                    if variant.overlap(exon.end + 1, exon.end + 5):
                        if not variant.isSubstitution():
                            if not (variant.pos == exon.end + 3 and len(variant.ref) == 2 and len(
                                    variant.alt) == 2): return True
                        else:
                            if variant.pos == exon.end + 5: return True
                return False
            else:
                for exon in self.exons:
                    if variant.overlap(exon.start - 4, exon.start):
                        if not variant.isSubstitution():
                            if not (variant.pos == exon.start - 3 and len(variant.ref) == 2 and len(
                                    variant.alt) == 2): return True
                        else:
                            if variant.pos == exon.start - 4: return True
                return False
        else:
            return False

    # Checking if the given variant affects the first or last 3 bases of an exon
    def isInFirstOrLast3BaseOfExon(self, variant):
        if not self.isOutsideTranslatedRegionPlus3(variant):
            for exon in self.exons:
                if variant.overlap(exon.start + 1, exon.start + 3): return True
                if variant.overlap(exon.end - 2, exon.end): return True
            return False
        else:
            return False

    # Checking where a given genomic position is located in the transcript
    def whereIsThisPosition(self, pos):
        # Iterating through exons and introns and checking if genomic position is located within

        for exon in self.exons:
            if exon.index > 1 and ((self.strand == 1 and prevexonend < pos <= exon.start) or (
                    self.strand == -1 and exon.end < pos <= prevexonend)):
                if self.intronLength(exon.index) > 5 or self.intronLength(exon.index) == 3:
                    return 'In' + str(exon.index - 1) + '/' + str(exon.index)
                else:
                    return 'fsIn' + str(exon.index - 1) + '/' + str(exon.index)
            if exon.start < pos <= exon.end:
                if (self.strand == 1 and pos < self.codingStartGenomic) or (
                        self.strand == -1 and pos > self.codingStartGenomic): return '5UTR'
                if (self.strand == 1 and pos > self.codingEndGenomic) or (
                        self.strand == -1 and pos < self.codingEndGenomic): return '3UTR'
                return 'Ex' + str(exon.index)
            prevexonend = exon.end if self.strand == 1 else exon.start

        return '.'

    # Checking where a given variant is located in the transcript
    def whereIsThisVariant(self, variant):
        # Getting the locations of both end points of the variant
        if variant.isInsertion():
            first = self.whereIsThisPosition(variant.pos - 1)
            second = self.whereIsThisPosition(variant.pos)
        else:
            first = self.whereIsThisPosition(variant.pos)
            second = self.whereIsThisPosition(variant.pos + len(variant.ref) - 1)
        if first == second:
            return first
        if self.strand == 1:
            return first + '-' + second
        else:
            return second + '-' + first

    # Getting the length of an intron, where idx is the index of the succeeding exon
    def intronLength(self, idx):
        for exon in self.exons:
            if exon.index == idx:
                if self.strand == 1:
                    return exon.start - prev
                else:
                    return prev - exon.end
            if self.strand == 1:
                prev = exon.end
            else:
                prev = exon.start


#######################################################################################################################

# Class representing a single exon
class Exon(object):
    # Constructor
    def __init__(self, index, start, end):
        self.index = index
        self.start = start
        self.end = end
        self.length = end - start

    def contains(self, pos):
        return (self.start + 1 <= pos <= self.end)


#######################################################################################################################

# Class representing a DNA sequence
class Sequence(str):
    # Translating to amino acid sequence
    def translate(self, letter):
        gencode = {}
        if letter == '1':
            gencode = {
                'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
                'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'}
        if letter == '3':
            gencode = {
                'ATA': 'Ile', 'ATC': 'Ile', 'ATT': 'Ile', 'ATG': 'Met',
                'ACA': 'Thr', 'ACC': 'Thr', 'ACG': 'Thr', 'ACT': 'Thr',
                'AAC': 'Asn', 'AAT': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
                'AGC': 'Ser', 'AGT': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
                'CTA': 'Leu', 'CTC': 'Leu', 'CTG': 'Leu', 'CTT': 'Leu',
                'CCA': 'Pro', 'CCC': 'Pro', 'CCG': 'Pro', 'CCT': 'Pro',
                'CAC': 'His', 'CAT': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
                'CGA': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGT': 'Arg',
                'GTA': 'Val', 'GTC': 'Val', 'GTG': 'Val', 'GTT': 'Val',
                'GCA': 'Ala', 'GCC': 'Ala', 'GCG': 'Ala', 'GCT': 'Ala',
                'GAC': 'Asp', 'GAT': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
                'GGA': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGT': 'Gly',
                'TCA': 'Ser', 'TCC': 'Ser', 'TCG': 'Ser', 'TCT': 'Ser',
                'TTC': 'Phe', 'TTT': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
                'TAC': 'Tyr', 'TAT': 'Tyr', 'TAA': 'X', 'TAG': 'X',
                'TGC': 'Cys', 'TGT': 'Cys', 'TGA': 'X', 'TGG': 'Trp'}
        if letter == '4':
            gencode = {
                'ATA': 'M', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGC': 'S', 'AGT': 'S', 'AGA': 'X', 'AGG': 'X',
                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
                'TGC': 'C', 'TGT': 'C', 'TGA': 'W', 'TGG': 'W'}
        if letter == '5':
            gencode = {
                'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
                'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
                'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
                'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
                'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
                'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                'TAA': 'X', 'TAC': 'Y', 'TAG': 'X', 'TAT': 'Y',
                'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                'TGA': 'X', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
                'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}
        ret = ''
        index = 0
        while index + 3 <= len(self):
            codon = self[index:index + 3].upper()
            if codon.replace('A', '').replace('C', '').replace('G', '').replace('T', ''):
                # if 'N' in codon:
                ret += '?'
                index += 3
                continue
            ret += gencode[codon]
            index += 3
        return ret

    # Getting reverse complement sequence
    def reverseComplement(self):
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "a": "t", "t": "a", "c": "g", "g": "c",
                      "n": "n", "*": "*"}
        ret = ''
        for base in self[::-1]: ret += complement[base]
        return ret


#######################################################################################################################

# Class representing the collection of options specified in the configuration file
class Options(object):
    # Constructor
    def __init__(self, configfn):

        self.defs = dict()
        self.configfn = configfn
        self.transcript2protein = dict()

        # Defining option flags
        self.defs['reference'] = ('string', '.')
        self.defs['ensembl'] = ('string', '.')
        self.defs['dbsnp'] = ('string', '.')
        self.defs['filter'] = ('boolean', False)
        self.defs['inputformat'] = ('string', 'VCF')
        self.defs['outputformat'] = ('string', 'VCF')
        self.defs['type'] = ('string', 'ALL')
        self.defs['logfile'] = ('boolean', False)
        self.defs['target'] = ('string', '.')
        self.defs['genelist'] = ('string', '.')
        self.defs['transcriptlist'] = ('string', '.')
        self.defs['snplist'] = ('string', '.')
        self.defs['nonannot'] = ('boolean', True)
        self.defs['givealt'] = ('boolean', True)
        self.defs['givealtflag'] = ('boolean', False)
        self.defs['ssrange'] = ('string', '8')
        self.defs['ontology'] = ('string', 'both')
        self.defs['impactdef'] = ('string', 'SG,ESS,FS|SS5,IM,SL,EE,IF,NSY|SY,SS,INT,5PU,3PU')
        self.defs['prefix'] = ('boolean', False)
        self.defs['codon_usage'] = ('string', '1')
        self.defs['normalized_mitochondrial_chrom'] = ('string', 'not_normalized')

        self.defs['transcript2protein'] = ('string', '.')

        # Reading options from file
        self.read()

        if (self.args['ensembl'] == '.' or self.args['ensembl'] == ''):
            d = os.path.dirname(os.path.realpath(__file__))
            dir = d[:d.rfind('env/lib')]
            self.args['ensembl'] = dir + '/defaultdb/ensembl75s.gz'

    # Reading options from configuration file
    def read(self):
        self.args = dict()
        for line in open(self.configfn):
            line = line.strip()
            if line.startswith('@'):
                key = line[1:line.index('=')].strip()
                if key in list(self.defs.keys()):
                    (typeofvar, default) = self.defs[key]
                    if typeofvar == 'string': self.args[key] = line[line.find('=') + 1:].strip()
                    if typeofvar == 'list': self.args[key] = line[line.find('=') + 1:].strip().split(',')
                    if typeofvar == 'boolean': self.args[key] = (line[line.find('=') + 1:].strip().upper() == 'TRUE')
        for key, (typeofvar, default) in self.defs.items():
            if not key in list(self.args.keys()): self.args[key] = default


#######################################################################################################################            

# Other basic utility functions 

# Reading gene, transcript or snp list from file
def readSet(options, tag):
    ret = set()
    if tag in list(options.args.keys()) and not (options.args[tag] == '' or options.args[tag] == '.'):
        for line in open(options.args[tag]):
            line = line.strip()
            if line == '' or line == '.': continue
            ret.add(line)
        if options.args['logfile']:
            txt = ''
            if tag == 'genelist': txt = 'Gene list'
            if tag == 'transcriptlist': txt = 'Transcript list'
            if tag == 'snplist': txt = 'SNP list'
            logging.info(txt + ' loaded.')
    return ret


# Reading dictionary (from transcript -> Protein)
def read_dict(options, tag):
    ret = dict()
    tx_to_prot_source_file = options.args['ensembl'].replace('db.gz', 'txt').replace('.gz','.txt')
    print("Reading option file dictionary : " + tx_to_prot_source_file + "\n")
    with open(tx_to_prot_source_file, 'r') as f:
        for line in f:
            if line == '' or line == '.' or line.startswith("#"): continue
            linedat = line.strip().split("\t")
            if len(linedat) == 3:
                ret[linedat[2]] = '.'
            else:
                ret[linedat[2]] = linedat[3]

    if options.args['logfile']:
        txt = options.args[tag]
        logging.info(txt + ' dictionary loaded from file ' + tx_to_prot_source_file)
    return ret


# Writing header information to output file
def writeHeader(options, header, outfile, stdout, version):
    if options.args['prefix']:
        prefix = 'CAVA_'
    else:
        prefix = ''
    headerinfo = '##INFO=<ID=' + prefix + 'TYPE,Number=.,Type=String,Description=\"Variant type: Substitution, Insertion, Deletion or Complex\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'GENE,Number=.,Type=String,Description=\"HGNC gene symbol\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'TRANSCRIPT,Number=.,Type=String,Description=\"Transcript identifier\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'GENEID,Number=.,Type=String,Description=\"Gene identifier\",Source=\"CAVA\",Version=\"1.2.4\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'TRINFO,Number=.,Type=String,Description=\"Transcript information: Strand/Length of transcript/Number of exons/Length of coding DNA + UTR/Protein length\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'LOC,Number=.,Type=String,Description=\"Location of variant in transcript\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'CSN,Number=.,Type=String,Description=\"CSN annotation\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'PROTPOS,Number=.,Type=String,Description=\"Protein position\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'PROTREF,Number=.,Type=String,Description=\"Reference amino acids\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'PROTALT,Number=.,Type=String,Description=\"Alternate amino acids\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'CLASS,Number=.,Type=String,Description=\"5PU: Variant in 5 prime untranslated region, 3PU: Variant in 3 prime untranslated region, INT: Intronic variant that does not alter splice site bases, SS: Intronic variant that alters a splice site base but not an ESS or SS5 base, ESS: Variant that alters essential splice site base (+1,+2,-1,-2), SS5: Variant that alters the +5 splice site base, but not an ESS base, SY: Synonymous change caused by a base substitution (i.e. does not alter amino acid), NSY: Nonsynonymous change (missense) caused by a base substitution (i.e. alters amino acid), IF: Inframe insertion and/or deletion (variant alters the length of coding sequence but not the frame), IM: Variant that alters the start codon, SG: Variant resulting in stop-gain (nonsense) mutation, SL: Variant resulting in stop-loss mutation, FS: Frameshifting insertion and/or deletion (variant alters the length and frame of coding sequence), EE: Inframe deletion, insertion or base substitution which affects the first or last three bases of the exon\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'SO,Number=.,Type=String,Description=\"Sequence Ontology term\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'ALTFLAG,Number=.,Type=String,Description=\"None: variant has the same CSN annotation regardless of its left or right-alignment, AnnNotClass/AnnNotSO/AnnNotClassNotSO: indel has an alternative CSN but the same CLASS and/or SO, AnnAndClass/AnnAndSO/AnnAndClassNotSO/AnnAndSONotClass/AnnAndClassAndSO: Multiple CSN with different CLASS and/or SO annotations\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'ALTANN,Number=.,Type=String,Description=\"Alternate CSN annotation\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'ALTCLASS,Number=.,Type=String,Description=\"Alternate CLASS annotation\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'ALTSO,Number=.,Type=String,Description=\"Alternate SO annotation\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'IMPACT,Number=.,Type=String,Description=\"Impact group the variant is stratified into\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'DBSNP,Number=.,Type=String,Description=\"rsID from dbSNP\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'HGVSc,Number=.,Type=String,Description=\"HGVS Nomenclature for cDNA changes\",Source=\"CAVA\",Version=\"'+ version +'\">\n'
    headerinfo += '##INFO=<ID=' + prefix + 'HGVSp,Number=.,Type=String,Description=\"HGVS Nomenclature for protein changes\",Source=\"CAVA\",Version=\"'+ version +'\">\n'

    dateline = '##fileDate=' + time.strftime("%Y-%m-%d")

    if options.args['outputformat'] == 'VCF':
        if header == '':
            if stdout:
                print(
                    '##fileformat=VCFv4.1\n' + dateline + '\n' + headerinfo + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
            else:
                outfile.write(
                    '##fileformat=VCFv4.1\n' + dateline + '\n' + headerinfo + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        else:
            headerv = header.split('\n')
            headerm = []
            startline = '##fileformat=VCFv4.1'
            for x in headerv:
                if x.startswith('##fileDate') or x.startswith('##filedate'): continue
                if x.startswith('##fileformat'):
                    startline = x
                    continue
                headerm.append(x)
            if stdout:
                print(startline + '\n' + dateline + '\n' + headerinfo + '\n'.join(headerm) + '\n')
            else:
                outfile.write(startline + '\n' + dateline + '\n' + headerinfo + '\n'.join(headerm) + '\n')

    if options.args['outputformat'] == 'TSV':
        str = 'ID\tCHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tTYPE'
        if (not options.args['ensembl'] == '.') and (not options.args['ensembl'] == ''):
            if options.args[
                'ontology'].upper() == 'CLASS': str += '\tTRANSCRIPT\tGENE\tGENEID\tTRINFO\tLOC\tCSN\tPROTPOS\tPROTREF\tPROTALT\tCLASS'
            if options.args[
                'ontology'].upper() == 'SO': str += '\tTRANSCRIPT\tGENE\tGENEID\tTRINFO\tLOC\tCSN\tPROTPOS\tPROTREF\tPROTALT\tSO'
            if options.args[
                'ontology'].upper() == 'BOTH': str += '\tTRANSCRIPT\tGENE\tGENEID\tTRINFO\tLOC\tCSN\tPROTPOS\tPROTREF\tPROTALT\tCLASS\tSO'

            if not (options.args['impactdef'] == '.' or options.args['impactdef'] == ''): str += '\tIMPACT'

            if options.args['givealt']:
                if options.args['ontology'].upper() == 'CLASS': str += '\tALTANN\tALTCLASS'
                if options.args['ontology'].upper() == 'SO': str += '\tALTANN\tALTSO'
                if options.args['ontology'].upper() == 'BOTH': str += '\tALTANN\tALTCLASS\tALTSO'

            if (not options.args['givealt']) or options.args['givealtflag']:
                str += '\tALTFLAG'

        if (not options.args['dbsnp'] == '.') and (not options.args['dbsnp'] == ''):
            str += '\tDBSNP'

        str += '\tHGVSC\tHGVSP'

        if stdout:
            print(str)
        else:
            outfile.write(str + '\n')


# Counting number of records in a file
def countRecords(filename):
    ret = 0
    if filename.endswith('.gz'):
        inputf = gzip.open(filename, 'rt')
    else:
        inputf = open(filename)
    for line in inputf:
        line = line.strip()
        if not (line.startswith("#") or line == ''): ret += 1
    return ret


# Checking if options are correct
def checkOptions(options):
    # Checking if @inputformat was given correct value
    str = options.args['inputformat'].upper()
    if not (str == 'VCF' or str == 'TXT'):
        print('ERROR: incorrect value of the tag @inputformat.')
        print('(Allowed values: \'VCF\' or \'TXT\')')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @inputformat.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @outputformat was given correct value
    str = options.args['outputformat'].upper()
    if not (str == 'VCF' or str == 'TSV'):
        print('ERROR: incorrect value of the tag @outputformat.')
        print('(Allowed values: \'VCF\' or \'TSV\')')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @outputformat.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @type was given correct value
    str = options.args['type'].upper()
    if not (
            str == 'ALL' or str == 'SUBSTITUTION' or str == 'INDEL' or str == 'INSERTION' or str == 'DELETION' or str == 'COMPLEX'):
        print('ERROR: incorrect value of the tag @type.')
        print('(Allowed values: \'all\', \'substitution\', \'indel\', \'insertion\', \'deletion\' or \'complex\')')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @type.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @ssrange was given correct value
    ssrange = int(options.args['ssrange'])
    if not ssrange >= 6:
        print('ERROR: incorrect value of the tag @ssrange.')
        print('(Minimum value allowed is 6.)')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @ssrange.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @ontology was given correct value
    str = options.args['ontology'].upper()
    if not (str == 'CLASS' or str == 'SO' or str == 'BOTH'):
        print('ERROR: incorrect value of the tag @ontology.')
        print('(Allowed values: \'CLASS\' or \'SO\' or \'both\')')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @ontology.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @reference file exists
    if not os.path.isfile(options.args['reference']):
        print('ERROR: the file given as @reference does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @reference does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @reference index file exists
    if not os.path.isfile(options.args['reference'] + '.fai'):
        print('ERROR: the .fa.fai index file for @reference is not found.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The .fa.fai index file for @reference is not found.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @ensembl file exists
    if not (options.args['ensembl'] == '.' or options.args['ensembl'] == '') and not os.path.isfile(
            options.args['ensembl']):
        print('ERROR: the file given as @ensembl does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @ensembl does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @ensembl index file exists
    if not (options.args['ensembl'] == '.' or options.args['ensembl'] == '') and not os.path.isfile(
            options.args['ensembl'] + '.tbi'):
        print('ERROR: the .gz.tbi index file for @ensembl is not found.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The .gz.tbi index file for @ensembl is not found.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @dbsnp file exists
    if not (options.args['dbsnp'] == '.' or options.args['dbsnp'] == '') and not os.path.isfile(options.args['dbsnp']):
        print('ERROR: the file given as @dbsnp does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @dbsnp does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @dbsnp index file exists
    if not (options.args['dbsnp'] == '.' or options.args['dbsnp'] == '') and not os.path.isfile(
            options.args['dbsnp'] + '.tbi'):
        print('ERROR: the .gz.tbi index file for @dbsnp is not found.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The .gz.tbi index file for @dbsnp is not found.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @target file exists
    if not (options.args['target'] == '.' or options.args['target'] == '') and not os.path.isfile(
            options.args['target']):
        print('ERROR: the file given as @target does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @target does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @target index file exists
    if not (options.args['target'] == '.' or options.args['target'] == '') and not os.path.isfile(
            options.args['target'] + '.tbi'):
        print('ERROR: the .bed.tbi index file for @target is not found.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The .bed.tbi index file for @target is not found.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @genelist file exists
    if not (options.args['genelist'] == '.' or options.args['genelist'] == '') and not os.path.isfile(
            options.args['genelist']):
        print('ERROR: the file given as @genelist does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @genelist does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @transcriptlist file exists
    if not (options.args['transcriptlist'] == '.' or options.args['transcriptlist'] == '') and not os.path.isfile(
            options.args['transcriptlist']):
        print('ERROR: the file given as @transcriptlist does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @transcriptlist does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @snplist file exists
    if not (options.args['snplist'] == '.' or options.args['snplist'] == '') and not os.path.isfile(
            options.args['snplist']):
        print('ERROR: the file given as @snplist does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('The file given as @snplist does not exist.')
            logging.info('No output file written. CAVA quit.')
        quit()

    # Checking if @transcript2protein file exists
    #  
    if not (('transcripts2protein' not in options.args) or options.args['transcript2protein'] == '.' or options.args[
        'transcript2protein'] == '') and not os.path.isfile(options.args['transcript2protein']):
        print('ERROR: the file given as @transcript2protein does not exist.')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        quit()

    # Checking if @normalized_mitochondrial_chrom was given correct value
    str = options.args['normalized_mitochondrial_chrom']
    if str != 'not_normalized' and not (str == 'M' or str == 'MT'):
        print('ERROR: incorrect value of the tag @normalized_mitochondrial_chrom.')
        print('(Allowed values: \'M\' or \'MT\')')
        print('\nNo output file written. CAVA quit.')
        print("--------------------------------------------------------------------\n")
        if options.args['logfile']:
            logging.error('Incorrect value of the tag @normalized_mitochondrial_chrom.')
            logging.info('No output file written. CAVA quit.')
        quit()

#######################################################################################################################
