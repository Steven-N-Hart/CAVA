#!/usr/bin/env python3


# CSN annotation
#######################################################################################################################
import sys
from . import core


# Class representing a CSN annotation
class CSNAnnot:
    # Constructor
    def __init__(self, coord1, intr1, coord2, intr2, dna, protein, coord1_ins, intr1_ins, coord2_ins, intr2_ins,
                 dna_ins):
        self.coord1 = coord1
        self.intr1 = intr1
        self.coord2 = coord2
        self.intr2 = intr2
        self.dna = dna
        self.protein = protein
        self.coord1_ins = coord1_ins
        self.intr1_ins = intr1_ins
        self.coord2_ins = coord2_ins
        self.intr2_ins = intr2_ins
        self.dna_ins = dna_ins
        self.nout1 = None
        self.nout2 = None

    # Getting annotation as a single String
    def getAsString(self):
        if (self.nout1 is not None and self.nout1>0) and (self.nout2 is not None and self.nout2>0):
            return ''   # Currently not used... HGVS does not support annotating variants  outside transcript
            # in this fashion. Rather assume start/end of mRNA are uncertain and use the c.- and c.* notation
        # Adding the first part of the csn annotation (coordinates)
        ret = 'c.' + str(self.coord1)
        if self.intr1 is not None and self.intr1 != 0:
            if self.intr1 > 0:
                ret += '+' + str(self.intr1)
            else:
                ret += str(self.intr1)
        if self.coord2 is not None and (self.coord1 != self.coord2 or (self.coord1 == self.coord2 and  self.intr1 != self.intr2)):
            ret += '_' + str(self.coord2)
            if self.intr2 != 0:
                if self.intr2 > 0:
                    ret += '+' + str(self.intr2)
                else:
                    ret += str(self.intr2)
        # Adding the second part of the csn annotation (DNA and protein level)
        ret += self.dna + self.protein
        return ret

    # Getting annotation as a list of separate values
    # This method is not used..
    def getAsFields(self):
        if not self.coord1_ins == '':
            return self.coord1_ins, self.intr1_ins, self.coord2_ins, self.intr2_ins, self.dna_ins, self.protein
        else:
            return self.coord1, self.intr1, self.coord2, self.intr2, self.dna, self.protein


#######################################################################################################################

# Getting CSN annotation of a given variant
# THIS only gets called with transcript.strand==1 ==> right shifted variants
#                        or transcript.strand==-1 ==> left shifted variants
def getAnnotation(variant, transcript, reference, prot, mutprot):
    # Creating csn annotation coordinates
    coord1, intr1, coord2, intr2, nout1, nout2 = calculateCSNCoordinates(variant, transcript)

    # Creating DNA level annotation
    if variant.alt.startswith("<") and variant.alt.endswith(">") and not ("," in variant.alt): # Tolerate gVCF or any symbolic ID (e.g. CNV)
        dna, dna_ins = 'X', 'X'
        protein, protchange = '', ('.', '.', '.')
        coord1_ins, intr1_ins, coord2_ins, intr2_ins = '', '', '', ''
        csn = CSNAnnot(coord1, intr1, coord2, intr2, dna, protein, coord1_ins, intr1_ins, coord2_ins, intr2_ins, dna_ins)
        csn.nout1 = nout1
        csn.nout2 = nout2
        return csn, protchange

    # Creating protein level annotation
    where = transcript.whereIsThisVariant(variant)
    if not '-' in where and where.startswith('Ex'):
        protein, protchange = makeProteinString(variant, prot, mutprot, coord1)
        skip_repeats = False
    else:  # large variant crossing intron/exon boundary .. or purely in intron
        protein, protchange = '', ('.', '.', '.')
        if (not ('-' in where)) and (where == '3UTR' or where == '5UTR' or where.startswith("In")):
            skip_repeats = False  #Can have repeat in single intron
                                  # Can have repeat in the UTR, even if the UTR spans multiple exons.
        else:
            skip_repeats = True  # Cannot have cDNA repeat across in intron/exon boundary


    try:
        dna, dna_ins = makeDNAannotation(variant, transcript, reference, coord1, intr1, coord2, intr2, nout1, nout2, skip_repeats)
    except TypeError:
        dna, dna_ins = 'X', 'X'

    # Transforming coordinates if the variant is a an insertion
    # of a multi-base repeat so that the range points to the repeat unit..
    # .. as opposed to the insertion site.
    skipped_repeat  = False
    if len(dna_ins)>0 and dna_ins.find('[') >=0: # repeated sequence,repeat_unit[repeat_len]
        coord1_ins, intr1_ins, coord2_ins, intr2_ins = coord1, intr1, coord2, intr2
        [left_result,right_result,full_result] = scan_for_repeat(variant, transcript, reference)
        if transcript.strand == 1:
            range_start = left_result[1]
            repeat_unit = left_result[4]
            n_repeat_ref = left_result[2]
            n_repeat_alt = left_result[3]
            range_end = range_start+len(repeat_unit)*n_repeat_ref -1
            # Check again to make sure that shifted range did not end up in coding region and repeat length is not multiple of 3.
            coord1new, intr1new, nout1new = transformToCSNCoordinate(range_start, transcript)
            coord2new, intr2new, nout2new = transformToCSNCoordinate(range_end, transcript)
            if len(repeat_unit) % 3 != 0 :
                if (range_start>= transcript.codingStartGenomic and range_start <= transcript.codingEndGenomic) or \
                        (range_end>= transcript.codingStartGenomic and range_end <= transcript.codingEndGenomic) :
                    if ((intr1new is None or intr1new ==0) or (intr2new is None or intr2new == 0)): # if not at least one end in CDS
                        skipped_repeat = True
            if skipped_repeat is False:
                coord1, intr1, nout1 = coord1new, intr1new, nout1new
                coord2, intr2, nout2 = coord2new, intr2new, nout2new
                dna = repeat_unit + '[' + str(n_repeat_ref) + ']%3B[' + str(n_repeat_alt) + ']'
            else: # Repeat is not allowed, try next best thing.
                try:
                    dna, dna_ins = makeDNAannotation(variant, transcript, reference, coord1, intr1, coord2, intr2,
                                                     nout1, nout2, True)
                except TypeError:
                    dna, dna_ins = 'X', 'X'
        else: # strand == -1
            range_end = right_result[1]
            repeat_unit = right_result[4]
            n_repeat_ref =right_result[2]
            n_repeat_alt = right_result[3]
            range_start = range_end - (len(repeat_unit)*n_repeat_ref -1)
            coord1new, intr1new, nout1new = transformToCSNCoordinate(range_end, transcript)
            coord2new, intr2new, nout2new = transformToCSNCoordinate(range_start, transcript)
            if len(repeat_unit) % 3 != 0 :
                if (range_start>= transcript.codingStartGenomic and range_start <= transcript.codingEndGenomic) or \
                        (range_end>= transcript.codingStartGenomic and range_end <= transcript.codingEndGenomic) :
                    if ((intr1new is None or intr1new ==0) or (intr2new is None or intr2new == 0)): # if not at least one end in CDS
                        skipped_repeat = True   # Cannot allow not multiple of 3 in CDS
            if skipped_repeat is False:
                dna = core.Sequence(repeat_unit).reverseComplement() + '[' + str(n_repeat_ref) + ']%3B['+ str(n_repeat_alt) + ']'
                coord1, intr1, nout1 = coord1new, intr1new, nout1new
                coord2, intr2, nout2 = coord2new, intr2new, nout2new
            else: # Repeat is not allowed, try next best thing.
                try:
                    dna, dna_ins = makeDNAannotation(variant, transcript, reference, coord1, intr1, coord2, intr2,
                                                     nout1, nout2, True)
                except TypeError:
                    dna, dna_ins = 'X', 'X'
    # shift coordinates if have to point to previous bases for dup or inv.
    if dna_ins == 'dup' or dna_ins == 'inv':
        # Shift position for insertions longer than 1 bp that are "dup"
        # if the repeat unit is "1", coord1 will already point to the bp before the insertion site.
        coord1_ins, intr1_ins, coord2_ins, intr2_ins = coord1, intr1, coord2, intr2
        # Point to range before insertion site, in cDNA coordinates, but do not skip over introns
        coord1, intr1, coord2, intr2, nout1, nout2 = duplicationCoordinates(variant, transcript)
    elif not ('[' in dna_ins): # del, ins*,delins*
        coord1_ins, intr1_ins, coord2_ins, intr2_ins = '', '', '', ''


    # Creating and returning csnAnnot object
    csn = CSNAnnot(coord1, intr1, coord2, intr2, dna, protein, coord1_ins, intr1_ins, coord2_ins, intr2_ins, dna_ins)
    csn.nout1 = nout1
    csn.nout2 = nout2

    return csn, protchange


# Calculating csn annotation coordinates
def calculateCSNCoordinates(variant, transcript):
    # Returning coordinates if variant is a base substitution
    if variant.is_substitution:
        x, y, nout = transformToCSNCoordinate(variant.pos, transcript)
        return x, y, None, None, nout, None

    # Returning coordinates if variant is an insertion
    if variant.is_insertion:
        startx, starty, start_nout = transformToCSNCoordinate(variant.pos - 1, transcript)  # For insertion, variant.pos is one after .. adjusting
        endx, endy, end_nout = transformToCSNCoordinate(variant.pos, transcript)
        if transcript.strand == 1:
            return startx, starty, endx, endy, start_nout, end_nout
        else:
            return endx, endy, startx, starty, end_nout, start_nout

    # Returning coordinates if variant is a deletion
    if variant.is_deletion:
        startx, starty, start_nout  = transformToCSNCoordinate(variant.pos, transcript)
        if len(variant.ref) == 1:
            return startx, starty, None, None, start_nout, None
        endx, endy, end_nout  = transformToCSNCoordinate(variant.pos + len(variant.ref) - 1, transcript)
        # Only return if at least one of the ends is inside transcript or cover the entire gene
        start_in_transcript = False
        if '*' in startx:
            afterlen = int(startx[(startx.index('*')+1):len(startx)])
            if afterlen <= transcript.three_prime_len:
                start_in_transcript = True
        elif '-' in startx:
            beforelen = int(startx[(startx.index('-')+1):len(startx)])
            if beforelen < transcript.codingStart:
                start_in_transcript = True
        elif startx != '0' and startx is not None:  # start inside coding region, it does not matter where endx is
            if transcript.strand == 1:
                return startx, starty, endx, endy, start_nout, end_nout
            else:
                return endx, endy, startx, starty, end_nout ,start_nout
        end_in_transcript = False
        if '*' in endx:
            afterlen = int(endx[(endx.index('*')+1):len(endx)])
            if afterlen <= transcript.three_prime_len:
                end_in_transcript = True
        elif '-' in endx:
            beforelen = int(endx[(1+endx.index('-')):len(endx)])
            if beforelen < transcript.codingStart:
                end_in_transcript = True
        elif endx != '0' and endx is not None: # endx in coding region, so it does not matter where startx is
            if transcript.strand == 1:
                return startx, starty, endx, endy, start_nout, end_nout
            else:
                return endx, endy, startx, starty, end_nout ,start_nout
        if start_in_transcript is True or end_in_transcript is True:  # At least One end in transcript
            if transcript.strand ==1:
                return startx, starty, endx, endy, start_nout, end_nout
            else:
                return endx, endy, startx, starty, end_nout ,start_nout
        if transcript.strand == 1:
            if '-' in startx and '*' in endx:  #Deletion encompassing whole transcript
                return startx, starty, endx, endy, start_nout, end_nout
        else:
            if '*' in startx and '-' in endx:  #Deletion encompassing whole transcript
                return endx, endy, startx, starty, end_nout ,start_nout
        return None, None, None, None, None, None  # Boths ends of variant outside transcript (one side of other)



    # Returning coordinates if variant is a complex indel
    if variant.is_complex:
        startx, starty, start_nout = transformToCSNCoordinate(variant.pos, transcript)
        if len(variant.ref) == 1:
            return startx, starty, None, None, start_nout, 0
        endx, endy, end_nout = transformToCSNCoordinate(variant.pos + len(variant.ref) - 1, transcript)
        if transcript.strand == 1:
            return startx, starty, endx, endy, start_nout, end_nout
        else:
            return endx, endy, startx, starty, end_nout, start_nout

    return None, None, None, None

# find if the alt-allele (or deletion) contains a repeated pattern (N full repeats )
#
# Can you shift a deletion not multiple of local unit (4 base del in a triplet repeat)? No! .. example
# ACGA[CGAC]GACG -> ACGAGACG
# ACG[ACGA]CGACG -> ACGCGACG  .. not same
# AC[GACG]ACGACG -> ACACGACG   .. not same
# A[CGAC]GACGACG -> AGACGACG   .. not same
# However, you can delete multiple copies of the repeat unit
# ACG[ACGACG]ACG -> ACGACG
# AC[GACGAC]GACG -> ACGACG
# A[CGACGA]CGACG -> ACGACG
# [ACGACG]ACGACG -> ACGACG

# This is why this function tries to split a deleted/inserted sequence into elementary repeat.
def find_repeat_unit(extraseq):
    if len(extraseq) == 0:
        return ["", 0]
    elif len(extraseq) == 1:
        return [extraseq, 1]
    else:
        irlen = 0
        nsegs = 0
        for irlen in range(1, min(len(extraseq) - 1, int((len(extraseq) + 1) / 2)) + 1):  # max repeat leg
            roll_seq = ""
            nsegs = int((len(extraseq)) / irlen)
            nsegsmatch = 0
            seg0 = extraseq[0:irlen]
            allsegsmatch = True
            for iseg in range(1, nsegs):
                segi = extraseq[(iseg * irlen):((iseg + 1) * irlen)]
                if seg0 != segi:
                    allsegsmatch = False
                    break
            if allsegsmatch is True:
                roll_seq = extraseq[(nsegs * irlen):]
                if len(roll_seq) > 0 and roll_seq != seg0:
                    allsegsmatch = False
                else:  # success, found full repeat,
                    return [seg0, nsegs]
        return [extraseq, 1]

    # An indel from NGS could be a repeat expansion..
    # Priority to the repeat uniq that leads to the most shifting
    #   inserted or deleted sequence is mono or exact repeat (no partial)
    #   .. if not, perfect repeat
    #              shift all the way right, then all the way left.. and
    #              scan for smallest possible repeat pattern that covers the interval
    #              Once have that pattern, scan for possible expansion of range
def scan_for_repeat(variant, transcript, reference):

    if not (variant.is_deletion or variant.is_insertion):
        return [None,None,None]
    else:
        if variant.is_insertion:
            rep0 = variant.alt
            pos_left_of_variant = variant.pos - 1  # where do you look for previous repeat
            pos_right_of_variant = variant.pos  # Where do you look for next repeat
        else:
            rep0 = variant.ref
            pos_left_of_variant = variant.pos - 1
            pos_right_of_variant = variant.pos + len(variant.ref)
        [rep, nrep] = find_repeat_unit(rep0)  # The insertion or deletion could be more than 1 repeat unit long.

        lseq = ""
        lrep = len(rep)
        left_context = ""
        nrep_left = 0
        # get reference operations are very expensive, better to get a big chunk.
        left_end = pos_left_of_variant + 1  # Just adding +1 to get loop initialized
        next_left = pos_left_of_variant - len(rep) + 1
        match_rep = True
        while match_rep is True:  # Scan for repeats and load bigger chunks of data as scan toward the end.
            left_right_end = left_end - 1
            left_end = left_end - 40 * lrep
            left_context = reference.getReference(variant.chrom, left_end, left_right_end) + left_context
            while match_rep is True and next_left > left_end + lrep:  # scan for repeat in current chunk
                pseq = left_context[next_left - left_end:(next_left + lrep - left_end)]
                if pseq != rep:
                    match_rep = False
                    lseq = pseq
                else:
                    next_left = next_left - lrep
                    nrep_left += 1

        # Now scan all the way to the right.
        rseq = ""
        nrep_right = 0
        # get reference operations are very expensive, better to get a big chunk.
        right_left_end = pos_right_of_variant  # left end of fasta sequence block
        right_end = right_left_end - 1  # initialize so loop works (will need right_left_end to be variant.pos+1)
        next_right = right_left_end  # start position of next repeat
        match_rep = True
        right_context = ""
        while match_rep is True:
            right_left_end = right_end + 1
            right_end = right_end + 40 * lrep
            right_context = right_context + reference.getReference(variant.chrom, right_left_end,
                                                                   right_end)  # positions are inclusive
            while match_rep is True and next_right + lrep - 1 <= right_end - lrep:
                pseq = right_context[next_right - right_left_end:(next_right + lrep - right_left_end)]
                if pseq != rep:
                    match_rep = False
                    rseq = pseq
                else:
                    next_right = next_right + lrep
                    nrep_right += 1
        # We need to check if changing the repeat unit might  allow a little more shifting
        # 3 cases arise
        # 1) extend on both sides, lead to more repeat units
        # e.g. CTATAT[AT]ATAC -> CTATATATAC was described as an AT deletion, CT[AT](3;2)AC or as a TA deletion C[TA](4;3)C
        # 2)  extend on left only, leading to more left-shifting, but same number of repeat units
        #  CTATAT[AT]ATC -> CTATATATC can be described as AT deletion CT[AT](4;3)C or TA deletion C[TA](4;3)C, the latter being more left-shifted
        # 3) extend to right only, leading to more right-shifting, but same number of repeat units.
        #   CATAT[AT]ATAC -> CATATATAC was described as an AT deletion, C[AT](4;3)AC or as CA[TA](4;3)C
        right_pad = 0
        left_pad = 0
        if lrep > 1 and len(lseq)== lrep and len(rseq) == lrep:
            # see if "rolling" the repeat pattern could extend the repeated section
            # We already know that we can only find a partial pattern on left and/or right

            for npad in range(1, lrep):
                # Check left side extension
                if rep[lrep - npad:] == lseq[lrep - npad:]:
                    left_pad = npad
                # Check right side extension of pattern
                if rep[0:npad] == rseq[0:npad]:
                    right_pad = npad

        extra_rep = 0
        if left_pad + right_pad >= lrep:
            extra_rep = 1
#        sys.stdout.write("left_pad="+str(left_pad)+"\n")
#        sys.stdout.write("right_pad="+str(right_pad)+"\n")

        # return fully left_padded variant with ref_repeats, alt_repeats, rep_unit, start_pos
        if variant.is_insertion:
            nrep_ref = 0
            nrep_alt = nrep
        elif variant.is_deletion:
            nrep_ref = nrep
            nrep_alt = 0
        newleft_rep = rep[lrep - left_pad:] + rep[0:lrep - left_pad]
        left_result = [next_left + lrep - left_pad,  # left shifted position of indel for variant obj
                       next_left + lrep - left_pad,  # HGVS of repeat is always
                       nrep_left + nrep_right + extra_rep + nrep_ref,  # Number of ref repeats
                       nrep_left + nrep_right + extra_rep + nrep_alt,  # number of alt repeats
                       newleft_rep,  # repeat_pattern
                       newleft_rep * nrep_ref,  # trimmed ref-allele
                       newleft_rep * nrep_alt]  # trimmed alt-allele
#
        newright_rep = rep[right_pad:] + rep[0:right_pad]
        right_result = [next_right - lrep * nrep + right_pad,  # right shifted position of indel for variant obj
                        next_left + lrep + right_pad,  #
                        nrep_left + nrep_right + extra_rep + nrep_ref,  # Number of ref repeats
                        nrep_left + nrep_right + extra_rep + nrep_alt ,  # number of alt repeats
                        newright_rep,  # repeat_pattern
                        newright_rep * nrep_ref,  # trimmed ref-allele
                        newright_rep * nrep_alt]  # trimmed alt-allele
        full_result = [left_result[0],  # left-most position
                       next_right +lrep-1 + right_pad,  # rightmost position, including repeat sequencue
                       left_context[next_left - left_end + lrep - left_pad:] +
                       variant.ref +
                       right_context[:(next_right - right_left_end + right_pad)],
                       left_context[next_left - left_end + lrep - left_pad:] +
                       variant.alt +
                       right_context[:(next_right - right_left_end + right_pad)]
                       ]  # sequence
        return [left_result,right_result,full_result]

#

                # return fully right_padded variant ref_repeats, alt_repeats, rep_unit, start_pos
                # return fully left+right variants
                # Repeats annotation does not Span introns.. so all this shifting is purely at the DNA level.
                # If the first repeat unit starts 5' (before) of the TSS, we want to shift the position 5'
                # If the first repeat unit starts after the TSS (3') or ends after the 3'end, we want to shift the repeat 3'
                # If the repeat is entirely within transcript, shift it 3'

    #
    #         if len(rep)<len(rep0):
    #         left_context=reference.getReference(variant.chrom, variant.pos - len(insert), variant.pos-1
    #         if transcript.strand == 1:
    #             # Full variant is already shifted 3' .. and should be
    #             before = reference.getReference(variant.chrom, variant.pos - len(insert), variant.pos-1)
    #                 # Checking if variant is a duplication, favor left shifting by 1 position.
    #                 if insert == before:
    #                     if len(insert) > 4:
    #                         return 'dup' + str(len(insert)), 'ins' + insert
    #                     return 'dup' + insert, 'ins' + insert
    #             # Detect dups, but do not shift left.
    #             if variant.pos + len(insert) <= transcript.transcriptEnd:
    #                 after =  reference.getReference(variant.chrom, variant.pos, variant.pos+len(insert)-1)
    #                 if insert == after:
    #                     return 'dup' + insert, ''  # uppercase INS insures position  will remain here.
    #                 return 'ins' + insert, ''
    #             return 'ins' + insert, ''
    #         else:  # transcript.strand == -1
    #             if (variant.pos + (len(insert)-1) <= transcript.transcriptEnd): #' 5' shift (right)
    #                 before = reference.getReference(variant.chrom, variant.pos, variant.pos + len(insert)- 1 )
    #             # Checking if variant is a duplication
    #                 if insert == before:
    #                     if len(insert) > 4: return 'dup' + str(len(insert)), 'ins' + insert.reverseComplement()
    #                     return 'dup' + insert.reverseComplement(), 'ins' + insert.reverseComplement()
    #             if variant.pos - len(insert) >= transcript.transcriptStart:
    #                 after = reference.getReference(variant.chrom, variant.pos - len(insert), variant.pos-1)
    #                 if insert == after:
    #                     return 'dup' +  insert.reverseComplement(), ''  # Uppercase insures variant is not shifted
    #             return 'ins' + insert.reverseComplement(), ''
    # #
    # # Returning DNA level annotation if variant is a deletion
    # if variant.is_deletion:
    #     if len(variant.ref) > 4:
    #         return 'del' + str(len(variant.ref)), ''
    #     if transcript.strand == 1:
    #         return 'del' + variant.ref, ''
    #     else:
    #         return 'del' + core.Sequence(variant.ref).reverseComplement(), ''
    #
    #     nDup = 0
    #     # Loop over potential repeat sizes (Smaller repeat sizes have priority)
    #     for SSRlen in range(1, len(trim_prot) + 1):
    #         # Make sure repeat size is a full multiple of the deletion
    #         nDups = int(len(trim_prot) / SSRlen)
    #         if nDups >= 1 and len(trim_prot) % SSRlen == 0:
    #             repeat_seq = trim_prot[0:SSRlen]
    #             nMatch_del = 0
    #             nMatch_ref = 0
    #             # Check if the repeat pattern matches every base of the trim_prot
    #             for iDup in range(0, nDups):
    #                 if trim_prot[(0 + SSRlen * iDup):(SSRlen * (iDup + 1))] == repeat_seq:
    #                     nMatch_del = nMatch_del + 1
    #                 else:
    #                     break
    #             if nMatch_del == nDups:  # Can have 0... nMatch_ref copies in reference sequence left after deletion
    #                 lowerlim = leftindex - SSRlen - 1
    #                 upperlim = lowerlim + SSRlen
    #                 while lowerlim >= 0 and protcopy[lowerlim:upperlim] == repeat_seq:
    #                     nMatch_ref = nMatch_ref + 1
    #                     lowerlim = lowerlim - SSRlen
    #                     upperlim = lowerlim + SSRlen
    #                 if not (
    #                         nMatch_del == 1 and nMatch_ref == 0 or nMatch_ref + nMatch_del <= 1):  # Deletion of a single lone copy is not repeat polymorphism
    #                     nDup = nDups
    #                     lowerlim = leftindex - SSRlen * nMatch_ref - 1
    #                     upperlim = lowerlim + SSRlen
    #                     break
    #
    #     if nDup == 0:  # No Repeats .. straight Deletion .. not frameshift .. not ins or complex
    #         if len(trim_prot) == 1:
    #             return '_p.' + changeTo3lettersTer(trim_prot[0]) + str(leftindex) + "del", (
    #             str(leftindex), trim_prot, '-')
    #         else:
    #             return '_p.' + changeTo3lettersTer(trim_prot[0]) + str(leftindex) + "_" + changeTo3lettersTer(
    #                 trim_prot[len(trim_prot) - 1]) + str(rightindex) + "del", (
    #                    str(leftindex) + '-' + str(leftindex + len(trim_prot) - 1), trim_prot, '-')
    #     else:
    #         if SSRlen == 1:  # Repeat unit of 1 base.
    #             if nMatch_ref == 0:  # Multiple >1 copies are deleted, leaving empty sequence. Guaranteed to be an interval of at least 2 AA positions to be called SSR
    #                 return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(lowerlim + 1) + "[" + str(
    #                     nMatch_ref + nMatch_del) + "]%3B[" + str(nMatch_ref) + "]", (
    #                        str(leftindex) + '-' + str(rightindex), trim_prot, '-')
    #             else:
    #                 if leftindex == rightindex:
    #                     return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(lowerlim + 1) + "[" + str(
    #                         nMatch_ref + nMatch_del) + "]%3B[" + str(nMatch_ref) + "]", (str(leftindex), trim_prot, '-')
    #                 else:
    #                     return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(lowerlim + 1) + "[" + str(
    #                         nMatch_ref + nMatch_del) + "]%3B[" + str(nMatch_ref) + "]", (
    #                            str(leftindex) + '-' + str(rightindex), trim_prot, '-')
    #
    #         else:  # Multi base repeat
    #             if nMatch_ref == 0:  # Deletion of repeats leaving no copies
    #                 return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(
    #                     lowerlim + 1) + "_" + changeTo3lettersTer(protcopy[upperlim - 1]) + str(upperlim) + "[" + str(
    #                     nMatch_ref + nMatch_del) + "]%3B[" + str(nMatch_ref) + "]", (
    #                        str(leftindex) + "-" + str(rightindex), trim_prot, '-')
    #             else:
    #                 return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(
    #                     lowerlim + 1) + "_" + changeTo3lettersTer(protcopy[upperlim - 1]) + str(upperlim) + "[" + str(
    #                     nMatch_ref + nMatch_del) + "]%3B[" + str(nMatch_ref) + "]", (
    #                        str(leftindex) + '-' + str(rightindex), trim_prot, '-')

#
# In original CAVA version, repeats can only be reported as a duplication.. even if it is part of a more complex repeat
# This is cDNA annotation, so the repeats should be detected within the confines of the transcript sequence
#   In version 20.05 of the HGVS nomenclature, the following note explains the limited scope of repeat annotation for cDNA
# exception: using a coding DNA reference sequence (“c.” description) a Repeated sequence variant description can be used
# only for repeat units with a length which is a multiple of 3, i.e. which can not affect the reading frame.
#   Consequently, use NM_024312.4:c.2692_2693dup and not NM_024312.4:c.2686A[10],
#   use NM_024312.4:c.1741_1742insTATATATA and not NM_024312.4:c.1738TA[6].
# so the original CAVA code is correct unless repeat units are *3
#
# The original CAVA version has one problem, it may call dups across exon/intron boundaries
#   .
#
#When the sequence is an inversion, we have to shift coordinates if the element is a duplication of inversion.
#
# New HGVS annotation, no longer cites the nucleotides AFTER the del, ins,delins,dup, or repeated
#
def makeDNAannotation(variant, transcript, reference, coord1, intr1, coord2, intr2, nout1, nout2, skip_repeats = False):
    if variant.alt.startswith('<') and variant.alt.endswith('>') and not ( "," in variant.alt):  # support <NON_REF>' <*> or SV
        return ''
    # Returning DNA level annotation if variant is a base substitution
    if variant.is_substitution:
        if transcript.strand == 1:
            return variant.ref + '>' + variant.alt, ''
        else:
            return variant.ref.reverseComplement() + '>' + variant.alt.reverseComplement(), ''

    # Returning genomic DNA level annotation if variant is an insertion
    rep_unit = variant.alt
    repn = 1
    if variant.is_insertion:
        [rep_unit, repn] = find_repeat_unit(variant.alt)
    elif variant.is_deletion:
        [rep_unit, repn] = find_repeat_unit(variant.ref)
    replen = len(rep_unit)

    if variant.is_insertion:
        if repn >1 and skip_repeats is False: # Should try to represent as repeat if allowable
            in_coding_region = variant.pos-1 >= transcript.codingStartGenomic and variant.pos<=transcript.codingEndGenomic
            fully_inside_exon = (coord1 is not None and (intr1 is None or intr1 == 0)) and (
                        coord2 is None or (intr2 is None or intr2 == 0))
            at_least_partly_inside_exon = (coord1 is not None and (intr1 is None or intr1 == 0)) or (
                coord2 is not None and (intr2 is None or intr2 == 0))
            # within cDNA, only multiple of 3 can be represented by repeat when fully or inside exon
            if in_coding_region is False or \
                    ((replen % 3 == 0 and fully_inside_exon is True) or  at_least_partly_inside_exon is False):
                if transcript.strand == 1:
                    rep_str = rep_unit+'['+str(repn)+']'
                    return rep_str,rep_str
                else:
                    rep_str = core.Sequence(rep_unit).reverseComplement()+ '[' + str(repn) + ']'
                    return rep_str,rep_str
            # else: Otherwise, try to represent as dup or  -- if all fails -- as dup.
        insert = core.Sequence(variant.alt)
        # the two cDNA coordinate point to 1=before the  insertion site .. and 1 after the insertion site

        # 2nd return value is set to !='' if variant is a dup and position needs to be altered.
        #  ... this will be done later.
        if transcript.strand == 1:
            # If get here, We are guaranteed that the variants is right shitfed if the transcript strand is +
            # so We only have to check for repeats before insertion point.
            # coord1 points to the repeated base i 1-base coordinates (does not point after insertion site like variant.pos)
            #if variant.pos - len(insert) >= transcript.transcriptStart and variant.pos - 1 <= transcript.transcriptEnd:
            if variant.pos - len(insert) >= 0:
                before = reference.getReference(variant.chrom, variant.pos - len(insert), variant.pos - 1)
            else: # right shifted variant insertion at beginning of chromosome, cannot be dup.
                return 'ins'+variant.alt,''
            # Checking if variant is a duplication, but no more
            if insert == before: # Dup
                return 'dup','dup'  # Dups, but coordinate of insertion will be adjusted to point to a range
            if len(insert)==1:
                return 'ins'+variant.alt,''  # cannot be an inversion
            rev = insert.reverseComplement()
            if rev == before:
                return 'inv','inv'
            return 'ins',''
        else:  # transcript.strand == -1
            # If get here, We are guaranteed that the variants is left shifted if the transcript strand is -1
            #if (variant.pos + (len(insert)-1) <= transcript.transcriptEnd): #' 5' shift (right)
            try:
                after = reference.getReference(variant.chrom, variant.pos, variant.pos + len(insert)- 1 )
            except: # Sequence right at the end, cannot be dup
                return 'ins' + variant.alt.reverseComplement(), ''
        # Checking if variant is a duplication
            if insert == after:
                return 'dup', 'dup'
            rev = insert.reverseComplement()
            if rev == after:
                if len(rev)>1:
                    return 'inv','inv'
                else:
                    return 'inv',''  # it's an inversion, but do not need to change coordinates
            return 'ins' + variant.alt.reverseComplement(), ''

    # Returning DNA level annotation if variant is a deletion
    if variant.is_deletion:
        return 'del',''
        #if len(variant.ref) > 4:
        #    return 'del' + str(len(variant.ref)), ''
        #if transcript.strand == 1:
        #    return 'del' + variant.ref, ''
        #else:
        #    return 'del' + core.Sequence(variant.ref).reverseComplement(), ''


    # Returning DNA level annotation if variant is a complex indel
    if variant.is_complex:
        if transcript.strand == 1:
            return 'delins' + variant.alt, ''
        else:
            return 'delins' + variant.alt.reverseComplement(), ''


# Calculating protein level annotation of the variant
# This code does not try to find new initiation sites upstream of Met1 .. if the the Methionine is destroyed
# This (new) code should mostly work even if reference protein is incomplete (e.g. not start with Met1 or end with X/ter).. except that we cannot call an extension unless last base is Ter.
#
# This (new) code has a new behavior for frameshift, different than the original CAVA, namely: the AAREF and AAMUT return values have the protein sequence.
# 5' Edge case: 
#     Anything deleting/mutating the initiating Methionine will be called a '?'
#     A deletion at the 5' end that leaves a Methionine will be called a deletion (choice between deleting first or second methionine
#     .. so will shift consequence 3' .. so initial Methionine will be considered not deleted.
#
#
#     The functions that call this program do not support 5' extensions  (getAnnotation/getCdsSequence).
#          To Support 5' extensions, they would have to scan/predict alternate Methionine start codons upstream
#          Variants that make the 5' end longer will be 
#                   - insertions (if initial Methionine is kept) or 
#                   - delins (if initial Methionine is deleted)
#          be called an insertion.(and the inserted sequence and the reference sequence have to start with the same sequence)
#
# Note: "reference" parameter is not actually used .. since HGVS doesn't allow looking at DNA to make a decision.
#
# HS: CSN returns  extX instead of extTer, that has to be fixed for HGVSp
# 
# Coord1 is 1-based position in the CDS .. and is ONLY used to report Synonymous variants .. it is not used to locate deletions or frameshifts
# 
def makeProteinString(variant, prot, mutprot, coord1_str):
    """

    :param variant:
    :param prot:
    :param mutprot:
    :param coord1_str:
    :return: '_p.Met1?', ('1', prot[0], 'X')
    """
    if isinstance(coord1_str,str):
        try:
            coord1 = int(coord1_str)
        except Exception as e:
            return '', ('.', '.', '.')   # "-" before coding or "*": after Stop codon
    else:
        coord1 = coord1_str
    if len(prot) == 0:
        return '', ('.', '.', '.')
    # this is used to distringuish frameshift from non-frameshifts.. though an Early Stop codon .. will be coded as nonsense
    #      even if it's a frameshift [HGVS NOTE in fs: the shortest frame shift variantis fsTer2, fsTer1  variants are by definition nonsense variants]
    is_not_frameshift = (len(variant.alt) - len(variant.ref)) % 3 == 0

    xindex = mutprot.find('X')
    # Edge case in HGVS, entire protein deleted .. starting
    if len(mutprot)==0:
        return '_p.Met1?', ('1', prot[0], '')
    if mutprot[0]=='X' and prot[0] == 'M':
        return '_p.Met1?', ('1', prot[0], 'X')
    if prot[0] == 'M' and mutprot[0] != 'M':
        return '_p.Met1?', ('1', prot[0], '')
    if mutprot[0] == 'X' and prot[0] != 'X':
        return '_p.?', ('1', prot[0], 'X')

        # Checking if there was no change in protein sequence
    if prot == mutprot:  # to reach this point, both must be len>0
        idx = int(coord1 / 3)
        if coord1 % 3 > 0:
            idx += 1
            if idx > len(prot):
                sys.stderr.write("Protein too short: prot="+prot+"\n in variant="+variant.id+" for coord1="+coord1_str+"\n")
        # Note old CAVA behavior of p.= is incorrect HGVS ... because that means there are no AA change across the entire protein for any variants.
        return '_p.' + changeTo3lettersTer(prot[idx - 1]) + str(idx) + '=', (str(idx), prot[idx - 1], prot[idx - 1])

    # Checking if the variant affects the initiating amino acid
    # ORIGINAL CODE   if prot[0] != mutprot[0]: return '_p.' + changeTo3letters(prot[0]) + '1?', ('1', prot[0], mutprot[0])
    # Change to match HGVS requirements.

    # Trimming common starting substring
    # leftindex is a 1-base position within prot
    protcopy = prot[:]
    mutprotcopy = mutprot[:]



    # XXX-HS rewrote because this consumed the most time.
    isame = -1
    while len(prot) > isame+1 and len(mutprot) > isame+1:
        if prot[isame+1] == mutprot[isame+1]:
            isame = isame + 1
        else:
            break
    if isame>=0:
        prot = prot[isame+1:]
        mutprot= mutprot[isame+1:]

    leftindex = isame +2


    # Any mutation where the first different AA is a Stop is a nonsense mut (wether insertion of deletion or frameshift)
    if len(mutprot)>0 and mutprot[0] == 'X' and len(prot)>0 and prot[0] != 'X':
        return '_p.' + changeTo3lettersTer(prot[0]) + str(leftindex) + changeTo3lettersTer(mutprot[0]), (
        str(leftindex), prot[0], 'X')

    # Trimming common ending substring
    trim_prot = prot
    trim_mutprot = mutprot
    rightindex = len(protcopy)


    ilast = 0
    while len(trim_prot)> ilast  and len(trim_mutprot) > ilast:
        if trim_prot[-(ilast + 1)] == trim_mutprot[-(ilast + 1)]:
            ilast = ilast + 1
        else:
            break
    if ilast > 0:
        trim_prot = trim_prot[:-ilast]
        trim_mutprot = trim_mutprot[:-ilast]
    rightindex -= ilast

    ####
    # Protein variants have to be considered in the following order ( no Inversion for Protein).
    ####
    # Variants Affecting Start Methionine
    # Variants Affecting End.
    # Subst
    # Del (not an out of frame frameshift .. and if long deletion, must not include the 3'end (Stop)
    # Dup are a type of insertion where previous sequence are identical.(In theory could come from a frameshift near the end)
    # VarSSR (Simple Sequence Repeat)
    #       a type of Ins repeated 1-6AA p.Ala2[10];[11] or p.Ala2_Val4[4];[5]
    #       if Del is n-AA, then try all repeat sizes that are multiple of n (e.g. n=6 .. TRY 2,3,6 .. )
    #       In theory could come from a frameshift near the end.
    # Ins (in-frame) - not a frameshift and not a dup or other repeated sequence SSR
    # FrameShift (out of frame del or ins .. and in frame deletions that include the 3'end)
    # Indel/delins (in frame) .. one or more amino acids are replaced with one or more other amino acids and which is not a substitution, frame shift or conversion.(Cannot tell conversion fro short read)
    # Insertion/Deletions/Frameshifts/SSR cannot occur on first AA of last AA

    # If the initial methionine changes, whether from a single or multiple base changes,or frameshift, the HGVS is still p.?
    #      ... We have no way to deal with 5' Extensions .. as Novel/alternate Methionine preciction upstream are not used
    #
    # Missense at first base or deletion or insertion (if initial protein is partial)
    # it doesn't matter if this is a frameshift as deletion of Initial Methionine dominates( also frameshifts cannot occur at first AA)
    # Deleting the Start Codon Dominates the HGVS function over SSR repeats, so don't even check for that either
    # out-of-frame delins at the first (or last) base are NOT frameshift, just p?

    if leftindex == 1 and len(prot) > 0 and len(mutprot) > 0:  # leftindex==1 means first base of ref is different
        if prot[0] == 'M':
            if len(trim_mutprot) == 0:  # Deletion of initial part of protein, including Methionine
                if rightindex != leftindex:
                    return '_p.?', ('1-' + str(rightindex), trim_prot, '-')
                else:
                    return '_p.?', ('1', trim_prot, '-')  # in-frame deletion
            else:
                if rightindex != leftindex:  # del/ins/complex
                    return '_p.?', ('1-' + str(rightindex), trim_prot, trim_mutprot)
                else:
                    return '_p.?', ('1', trim_prot, trim_mutprot)

        else:  # Incomplete reference protein without a Methionine at the Start
            if mutprot[0] == 'M':  # Incomplete reference protein without start becoming a start
                if len(trim_prot) == 1:
                    return '_p.' + changeTo3lettersTer(prot[0]) + '1' + changeTo3lettersTer(mutprot[0]), (
                    '1', trim_prot, trim_mutprot)
                else:
                    return '_p.' + changeTo3lettersTer(prot[0]) + '1' + changeTo3lettersTer(mutprot[0]), (
                    '1-' + str(len(trim_prot)), trim_prot, trim_mutprot)
            else:
                if len(trim_prot) == 1:
                    return '_p.?', ('1', trim_prot, trim_mutprot)
                else:
                    return '_p.?', ('1-' + str(len(trim_prot)), trim_prot, trim_mutprot)

    # Note .. if a deletion is shifted 3' because of a repeat sequence .. it has to be called at this 3' position.
    #    The particular scanning algorithm insures that AA-level variants are also shifted right.
    #
    # Single AA Substitution has Precedence over frameshift, deleting, or insertion or complex (Except for Start Stop Subs)
    #               For the case of Ins of fs.. if an early Stop occur.. and is encodable as a substitution

    # Any mutation where the first different AA is a new Stop ==> priority as nonsense mut (wether ins, del or fs, or complex )
    if len(mutprot)>0 and mutprot[0] == 'X' and len(prot)>0 and prot[0] != 'X':
        return '_p.' + changeTo3lettersTer(prot[0]) + str(leftindex) + changeTo3lettersTer(mutprot[0]), (
        str(leftindex), prot[0], 'X')

    # Checking if the first base altered results in a stop lost mutation
    if len(prot)>0 and prot[0] == 'X':  # Don't have to check for frameshift or SSR if stop codon is first mutated.
        if len(mutprot) == 0:  # 3' Extension of unknown length on reference sequence without reference sequence that can lead to extension
            return '_p.Ter' + str(leftindex) + '?ext*?', (str(leftindex), 'X', '?')
        else:  # by definition or "prot" and left-scanning ==> mutprot[0] != 'X' : # 3' Extension
            nextstop = mutprot.find('X')  # mutprot starts at mutated Stop codon.
            # HS: Changed CSN extX to extTer
            if nextstop != -1: # ext1 is moving the stop codon by 1 position(nextstop == 1) e.g. REF
                return '_p.Ter' + str(leftindex) + changeTo3lettersTer(mutprot[0]) + 'extTer' + str(nextstop), (
                    str(leftindex), 'X', mutprot[0:(nextstop + 1)])
            else:
                return '_p.Ter' + str(leftindex) + changeTo3lettersTer(mutprot[0]) + 'ext*?', (
                str(leftindex), 'X', mutprot)

    #
    # note that HGVS says "deletions starting N-terminal of and including the translation termination (stop) codon are described as Frame shift. .. "
    #  ... irrespective of wether it is an actual frameshift (at the dna level)
    # .. so we have to identify deletion that affect the stop codon (rightindex<len(protcopy) and treat them as frameshift
    # ...Because we need DNA to identify actual frameshift, the standard made some calls.
    #
    # 'deletions starting N-terminal of and including the translation termination (stop) codon are described as Frame shift'
    # 'delins starting N-terminal of and including the translation termination (stop) codon are described as Frame shift'
    # 'insertions extending the full-length amino acid sequence at the C-terminal end with one or more amino acids are described as Extension.'
    # Deletions of Stop Codon can result in finding another Ter in 3'UTR ... and it's impossible to tell a delins from a frameshift .. if they reach the end.
    #
    # Any Variant where the prot and mut_prot end in Methionine.. and it's not a simple 1 bp point mutation
    # is_not_frameshift and rightindex==len(protcopy)-1 and not (len(trim_prot)==1 and len(trim_mutprot)==1)
    # Variant where the last AA in the prot and mutprot are not the same (keep in mind, we already covered the case of the simple 1 bp extension.
    # is_not_frameshift and rightindex==len(protcopy)   and not (len(trim_prot)==1 and len(trim_mutprot)==1) .... no restriction on left index.
    #
    # Mismatch region cannot start at stop codon .. because otherwise it will be an extension (which we covered above)
    # .. if leftindex==len(protcopy) ==> extension

    # Exclude pure insertions (e.g. not delins) that are not frameshift .. unless they  extend AA at 5' end
    # e.g.len(trim_prot)==0 and is_not_frameshift and len(protcopy)<len(mutprot) and
    # Deletions that are not frameshift should still be treated as frameshift if they affect last exon (e.g. len(trim_mutprot)==0 .. OK)
    couldbe_delstop_or_ins_extendingstop = is_not_frameshift and rightindex >= len(protcopy) - 1 and leftindex < len(
        protcopy) and not ((len(trim_prot) == 1 and len(trim_mutprot) == 1) or (
                len(trim_prot) == 0 and len(protcopy) > len(mutprotcopy)))

    # protein Delins are by definition not allowed to impact the stop codon ( and don't look at DNA to check if it's complex.)
    #         if first base mutated is stop(leftindex==len(protcopy), then it's an extension .. and it was caught earlier.
    # if both alleles end with a Ter ... then Ter will be trimmed.. and this will be an internal delins (not to be treated as FS) .. require rightindex==len(protcopy)
    # 1 BP changes are single-base substitutions .. so exclude those
    # Delins that may affect Stop are already included in couldbe_delstop
    #    delins_tobetreatedasFS = is_not_frameshift and len(trim_mutprot)>=1 and len(trim_prot)>=1 and leftindex<len(protcopy) and rightindex>=len(protcopy)-1 and not (len(trim_prot)==1 and len(trim_mutprot)==1)

    #
    # Finally, we can call single-base substitutions ( after dealing with start and stops changes)
    #   Note that if (somehow?) a frameshift result in a single-base substitution .. it will be called a substitutions
    if len(trim_prot) == 1 and len(trim_mutprot) == 1:
        return '_p.' + changeTo3lettersTer(trim_prot) + str(leftindex) + changeTo3lettersTer(trim_mutprot), (
        str(leftindex), trim_prot, trim_mutprot)
    #
    # For a Variant to be called a repeat change, the entire event must be describable as a repeat variation.. There is no repeat+1 base change
    #
    # Pure Deletion Repeat Event (not a frameshift) that does not include the beginning or end (because such highly impactfull changes were dealth earlier)
    #

    if rightindex < len(protcopy) and len(trim_mutprot) == 0 and is_not_frameshift:
        nDup = 0
        # Loop over potential repeat sizes (Smaller repeat sizes have priority)
        for SSRlen in range(1, len(trim_prot) + 1):
            # Make sure repeat size is a full multiple of the deletion
            nDups = int(len(trim_prot) / SSRlen)
            if nDups >= 1 and len(trim_prot) % SSRlen == 0:
                repeat_seq = trim_prot[0:SSRlen]
                nMatch_del = 0
                nMatch_ref = 0
                # Check if the repeat pattern matches every base of the trim_prot
                for iDup in range(0, nDups):
                    if trim_prot[(0 + SSRlen * iDup):(SSRlen * (iDup + 1))] == repeat_seq:
                        nMatch_del = nMatch_del + 1
                    else:
                        break
                if nMatch_del == nDups:  # Can have 0... nMatch_ref copies in reference sequence left after deletion
                    lowerlim = leftindex - SSRlen - 1
                    upperlim = lowerlim + SSRlen
                    while lowerlim >= 0 and protcopy[lowerlim:upperlim] == repeat_seq:
                        nMatch_ref = nMatch_ref + 1
                        lowerlim = lowerlim - SSRlen
                        upperlim = lowerlim + SSRlen
                    if not (
                            nMatch_del == 1 and nMatch_ref == 0 or nMatch_ref + nMatch_del <= 1):  # Deletion of a single lone copy is not repeat polymorphism
                        nDup = nDups
                        lowerlim = leftindex - SSRlen * nMatch_ref - 1
                        upperlim = lowerlim + SSRlen
                        break

        if nDup == 0:  # No Repeats .. straight Deletion .. not frameshift .. not ins or complex
            if len(trim_prot) == 1:
                return '_p.' + changeTo3lettersTer(trim_prot[0]) + str(leftindex) + "del", (
                str(leftindex), trim_prot, '-')
            else:
                return '_p.' + changeTo3lettersTer(trim_prot[0]) + str(leftindex) + "_" + changeTo3lettersTer(
                    trim_prot[len(trim_prot) - 1]) + str(rightindex) + "del", (
                       str(leftindex) + '-' + str(leftindex + len(trim_prot) - 1), trim_prot, '-')
        else:
            if SSRlen == 1:  # Repeat unit of 1 base.
                if nMatch_ref == 0:  # Multiple >1 copies are deleted, leaving empty sequence. Guaranteed to be an interval of at least 2 AA positions to be called SSR
                    return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(lowerlim + 1) + "[" + str(
                        nMatch_ref + nMatch_del) + "]%3B[" + str(nMatch_ref) + "]", (
                           str(leftindex) + '-' + str(rightindex), trim_prot, '-')
                else:
                    if leftindex == rightindex:
                        return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(lowerlim + 1) + "[" + str(
                            nMatch_ref + nMatch_del) + "]%3B[" + str(nMatch_ref) + "]", (str(leftindex), trim_prot, '-')
                    else:
                        return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(lowerlim + 1) + "[" + str(
                            nMatch_ref + nMatch_del) + "]%3B[" + str(nMatch_ref) + "]", (
                               str(leftindex) + '-' + str(rightindex), trim_prot, '-')

            else:  # Multi base repeat
                if nMatch_ref == 0:  # Deletion of repeats leaving no copies
                    return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(
                        lowerlim + 1) + "_" + changeTo3lettersTer(protcopy[upperlim - 1]) + str(upperlim) + "[" + str(
                        nMatch_ref + nMatch_del) + "]%3B[" + str(nMatch_ref) + "]", (
                           str(leftindex) + "-" + str(rightindex), trim_prot, '-')
                else:
                    return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(
                        lowerlim + 1) + "_" + changeTo3lettersTer(protcopy[upperlim - 1]) + str(upperlim) + "[" + str(
                        nMatch_ref + nMatch_del) + "]%3B[" + str(nMatch_ref) + "]", (
                           str(leftindex) + '-' + str(rightindex), trim_prot, '-')

    #
    # Pure insertion/Repeat/Dup .. not frameshift Extension
    #
    # Insertions Repeats (or Dups) require that "a sequence where, compared to a reference sequence, a segment of one or more amino acids (the repeat unit) is present several times, one after the other.."
    #   so we cannot have nMatch_ref==0
    #   Make sure insertion does not include multiple "Ter"
    xindex = trim_mutprot.find('X')
    if leftindex > 1 and rightindex < len(protcopy) and len(trim_prot) == 0 and len(
            trim_mutprot) > 0 and is_not_frameshift:
        nDup = 0
        # Loop over potential repeat sizes (Smaller repeat sizes have priority)
        for SSRlen in range(1, len(trim_mutprot) + 1):
            # Make sure repeat size is a full multiple of the deletion
            nDups = int(len(trim_mutprot) / SSRlen)
            if nDups >= 1 and len(trim_mutprot) % SSRlen == 0:
                repeat_seq = trim_mutprot[0:SSRlen]
                nMatch_ins = 0
                nMatch_ref = 0
                # Check if the repeat pattern matches every base of the trim_prot
                for iDup in range(0, nDups):
                    if trim_mutprot[(0 + SSRlen * iDup):(SSRlen * (iDup + 1))] == repeat_seq:
                        nMatch_ins = nMatch_ins + 1
                    else:
                        break
                if nMatch_ins == nDups:  # Can have 1... nMatch_ref copies in reference sequence (at least 1 copy is requires for insertions)
                    lowerlim = leftindex - SSRlen - 1
                    upperlim = lowerlim + SSRlen
                    #                         Assume that Variant was right shifted to begin with
                    while lowerlim >= 0 and protcopy[lowerlim:upperlim] == repeat_seq:
                        nMatch_ref = nMatch_ref + 1
                        lowerlim = lowerlim - SSRlen
                        lowerlim + SSRlen
                    if nMatch_ref > 0:  # We require a copy to be present on the reference "protein" .. as per HGVS
                        nDup = nDups
                        lowerlim = leftindex - SSRlen * nMatch_ref - 1
                        upperlim = lowerlim + SSRlen
                        break  # break out of the for iDup loop

        if xindex != -1 or nDup == 0:  # Insertion contain Ter (don't allow repeats with a Ter-containing pattern)
            # No Repeats .. straight Insertion .. not frameshift .. not Deletion or complex
            # insertion sequence include 1 AA  before.(leftindex is guaranteed to be 2 or more)
            # Make sure inserted sequence does not include a Stop codon.
            if xindex == -1:
                return '_p.' + changeTo3lettersTer(protcopy[leftindex - 2]) + str(
                    leftindex - 1) + '_' + changeTo3lettersTer(protcopy[leftindex - 1]) + str(
                    leftindex) + "ins" + changeTo3lettersTer(trim_mutprot), (
                       str(leftindex - 1) + "-" + str(rightindex + 1), '-', trim_mutprot)
            else:  # Inserted protein includes a Stop (but not as first sequence ... that would have been caught earlier)
                # also, because rightindex<len(protcopy) .. we know last Ter matches.
                return '_p.' + changeTo3lettersTer(protcopy[leftindex - 2]) + str(
                    leftindex - 1) + '_' + changeTo3lettersTer(protcopy[rightindex]) + str(
                    rightindex + 1) + "ins" + changeTo3lettersTer(trim_mutprot[0:(xindex + 1)]), (
                       str(leftindex - 1) + '-' + str(rightindex + 1), '-', trim_mutprot[0:(xindex + 1)])
        else:  # Single and Multi-base repeat insertion
            if nDup == 1 and nMatch_ref == 1:  # Duplications, special treatment
                if SSRlen == 1:
                    return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(lowerlim + 1) + "dup", (
                    str(leftindex - 1) + '-' + str(rightindex + 1), '-', trim_mutprot)
                else:
                    return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(
                        lowerlim + 1) + '_' + changeTo3lettersTer(protcopy[upperlim - 1]) + str(upperlim) + "dup", (
                           str(leftindex - 1) + '-' + str(rightindex + 1), '-', trim_mutprot)
            else:
                if SSRlen == 1:
                    return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(lowerlim + 1) + "[" + str(
                        nMatch_ref) + "]%3B[" + str(nMatch_ref + nMatch_ins) + "]", (
                           str(leftindex - 1) + '-' + str(rightindex + 1), '-', trim_mutprot)
                else:
                    return '_p.' + changeTo3lettersTer(protcopy[lowerlim]) + str(
                        lowerlim + 1) + '_' + changeTo3lettersTer(protcopy[upperlim - 1]) + str(upperlim) + "[" + str(
                        nMatch_ref) + "]%3B[" + str(nMatch_ref + nMatch_ins) + "]", (
                           str(leftindex - 1) + '-' + str(rightindex + 1), '-', trim_mutprot)

    # Frameshift mutations (assume len(prot)>0 from now on)
    # New behavior for Frameshift , PROTALT will be the sequence (not just a ".")
    #
    # Also include Some "In-phase deletions or delins:
    #       HGVS Deletion definition: "deletions starting N-terminal of and including the translation termination (stop) codon are described as Frame shift. .. "
    #  ... irrespective of wether it is an actual frameshift (at the dna level) .. this includes Complex Variants that affect Stop too.
    #
    #       HGVS delins definition: deletion-insertion variants starting N-terminal () of and including the translation termination (stop) codon are described as frame shift.
    #       Note Mutated protein can have multiple 'X' (0 or more) after the frameshift location.
    if (len(variant.alt) - len(variant.ref)) % 3 > 0 or couldbe_delstop_or_ins_extendingstop:
        if len(mutprot) == 0:  # Deletion until the end of the protein and whole UTR (since no extra AA picked up from UTR)
            #         Since len(prot)>0, then the variant is a deletion (could NOT be a frameshift that causes an early Stop .. because len(mutprot)==0
            #         It is possible that the DNA variant is a complex variant and not a del, but that is OK (protein annotations should not consider DNA).
            if len(prot) == 1:  # Unless last AA is not Ter .. then this should have been dealth with by Code above
                return '_p.', changeTo3lettersTer(prot[0]) + str(leftindex) + "del", (str(leftindex), prot, '-')
            else:
                return '_p.', changeTo3lettersTer(prot[0]) + str(leftindex) + "_" + changeTo3lettersTer(
                    prot[len(prot) - 1]) + str(rightindex) + "del", (str(leftindex) + '-' + str(rightindex), prot, '-')
        else:
            xindex = mutprot.find("X")
            if xindex != -1:  # p.(Arg123LysfsTer34)
                # Special Rules for when frameshift causes the Initial Met to become Stop
                if xindex == 0 and leftindex == 1 and protcopy[0] == "M":
                    return '_p.(Met1Ter)', ('1', 'M', 'X')
                # Check to see if new Stop is a 3' extension.. can only occur if last AA is the first one changed.
                if leftindex == len(protcopy):  # last AA of protein is first mutated ==> Extensions
                    if xindex + 1 > len(prot):  # extension, at or past end of current protein.
                        # do not check to see if last original codon is X,x,* .. to allow partial reference.
                        return '_p.' + changeTo3lettersTer(prot[0]) + str(len(protcopy)) + changeTo3lettersTer(
                            mutprot[0]) + 'extTer' + str(xindex + 1), (str(leftindex), prot[0], mutprot[0:(xindex + 1)])
                    # else xindex+1<=len(prot) .. only xindex+1==len(prot) possible .. which would mean that the new Ter is at the end of the protein.(only possible is original protein is partial .. considered below)
                if xindex + 1 == len(
                        prot):  # Stop codon in mutprot at end of original protein .. either no change to original stop codon, or original prot did not have a stop codon there.
                    if len(prot) == 1:  # Extension original prot did not have Stop there.. but mutated one does
                        return '_p.' + changeTo3lettersTer(prot[0]) + str(len(protcopy)) + changeTo3lettersTer(
                            mutprot[0]) + 'extTer' + str(xindex + 1), (str(leftindex), prot[0], mutprot[0:(xindex + 1)])
                    # else Frameshift occuring before the end of the protein with the Ter occuring a few AA later at the original prot position (will be dealt below)
                if xindex == 0 and len(
                        prot) > 1:  # First Modified Base is a stop codon before end (Independent wether mutation is deletion or insertion
                    return '_p.' + changeTo3lettersTer(prot[0]) + str(leftindex) + 'Ter', (
                    str(leftindex), prot[0], mutprot[0])
                if leftindex < len(protcopy) and xindex + 1 <= len(
                        prot):  # First Modified Base before the end of the protein .. and new Stop not past original Stop (e.g. not an extension)
                    # (btw Minimal length of prot==2 for protein ending in Ter)
                    if len(trim_prot) == 0:  # insertion
                        #                       Duplications notation has priority over insertion (except if the Stop is within the new allele)
                        #       Note that because len(trim_mutprot)>0 ==> leftindex>1 is required to be an insertion
                        #         also note that anything more than a Dup .. was dealt before
                        if len(trim_mutprot) > 0 and (leftindex - len(trim_mutprot)) > 0 and protcopy[leftindex - len(
                                trim_mutprot) - 1:leftindex - 1] == trim_mutprot:
                            if len(trim_mutprot) == 1:
                                return '_p.' + changeTo3lettersTer(protcopy[leftindex - len(trim_mutprot) - 1]) + str(
                                    leftindex - len(trim_mutprot)) + 'dup', (
                                       str(leftindex - 1) + '-' + str(rightindex + 1), '-', trim_mutprot)
                            else:
                                return '_p.' + changeTo3lettersTer(protcopy[leftindex - len(trim_mutprot) - 1]) + str(
                                    leftindex - len(trim_mutprot)) + '_' + changeTo3lettersTer(
                                    protcopy[leftindex - 2]) + str(
                                    leftindex - 1) + 'dup', (
                                       str(leftindex - 1) + '-' + str(rightindex + 1), '-', trim_mutprot)

                    # insertion or deletion frameshift starting and ending before the original protein end.
                    return '_p.' + changeTo3lettersTer(prot[0]) + str(leftindex) + changeTo3lettersTer(
                        mutprot[0]) + 'fsTer' + str(xindex + 1), (str(leftindex), prot[0], mutprot[0:(xindex + 1)])
                else:  # xindex+1 > len(prot) .. so new protein past end .. and since  extension already returned above, so leftindex<len(protcopy).
                    return '_p.' + changeTo3lettersTer(prot[0]) + str(leftindex) + changeTo3lettersTer(
                        mutprot[0]) + 'fsTer' + str(xindex + 1), (str(leftindex), prot[0], mutprot[0:(xindex + 1)])
            else:  # No Stop codon until the end of the transcript.. Nevertheless Need to supply the protein, so we can test if nonsense Mediated Decay occurs.
                if leftindex == len(protcopy):  # Extension without predicted Ter
                    return '_p.(' + changeTo3lettersTer(prot[0]) + str(leftindex) + changeTo3lettersTer(
                        mutprot[0]) + 'ext*?', (str(leftindex), prot[0], mutprot)
                else:
                    return '_p.(' + changeTo3lettersTer(prot[0]) + str(leftindex) + changeTo3lettersTer(
                        mutprot[0]) + 'fs*?)', (str(leftindex), prot[0], mutprot)

    # Past this point, only Complex delins not including the stop codon, 5' extension insertions


    #  (len(prot)==0 and len(mutprot)==0) should have been covered by  earlier cases

    # deletion at the ends (Already checked for "Internal deletions" and deleted first base Methionine earlier)
    if len(trim_mutprot) == 0:
        if len(trim_prot) == 1:
            return '_p.' + changeTo3lettersTer(trim_prot) + str(leftindex) + 'del', (str(leftindex), trim_prot, '-')
        else:
            protpos = str(leftindex) + '-' + str(rightindex)
            return '_p.' + changeTo3lettersTer(trim_prot[0]) + str(leftindex) + '_' + changeTo3lettersTer(
                trim_prot[-1]) + str(
                rightindex) + 'del', (protpos, trim_prot, '-')

    # Insertion-like, but affecting first or last.
    if len(prot) == 0:
        # New functionality (includes some 5' extension) when leftindex==1 and len(prot)==0
        if leftindex == 1:  # Inserted sequences before the Original Methionine
            if trim_mutprot[0] == 'M':  # Extension
                if protcopy[
                    0] == 'M':  # insertion e.g. MA-> MRMA : M1_A2insRM (prot="" -> trim_mutprot=MR)  .. Only this one is realistic.
                    return '_p.' + changeTo3lettersTer(protcopy[0]) + '1_' + changeTo3lettersTer(
                        protcopy[1]) + '2ins' + changeTo3lettersTer(trim_mutprot[1:]) + changeTo3lettersTer('M'), (
                               '1-2', 'M', trim_mutprot + 'M')
                else:  # e.g. A->MRA (prot="" --> trim_mutprot="MR")
                    return '_p.?', ('1', protcopy[0], trim_mutprot + protcopy[0])
            else:  # MAR->YMAR (""->A) or AR->GAR (""->G)
                return '_p.?', ('1', '-', trim_mutprot)
        if leftindex == len(
                protcopy):  # Only case is if the original protein did not have a Ter e.g. ARG->ARGYX -> prot="", trim_mutprot="YX"
            # Not covered by the HGVS standard, unless it's a dup .. and that was covered above
            return '_p.?', (leftindex, '-', trim_mutprot)

    # Checking if variant results in a complex change
    if len(trim_prot) > 0 and len(trim_mutprot) > 0:
        xindex = trim_mutprot.find('X')
        if xindex != -1:
            trim_mutprot = trim_mutprot[0:(xindex + 1)]
        ret = '_p.'
        if len(trim_prot) == 1:
            ret += changeTo3lettersTer(trim_prot) + str(leftindex)
        else:
            ret += changeTo3lettersTer(trim_prot[0]) + str(leftindex) + '_' + changeTo3lettersTer(trim_prot[-1]) + str(
                rightindex)
        ret += 'delins'

        ret += changeTo3lettersTer(trim_mutprot)

        if leftindex == rightindex:
            protpos = str(leftindex)
        else:
            protpos = str(leftindex) + '-' + str(rightindex)
        return ret, (protpos, trim_prot, trim_mutprot)
    sys.stderr.write(
        "\nBUG: Cannot compute CSN for : original prot=" + protcopy + "\n   : Cannot compute CSN for : mutated  prot=" + mutprotcopy + "\n")

    return "", ()


# Transforming a genomic position to csn coordinate
def transformToCSNCoordinate(pos, transcript):
    prevExonEnd = 99999999

    # Checking if genomic position is not outside CDS
    if not transcript.isPositionOutsideCDS(pos):
        sumOfExonLengths = -transcript.codingStart + 1
        # Iterating through exons
        for i in range(len(transcript.exons)):
            exon = transcript.exons[i]
            if i > 0:
                if transcript.strand == 1:
                    # Checking if genomic position is within intron
                    if prevExonEnd < pos < exon.start + 1:
                        if pos <= int((exon.start + 1 - prevExonEnd) / 2) + prevExonEnd:  # count from previous exon
                            x, y, nout = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, pos - prevExonEnd, nout
                        else:   # count from exon coming up
                            x, y, nout = transformToCSNCoordinate(exon.start + 1, transcript)
                            return x, pos - exon.start - 1, nout
                else:
                    # Checking if genomic position is within intron
                    if exon.end < pos < prevExonEnd:
                        if pos >= int((prevExonEnd - exon.end + 1) / 2) + exon.end:
                            x, y, nout = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, prevExonEnd - pos, nout
                        else:
                            x, y, nout = transformToCSNCoordinate(exon.end, transcript)
                            return x, exon.end - pos, nout
            # Checking if genomic position is within exon
            if exon.start + 1 <= pos <= exon.end:
                if transcript.strand == 1:
                    return str(sumOfExonLengths + pos - exon.start), 0, 0
                else:
                    return str(sumOfExonLengths + exon.end - pos + 1), 0, 0
            # Calculating sum of exon lengths up to this point
            sumOfExonLengths += exon.length
            if transcript.strand == 1:
                prevExonEnd = exon.end
            else:
                prevExonEnd = exon.start + 1
    # If genomic position is outside CDS
    else:
        sumpos = 0
        lastexon  = len(transcript.exons)-1
        noutside = 0
        for i in range(lastexon+1):
            exon = transcript.exons[i]
            if i > 0:  # checking is position is within intron before current exon
                if transcript.strand == 1:
                    # Checking if genomic position is within intron
                    if prevExonEnd < pos < exon.start + 1:
                        if pos <= int((exon.start + 1 - prevExonEnd) / 2) + prevExonEnd:
                            x, y, nout = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, pos - prevExonEnd, nout
                        else:
                            x, y, nout = transformToCSNCoordinate(exon.start + 1, transcript)
                            return x, pos - exon.start - 1, nout
                else:
                    # Checking if genomic position is within intron
                    if exon.end < pos < prevExonEnd:
                        if pos >= int((prevExonEnd - exon.end + 1) / 2) + exon.end:
                            x, y, nout = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, prevExonEnd - pos, nout
                        else:
                            x, y, nout  = transformToCSNCoordinate(exon.end, transcript)
                            return x, exon.end - pos, nout

            if transcript.strand == 1:
                if pos > transcript.codingEndGenomic:
                    if transcript.codingEndGenomic < exon.start + 1 and exon.end < pos:  # coding end is in an earlier exon
                        sumpos += exon.length
                        if i==lastexon: # Variant past end of transcript
                            noutside += pos-exon.end
                    elif exon.contains(transcript.codingEndGenomic) and exon.end < pos: # coding end is in this exon
                        sumpos += exon.end - transcript.codingEndGenomic + 1
                        if i==lastexon: # Variant past end of transcript
                            noutside += pos-exon.end
                    elif exon.contains(transcript.codingEndGenomic) and exon.contains(pos):
                        sumpos += pos - transcript.codingEndGenomic
                    elif transcript.codingEndGenomic < exon.start + 1 and exon.contains(pos):
                        sumpos += pos - exon.start - 1
                    # no else, all cases covered
                if pos < transcript.codingStartGenomic:
                    if pos < exon.start + 1 and exon.end < transcript.codingStartGenomic:  # coding start is in a later exon
                        sumpos += exon.length
                        if i == 0:  # pos 5' of TSS
                            noutside += (exon.start + 1 - pos)
                    elif pos < exon.start + 1 and exon.contains(transcript.codingStartGenomic):
                        sumpos += transcript.codingStartGenomic - exon.start
                        if i==0:  # pos 5' of TSS
                            noutside += (exon.start + 1 - pos)
                    elif exon.contains(pos) and exon.contains(transcript.codingStartGenomic):
                        sumpos += transcript.codingStartGenomic - pos
                    elif exon.contains(pos) and exon.end < transcript.codingStartGenomic:
                        sumpos += exon.end - pos

            if transcript.strand == -1:  # exons are from higher coordinate to lower, but still start<end
                if pos < transcript.codingEndGenomic:
                    if pos < exon.start + 1 and exon.end < transcript.codingEndGenomic:  # coding end is in an earlier (higher coodginate) exon
                        sumpos += exon.length
                        if i==lastexon:  # variant is downstream of transcript
                            noutside += (exon.start+1 - pos)
                    elif pos < exon.start + 1 and exon.contains(transcript.codingEndGenomic): #
                        sumpos += transcript.codingEndGenomic - exon.start
                        if i==lastexon: # variant is downstream of transcript
                            noutside += (exon.start+1 - pos)
                    elif exon.contains(pos) and exon.contains(transcript.codingEndGenomic):  #pos < transcript.codingEndGenomic by default
                        sumpos += transcript.codingEndGenomic - pos
                    elif exon.contains(pos) and exon.end < transcript.codingEndGenomic:
                        sumpos += exon.end - pos
                if transcript.codingStartGenomic < pos:
                    if transcript.codingStartGenomic < exon.start + 1 and exon.end < pos:  # coding start is in a later (downstream) exon
                        sumpos += exon.length
                        if i==0: # variant is 5' of TSS
                            noutside += (pos - exon.end)
                    elif exon.contains(transcript.codingStartGenomic) and exon.end < pos:
                        sumpos += exon.end - transcript.codingStartGenomic + 1
                        if i==0: # variant is 5' of TSS
                            noutside += (pos - exon.end)
                    elif exon.contains(transcript.codingStartGenomic) and exon.contains(pos):
                        sumpos += pos - transcript.codingStartGenomic
                    elif transcript.codingStartGenomic < exon.start + 1 and exon.contains(pos):
                        sumpos += pos - exon.start - 1

            if transcript.strand == 1:
                prevExonEnd = exon.end
            else:
                prevExonEnd = exon.start + 1

        #  .. If decide to not annotate variants outside transcript, use this
        #  .. but according to HGVS, it is OK to refere to variants outside transcripts
        #  ... as long as it's relative to the CDS start/end (and not relative to transcripts annotated ends which can
        #  ...   be uncertain)
        # if transcript.strand == 1:
        #     if pos > transcript.codingEndGenomic:  return '*' + str(sumpos), 0, noutside
        #     if pos < transcript.codingStartGenomic: return '-' + str(sumpos), 0, noutside
        # else:
        #     if pos < transcript.codingEndGenomic: return '*' + str(sumpos), 0, noutside
        #     if pos > transcript.codingStartGenomic: return '-' + str(sumpos), 0, noutside
        #
        # return str(sumpos), 0, noutside
        if transcript.strand == 1:
            if pos > transcript.codingEndGenomic:
                return '*' + str(sumpos+noutside), 0, 0
            if pos < transcript.codingStartGenomic:
                return '-' + str(sumpos+noutside), 0, 0
        else:
            if pos < transcript.codingEndGenomic:
                return '*' + str(sumpos+noutside), 0, 0
            if pos > transcript.codingStartGenomic:
                return '-' + str(sumpos+noutside), 0, 0

        return str(sumpos), 0, 0


# Calculating csn coordinates for duplications.. by 3' shifting by
def duplicationCoordinates(variant, transcript):
    if transcript.strand == 1:
        coord1, intr1, nout1 = transformToCSNCoordinate(variant.pos - len(variant.alt), transcript)
        if len(variant.alt) == 1:
            coord2, intr2, nout2 = None, None, None
        else:
            coord2, intr2, nout2 = transformToCSNCoordinate(variant.pos - 1, transcript)
    else:
        coord1, intr1, nout1 = transformToCSNCoordinate(variant.pos + len(variant.alt) - 1, transcript)
        if len(variant.alt) == 1:
            coord2, intr2, nout2 = None, None, None
        else:
            coord2, intr2, nout2 = transformToCSNCoordinate(variant.pos, transcript)
    return coord1, intr1, coord2, intr2, nout1, nout2




# Changing protein sequence of 1-letter amino acid code to 3-letter code
def changeTo3letters(aas):
    ret = ''
    codes = {
        'I': 'Ile', 'M': 'Met', 'T': 'Thr', 'N': 'Asn',
        'K': 'Lys', 'S': 'Ser', 'R': 'Arg', 'L': 'Leu',
        'P': 'Pro', 'H': 'His', 'Q': 'Gln', 'V': 'Val',
        'A': 'Ala', 'D': 'Asp', 'E': 'Glu', 'G': 'Gly',
        'F': 'Phe', 'Y': 'Tyr', 'C': 'Cys', 'W': 'Trp',
        '*': 'X', 'X': 'X', 'x': 'X', '?': '?', 'U' : 'Sel'}
    for aa in aas: ret += codes[aa]
    return ret


# Changing protein sequence of 1-letter amino acid code to 3-letter code, but with Stops changed to "Ter"
def changeTo3lettersTer(aas):
    ret = ''
    # 5/13/21 Updated to account for "?"
    codes = {
        'I': 'Ile', 'M': 'Met', 'T': 'Thr', 'N': 'Asn',
        'K': 'Lys', 'S': 'Ser', 'R': 'Arg', 'L': 'Leu',
        'P': 'Pro', 'H': 'His', 'Q': 'Gln', 'V': 'Val',
        'A': 'Ala', 'D': 'Asp', 'E': 'Glu', 'G': 'Gly',
        'F': 'Phe', 'Y': 'Tyr', 'C': 'Cys', 'W': 'Trp',
        '*': 'Ter', 'X': 'Ter', 'x': 'Ter', '?': '?', 'U' : 'Sel'}
    for aa in aas: ret += codes[aa]
    return ret

#######################################################################################################################
