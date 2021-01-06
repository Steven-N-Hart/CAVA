#!/usr/bin/env python3


# CSN annotation
#######################################################################################################################
import cava_.core as core



# Class representing a CSN annotation
class CSNAnnot():
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

    # Getting annotation as a single String
    def getAsString(self):
        # Adding the first part of the csn annotation (coordinates)
        ret = 'c.' + str(self.coord1)
        if self.intr1 != 0:
            if self.intr1 > 0:
                ret += '+' + str(self.intr1)
            else:
                ret += str(self.intr1)
        if self.coord2 is not None:
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
    def getAsFields(self):
        if not self.coord1_ins == '':
            return self.coord1_ins, self.intr1_ins, self.coord2_ins, self.intr2_ins, self.dna_ins, self.protein
        else:
            return self.coord1, self.intr1, self.coord2, self.intr2, self.dna, self.protein


#######################################################################################################################

# Getting CSN annotation of a given variant
def getAnnotation(variant, transcript, reference, prot, mutprot):
    # Creating csn annotation coordinates
    coord1, intr1, coord2, intr2 = calculateCSNCoordinates(variant, transcript)

    # Creating DNA level annotation
    try:
        dna, dna_ins = makeDNAannotation(variant, transcript, reference)
    except TypeError:
        dna, dna_ins = 'X', 'X'

    # Creating protein level annotation
    where = transcript.whereIsThisVariant(variant)
    if not '-' in where and where.startswith('Ex'):
        protein, protchange = makeProteinString(variant, transcript, reference, prot, mutprot, coord1)
    else:
        protein, protchange = '', ('.', '.', '.')

    # Transforming coordinates if the variant is a duplication
    if not dna_ins == '':
        coord1_ins, intr1_ins, coord2_ins, intr2_ins = coord1, intr1, coord2, intr2
        coord1, intr1, coord2, intr2 = duplicationCoordinates(variant, transcript)
    else:
        coord1_ins, intr1_ins, coord2_ins, intr2_ins = '', '', '', ''

    # Creating and returning csnAnnot object
    csn = CSNAnnot(coord1, intr1, coord2, intr2, dna, protein, coord1_ins, intr1_ins, coord2_ins, intr2_ins, dna_ins)
    return csn, protchange


# Getting CSN annotation of a given variant
def getAnnotationNew(variant, transcript, reference, prot, mutprot):
    # Creating csn annotation coordinates
    coord1, intr1, coord2, intr2 = calculateCSNCoordinates(variant, transcript)

    # Creating DNA level annotation
    try:
        dna, dna_ins = makeDNAannotation(variant, transcript, reference)
#        dna, dna_ins = makeDNAannotationNew(variant, transcript, reference)
    except TypeError:
        dna, dna_ins = 'X', 'X'

    # Creating protein level annotation
    where = transcript.whereIsThisVariant(variant)
    if not '-' in where and where.startswith('Ex'):
        protein, protchange = makeProteinStringNew(variant, transcript, reference, prot, mutprot, coord1)
    else:
        protein, protchange = '', ('.', '.', '.')

    # Transforming coordinates if the variant is a duplication
    if not dna_ins == '':
        coord1_ins, intr1_ins, coord2_ins, intr2_ins = coord1, intr1, coord2, intr2
        coord1, intr1, coord2, intr2 = duplicationCoordinates(variant, transcript)
    else:
        coord1_ins, intr1_ins, coord2_ins, intr2_ins = '', '', '', ''

    # Creating and returning csnAnnot object
    csn = CSNAnnot(coord1, intr1, coord2, intr2, dna, protein, coord1_ins, intr1_ins, coord2_ins, intr2_ins, dna_ins)
    return csn, protchange


# Calculating csn annotation coordinates
def calculateCSNCoordinates(variant, transcript):
    # Returning coordinates if variant is a base substitution
    if variant.isSubstitution():
        x, y = transformToCSNCoordinate(variant.pos, transcript)
        return x, y, None, None

    # Returning coordinates if variant is an insertion
    if variant.isInsertion():
        startx, starty = transformToCSNCoordinate(variant.pos - 1, transcript)
        endx, endy = transformToCSNCoordinate(variant.pos, transcript)
        if transcript.strand == 1:
            return startx, starty, endx, endy
        else:
            return endx, endy, startx, starty

    # Returning coordinates if variant is a deletion
    if variant.isDeletion():
        startx, starty = transformToCSNCoordinate(variant.pos, transcript)
        if len(variant.ref) == 1:
            return startx, starty, None, None
        endx, endy = transformToCSNCoordinate(variant.pos + len(variant.ref) - 1, transcript)
        if transcript.strand == 1:
            return startx, starty, endx, endy
        else:
            return endx, endy, startx, starty

    # Returning coordinates if variant is a complex indel
    if variant.isComplex():
        startx, starty = transformToCSNCoordinate(variant.pos, transcript)
        if len(variant.ref) == 1:
            return startx, starty, None, None
        endx, endy = transformToCSNCoordinate(variant.pos + len(variant.ref) - 1, transcript)
        if transcript.strand == 1:
            return startx, starty, endx, endy
        else:
            return endx, endy, startx, starty

    return None, None, None, None


# Calculating DNA level annotation of the variant
def makeDNAannotation(variant, transcript, reference):
    # Returning DNA level annotation if variant is a base substitution
    if variant.isSubstitution():
        if transcript.strand == 1:
            return variant.ref + '>' + variant.alt, ''
        else:
            return variant.ref.reverseComplement() + '>' + variant.alt.reverseComplement(), ''

    # Returning DNA level annotation if variant is an insertion
    if variant.isInsertion():
        insert = core.Sequence(variant.alt)
        if transcript.strand == 1:
            before = reference.getReference(variant.chrom, variant.pos - len(insert), variant.pos - 1)
            # Checking if variant is a duplication
            if insert == before and (variant.pos - len(insert) >= transcript.transcriptStart):
                if len(insert) > 4: return 'dup' + str(len(insert)), 'ins' + insert
                return 'dup' + insert, 'ins' + insert
            else:
                return 'ins' + insert, ''
        else:
            before = reference.getReference(variant.chrom, variant.pos, variant.pos + len(insert) - 1)
            # Checking if variant is a duplication
            if insert == before and (variant.pos + len(insert) - 1 <= transcript.transcriptEnd):
                if len(insert) > 4: return 'dup' + str(len(insert)), 'ins' + insert.reverseComplement()
                return 'dup' + insert.reverseComplement(), 'ins' + insert.reverseComplement()
            else:
                return 'ins' + insert.reverseComplement(), ''

    # Returning DNA level annotation if variant is a deletion
    if variant.isDeletion():
        if len(variant.ref) > 4:
            return 'del' + str(len(variant.ref)), ''
        if transcript.strand == 1:
            return 'del' + variant.ref, ''
        else:
            return 'del' + core.Sequence(variant.ref).reverseComplement(), ''

    # Returning DNA level annotation if variant is a complex indel
    if variant.isComplex():
        if transcript.strand == 1:
            return 'delins' + variant.alt, ''
        else:
            return 'delins' + variant.alt.reverseComplement(), ''


# Calculating DNA level annotation of the variant
# this is just a placeholder routine for new functionality.
def makeDNAannotationNew(variant, transcript, reference):
    # Returning DNA level annotation if variant is a base substitution
    if variant.isSubstitution():
        if transcript.strand == 1:
            return variant.ref + '>' + variant.alt, ''
        else:
            return variant.ref.reverseComplement() + '>' + variant.alt.reverseComplement(), ''

# XXX-HS . To support repeat annotation (only applied is (len(variant.alt) %3==0, this code must also be refactored to return the genomic annotation and the cDNA annotation
#      the '[]' notation is to be made at the genomic level... but with shifting that depends on the strand orientation
#      .. also may need to support repeats that are have single sequencing errors... hybrid repeats
# 
    if variant.isInsertion():
        insert = core.Sequence(variant.alt)
        if transcript.strand == 1:
            if len(variant.alt) %3==0:
                a=1 # Placeholder for code to check if there is a repeat at the DNA level.

            before = reference.getReference(variant.chrom, variant.pos - len(insert), variant.pos - 1)
            # Checking if variant is a duplication
            if insert == before and (variant.pos - len(insert) >= transcript.transcriptStart):
                if len(insert) > 4: return 'dup' + str(len(insert)), 'ins' + insert
                return 'dup' + insert, 'ins' + insert
            else:
                return 'ins' + insert, ''
        else:
            if len(variant.alt) %3==0:
                a=1 # Placeholder for code to check if there is a repeat at the DNA level.
            before = reference.getReference(variant.chrom, variant.pos, variant.pos + len(insert) - 1)
            # Checking if variant is a duplication
            if insert == before and (variant.pos + len(insert) - 1 <= transcript.transcriptEnd):
                if len(insert) > 4: 
                    return 'dup' + str(len(insert)), 'ins' + insert.reverseComplement()
                return 'dup' + insert.reverseComplement(), 'ins' + insert.reverseComplement()
            else:
                return 'ins' + insert.reverseComplement(), ''


    # Returning DNA level annotation if variant is a deletion
    if variant.isDeletion():
        if len(variant.ref) %3==0:
            a=1 # Placeholder for code to check if there is a repeat at the DNA level.. Strand-Dependent.
        if len(variant.ref) > 4:
            return 'del' + str(len(variant.ref)), ''
        if transcript.strand == 1:
            return 'del' + variant.ref, ''
        else:
            return 'del' + core.Sequence(variant.ref).reverseComplement(), ''


    # Returning DNA level annotation if variant is a complex indel
    if variant.isComplex():
        if len(variant.alt) %3==0 and len(variant.ref) %3==0:
                a=1 # Placeholder for code to check if there is a repeat at the DNA level. 
                    # This code must be able to handle repeats with sequence variation.. and be strand-dependent.
        if transcript.strand == 1:
            return 'delins' + variant.alt, ''
        else:
            return 'delins' + variant.alt.reverseComplement(), ''



# Calculating protein level annotation of the variant
def makeProteinString(variant, transcript, reference, prot, mutprot, coord1):

    protcopy = prot[:]

    # Checking if there was no change in protein sequence
    if prot == mutprot:
        idx = coord1 / 3
        if coord1 % 3 > 0: idx += 1
        return '_p.=', (str(idx), prot[idx-1], prot[idx-1])

    # Checking if the variant affects the initiating amino acid
    if prot[0] != mutprot[0]: return '_p.' + changeTo3letters(prot[0]) + '1?', ('1', prot[0], mutprot[0])

    # Trimming common starting substring
    leftindex = 1
    rightindex = len(prot)
    while len(prot) > 0 and len(mutprot) > 0:
        if prot[0] == mutprot[0]:
            prot = prot[1:]
            mutprot = mutprot[1:]
            leftindex += 1
        else:
            break

    if prot == '': return '', ('.','.','.')

    # Frameshift mutations
    if (len(variant.alt) - len(variant.ref)) % 3 > 0:
        if len(prot) > 0: return '', (str(leftindex), prot[0], '.')
    else: return '', (str(leftindex), '.', '.')

    # Checking if variant results in a stop lost mutation
    if prot[0] == 'X' and len(mutprot) == 0:
        return '_p.X' + str(leftindex) + '?extX?', (str(leftindex), 'X', '?')

    if prot[0] == 'X' and mutprot[0] != 'X':
        nextstop = mutprot.find('X')
        if nextstop != -1:
            return '_p.X' + str(leftindex) + changeTo3letters(mutprot[0]) + 'extX' + str(nextstop), (str(leftindex), 'X', mutprot[0])
        else:
            return '_p.X' + str(leftindex) + changeTo3letters(mutprot[0]) + 'extX?', (str(leftindex), 'X', mutprot[0])

    # Trimming common ending substring
    while len(prot) > 0 and len(mutprot) > 0:
        if prot[-1] == mutprot[-1]:
            prot = prot[:-1]
            mutprot = mutprot[:-1]
            rightindex -= 1
        else:
            break

    # Checking if variant results in an amino acid change
    if len(prot) == 1 and len(mutprot) == 1: return '_p.' + changeTo3letters(prot) + str(leftindex) + changeTo3letters(
        mutprot), (str(leftindex), prot, mutprot)

    # Checking if variant results in a deletion
    if len(mutprot) == 0:
        if len(prot) == 1:
            return '_p.' + changeTo3letters(prot) + str(leftindex) + 'del', (str(leftindex), prot, '-')
        else:
            if leftindex == rightindex: protpos = str(leftindex)
            else: protpos = str(leftindex) + '-' + str(rightindex)
            return '_p.' + changeTo3letters(prot[0]) + str(leftindex) + '_' + changeTo3letters(prot[-1]) + str(
                rightindex) + 'del',  (protpos, prot, '-')

    # Checking if variant results in an insertion or duplication
    if len(prot) == 0:
        if protcopy[leftindex - len(mutprot) - 1:leftindex - 1] == mutprot:
            if len(mutprot) == 1:
                return '_p.' + changeTo3letters(protcopy[leftindex - len(mutprot) - 1]) + str(
                    leftindex - len(mutprot)) + 'dup',  (str(leftindex - 1)+'-'+str(rightindex + 1), '-', mutprot)
            else:
                return '_p.' + changeTo3letters(protcopy[leftindex - len(mutprot) - 1]) + str(
                    leftindex - len(mutprot)) + '_' + changeTo3letters(protcopy[leftindex - 2]) + str(
                    leftindex - 1) + 'dup',  (str(leftindex - 1)+'-'+str(rightindex + 1), '-', mutprot)
        if 'X' in mutprot: mutprot = mutprot[:mutprot.find('X') + 1]
        return '_p.' + changeTo3letters(protcopy[leftindex - 2]) + str(leftindex - 1) + '_' + changeTo3letters(
            protcopy[rightindex]) + str(rightindex + 1) + 'ins' + changeTo3letters(mutprot), (str(leftindex - 1)+'-'+str(rightindex + 1), '-', mutprot)

    # Checking if variant results in a complex change
    if len(prot) > 0 and len(mutprot) > 0:
        ret = '_p.'
        if len(prot) == 1:
            ret += changeTo3letters(prot) + str(leftindex)
        else:
            ret += changeTo3letters(prot[0]) + str(leftindex) + '_' + changeTo3letters(prot[-1]) + str(rightindex)
        ret += 'delins'
        if 'X' in mutprot: mutprot = mutprot[:mutprot.find('X') + 1]
        ret += changeTo3letters(mutprot)

        if leftindex == rightindex: protpos = str(leftindex)
        else: protpos = str(leftindex) + '-' + str(rightindex)
        return ret, (protpos, prot, mutprot)

    return "", ()


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
def makeProteinStringNew(variant, transcript, reference, prot, mutprot, coord1):
    if prot == '': return '', ('.', '.', '.')
# this is used to distringuish frameshift from non-frameshifts.. though an Early Stop codon .. will be coded as nonsense
#      even if it's a frameshift [HGVS NOTE in fs: the shortest frame shift variantis fsTer2, fsTer1  variants are by definition nonsense variants]
    is_not_frameshift = (len(variant.alt) - len(variant.ref)) % 3 == 0

    xindex=mutprot.find('X')
# Edge case in HGVS, entire protein deleted .. starting 
    if mutprot=='':
        return '_p.Met1?',('1',prot[0],'') 
    if xindex==0:
        return '_p.Met1?',('1',prot[0],'X') 

    # Checking if there was no change in protein sequence
    if prot == mutprot:
        idx = int(coord1 / 3)
        if coord1 % 3 > 0: idx += 1
# Note old CAVA behavior of p.= is incorrect HGVS ... because that means there are no AA change across the entire protein for any variants.
        return '_p.'+changeTo3lettersTer(prot[idx-1])+str(idx)+'=', (str(idx), prot[idx - 1], prot[idx - 1])

    # Checking if the variant affects the initiating amino acid
# ORIGINAL CODE   if prot[0] != mutprot[0]: return '_p.' + changeTo3letters(prot[0]) + '1?', ('1', prot[0], mutprot[0])
# Change to match HGVS requirements.


    # Trimming common starting substring
    # leftindex is a 1-base position within prot
    protcopy = prot[:]
    mutprotcopy = mutprot[:]

    leftindex = 1
    rightindex = len(protcopy)
    while len(prot) > 0 and len(mutprot) > 0:
        if prot[0] == mutprot[0]:
            prot = prot[1:]
            mutprot = mutprot[1:]
            leftindex += 1
        else:
            break

# Any mutation where the first different AA is a Stop is a nonsense mut (wether insertion of deletion or frameshift)
    if mutprot[0]=='X' and prot[0]!='X':
        return '_p.' + changeTo3lettersTer(prot[0]) + str(leftindex) + changeTo3lettersTer(mutprot[0]), (str(leftindex), prot[0], 'X')

    # Trimming common ending substring
    trim_prot=prot
    trim_mutprot=mutprot

    while len(trim_prot) > 0 and len(trim_mutprot) > 0:
        if trim_prot[-1] == trim_mutprot[-1]:
            trim_prot = trim_prot[:-1]
            trim_mutprot = trim_mutprot[:-1]
            rightindex -= 1
        else:
            break

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

    if leftindex==1 and len(prot)>0 and len(mutprot)>0:# leftindex==1 means first base of ref is different
        if prot[0]=='M':
            if len(trim_mutprot)==0: # Deletion of initial part of protein, including Methionine
                if rightindex!=leftindex:
                    return '_p.?', ('1-'+str(rightindex), trim_prot, '-')
                else:
                    return '_p.?', ('1', trim_prot, '-')
            else:
                if rightindex!=leftindex: # del/ins/complex
                    return '_p.?', ('1-'+str(rightindex), trim_prot, trim_mutprot)
                else:
                    return '_p.?', ('1', trim_prot, trim_mutprot)

        else: # Incomplete reference protein without a Methionine at the Start
            if mutprot[0]=='M': # Incomplete reference protein without start becoming a start
               if len(trim_prot)==1:
                    return '_p.' + changeTo3lettersTer(prot[0]) + '1' + changeTo3lettersTer(mutprot[0]) , ('1', trim_prot, trim_mutprot)
               else:
                    return '_p.' + changeTo3lettersTer(prot[0]) + '1' + changeTo3lettersTer(mutprot[0]) , ('1-'+str(len(trim_prot)), trim_prot, trim_mutprot)
            else: 
                if len(trim_prot)==1:
                    return '_p.?', ('1', trim_prot, trim_mutprot)
                else:
                    return '_p.?', ('1-'+str(len(trim_prot)), trim_prot, trim_mutprot)

# Note .. if a deletion is shifted 3' because of a repeat sequence .. it has to be called at this 3' position.
#    The particular scanning algorithm insures that AA-level variants are also shifted right.
#
# Single AA Substitution has Precedence over frameshift, deleting, or insertion or complex (Except for Start Stop Subs)
#               For the case of Ins of fs.. if an early Stop occur.. and is encodable as a substitution

# Any mutation where the first different AA is a new Stop ==> priority as nonsense mut (wether ins, del or fs, or complex )
    if mutprot[0]=='X' and prot[0]!='X':
        return '_p.' + changeTo3lettersTer(prot[0]) + str(leftindex) + changeTo3lettersTer(mutprot[0]), (str(leftindex), prot[0], 'X')

# Checking if the first base altered results in a stop lost mutation
    if prot[0] == 'X': # Don't have to check for frameshift or SSR if stop codon is first mutated.
        if  len(mutprot)==0: # 3' Extension of unknown length on reference sequence without reference sequence that can lead to extension
            return '_p.Ter' + str(leftindex) + '?extX?', (str(leftindex), 'X', '?')
        else: #  by definition or "prot" and left-scanning ==> mutprot[0] != 'X' : # 3' Extension
            nextstop = mutprot.find('X')
# HS: Changed CSN extX to ext*
            if nextstop != -1:
                return '_p.Ter' + str(leftindex) + changeTo3lettersTer(mutprot[0]) + 'extX' + str(nextstop+1), (
                    str(leftindex), 'X', mutprot[0:(nextstop+1)])
            else:
                return '_p.Ter' + str(leftindex) + changeTo3lettersTer(mutprot[0]) + 'ext*?', (str(leftindex), 'X', mutprot)

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
    couldbe_delstop_or_ins_extendingstop=is_not_frameshift and rightindex>=len(protcopy)-1 and leftindex<len(protcopy) and not ( (len(trim_prot)==1 and len(trim_mutprot)==1) or (len(trim_prot)==0 and len(protcopy)>len(mutprotcopy)))


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
        return '_p.' + changeTo3lettersTer(trim_prot) + str(leftindex) + changeTo3lettersTer(trim_mutprot), (str(leftindex), trim_prot, trim_mutprot)
#
# For a Variant to be called a repeat change, the entire event must be describable as a repeat variation.. There is no repeat+1 base change
#
# Pure Deletion Repeat Event (not a frameshift) that does not include the beginning or end (because such highly impactfull changes were dealth earlier)
#

    if (rightindex<len(protcopy) and len(trim_mutprot)==0 and is_not_frameshift):
        nDup=0
# Loop over potential repeat sizes (Smaller repeat sizes have priority)
        for SSRlen in range(1,len(trim_prot)+1):
# Make sure repeat size is a full multiple of the deletion
            nDups=int(len(trim_prot)/SSRlen)
            if nDups>=1 and len(trim_prot) % SSRlen ==0 :
                repeat_seq=trim_prot[0:SSRlen]
                nMatch_del=0
                nMatch_ref=0
# Check if the repeat pattern matches every base of the trim_prot
                for iDup in range(0,nDups):
                    if trim_prot[(0+SSRlen*iDup):(SSRlen*(iDup+1))]==repeat_seq:
                        nMatch_del=nMatch_del+1
                    else:
                        break
                if nMatch_del==nDups: # Can have 0... nMatch_ref copies in reference sequence left after deletion 
                    lowerlim=leftindex - SSRlen - 1
                    upperlim=lowerlim+SSRlen
                    while lowerlim>=0 and protcopy[lowerlim:upperlim]==repeat_seq:
                        nMatch_ref=nMatch_ref+1
                        lowerlim=lowerlim-SSRlen
                        upperlim=lowerlim+SSRlen
                    if not (nMatch_del==1 and nMatch_ref==0 or nMatch_ref+nMatch_del<=1): # Deletion of a single lone copy is not repeat polymorphism 
                        nDup=nDups
                        lowerlim=leftindex-SSRlen*(nMatch_ref)-1
                        upperlim=lowerlim+SSRlen
                        break

        if nDup==0: # No Repeats .. straight Deletion .. not frameshift .. not ins or complex
            if len(trim_prot)==1:
                return '_p.'+changeTo3lettersTer(trim_prot[0])+str(leftindex)+"del", (str(leftindex), trim_prot, '-')
            else:
                return '_p.'+changeTo3lettersTer(trim_prot[0])+str(leftindex)+"_"+changeTo3lettersTer(trim_prot[len(trim_prot)-1])+str(rightindex)+"del", (str(leftindex)+'-'+str(leftindex+len(trim_prot)-1),trim_prot,'-')
        else:
            if SSRlen==1: # Repeat unit of 1 base.
                if nMatch_ref==0: # Multiple >1 copies are deleted, leaving empty sequence. Guaranteed to be an interval of at least 2 AA positions to be called SSR
                    return '_p.'+changeTo3lettersTer(protcopy[lowerlim])+str(lowerlim+1)+"["+str(nMatch_ref+nMatch_del)+"]%3B["+str(nMatch_ref) + "]",(str(leftindex)+'-'+str(rightindex),trim_prot,'-')
                else:
                    if leftindex==rightindex:
                        return '_p.'+changeTo3lettersTer(protcopy[lowerlim])+str(lowerlim+1)+"["+str(nMatch_ref+nMatch_del)+"]%3B["+str(nMatch_ref) + "]",(str(leftindex),trim_prot,'-')
                    else:
                        return '_p.'+changeTo3lettersTer(protcopy[lowerlim])+str(lowerlim+1)+"["+str(nMatch_ref+nMatch_del)+"]%3B["+str(nMatch_ref) + "]",(str(leftindex)+'-'+str(rightindex),trim_prot,'-')

            else: # Multi base repeat
                if nMatch_ref==0: # Deletion of repeats leaving no copies
                    return '_p.'+changeTo3lettersTer(protcopy[lowerlim])+str(lowerlim+1)+"_"+changeTo3lettersTer(protcopy[upperlim-1])+str(upperlim)+"["+str(nMatch_ref+nMatch_del)+"]%3B["+str(nMatch_ref) + "]",(str(leftindex)+"-"+str(rightindex),trim_prot,'-')
                else:
                    return '_p.'+changeTo3lettersTer(protcopy[lowerlim])+str(lowerlim+1)+"_"+changeTo3lettersTer(protcopy[upperlim-1])+str(upperlim)+"["+str(nMatch_ref+nMatch_del)+"]%3B["+str(nMatch_ref) + "]",(str(leftindex)+'-'+str(rightindex),trim_prot,'-')

#
# Pure insertion/Repeat/Dup .. not frameshift Extension
#
# Insertions Repeats (or Dups) require that "a sequence where, compared to a reference sequence, a segment of one or more amino acids (the repeat unit) is present several times, one after the other.."
#   so we cannot have nMatch_ref==0
#   Make sure insertion does not include multiple "Ter"
    xindex=trim_mutprot.find('X')
    if leftindex>1 and rightindex<len(protcopy) and len(trim_prot)==0 and len(trim_mutprot)>0 and is_not_frameshift:
        nDup=0
# Loop over potential repeat sizes (Smaller repeat sizes have priority)
        for SSRlen in range(1,len(trim_mutprot)+1):
# Make sure repeat size is a full multiple of the deletion
            nDups=int(len(trim_mutprot)/SSRlen)
            if nDups>=1 and len(trim_mutprot) % SSRlen ==0 :
                repeat_seq=trim_mutprot[0:SSRlen]
                nMatch_ins=0
                nMatch_ref=0
# Check if the repeat pattern matches every base of the trim_prot
                for iDup in range(0,nDups):
                    if trim_mutprot[(0+SSRlen*iDup):(SSRlen*(iDup+1))]==repeat_seq:
                        nMatch_ins=nMatch_ins+1
                    else:
                        break
                if nMatch_ins==nDups: # Can have 1... nMatch_ref copies in reference sequence (at least 1 copy is requires for insertions)
                    lowerlim=leftindex - SSRlen-1
                    upperlim=lowerlim+SSRlen
#                         Assume that Variant was right shifted to begin with
                    while lowerlim>=0 and protcopy[lowerlim:upperlim]==repeat_seq:
                        nMatch_ref=nMatch_ref+1
                        lowerlim=lowerlim-SSRlen
                        upperlime=lowerlim+SSRlen
                    if nMatch_ref>0: # We require a copy to be present on the reference "protein" .. as per HGVS
                        nDup=nDups
                        lowerlim=leftindex-SSRlen*(nMatch_ref)-1
                        upperlim=lowerlim+SSRlen
                        break   # break out of the for iDup loop
                    
        if xindex!=-1 or nDup==0: # Insertion contain Ter (don't allow repeats with a Ter-containing pattern)
                                  # No Repeats .. straight Insertion .. not frameshift .. not Deletion or complex
                                  # insertion sequence include 1 AA  before.(leftindex is guaranteed to be 2 or more)
                                  # Make sure inserted sequence does not include a Stop codon.
            if xindex==-1:
                return '_p.'+changeTo3lettersTer(protcopy[leftindex-2])+str(leftindex-1)+'_'+changeTo3lettersTer(protcopy[leftindex-1])+str(leftindex)+"ins"+changeTo3lettersTer(trim_mutprot), (str(leftindex-1)+"-"+str(rightindex+1),'-',trim_mutprot)
            else: # Inserted protein includes a Stop (but not as first sequence ... that would have been caught earlier)
                  # also, because rightindex<len(protcopy) .. we know last Ter matches.
                return '_p.'+changeTo3lettersTer(protcopy[leftindex-2])+str(leftindex-1)+'_'+changeTo3lettersTer(protcopy[rightindex])+str(rightindex+1)+"ins"+changeTo3lettersTer(trim_mutprot[0:(xindex+1)]), (str(leftindex-1)+'-'+str(rightindex+1),'-',trim_mutprot[0:(xindex+1)])
        else: # Single and Multi-base repeat insertion
            if nDup==1 and nMatch_ref==1: # Duplications, special treatment
                if SSRlen==1:
                    return '_p.'+changeTo3lettersTer(protcopy[lowerlim])+str(lowerlim+1)+"dup",(str(leftindex-1)+'-'+str(rightindex+1),'-',trim_mutprot)
                else:
                    return '_p.'+changeTo3lettersTer(protcopy[lowerlim])+str(lowerlim+1)+'_'+changeTo3lettersTer(protcopy[upperlim-1])+str(upperlim)+"dup",(str(leftindex-1)+'-'+str(rightindex+1),'-',trim_mutprot)
            else:
                if SSRlen==1:
                    return '_p.'+changeTo3lettersTer(protcopy[lowerlim])+str(lowerlim+1)+"["+str(nMatch_ref)+"]%3B["+str(nMatch_ref+nMatch_ins) + "]",(str(leftindex-1)+'-'+str(rightindex+1),'-',trim_mutprot)
                else:
                    return '_p.'+changeTo3lettersTer(protcopy[lowerlim])+str(lowerlim+1)+'_'+changeTo3lettersTer(protcopy[upperlim-1])+str(upperlim)+"["+str(nMatch_ref)+"]%3B["+str(nMatch_ref+nMatch_ins) + "]",(str(leftindex-1)+'-'+str(rightindex+1),'-',trim_mutprot)
            


# Frameshift mutations (assume len(prot)>0 from now on)
# New behavior for Frameshift , PROTALT will be the sequence (not just a ".")
#
# Also include Some "In-phase deletions or delins:
#       HGVS Deletion definition: "deletions starting N-terminal of and including the translation termination (stop) codon are described as Frame shift. .. "
#  ... irrespective of wether it is an actual frameshift (at the dna level) .. this includes Complex Variants that affect Stop too.
#   
#       HGVS delins definition: deletion-insertion variants starting N-terminal () of and including the translation termination (stop) codon are described as frame shift. 
#       Note Mutated protein can have multiple 'X' (0 or more) after the frameshift location.
    if (len(variant.alt) - len(variant.ref)) % 3 > 0  or couldbe_delstop_or_ins_extendingstop:
        if len(mutprot)==0: # Deletion until the end of the protein and whole UTR (since no extra AA picked up from UTR)
#         Since len(prot)>0, then the variant is a deletion (could NOT be a frameshift that causes an early Stop .. because len(mutprot)==0
#         It is possible that the DNA variant is a complex variant and not a del, but that is OK (protein annotations should not consider DNA).
            if len(prot)==1: # Unless last AA is not Ter .. then this should have been dealth with by Code above
                return '_p.',changeTo3lettersTer(prot[0])+str(leftindex)+"del", (str(leftindex), prot, '-')
            else:
                return '_p.',changeTo3lettersTer(prot[0])+str(leftindex)+"_"+changeTo3lettersTer(prot[len(prot)-1])+str(rightindex)+"del", (str(leftindex)+'-'+str(rightindex), prot, '-')
        else:
            xindex=mutprot.find("X")
            if xindex!=-1: #  p.(Arg123LysfsTer34)
# Special Rules for when frameshift causes the Initial Met to become Stop
                if xindex==0 and leftindex==1 and protcopy[0]=="M" : 
                    return '_p.(Met1Ter)',('1','M','X')
# Check to see if new Stop is a 3' extension.. can only occur if last AA is the first one changed.
                if leftindex==len(protcopy): # last AA of protein is first mutated ==> Extensions
                    if xindex+1>len(prot) : # extension, at or past end of current protein. do not check to see if last original codon is X,x,* .. to allow partial reference.
                        return '_p.'+changeTo3lettersTer(prot[0])+str(len(protcopy))+changeTo3lettersTer(mutprot[0])+'extX'+str(xindex+1),(str(leftindex),prot[0],mutprot[0:(xindex+1)])
                    # else xindex+1<=len(prot) .. only xindex+1==len(prot) possible .. which would mean that the new Ter is at the end of the protein.(only possible is original protein is partial .. considered below)
                if xindex+1== len(prot): # Stop codon in mutprot at end of original protein .. either no change to original stop codon, or original prot did not have a stop codon there.
                    if len(prot)==1: # Extension original prot did not have Stop there.. but mutated one does
                        return '_p.'+changeTo3lettersTer(prot[0])+str(len(protcopy))+changeTo3lettersTer(mutprot[0])+'extX'+str(xindex+1),(str(leftindex),prot[0],mutprot[0:(xindex+1)])
                    # else Frameshift occuring before the end of the protein with the Ter occuring a few AA later at the original prot position (will be dealt below)
                if xindex==0 and len(prot)>1  : #First Modified Base is a stop codon before end (Independent wether mutation is deletion or insertion
                    return '_p.'+changeTo3lettersTer(prot[0])+str(leftindex)+'Ter',(str(leftindex),prot[0],mut[0])
                if leftindex<len(protcopy) and xindex+1<= len(prot): #First Modified Base before the end of the protein .. and new Stop not past original Stop (e.g. not an extension)
                                                                     # (btw Minimal length of prot==2 for protein ending in Ter)
                    if len(trim_prot) == 0: # insertion
#                       Duplications notation has priority over insertion (except if the Stop is within the new allele)
#       Note that because len(trim_mutprot)>0 ==> leftindex>1 is required to be an insertion
#         also note that anything more than a Dup .. was dealt before
                        if len(trim_mutprot)>0  and (leftindex - len(trim_mutprot)) >0 and protcopy[leftindex - len(trim_mutprot) - 1:leftindex - 1] == trim_mutprot:
                            if len(trim_mutprot) == 1:
                                return '_p.' + changeTo3lettersTer(protcopy[leftindex - len(trim_mutprot) - 1]) + str(
                                        leftindex - len(trim_mutprot)) + 'dup', (str(leftindex - 1) + '-' + str(rightindex + 1), '-', trim_mutprot)
                            else:
                                return '_p.' + changeTo3lettersTer(protcopy[leftindex - len(trim_mutprot) - 1]) + str(
                                        leftindex - len(trim_mutprot)) + '_' + changeTo3lettersTer(protcopy[leftindex - 2]) + str(
                                        leftindex - 1) + 'dup', (str(leftindex - 1) + '-' + str(rightindex + 1), '-', trim_mutprot)

                    # insertion or deletion frameshift starting and ending before the original protein end.
                    return '_p.' + changeTo3lettersTer(prot[0])+str(leftindex)+changeTo3lettersTer(mutprot[0])+'fsTer'+str(xindex+1),(str(leftindex),prot[0],mutprot[0:(xindex+1)])
                else: # xindex+1 > len(prot) .. so new protein past end .. and since  extension already returned above, so leftindex<len(protcopy).
                    return '_p.' + changeTo3lettersTer(prot[0])+str(leftindex)+changeTo3lettersTer(mutprot[0])+'fsTer'+str(xindex+1),(str(leftindex),prot[0],mutprot[0:(xindex+1)])
            else: # No Stop codon until the end of the transcript.. Nevertheless Need to supply the protein, so we can test if nonsense Mediated Decay occurs.
                if leftindex==len(protcopy) : # Extension without predicted Ter
                    return '_p.(' + changeTo3lettersTer(prot[0])+str(leftindex)+changeTo3lettersTer(mutprot[0])+'ext*?', (str(leftindex),prot[0],mutprot)
                else: 
                    return '_p.(' + changeTo3lettersTer(prot[0])+str(leftindex)+changeTo3lettersTer(mutprot[0])+'fs*?)',(str(leftindex),prot[0],mutprot)


# Past this point, only Complex delins not including the stop codon, 5' extension insertions



    # Trimming common ending substring
    prot=trim_prot
    mutprot=trim_mutprot

    prot=trim_prot
    mutprot=trim_mutprot

    #  (len(prot)==0 and len(mutprot)==0) should have been covered by  earlier cases

    # deletion at the ends (Already checked for "Internal deletions" and deleted first base Methionine earlier)
    if len(trim_mutprot) == 0: 
        if len(trim_prot) == 1: 
            return '_p.' + changeTo3lettersTer(trim_prot) + str(leftindex) + 'del', (str(leftindex), trim_prot, '-')
        else:
            protpos = str(leftindex) + '-' + str(rightindex)
            return '_p.' + changeTo3lettersTer(trim_prot[0]) + str(leftindex) + '_' + changeTo3lettersTer(trim_prot[-1]) + str(
                rightindex) + 'del', (protpos, trim_prot, '-')

# Insertion-like, but affecting first or last.
    if len(prot) == 0:
# New functionality (includes some 5' extension) when leftindex==1 and len(prot)==0
        if leftindex==1: # Inserted sequences before the Original Methionine
            if trim_mutprot[0]=='M': # Extension 
                if protcopy[0]=='M': # insertion e.g. MA-> MRMA : M1_A2insRM (prot="" -> trim_mutprot=MR)  .. Only this one is realistic.
                    return '_p.' + changeTo3lettersTer(protcopy[0]) + '1_' + changeTo3lettersTer(
                        protcopy[1]) + '2ins' + changeTo3lettersTer(trim_mutprot[1:])+changeTo3lettersTer('M'), (
                            '1-2', 'M', trim_mutprot+'M')
                else: # e.g. A->MRA (prot="" --> trim_mutprot="MR")
                    return '_p.?',('1',protcopy[0],trim_mutprot+protcopy[0])
            else: # MAR->YMAR (""->A) or AR->GAR (""->G)
                return '_p.?',('1','-',trim_mutprot)
        if leftindex==len(protcopy): # Only case is if the original protein did not have a Ter e.g. ARG->ARGYX -> prot="", trim_mutprot="YX"
# Not covered by the HGVS standard, unless it's a dup .. and that was covered above
            return '_p.?',(leftindex,'-',trim_mutprot)

    # Checking if variant results in a complex change
    if len(trim_prot) > 0 and len(trim_mutprot) > 0:
        xindex=trim_mutprot.find('X')
        if xindex!=-1:
            trim_mutprot=trim_mutprot[0:(xindex+1)]
        ret = '_p.'
        if len(trim_prot) == 1:
            ret += changeTo3lettersTer(trim_prot) + str(leftindex)
        else:
            ret += changeTo3lettersTer(trim_prot[0]) + str(leftindex) + '_' + changeTo3lettersTer(trim_prot[-1]) + str(rightindex)
        ret += 'delins'
        
        ret += changeTo3lettersTer(trim_mutprot)

        if leftindex == rightindex:
            protpos = str(leftindex)
        else:
            protpos = str(leftindex) + '-' + str(rightindex)
        return ret, (protpos, trim_prot, trim_mutprot)
    sys.stderr.write("\nBUG: Cannot compute CSN for : original prot="+protcopy+"\n   : Cannot compute CSN for : mutated  prot="+mutprotcopy+"\n")


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
                        if pos <= int((exon.start + 1 - prevExonEnd) / 2) + prevExonEnd:
                            x, y = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, pos - prevExonEnd
                        else:
                            x, y = transformToCSNCoordinate(exon.start + 1, transcript)
                            return x, pos - exon.start - 1
                else:
                    # Checking if genomic position is within intron
                    if exon.end < pos < prevExonEnd:
                        if pos >= int((prevExonEnd - exon.end + 1) / 2) + exon.end:
                            x, y = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, prevExonEnd - pos
                        else:
                            x, y = transformToCSNCoordinate(exon.end, transcript)
                            return x, exon.end - pos
            # Checking if genomic position is within exon
            if exon.start + 1 <= pos <= exon.end:
                if transcript.strand == 1:
                    return sumOfExonLengths + pos - exon.start, 0
                else:
                    return sumOfExonLengths + exon.end - pos + 1, 0
            # Calculating sum of exon lengths up to this point
            sumOfExonLengths += exon.length
            if transcript.strand == 1:
                prevExonEnd = exon.end
            else:
                prevExonEnd = exon.start + 1

    # If genomic position is outside CDS
    else:
        sumpos = 0
        for i in range(len(transcript.exons)):
            exon = transcript.exons[i]

            if i > 0:
                if transcript.strand == 1:
                    # Checking if genomic position is within intron
                    if prevExonEnd < pos < exon.start + 1:
                        if pos <= int((exon.start + 1 - prevExonEnd) / 2) + prevExonEnd:
                            x, y = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, pos - prevExonEnd
                        else:
                            x, y = transformToCSNCoordinate(exon.start + 1, transcript)
                            return x, pos - exon.start - 1
                else:
                    # Checking if genomic position is within intron
                    if exon.end < pos < prevExonEnd:
                        if pos >= int((prevExonEnd - exon.end + 1) / 2) + exon.end:
                            x, y = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, prevExonEnd - pos
                        else:
                            x, y = transformToCSNCoordinate(exon.end, transcript)
                            return x, exon.end - pos

            if transcript.strand == 1:
                if pos > transcript.codingEndGenomic:
                    if transcript.codingEndGenomic < exon.start + 1 and exon.end < pos:
                        sumpos += exon.length
                    elif exon.contains(transcript.codingEndGenomic) and exon.end < pos:
                        sumpos += exon.end - transcript.codingEndGenomic + 1
                    elif exon.contains(transcript.codingEndGenomic) and exon.contains(pos):
                        sumpos += pos - transcript.codingEndGenomic
                    elif transcript.codingEndGenomic < exon.start + 1 and exon.contains(pos):
                        sumpos += pos - exon.start - 1
                if pos < transcript.codingStartGenomic:
                    if pos < exon.start + 1 and exon.end < transcript.codingStartGenomic:
                        sumpos += exon.length
                    elif pos < exon.start + 1 and exon.contains(transcript.codingStartGenomic):
                        sumpos += transcript.codingStartGenomic - exon.start
                    elif exon.contains(pos) and exon.contains(transcript.codingStartGenomic):
                        sumpos += transcript.codingStartGenomic - pos
                    elif exon.contains(pos) and exon.end < transcript.codingStartGenomic:
                        sumpos += exon.end - pos

            if transcript.strand == -1:
                if pos < transcript.codingEndGenomic:
                    if pos < exon.start + 1 and exon.end < transcript.codingEndGenomic:
                        sumpos += exon.length
                    elif pos < exon.start + 1 and exon.contains(transcript.codingEndGenomic):
                        sumpos += transcript.codingEndGenomic - exon.start
                    elif exon.contains(pos) and exon.contains(transcript.codingEndGenomic):
                        sumpos += transcript.codingEndGenomic - pos
                    elif exon.contains(pos) and exon.end < transcript.codingEndGenomic:
                        sumpos += exon.end - pos
                if transcript.codingStartGenomic < pos:
                    if transcript.codingStartGenomic < exon.start + 1 and exon.end < pos:
                        sumpos += exon.length
                    elif exon.contains(transcript.codingStartGenomic) and exon.end < pos:
                        sumpos += exon.end - transcript.codingStartGenomic + 1
                    elif exon.contains(transcript.codingStartGenomic) and exon.contains(pos):
                        sumpos += pos - transcript.codingStartGenomic
                    elif transcript.codingStartGenomic < exon.start + 1 and exon.contains(pos):
                        sumpos += pos - exon.start - 1

            if transcript.strand == 1:
                prevExonEnd = exon.end
            else:
                prevExonEnd = exon.start + 1

        if transcript.strand == 1:
            if pos > transcript.codingEndGenomic:  return '+' + str(sumpos), 0
            if pos < transcript.codingStartGenomic: return '-' + str(sumpos), 0
        else:
            if pos < transcript.codingEndGenomic: return '+' + str(sumpos), 0
            if pos > transcript.codingStartGenomic: return '-' + str(sumpos), 0

        return str(sumpos), 0


# Calculating csn coordinates for duplications
def duplicationCoordinates(variant, transcript):
    if transcript.strand == 1:
        coord1, intr1 = transformToCSNCoordinate(variant.pos - len(variant.alt), transcript)
        if len(variant.alt) == 1:
            coord2, intr2 = None, None
        else:
            coord2, intr2 = transformToCSNCoordinate(variant.pos - 1, transcript)
    else:
        coord1, intr1 = transformToCSNCoordinate(variant.pos + len(variant.alt) - 1, transcript)
        if len(variant.alt) == 1:
            coord2, intr2 = None, None
        else:
            coord2, intr2 = transformToCSNCoordinate(variant.pos, transcript)
    return coord1, intr1, coord2, intr2


# Changing protein sequence of 1-letter amino acid code to 3-letter code
def changeTo3letters(aas):
    ret = ''
    codes = {
        'I': 'Ile', 'M': 'Met', 'T': 'Thr', 'N': 'Asn',
        'K': 'Lys', 'S': 'Ser', 'R': 'Arg', 'L': 'Leu',
        'P': 'Pro', 'H': 'His', 'Q': 'Gln', 'V': 'Val',
        'A': 'Ala', 'D': 'Asp', 'E': 'Glu', 'G': 'Gly',
        'F': 'Phe', 'Y': 'Tyr', 'C': 'Cys', 'W': 'Trp',
        '*': 'X', 'X': 'X', 'x': 'X', '?': '?'}
    for aa in aas: ret += codes[aa]
    return ret


# Changing protein sequence of 1-letter amino acid code to 3-letter code, but with Stops changed to "Ter"
def changeTo3lettersTer(aas):
    ret = ''
    codes = {
        'I': 'Ile', 'M': 'Met', 'T': 'Thr', 'N': 'Asn',
        'K': 'Lys', 'S': 'Ser', 'R': 'Arg', 'L': 'Leu',
        'P': 'Pro', 'H': 'His', 'Q': 'Gln', 'V': 'Val',
        'A': 'Ala', 'D': 'Asp', 'E': 'Glu', 'G': 'Gly',
        'F': 'Phe', 'Y': 'Tyr', 'C': 'Cys', 'W': 'Trp',
        '*': 'Ter', 'X': 'Ter', 'x': 'Ter'}
    for aa in aas: ret += codes[aa]
    return ret

#######################################################################################################################
