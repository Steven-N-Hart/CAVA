#!/usr/bin/env python3


# CLASS annotation
#######################################################################################################################

# Getting CLASS annotation of a given variant
import sys
def getClassAnnotation(variant, transcript, protein, mutprotein, loc, ssrange):
    # Variants in UTR, if "Ex" present, variant overlaps CDS
    chkutr = checkUTR(transcript, variant)
    if chkutr == 'UTR5': # Even if a variant covers more than UTR5, UTR5 has priority. .. also includes UTR5-Ex* .. or Ex*-UTR3
        if 'Ex' in loc and (variant.is_deletion or variant.is_complex): # "Ex" actually means coding region
            return 'IM'
        return '5PU'
    elif chkutr == 'UTR3':
        if 'Ex' in loc and (variant.is_deletion or variant.is_complex):
            return 'SL'
        return '3PU'

    # Variants crossing exon-intron boundaries
    if '-' in loc:
        return 'ESS'

    # Intronic, splice site and essential splice site variants
    if 'In' in loc:
        if transcript.isInEssentialSpliceSite(variant): return 'ESS'
        if transcript.intronLength(int(loc[loc.find('/') + 1:])) >= 9:
            if transcript.isIn_SS5_Site(variant): return 'SS5'
        if transcript.isInSplicingRegion(variant, ssrange): return 'SS'
        return 'INT'

    # Checking if variant affect the first or last three bases of an exon
    potSS = transcript.isInFirstOrLast3BaseOfExon(variant)

    if mutprotein is None:
        return 'INT'
    protL = len(protein)
    mutprotL = len(mutprotein)

    # Synonymous coding variants
    if protein == mutprotein:
        if potSS:
            return 'EE'
        return 'SY'

    # Variants affecting the initiation amino acid
    if (mutprotL>0 and protL>0 and protein[0] != mutprotein[0]):
        return 'IM'
    if mutprotL==0 and protL>0: # if it was a SL or IM, it should have been caught earlier. This is a variant overlapping the ends
        return ''

    # Stop gain and stop lost variants
    # XXX-HS rewrote because this was the 3rd sink of time.
    # Remove AA that are identical between prot and mutprot (beginning)
    isame = -1
    while protL > isame+1 and mutprotL > isame+1:
        if protein[isame+1] == mutprotein[isame+1]:
            isame = isame + 1
        else:
            break
    if isame>=0:
        protein = protein[isame+1:]
        mutprotein = mutprotein[isame+1:]

    # By this point, variant is guaranteed to be inside CDS
    if protein == '':  #mutprotein cannot be len(0) .. would have returnes SY/potSS.. so this must be an extension
        # Can only get here if CDS has no stop codon on reference (annotation issue)
        sys.stderr.write("ERROR:Inconsistent CAVA, looks like CDS annotation does not end in a Stop codon for "+variant.id+"\n")
        return 'SL'

    if protein[0] == 'X' and len(mutprotein) == 0: return 'SL'
    if protein[0] == 'X' and mutprotein[0] != 'X': return 'SL'

    if len(mutprotein) > 0 and  mutprotein[0] == 'X' and \
            len(protein) > 0 and \
            mutprotein[0] != protein[0] : # as Per HGVS SG , since shortest FS is Ter2 (7-26-2022)
        return 'SG'

    # Frame-shift coding variants
    if variant.is_in_frame is False:
        return 'FS'

    # XXX-HS rewrote this section because it was 3rd slowest
    # Trim the end of the protein and mutprotein that are identical
    # At this point, any change is either a point mutation of insertion/deletion of multiple AA
    ilast  = 0
    while ilast<len(protein) and ilast<len(mutprotein) :
        if protein[-(ilast+1)] == mutprotein[-(ilast+1)]:
            ilast = ilast +1
        else:
            break
    if ilast >0:
        protein = protein[:-ilast]
        mutprotein = mutprotein[:-ilast]

    #if 'X' in mutprotein:  # In frame insertion of a Stop Codo, but not first AA change
    # According to HGVS, this is described as an IF (deletion/insertion)
    #    return 'SG'

    # Non-synonymous coding variants
    if protL == mutprotL and len(protein)==1:
        if potSS:
            return 'EE'
        return 'NSY'
    # IF variants from here on.
    if potSS:
        return 'EE'
    # In-frame variants
    if potSS:
        return 'EE'
    return 'IF'


#######################################################################################################################


# Sequence Ontology (SO)
#######################################################################################################################


def getSequenceOntologyAnnotation(variant, transcript, protein, mutprotein, loc):

    # Variants in UTR
    chkutr = checkUTR(transcript, variant)
    if chkutr == 'UTR5':
        if 'Ex' in loc and (variant.is_deletion or variant.is_complex): return 'initiator_codon_variant'
        return '5_prime_UTR_variant'
    elif chkutr == 'UTR3':
        if 'Ex' in loc and (variant.is_deletion or variant.is_complex): return 'stop_lost'
        return '3_prime_UTR_variant'

    where = transcript.whereIsThisVariant(variant)

    if '-' in where:
        first = where[:where.find('-')]
        if first.startswith('In'): return 'splice_acceptor_variant'
        if first.startswith('Ex'): return 'splice_donor_variant'
        if first.startswith('fsIn'): return 'splice_acceptor_variant'

    if 'In' in where:

        if isInSpliceDonor(transcript, variant): return 'splice_donor_variant'
        if isInSpliceAcceptor(transcript, variant): return 'splice_acceptor_variant'

        if transcript.intronLength(int(where[where.find('/') + 1:])) >= 9:
            if transcript.isIn_SS5_Site(variant): return 'splice_donor_5th_base_variant'

    out = []

    if isInSplicingRegion(transcript, variant):
        if where.startswith('In') or where.startswith('fsIn'):
            return 'intron_variant|splice_region_variant'
        else:
            out.append('splice_region_variant')

    if where.startswith('In') or where.startswith('fsIn'): return 'intron_variant'

    if variant.is_in_frame is True:
        if variant.is_deletion:
            out.append('inframe_deletion')
        if variant.is_insertion:
            out.append('inframe_insertion')
        if variant.is_complex:
            out.append('inframe_indel')
    else:
        out.append('frameshift_variant')
        return '|'.join(out)

    if mutprotein is None:
        return '.'
    if len(protein) ==0: # both ends of the variant completely outside transcript.
        return '.'
    if len(protein)>0 and len(mutprotein)==0: # This can only happen when beginning of protein is deleted
        return 'initiator_codon_variant'

    if protein == mutprotein:
        out.append('synonymous_variant')
        return '|'.join(out)

    if protein[0] != mutprotein[0]:
        out.append('initiator_codon_variant')
        return '|'.join(out)

    if (not protein == mutprotein) and len(protein) == len(mutprotein):
        out.append('missense_variant')

    # XXX-HS rewrote because this was the 3rd sink of time.
    isame = -1
    while len(protein) > isame+1 and len(mutprotein) > isame+1:
        if protein[isame+1] == mutprotein[isame+1]:
            isame = isame + 1
        else:
            break
    if isame>=0:
        protein = protein[isame+1:]
        mutprotein = mutprotein[isame+1:]

    if protein == '': return '3_prime_UTR_variant'

    if protein[0] == 'X' and len(mutprotein) == 0:
        out.append('stop_lost')
    else:
        if protein[0] == 'X' and mutprotein[0] != 'X': out.append('stop_lost')

# XXX rewrote because it was among top slowest parts of CAVA
#


    ilast  = 0
    while ilast<len(protein) and ilast<len(mutprotein) :
        if protein[-(ilast+1)] == mutprotein[-(ilast+1)]:
            ilast = ilast +1
        else:
            break
    if ilast >0:
        protein = protein[:-ilast]
        mutprotein = mutprotein[:-ilast]

    if 'X' in mutprotein: out.append('stop_gained')

    if ('stop_gained' in out) or ('stop_lost' in out):
        if 'missense_variant' in out: out.remove('missense_variant')

    return '|'.join(out)


def isInSpliceDonor(transcript, variant):
    if transcript.strand == 1:
        for exon in transcript.exons:
            isLastExon = (exon.index == len(transcript.exons))
            if not isLastExon and variant.overlap(exon.end + 1, exon.end + 2):
                return True
        return False
    else:
        for exon in transcript.exons:
            isLastExon = (exon.index == len(transcript.exons))
            if not isLastExon and variant.overlap(exon.start - 1, exon.start):
                return True
        return False


def isInSpliceAcceptor(transcript, variant):
    if transcript.strand == 1:
        for exon in transcript.exons:
            isFirstExon = (exon.index == 1)
            if not isFirstExon and variant.overlap(exon.start - 1, exon.start):
                return True
        return False
    else:
        for exon in transcript.exons:
            isFirstExon = (exon.index == 1)
            if not isFirstExon and variant.overlap(exon.end + 1, exon.end + 2):
                return True

        return False


def isInSplicingRegion(transcript, variant):
    if transcript.strand == 1:
        for exon in transcript.exons:
            isFirstExon = (exon.index == 1)
            isLastExon = (exon.index == len(transcript.exons))
            if not isLastExon and variant.overlap(exon.end - 2, exon.end + 8): return True
            if not isFirstExon and variant.overlap(exon.start - 7, exon.start + 3): return True
        return False
    else:
        for exon in transcript.exons:
            isFirstExon = (exon.index == 1)
            isLastExon = (exon.index == len(transcript.exons))
            if not isLastExon and variant.overlap(exon.start - 7, exon.start + 3): return True
            if not isFirstExon and variant.overlap(exon.end - 2, exon.end + 8): return True
        return False


#######################################################################################################################
#
# returns UTR3/5 if at least one end of the variant is in the UTR
#
def checkUTR(transcript, variant):
    if variant.is_insertion:
        x = variant.pos - 1
        y = variant.pos
    else:
        x = variant.pos
        y = variant.pos + len(variant.ref) - 1
    if transcript.isPositionOutsideCDS_5prime(x) or transcript.isPositionOutsideCDS_5prime(y): return 'UTR5'
    if transcript.isPositionOutsideCDS_3prime(x) or transcript.isPositionOutsideCDS_3prime(y): return 'UTR3'
    return ''
