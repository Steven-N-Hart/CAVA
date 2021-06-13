import unittest

from cava_.core import Variant
from cava_.csn import makeProteinStringHGVS as makeProteinString


class TestmakeProteinString(unittest.TestCase):

    def test_makeProteinString_emptyProt(self):
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "", "MLX", 1)
        expected = ('', ('.', '.', '.'))
        print("Testing prot= empty\n")
        self.assertEqual(actual, expected)

    def test_makeProteinString_Syn(self):
        variant = Variant("chr1", 1000, "C", "T")
        #        print("Testing Syn Change\n")
        actual = makeProteinString(variant, "NM", "chr1", "MRX", "MLX", 6)
        expected = ('_p.Arg2Leu', ('2', 'R', 'L'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_EarlyStopInMiddleInPhase(self):
        variant = Variant("chr1", 1000, "C", "T")
        print("Early Stop - 3rd pos of 6\n")
        actual = makeProteinString(variant, "NM", "chr1", "MYLRGX", "MYX", 9)
        expected = ('_p.Leu3Ter', ('3', 'L', 'X'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_EarlyStopInMiddleOutofPhase(self):
        variant = Variant("chr1", 1000, "C", "T")
        print("Early Stop - 3rd pos of 6\n")
        actual = makeProteinString(variant, "NM", "chr1", "MYLRGX", "MYX", 8)
        expected = ('_p.Leu3Ter', ('3', 'L', 'X'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_EarlyStopNearEnd(self):
        variant = Variant("chr1", 1000, "C", "T")
        print("Early Stop - 2nd pos of 3\n")
        actual = makeProteinString(variant, "NM", "chr1", "MRX", "MX", 4)
        expected = ('_p.Arg2Ter', ('2', 'R', 'X'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_delMet1(self):
        print("Testing deletion of Initial Methionine")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRX", "LRX", 1)
        expected = ('_p.?', ('1', 'M', '-'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_del2(self):
        print("Testing deletion 2AA")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRYX", "MYX", 1)
        expected = ('_p.Leu2_Arg3del', ('2-3', 'LR', '-'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_del1(self):
        print("Testing deletion of 1AA")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRX", "MRX", 1)
        expected = ('_p.Leu2del', ('2', 'L', '-'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_checkMetIns4(self):
        print("Testing checking insertion after Met")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRX", "MQYMLRX", 3)
        expected = ('_p.Met1_Leu2insGlnTyrMet', ('1-2', '-', 'QYM'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_checkMetIns1(self):
        print("Testing checking 1 AA insertion after Met")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRX", "MQLRX", 3)
        expected = ('_p.Met1_Leu2insGln', ('1-2', '-', 'Q'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_checkInsTer(self):
        print("Testing checking 3AA Insertion that contains a Ter in the middle")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRQX", "MLLXLRQX", 3)
        expected = ('_p.Leu2_Arg3insLeuTer', ('2-3', '-', 'LX'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_checkShortExt(self):
        print("Testing checking extension")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRX", "MLRQX", 12)
        expected = ('_p.Ter4GlnextX2', ('4', 'X', 'QX'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_checkLongExtInPhase(self):
        print("Testing checking long extension in phase")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRX", "MLRQLVYX", 12)
        expected = ('_p.Ter4GlnextX5', ('4', 'X', 'QLVYX'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_checkLongExtnotPhase(self):
        print("Testing checking long extension not in phase")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRX", "MLRQLVYX", 11)
        expected = ('_p.Ter4GlnextX5', ('4', 'X', 'QLVYX'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_checkShortFS(self):
        print("Testing checking frameshift with mutprot<prot")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRYQVRX", "MLSVX", 8)
        expected = ('_p.Arg3SerfsTer3', ('3', 'R', 'SVX'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_checkDellookLikeShortFS(self):
        print("Testing checking Deletion")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRYQVQX", "MLQX", 7)
        expected = ('_p.Arg3_Val6del', ('3-6', 'RYQV', '-'))
        self.assertEqual(actual, expected)

    # Deletions near 3' end are deemed to include the Stop codon
    def test_makeProteinString_checkDelasFS(self):
        print("Testing checking check frameshift")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRYQVRX", "MLRVX", 6)
        expected = ('_p.Tyr4ValfsTer2', ('4', 'Y', 'VX'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_checkLongFS(self):
        print("Testing checking long frameshift")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRYQVQX", "MLQYLVMSNX", 7)
        expected = ('_p.Arg3GlnfsTer8', ('3', 'R', 'QYLVMSNX'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_delins(self):
        print("Testing checking delins")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRYQVQX", "MLRLVISVQX", 7)
        expected = ('_p.Tyr4_Gln5delinsLeuValIleSer', ('4-5', 'YQ', 'LVIS'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_delinswithTerAsFS(self):
        print("Testing checking delins 1 AA becoming 3 (with Ter in middle) ")
        variant = Variant("chr1", 1000, "C", "T")
        # Q --> LXR
        actual = makeProteinString(variant, "NM", "chr1", "MLRYQVQX", "MLRYLXRVQX", 7)
        expected = ('_p.Gln5delinsLeuTer', ('5', 'Q', 'LX'))
        self.assertEqual(actual, expected)

    # deletion-insertion variants starting N-terminal () of and including the translation termination (stop) codon are described as frame shift.

    def test_makeProteinString_delinsasFS(self):
        print("Testing checking delins affecting Stop Should be an FS ")
        variant = Variant("chr1", 1000, "C", "T")
        # Q --> LXR
        actual = makeProteinString(variant, "NM", "chr1", "MLRYQX", "MLRLRX", 7)
        expected = ('_p.Tyr4LeufsTer3', ('4', 'Y', 'LRX'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_extNoTer(self):
        print("Testing checking delins")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRX", "MLRLVI", 7)
        expected = ('_p.Ter4Leuext*?', ('4', 'X', 'LVI'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_dup1(self):
        print("Testing checking dup of 1 base")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRX", "MLLRX", 6)
        expected = ('_p.Leu2dup', ('2-3', '-', 'L'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_dup_AA(self):
        print("Testing checking dup of pattern 2 AA long")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLYRX", "MLYLYRX", 6)
        expected = ('_p.Leu2_Tyr3dup', ('3-4', '-', 'LY'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_ssr2_gain1to3(self):
        print("Testing checking ssr 1 to 3 copies of 2-long SSR")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLYRX", "MLYLYLYRX", 6)
        expected = ('_p.Leu2_Tyr3[1]%3B[3]', ('3-4', '-', 'LYLY'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_ssr2_gain0to3(self):
        print("Testing checking ssr 0 to 3 copies of 2-long SSR .. should not be an SSR")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MRX", "MLYLYLYRX", 6)
        expected = ('_p.Met1_Arg2insLeuTyrLeuTyrLeuTyr', ('1-2', '-', 'LYLYLY'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_ssr2_loss2to1(self):
        print("Testing checking ssr 2-long deletion")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLYLYRX", "MLYRX", 6)
        expected = ('_p.Leu2_Tyr3[2]%3B[1]', ('4-5', 'LY', '-'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_ssr2_loss2to0(self):
        print("Testing checking ssr 2-long deletion to 0")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLYLYRX", "MRX", 6)
        expected = ('_p.Leu2_Tyr3[2]%3B[0]', ('2-5', 'LYLY', '-'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_ssr1_gain1to3(self):
        print("Testing checking 1 to 3 copies of 1-long SSR")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLRX", "MLLLRX", 6)
        expected = ('_p.Leu2[1]%3B[3]', ('2-3', '-', 'LL'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_ssr1_loss2to1(self):
        print("Testing checking ssr 2-long deletion")
        variant = Variant("chr1", 1000, "C", "T")
        actual = makeProteinString(variant, "NM", "chr1", "MLLRX", "MLRX", 6)
        expected = ('_p.Leu2[2]%3B[1]', ('3', 'L', '-'))
        self.assertEqual(actual, expected)

    def test_makeProteinString_basic(self):
        print("Testing checking basic change")
        variant = Variant("chr17", 43045712, "T", "C")
        REF = "MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQETSLSVQLSNLGTVRTLRTKQRIQPQKTSVYIELGSDSSEDTVNKATYCSVGDQELLQITPQGTRDEISLDSAKKAACEFSETDVTNTEHHQPSNNDLNTTEKRAAERHPEKYQGSSVSNLHVEPCGTNTHASSLQHENSSLLLTKDRMNVEKAEFCNKSKQPGLARSQHNRWAGSKETCNDRRTPSTEKKVDLNADPLCERKEWNKQKLPCSENPRDTEDVPWITLNSSIQKVNEWFSRSDELLGSDDSHDGESESNAKVADVLDVLNEVDEYSGSSEKIDLLASDPHEALICKSERVHSKSVESNIEDKIFGKTYRKKASLPNLSHVTENLIIGAFVTEPQIIQERPLTNKLKRKRRPTSGLHPEDFIKKADLAVQKTPEMINQGTNQTEQNGQVMNITNSGHENKTKGDSIQNEKNPNPIESLEKESAFKTKAEPISSSISNMELELNIHNSKAPKKNRLRRKSSTRHIHALELVVSRNLSPPNCTELQIDSCSSSEEIKKKKYNQMPVRHSRNLQLMEGKEPATGAKKSNKPNEQTSKRHDSDTFPELKLTNAPGSFTKCSNTSELKEFVNPSLPREEKEEKLETVKVSNNAEDPKDLMLSGERVLQTERSVESSSISLVPGTDYGTQESISLLEVSTLGKAKTEPNKCVSQCAAFENPKGLIHGCSKDNRNDTEGFKYPLGHEVNHSRETSIEMEESELDAQYLQNTFKVSKRQSFAPFSNPGNAEEECATFSAHSGSLKKQSPKVTFECEQKEENQGKNESNIKPVQTVNITAGFPVVGQKDKPVDNAKCSIKGGSRFCLSSQFRGNETGLITPNKHGLLQNPYRPPLFPIKSFVKTKCKKNLLEENFEEHSMSPEREMGNENIPSTVSTISRNNIRENVFKEASSSNINEVGSSTNEVGSSINEIGSSDENIQAELGRNRGPKLNAMLRLGVLQPEVYKQSLPGSNCKHPEIKKQEYEEVVQTVNTDFSPYLISDNLEQPMGSSHASQVCSETPDDLLDDGEIKEDTSFAENDIKESSAVFSKSVQKGELSRSPSPFTHTHLAQGYRRGAKKLESSEENLSSEDEELPCFQHLLFGKVNNIPSQSTRHSTVATECLSKNTEENLLSLKNSLNDCSNQVILAKASQEHHLSEETKCSASLFSSQCSELEDLTANTNTQDPFLIGSSKQMRHQSESQGVGLSDKELVSDDEERGTGLEENNQEEQSMDSNLGEAASGCESETSVSEDCSGLSSQSDILTTQQRDTMQHNLIKLQQEMAELEAVLEQHGSQPSNSYPSIISDSSALEDLRNPEQSTSEKAVLTSQKSSEYPISQNPEGLSADKFEVSADSSTSKNKEPGVERSSPSKCPSLDDRWYMHSCSGSLQNRNYPSQEELIKVVDVEEQQLEESGPHDLTETSYLPRQDLEGTPYLESGISLFSDDPESDPSEDRAPESARVGNIPSSTSALKVPQLKVAESAQSPAAAHTTDTAGYNAMEESVSREKPELTASTERVNKRMSMVVSGLTPEEFMLVYKFARKHHITLTNLITEETTHVVMKTDAEFVCERTLKYFLGIAGGKWVVSYFWVTQSIKERKMLNEHDFEVRGDVVNGRNHQGPKRARESQDRKIFRGLEICCYGPFTNMPTDQLEWMVQLCGASVVKELSSFTLGTGVHPIVVVQPDAWTEDNGFHAIGQMCEAPVVTREWVLDSVALYQCQELDTYLIPQIPHSHY"
        ALT = "MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQETSLSVQLSNLGTVRTLRTKQRIQPQKTSVYIELGSDSSEDTVNKATYCSVGDQELLQITPQGTRDEISLDSAKKAACEFSETDVTNTEHHQPSNNDLNTTEKRAAERHPEKYQGSSVSNLHVEPCGTNTHASSLQHENSSLLLTKDRMNVEKAEFCNKSKQPGLARSQHNRWAGSKETCNDRRTPSTEKKVDLNADPLCERKEWNKQKLPCSENPRDTEDVPWITLNSSIQKVNEWFSRSDELLGSDDSHDGESESNAKVADVLDVLNEVDEYSGSSEKIDLLASDPHEALICKSERVHSKSVESNIEDKIFGKTYRKKASLPNLSHVTENLIIGAFVTEPQIIQERPLTNKLKRKRRPTSGLHPEDFIKKADLAVQKTPEMINQGTNQTEQNGQVMNITNSGHENKTKGDSIQNEKNPNPIESLEKESAFKTKAEPISSSISNMELELNIHNSKAPKKNRLRRKSSTRHIHALELVVSRNLSPPNCTELQIDSCSSSEEIKKKKYNQMPVRHSRNLQLMEGKEPATGAKKSNKPNEQTSKRHDSDTFPELKLTNAPGSFTKCSNTSELKEFVNPSLPREEKEEKLETVKVSNNAEDPKDLMLSGERVLQTERSVESSSISLVPGTDYGTQESISLLEVSTLGKAKTEPNKCVSQCAAFENPKGLIHGCSKDNRNDTEGFKYPLGHEVNHSRETSIEMEESELDAQYLQNTFKVSKRQSFAPFSNPGNAEEECATFSAHSGSLKKQSPKVTFECEQKEENQGKNESNIKPVQTVNITAGFPVVGQKDKPVDNAKCSIKGGSRFCLSSQFRGNETGLITPNKHGLLQNPYRPPLFPIKSFVKTKCKKNLLEENFEEHSMSPEREMGNENIPSTVSTISRNNIRENVFKEASSSNINEVGSSTNEVGSSINEIGSSDENIQAELGRNRGPKLNAMLRLGVLQPEVYKQSLPGSNCKHPEIKKQEYEEVVQTVNTDFSPYLISDNLEQPMGSSHASQVCSETPDDLLDDGEIKEDTSFAENDIKESSAVFSKSVQKGELSRSPSPFTHTHLAQGYRRGAKKLESSEENLSSEDEELPCFQHLLFGKVNNIPSQSTRHSTVATECLSKNTEENLLSLKNSLNDCSNQVILAKASQEHHLSEETKCSASLFSSQCSELEDLTANTNTQDPFLIGSSKQMRHQSESQGVGLSDKELVSDDEERGTGLEENNQEEQSMDSNLGEAASGCESETSVSEDCSGLSSQSDILTTQQRDTMQHNLIKLQQEMAELEAVLEQHGSQPSNSYPSIISDSSALEDLRNPEQSTSEKAVLTSQKSSEYPISQNPEGLSADKFEVSADSSTSKNKEPGVERSSPSKCPSLDDRWYMHSCSGSLQNRNYPSQEELIKVVDVEEQQLEESGPHDLTETSYLPRQDLEGTPYLESGISLFSDDPESDPSEDRAPESARVGNIPSSTSALKVPQLKVAESAQSPAAAHTTDTAGYNAMEESVSREKPELTASTERVNKRMSMVVSGLTPEEFMLVYKFARKHHITLTNLITEETTHVVMKTDAEFVCERTLKYFLGIAGGKWVVSYFWVTQSIKERKMLNEHDFEVRGDVVNGRNHQGPKRARESQDRKIFRGLEICCYGPFTNMPTDQLEWMVQLCGASVVKELSSFTLGTGVHPIVVVQPDAWTEDNGFHAIGQMCEAPVVTREWVLDSVALYQCQELDTCLIPQIPHSHY"
        actual = makeProteinString(variant, "ENST00000357654.9", "chr17", REF, ALT, 5558)
        expected = ('_p.Tyr1852Cys', ('1852', 'Y', 'C'))
        self.assertEqual(actual, expected)
