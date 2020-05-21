
from utility.IUPACData import * # variable import
from Bio.Seq import translate,transcribe




def validateseq(dna_seq):
    """ Program to check if Nucleotide sequence is Valid """
    tempseq = dna_seq.upper()
    tmpfreqdict = { "A":0, "C":0, "G":0, "T":0 }
    for nuc in tempseq:
        if nuc not in nucleotides:
            return "error"
        
    for nuc in tempseq:
        tmpfreqdict[nuc] += 1
        
    a = tmpfreqdict["A"]
    c = tmpfreqdict["C"]
    g = tmpfreqdict["G"]
    t = tmpfreqdict["T"]
        
    return tempseq,a,c,g,t


def validaterna(rna_seq):
    """ Program to check if Nucleotide RNA sequence is Valid """
    tempseq = rna_seq.upper()
    UnambiguousRNA = ["A","U","G","C"]
    for nuc in tempseq:
        if nuc not in UnambiguousRNA:
            return "error"
        
    return tempseq


def validateprotein(protein_seq):
    """ Program to check if Protein sequence is valid """
    tempseq = protein_seq.upper()
    Unambiguousprotein = ["A","C","D","E","F","G","H","I","K","L","M","N","O","P","Q","R","S","T","U","V","W","Y"]
    for pro in tempseq:
        if pro not in Unambiguousprotein:
            return "error"
    
    return tempseq


def countrna(seq):
    """Program to count RNA nucleobases """
    tmpfreqdictrna = { "A":0, "C":0, "G":0, "U":0 }
    for nuc in seq:
        tmpfreqdictrna[nuc] += 1

    a = tmpfreqdictrna["A"]
    c = tmpfreqdictrna["C"]
    g = tmpfreqdictrna["G"]
    u = tmpfreqdictrna["U"]

    return a,c,g,u


def countprotein(seq):
    tmpfreqdictprotein = {
    "A": 0,
    "C": 0,
    "D": 0,
    "E": 0,
    "F": 0,
    "G": 0,
    "H": 0,
    "I": 0,
    "K": 0,
    "L": 0,
    "M": 0,
    "N": 0,
    "O": 0,
    "P": 0,
    "Q": 0,
    "R": 0,
    "S": 0,
    "T": 0,
    "U": 0,
    "V": 0,
    "W": 0,
    "Y": 0
    }
    for pro in seq:
        tmpfreqdictprotein[pro] += 1

    return tmpfreqdictprotein["A"], tmpfreqdictprotein["C"], tmpfreqdictprotein["D"], tmpfreqdictprotein["E"], tmpfreqdictprotein["F"], tmpfreqdictprotein["G"], tmpfreqdictprotein["H"], tmpfreqdictprotein["I"], tmpfreqdictprotein["K"], tmpfreqdictprotein["L"], tmpfreqdictprotein["M"], tmpfreqdictprotein["N"], tmpfreqdictprotein["O"], tmpfreqdictprotein["P"], tmpfreqdictprotein["Q"], tmpfreqdictprotein["R"], tmpfreqdictprotein["S"], tmpfreqdictprotein["T"], tmpfreqdictprotein["T"], tmpfreqdictprotein["U"], tmpfreqdictprotein["V"], tmpfreqdictprotein["W"], tmpfreqdictprotein["Y"]


def gc_content(seq):
    """ Returns percentage of G and C nucleotides in a DNA/RNA sequence. """
    return ("%.2f" % float((seq.count("G") + seq.count("C")) / len(seq) * 100 ))

def GC(seq):
    """Calculate G+C content, return percentage (as float between 0 and 100).

    Copes mixed case sequences, and with the ambiguous nucleotide S (G or C)
    when counting the G and C content.  The percentage is calculated against
    the full length, e.g.:

    >>> from Bio.SeqUtils import GC
    >>> GC("ACTGN")
    40.0

    Note that this will return zero for an empty sequence.
    """
    gc = sum(seq.count(x) for x in ["G", "C", "g", "c", "S", "s"])
    try:
        return gc * 100.0 / len(seq)
    except ZeroDivisionError:
        return 0.0

def gc_content_subsec(seq,k=22):
    """ GC content in a DNA/RNA subsequence having length k , k = 22 by default """
    res = []
    for i in range(0,len(seq) - k + 1, k ):
        subseq = seq[i:i+k]
        res.append(gc_content(subseq))
        return res

def reverse(seq):
    """ Return reverse of input string """
    return seq[::-1]

def complement(seq):
    """Return the complement sequence """
    mapping = str.maketrans('ATGC', 'TACG')
    return seq.translate(mapping)


def reverse_complement(seq):
    """ Computes the reverse complement of the DNA sequence Computes the 
    reverse complement of the DNA sequence"""

    return complement(seq)[::-1]

 
def molecular_weight(seq, seq_type=None, double_stranded=False, circular=False, monoisotopic=False):
    """Calculate the molecular mass of DNA, RNA or protein sequences as float.
    Only unambiguous letters are allowed. Nucleotide sequences are assumed to
    have a 5' phosphate. """

    if seq_type == "DNA":
        if monoisotopic:
            weight_table = monoisotopic_unambiguous_dna_weights
        else:
            weight_table = unambiguous_dna_weights
    elif seq_type == "RNA":
        if monoisotopic:
            weight_table = monoisotopic_unambiguous_rna_weights
        else:
            weight_table = unambiguous_rna_weights
    elif seq_type == "protein":
        if monoisotopic:
            weight_table = monoisotopic_protein_weights
        else:
            weight_table = protein_weights
    else:
        return("Allowed seq_types are DNA, RNA or protein, not %r" % seq_type)

    if monoisotopic:
        water = 18.010565
    else:
        water = 18.0153

    try:
        weight = sum(weight_table[x] for x in seq) - (len(seq) - 1) * water
        if circular:
            weight -= water
    except KeyError as e:
        raise ValueError(
            "%s is not a valid unambiguous letter for %s" % (e, seq_type)
        ) from None

    if seq_type in ("DNA", "RNA") and double_stranded:
        seq = str(seq.complement())
        weight += sum(weight_table[x] for x in seq) - (len(seq) - 1) * water
        if circular:
            weight -= water
    elif seq_type == "protein" and double_stranded:
        raise ValueError("double-stranded proteins await their discovery")

    return weight

def translate_seq(self, init_pos=0):
        """Translates a DNA sequence into an aminoacid sequence"""
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]


def six_frame_translations(seq, genetic_code=1):
    """Return pretty string showing the 6 frame translations and GC content.
    coded and written by casesagar aka sagar saini

    >>> from Bio.SeqUtils import six_frame_translations
    >>> print(six_frame_translations("AUGGCCAUUGUAAUGGGCCGCUGA"))
    GC_Frame: a:5 t:0 g:8 c:5 
    Sequence: auggccauug ... gggccgcuga, 24 nt, 54.17 %GC
    <BLANKLINE>
    <BLANKLINE>
    1/1
      G  H  C  N  G  P  L
     W  P  L  *  W  A  A
    M  A  I  V  M  G  R  *
    auggccauuguaaugggccgcuga 
    uaccgguaacauuacccggcgacu
    A  M  T  I  P  R  Q 
     H  G  N  Y  H  A  A  S
      P  W  Q  L  P  G  S
    <BLANKLINE>
    <BLANKLINE>

    """  # noqa for pep8 W291 trailing whitespace
    from Bio.Seq import reverse_complement, translate

    anti = reverse_complement(seq)
    comp = anti[::-1]
    length = len(seq)
    frames = {}
    for i in range(0, 3):
        fragment_length = 3 * ((length - i) // 3)
        frames[i + 1] = translate(seq[i : i + fragment_length], genetic_code,)
        frames[-(i + 1)] = translate(anti[i : i + fragment_length], genetic_code,)[::-1]

    # create header
    if length > 20:
        short = "%s ... %s" % (seq[:10], seq[-10:])
    else:
        short = seq
    header = "GC_Frame: "
    for nt in ["a", "t", "g", "c"]:
        header += "%s:%d " % (nt, seq.count(nt.upper()))

    header += "\nSequence: %s, %d nt, %0.2f %%GC\n\n\n" % (
        short.lower(),
        length,
        GC(seq),
    )
    res = header

    

        
    
    frame_3 = frames[3]
    frame_2 = frames[2]
    frame_1 = frames[1]
    # seq

    # - frames
    framecomp_2 = frames[-2]
    framecomp_1 = frames[-1]
    framecomp_3 = frames[-3]
    return res, frame_1, frame_2, frame_3, seq, comp, framecomp_2, framecomp_1, framecomp_3

