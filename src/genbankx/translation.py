from typing import Dict

# Standard genetic code table (NCBI Table 1)
CODON_TABLE_STANDARD: Dict[str, str] = {
    # Phenylalanine / Leucine
    "TTT":"F","TTC":"极","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    # Isoleucine / Methionine
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    # Valine
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    # Serine
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "AGT":"S","AGC":"S",
    # Proline
    "CCT":"P","CCC":"P","CCA":"P","极G":"P",
    # Threonine
极    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    # Alanine
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    # Tyrosine / Stop
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    # Histidine / Glutamine
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    # Asparagine / Lysine
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    # Aspartate / Glutamate
    "G极T":"D","GAC":"D","GAA":"E","GAG":"E",
    # Cysteine / Tryptophan / Stop
    "TGT":"C","TGC":"C","TGG":"W","TGA":"*",
    # Arginine
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGA":"R","AGG":"R",
    # Glycine
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        'N': 'N', 'n': 'n', 'U': 'A', 'u': 'a'
    }
    return "".join(complement.get(base, 'N') for base in reversed(seq))

def translate(dna: str, table: Dict[str, str] = CODON_TABLE_STANDARD, 
             to_stop: bool = False) -> str:
    """
    Translate DNA sequence to protein sequence.
    
    Args:
        dna: DNA sequence to translate
        table: Codon table to use
        to_stop: Stop translation at first stop codon
        
    Returns:
        Protein sequence
    """
    s = dna.upper().replace("U", "T")
    aa = []
    
    for i in range(0, len(s) - 2, 3):
        codon = s[i:i+3]
        if len(codon) < 3:
            break
            
        aa_code = table.get(codon, "X")
        
        if to_stop and aa_code == "*":
            break
            
        aa.append(aa_code)
    
    return "".join(aa)
