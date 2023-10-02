codon_table = {
    "UUU": "Phe - (F)", "UUC": "Phe - (F)", "UUA": "Leu - (L)", "UUG": "Leu - (L)",
    "UCU": "Ser - (S)", "UCC": "Ser - (S)", "UCA": "Ser - (S)", "UCG": "Ser - (S)",
    "UAU": "Tyr - (Y)", "UAC": "Tyr - (Y)", "UAA": "_", "UAG": "_",
    "UGU": "Cys - (C)", "UGC": "Cys - (C)", "UGA": "_", "UGG": "Trp - (W)",

    "CUU": "Leu - (L)", "CUC": "Leu - (L)", "CUA": "Leu - (L)", "CUG": "Leu - (L)",
    "CCU": "Pro - (P)", "CCC": "Pro - (P)", "CCA": "Pro - (P)", "CCG": "Pro - (P)",
    "CAU": "His - (H)", "CAC": "His - (H)", "CAA": "Gln - (Q)", "CAG": "Gln - (Q)",
    "CGU": "Arg - (R)", "CGC": "Arg - (R)", "CGA": "Arg - (R)", "CGG": "Arg - (R)",

    "AUU": "Ile - (I)", "AUC": "Ile - (I)", "AUA": "Ile - (I)", "AUG": "Met - (M)",
    "ACU": "Thr - (T)", "ACC": "Thr - (T)", "ACA": "Thr - (T)", "ACG": "Thr - (T)",
    "AAU": "Asn - (N)", "AAC": "Asn - (N)", "AAA": "Lys - (K)", "AAG": "Lys - (K)",
    "AGU": "Ser - (S)", "AGC": "Ser - (S)", "AGA": "Arg - (R)", "AGG": "Arg - (R)",

    "GUU": "Val - (V)", "GUC": "Val - (V)", "GUA": "Val - (V)", "GUG": "Val - (V)",
    "GCU": "Ala - (A)", "GCC": "Ala - (A)", "GCA": "Ala - (A)", "GCG": "Ala - (A)",
    "GAU": "Asp - (D)", "GAC": "Asp - (D)", "GAA": "Glu - (E)", "GAG": "Glu - (E)",
    "GGU": "Gly - (G)", "GGC": "Gly - (G)", "GGA": "Gly - (G)", "GGG": "Gly - (G)",
}

nucleotides_dic = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def check(dna):
    valid = dna.count("A") + dna.count("C") + dna.count("G") + dna.count("T")
    if valid == len(dna):
        return True
    else:
        return False


def complement(dna_sequence):
    dna_sequence = dna_sequence.upper()
    if check(dna_sequence):
        complement_sequence = "".join(nucleotides_dic[nucleotide] for nucleotide in dna_sequence)
        print("Complement = ", complement_sequence)
        return complement_sequence
    else:
        print("Invalid DNA sequence")
        return False


def complement_to_mrna(complement_sequence):
    mrna = complement_sequence.replace("T", "U")
    print("mRNA = ", mrna)
    return mrna


def mrna_to_protein(mrna_sequence):
    protein_sequence = ""
    for i in range(0, len(mrna_sequence), 3):  # divide the mRNA sequence into groups of 3
        codon = mrna_sequence[i:i + 3]
        if codon in codon_table:
            amino_acid = codon_table[codon]
            if amino_acid == "Stop":
                break  # Stop
            protein_sequence += amino_acid + " "
        else:
            protein_sequence += "U "  # Unknown codon = "U"

    print("Aminoacid = ", protein_sequence.strip())


# main code
DNA = input("Input DNA = ")  # Ask the user for the DNA sequence
complement = complement(DNA)
if complement:
    mRNA = complement_to_mrna(complement)
    mrna_to_protein(mRNA)
