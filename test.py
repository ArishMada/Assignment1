# Provided codon table
table = {
    'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
    'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
    'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
    'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_',
    'UGC': 'C', 'UGU': 'C', 'UGA': '_', 'UGG': 'W',
}


# Function to translate an amino acid sequence to RNA using the provided codon table
def translate_amino_acid_sequence_to_rna(table, amino_acid_sequence):
    rna_seq= ""
    for amino_acid in amino_acid_sequence:
        found = False
        for codon, codon_amino_acid in table.items():
            if codon_amino_acid == amino_acid:
                rna_seq += codon.replace("T", "U")
                found = True
                break
        if not found:
            return None  # Amino acid not found in the codon table
    return rna_seq


def codon_frequency(rna_seq):
    frequency = {}

    for i in range(0, len(rna_seq), 3):
        cod = rna_seq[i:i + 3]
        frequency[cod] = frequency.get(cod, 0) + 1

    return frequency


# Given amino acid sequence
input_aminoacid = input("Input aminoacid = ").upper()

# Translate the given amino acid sequence to RNA
rna_sequence = translate_amino_acid_sequence_to_rna(table, input_aminoacid)

if rna_sequence:
    print(f"RNA Sequence for {input_aminoacid}: {rna_sequence}")

    for codon, count in codon_frequency(rna_sequence).items():
        print(f"{codon}: {count}")
else:
    print(f"No RNA sequence found for the amino acid sequence {input_aminoacid}.")
