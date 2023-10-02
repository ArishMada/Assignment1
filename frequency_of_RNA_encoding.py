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


def translate_amino_acid_sequence_to_rna(amino_acid_sequence):
    possible_rna_sequences = [""]

    for amino_acid in amino_acid_sequence:
        possible_codons = []
        for cod, codon_amino_acid in table.items():
            if codon_amino_acid == amino_acid:
                possible_codons.append(cod)

        if not possible_codons:
            return None  # Amino acid not found in the codon table

        possible_sequences = []
        for rna_seq in possible_rna_sequences:
            for cod in possible_codons:
                possible_sequences.append(rna_seq + cod)

        possible_rna_sequences = possible_sequences

    return possible_rna_sequences


def codon_frequency(rna_seq):
    frequency = {}

    for i in range(0, len(rna_seq), 3):
        cod = rna_seq[i:i + 3]
        frequency[cod] = frequency.get(cod, 0) + 1

    return frequency


# main code
input_aminoacid = input("Amino Acid = ").upper()  # Ask the user for the aminoacid

rna_sequences = translate_amino_acid_sequence_to_rna(input_aminoacid)

if rna_sequences:
    print(f"All Possible RNA Sequences for {input_aminoacid}:")
    for rna_sequence in rna_sequences:
        print(f"RNA Sequence: {rna_sequence}")

        for codon, count in codon_frequency(rna_sequence).items():
            print(f"{codon}: {count}")
        print()
else:
    print(f"No RNA sequences found {input_aminoacid}.")
