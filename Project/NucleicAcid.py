'''
Name: Jaron Bialecki
Project: NucleicAcid Class
Date: 11/23/2018
Copyrighted Year: 2018
'''


from Bio.Seq import Seq
import random


class NucleicAcid:
    def __init__(self, seq):
        self.__sequence = seq
        self.__type_of_Nacid = self.find_base_type()
        self.__base_pool = "ATCCGG"

    def find_base_type(self):  # Checks to see if it is DNA or RNA
        if "T" in self.__sequence:
            return "DNA"
        return "RNA"

    def generate_seq(self, length=None): # Function for testing purposes
        new_seq = ''
        for i in range(length):
            new_seq += random.choice(self.__base_pool)
        return new_seq

    def ratio_of_base(self, chemical_base):  # Finds ratio of A T U C or G ratio w/in Sequence
        return self.__sequence.count(chemical_base) / len(self.__sequence)  # Occurence / Size


class DNA(NucleicAcid):
    def __init__(self, seq, coding_seq=None):
        NucleicAcid.__init__(self, seq)
        self.__DNA_sequence = seq
        self.__DNA_coding = coding_seq
        self.__DNA_template = self.__generate_dna_template()

    def __generate_dna_template(self):
        if self.__DNA_coding is None:
            return None
        return self.__DNA_coding.complement()

    def generate_mRNA(self):
        if self.__DNA_coding is None:
            return None
        return self.__DNA_sequence.transcribe()

    def generate_protein_dna(self):
        if self.__DNA_coding is None:
            return None
        return self.__DNA_coding.translate()


class RNA(NucleicAcid): # Current version of RNA class represents mRNA
    def __init__(self, seq):
        NucleicAcid.__init__(self, seq)
        self.__RNA_sequence = seq
        self.__typeOfRNA = ''  # t, m and r types of RNA
        self.amino_acids = ["Phe", "Leu", "Ile", "Met", "Val", "Ser", "Pro", "Thr", "Ala", "Tyr"
                            , "stop", "His", "Gln", "Asn", "Lys", "Asp", "Glu", "Cys", "Trp", "Arg"
                            , "Gly"]
        self.codon_encoding = {"UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
                               "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
                               "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
                               "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
                               "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
                               "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
                               "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
                               "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
                               "UAU": "Tyr", "UAC": "Tyr", "UAA": "stop", "UAG": "stop",
                               "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
                               "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
                               "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
                               "UGU": "Cys", "UGC": "Cys", "UGA": "stop", "UGG": "Trp",
                               "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
                               "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
                               "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"}

    def generate_protein_mRNA(self):
        return self.__RNA_sequence.translate()

    def find_amino_acid(self, codon_base):  # Returns amino acid associated with codon_base
        for key in self.codon_encoding.keys():
            if codon_base == key:
                return self.codon_encoding[key]
        return None  # If codon_base was non-existent it returns None

    '''
    Returns positions where base codon has been mutated. If position is 0 no mutation 
    had occured at that point. If positive that means the letter value shifted downward
    and if negative the letter shifted upwards.
    '''

    def compare_for_mutations(self, mutation):
        og_seq = str(self.__RNA_sequence)
        mu_seq = str(mutation)
        mutation_list = {}
        start = 0
        stop = 2

        def find_mutations(original_base, mutation_base):
            mu_positions = []  # Positions at which mutations occured
            original_base = str(original_base)
            mutation_base = str(mutation_base)
            for index in range(0, len(original_base)):  # Negative numbers indicate shift up & positive means downward
                mutation_shift = ord(original_base[index]) - ord(mutation_base[index])  # Calculates mutation
                mu_positions.append(mutation_shift)
            return mu_positions

        while stop < len(self.__RNA_sequence):  # Calculates at which position mutations occured
            og_base = og_seq[start:stop+1]
            mu_base = mu_seq[start:stop+1]
            if mu_base != og_base:
                codon_base_position = start
                mutation_positions = find_mutations(og_base, mu_base)
                mutation_list.update({codon_base_position: mutation_positions})
            start += 3
            stop += 3

        def mutations(mutations_dict):  # Tell you the original base and mutated base at a position
            original_to_base = {}
            for key in mutations_dict:
                calculate_mutant = ''
                for m_index, elem in enumerate(mutations_dict[key]):
                    calculate_mutant = calculate_mutant + chr(ord(self.__RNA_sequence[key + m_index]) - elem)
                original_to_base.update({key: [str(self.__RNA_sequence[key:key+3]), calculate_mutant]})
            return original_to_base

        return mutations(mutation_list)


'''
x = Seq("ATCGATCG")
#dna_unknown_seq = NucleicAcid(x)
#print(dna_unknown_seq.ratio_of_base("T"))
dna_unknown_seq = DNA(x)
print(dna_unknown_seq.generate_mRNA())
print(dna_unknown_seq.generate_protein_dna())
'''
'''
x = Seq("AUCGAUCG")
rna_unknown_seq = RNA(x)
print(rna_unknown_seq.find_amino_acid("UUU"))
print(rna_unknown_seq.find_mutations("ATCGATCG")) # Not sure if it works in this scenario
'''
'''
y = Seq("AUCGAUCGUAACGUUGGG")
rna_unknown_seq = RNA(y)
mutations = rna_unknown_seq.compare_for_mutations("GUCGAUCGUAACGUUGGA")
print(rna_unknown_seq.compare_for_mutations("GUCGAUCGUAACGUUGGA"))

y = Seq("CCACTTTCA")
rna_unknown_seq = RNA(y)
mutations = rna_unknown_seq.compare_for_mutations("CGGCATACA")
print(mutations)
'''