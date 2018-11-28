'''
Name: George Tsai
Project: PrimerDesign
Date: 11/23/2018
Copyrighted Year: 2018
'''

from Bio.Seq import Seq
import Project.NucleicAcid as NC


class Pcr_Primer:
    # length of 18-24 bases
    # 40-60% G/C content
    # start and end with 1-2 G/C pairs
    # melting temperature (Tm) of 50-60 degree
    # primer pairs should have a Tm within 5 degree of each other
    # primer pairs should not have complementary regions
    def __init__(self, seq_file):
        self.seq_str = ''
        for i in range(len(seq_file)):
            if seq_file[i] is not "\n":
                self.seq_str+=seq_file[i]

    def compltemented_seq(self, seq = None):
        compltemented_seq = ""
        for i in range(len(seq)):
            if seq[i] is "A":
                compltemented_seq += "T"
            elif seq[i] is "T":
                compltemented_seq += "A"
            elif seq[i] is "C":
                compltemented_seq += "G"
            elif seq[i] is "G":
                compltemented_seq += "C"
        return compltemented_seq

    def Tm_calculator(self, seq=None):
        if len(seq) is not 20:
            print("The primer is not with a right length of 20")
        else:
            Tm = 0
            Tm = seq.count("C") * 4 + seq.count("G") * 4 + seq.count("A") * 2 + seq.count("T") * 2
            return Tm

    def get_GC_ratio(self, seq):
        count = seq.count("C") + seq.count("G")
        ratio = count / len(seq)
        return ratio

    def design(self, product_len = 200):
        if len(self.seq_str) < 1000:    # Make it divisible by 3 #####################################################
            print("The sequence is not enough long of 1000 base pairs")
        else:
            if product_len < 200 and product_len > 400:
                print("The product length range should be in between 200 and 400")
            else:
                select_list = []
                count = 0
                for i in range(0, len(self.seq_str) - product_len, 1):
                    if self.seq_str[i] is "C" or self.seq_str[i] is "G":
                        if self.seq_str[i + 1] is "C" or self.seq_str[i + 1] is "G":
                            if self.seq_str[i + product_len] is "C" or self.seq_str[i + product_len] is "G":
                                if self.seq_str[i + product_len - 1] is "C" or self.seq_str[i + product_len -1] is "G":
                                    if self.get_GC_ratio(self.seq_str[i:i + product_len]) > 0.4 and self.get_GC_ratio(self.seq_str[i:i + product_len]) < 0.6:
                                        if self.Tm_calculator(self.seq_str[i : i+20]) in range(50, 60):
                                            if abs(self.Tm_calculator(self.seq_str[i : i+20]) - self.Tm_calculator(self.seq_str[i+ product_len-20 : i+product_len])) < 5:
                                                if self.seq_str[i:i + 20] is not self.compltemented_seq(self.seq_str[i+product_len-20 : i + product_len]):
                                                    count+=1
                                                    print(count)
                                                    print("primer:", self.Tm_calculator(self.seq_str[i:i + 20]))
                                                    print("anti-primer:", self.Tm_calculator(self.seq_str[i+product_len-20 : i+product_len]))
                                                    select = ["position:", i, "end_position:", i+product_len, "primer_part:", self.seq_str[i:i + 20], "whole_sequence:",
                                                              self.seq_str[i:i + product_len]]
                                                    print(select)
                                                    select_list.append(select)
            return select_list

def str_format(seq=None):
    str = ''
    for i in range(len(seq)):
        if seq[i] is not "\n":
            str += seq[i]
    return str

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
'''

handle = open("Homo sapiens coagulation factor IX (F9) mRNA.txt", "r")
o = Seq(handle.read())
o = str_format(o)
pcr_design_result = Pcr_Primer(o).design(200)
print(pcr_design_result[0][7])

handle = open("Homo sapiens coagulation factor IX (F9) mutated mRNA.txt", "r")
m = Seq(handle.read())
m = str_format(m)
rna_unknown_seq = NC.RNA(m[0:217])
mutation = rna_unknown_seq.compare_for_mutations(o[0:217])
print(mutation)
