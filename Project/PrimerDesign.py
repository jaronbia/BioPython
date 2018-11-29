from Bio.Seq import Seq
from NucleicAcid import RNA

class Pcr_Primer:
    # PCR chemical reaction requirements as bellow:
    # length of 20 bases
    # 40-60% G/C content
    # start and end with 1 G/C pairs
    # melting temperature (Tm) of 50-60 degree.
    # primer pairs should have a Tm within 5 degree of each other
    # primer pairs should not have complementary regions
    def __init__(self, seq):
        self.seq_str = ''
        for i in range(len(seq)):
            if seq[i] is not "\n":
                self.seq_str += seq[i]

    def compltemented_seq(self, seq=None):
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

    # The F(Tm) = the total number of C and G * 4 + the total number of A and T * 2
    def Tm_calculator(self, seq=None):
        if len(seq) is not 20:
            print("The primer is not with a right length of 20")
        else:
            Tm = 0
            Tm = seq.count("C") * 4 + seq.count("G") * 4 + seq.count("A") * 2 + seq.count("T") * 2
            return Tm

    # a quick GC ration return to test all possible pcr products
    def get_GC_ratio(self, seq=None):
        count = seq.count("C") + seq.count("G")
        ratio = (count) / len(seq)
        return ratio

    # All if conditions are used to filter and return all possible pcr products by pcr primer rule
    # allow to enter a flexible pcr product length
    def design(self, product_len=200):
        if len(self.seq_str) < 500:
            print("The sequence is not enough long of 500 base pairs")
        else:
            if product_len < 200 or product_len > 500:
                print("The product length range should be in between 100 and 500")
            else:
                select_list = []
                count = 0
                for i in range(0, len(self.seq_str) - product_len, 1):
                    if self.seq_str[i] is "C" or self.seq_str[i] is "G":
                        if self.seq_str[i + product_len] is "C" or self.seq_str[i + product_len] is "G":
                            if self.get_GC_ratio(self.seq_str[i:i + product_len]) > 0.4 and self.get_GC_ratio(
                                    self.seq_str[i:i + product_len]) < 0.6:
                                if self.Tm_calculator(self.seq_str[i: i + 20]) in range(50, 60):
                                    if abs(self.Tm_calculator(self.seq_str[i: i + 20]) - self.Tm_calculator(
                                            self.seq_str[i + product_len - 20: i + product_len])) < 5:
                                        if self.seq_str[i:i + 20] is not self.compltemented_seq(
                                                self.seq_str[i + product_len - 20: i + product_len]):
                                            count += 1
                                            print(count)
                                            print("primer Tm degree:", self.Tm_calculator(self.seq_str[i:i + 20]))
                                            print("anti-primer Tm degree:", self.Tm_calculator(
                                                self.seq_str[i + product_len - 20: i + product_len]))
                                            select = ["position:", i, "end_position:", i + product_len,
                                                      "primer_part:", self.seq_str[i:i + 20], "whole_sequence:",
                                                      self.seq_str[i:i + product_len]]
                                            print(select)
                                            select_list.append(select)
                return select_list
    # allow to enter a specific spot to check the point mutation
    # if there is no filtered result, it returns the list is empty
    # else it will collect all suitable selections
    def filter(select_list = None, position = None):
        filter_list = []
        for i in range(len(select_list)):
            if select_list[i][1] < position-1 and select_list[i][3] > position-1:
                print(select_list[i])
                filter_list.append(select_list[i])
        return filter_list

# to remove the new line from the text file
def str_format(seq=None):
    str = ''
    for i in range(len(seq)):
        if seq[i] is not "\n":
            str += seq[i]
    return str

# to adjust the pcr product output is divisible by 3 and translate rna to amino acid
def adjust_frame(string = None, start = None, end = None):
    new_start = start - start%3
    new_end = end + end%3
    i = 0
    if new_start > 9:
        i = 5
        print(new_start, " " * ((new_end - new_start) - i), new_end)
    elif new_start > 99:
        i = 6
        print(new_start, " " * ((new_end - new_start) - i), new_end)
    else :
        i = 4
        print(new_start, " " * ((new_end - new_start) - i), new_end)
    print("|", " "*((new_end-new_start)-3), "|")
    print(string[new_start:new_end+1])
    return string[new_start:new_end+1]

# receive the point mutation data dictionary to analyze the positions of mutations, the change of PCR original product
# to mutated product, the change of amino acid from original to mutated in protein sequence for disease research
def amino_aicd_mutant(dic = None, new_start = 0):
    analyze_list=[]
    analyze_line=[]
    for x in dic.keys():
        analyze_line = ["point_mutation_position:", x+new_start, "code:", dic[x][0], "->", dic[x][1], "amino_acid:", original_pcr.find_amino_acid(DnaToRna(dic[x][0])), "->", original_pcr.find_amino_acid(DnaToRna(dic[x][1]))]
        print(analyze_line)
        analyze_list.append(analyze_line)
    if len(analyze_list) is 0:
        print("There is no point mutation")
    else:
        return analyze_list

# turn a cDNA (PCR product) back to mRNA coding
def DnaToRna(seq=None):
    return seq.replace("T", "U")


if __name__ == "__main__":
    print("\n")
    print("\tWelcome to Point Mutation Diagnosis and PCR Primer Design tool")
    opened = False
    while opened is False:
        try:
            filename = input("Please enter the original sequence fasta file:")
            handle = open(filename, "r")
            o = Seq(handle.read())
            o = str_format(o)
            opened = True
        except:
            opened = False
    found = False
    while found is False:
        print("enter the expected pcr length:")
        length = int(input())
        pcr_design_result = Pcr_Primer(o).design(length)
        print("enter the target spot position:")
        spot = int(input())
        filtered = Pcr_Primer.filter(pcr_design_result, spot)
        if len(filtered) is 0:
            print("There is no filtered result for the pcr")
            found = False
        else:
            found = True
    print("enter new start point for pcr:")
    new_start = int(input())
    print("enter new end point for pcr:")
    new_end = int(input())
    original_pcr = RNA(adjust_frame(o, new_start, new_end))
    handle.close()
    opened = False
    while opened is False:
        try:
            print("Please enter the mutated sequence fasta file:")
            filename = input()
            handle = open(filename, "r")
            m = Seq(handle.read())
            m = str_format(m)
            opened = True
        except:
            opened = False
    mutation = original_pcr.compare_for_mutations(adjust_frame(m, new_start, new_end))
    amino_aicd_mutant(mutation, new_start)
    handle.close()
    print("Thank you for using the tool")