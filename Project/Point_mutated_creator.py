'''
Name: George Tsai
Project: NucleicAcid Class
Date: 11/23/2018
Copyrighted Year: 2018
'''


'''
for i in range(5):
    for j in range(5):
        print(i, j)

l=[]
for i in range(5):
    for j in range(5):
        subl=[i, j]
        l.append(subl)
print(l)

for i in range(len(l)):
    l[i].append(i)

print(l)
print(len(l))
print(min(l))

l=[["a", "b", 20], ["c", "d", 10], ["e", "f", 5]]
for i in range(len(l)):
    l[i].insert(0, i)
print(l)
print(max(l))
print(min(l))
'''

def str_format(seq=None):
    str = ""
    for i in range(len(seq)):
        if seq[i] is not "\n":
            str += seq[i]
    return str

print("----------"*5)

print("Hemophilia - F9")
handle = open("Homo sapiens coagulation factor IX (F9) mRNA.txt", "r")
mRNA_string = handle.read()
mRNA_string = str_format(mRNA_string)
print(len(mRNA_string))
print(mRNA_string[105:108])

patient_mRNA = []
for i in range(len(mRNA_string)):
    patient_mRNA.append(mRNA_string[i])
patient_mRNA[106] = "A"
print(patient_mRNA[105:108])

patient_fasta = "".join(patient_mRNA)
patient_file = open("Homo sapiens coagulation factor IX (F9) mutated mRNA.txt", "w")
patient_file.write(patient_fasta)

handle.close()
patient_file.close()

print("----------"*5)
print("Cancer - p53")

handle = open("Homo sapiens tumor protein p53 (TP53) mRNA.txt", "r")
mRNA_string = handle.read()
mRNA_string = str_format(mRNA_string)
print(len(mRNA_string))
print(mRNA_string[208:246])

patient_mRNA = []
for i in range(len(mRNA_string)):
    patient_mRNA.append(mRNA_string[i])

patient_mRNA[200] = "A"
patient_mRNA[208] = "C"
patient_mRNA[218] = "G"
patient_mRNA[245] = "T"
print(patient_mRNA[243:246])

patient_fasta = "".join(patient_mRNA)
patient_file = open("Homo sapiens tumor protein p53 (TP53) mutated mRNA.txt", "w")
patient_file.write(patient_fasta)

handle.close()
patient_file.close()

print("----------"*5)
print("Sickle cell anemia - HBB")

handle = open("Homo sapiens hemoglobin subunit beta (HBB) mRNA.txt", "r")
mRNA_string = handle.read()
mRNA_string = str_format(mRNA_string)
print(len(mRNA_string))
print(mRNA_string[159:162])

for i in range(0, len(mRNA_string)-3, 3):
    if mRNA_string[i]=="C" and mRNA_string[i+1]=="T" and mRNA_string[i+2]=="T":
        print(i, mRNA_string[i:i+3])

patient_mRNA = []
for i in range(len(mRNA_string)):
    patient_mRNA.append(mRNA_string[i])
patient_mRNA[160] = "A"
print(patient_mRNA[159:162])

patient_fasta = "".join(patient_mRNA)
print(patient_fasta[159:162])
patient_file = open("Homo sapiens hemoglobin subunit beta (HBB) mutated mRNA.txt", "w")
patient_file.write(patient_fasta)

handle.close()
patient_file.close()

handle = open("Homo sapiens hemoglobin subunit beta (HBB) mutated mRNA.txt", "r")
mRNA_string = handle.read()
mRNA_string = str_format(mRNA_string)
print(len(mRNA_string))
print(mRNA_string[159:162])