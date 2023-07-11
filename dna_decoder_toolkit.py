# importing the random module 
import random

#load sequence data files
filename = 'data/dna_data.txt' # your file name goes here
#filename = ' ' # your file name goes here

# another way you can do this is: 
#filename = input("Please enter your dna file here: ") # This way you can just enter your file name here as input  

#read sequence data files
class grab_data:
    def read_sequence_dna(self, filename):
        with open(filename, "r") as file:
            dna_seq_data = file.read()
        dna_seq_data = dna_seq_data.replace("\n", "")
        return dna_seq_data

grab = grab_data()
dna = grab.read_sequence_dna(filename)

# Analysis class tool 
class dna_analysis_tool:
    def __init__(self): 
        self.count_A = dna.count("A")
        self.count_T = dna.count("T")
        self.count_G = dna.count("G")
        self.count_C = dna.count("C")
        print("Dna Decoder Toolkit Initiated\n") # This is the first line that should print in the terminal 

    # gives the version  number 
    def get_version(self):
        return "version 2.1.1"
    
    # return the seq as a string 
    def Seq(self, dna):
        return str(dna)
    
    # returns the reverse of the original seq 
    def reverse_seq(self, dna):
        return dna[::-1]
    
    # converting thymine to uracil 
    def transcription(self, dna):
        output = dna.replace("T", "U")
        return output
    
    # count the four base pair necleotides in the seq 
    def count_nucleotides(self, dna):
        output = f"A: {self.count_A}, T: {self.count_T}, G: {self.count_G}, C: {self.count_C}"
        return output 
    
    # index the sequence
    def indexing_seq(self):
        for index, letter in enumerate(dna):
            print("%i %s" % (index, letter))

    # generate a random neclueotides of sequence 
    def generate_random_neclueotides(self, strand_size):
        neclueotides = ["A", "T", "G", "C"]
        for x in range(strand_size): # how long you want your strand to be 
            random_neclueotides = random.choice(neclueotides)
            list_neclueotides = random_neclueotides
            #print(list_neclueotides + "\n", end='')
            print(list_neclueotides, end='')

    # returns a gc content 
    def gc_content(self, dna):
        output = f"GC Content is {round(float(self.count_C + self.count_G) / len(dna) * 100)}%"
        return output
     
    # returns a at content   
    def at_content(self, dna):
        output = f"AT Content is {round(float(self.count_A + self.count_T) / len(dna) * 100)}%"
        return output
    
    # complete the sequence
    def complement_strand(self, dna):
        complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
        bases = list(dna)
        bases = [complement[base] for base in bases]
        return ''.join(bases)
    

    # complete the  rna sequence
    def complement_rna_strand(self, rna):
        complement_rna = {"A": "U", "C": "G", "G": "C", "U": "A"}
        bases = list(rna)
        bases = [complement_rna[base] for base in bases]
        return ''.join(bases)
    

    # translation --->  Sequence ---> to Protein 
    def translate(self, dna):
        table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
            'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
        }

        # logic portion of the code
        protein = ""
        dna = dna.upper()

        for i in range(0, len(dna), 3):
            codon = dna[i:i+3]
            if codon in table:
                protein += table[codon]
            else:
                protein += 'X' # Represents unknown amino acid
        return protein


    # Fasta format file 
    def fasta_format(self, name, dna):
        return f"{name}\n%s\n" % dna

    ''' dna alignment sequence example: ATGC 
                                        ||||
                                        TACG
    '''
    

    def alignment_dna_sequence(self, dna):
        result = dna + "\n" + f"{''.join(['|' for c in range(len(dna))])} ""\n" + analyze_tools.complement_strand(dna)
        return result


    # translation ---> RNA Sequence ---> to Protein 
    def translate_rna(self, rna):
        rna_table = {
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
            'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*',
            'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W',
        }

        protein = ""
        rna.upper()

        for i in range(0, len(rna), 3):
            codon = rna[i:i+3]
            if codon in rna_table:
                protein += rna_table[codon]
            else:
                protein += 'X'
        return protein


    ''' alignment the sequence enample: AUGC
                                        ||||
                                        UACG 
    '''

    def rna_alignment(self, rna):
        return rna + "\n" + f"{''.join(['|' for c in range(len(dna))])} ""\n" + analyze_tools.complement_rna_strand(rna)
    

    # generates a random mutation sequence
    def random_mutation(self, dna):
        dna_list = list(dna)
        mutation_site = random.randint(0, len(dna_list) - 1)
        dna_list[mutation_site] = random.choice(list('ATCG'))
        return ''.join(dna_list)
    
    ''''
    This line of code looks for any mutation in the sequence any sequence that is identical example A == A or G == G
    
    '''
    # count all the mutation seq 
    def count_mutations(self, dna, complementary):
        mutation_count = 0
        if len(dna) != len(complementary):
            return "Error! DNA Sequences must have the same length."
        for dna_sequence1, dna_sequence2 in zip(dna, complementary):
            if dna_sequence1 == dna_sequence2: # if dna_sequrnce1 != dna_sequrnce2:
                mutation_count += 1
        return mutation_count


analyze_tools = dna_analysis_tool()

# This line need to be in the main.py file 
# This class save all of the data and output and save it into a file: 
class save_data_function:
    def save_translate_data(self, save_filename):
        with open(save_filename, "w") as save_file:
            save_file.write(analyze_tools.translate(dna))

save_file_data = save_data_function()


