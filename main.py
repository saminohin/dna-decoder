#from dna_decoder_toolkit import * # import everything for the dna_toolkit module 
# imports everything from the dna analysis tool class
from dna_decoder_toolkit import analyze_tools
# import data for the grab data class 
from dna_decoder_toolkit import dna 
#imports everything for the save_data_as_file class 
from dna_decoder_toolkit import save_data_function
from dna_decoder_toolkit import save_file_data 


random_dna = "ATGCGTA"
rna = analyze_tools.transcription(dna)
translate_rna = analyze_tools.translate_rna(rna)

#This function outputs results on the terminal: 
def output_in_the_terminal():
    print(f'Dna toolkit module {analyze_tools.get_version()}') # gets the code version number 

    print(f'The length of the sequence is {len(dna)}') # prints out the length of the sequence

    print(analyze_tools.gc_content(dna)) # prints out gc content 
    print(analyze_tools.at_content(dna)) # prints out at content 

    print(analyze_tools.count_nucleotides(dna)) # counts A, T, G, C 

    print(f'Dna sequence:{dna}') #print out dna sequence

    print(f'complement_strand: {analyze_tools.complement_strand(dna)}') # prints out the dna's complementary strand

    print(f'Rna sequence:   {analyze_tools.transcription(dna)}') #  transcription --> thymine to uracil

    print(f'Rna compleated: {analyze_tools.complement_rna_strand(rna)}') # print out rna and its complementary sequence  

    print(f'translate: {analyze_tools.translate(dna)}')  # translation --->  dna Sequence ---> to Protein 
    print(f'translate_rna: {analyze_tools.translate_rna(dna)}\n') # translation --->  rna Sequence ---> to Protein 

    print("alignment_dna_sequence:\n")
    print(analyze_tools.alignment_dna_sequence(dna),"\n") # print out dna and its complementary sequence
    
    print("rna_alignment:\n")
    '''This two line does the same thing: '''
     
    #print(rna + "\n" + f"{''.join(['|' for c in range(len(dna))])} ""\n" + analyze.complement_rna_strand(rna))
    print(analyze_tools.rna_alignment(rna)) # alignment for rna and complementary strand
    
    print(f'original: {random_dna} , random mutation: {analyze_tools.random_mutation(random_dna)}') 
    print(f'Mutation count:{analyze_tools.count_mutations("TGC", "TCC")}') # counts all the mutation in the sequence 
    
    save_file_data.save_translate_data('output/translate_data.txt') # you can also put this string of text as a variable
output_in_the_terminal() 

