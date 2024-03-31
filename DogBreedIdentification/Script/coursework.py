from Bio import SeqIO, Phylo
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def load_sequences_from_fasta(file_path):
#Reads and returns sequences from a specified FASTA file for analysis.
    return list(SeqIO.parse(file_path, "fasta"))

def identify_similar_sequence(known_sequences, query_sequence):
    #Finds the sequence that most closely matches the query sequence among the known sequences using PairwiseAligner.
    sequence_aligner = PairwiseAligner()
    highest_score = float('-inf')
    most_similar_seq = None
    for sequence in known_sequences:
        #Alignment score is calculated for each sequence and the highest score is stored.
        alignment_score = sequence_aligner.score(query_sequence.seq, sequence.seq)
        if alignment_score > highest_score:
            highest_score = alignment_score
            most_similar_seq = sequence
    return most_similar_seq, highest_score

def build_phylo_tree(input_sequences):
    # Creates a phylogenetic tree from the input sequences.
    seq_alignment = MultipleSeqAlignment(input_sequences) #Sequence alignment is performed.
    dist_calc = DistanceCalculator('identity')
    dist_matrix = dist_calc.get_distance(seq_alignment)
    tree_builder = DistanceTreeConstructor()
    phylo_tree = tree_builder.nj(dist_matrix) #Construction of the phylogenetic tree.
    return phylo_tree

def execute():
    path_to_database_sequences = '/Users/alexandroskaikkis/Desktop/project_dog_dna 2/dog_breeds.fa'
    path_to_query_sequence = '/Users/alexandroskaikkis/Desktop/project_dog_dna 2/mystery.fa'
    # Sequences are loaded from the specified FASTA files.
    known_sequences = load_sequences_from_fasta(path_to_database_sequences)
    query_sequence = load_sequences_from_fasta(path_to_query_sequence)[0]
    #Identtify the closest match for the query sequence.
    similar_sequence, alignment_score = identify_similar_sequence(known_sequences, query_sequence)
    print(f"Closest match: {similar_sequence.id} with alignment score: {alignment_score}")
    # Includes all sequences for the phylogenetic analysis.
    all_sequences_for_analysis = known_sequences + [query_sequence]
    phylogenetic_tree = build_phylo_tree(all_sequences_for_analysis)
    Phylo.write(phylogenetic_tree, "phylogenetic_tree.xml", "phyloxml") #The phylogenetic tree is saved to a file.
    print("Phylogenetic tree has been successfully created and stored as 'phylogenetic_tree.xml'.")

if __name__ == "__main__":
    execute()
