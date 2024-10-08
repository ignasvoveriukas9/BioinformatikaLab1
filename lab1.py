from Bio import SeqIO
import numpy as np

def parse_fasta(file):
    for record in SeqIO.parse(file, "fasta"):
        return ( record.id, str(record.seq).upper() )

def split_sequence(sequence):
    sequences = []

    for i in range ( 0, 3 ):
        sequences.append ( split_into_triplets ( sequence [ i: ] ) )

    for i in range ( 0, 3 ):
        sequences.append ( split_into_triplets ( reverse_compliment (sequence) [ i: ] ) )

    return sequences

def split_into_triplets(sequence):
    triplets = []
    for i in range ( 0 , len(sequence), 3 ):
        triplets.append ( sequence [ i : i + 3 ] )

    if len ( triplets [ -1 ] ) < 3:
        del triplets [ -1 ]

    return triplets

def reverse_compliment ( sequence ):
    if set(sequence) == {"A", "T", "C", "G"}.union(set(sequence)):
        sequence = sequence.replace( "A", "t" ).replace( "C", "g" ).replace("T", "a").replace("G", "c")
        sequence = sequence.upper()

        sequence = sequence[::-1]
        return sequence
    elif set(sequence) == {"A", "U", "C", "G"}.union(set(sequence)):
        sequence = sequence.replace( "A", "u" ).replace( "C", "g" ).replace("U", "a").replace("G", "c")
        sequence = sequence.upper()

        sequence = sequence[::-1]
        return sequence
    else:
        return "Invalid sequence"

codontab = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}

def turn_into_codon ( sequences ):
    codons = [[],[],[],[],[],[]]
    for idx, seq in enumerate ( sequences ):
        for triplet in seq:
            codons [idx].append ( codontab [ triplet ] )

    return codons

def find_proteins ( sequence ):
    proteins = []

    proteinStart = -1
    for i in range ( 0, len ( sequence ) ):
        if ( sequence [ i ] == 'M' and proteinStart == -1 ):
            proteinStart = i
        elif ( sequence [ i ] == '*' and proteinStart != -1 ):
            if len ( sequence [ proteinStart  : i + 1 ] ) >= 33:
                proteins.append ( sequence [ proteinStart + 1  : i ] )
            proteinStart = -1
    return proteins

codonList = [ 'S', 'F', 'L', 'Y', 'C', 'W', 'P', 'H', 'Q', 'R', 'I', 'M', 'T', 'N', 'K', 'V', 'A', 'D', 'E', 'G' ]

def make_codon_dic ( codonList ):
    codonDic = {}
    for codon in codonList:
        codonDic [ codon ] = 0

    return codonDic

def make_dicodon_dic (codonList ):
    dicodonDic = {}
    for codon1 in codonList:
        for codon2 in codonList:
            dicodonDic [ codon1 + codon2 ] = 0

    return dicodonDic

def calc_codon_freq ( proteinsForAllSeq ):
    codonDic = make_codon_dic ( codonList )
    codonCount = 0

    for seq in proteinsForAllSeq:
        for protein in seq:
            for codon in protein:
                codonDic [ codon ] += 1
                codonCount += 1

    for codon in codonList:
        codonDic [ codon ] = codonDic [ codon ] / codonCount

    return codonDic

def calc_dicodon_freq ( proteinsForAllSeq ):
    dicodonDic = make_dicodon_dic ( codonList )
    dicodonCount = 0

    for seq in proteinsForAllSeq:
        for protein in seq:
            for idx, codon in enumerate ( protein ):
                if idx != ( len ( protein ) - 1 ) :
                    dicodonDic [ codon + protein [ idx + 1] ] += 1
                    dicodonCount += 1

    for ( key, value ) in dicodonDic.items():
        dicodonDic [ key ] = value / dicodonCount

    return dicodonDic

def euclidean_distance ( dic1, dic2 ):
    sum_sq_diff = 0.0

    for key in dic1.keys():
        sum_sq_diff += ( dic1 [ key ] - dic2 [ key ] ) ** 2

    return np.sqrt ( sum_sq_diff )

def distance_matrixes ( seqs ):
    n = len ( seqs.keys() )
    codonMatrix = np.zeros ( ( n, n ) )
    dicodonMatrix = np.zeros ( ( n, n ) )
    nameList = []

    for i, key1 in enumerate ( seqs.keys() ):

        ( codonFreq1, dicodonFreq1 ) = seqs [ key1 ]

        nameList.append ( key1 )

        for j, key2 in enumerate ( seqs.keys() ):

            ( codonFreq2, dicodonFreq2 ) = seqs [ key2 ]

            codonDist = euclidean_distance ( codonFreq1, codonFreq2 )
            dicodonDist = euclidean_distance ( dicodonFreq1, dicodonFreq2 )
            
            codonMatrix [ i ] [ j ] = codonDist
            dicodonMatrix [ i ] [ j ] = dicodonDist

    return ( nameList, codonMatrix, dicodonMatrix )

def save_phylip_format ( nameList, matrix, filename ):
    with open ( filename, 'w') as f:
        f.write ( f"{len(nameList)}\n" )
        for i, name in enumerate ( nameList ):
            f.write ( f"{name} " + " ".join ( f"{val}" for val in matrix [ i ] ) + "\n" )

files = [ 'bacterial1.fasta', 'bacterial2.fasta', 'bacterial3.fasta', 'bacterial4.fasta', 'mamalian1.fasta', 'mamalian2.fasta', 'mamalian3.fasta', 'mamalian4.fasta' ]

seqs = {}

for file in files:
    (name, seq) = parse_fasta ( file )

    splitSeqs = split_sequence ( seq )

    codons = turn_into_codon ( splitSeqs )

    proteinsForAllSplits = []
    for split in codons:
        proteins = find_proteins ( split )
        proteinsForAllSplits.append ( proteins )

    codonFreq = calc_codon_freq ( proteinsForAllSplits )
    dicodonFreq = calc_dicodon_freq ( proteinsForAllSplits )

    seqs [ name ] = ( codonFreq, dicodonFreq )

( nameList, codonMatrix, dicodonMatrix ) = distance_matrixes ( seqs )

save_phylip_format ( nameList, codonMatrix, 'codonMatrix.txt' )
save_phylip_format ( nameList, dicodonMatrix, 'dicodonMatrix.txt' )

#print ( distance_matrixes ( seqs ) )

#print ( seqs )


#(idd, rez) = parse_fasta('bacterial1.fasta')

#triplets = split_into_triplets(rez)

#print (type(rez))
#print ( reverse_compliment(rez) )

#print (triplets)

#finalSeq = split_sequence ( rez )

#for seq in finalSeq:
#    print (seq)
#    print ("\r\n")

#codons = turn_into_codon ( finalSeq )

#print ( codons )
#proteinsForAllSeq = []
#for seq in codons:
    #proteins = find_proteins ( seq )
    #proteinsForAllSeq.append ( proteins )


#print ( proteinsForAllSeq )
#dicodonDic = make_dicodon_dic ( codonList )

#codonFreq = calc_codon_freq ( proteinsForAllSeq )
#dicodonFreq = calc_dicodon_freq ( proteinsForAllSeq )

#print ( dicodonFreq )

#print ( euclidean_distance ( codonFreq, codonFreq ) )

