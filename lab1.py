from Bio import SeqIO

def parse_fasta(file):
    for record in SeqIO.parse(file, "fasta"):
        return str(record.seq).upper()

def split_sequence(sequence):
    sequences = []

    for i in range ( 0, 3 ):
        sequences.append ( split_into_triplets ( sequence [ i: ] ) )

    #for i in range ( 0, 3 ):
     #   sequneces.add ( split_intotriplets (  )

    return sequences

def split_into_triplets(sequence):
    triplets = []
    for i in range ( 0 , len(sequence), 3 ):
        triplets.append ( sequence [ i : i + 3 ] )
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



rez = parse_fasta('/home/ignasvoveriukas/VU/Bioinformatika/lab1/BioinformatikaLab1/bacterial1.fasta')

triplets = split_into_triplets(rez)

#print (rez)
#print (type(rez))
#print ( reverse_compliment(rez) )

print (triplets)
