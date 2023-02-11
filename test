import sys

program_name = sys.argv[0]
human_fasta = sys.argv[1]
mouse_fasta = sys.argv[2]
kmer_size = int(sys.argv[3])
transcript_num = int(sys.argv[4])


def parse_fasta_file(file_name):
    transcripts = {}
    transcript = ""
    transcript_name = ""
    with open(file_name, "r") as f:
        for line in f:
            if line[0] == ">":
                if transcript_name:
                    transcripts[transcript_name] = transcript
                    transcript = ""
                transcript_name = line.strip()[1:]
            else:
                transcript += line.strip()
    transcripts[transcript_name] = transcript
    return transcripts


def kmer_counts_fasta(human_fasta_file, mouse_fasta_file, k, transcript_num):
    human_kmers = []
    mouse_kmers = []
    shared_kmers = []
    shared_kmer_positions = []

    human_transcripts = parse_fasta_file(human_fasta_file)
    mouse_transcripts = parse_fasta_file(mouse_fasta_file)

    human_transcript = human_transcripts[list(human_transcripts.keys())[transcript_num]]
    mouse_transcript = mouse_transcripts[list(mouse_transcripts.keys())[transcript_num]]

    human_kmer_set = set()
    for i in range(len(human_transcript) - k + 1):
        kmer = human_transcript[i:i + k]
        human_kmer_set.add(kmer)
    human_kmers.append(human_kmer_set)

    mouse_kmer_set = set()
    for i in range(len(mouse_transcript) - k + 1):
        kmer = mouse_transcript[i:i + k]
        mouse_kmer_set.add(kmer)
    mouse_kmers.append(mouse_kmer_set)

    shared_kmer_set = human_kmer_set.intersection(mouse_kmer_set)
    shared_kmers.append(len(shared_kmer_set))

    for kmer in shared_kmer_set:
        human_positions = [i for i, x in enumerate(human_transcript) if human_transcript.startswith(kmer, i)]
        mouse_positions = [i for i, x in enumerate(mouse_transcript) if mouse_transcript.startswith(kmer, i)]
        shared_kmer_positions.append((kmer, human_positions, mouse_positions))

    return [len(kmer_set) for kmer_set in human_kmers], [len(kmer_set) for kmer_set in
                                                         mouse_kmers], shared_kmers, shared_kmer_positions


num_human_kmer, num_mouse_kmer, num_shared_kmer, shared_kmer_pos = kmer_counts_fasta(human_fasta, mouse_fasta, kmer_size, transcript_num)
print(num_human_kmer, num_mouse_kmer, num_shared_kmer) ## add "shared_kmer_pos" to see the kmers and their respective match positions in the human and mouse transcripts

