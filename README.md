## Kmer_fouder

"kmer_founder" is a small Python script that searches for the most over~~dosed~~represented kmers in a .fasta file.
kmer class has the following attributes:
1. sequence
2. counter
3. list of its coordinates (start and stop)

and the following methods:
1. increase counter
2. add new coordinate
3. show information about itself

## Kmer spectrum

"kmer_spectrum" is small python script that analyzes fastq file and builds kmer spectrum of defined k.
Special class "kmer_spectrum_builder" does all the operations of spectrum building step by step. It has the following methods:
1. analyse_file(k_length , quality_tres)
2. spectrum_data
3. draw_spectrum
4. assess_genome_size
5. clear_RAM

During "spectrum data" procedure it automatically establishes cut-off for noisy part of the spectrum, which is later drawn on the figure as a vertical line.


kmer_spectrum.py was tested on the "test_kmer.fastq".
Spectrums for k=17 are presented in the current repository.

In case of k=17 genome size estimated values were:
a) without quality filter

b)with quality filter (q>20)
Genome size without noise cut off: about 2,8 Mb
Corrected genome size with notoriuous noise cutted off 2,1 Mb

In case of k=31 with quality filter (q>20) the following values were obtained:
Genome size without noise cut off 3,2 Mb
Corrected genome size with notoriuous noise cutted off 2,0 Mb

So it can be seen that obtained estimations of genome size are quite small and considering it we can assume that it was some prokaryotic organism sequenced.








