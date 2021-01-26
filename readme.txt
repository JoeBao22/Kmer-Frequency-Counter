=========================================
= KmerFrequencyCalculator version 0.3.0 =
=========================================

by Xiaoyi Bao in SJTU & Qi WU* in IM, CAS.

Bao: aaa@sjtu.edu.cn
WU: wuq@im.ac.cn

*: WU is the correspodence author.


LICENCE
=======

KFC is free software which may be redistributed and modified 
under the terms of the GNU General Public License 
...


INTRODUCTION
============

KmerFrequencyCalculator (KFC) is a software package for alignment-free
sequence analyses. It calculate the copy number or frequencies of oligo-
nucleotides with given length from target DNA sequence as the global
statistical feature. The frequencies or copy numbers can be used for further
comparisons or other analyses between target sequences.


DESCRIPTION
===========

kfc.py:

Takes a fasta sequence file or a number of fasta files in one directory to
perform certain kmer computation defined by a series of parameters.

Options:
      -h                    prints help information.
      --name_file           input list file name, contains sequences' name
      --file_dir            input sequences' dir 
      --core                multiple processes or CPU cores to be used    
      --length              kmer length
      --new_encoder         using new codes for nucleotides with 1, 2, 3, 4
                            representing A, C, G, T respectively. Default
                            setting using codes in cvtree(Qi et al. 2004), in
                            which A, C, G and T are represented as 0, 1, 2
                            and 3. when using new codes, the kmer length need
                            not be given for decoding process. 
      --subtraction         Using random background subtraction which is the
                            same as that in cvtree.
      --frequency           Computing the frequency instead of counting the
                            copy number. Default is just counting the copy
                            number.
      --combine             combining the kmer with the same position in
                            reverse-complement DNA chain. Note that this
                            parameter is necessary for the usage of spaced-
                            word with uneven pattern .
      --space               Using spaced-word method (Leimester et al. 2014).
      --interval INTERVAL   Defining the interval length for spaced-word with
                            even pattern, in which all the "match" positions
                            are divided by intervals with same length, which
                            is consist of "don't care" positions.
      --direction           This parameter and the following --position
                            parameter define the uneven pattern of spaced-
                            word. The uneven patterns are those spaced
                            patterns most similar with the even pattern. That
                            means, of all the match position excluding the
                            begin and end position, only one position may
                            deviate left or right one base position. This
                            parameter defines the deviation direction, -1
                            means left direction of the kmer (or 5-'
                            direction in DNA sequence. While 1 means right
                            direction of the kmer or 3-' direction in DNA
                            sequence.
      --position POSITION   This parameter defines which of the "match
                            postion" is chosed to deviate to form an uneven
                            pattern.

Programs in the utils_kfc directory are provided for developers or advanced
users.

utils_kfc/encoder_decoder.py

Encoding a kmer of sequence to an integral number, or decoding a number to
the sequence of that kmer using number of 1, 2, 3 and 4 to represent A, C, G
and T as described in the paramter of --new_encoder in kfc.py program.

Methods:
      --encoder_decoder.EncoderDecoder.encoder(target_kmer_string)
      --encoder_decoder.EncoderDecoder.decoder(target_kmer_code)

utils_kfc/encoder_decoder_old.py

Encoding/decoding a kmer using number of 0, 1, 2 and 3 to represent A, C, G
and T, which is the same as that used in cvtree.

Methods:
      --encoder_decoder.EncoderDecoder.encoder(kmer_string)
      --encoder_decoder.EncoderDecoder.decoder(kmer_code, kmer_length)


utils_kfc/io_manager.py

There are some classes and methods used for input/output treatment.

utils_kfc/kmer_statistics.py

There are a series of classes and functions in this program, which
are used in the main program kfc.py.

References:
-----------
Leimeister, C.-A., Boden, M., Horwege, S., Lindner, S., & Morgenstern, B.
(2014). Fast alignment-free sequence comparison using spaced-word
frequencies. BIOINFORMATICS, 30(14), 1991–1999.
Qi, J., Luo, H., & Hao, B. (2004). CVTree: A phylogenetic tree
reconstruction tool based on whole genomes. JOURNAL OF MOLECULAR EVOLUTION,
58(1), 1–11.


USAGES AND EXAMPLES
===================

1: basic
python kfc.py --file_dir /input/folder/ --core 2 --length 4 --name_file \
names.txt

2: subtraction with frequency
python kfc.py --file_dir /input/folder/  --core 2 --length 4 --name_file \
names.txt --subtraction --frequency

3: combine reverse-complement kmers
python kfc.py --file_dir /input/folder/  --core 2 --length 2 --name_file \
names.txt --combine

4: spaced-word with pattern of even pattern 
python kfc.py --file_dir /input/folder/  --core 2 --length 3 --name_file \
names.txt --space --interval 1

5: spaced-word direction with pattern of uneven pattern
python kfc.py --file_dir /input/folder/  --core 2 --length 3 --name_file \
names.txt --space --interval 1 --direction 1

python kfc.py --file_dir /input/folder/  --core 2 --length 3 --name_file \
names.txt --space --interval 1 --direction -1


CITATIONS
=========

Jiuhong Dong#, Shuai Liu#, Yaran Zhang, Yi Dai and Qi Wu (2020) A New
Alignment-Free Whole Metagenome Comparison Tool and Its Application on Gut
Microbiomes of Wild Giant Pandas. Frontiers in Microbiology. 11:1061.


BUGS
====

please email a@a.edu.cn or wuq@im.ac.cn if you find one!


CHANGE LOG
==========

0.3.0:
* Objective Oriented programming.
* revised random background algrithum.
* uneven pattern of spaced-word. In this pattern the reverse-complement kmers
  should not be combined.
  This is the reason of adding --combine parameter.
* integrity of kmer copy number computation and kmer frequency computation. 

0.2.1:
Distributed in github (https://github.com/LDS-axel/KmerFreqCalc) written by
Shuai Liu supervised by Qi WU.

