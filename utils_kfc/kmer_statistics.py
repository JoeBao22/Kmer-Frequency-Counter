from .encoder_decoder_old import EncoderDecoderOld
from .encoder_decoder import EncoderDecoder


def item_in_dict(d):
    counter = 0
    for v in d.values():
        counter += v
    return counter


def reverse_complement(s):
    """
    Obtain Inverted complement string of a DNA sequence
    :param s: string, sequence
    :return: the reversed complement string of the original DNA sequence
    """
    base_complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    letters = [base_complement[character] if character in base_complement.keys() else "N" for character in s]
    return ''.join(letters)[::-1]


def get_smaller(s):
    """
    Return the smaller one of 1) original sequence 2) its reverse complement
    :param s: string, sequence
    :return: return the smaller one of the original sequence and its reverse complement
    """
    reverse_complement_s = reverse_complement(s)
    return min(s, reverse_complement_s)


def modify_sequence(s):
    """
    concatenate the reverse-complement string with the original one with character N.
    N can natually separate the two string when we count
    :param s: str, the original sequence
    :return: str
    """
    modified_s = s + "N" + reverse_complement(s)
    return modified_s


def no_deviation(kmer):
    """
    determine whether the sequence, kmer, contains any characters other than A, T, C, G
    :param kmer: str
    :return: bool, True if the sequence contains no characters other than A, T, C, G
    """
    return kmer.count('A') + kmer.count('T') + kmer.count('C') + kmer.count('G') == len(kmer)


def vector_calculator(gross_f, kmer_statistics):
    """
    given the original copy number, generate a list of frequency or copy number. Use subtraction or not
    :param gross_f: {key: value}
    :return: [(key:str, encoded_key: int, value: number)]
    """
    if kmer_statistics.subtraction:
        # for str s, we get substrings:
        # s1 = s[:-1] (remove the last character)
        # s2: s[1:] (remove the first character)
        # s3: s[1:-1] (remove the first and the last character)
        # and we've got the frequency dicts for pattern s, s1, s2, s3 in freq, freq1, freq2, freq3, respectively
        # the frequency dicts are grouped according to their patterns (positions being chosen in each substring)
        # the background for key s is defined as :
        # bgd(s) = freq_1(s1) * freq_2(s2) / freq_3(s3)
        # we define the frequency with respect to background as :
        # del_bgd(s) = freq(s) / bgd(s) - 1
        # we will record the frequency of s as long as: s[1: ] in mid_dict
        probability = [{} for _ in range(4)]
        bgd = {}
        del_bgd = {}
        if kmer_statistics.frequency:
            for index, freq_count in enumerate(gross_f):
                num_of_items = item_in_dict(freq_count)
                for (k, v) in freq_count.items():
                    probability[index][kmer_statistics.e_d.encoder(k)] = v / num_of_items
        else:
            # it will not the branch, because we require the use of frequency when using subtraction
            for index, freq_count in enumerate(gross_f):
                for (k, v) in freq_count.items():
                    probability[index][kmer_statistics.e_d.encoder(k)] = v
        for (k, v) in gross_f[3].items():
            for pre_ch in ["A", "C", "G", "T"]:
                for after_ch in ["A", "C", "G", "T"]:
                    full_k = get_smaller(pre_ch + k + after_ch)
                    k_sub1 = get_smaller(pre_ch + k)
                    k_sub2 = get_smaller(k + after_ch)
                    k_sub3 = k
                    if k_sub1 in gross_f[1].keys() and  \
                        k_sub2 in gross_f[2].keys():
                        encoded_0, encoded_1, encoded_2, encoded_3 = kmer_statistics.e_d.encoder(full_k), kmer_statistics.e_d.encoder(
                            k_sub1), kmer_statistics.e_d.encoder(k_sub2), kmer_statistics.e_d.encoder(k_sub3)
                        bgd[encoded_0] = probability[1][encoded_1] * probability[2][encoded_2] / probability[3][
                            encoded_3]
                        del_bgd[full_k] = ((full_k, encoded_0, probability[0].get(encoded_0, 0) / bgd[encoded_0] - 1))
        del_bgd_list = sorted(list(del_bgd.values()))
        return del_bgd_list

    else:
        probability = []
        num_of_items = item_in_dict(gross_f[0])  # rename
        if kmer_statistics.frequency:
            for (k, v) in gross_f[0].items():
                probability.append((k, kmer_statistics.e_d.encoder(k), v / num_of_items))
        else:
            for (k, v) in gross_f[0].items():
                probability.append((k, kmer_statistics.e_d.encoder(k), v))
        probability.sort()
        return probability


class ArgumentManager:
    """
    Some functions related to the calculation of Kmer Statistics
    """
    def __init__(self, args):
        self.combine = args.combine
        self.core = args.core
        self.d = args.interval  # sample every (d+1) positions
        self.direction = args.direction
        self.file_dir = args.file_dir
        self.frequency = args.frequency
        self.n = args.length  # length of the spaced sequence
        self.new_encoder = args.new_encoder
        self.position = args.position - 1
        self.space = args.space
        self.subtraction = args.subtraction
        
        self.k = self.n + (self.n - 1) * self.d
        if self.new_encoder:
            self.e_d = EncoderDecoder()
        else:
            self.e_d = EncoderDecoderOld()
        self.loc = list(range(0, self.k, self.d + 1))  # the chosen positions when spaced
        self.loc[self.position] += self.direction






