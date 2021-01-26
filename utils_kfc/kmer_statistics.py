from .encoder_decoder_old import EncoderDecoderOld
from .encoder_decoder import EncoderDecoder


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

class KmerCount:
    """
    Count the frequency of each subsequence, and save in a dict
    """
    def __init__(self, combine=True):
        """
        :param combine: bool.
        If combine is True, then s and reverse_complement(s) share the same key when the frequency is counted
        """
        self.freq = {}
        self.combine = combine

    def freq_update(self, new_kmer):
        """
        update the kmer counter with current kmer, a (maybe sampled) substring of the original sequence
        the smaller one of two keys new_kmer, reverse_complement(new_kmer) is used as the key in the dict
        :param new_kmer: str
        :return: None
        """
        if self.combine:
            new_kmer = get_smaller(new_kmer)
        self.freq[new_kmer] = self.freq.get(new_kmer, 0) + 1

    def get_freq(self):
        """
        return the frequency dict
        :return: dict{str: number}
        """
        return self.freq


class FrequencyDict:
    """
    Combine the result of several KmerCount, and return a gross dict
    """
    def __init__(self):
        self.frequency_dict = {}
        self.num_of_items = 0

    def integrate(self, new_counter):
        """
        integrate the new counter into the temperate frequency dict
        :param new_counter: KmerCount
        :return: None
        """
        new_dict = new_counter.get_freq()
        for k, v in new_dict.items():
            self.num_of_items += new_dict[k]
            self.frequency_dict[k] = self.frequency_dict.get(k, 0) + new_dict[k]

    def get_frequency_dict(self):
        """
        return the frequency dict stored within
        :return: dict{str: number}
        """
        return self.frequency_dict


class KmerStatistics:
    """
    Some functions related to the calculation of Kmer Statistics
    """
    def __init__(self, args, logging_fd, error_fd):
        self.d = args.interval  # sample every (d+1) positions
        self.n = args.length  # length of the spaced sequence
        self.subtraction = args.subtraction
        self.space = args.space
        self.frequency = args.frequency
        self.new_encoder = args.new_encoder
        if self.new_encoder:
            self.e_d = EncoderDecoder()
        else:
            self.e_d = EncoderDecoderOld()
        self.combine = args.combine
        self.k = self.n + (self.n - 1) * self.d
        self.direction = args.direction
        self.position = args.position - 1
        self.loc = list(range(0, self.k, self.d + 1))  # the chosen positions when spaced
        self.loc[self.position] += self.direction

    def get_freq_function(self, seq_dict):
        """
        return a function to get frequency dict from the original sequence dict
        :param seq_dict:
        :return: a function used to build kmer frequency dict from the original sequence
        """
        if self.subtraction:
            return self.__kmer_freq_subtraction(seq_dict, spaced=self.space)
        else:
            return self.__kmer_freq(seq_dict, spaced=self.space)

    @staticmethod
    def __modify_sequence(s):
        """
        concatenate the reverse-complement string with the original one with character N.
        N can natually separate the two string when we count
        :param s: str, the original sequence
        :return: str
        """
        modified_s = s + "N" + reverse_complement(s)
        return modified_s

    @staticmethod
    def __no_deviation(kmer):
        """
        determine whether the sequence, kmer, contains any characters other than A, T, C, G
        :param kmer: str
        :return: bool, True if the sequence contains no characters other than A, T, C, G
        """
        return kmer.count('A') + kmer.count('T') + kmer.count('C') + kmer.count('G') == len(kmer)

    def gross_frequency(self, frequency_count):
        """
        combine the result on different threads
        :param frequency_count: [[KmerCount]]
        :return: [FrequencyDict]
        """
        if self.subtraction:
            freq_list_return = [FrequencyDict() for _ in range(4)]  # long, mid1, mid2, short
            for thread in frequency_count:
                for index, type_item in enumerate(thread):
                    freq_list_return[index].integrate(type_item)
        else:
            freq_list_return = [FrequencyDict()]
            for thread in frequency_count:
                item = thread[0]
                freq_list_return[0].integrate(item)
        return freq_list_return

    def __kmer_freq(self, seq_dict, spaced=False):
        """
        Designed for Kmer counting without subtraction of background
        :param seq_dict: dict
        :param spaced:  whether to sample with certain pattern or not
        :return: [KmerCount]
        """
        k_count = KmerCount(self.combine)
        for (seq_key, seq_value) in seq_dict.items():
            modified_seq_value = self.__modify_sequence(seq_value)
            for start_index in range(len(modified_seq_value) - self.k + 1):
                sub_seq = modified_seq_value[start_index: start_index + self.k]
                if spaced:
                    sub_seq = ''.join([sub_seq[x] for x in self.loc])
                if self.__no_deviation(sub_seq):
                    k_count.freq_update(sub_seq)  # update the counter only when the sub sequence is valid
        return [k_count]

    def __kmer_freq_subtraction(self, seq_dict, spaced=False):
        """
        Designed for Kmer counting with subtraction of background
        :param seq_dict:
        :param spaced: whether to sample with certain pattern or not
        :return: [KmerCount]
        """
        k_count_long, k_count_mid1, k_count_mid2, k_count_short = \
            KmerCount(self.combine), KmerCount(self.combine), KmerCount(self.combine), KmerCount(self.combine)
        for seq_key, seq_value in seq_dict.items():
            modified_seq_value = self.__modify_sequence(seq_value)
            for start_index in range(len(modified_seq_value) - self.k + 1):
                sub_seq = modified_seq_value[start_index: start_index + self.k]
                if spaced:
                    sub_seq = ''.join([sub_seq[x] for x in self.loc])
                if self.__no_deviation(sub_seq):
                    # subsequence of length k
                    k_count_long.freq_update(sub_seq)
                    mid_seq1, mid_seq2 = sub_seq[: -1], sub_seq[1: ]
                    k_count_mid1.freq_update(mid_seq1)
                    k_count_mid2.freq_update(mid_seq2)
                    short_seq = sub_seq[1:-1]
                    k_count_short.freq_update(short_seq)

        return [k_count_long, k_count_mid1, k_count_mid2, k_count_short]

    def vector_calculator(self, gross_f):
        """
        given the original copy number, generate a list of frequency or copy number. Use subtraction or not
        :param gross_f: {key: value}
        :return: [(key:str, encoded_key: int, value: number)]
        """
        if self.subtraction:
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
            if self.frequency:
                for index, freq_count in enumerate(gross_f):
                    num_of_items = freq_count.num_of_items
                    for (k, v) in freq_count.get_frequency_dict().items():
                        probability[index][self.e_d.encoder(k)] = v / num_of_items
            else:
                # it will not the branch, because we require the use of frequency when using subtraction
                for index, freq_count in enumerate(gross_f):
                    for (k, v) in freq_count.get_frequency_dict().items():
                        probability[index][self.e_d.encoder(k)] = v
            for (k, v) in gross_f[3].get_frequency_dict().items():
                for pre_ch in ["A", "C", "G", "T"]:
                    for after_ch in ["A", "C", "G", "T"]:
                        full_k = get_smaller(pre_ch + k + after_ch)
                        k_sub1 = get_smaller(pre_ch + k)
                        k_sub2 = get_smaller(k + after_ch)
                        k_sub3 = k
                        if k_sub1 in gross_f[1].get_frequency_dict().keys() and  \
                            k_sub2 in gross_f[2].get_frequency_dict().keys():
                            encoded_0, encoded_1, encoded_2, encoded_3 = self.e_d.encoder(full_k), self.e_d.encoder(
                                k_sub1), self.e_d.encoder(k_sub2), self.e_d.encoder(k_sub3)
                            bgd[encoded_0] = probability[1][encoded_1] * probability[2][encoded_2] / probability[3][
                                encoded_3]
                            del_bgd[full_k] = ((full_k, encoded_0, probability[0].get(encoded_0, 0) / bgd[encoded_0] - 1))
            del_bgd_list = sorted(list(del_bgd.values()))
            return del_bgd_list

        else:
            probability = []
            num_of_items = gross_f[0].num_of_items  # rename
            if self.frequency:
                for (k, v) in gross_f[0].get_frequency_dict().items():
                    probability.append((k, self.e_d.encoder(k), v / num_of_items))
            else:
                for (k, v) in gross_f[0].get_frequency_dict().items():
                    probability.append((k, self.e_d.encoder(k), v))
            probability.sort()
            return probability







