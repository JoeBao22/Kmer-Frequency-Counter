from utils_kfc.kmer_statistics import get_smaller, item_in_dict
from multiprocessing import Pool
from functools import partial


def encoded_list(item, encoder):
    k, v = item
    return (k, encoder(k), v)


def encoded_list_freq(item, num_of_items, encoder):
    k, v = item
    return (k, encoder(k), v / num_of_items)
    
def feature_vector(gross_f, kmer_statistics):
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
        assert kmer_statistics.frequency
        for index, freq_count in enumerate(gross_f):
            num_of_items = item_in_dict(freq_count)
            for (k, v) in freq_count.items():
                probability[index][kmer_statistics.encoder(k)] = v / num_of_items
        for (k, v) in gross_f[3].items():
            for pre_ch in ["A", "C", "G", "T"]:
                for after_ch in ["A", "C", "G", "T"]:
                    full_k = get_smaller(pre_ch + k + after_ch)
                    k_sub1 = get_smaller(pre_ch + k)
                    k_sub2 = get_smaller(k + after_ch)
                    k_sub3 = k
                    if gross_f[1][k_sub1] and gross_f[2][k_sub2]:
                        encoded_0, encoded_1, encoded_2, encoded_3 = kmer_statistics.encoder(full_k), kmer_statistics.encoder(
                            k_sub1), kmer_statistics.encoder(k_sub2), kmer_statistics.encoder(k_sub3)
                        bgd[encoded_0] = probability[1][encoded_1] * probability[2][encoded_2] / probability[3][
                            encoded_3]
                        del_bgd[full_k] = ((full_k, encoded_0, probability[0].get(encoded_0, 0) / bgd[encoded_0] - 1))
        del_bgd_list = sorted(list(del_bgd.values()))
        return del_bgd_list

    else:
        num_of_items = item_in_dict(gross_f[0])  # rename
        with Pool(kmer_statistics.core) as pool:
            if kmer_statistics.frequency:
                partial_res = partial(encoded_list,  encoder=kmer_statistics.encoder)
                iterator_p = pool.imap_unordered(partial_res, gross_f[0].items(), chunksize=5000)
            else:
                partial_res = partial(encoded_list_freq, num_of_items=num_of_items, encoder=kmer_statistics.encoder)
                iterator_p = pool.imap_unordered(partial_res, gross_f[0].items(), chunksize=5000)
            probability = sorted(list(iterator_p))
            return probability
