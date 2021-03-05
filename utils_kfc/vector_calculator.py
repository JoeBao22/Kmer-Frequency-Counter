from utils_kfc.kmer_statistics import get_smaller, item_in_dict
from multiprocessing import Pool, Manager
from functools import partial
import itertools

def encoded_list(item, encoder):
    k, v = item
    return (k, encoder(k), v)


def encoded_list_freq(item, num_of_items, encoder):
    k, v = item
    return (k, encoder(k), v / num_of_items)
    

def subtraction_counter_task(k, all_counter, kmer_statistics):
    del_bgd = []
    for pre_ch in ["A", "C", "G", "T"]:
        for after_ch in ["A", "C", "G", "T"]:
            if kmer_statistics.combine:
                full_k = get_smaller(pre_ch + k + after_ch)
            else:
                full_k = pre_ch + k + after_ch
            k_sub1 = get_smaller(full_k[:-1])
            k_sub2 = get_smaller(full_k[1:])
            if all_counter[1][k_sub1] and all_counter[2][k_sub2]:
                prob = (all_counter[0][full_k] * all_counter[3][k]) / \
                (all_counter[1][k_sub1] * all_counter[2][k_sub2]) - 1
                encoded = kmer_statistics.encoder(full_k)
                del_bgd.append((full_k, encoded, prob))
    return del_bgd
    
def feature_vector(all_counter, kmer_statistics):
    """
    given the original copy number, generate a list of frequency or copy number. Use subtraction or not
    :param gross_f: a list of counters
    :return: [(key:str, encoded_key: int, value: number)]
    """
    if kmer_statistics.subtraction:
        """
        for str s, we get substrings:
        s1 = s[:-1] (remove the last character)
        s2: s[1:] (remove the first character)
        s3: s[1:-1] (remove the first and the last character)
        and we've got the frequency dicts for pattern s, s1, s2, s3 in freq, freq1, freq2, freq3, respectively
        the frequency dicts are grouped according to their patterns (positions being chosen in each substring)
        the background for key s is defined as :
        bgd(s) = freq_1(s1) * freq_2(s2) / freq_3(s3)
        we define the frequency with respect to background as :
        del_bgd(s) = freq(s) / bgd(s) - 1
        as the total number of items in each frequency dict are always the same, 
        del_bgd(s) = (counter(s) * counter3(s3)) / (counter1(s1) * counter2(s2)) - 1
        we will record the frequency of s as long as: s[1: ] in mid_dict
        """
        assert kmer_statistics.frequency
        with Pool(kmer_statistics.core) as pool: 
            task = partial(subtraction_counter_task, all_counter=all_counter, kmer_statistics=kmer_statistics)
            iterator_p = pool.imap_unordered(task, all_counter[3].keys(), chunksize=5000)
            del_bgd = list(itertools.chain(*iterator_p))  # flatten a list of list
            return del_bgd

    else:
        num_of_items = item_in_dict(all_counter[0])  # rename
        with Pool(kmer_statistics.core) as pool:
            if kmer_statistics.frequency:
                partial_res = partial(encoded_list,  encoder=kmer_statistics.encoder)
                iterator_p = pool.imap_unordered(partial_res, all_counter[0].items(), chunksize=5000)
            else:
                partial_res = partial(encoded_list_freq, num_of_items=num_of_items, encoder=kmer_statistics.encoder)
                iterator_p = pool.imap_unordered(partial_res, all_counter[0].items(), chunksize=5000)
            probability = sorted(list(iterator_p))
            return probability
