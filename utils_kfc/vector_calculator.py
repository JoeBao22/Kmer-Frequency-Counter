from utils_kfc.kmer_statistics import get_smaller, item_in_dict, reverse_complement
from multiprocessing import Pool, Manager, Process, Queue
from functools import partial
import itertools, time

def encoded_list(item, encoder):
    k, v = item
    return (k, encoder(k), v)


def encoded_list_freq(item, num_of_items, encoder):
    k, v = item
    return (k, encoder(k), v / num_of_items)
    

def encoded_list_use_fastq(item, encoder):
    k, v = item
    grouped_value_list = ["{}-{}-{}".format(x[0], x[1], x[2]) for x in v]
    # grouped_value = ",".join(map(lambda x: "{}-{}-{}".format(x[0], x[1], x[2]), v))
    return (k, encoder(k), len(v), ",".join(grouped_value_list))

def subtraction_counter_task(ks, all_counter, kmer_statistics, core_index):
    del_bgd = []
    filters = [("T", "A"), ("T", "C"), ("T", "G"), ("T", "T"),
            ("G", "T"), ("G", "G"), ("G", "C"), ("C", "T")]
    for k in ks:
        for pre_ch in ["A", "C", "G", "T"]:
            for after_ch in ["A", "C", "G", "T"]:
                if k == reverse_complement(k) and (pre_ch, after_ch) in filters:
                    continue
                if kmer_statistics.combine:
                    full_k = get_smaller(pre_ch + k + after_ch)
                    k_sub1 = get_smaller(full_k[:-1])
                    if not all_counter[1][k_sub1]:
                        continue
                    k_sub2 = get_smaller(full_k[1:])
                    if not all_counter[2][k_sub2]:
                        continue
                else:
                    full_k = pre_ch + k + after_ch
                    k_sub1 = full_k[:-1]
                    if not all_counter[1][k_sub1]:
                        continue
                    k_sub2 = full_k[1:]
                    if not all_counter[2][k_sub2]:
                        continue
                if not all_counter[0][full_k]:
                    prob = -1
                else:
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
        # with Pool(kmer_statistics.core) as pool: 
        #     task = partial(subtraction_counter_task, all_counter=all_counter, kmer_statistics=kmer_statistics)
        #     iterator_p = pool.imap_unordered(task, all_counter[3].keys(), chunksize=5000)
        #     del_bgd = list(itertools.chain(*iterator_p))  # flatten a list of list
        #     return del_bgd
        pool = Pool(kmer_statistics.core)
        all_counter_2 = []
        all_keys = list(all_counter[3].keys())
        keys_each = len(all_keys) // (kmer_statistics.core)
        all_tasks = []
        for core_index in range(kmer_statistics.core):
            if core_index == kmer_statistics.core - 1:
                temp_keys = all_keys[core_index * keys_each:]
            else:
                temp_keys = all_keys[core_index * keys_each: (core_index+ 1) * keys_each]
            all_tasks.append(temp_keys)
        for task in all_tasks:
            all_counter_2.append(pool.apply_async(subtraction_counter_task, args=(task, all_counter, kmer_statistics, core_index)) )
        pool.close()
        pool.join()

        del_bgd = []
        for i in all_counter_2:
            del_bgd = del_bgd + i.get()
        del_bgd.sort()
        
        return del_bgd
    else:
        num_of_items = item_in_dict(all_counter[0])  # rename
        with Pool(kmer_statistics.core) as pool:
            if kmer_statistics.frequency:
                partial_res = partial(encoded_list_freq, num_of_items=num_of_items, encoder=kmer_statistics.encoder)
                iterator_p = pool.imap_unordered(partial_res, all_counter[0].items(), chunksize=5000)
            else:
                partial_res = partial(encoded_list,  encoder=kmer_statistics.encoder)
                iterator_p = pool.imap_unordered(partial_res, all_counter[0].items(), chunksize=5000)
            probability = sorted(list(iterator_p))
            return probability
