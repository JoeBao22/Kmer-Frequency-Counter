import argparse,time
import os
import multiprocessing
from collections import Counter
from utils_kfc import io_manager
from utils_kfc.kmer_statistics import ArgumentManager, item_in_dict, reverse_complement, get_smaller, modify_sequence, no_deviation
from utils_kfc.vector_calculator import feature_vector


def check_parameters(args):
    """
    check whether the parameters are valid. If not, throw out warning and exit.
    :return: None
    """
    global logging
    global error

    if args.space and args.interval == 0:
        print("set space == True but interval not set. ")
        print("Please specify interval by --interval INTERVAL, where INTERVAL is an integer number > 0")
        error.write("set space == True but interval not set. \n")
        error.write(
            "Please specify interval by --interval INTERVAL, where INTERVAL is an integer number > 0\n")
        exit(1)
    if args.subtraction and not args.frequency:
        print("must use frequency instead of copy number when using subtraction")
        error.write("must use frequency instead of copy number when using subtraction\n")
        exit(1)
    if args.direction != 0 and args.position == 0:
        print("set shift direction == {} but position not chosen.".format(args.direction))
        error.write("set shift direction == {} but position not chosen.".format(args.direction))
        exit(1)
    if args.direction != 0 and args.position < 0:
        print("chosen position is exactly on or out of left boundary 1")
        error.write("chosen position is exactly on or out of left boundary 1\n")
        exit(1)
    if args.direction != 0 and args.position >= args.length-1:
        print("chosen position {} is exactly on or out of the right boundary {}".format(args.position+1, args.length))
        error.write("chosen position {} is exactly on or out of the right boundary {}\n".format(args.position+1, args.length))
        exit(1)
    if args.direction != 0 and args.combine:
        error.write("can't use position shift and combination at the same time\n")
        print("can't use position shift and combination at the same time")
        exit(1)
    if args.subtraction and args.length <= 2:
        error.write("length must be greater or equal to 3 when using subtraction\n")
        print("length must be greater or equal to 3 when using subtraction")
        exit(1)
    if args.use_fastq and (args.frequency):
        error.write("can not use subtraction when the input is *.fastq\n")
        print("can not use subtraction when the input is *.fastq")
        exit(1)


# def kmer_freq_update_fastq(seq_info_list, k):
#     record_dict = {}
#     for seq_info in seq_info_list:
#         file_num, line_num, seq = seq_info
#         for start_index in range(len(seq) - k + 1):
#             sub_seq = seq[start_index: start_index + k]
#             if no_deviation(sub_seq):
#                 record_dict[sub_seq] = record_dict.get(sub_seq, []) + [(file_num, line_num, start_index)]
#     return record_dict


def kmer_freq_update(seq, k, space, combine, loc):
    counter_list = [Counter()]
    modified_seq_value = modify_sequence(seq)
    for start_index in range(len(modified_seq_value) - k + 1):
        sub_seq = modified_seq_value[start_index: start_index + k]
        if space:
            sub_seq = ''.join([sub_seq[x] for x in loc])
        if combine:
            sub_seq = get_smaller(sub_seq)
        if no_deviation(sub_seq):
            counter_list[0][sub_seq] += 1
    return counter_list


def kmer_freq_subtraction_update(seq, k, space, combine, loc):
    counter_list = [Counter() for _ in range(4)]
    modified_seq_value = modify_sequence(seq)
    for start_index in range(len(modified_seq_value) - k + 1):
        sub_seq = modified_seq_value[start_index: start_index + k]
        if space:
            sub_seq = ''.join([sub_seq[x] for x in loc])
        if no_deviation(sub_seq):
            if combine:
                sub_seq = get_smaller(sub_seq)
            mid_seq1, mid_seq2 = sub_seq[: -1], sub_seq[1: ]
            short_seq = sub_seq[1:-1]
            if combine:
                mid_seq1 = get_smaller(mid_seq1)
                mid_seq2 = get_smaller(mid_seq2)
                short_seq = get_smaller(short_seq)
            seq_list = [sub_seq, mid_seq1, mid_seq2, short_seq]
            for i in range(4):
                counter_list[i][seq_list[i]] += 1
    return counter_list
    
    
def sum_counter(counter_list):
    summed = Counter()
    for counter in counter_list:
        summed += counter
    return summed

def merge_counter(all_res_list, core):
    all_counter_list = [res.get() for res in all_res_list]
    print(len(all_counter_list))
    num_of_counter = len(all_counter_list[0])
    print("num of counter", num_of_counter)
    merged_counter = []
    for i in range(num_of_counter):
        res = []
        pool =  multiprocessing.Pool(core)
        for c in range(core):
            chosen_counter = [all_counter_list[j][i] for j in range(c, len(all_counter_list), core)]
            # start_index = c * (len(all_counter_list) // core)
            # end_index = c *  min(len(all_counter_list) // core + 1, len(all_counter_list))
            # chosen_counter = [all_counter_list[j][i] for j in range(start_index, end_index)]
            res.append(pool.apply_async(sum_counter, args=(chosen_counter, ))) 
        pool.close()
        pool.join()
        merged_counter.append(sum([item.get() for item in res], Counter()))
    return merged_counter
    
def merge_record(all_record_list):
    print(all_record_list)
    all_record_dict = [res.get() for res in all_record_list]
    merged_record = {}
    for record in all_record_dict:
        # print(record)
        for k, v in record.items():
            merged_record[k] = merged_record.get(k, []) + v
    print(merged_record)
    return merged_record
    

def main():
    parser = argparse.ArgumentParser(description="A parser for KFC that accepts input from the user")
    parser.add_argument('--file_dir', type=str, required=True,
                        help="the directory that contains all the dna-sequence files")
    parser.add_argument('--name_file', type=str, required=True,
                        help="the path to the file that contains all the names. The names don't have any suffix like .fa")
    parser.add_argument('--core', type=int, required=True,
                        help="the number of CPU be used")
    parser.add_argument('--length', type=int, required=True,
                        help="the length of spaced word, N>1")
    parser.add_argument('--new_encoder', action="store_true",
                        help="whether to use the new encoder-decoder with base 5 or " \
                        "the original one with base 4")
    parser.add_argument("--subtraction", action="store_true",
                        help="whether to do background subtraction or not")
    parser.add_argument("--frequency", action="store_true",
                        help="whether to use the original copy number or the frequency")
    parser.add_argument('--combine', action="store_true",
                        help="whether to combine the result for s with that for reverse_complement(s)")
    parser.add_argument("--space", action="store_true",
                        help="whether the original sequence is spaced or not")
    parser.add_argument('--interval', type=int, default=0,
                        help="the interval length of degenerated Kmer, sample every (d+1) positions. d>0 if spaced")
    parser.add_argument('--direction', type=int, default=0,
                        help="shift of one position in pattern. "
                            "0: None, 1: to the right of the original sample position, -1: to the left")
    parser.add_argument('--position', type=int, default=0,
                        help="the position of shift, count from 1. "
                            "eg: 1, 2, 3, 4...9 ; sampled: 1, 3, 5, 7, 9"
                            "if position = 5 and direction = 1:"
                            "shifted: 1, 3, 6, 7, 9")
    parser.add_argument('--use_fastq', action="store_true",
                        help="whether the input files are in the format of fastq")

    t = time.time()
    args = parser.parse_args()
    names, output_folder, logging, error = io_manager.output_error_logging_config(args)
    check_parameters(args)
    print("------Arguments:------\n{}".format(args))
    logging.write("------ARGUMENTS:------\n{}\n".format(args))
    kmer_statistics = ArgumentManager(args)
    for name in names:
        print("calculating {}".format(name))
        logging.write("CALCULATING {}\n".format(name.upper()))
        if kmer_statistics.use_fastq:
            suffix = "fq"
            file_path_upstream = os.path.join(kmer_statistics.file_dir, "{}_{}.{}".format(name, "R1", suffix))
            if not os.path.exists(file_path_upstream):
                error.write("------upstream file not found for {} ------\n".format(name))
                error.write("the upstream file should be like 'apple_R1.fq', if the species name is 'apple'")
            file_path_downstream = os.path.join(kmer_statistics.file_dir, "{}_{}.{}".format(name, "R2", suffix))
            if not os.path.exists(file_path_downstream):
                error.write("------downstream file not found for {} ------\n".format(name))
                error.write("the downstream file should be like 'apple_R2.fq', if the species name is 'apple'")
            logging.write("------extracting sequences from path {} and path {} ------\n".format(
                file_path_upstream, file_path_downstream))
            seq_list = io_manager.extract_fastq(file_path_upstream, file_path_downstream)
        else:
            all_file_path = [os.path.join(kmer_statistics.file_dir, "{}.{}".format(name, suffix)) for suffix in ["fa","fna","ffn","fsa","fasta"]]
            chosen_file_path = None  # choose the first file that satisfies certain format
            for file_path in all_file_path:
                if os.path.exists(file_path):
                    chosen_file_path = file_path
                    break
            if not chosen_file_path:
                error.write("------no file found for {} ------\n".format(name))
                exit(1)
            logging.write("------extracting sequences from the given path {} ------\n".format(chosen_file_path))

            seq_list = io_manager.extract_fasta(chosen_file_path)
            
        # print("seq_list", seq_list)
        pool = multiprocessing.Pool(kmer_statistics.core)
        if kmer_statistics.subtraction:
            function_to_call = kmer_freq_subtraction_update
        else:
            function_to_call = kmer_freq_update
        segmentation_len = 5000
        # the ACTUAL segmentLen will be segmentLen + k - 1 because of overlapping.
        all_tasks = []
        all_counter = []
        for seq_value in seq_list:
            for start_index in range(0, len(seq_value), segmentation_len):
                temp_seq = seq_value[max(0, start_index - kmer_statistics.k + 1): \
                    min(start_index + segmentation_len, len(seq_value))]
                all_tasks.append(temp_seq)
        # print("all tasks:", all_tasks)
        for task in all_tasks:
            all_counter.append(
                    pool.apply_async(function_to_call, 
                                args=(task,
                                    kmer_statistics.k, kmer_statistics.space, 
                                    kmer_statistics.combine, kmer_statistics.loc))
                )
        pool.close()
        pool.join()
        print('input kmer dict finished!',time.time()-t)
        logging.write('input kmer dict finished!' + str(time.time()-t) + "\n")
        t=time.time()

        freq_list_return = merge_counter(all_counter, kmer_statistics.core)
        print('merge counter finished!',time.time()-t)
        
        logging.write('merge counter finished!' + str(time.time()-t) + "\n")
        t=time.time()
        
        vector_frequency = feature_vector(freq_list_return, kmer_statistics)
        print('calculate vector frequency finished!' + str(time.time()-t) + "\n")
        t=time.time()

        if kmer_statistics.subtraction:
            logging.write("long dict: \n")
            logging.write("{:<.50}...\n".format(str(freq_list_return[0])))
            logging.write("mid dict 1: \n")
            logging.write("{:<.50}...\n".format(str(freq_list_return[1])))
            logging.write("mid dict 2: \n")
            logging.write("{:<.50}...\n".format(str(freq_list_return[2])))
            logging.write("short dict: \n")
            logging.write("{:<.50}...\n".format(str(freq_list_return[3])))
        else:
            logging.write("{:<.50}...\n".format(str(freq_list_return[0])))

        for item in vector_frequency[:min(4, len(vector_frequency))]:
            logging.write("key: {0}\tencoded key: {1}\tvalue: {2}\n".format(
                item[0], item[1], item[2]))

        logging.write("------saving result for {}: ------\n".format(name))
        io_manager.output_content(kmer_statistics, output_folder, vector_frequency, name)

        logging.write("\n")
        print('all finished!',time.time()-t)
        logging.write('all finished!'+ str(time.time()-t) + "\n")
    logging.close()

    exit(time.time()-t)
    
    
if __name__ == '__main__':
    main()
