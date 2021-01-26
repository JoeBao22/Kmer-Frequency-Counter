import argparse
import os
import multiprocessing
from utils_kfc import io_manager, kmer_statistics


# update 0920:
# rename for output folder
# check parameter before calculating


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


def calculate(args, name, output_folder, kmer_statistics):
    global logging
    global error
    all_file_path = [os.path.join(args.file_dir, "{}.{}".format(name, suffix)) for suffix in ["fa","fna","ffn","fsa","fasta"]]
    chosen_file_path = None  # choose the first file that satisfies certain format
    for file_path in all_file_path:
        if os.path.exists(file_path):
            chosen_file_path = file_path
            break
    if not chosen_file_path:
        error.write("------no file found for {} ------\n".format(name))
        return
    logging.write("------extracting sequences from the given path {} ------\n".format(chosen_file_path))

    seq_dict = io_manager.extract_dict(chosen_file_path)
    queue = [{} for _ in range(args.core)]  # allocate sequences to different cores
    thread_pool = multiprocessing.Pool(args.core)
    thread_counter = 0
    for seq_key, seq_value in seq_dict.items():
        chosen_core = thread_counter % args.core
        queue[chosen_core][seq_key] = seq_value
        thread_counter += 1
    # logging.write("------tasks on each core:------\n")
    # for i in range(args.core):
    #     logging.write("core #{}\n".format(i + 1))
    #     logging.write(str(queue[i]) + "\n")
    kmer_statistics_function = kmer_statistics.get_freq_function
    freq_count = thread_pool.map(kmer_statistics_function, queue)
    logging.write("------counting result on each core------\n")
    for i in range(args.core):
        logging.write("core #{}\n".format(i + 1))
        if args.subtraction:
            logging.write("long dict: \n")
            logging.write("{:<.50}...\n".format(str(freq_count[i][0].get_freq())))
            logging.write("mid dict 1: \n")
            logging.write("{:<.50}...\n".format(str(freq_count[i][1].get_freq())))
            logging.write("mid dict 2: \n")
            logging.write("{:<.50}...\n".format(str(freq_count[i][2].get_freq())))
            logging.write("short dict : \n")
            logging.write("{:<.50}...\n".format(str(freq_count[i][3].get_freq())))
        else:
            logging.write("{:<.50}...\n".format(str(freq_count[i][0].get_freq())))
    thread_pool.close()
    thread_pool.join()

    logging.write("------gross frequency of all sequences and vector representation: ------\n")
    gross_frequency = kmer_statistics.gross_frequency(freq_count)  # returns a list
    vector_frequency = kmer_statistics.vector_calculator(gross_frequency)

    if args.subtraction:
        logging.write("long dict: \n")
        logging.write("{:<.50}...\n".format(str(gross_frequency[0].get_frequency_dict())))
        logging.write("mid dict 1: \n")
        logging.write("{:<.50}...\n".format(str(gross_frequency[1].get_frequency_dict())))
        logging.write("mid dict 2: \n")
        logging.write("{:<.50}...\n".format(str(gross_frequency[2].get_frequency_dict())))
        logging.write("short dict: \n")
        logging.write("{:<.50}...\n".format(str(gross_frequency[3].get_frequency_dict())))
    else:
        logging.write("{:<.50}...\n".format(str(gross_frequency[0].get_frequency_dict())))

    for item in vector_frequency[:min(4, len(vector_frequency))]:
        logging.write("key: {0}\tencoded key: {1}\tvalue: {2}\n".format(
            item[0], item[1], item[2]))

    logging.write("------saving result for {}: ------\n".format(name))
    io_manager.output_content(args, output_folder, vector_frequency, name)


if __name__ == '__main__':
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

    args = parser.parse_args()
    names, output_folder, logging, error = io_manager.output_error_logging_config(args)
    check_parameters(args)
    print("------Arguments:------\n{}".format(args))
    logging.write("------ARGUMENTS:------\n{}\n".format(args))
    kmer_tool = kmer_statistics.KmerStatistics(args, logging, error)
    for name in names:
        print("calculating {}".format(name))
        logging.write("CALCULATING {}\n".format(name.upper()))
        calculate(args, name, output_folder, kmer_tool)
        logging.write("\n")
    logging.close()

