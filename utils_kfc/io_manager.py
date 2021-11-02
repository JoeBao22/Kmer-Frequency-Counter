import os
import time


def create_dir(full_path):
    """
    create directory given path
    :param full_path: full path to the directory
    :return: None
    """
    if not os.path.exists(full_path):
        os.mkdir(full_path)


def extract_fasta(filename):
    """
    extract a list from *.fa
    :param filename: path to the *.fa file
    :return: all the items after the key
    """
    seq_list = []
    temp_seq = None
    for row in open(filename, 'r'):
        # in the file, each key is written in a line beginning with >
        # after that, each
        if row[0] == '>':
            if temp_seq is not None:
                seq_list.append(temp_seq)
            temp_seq = ""
        else:
            temp_seq += row.strip()
    if temp_seq is not None:
        seq_list.append(temp_seq)
    return seq_list

def extract_fastq(upstream_filename, downstream_filename):
    # use 1 to represent upstream_file, and 2 to represent downstream file
    seq_list = []
    for filename in [upstream_filename, downstream_filename]:
        temp_seq = ""
        for row_index, row in enumerate(open(filename, 'r')):
            if row_index % 4 == 1:
                if len(temp_seq) == 0 :
                    temp_seq = row.strip()
                else:
                    temp_seq = temp_seq + "N" + row.strip()
        seq_list.append(temp_seq)
    return seq_list   

# def extract_fastq(upstream_filename, downstream_filename):
#     # use 1 to represent upstream_file, and 2 to represent downstream file
#     seq_list = []
#     is_next_line_seq = False
#     for row_index, row in enumerate(open(upstream_filename, 'r')):
#         if row_index % 4 == 1:
#             seq_list.append((1, row_index+1, row.strip()))
#     for row_index, row in enumerate(open(downstream_filename, 'r')):
#         if row_index % 4 == 1:
#             seq_list.append((2, row_index+1, row.strip()))
#     return seq_list     
    
    
def fetch_names(name_path):
    if not os.path.exists(name_path):
        print("Name path {} doesn't exist.")
        exit(1)
    name_fd = open(name_path, "r")
    names = [name.strip() for name in name_fd.readlines() if name.strip() ]
    return names


def get_output_folder(args):
    if args.space:
        k = args.length + (args.length - 1) * args.interval
        suffix = "_space_N{}_{}".format(args.length, args.interval)
        if args.direction:
            suffix += "_{}_{}".format(args.position, args.direction)
        else:
            suffix += "_-0"
    else:
        k = args.length
        suffix = ""
    folder_path = "K{}".format(str(k))
    if args.combine:
        folder_path += "_comb"
    if args.frequency:
        folder_path += "_freq"
    if args.subtraction:
        folder_path += "_sub"
    if args.new_encoder:
        folder_path += "_new"
    else:
        folder_path += "_old"
    if args.use_fastq:
        folder_path += "_fastq"
    folder_path += suffix
    create_dir(folder_path)
    return folder_path


def output_error_logging_config(args):
    temp_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    config_filename = "config_{}.txt".format(temp_time)
    logging_filename = "logging_{}.txt".format(temp_time)
    error_filename = "error_{}.txt".format(temp_time)
    logging_fd = open(os.path.join("\\".join(args.name_file.split("\\")[: -1]), logging_filename), "w")
    error_fd = open(os.path.join("\\".join(args.name_file.split("\\")[: -1]), error_filename), "w")
    if not os.path.exists(args.name_file):
        print("Name path {} doesn't exist.".format(args.name_file))
        error_fd.write("Name path {} doesn't exist\n.".format(args.name_file))
        exit(1)
    name_fd = open(args.name_file, "r")
    names = [name.strip() for name in name_fd.readlines() if name.strip()]
    config = ""
    for (k, v) in sorted(vars(args).items()):
        config += "{key}: {value}\n".format(key=k, value=v)
    output_folder = get_output_folder(args)
    open(os.path.join(output_folder, config_filename), "w").write(config)
    return names, output_folder, logging_fd, error_fd


def output_content(args, output_folder, frequency, species_name):
    # Output as a format of cvtree
    count_none_zero = 0
    inner_product = sum([item[2] ** 2 for item in frequency])  # item[2]: value
    content = ""
    for item in frequency:
        if item[2] != 0.:  
            count_none_zero += 1
            content += "{} {}\n".format(item[1], item[2])  # item[1]: key, item[2]: value
    title = "{}\n{}\n{}\n".format(args.n, inner_product, count_none_zero)
    create_dir(output_folder)
    filename_species = "{}.cv.txt".format(species_name)
    open(os.path.join(output_folder, filename_species), 'w').write(title + content)


def output_content_fastq(args, output_folder, record, species_name):
    create_dir(output_folder)
    content = ""
    for item in record:
        content += "{}\t{}\t{}\n".format(item[1], item[2], item[3])
    filename_species = "{}.cv.txt".format(species_name)
    open(os.path.join(output_folder, filename_species), 'w').write(content)