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


def extract_dict(filename):
    """
    extract a dict from *.fa
    :param filename: path to the *.fa file
    :return: dict. key -> all the items after the key
    """
    seq_dict = {}
    temp_key = None
    for row in open(filename, 'r'):
        # in the file, each key is written in a line beginning with >
        # after that, each
        if row[0] == '>':
            temp_key = row[1:-1]
            seq_dict[temp_key] = []
        else:
            seq_dict[temp_key].append(row.strip())
    for key, value in seq_dict.items():
        seq_dict[key] = ''.join(value).upper()
    return seq_dict


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
