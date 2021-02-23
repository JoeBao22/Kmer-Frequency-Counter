
def item_in_dict(d):
    return sum(d.values())


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
            from .encoder_decoder_old import encoder, decoder
        else:
            from .encoder_decoder import encoder, decoder
        self.encoder = encoder
        self.decoder = decoder
        self.loc = list(range(0, self.k, self.d + 1))  # the chosen positions when spaced
        self.loc[self.position] += self.direction






