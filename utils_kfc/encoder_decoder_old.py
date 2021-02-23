
encoder_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
decoder_dict = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}


def encoder(kmer_str):
    """
    encoded = sum(encoder_dict[i] * (4 ** i)), where i is the index for each character in the reversed string
    :param kmer_str: a string that is composed of A, C, G, T
    :return: the numerical representation for the string
    """
    global encoder_dict
    numerical_representation = 0
    for index, character in enumerate(kmer_str[::-1]):  # reverse string, the last character has the smallest index
        numerical_representation += encoder_dict[character] * 4 ** index
    return numerical_representation


def decoder(kmer_int, length):
        """
        find the original string from the encoded one(numerical representation)
        :param kmer_int: the numerical representation for the string
        :param length: the length of the decoded str
        :return: a string that is composed of A, C, G, T
        """
        global decoder_dict
        string_representation = ""
        while length > 0:
            temp_remainder = kmer_int % 4  # consider the last character that has not been considered
            kmer_int = kmer_int // 4
            string_representation = decoder_dict[temp_remainder] + string_representation
            length -= 1
        return string_representation
    