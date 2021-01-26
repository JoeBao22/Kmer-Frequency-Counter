class EncoderDecoderOld:
    def __init__(self):
        self.encoder_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        self.decoder_dict = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

    def encoder(self, kmer_str):
        """
        encoded = sum(encoder_dict[i] * (4 ** i)), where i is the index for each character in the reversed string
        :param kmer_str: a string that is composed of A, C, G, T
        :return: the numerical representation for the string
        """
        numerical_representation = 0
        for index, character in enumerate(kmer_str[::-1]):  # reverse string, the last character has the smallest index
            numerical_representation += self.encoder_dict[character] * 4 ** index
        return numerical_representation

    def decoder(self, kmer_int, length):
        """
        find the original string from the encoded one(numerical representation)
        :param kmer_int: the numerical representation for the string
        :param length: the length of the decoded str
        :return: a string that is composed of A, C, G, T
        """
        string_representation = ""
        while length > 0:
            temp_remainder = kmer_int % 4  # consider the last character that has not been considered
            kmer_int = kmer_int // 4
            string_representation = self.decoder_dict[temp_remainder] + string_representation
            length -= 1
        return string_representation


if __name__ == "__main__":
    import random


    def str2int(Kmer):
        nucleotide = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        Kmer_list = list(Kmer)
        location = len(Kmer) - 1
        Kmer_value = 0
        for i in Kmer_list:
            Kmer_value += nucleotide[i] * 4 ** location
            location -= 1
        return Kmer_value


    def generate_string(length):
        choices = ["A", "C", "T", "G"]
        result_str = "".join(random.choices(choices, k=length))
        return result_str

    e_d = EncoderDecoderOld()
    repeat_time = 100
    n = 10

    for _ in range(repeat_time):
        s = generate_string(n)
        # print(s)
        encoded = e_d.encoder(s)
        # print("my:", encoded)
        # print("original:", str2int(s))
        assert encoded == str2int(s), "error encoding when stirng = {}".format(s)
        # print(encoded)
        recovered = e_d.decoder(encoded, n)
        assert recovered == s, "error when string = {}".format(s)
    print("All tests passed")
