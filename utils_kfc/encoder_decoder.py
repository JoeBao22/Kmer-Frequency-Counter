class EncoderDecoder:
    def __init__(self):
        self.encoder_dict = {'A': 1, 'C': 2, 'G': 3, 'T': 4}  # zero: Nothing left
        self.decoder_dict = {1: 'A', 2: 'C', 3: 'G', 4: 'T'}

    def encoder(self, kmer_str):
        """
        encoded = sum(encoder_dict[i] * (5 ** i)), where i is the index for each character in the reversed string
        :param kmer_str: a string that is composed of A, C, G, T
        :return: the numerical representation for the string
        """
        numerical_representation = 0
        for index, character in enumerate(kmer_str[::-1]):  # reverse string, the last character has the smallest index
            numerical_representation += self.encoder_dict[character] * 5 ** index
        return numerical_representation

    def decoder(self, kmer_int):
        """
        find the original string from the encoded one(numerical representation)
        :param kmer_int: the numerical representation for the string
        :return: a string that is composed of A, C, G, T
        """
        string_representation = ""
        while kmer_int != 0:
            temp_remainder = kmer_int % 5  # consider the last character that has not been considered
            kmer_int = kmer_int // 5
            string_representation = self.decoder_dict[temp_remainder] + string_representation
        return string_representation


if __name__ == "__main__":
    import random

    def generate_string(length):
        choices = ["A", "C", "T", "G"]
        result_str = "".join(random.choices(choices, k=length))
        return result_str

    e_d = EncoderDecoder()
    repeat_time = 100
    n = 10
    for _ in range(repeat_time):
        s = generate_string(n)
        encoded = e_d.encoder(s)
        recovered = e_d.decoder(encoded)
        assert recovered == s, "error when string = {}".format(s)
    print("All tests passed")
