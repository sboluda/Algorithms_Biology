import re

def indexing(string):
    """
    Write a function indexing(string) that returns a sorted list of tuples containing all possible
    offsets of the substrings in the string. Each tuple has two elements: the substring and the
    list of offsets for this substring in the string.
    
    >>> indexing("THEFASTCAT")
    [('A', [4, 8]), ('AS', [4]), ('AST', [4]), ('ASTC', [4]), ('ASTCA', [4]), ('ASTCAT', [4]), ('AT', [8]), ('C', [7]), ('CA', [7]), ('CAT', [7]), ('E', [2]), ('EF', [2]), ('EFA', [2]), ('EFAS', [2]), ('EFAST', [2]), ('EFASTC', [2]), ('EFASTCA', [2]), ('EFASTCAT', [2]), ('F', [3]), ('FA', [3]), ('FAS', [3]), ('FAST', [3]), ('FASTC', [3]), ('FASTCA', [3]), ('FASTCAT', [3]), ('H', [1]), ('HE', [1]), ('HEF', [1]), ('HEFA', [1]), ('HEFAS', [1]), ('HEFAST', [1]), ('HEFASTC', [1]), ('HEFASTCA', [1]), ('HEFASTCAT', [1]), ('S', [5]), ('ST', [5]), ('STC', [5]), ('STCA', [5]), ('STCAT', [5]), ('T', [0, 6, 9]), ('TC', [6]), ('TCA', [6]), ('TCAT', [6]), ('TH', [0]), ('THE', [0]), ('THEF', [0]), ('THEFA', [0]), ('THEFAS', [0]), ('THEFAST', [0]), ('THEFASTC', [0]), ('THEFASTCA', [0]), ('THEFASTCAT', [0])]
    
    """
    indexes = []
    substrings = {}
    for i in range(1, len(string) + 1):
        for j in range(len(string) - i + 1):
            # This way we obtain all the possible substrings
            substring = string[j : j+i] 
            if substring not in substrings: # We check if we already computeted the substring
                # Using regex and re.finditer we obtain a fast way to get all starting indexes (m.start()) for all substrings
                index_positions = [m.start() for m in re.finditer(substring, string)]

                # We combine both in a tuple
                substring_index = (substring, index_positions)

                # We add it to the main list
                indexes.append(substring_index)

                # We store the substring to avoid future repetitions
                substrings[substring] = "X"

    return sorted(indexes, key=lambda x: x[0])