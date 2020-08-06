import collections



def findOric(kstr,kmer):
    count = 0
    words = collections.defaultdict(list)
    mwords = []
    for i in range(0, len(kstr) - kmer + 1):
        word = kstr[i:i + kmer]
        if len(word) == kmer:
            if word in words:
                words[word] = words[word] + 1
                if words[word] > count:
                    count = words[word]
                    mwords = []
                    mwords.append(word)
                elif words[word] == count:
                    mwords.append(word)
            else:
                words[word] = 1

    complete=""
    for w in mwords:
        complete=complete+w


    oric_dict={}
    oric_dict = {"OriginOfReplication": str(complete),"totalOccurences":count,"LengthOfKmer":kmer}
    return oric_dict


