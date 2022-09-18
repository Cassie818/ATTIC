import os
import re
import numpy as np

datapath = os.path.abspath('.')

class Sequence(object):
    def __init__(self, file):
        self.file = file
        self.min_length = 51

    def read_fasta(self):
        """
        read fasta sequence
        :param file:
        :return:
        """
        msg = ''
        if not os.path.exists(self.file):
            msg = 'Error: file %s does not exist.' % self.file
            return [], None, msg
        with open(self.file) as f:
            records = f.read()

        if re.search('>', records) == None:
            print("The input file seems not in fasta format!")

        records = records.split('>')[1:]
        myFasta = []
        for fasta in records:
            array = fasta.split('\n')
            name, sequence = array[0].split()[0], re.sub('[^AUCG-]', '-', ''.join(array[1:]).upper())
            myFasta.append([name, sequence])

        return myFasta

    def output_fasta(self, fastafile,species):
        single_seq = ''
        sample = []
        for i in fastafile:
            length = len(i[1])
            if species == "M.musculus":
                for j in range(20, length-20):
                    if i[1][j] == 'A':
                        single_seq += '>' + i[0] + '_' + str(j) + '\n' + i[1][j-20:j+21] + '\n'
                        sample.append((i[0] + '_' + str(j)))
                with open(datapath +'/data/M/M.txt', 'w') as f:
                    f.write(single_seq)
            else:
                for j in range(25, length-25):
                    if i[1][j] == 'A':
                        single_seq += '>' + i[0] + '_' + str(j) + '\n' + i[1][j-25:j+26] + '\n'
                        sample.append((i[0] + '_' + str(j)))
                if species == "D.melanogaster":
                    with open(datapath +'/data/D/D.txt', 'w') as f:
                        f.write(single_seq)
                else:
                    with open(datapath + '/data/H/H.txt', 'w') as f:
                        f.write(single_seq)
        return sample



if __name__ == '__main__':
    sequence = Sequence('/Users/cassie/PycharmProjects/Oner/data/D/D.fasta')
    myfasta = sequence.read_fasta()
    samplename = sequence.output_fasta(myfasta,"D.melanogaster")
    print(samplename)




