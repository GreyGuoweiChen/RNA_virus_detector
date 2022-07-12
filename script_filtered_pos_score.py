
import subprocess
from collections import Counter
global positive_cluster
positive_cluster = []

def read_thresholding():
    filename = "data/RNA_virus_4/thresholding_adaptive_T1_neg.txt"
    threshold = {}
    with open(filename, 'r') as f:
        for line in f:
            t = line.strip().split()
            threshold[int(t[0])] = float(t[1])
    return threshold

def read_rv_acc():
    filename = "data/RNA_virus_4/sequence_all_acc.txt"
    db_rc_acc = set()
    with open(filename, 'r') as f:
        for line in f:
            db_rc_acc.add(line.strip())
    return db_rc_acc


class contig:
    def __init__(self,fullname):
        self.fullname = fullname
        self.name = fullname.strip('>').split()[0]
        # self.length = int(fullname.strip().split()[-1].strip('len='))
        self.seq = None
        self.proteins = {}
        self.rnaviralness = 0

    def add_protein(self,prot_name,prot):
        self.proteins[prot_name] = prot

    def calculate_rnaviralness(self):
        t = 0
        for key, value in  self.proteins.items():
            if value.rnavaralness:
                t += 1
        if len(self.proteins):
            self.rnaviralness = t / len(self.proteins)
#####################################################################################
class protein:

    db_threshold = read_thresholding()
    db_rv_acc = read_rv_acc()

    def __init__(self,fullname):
        self.fullname = fullname
        self.contig_name = self._set_contig_name(fullname)
        self.seq=''

        self.best_hit = 0
        self.e_value = 1
        self.score = 0
        self.rnavaralness = 0
        self.diamond = None

    def _set_contig_name(self, fullname):
        t = fullname.strip('>').split()[0]
        t = t.rsplit('_',1)[0]
        return t

    def _set_protein_name(self,fullname):
        t=fullname.split('#')[0][1:-1]
        return t

    def filter2neg(self):
        # only keep the best hit; If it is not pos, then block other hit
        # self.best_hit = 0
        self.e_value = 1
        self.rnavaralness = 0
        self.score = 0

    def parse_match_search(self, result):
        t = result.split()
        hit_clustet_index, e_value, score  \
            = int(t[2].strip("cluster_")), float(t[4]), float(t[5])

        if score > self.score:
            self.score, self.e_value, self.best_hit \
                = score, e_value, hit_clustet_index

    def rnavaralness_for_search(self):
        if self.best_hit and self.score >= protein.db_threshold[self.best_hit] and self.e_value< 1e-3:
            self.rnavaralness = 1
            positive_cluster.append(self.best_hit)
            print(self.fullname.split('#')[0][1:-1], '\t', self.best_hit, '\t', self.score, '>',
                  protein.db_threshold[self.best_hit], '\te_value:', self.e_value)


    def parse_match_diamond_blastx(self, result):
        t = result.split()
        hit_acc = t[1]
        if hit_acc in protein.db_rv_acc:
            self.rnavaralness, self.diamond_acc = \
                1, hit_acc
            print(self.fullname.split('#')[0][1:-1], '\t', hit_acc,
                  '\te_value:',t[-2],'\tscore:',t[-1])

#####################################################################################
def predict(filepath,
            filename1, filename2, filename3,
            filename4, filename5, prot_pred=False):
    """

    :param filename1: proteins file for contigs need to be predicted
    :param filename2: hmmscan result for filename1 (tbl_report)
    :param filename3: contigs need to be predicted
    :param filename4: positive contigs predicted by our strategy
    :param filename5: diamond result for filename1 (diamond)

    :return: no output, but write down the predicted fasta
    """
    def protein_predict(filepath, filename1, filename2):
        subprocess.call("source activate PHD1;"
                        "cd %s;"
                        "prodigal -p meta -a %s -i %s;"
                        "conda deactivate "
                        % (filepath, filename2, filename1), shell=True)

    def parse_protein(filename):
        proteins={}
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('>'):
                    prot_name = line.strip('>').split()[0]
                    proteins[prot_name]=protein(line)
        return proteins

    def parse_hmmsearch(filename,proteins):
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    t = line.strip().split()
                    prot_name = t[0]
                    if prot_name in proteins:
                        proteins[prot_name].parse_match_search(line)

        for prot_name, protein in proteins.items():
            protein.rnavaralness_for_search()

        return proteins

    def parse_diamond_blastx(filename, proteins):
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    t = line.strip().split()
                    prot_name = t[0]
                    if prot_name in proteins:
                        proteins[prot_name].parse_match_diamond_blastx(line)

        return proteins

    def parse_contig(filename,proteins):
        contigs={}
        with open(filename,'r') as f:
            for line in f:
                if line.startswith('>'):
                    contig_name = line.strip('>').split()[0]
                    contigs[contig_name] = contig(line)

        for prot_name, prot in proteins.items():
            if prot.contig_name in contigs:
                contigs[prot.contig_name].add_protein(prot_name,prot)

        for c_name, c in contigs.items():
            c.calculate_rnaviralness()

        return contigs

    def output_rnaviralness(filename, contigs):
        i = 0
        print('contig_name\tRNAviralness\tprotein_num')
        positive_contigs = {}
        # 如果contig呈阳性，则从contigs中记录下该contig(不包含seq)
        for c_name,c in contigs.items():
            # 此处可以对rnaviralness进行阈值的划分，
            # 考虑到该数据组装出来的片段较小，暂时不划分。
            if c.rnaviralness >= 0.0625:
                print(c_name,'\t',c.rnaviralness,'\t',len(c.proteins))
                positive_contigs[c_name] = c
            if c.proteins.__len__()>0:
                i+=1

        print("total num of contigs:",contigs.__len__())
        print("num of contigs containing gene(s):",i)
        print("num of contigs lacking gene(s)",contigs.__len__()-i)
        print("num of positive_contigs:",positive_contigs.__len__())
        print("Counter of cluster hit",Counter(positive_cluster))

        # 对positive contig 补充记录seq
        with open(filename,'r') as f:

            significant =False
            seq = ''
            contig_name = None

            for line in f:
                if line.startswith('>'):
                    #在positive_contigs中，存下对应的seq
                    if significant:
                        positive_contigs[contig_name].seq= seq

                    significant = False
                    contig_name = line.strip('>').split()[0]
                    if contig_name in positive_contigs:
                        significant = True
                        seq=''
                elif significant:
                    seq += line
            #do not ignore the last contig
            if significant:
                positive_contigs[contig_name].seq = seq

        return positive_contigs

    def write_positive_file(filename,positive_contigs):
        with open(filename,'w') as f:
            for c_name,c in positive_contigs.items():
                f.write(c.fullname)
                f.write(c.seq)

#####################################################################################

    # add prodigal module
    # gene prediction
    if prot_pred:
        protein_predict(filepath,filename3,filename1)

    # store all the predicted proteins into a ditc {prot_name: proteins} without the seq,
    proteins = parse_protein(filepath + filename1)
    print("num of proteins",proteins.__len__())
    # add hmmsearch module for prot

    # parse the match result for all proteins
    # proteins = parse_hmmsearch(filepath + filename2, proteins)
    # print("parsing of protein HMM-match result finished")

    proteins = parse_diamond_blastx(filepath + filename5, proteins)
    print("parsing of protein DIAMOND-match result finished")

    #统计hmmscan结果，并统计到contig里面去。
    contigs = parse_contig(filepath + filename3, proteins)

    #对每条contig进行评分，返回包含seq 的 positive contigs
    positive_contigs = output_rnaviralness(filepath + filename3,contigs)

    #所有postitive contig写文件
    write_positive_file(filepath + filename4,positive_contigs)

#####################################################################################
if __name__ == "__main__":

    # filepath = "data/other/SPAV1/"
    # filename1 = "contigs_500.faa"
    # filename2 = "tbl4/tbl_search_contig_500"
    # filename3 = "contigs_500.fasta"
    # filename4 = "tbl4/pos_contig.fasta"
    # predict(filepath,
    #         filename1,
    #         filename2,
    #         filename3,
    #         filename4,
    #         prot_pred = False)

    # filepath = "data/other/Salmon_meta/"
    # filename1 = "contigs_1000.faa"
    # filename2 = "tbl4/tbl_report_non_salmon_nr"
    # filename3 = "tbl4/non_salmon_1000.contigs.fa"
    # filename4 = "pos_contig.fasta"
    # predict(filepath,
    #         filename1,
    #         filename2,
    #         filename3,
    #         filename4,
    #         prot_pred = False)

    # filepath = "data/other/PRJNA605028/"
    # filename1 = "contigsPRJNA605028_RNA.over500.faa"
    # filename2 = "tbl4/tbl_search_PRJNA_500"
    # filename3 = "contigsPRJNA605028_RNA.over500.fasta"
    # filename4 = "tbl4/PRJNA_500.rvd.fasta"
    # filename5 = "tbl4/pos_contig_500.diamond.rvd.out"
    # F = False
    # predict(filepath,
    #         filename1,
    #         filename2,
    #         filename3,
    #         filename4,
    #         filename5,
    #         prot_pred = False)

    # filepath = "data/other/marine_virus/"
    # filename1 = "contigs_1000.faa"
    # filename2 = "tbl4/tbl_search_marine_1000"
    # filename3 = "tbl4/contigsMarine.over1000.fasta"
    # filename4 = "pos_contig.fasta"
    # predict(filepath,
    #         filename1,
    #         filename2,
    #         filename3,
    #         filename4,
    #         prot_pred = False)

    ERR_num = 4
    filepath = "data/other/ERR%d/"%(ERR_num)
    filename1 = "ERR%d_500.contig.faa"%(ERR_num)
    filename2 = "tbl4/tbl_search_ERR%d_500"%(ERR_num)
    filename3 = "ERR%d_500.contig.fasta"%(ERR_num)
    filename4 = "tbl4/ERR%d_500.rvd.fasta"%(ERR_num)
    filename5 = "tbl4/ERR%d_500_hit.diamond.rvd.out"%(ERR_num)
    predict(filepath,
            filename1,
            filename2,
            filename3,
            filename4,
            filename5,
            prot_pred = False)

    # filepath = "data/other/aay5981_data_s2/ssRNAphage_sequences/allNew/"
    # filename1 = "allNew_all.faa"
    # filename2 = "tbl4/tbl_search_allNew_500"
    # filename3 = "allNew_500.fasta"
    # filename4 = "tbl4/allNew_500.rvd.fasta"
    # filename5 = "tbl4/pos_contig_500.diamond.rvd.out"
    # predict(filepath,
    #         filename1,
    #         filename2,
    #         filename3,
    #         filename4,
    #         filename5,
    #         prot_pred = False)

    # filepath = "data/other/gut_virome/"
    # filename1 = "final_500.contigs.faa"
    # filename2 = "tbl/tbl_search_500"
    # filename3 = "final.contigs.fasta"
    # filename4 = "tbl/pos_contig_500.rvd.fasta"
    # predict(filepath,
    #         filename1,
    #         filename2,
    #         filename3,
    #         filename4,
    #         prot_pred = False)

    # filepath = "data/other/human_CYRN/"
    # filename1 = "spades/spades_500.faa"
    # filename2 = "tbl4/tbl_report_human_spades"
    # filename3 = "spades/spades_500.fasta"
    # filename4 = "spades/pos_contig_spades_500.fasta"
    # predict(filepath,
    #         filename1,
    #         filename2,
    #         filename3,
    #         filename4,
    #         prot_pred = False)