import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import seaborn as sns

def read_distribution(N):

    def read_pos(filename, i):

        score = 5000
        record_stat = 0
        score_list = []
        with open(filename + 'report_%d.out' % i, 'r') as f:
            # #statistic e-value in one clusters
            for line in f:
                if line.startswith('    -------'):
                    record_stat = 1
                    continue
                if record_stat:
                    if line == '\n' or line.startswith('  ------'):
                        break
                    #score_list.append(float(line.split()[1]))
                    score_list.append(float(line.split()[1]))


        index = -1
        #index = score_list.__len__()//10-1
        #index = -2
        score = score_list[index]
        return score

    def read_neg(filename, neg_bound):

        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                t = line.split()
                cluster, score = int(t[2].split('_')[1]), float(t[5])
                neg_bound[cluster] = max(neg_bound[cluster], score)
        return neg_bound

    pos_bound = defaultdict(int)
    for cluster in range(1, N+1):
        filename = "data/RNA_virus_4/clustering/selected_clusters/mapback_report/"
        pos_bound[cluster] = read_pos(filename, cluster)

    neg_bound = defaultdict(int)
    filename = "data/RNA_virus_4/clustering/selected_clusters/threshold/tbl_report_SP_cell_T1"
    neg_bound = read_neg(filename, neg_bound)

    filename = "data/RNA_virus_4/clustering/selected_clusters/threshold/tbl_report_SP_DNAV_T1"
    neg_bound = read_neg(filename, neg_bound)

    filename = "data/RNA_virus_4/clustering/selected_clusters/threshold/tbl_rep_bac"
    neg_bound = read_neg(filename, neg_bound)


    return  pos_bound, neg_bound

def drew_distribution(diff):
    bin = []
    for i in range(-1000,2600,100):
        bin.append(i)
    plt.style.use('ggplot')
    # ax1.grid(True)
    ax1 = sns.histplot(diff, bins=bin, kde=False, line_kws={'color': "red"}, label='diff',stat="probability") # 'steelblue'
    ax1.set_xlim(-1000, 2500)
    # ax1.grid(True)

    # ax1.set_xlabel('Score')
    # plt.title('Bit score difference histogram')
    plt.xlabel('bit score difference')
    plt.ylabel('probability')
    # plt.ylabel('Frequency')
    # plt.legend()
    # plt.savefig('bar-5.png')
    plt.show()


if __name__ == "__main__":

    N = 1384
    pos_bound, neg_bound =  read_distribution(N)
    diff = []
    min, max = float("inf"), float("-inf")
    min_clu , max_clu =0, 0
    num_range = 0
    for i in range(1,N+1):
        diff.append(pos_bound[i]-neg_bound[i])
        if pos_bound[i]-neg_bound[i] <= 0:
            #print(i)
            if pos_bound[i]-neg_bound[i] < min :
                min =  pos_bound[i]-neg_bound[i]
                min_clu = i
        if pos_bound[i]-neg_bound[i] > max :
            max = pos_bound[i]-neg_bound[i]
            max_clu = i
        if pos_bound[i]-neg_bound[i]>=0 and pos_bound[i]-neg_bound[i]<= 200:
            num_range += 1
        if pos_bound[i]-neg_bound[i]>=0 and pos_bound[i]-neg_bound[i]<= 100:
            print("!!!!!",i)
    drew_distribution(diff)
    print(max_clu,"MAX:", max)
    print(min_clu,"MIN:", min)
    print(num_range, num_range/N)