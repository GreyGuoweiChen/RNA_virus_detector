import numpy as np
import matplotlib.pyplot as plt
import scipy
import seaborn as sns


def read_distribution(cluster):

    def read_pos(filename, i):

        score_list = []
        record_stat = 0
        with open(filename + 'report_%d.out' % i, 'r') as f:
            # #statistic e-value in one clusters
            for line in f:
                if line.startswith('    -------'):
                    record_stat = 1
                    continue
                if record_stat:
                    if line == '\n' or line.startswith('  ------'):
                        break
                    score_list.append(float(line.split()[1]))
        return score_list

    def read_neg(filename, i):

        score_list = []
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                t = line.split()
                cluster, score = int(t[2].split('_')[1]), float(t[5])
                if cluster == i:
                    score_list.append(score)
        return score_list

    filename = "data/RNA_virus_4/clustering/selected_clusters/mapback_report/"
    pos = read_pos(filename, cluster)

    filename = "data/RNA_virus_4/clustering/selected_clusters/threshold/tbl_report_SP_cell_T1"
    neg = read_neg(filename, cluster)

    filename = "data/RNA_virus_4/clustering/selected_clusters/threshold/tbl_report_SP_DNAV_T1"
    neg += read_neg(filename, cluster)

    filename = "data/RNA_virus_4/clustering/selected_clusters/threshold/tbl_rep_bac"
    neg += read_neg(filename, cluster)

    print("Pos_distribution of cluster_%d"%cluster,pos)
    print("Neg_distribution of cluster_%d"%cluster,neg)
    return  neg, pos

def drew_distribution_1(cluster,neg, pos):
    plt.hist(neg, bins=1000, range=(0,5000), # normed=True,
                weights=None, cumulative=False, bottom=None,
                histtype=u'bar', align=u'left', orientation=u'vertical',
                rwidth=0.8, log=False, color='r', label=None, stacked=False)
                #hold=None)
    plt.hist(pos, bins=50, range=(0,5000), # normed=True,
                weights=None, cumulative=False, bottom=None,
                histtype=u'bar', align=u'right', orientation=u'vertical',
                rwidth=0.8, log=False, color='b', label=None, stacked=False)
                #hold=None)

    plt.show()





def drew_distribution_2(cluster,neg, pos):
    plt.style.use('ggplot')
    #处理中文乱码
    #plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
    #plt.rcParams['axes.unicode_minus']=False
    # seaborn module drew the histogram and curve;
    # 绘制男女乘客年龄的直方图
    ax1 = sns.histplot(pos, bins = 20, kde = False, color="steelblue", label = 'Pos',stat="probability")
    # ax1.grid(True)
    # ax1.set_xlabel('Score')
    # 绘制女性年龄的直方图
    # ax2 = ax1.twinx()
    # ax2.set_ylabel('Neg_hist')
    # ax2.set_xlabel('Score')   , ax=ax2
    sns.histplot(neg, bins = 4, kde = False, color="purple", label = 'Neg',stat="probability")
    # plt.title('Bit score histogram of C_%d'%cluster)
    plt.xlabel('Bit score')
    plt.ylabel('probability')
    # plt.ylabel('Frequency')
    # 显示图例
    plt.legend()
    # 显示图形
    #plt.savefig('bar-6-1.png')
    plt.show()

    sns.histplot(neg, bins=4, kde=False, color="purple", label='Neg', stat="probability")
    plt.xlabel('Bit score')
    plt.ylabel('probability')
    #plt.savefig('bar-6-2.png')
    plt.show()

    # kernel density curve graph

    ax1 = sns.distplot(pos, hist = False, kde_kws = {'color':'red', 'linestyle':'-'},
                 norm_hist = True, label = 'Pos')
    ax1.set_ylabel('density')
    ax1.set_xlabel('Bit score')
    plt.yticks([])
    #plt.axis('off')
    # 绘制女性年龄的核密度图
    ax2 = ax1.twinx()
    sns.distplot(neg, ax=ax2, hist = False, kde_kws = {'color':'black', 'linestyle':'--'},
                 norm_hist = True, label = 'Neg')
    ax2.set_ylabel('')
    plt.yticks([])
    plt.xticks([])
    #ax2.set_xlabel('Bit score')
    #plt.title('Bit score distribution')
    #plt.legend()
    # plt.axis('off')
    ax1.set_facecolor('w')
    plt.show()



if __name__ == "__main__":

    cluster = 133
    neg, pos = read_distribution(cluster)

    drew_distribution_2(cluster,neg, pos)


