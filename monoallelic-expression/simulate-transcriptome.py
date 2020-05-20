#!/usr/bin/env python3

import random
from scipy.stats import poisson, gamma, nbinom, beta
from tqdm import tqdm as tqdm

random.seed(34)

NS = "ATGC"
P_MUT = 0.01
READ_LEN = 100
N_GENES = 10000


def gen_seq(l):
    return "".join([random.choice(NS) for _ in range(l)])


def gen_mutated_seq(seq):
    mseq = ""
    for i in range(len(seq)):
        n, new_n = seq[i], seq[i]
        if random.random() <= P_MUT:
            while new_n == n:
                new_n = random.choice(NS)
        else:
            new_n = n
        mseq += new_n
    return mseq


def get_gene_len():
    return nbinom.rvs(4, 0.002)


def get_lambda():
    ## expression per nucleotide
    return nbinom.rvs(1, 0.1)


def get_ai():
    if random.random() <= 0.2:
        return beta.rvs(1 / 5, 1 / 5)
    else:
        return beta.rvs(30, 30)


def gen_files(dir_name, n):
    with open(dir_name + "/" + str(n) + "_mat_transcriptome.fa", "w") as mw, open(
        dir_name + "/" + str(n) + "_pat_transcriptome.fa", "w"
    ) as pw, open(dir_name + "/" + str(n) + "_reads.fa", "w") as rw, open(
        dir_name + "/" + str(n) + "_gene_info.tsv", "w"
    ) as iw:
        for i in range(N_GENES):
            l = get_gene_len()
            if l < READ_LEN * 3:
                continue
            p_seq = gen_seq(l)
            m_seq = gen_mutated_seq(p_seq)
            if p_seq == m_seq:
                continue

            pw.write(">pat_gene_" + str(i) + "\n")
            pw.write(p_seq + "\n")
            mw.write(">mat_gene_" + str(i) + "\n")
            mw.write(m_seq + "\n")

            lam = get_lambda()
            ai = get_ai()

            iw.write("\t".join([str(x) for x in ["gene_" + str(i), ai, lam]]) + "\n")

            read_count = lam * l // READ_LEN
            for i in range(read_count):
                pos = random.randint(0, l - READ_LEN - 1)
                sourse = p_seq if random.random() < ai else m_seq
                read = sourse[pos : pos + READ_LEN]
                rw.write(">read_" + str(i) + "\n")
                rw.write(read + "\n")


if __name__ == "__main__":
    for i in tqdm(range(50)):
        gen_files("sim_data", i)
