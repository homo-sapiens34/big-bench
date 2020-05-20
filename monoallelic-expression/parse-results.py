#!/usr/bin/env python3

from os import listdir
from os.path import isfile, join
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt


black_list = set(
    [0, 4, 5, 6, 9, 12, 15, 42, 45, 49]
)  # smth went wrong with star mapping


def get_files(mypath):
    return [f for f in listdir(mypath) if isfile(join(mypath, f))]


def parse_kallisto_file(fname):
    df = pd.read_csv(fname, sep="\t")
    gene2counts = {}
    for i in range(len(df)):
        s = df.iloc[i]
        gene, counts = s["target_id"], s["est_counts"]
        gene_id = gene[4:]
        if not gene_id in gene2counts:
            gene2counts[gene_id] = {}
        gene2counts[gene_id][gene[:3]] = counts
    gene2ans = {}
    for gene in gene2counts:
        su = gene2counts[gene]["pat"] + gene2counts[gene]["mat"]
        if su == 0:
            gene2ans[gene] = None
        else:
            gene2ans[gene] = gene2counts[gene]["pat"] / su
    return gene2ans


def calc_metrics(dict_ans, dict_true):
    nans = 0
    n = 0
    su = 0
    v1, v2 = [], []
    for gene in dict_true:
        ans, true = dict_ans[gene], dict_true[gene]
        if ans == None:
            nans += 1
            continue
        su += (ans - true) ** 2
        n += 1
        v1.append(ans)
        v2.append(true)
    r = pearsonr(v1, v2)[0]
    return {"MSE": su / n, "r": r, "nans_fraction": nans / len(dict_true)}


def parse_sim_info_file(fname):
    df = pd.read_csv(fname, sep="\t", header=None)
    gene2ans = {}
    for i in range(len(df)):
        s = df.iloc[i]
        gene, ans = s[0], s[1]
        gene2ans[gene] = ans
    return gene2ans


def parse_rsem_file(fname):
    df = pd.read_csv(fname, sep="\t")
    gene2counts = {}
    for i in range(len(df)):
        s = df.iloc[i]
        gene, counts = s["gene_id"], s["expected_count"]
        gene_id = gene[4:]
        if not gene_id in gene2counts:
            gene2counts[gene_id] = {}
        gene2counts[gene_id][gene[:3]] = counts
    gene2ans = {}
    for gene in gene2counts:
        su = gene2counts[gene]["pat"] + gene2counts[gene]["mat"]
        if su == 0:
            gene2ans[gene] = None
        else:
            gene2ans[gene] = gene2counts[gene]["pat"] / su
    return gene2ans


def parse_dir(dir_name):
    files = get_files(dir_name)
    num2metrics = {}
    for fi in files:
        num = int(fi.split(".")[0].split("_")[1])
        if num in black_list:
            continue
        ans = (
            parse_kallisto_file(join(dir_name, fi))
            if "kallisto" in dir_name
            else parse_rsem_file(join(dir_name, fi))
        )
        true = parse_sim_info_file("sim_data/" + str(num) + "_gene_info.tsv")
        res = calc_metrics(ans, true)
        num2metrics[num] = res
    return num2metrics


if __name__ == "__main__":
    num2kal_metrics = parse_dir("kallisto_res")
    num2rsem_metrics = parse_dir("rsem_res")

    plt.figure(figsize=(5, 5))

    kal, rsem = [], []

    for n in num2kal_metrics:
        kal.append(num2kal_metrics[n]["r"])
        rsem.append(num2rsem_metrics[n]["r"])

    diff = [kal[i] - rsem[i] for i in range(len(kal))]
    plt.plot([0.924, 0.932], [0.924, 0.932], "r")
    plt.plot(kal, rsem, "bo")
    plt.xlabel("Kallisto correlation coefficient")
    plt.ylabel("RSEM correlation coefficient")
    plt.show()

    plt.pie([len([x for x in diff if x > 0]), len([x for x in diff if x < 0])])
    plt.legend(["Kallisto is better", "RSEM is better"], loc="lower right")
    plt.show()

