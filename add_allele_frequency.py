import requests
import re
import argparse

def get_allele_freq(vcf_line):
    '''this function retrieve exac database to get allele frequency'''
    items = vcf_line.strip().split('\t')
    web = 'http://exac.broadinstitute.org/variant/'
    web += '-'.join([items[0], items[1], items[3], items[4]])
    r = requests.get(web)
    if '<dt>Allele Frequency</dt>' in r.text:
        res = r.text.split('\n')
        for i in range(len(res)):
            if '<dt>Allele Frequency</dt>' in res[i]:
                allele_info = res[i+1]
                break
        allele_freq = re.search('(?<=\>).+?(?=\<)',allele_info).group(0)
    else:
        allele_freq = '-'
    return allele_freq


def split_multi_mutations(vcf_line):
    '''this function split vcf lines that have more than one mutations, it creates a line for each mutation'''
    items = vcf_line.strip().split('\t')
    n = len(items[4].strip().split(','))
    lines = [[''] * len(items) for i in range(n)]
    single_index = [0, 1, 2, 3, 5, 6] + list(range(8, 19))
    # first fill in items that the same values for different mutations
    for i in single_index:
        for line in lines:
            line[i] = items[i]
    # second fill in items that have multiple values correspond to different mutations
    for item in items[4].split(','):
        for line in lines:
            line[4] = item
    for item in  items[7].split(';'):
        index = item.index('=')
        pre = item[:index+1]
        post = item[index+1:].split(',')
        if len(post) == 1:
            for i in range(n):
                lines[i][7] += pre + post[0] + ';'
        else:
            for i in range(n):
                lines[i][7] += pre + post[i] + ';'
    res = ['\t'.join(l) for l in lines]
    return res



parser = argparse.ArgumentParser(description='add allele frequency information to annotated vcf file')
parser.add_argument('-i','--input',action='store',dest='input',help='input annotated vcf file')
parser.add_argument('-o','--out',action='store',dest='out',help='output vcf file')

args = parser.parse_args()

if __name__ == '__main__':
    anno_vcf_fn = args.input
    anno_vcf_freq = args.out
    # anno_vcf_fn = '/data/shangzhong/temp/anno_vcf_anno.vcf'
    # anno_vcf_freq = '/data/shangzhong/temp/anno_vcf_anno_freq.vcf'
    with open(anno_vcf_fn) as f, open(anno_vcf_freq, 'w') as out:
        for line in f:
            if line.startswith('#'):
                out.write(line.strip() + '\t' + 'allele_frequency' + '\n')
                continue
            lines = split_multi_mutations(line)
            allele = []
            for l in lines:
                freq = get_allele_freq(l)
                allele.append(freq)
            out.write(line.strip() + '\t' + ','.join(allele) + '\n')   