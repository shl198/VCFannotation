# build interval trees
from intervaltree import Interval, IntervalTree
from Bio import SeqIO
from Bio.Seq import Seq
import re
import argparse

def build_interval_tree(pred_fn):
    '''this functions build interval trees using genepred file
    the intervals are start and end positions of all the mRNAs.
    With this tree we can quickly decide if a mutation locates in 
    any rna
    * pred_fn: genepred file name
    return two objects, trees: all rna intervals.  '''
    trees = {}
    rna_info_dic = {}
    with open(pred_fn) as f:
        for line in f:
            items = line.strip().split('\t')
            rna, chrom, start, end = items[0], items[1], int(items[3]), int(items[4])
            trees.setdefault(chrom,IntervalTree())[start:end] = rna
            rna_info_dic[rna] = items
    return trees, rna_info_dic


def build_rna_tree(rna_info):
    '''this function builds a tree for a single rna. It has intervals of UTR5, CDS, UTR3 and intron intervals.
    This can be used to check if a mutation overlap with any of the regions.
    * rna_info: value in rna_info_dict. Which is a line in genepred file'''
    rna_name, chrom, strand, start, end, cds_s, cds_e, exn_n, \
                exn_starts, exn_ends, score, name, cds_s_stat, cds_e_stat,frame = rna_info
    cds_s = int(cds_s)
    cds_e = int(cds_e)
    exn_starts = [int(p) for p in exn_starts.split(',')[:-1]]
    exn_ends =   [int(p) for p in exn_ends.split(',')[:-1]]

    inter_rna = IntervalTree()
    for s, e in zip(exn_starts, exn_ends):
        inter_rna[s:e] = 'exon'
    if cds_s + 1 == cds_e:  # if it's non coding rna, then the regions would be annotated as exon, we return it.
        return inter_rna
    
    # if it's protein coding rna, continue further slice the interval
    inter_rna.slice(cds_s)
    inter_rna.slice(cds_e)
    if strand == '+':
        start_anno = 'UTR5'
        end_anno = 'UTR3'
    else:
        start_anno = 'UTR3'
        end_anno = 'UTR5'
    anno = start_anno
    new_inter_rna = IntervalTree()
    # add functional annotation to interval: eg. UTR5, CDS, UTR3
    for inter in sorted(inter_rna):
        if inter[0] == cds_s:
            anno = 'CDS'
        elif inter[0] == cds_e:
            anno = end_anno
        new_inter_rna[inter[0]:inter[1]] = anno
    new_intervals = sorted(new_inter_rna)
    for i in range(len(new_intervals)-1):
        if new_intervals[i][1] < new_intervals[i+1][0]:
            new_inter_rna[new_intervals[i][1]:new_intervals[i+1][0]] = 'intron'
    return new_inter_rna


def indel_vcf_cigar(cigar):
    '''this function process cigar string of insertion and deletion into a list of integers
    the list should have 3 numbers. The middle number represents the indel nucleotides, the
    other represents the matched nucleotides that are the same with reference
    * cigar: eg: 1M3I5M, this meanns the middle 3 base pairs are inserted'''
    match = re.findall(r'(\d+)([A-Z]{1})', cigar)
    if len(match) == 3:
        res = [int(i[0]) for i in match]
    elif len(match) == 2:
        if match[0][0] in ['I', 'D']:
            res = [0] + [int(i[0]) for i in match]
        else:
            res = [int(i[0]) for i in match] + [0]
    else:
        res = [0,0,0]
    return res


def get_vcf_type_in_cds(rna_tree, rna_info, rna_seq_index, vcf_info):
    '''this function analyze mutations locate in CDS regions. It will get the type of mutation (eg, synonymous, nonsynonymous)
    and provide the position, nucleotide and amino acide change
    * rna_tree: all intervals of a rna
    * rna_info: a line in genepred file that has all position information of an rna
    * rna_seq_index: SeqIO object that have access to the full sequence of all the rnas'''
    for inter in sorted(rna_tree):
        if inter[2] == 'intron':
            rna_tree.remove(inter)
    rna_name, chrom, strand, start, end, cds_s, cds_e, exn_n, \
            exn_s, exn_e, score, name, cds_s_stat, cds_e_stat,frame = rna_info # rna_info_dic[rna]
    
    vcf_start = int(vcf_info[1]) - 1
    vcf_end = vcf_start + len(vcf_info[3])
    cigar = indel_vcf_cigar(re.search('(?<=CIGAR=).+?(;)', vcf_info[7]).group(0))
    if strand == '-':
        cigar = cigar[::-1]
    rna_tree.slice(vcf_start)
    rna_tree.slice(vcf_end)
    # calculate positions of CDS and mutation relative to TSS of mRNA sequence
    vcfStartInRNA = 0
    vcfEndInRNA = 0
    cdsStartInRNA = 0
    cdsEndInRNA = 0
    if strand == '-':
        for inter in sorted(rna_tree)[::-1]:
            if inter[0] >= vcf_start:
                vcfEndInRNA += inter[1] - inter[0]
            if inter[0] >= vcf_end:
                vcfStartInRNA += inter[1] - inter[0]
            if inter[2] == 'UTR5':
                cdsStartInRNA += inter[1] - inter[0]
    else:
        for inter in sorted(rna_tree):
            if inter[1] <= vcf_start:
                vcfStartInRNA += inter[1] - inter[0]
            if inter[1] <= vcf_end:
                vcfEndInRNA += inter[1] - inter[0]
            if inter[2] == 'UTR5':
                cdsStartInRNA += inter[1] - inter[0]
    vcfStartInCDS = vcfStartInRNA - cdsStartInRNA
    # ----------------------  1). process SNP
    if len(vcf_info[3]) == len(vcf_info[4]) == 1:
        snp_reverse_dic = {'A':'T','C':'G','G':'C','T':'A'}
        ref = vcf_info[3] # reference dna sequence
        alt = vcf_info[4] # alternative dna sequence
        if strand == '-': # reverse strand need to get complement of the nucleotides
            ref = snp_reverse_dic[ref]
            alt = snp_reverse_dic[alt]
        frame = vcfStartInCDS % 3 # calculate which frame does the snp locates
        aaStartInCDS = vcfStartInCDS // 3
        if frame == 0:
            ref_seq = rna_seq_index[rna_name].seq[vcfStartInRNA:vcfStartInRNA + 3]
            alt_seq = alt + ref_seq[1:]
        elif frame == 1:
            ref_seq = rna_seq_index[rna_name].seq[vcfStartInRNA-1:vcfStartInRNA + 2]
            alt_seq = ref_seq[0] + alt + ref_seq[2]
        elif frame == 2:
            ref_seq = rna_seq_index[rna_name].seq[vcfStartInRNA-2:vcfStartInRNA + 1]
            alt_seq = ref_seq[:2] + alt
        if isinstance(ref_seq, str):
            ref_seq = Seq(ref_seq)
        if isinstance(alt_seq, str):
            alt_seq = Seq(alt_seq)
        ref_aa = str(ref_seq.translate())  # get reference amino acid
        alt_aa = str(alt_seq.translate())
        vcf_detail = 'c.' + str(ref) + str(vcfStartInCDS + 1) + str(alt) + ':p.' + str(ref_aa) + str(aaStartInCDS + 1) + str(alt_aa)
        if ref_aa == alt_aa:
            vcf_type = 'synonymous SNV'
        else:
            vcf_type = 'non synonymous SNV'
        if ref_aa == '*':
            vcf_type = 'stop loss'
        if alt_aa == '*':
            vcf_type = 'stop gain'
    #-------------------------- 2).process insertion
    elif len(vcf_info[3]) < len(vcf_info[4]):
        insStartInCDS = vcfStartInCDS + cigar[0] #len(vcf_info[3])
        rna_full_seq = rna_seq_index[rna_name].seq
        ins_seq = vcf_info[4][cigar[0]:cigar[0] + cigar[1]]
        # if in reverse strand, get reverse complement of the sequence
        if strand == '-':
            ins_seq = str(Seq(ins_seq).reverse_complement())
        # include 3 nt before insertion and 3 nt after insertion
        insStartInRNA = vcfStartInRNA + cigar[0]
        # insertions can be inserted in any frame', so need to consider them separately
        if len(ins_seq) % 3 == 0: # it's non frame shifting insertion
            frame = (insStartInCDS - 1)  % 3
            if frame == 0:
                ref_dna = rna_full_seq[insStartInRNA-1 : insStartInRNA+2]
                alt_dna = rna_full_seq[insStartInRNA-1 : insStartInRNA] + ins_seq + rna_full_seq[insStartInRNA : insStartInRNA+2]
            elif frame == 1:
                ref_dna = rna_full_seq[insStartInRNA-2 : insStartInRNA+1]
                alt_dna = rna_full_seq[insStartInRNA-2 : insStartInRNA] + ins_seq + rna_full_seq[insStartInRNA : insStartInRNA+1]
            elif frame == 2:
                ref_dna = rna_full_seq[insStartInRNA-3 : insStartInRNA]
                alt_dna = rna_full_seq[insStartInRNA-3 : insStartInRNA] + ins_seq
            ref_aa = ref_dna.translate()
            alt_aa = alt_dna.translate()
        
            aaStartCDS = insStartInCDS // 3 + 1
            vcf_detail = 'c.' + str(insStartInCDS+1) + '_' + str(insStartInCDS+2) + ':p.' + str(ref_aa) + str(aaStartCDS) + 'delins' + str(alt_aa)
            vcf_type = 'non frame shift insertion'
        else:  # frame shifting insertion
            frame = (insStartInCDS - 1) % 3
            if frame == 0:
                ref_dna = rna_full_seq[insStartInRNA-1 : insStartInRNA+2]
            elif frame == 1:
                ref_dna = rna_full_seq[insStartInRNA-2 : insStartInRNA+1]
            elif frame == 2:
                ref_dna = rna_full_seq[insStartInRNA-3 : insStartInRNA]
            ref_aa = str(ref_dna.translate())
            aaStartCDS = insStartInCDS // 3 + 1
            vcf_detail = 'c.' + str(insStartInCDS+1) + '_' + str(insStartInCDS+2) + 'ins' + str(ins_seq) + \
                        ':p.' + str(ref_aa) + str(aaStartCDS) + 'fs'
            vcf_type = 'frame shift insertion'
    #--------------------------- 3). process deletion
    elif len(vcf_info[3]) > len(vcf_info[4]):
        delStartInCDS = vcfStartInCDS + cigar[0]
        delEndInCDS = vcfStartInCDS + cigar[0] + cigar[1]
        rna_full_seq = rna_seq_index[rna_name].seq
        del_seq = vcf_info[3][cigar[0]:cigar[0]+cigar[1]]
        if strand == '-':
            del_seq = str(Seq(del_seq).reverse_complement())
        delStartInRNA = vcfStartInRNA + cigar[0]
        delEndInRNA = vcfStartInRNA + cigar[0] + cigar[1]
        if len(del_seq) % 3 == 0: # non frame shift deletion
            frame = delStartInCDS % 3 
            if frame == 0:
                ref_dna = del_seq
                alt_dna = ''
            elif frame == 1:
                ref_dna = rna_full_seq[delStartInRNA-1:delStartInRNA]+del_seq+rna_full_seq[delEndInRNA:delEndInRNA+2]
                alt_dna = rna_full_seq[delStartInRNA-1:delStartInRNA] + rna_full_seq[delEndInRNA:delEndInRNA+2]
            elif frame == 2:
                ref_dna = rna_full_seq[delStartInRNA-2:delStartInRNA]+del_seq+rna_full_seq[delEndInRNA:delEndInRNA+1]
                alt_dna = rna_full_seq[delStartInRNA-2:delStartInRNA] + rna_full_seq[delEndInRNA:delEndInRNA+1]
            if isinstance(ref_dna, str):
                ref_dna = Seq(ref_dna)
            if isinstance(alt_dna, str):
                alt_dna = Seq(alt_dna)
            ref_aa = ref_dna.translate()
            alt_aa = alt_dna.translate()
            aaStartCDS = delStartInCDS // 3 + 1
            aaEndCDS = aaStartCDS + len(ref_aa)
            vcf_detail = 'c.' + str(delStartInCDS+1) + '_' + str(delEndInCDS) + 'del:p.' + str(aaStartCDS) + '_' + \
                        str(aaEndCDS) + str(ref_aa) + 'del'
            if str(alt_aa) != '':
                vcf_info += 'ins' + str(alt_aa)
            vcf_type = 'non frame shift deletion'
        else: # frame shift deletion
            frame = delStartInCDS % 3
            offset = (delEndInCDS - delStartInCDS) % 3
            if frame == 0:
                ref_dna = del_seq + rna_full_seq[delEndInRNA:delEndInRNA+3-offset]
            elif frame == 1:
                ref_dna = rna_full_seq[delStartInRNA-1:delStartInRNA]+del_seq+rna_full_seq[delEndInRNA:delEndInRNA+2]
            elif frame == 2:
                ref_dna = rna_full_seq[delStartInRNA-2:delStartInRNA]+del_seq+rna_full_seq[delEndInRNA:delEndInRNA+1]
            ref_aa = ref_dna.translate()
            aaStartCDS = delStartInCDS // 3 + 1
            aaEndCDS = aaStartCDS + len(ref_aa) - 1 
            vcf_detail = 'c.'+ str(delStartInCDS+1) +'_' + str(delEndInCDS) + 'del:p.' + str(aaStartCDS) + '_' + \
                        str(aaEndCDS) + str(ref_aa) + 'delfs'  
            vcf_type = 'frame shift deletion'
    #------------------------------- 4) equal length replacement
    elif len(vcf_info[3]) == len(vcf_info[4]) and len(vcf_info[3]) > 1 and len(vcf_info[4]) > 1:
        rna_full_seq = rna_seq_index[rna_name].seq
        ref_seq = Seq(vcf_info[3])
        alt_seq = Seq(vcf_info[4])
        if strand == '-':
            ref_seq = ref_seq.reverse_complement()
            alt_seq = alt_seq.reverse_complement()
        frame = vcfStartInCDS % 3
        offset = len(ref_seq) % 3
        vcfEndInCDS = vcfStartInCDS + len(ref_seq)
        if frame == 0:
            ref_dna = ref_seq + rna_full_seq[vcfEndInRNA:vcfEndInRNA+3-offset]
            alt_dna = alt_seq + rna_full_seq[vcfEndInRNA:vcfEndInRNA+3-offset]
        elif frame == 1:
            ref_dna = rna_full_seq[vcfStartInRNA-1:vcfStartInRNA] + ref_seq + rna_full_seq[vcfEndInRNA:vcfEndInRNA+3-offset-1]
            alt_dna = rna_full_seq[vcfStartInRNA-1:vcfStartInRNA] + alt_seq + rna_full_seq[vcfEndInRNA:vcfEndInRNA+3-offset-1]
        elif frame == 2:
            ref_dna = rna_full_seq[vcfStartInRNA-2:vcfStartInRNA] + ref_seq + rna_full_seq[vcfEndInRNA:vcfEndInRNA+3-offset+1]
            alt_dna = rna_full_seq[vcfStartInRNA-2:vcfStartInRNA] + alt_seq + rna_full_seq[vcfEndInRNA:vcfEndInRNA+3-offset+1]
        ref_aa = ref_dna.translate()
        alt_aa = alt_dna.translate()
        aaStartCDS = vcfStartInCDS // 3 + 1
        aaEndCDS = aaStartCDS + len(ref_aa) - 1
        vcf_detail = 'c.' + str(ref_seq) + str(vcfStartInCDS + 1) + '_' + str(vcfEndInCDS) + str(alt_seq) + ':p.' + \
                        str(ref_aa) + str(aaStartCDS) + '_' + str(aaEndCDS) + str(alt_aa)
        if ref_aa == alt_aa:
            vcf_type = 'synonymous SNV'
        else:
            vcf_type = 'non synonymous SNV'
    return vcf_detail, vcf_type


def deleterious_level(vcf_type):
    '''this function assigns each variation a level, smaller number means more deleterious'''
    if vcf_type in ['non synonymous SNV', 'frame shift insertion', 'frame shift deletion','stop gain','stop loss']:
        return 1
    elif vcf_type in ['non frame shift insertion', 'non frame shift deletion']:
        return 2
    elif vcf_type in ['synonymous SNV','UTR5','UTR3','intron','intergenic']:
        return 3
    else:
        raise 'vcf type not recognized'


def get_vcf_type_for_one_mutation(vcf_line, trees, rna_seq_index):
    '''most of the time each line in vcf file only has one mutation, but sometimes there may exist more than one
    In the latter case, there is another function to split them, then pass to this function to parse one mutation at a time'''
    # read througheach line of vcf file
    items = vcf_line.strip().split()
    vcf_chrom, vcf_start = items[0], int(items[1]) - 1
    vcf_end = vcf_start + len(items[3])
    # get rnas intervals that overlap with the mutation
    rna_inters = sorted(trees[vcf_chrom][vcf_start:vcf_end])
    vcf_anno = []
    if rna_inters == []:  # intergenic regions
        vcf_type = 'intergenic'
        vcf_detail = '.'
        final_vcf_anno = [['.','.','intergenic']]
    else:
        # process for each rna that has that mutation
        for rna_inter in rna_inters:
            # 1). build rna tree
            rna = rna_inter[2]
            rna_tree = build_rna_tree(rna_info_dic[rna])
            rna_name, chrom, strand, start, end, cds_s, cds_e, exn_n, \
                exn_s, exn_e, score, name, cds_s_stat, cds_e_stat,frame = rna_info_dic[rna]
            # 2). check which region of rna 
            rna_regions_inters = sorted(rna_tree[vcf_start:vcf_end])
            if len(rna_regions_inters) == 1:
                vcf_type = rna_regions_inters[0][2]
                vcf_detail = ''
            else:
                # if it maps to more than one regions, it must be a deletion
                # the severity order is CDS > UTR > intron
                types = [p[2] for p in rna_regions_inters]
                if 'CDS' in types:
                    vcf_type = 'CDS'
                elif 'UTR5' in types:
                    vcf_type = 'UTR5'
                elif 'UTR3' in types:
                    vcf_type = 'UTR3'
                else:
                    vcf_type = 'intron'
                vcf_detail = ''
            # if locate in CDS, further parse the results
            if vcf_type == 'CDS': 
                vcf_detail, vcf_type = get_vcf_type_in_cds(rna_tree, rna_info_dic[rna], rna_seq_index, items)
            vcf_anno.append([rna_name, vcf_detail, vcf_type])
        for anno in vcf_anno:
            anno.append(deleterious_level(anno[2]))
        vcf_type_best_score = min([deleterious_level(s[2]) for s in vcf_anno])
        final_vcf_anno = [v for v in vcf_anno if v[3] == vcf_type_best_score]
    return final_vcf_anno


def split_multi_mutations(vcf_line):
    '''this function split vcf lines that have more than one mutations, it creates a line for each mutation'''
    items = vcf_line.strip().split('\t')
    n = len(items[4].strip().split(','))
    lines = [[''] * len(items) for i in range(n)]
    single_index = [0, 1, 2, 3, 5, 6, 8, 9, 10]
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


def coverage_at_mutation_site(line):
    '''this function gets some statistics about coverage for the mutations
    it will return 6 elements: [normal sample total read coverage,      patient sample total read coverage,
                                normal sample mutation read coverage,   patient sample mutation read coverage,
                                normal sample mutation read percentage, patient sample mutation read percentage]'''
    items = line.strip().split('\t')
    normal = items[-2].split(':')
    vaf5 = items[-1].split(':')
    # total coverage at mutation loci
    normal_cov = normal[2]
    vaf5_cov = vaf5[2]
    # reads support mutation
    normal_mut_cov = normal[6]
    vaf5_mut_cov = vaf5[6]
    # percentage of reads supporting the variants vS reads supporting reference
    normal_mut_covs = [int(n) for n in normal_mut_cov.split(',')]
    vaf5_mut_covs = [int(n) for n in vaf5_mut_cov.split(',')]
    normal_mut_per = ','.join([str(round(n / int(normal_cov) * 100, 2)) + '%' for n in normal_mut_covs])
    normal_ref_per = str(round(100 - sum(normal_mut_covs) / int(normal_cov) * 100, 2)) + '%'
    vaf5_mut_per = ','.join([str(round(n / int(vaf5_cov) * 100, 2)) + '%' for n in vaf5_mut_covs])
    vaf5_ref_per = str(round(100 - sum(vaf5_mut_covs) / int(vaf5_cov) * 100, 2)) + '%'
    return [normal_cov, vaf5_cov, normal_mut_cov, vaf5_mut_cov, 
            normal_ref_per+','+normal_mut_per, vaf5_ref_per+','+vaf5_mut_per]
   



parser = argparse.ArgumentParser(description='annotate vcf file')
parser.add_argument('-i','--input',action='store',dest='input',help='input vcf file')
parser.add_argument('-p','--pred',action='store',dest='pred',help='pred file format of gff file')
parser.add_argument('-f','--rnaseq',action='store',dest='rnaseq',help='full sequence of mRNA')
parser.add_argument('-o','--out',action='store',dest='out',help='out file name')

args = parser.parse_args()

if __name__ == '__main__':
    # anno_vcf_fn = '/data/shangzhong/temp/anno_vcf_anno.vcf'
    # pred_fn = '/data/shangzhong/temp/human_g1k_v37.pred'
    # vcf_fn = '/data/shangzhong/temp/Challenge_data.vcf'
    # rna_seq_fn = '/data/shangzhong/temp/human_g1k_v37.rna.fa'
    # build tree and rna information
    vcf_fn = args.input
    pred_fn = args.pred
    rna_seq_fn = args.rnaseq
    anno_vcf_fn = args.out
    

    trees, rna_info_dic = build_interval_tree(pred_fn)
    # index mRNA seqeunce fasta file
    rna_seq_index = SeqIO.index(rna_seq_fn,'fasta')
    # annotation the file
    with open(vcf_fn) as f, open(anno_vcf_fn, 'w') as out:
        out.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'normal', 'vaf5',
                             'variat_type','variant_detail','normal_cov', 'vaf5_cov','normal_mut_cov',
                             'vaf5_mut_cov','normal_mut_cov_percent','vaf5_mut_cov_percent']) + '\n')
        for line in f:
            if line.startswith('#'):
                continue
            #------------------ Step 1: get mutation type and mutation detials
            # make sure each mutation is in a separate line
            lines = split_multi_mutations(line)
            vcf_type = []
            full_annotation = []
            for l in lines:
                vcf_annos = get_vcf_type_for_one_mutation(l, trees, rna_seq_index)
                types = []
                annos = []
                for a in vcf_annos:
                    types.append(a[2])
                    annos.append(':'.join(a[:2]))
                    if annos[-1] == '.:.': annos[-1] = '.'
                vcf_type.extend(types)
                full_annotation.append(';'.join(annos))
            #------------------ Step 2 & 3 & 4: get depth of sequence coverage at he site of variation and reads support mutation
            coverage = coverage_at_mutation_site(line)
            # write to file
            out.write(line.strip() + '\t' + ','.join(list(set(vcf_type))) + '\t' + ','.join(full_annotation)  + '\t' + \
                      '\t'.join(coverage) + '\n')
