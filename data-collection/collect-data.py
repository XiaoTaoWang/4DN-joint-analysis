"""
This script is able to download data from ENCODE, taking into account
the presence of an unstable internet connection.

Both input files metadata.tsv and experiment_report_2023_4_10_0h_40m.tsv were downloaded from ENCODE data portal.

"""

import subprocess, joblib, os, time, sys
from collections import Counter

def start_download(command):

    download_process = subprocess.run(' '.join(command), shell=True)

    return download_process

def run_aria2c_recursive(dest, url):

    try:
        command = [
                    'aria2c', '-x 5', '-c',
                    '-d', os.path.abspath(os.path.split(dest)[0]),
                    '-o', os.path.split(dest)[1], url
                ]
        subprocess.check_call(' '.join(command), shell=True)
    except:
        print('Error occurred, try to re-connect ...')
        run_aria2c_recursive(dest, url)

def run_aria2c(dest, url):

    while True:
        command = [
                    'aria2c', '-x 5', '-c',
                    '-d', os.path.abspath(os.path.split(dest)[0]),
                    '-o', os.path.split(dest)[1], url
                ]
        download_process = start_download(command)
        if download_process.returncode == 0:
            break


def read_header(metadata_fil):

    with open(metadata_fil, 'r') as source:
        header = source.readline().rstrip().split('\t')
    
    header_map = dict(zip(header, range(len(header))))

    return header_map

def check_genetic_info(line, header_map):

    parse = line.rstrip().split('\t')
    if parse[header_map['Biosample genetic modifications methods']] or \
       parse[header_map['Biosample genetic modifications categories']] or \
       parse[header_map['Biosample genetic modifications targets']] or \
       parse[header_map['Biosample genetic modifications gene targets']] or \
       parse[header_map['Biosample genetic modifications site coordinates']] or \
       parse[header_map['Biosample genetic modifications zygosity']]:
        return True
    else:
        return False   

def encode_treatment_duration(string):

    tmp = string.split(', ')
    fields = []
    for t in tmp:
        parse = t.split()
        if len(parse)==1:
            fields.append(parse[0])
        elif parse[1]=='day':
            fields.append('{0}d'.format(parse[0]))
        elif parse[1]=='hour':
            fields.append('{0}h'.format(parse[0]))
    
    return '-'.join(fields)

def read_experiment_info(fil):

    D = {}
    with open(fil, 'r') as source:
        source.readline()
        header = source.readline().rstrip().split('\t')
        idx1 = header.index('Accession')
        idx2 = header.index('Cellular component')
        for line in source:
            parse = line.split('\t')
            ID = parse[idx1]
            cell_component = parse[idx2]
            D[ID] = cell_component
    
    return D

experiment_info = read_experiment_info('experiment_report_2023_4_10_0h_40m.tsv')

included_assays = {
    'ATAC-seq': 'ATAC-seq',
    'DNase-seq': 'DNase-seq',
    'Bru-seq': 'Bru-seq', # nascent transcripts
    'ChIA-PET': 'ChIA-PET',
    'Histone ChIP-seq': 'ChIP-seq',
    'Mint-ChIP-seq': 'Mint-ChIP-seq',
    'TF ChIP-seq': 'ChIP-seq',
    'WGBS': 'WGBS',
    'in situ Hi-C': 'in-situ-HiC',
    'intact Hi-C': 'INT-HiC', # omit the biotin incorporation, only need a few hundred cells
    'polyA plus RNA-seq': 'polyA-RNASeq',
    'total RNA-seq': 'total-RNASeq',
    'snATAC-seq': 'snATAC-seq'
}

cell_name_map = {
    'B cell': 'B-cell',
    'CD4-positive, alpha-beta T cell': 'T-cell-CD4+',
    'CD4-positive, alpha-beta memory T cell': 'memory-T-cell-CD4+',
    'CD8-positive, alpha-beta T cell': 'T-cell-CD8+',
    'CD8-positive, alpha-beta memory T cell': 'memory-T-cell-CD8+',
    'MCF 10A': 'MCF10A',
    'activated B cell': 'B-cell',
    'activated CD4-positive, alpha-beta T cell': 'T-cell-CD4+',
    'activated CD4-positive, alpha-beta memory T cell': 'memory-T-cell-CD4+',
    'activated CD8-positive, alpha-beta T cell': 'T-cell-CD8+',
    'activated CD8-positive, alpha-beta memory T cell': 'memory-T-cell-CD8+',
    'activated T-cell': 'T-cell',
    'activated naive CD4-positive, alpha-beta T cell': 'naive-T-cell-CD4+',
    'activated naive CD8-positive, alpha-beta T cell': 'naive-T-cell-CD8+',
    'common myeloid progenitor, CD34-positive': 'myeloid-progenitor-CD34+',
    'endothelial cell of umbilical vein': 'HUVEC',
    'mammary epithelial cell': 'HMEC',
    'naive thymus-derived CD4-positive, alpha-beta T cell': 'naive-T-cell-CD4+',
    'naive thymus-derived CD8-positive, alpha-beta T cell': 'naive-T-cell-CD8+'
}

treatment_map = {
    'CpG ODN': 'CpG-ODN',
    'anti-CD3 and anti-CD28 coated beads, Interleukin-2': 'IL2-TCR',
    'anti-CD3 and anti-CD28 coated beads': 'TCR',
    'DMSO': 'DMSO'
}

header_map = read_header('metadata.tsv')

collect = {}
biotype_map = {}
with open('metadata.tsv', 'r') as source:
    source.readline()
    for line in source:
        if check_genetic_info(line, header_map):
            continue
        parse = line.rstrip().split('\t')
        genome_assembly = parse[header_map['File assembly']]
        ori_cell_name = parse[header_map['Biosample term name']]
        ID = parse[header_map['File accession']]
        output_type = parse[header_map['Output type']]
        ftype = parse[header_map['File format']]
        experiment_ID = parse[header_map['Experiment accession']]
        assay = parse[header_map['Assay']]
        target = parse[header_map['Experiment target']].split('-')[0]
        reps = parse[header_map['Biological replicate(s)']]
        url = parse[header_map['File download URL']]
        #lab = parse[header_map['Lab']]
        bio_type = parse[header_map['Biosample type']]
        bio_treatment = parse[header_map['Biosample treatments']]
        treatment_duration = parse[header_map['Biosample treatments duration']]
        
        if experiment_info[experiment_ID]:
            if not assay in ['snATAC-seq']:
                continue # skip RNA-seq data that are from sub-nuclear fractions
        
        if bio_treatment:
            if not bio_treatment in treatment_map:
                continue
        
        rep_label = '_'.join(reps.split(', '))
        if assay in included_assays:
            if ori_cell_name in cell_name_map:
                cell_name = cell_name_map[ori_cell_name]
            else:
                cell_name = ori_cell_name
            
            if not bio_treatment:
                sample_name = cell_name
            else:
                sample_name = '_'.join([cell_name, treatment_map[bio_treatment], encode_treatment_duration(treatment_duration)])
            if not sample_name in collect:
                collect[sample_name] = {}
                biotype_map[sample_name] = bio_type
            
            if assay in ['ATAC-seq']:
                # will download signal tracks and peaks
                if (ftype == 'bigWig') and (output_type == 'signal p-value'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.pvalue_signal{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)
                elif (ftype == 'bigWig') and (output_type == 'fold change over control'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.fc_signal{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)
                elif (ftype == 'bed narrowPeak') and (output_type == 'pseudoreplicated peaks'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.peaks{1}.{2}'.format(ID, rep_label, 'bed.gz')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)

            elif assay in ['DNase-seq']:
                # will download signal tracks, peaks, and footprints
                if (ftype == 'bigWig') and (output_type == 'read-depth normalized signal'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.normalized_signal{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)
                elif (ftype == 'bed narrowPeak') and (output_type == 'peaks'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.peaks{1}.{2}'.format(ID, rep_label, 'bed.gz')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)
                elif (ftype == 'bed bed3+') and (output_type == 'footprints'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.footprints{1}.{2}'.format(ID, rep_label, 'bed.gz')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)
            
            elif assay in ['Bru-seq', 'polyA plus RNA-seq', 'total RNA-seq']:
                # will download the plus/minus signals in bigwig, and gene quantifications in tsv
                if (ftype == 'bigWig') and (output_type == 'plus strand signal of unique reads'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.plus_strand{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)
                elif (ftype == 'bigWig') and (output_type == 'minus strand signal of unique reads'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.minus_strand{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)
                elif (ftype == 'tsv') and (output_type == 'gene quantifications'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.gene_quant{1}.{2}'.format(ID, rep_label, 'tsv')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)
            
            elif assay in ['ChIA-PET']:
                # will download the contact matrix, 1D signals in bigwig, loops, and peaks
                assay_ = '{0} {1}'.format(target, assay)
                if (ftype == 'bigWig') and (output_type == 'signal of unique reads'):
                    if not assay_ in collect[sample_name]:
                        collect[sample_name][assay_] = {}
                    if not experiment_ID in collect[sample_name][assay_]:
                        collect[sample_name][assay_][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.signal{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay_][experiment_ID][key] = (url, dest)
                elif (ftype == 'hic') and (output_type == 'contact matrix'):
                    if not assay_ in collect[sample_name]:
                        collect[sample_name][assay_] = {}
                    if not experiment_ID in collect[sample_name][assay_]:
                        collect[sample_name][assay_][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.contact_matrix{1}.{2}'.format(ID, rep_label, 'hic')
                    collect[sample_name][assay_][experiment_ID][key] = (url, dest)
                elif (ftype == 'bed broadPeak') and (output_type == 'peaks'):
                    if not assay_ in collect[sample_name]:
                        collect[sample_name][assay_] = {}
                    if not experiment_ID in collect[sample_name][assay_]:
                        collect[sample_name][assay_][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.peaks{1}.{2}'.format(ID, rep_label, 'bed.gz')
                    collect[sample_name][assay_][experiment_ID][key] = (url, dest)
                elif (ftype == 'bed narrowPeak') and (output_type == 'peaks'):
                    if not assay_ in collect[sample_name]:
                        collect[sample_name][assay_] = {}
                    if not experiment_ID in collect[sample_name][assay_]:
                        collect[sample_name][assay_][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.peaks{1}.{2}'.format(ID, rep_label, 'bed.gz')
                    collect[sample_name][assay_][experiment_ID][key] = (url, dest)
                elif (ftype == 'bedpe') and (output_type == 'loops'):
                    if not assay_ in collect[sample_name]:
                        collect[sample_name][assay_] = {}
                    if not experiment_ID in collect[sample_name][assay_]:
                        collect[sample_name][assay_][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.loops{1}.{2}'.format(ID, rep_label, 'bedpe.gz')
                    collect[sample_name][assay_][experiment_ID][key] = (url, dest)

            elif assay in ['Histone ChIP-seq', 'Mint-ChIP-seq']:
                # will download signal tracks and peaks
                if assay == 'Histone ChIP-seq':
                    assay_ = '{0} {1}'.format(target, 'ChIP-seq')
                else:
                    assay_ = '{0} {1}'.format(target, 'Mint-ChIP-seq')
                if (ftype == 'bigWig') and (output_type == 'signal p-value'):
                    if not assay_ in collect[sample_name]:
                        collect[sample_name][assay_] = {}
                    if not experiment_ID in collect[sample_name][assay_]:
                        collect[sample_name][assay_][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.pvalue_signal{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay_][experiment_ID][key] = (url, dest)
                elif (ftype == 'bigWig') and (output_type == 'fold change over control'):
                    if not assay_ in collect[sample_name]:
                        collect[sample_name][assay_] = {}
                    if not experiment_ID in collect[sample_name][assay_]:
                        collect[sample_name][assay_][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.fc_signal{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay_][experiment_ID][key] = (url, dest)
                elif (ftype == 'bed narrowPeak') and (output_type == 'pseudoreplicated peaks'):
                    if not assay_ in collect[sample_name]:
                        collect[sample_name][assay_] = {}
                    if not experiment_ID in collect[sample_name][assay_]:
                        collect[sample_name][assay_][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.peaks{1}.{2}'.format(ID, rep_label, 'bed.gz')
                    collect[sample_name][assay_][experiment_ID][key] = (url, dest)
            
            elif assay in ['TF ChIP-seq']:
                # will download signal tracks and peaks
                assay_ = '{0} {1}'.format(target, 'ChIP-seq')
                if (ftype == 'bigWig') and (output_type == 'signal p-value'):
                    if not assay_ in collect[sample_name]:
                        collect[sample_name][assay_] = {}
                    if not experiment_ID in collect[sample_name][assay_]:
                        collect[sample_name][assay_][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.pvalue_signal{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay_][experiment_ID][key] = (url, dest)
                elif (ftype == 'bigWig') and (output_type == 'fold change over control'):
                    if not assay_ in collect[sample_name]:
                        collect[sample_name][assay_] = {}
                    if not experiment_ID in collect[sample_name][assay_]:
                        collect[sample_name][assay_][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.fc_signal{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay_][experiment_ID][key] = (url, dest)
                elif (ftype == 'bed narrowPeak') and (output_type == 'IDR thresholded peaks'):
                    if not assay_ in collect[sample_name]:
                        collect[sample_name][assay_] = {}
                    if not experiment_ID in collect[sample_name][assay_]:
                        collect[sample_name][assay_][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.peaks{1}.{2}'.format(ID, rep_label, 'bed.gz')
                    collect[sample_name][assay_][experiment_ID][key] = (url, dest)
            
            elif assay in ['WGBS']:
                # will download the CpG coverage and methylation state in bigwig
                if (ftype == 'bigWig') and (output_type == 'plus strand methylation state at CpG'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.plus_CpG_methyl{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)
                elif (ftype == 'bigWig') and (output_type == 'minus strand methylation state at CpG'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.minus_CpG_methyl{1}.{2}'.format(ID, rep_label, 'bw')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)
                elif (ftype == 'bigWig') and (output_type == 'CpG sites coverage'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.CpG_cov{1}.{2}'.format(ID, rep_label, 'bw')
            
            elif assay in ['snATAC-seq']:
                if (ftype == 'tar') and (output_type == 'fragments'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.fragments{1}.{2}'.format(ID, rep_label, 'tar.gz')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)
            
            elif assay in ['intact Hi-C', 'in situ Hi-C']:
                # download contact matrices
                if (ftype == 'hic') and (output_type == 'mapping quality thresholded contact matrix'):
                    if not assay in collect[sample_name]:
                        collect[sample_name][assay] = {}
                    if not experiment_ID in collect[sample_name][assay]:
                        collect[sample_name][assay][experiment_ID] = {}
                    key = (output_type, rep_label)
                    dest = '{0}.contact_matrix{1}.{2}'.format(ID, rep_label, 'hic')
                    collect[sample_name][assay][experiment_ID][key] = (url, dest)


# skip data of individual replicates if multiple replicates exist
for sample in collect:
    for assay in collect[sample]:
        for experiment_ID in collect[sample][assay]:
            keys = list(collect[sample][assay][experiment_ID])
            for key in keys:
                output_type, rep_id = key
                reps = rep_id.split('_')
                if len(reps) > 1:
                    for k in keys:
                        if k in collect[sample][assay][experiment_ID]:
                            if (k[1] in reps) and (k[0] == output_type):
                                del collect[sample][assay][experiment_ID][k]

# manually remove datasets that are irrelevant
samples_to_remove = [
    'K562_DMSO_4h',
    'K562_DMSO_48h',
    'K562_DMSO_24h',
    'K562_DMSO_12h',
    'memory-T-cell-CD8+_TCR_4h',
    'memory-T-cell-CD8+_TCR_1h',
    'memory-T-cell-CD8+_IL2-TCR_4h',
    'memory-T-cell-CD8+_IL2-TCR_48h',
    'memory-T-cell-CD8+_IL2-TCR_24h',
    'memory-T-cell-CD8+_IL2-TCR_1h',
    'T-cell_TCR_4h',
    'T-cell_TCR_48h',
    'T-cell_TCR_24h',
    'T-cell_TCR_1h',
    'T-cell_IL2-TCR_48h',
    'T-cell_IL2-TCR_1h',
    'T-cell-CD4+_TCR_4h',
    'T-cell-CD4+_TCR_48h',
    'T-cell-CD4+_TCR_24h',
    'T-cell-CD4+_TCR_1h',
    'T-cell-CD4+_IL2-TCR_4h',
    'T-cell-CD4+_IL2-TCR_48h',
    'T-cell-CD4+_IL2-TCR_24h',
    'T-cell-CD4+_IL2-TCR_1h',
    'T-cell-CD4+_IL2-TCR_16h',
    'naive-T-cell-CD4+_IL2-TCR_72h',
    'naive-T-cell-CD4+_IL2-TCR_24h',
    'naive-T-cell-CD4+_IL2-TCR_5h-36h',
    'naive-T-cell-CD8+_IL2-TCR_24h',
    'naive-T-cell-CD8+_IL2-TCR_16h',
    'T-cell-CD4+_IL2-TCR_24d-14d',
    'naive-T-cell-CD8+_TCR_36h',
    'T-cell_IL2-TCR_4h',
    'T-cell_IL2-TCR_24h',
    'naive-T-cell-CD4+_IL2-TCR_16h'
]

for sample in samples_to_remove:
    del collect[sample]

# sort the rows
sample_stats = []
for sample in collect:
    sample_stats.append((biotype_map[sample], len(collect[sample]), sample))
sample_stats.sort(reverse=True)

# sort the columns
assay_list = []
assay_primary_list = []
for sample in collect:
    for assay in collect[sample]:
        assay_list.append(assay)
        if biotype_map[sample]=='primary cell':
            assay_primary_list.append(assay)

assay_counts = Counter(assay_list).most_common()
assay_primary_counts = Counter(assay_primary_list)
assay_order = [
    'CTCF ChIA-PET',
    'POLR2A ChIA-PET',
    'intact Hi-C',
    'in situ Hi-C',
    'DNase-seq',
    'ATAC-seq',
    'snATAC-seq',
    'total RNA-seq',
    'polyA plus RNA-seq',
    'Bru-seq',
    'WGBS'
]
for assay, count in assay_counts:
    if (assay.startswith('H3') or assay.startswith('H2') or assay.startswith('H4')) and (assay.split()[1]=='ChIP-seq'):
        if (not assay in assay_order) and (count > 1) and (assay_primary_counts[assay] > 0):
            assay_order.append(assay)

for assay, count in assay_counts:
    if assay.endswith('Mint-ChIP-seq'):
        if (not assay in assay_order) and (count > 1) and (assay_primary_counts[assay] > 0):
            assay_order.append(assay)

for assay, count in assay_counts:
    if assay.endswith('ChIP-seq'):
        if (not assay in assay_order) and (count > 1) and (assay_primary_counts[assay] > 0):
            assay_order.append(assay)

# write the table
with open('ENCODE-data-collection.tsv', 'w') as out:
    header = ['', ''] + assay_order
    out.write('\t'.join(header)+'\n')
    for biotype, _, sample in sample_stats:
        line = [biotype, sample]
        for assay in assay_order:
            if assay in collect[sample]:
                #line.append(str(len(collect[sample][assay])))
                line.append('\u2714')
            else:
                line.append('')
        out.write('\t'.join(line)+'\n')

# download the data
log_folder = 'log'
if not os.path.exists(log_folder):
    os.mkdir(log_folder)

for assay in assay_order:
    for sample in collect:
        if assay in collect[sample]:
            for experiment_ID in collect[sample][assay]:
                tmp = assay.split()
                assay_ID = '_'.join(tmp)
                outdir = '{0}/{1}/{2}'.format(sample, assay_ID, experiment_ID)
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                for k in collect[sample][assay][experiment_ID]:
                    url, dest = collect[sample][assay][experiment_ID][k]
                    Indicator = os.path.join(log_folder, '{0}.completed'.format(dest))
                    dest = os.path.join(outdir, dest)
                    
                    if os.path.exists(Indicator) and os.path.exists(dest):
                        continue
                    
                    run_aria2c(dest, url)

                    completed = open(Indicator, 'wb')
                    completed.close()