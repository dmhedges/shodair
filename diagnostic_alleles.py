import pickle
import pandas as pd
from pprint import pprint
from collections import defaultdict


def get_all_alleles(build=37):
    if build == 37:
        uri = '/Users/david.hedges/referenceAndData/pharmvar-4.0.2/GRCh37_all/combined37.tsv'
    if build == 38:
        uri = '/Users/david.hedges/referenceAndData/pharmvar-4.0.2/GRCh38_all/combined38.tsv'
    all_alleles = defaultdict(list)

    df = pd.read_csv(uri, sep='\t', na_filter=None)
    df = df.drop(['ReferenceSequence'], axis=1)
    df.rename(columns={'Haplotype Name': 'HaplotypeName', 'Variant Start': 'VariantStart', 
                   'Variant Stop':'VariantStop','Reference Allele':'ReferenceAllele',
                   'Variant Allele':'VariantAllele'}, inplace=True)
    #df = df[~df.HaplotypeName.str.contains('.0')]
    df = df.reset_index().drop(['index'], axis=1)
    allel_dict = df.to_dict('index')

    var_list = []
    for k,v in allel_dict.items():
        for key,val in v.items():
            gene = v['Gene']
            allele = v['HaplotypeName']
            rsID = v['rsID']
            ref = v['ReferenceAllele']
            var = v['VariantAllele']
            start = v['VariantStart']
            stop = v['VariantStop']
            vtype = v['Type']
            var_list.append((allele,gene,rsID,ref,var,start,stop,vtype))
    var_set = set(var_list)
    for j in var_set:
        all_alleles[j[0]].append(j[1:])
    return all_alleles


def remove_suballeles(build=37):
    alleles = defaultdict(list)
    all_alleles = get_all_alleles(build=build)
    for i in list(all_alleles.keys()):
        try:
            allele = i.split('.')[0]
            suballele = i.split('.')[1]
        except:
            allele = i
            suballele = ''
        variant_info = all_alleles[i]
        alleles[allele].append(variant_info)
    return alleles


def get_diagnostic_alleles(build=37):
    final_diagnostic = defaultdict(list)
    alleles = remove_suballeles(build=build)
    for k,v in alleles.items():
        counter=0
        intersected = v[0]
        if len(v) >= 2:
            for i in v:
                a = intersected
                b = v[counter]
                counter+=1        
                intersected = set(a).intersection(b)
        final_diagnostic[k].append(intersected)
    return final_diagnostic