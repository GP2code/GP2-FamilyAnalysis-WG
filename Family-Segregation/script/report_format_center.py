#/usr/bin/env python
import pandas as pd
import openpyxl
import re
import os
import numpy as np

df=pd.read_csv(snakemake.input.tsv, sep="\t",dtype=str)
ped = pd.read_csv(snakemake.input.ped,sep='\t',dtype=str)

df=df.drop(['depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)'],axis=1)
#aggregate sample-level fields for non-comphet mode
sample_cols=['sample_id','genotype(sample,dad,mom)']
others=df.loc[df['#mode']!='comphet']
sample_count=others.groupby(['chr:pos:ref:alt','#mode','family_id'])[sample_cols].agg('|'.join).reset_index()
#deal with duplicated variants (x-linked) until find a way with slivar expr
sample_count= sample_count.drop_duplicates(subset=['chr:pos:ref:alt','family_id'],keep='last')

# de_al with comphet duplicated position per individaul (cuz slivar report the pair), drop not segregating comphet by counting samples in sample_id column
comphet=df.loc[df['#mode']=='comphet'].drop_duplicates()
li=[]
# find segregating comphet
for fam, dat in comphet.groupby('family_id'):
    # Filter the affected samples for this family
    matching_samples = ped.loc[(ped['family_id'] == fam)&(ped['phenotype']=='2'), 'id']
    # Check if we actually found any samples
    if not matching_samples.empty:
        samples = matching_samples.iloc[0]
        # grouping and aggregation
        dat = dat.groupby(['chr:pos:ref:alt', '#mode', 'family_id', 'gene'])[sample_cols].agg('|'.join).reset_index()
        # Filter dat based on the samples found
        count = dat.loc[(dat['sample_id'] == samples)]
        li.append(count)
    else:
        print(f"Warning: Family {fam} not found in ped. Skipping...")

# Check if the list contains anything
if len(li) > 0:
    comphet_count = pd.concat(li, axis=0)
    # Continue with your script (e.g., saving to CSV)
else:
    print("Warning: No data was collected. 'li' is empty.")
    # Create an empty dataframe with expected columns to prevent downstream errors
    comphet_count = pd.DataFrame()

sample_df=pd.concat([sample_count,comphet_count],axis=0).sort_values(['#mode', 'chr:pos:ref:alt'])
sample_df = sample_df.drop('gene', axis=1, errors='ignore')

drop_cols=['family_id']+sample_cols
anno=df.drop(drop_cols,axis=1).drop_duplicates(subset=['#mode','chr:pos:ref:alt'],keep='last').sort_values(['#mode', 'chr:pos:ref:alt'])

#merge back aggregate and anno fields
out=pd.merge(sample_df,anno,on=["chr:pos:ref:alt","#mode"], how='inner')

def get_max_str(lst):
    return max(lst, key=len)

# move vep annotation to the last position
column_to_move = get_max_str(out.columns)
out[column_to_move] = out.pop(column_to_move)

anno= out[get_max_str(out.columns)].str.split(';', expand=True).add_prefix('ann')
tmp=pd.concat([out,anno],axis=1)
csq_column=get_max_str(out.columns).split(';')
other_column=[ c for c in tmp.columns if not(c.startswith('ann'))]
tmp.columns=other_column+csq_column
tmp=tmp.loc[:,~tmp.columns.duplicated()]

##split var_synmoous
def make_OMIMlink(value):
    url = "https://www.omim.org/entry/{}"
    return '=HYPERLINK("%s", "%s")' % (url.format(value), value)

def make_Clinvarlink(value):
    url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/{}"
    return '=HYPERLINK("%s", "%s")' % (url.format(value), value)

##split var_synmoous (get only clinvar accession and OMIM)
tmp['ClinVar']=tmp["VAR_SYNONYMS"].str.extract(r'(VCV\d*)')

#columns to drop
allvars= tmp.drop([get_max_str(out.columns),'VAR_SYNONYMS'],axis=1)

##need to remove amp-pd once it's fixed
allvars.rename(columns = {'Gene':'Ensembl_geneID', '#mode':'mode',
                          'NEAREST':'NEAREST_gene',
                          'gnomADg_Popmax_AF_filter':'gnomADg_filter',
                         }, inplace = True)

allvars['gnomADg_filter'] = allvars['gnomADg_filter'].str.replace('gnomADg_Popmax_AF_filter','Failed')

#move sample-level column to the last
all_cols=allvars.columns.to_list()
cols_to_move=['sample_id','genotype(sample,dad,mom)']
other_cols = [x for x in all_cols if x not in cols_to_move]
new_cols=other_cols+cols_to_move
allvars=allvars[new_cols]

#sort dataframe
allvars = allvars.sort_values(['mode','chr:pos:ref:alt','BIOTYPE','highest_impact'], ascending = (True, True, False, True))
allvars[['gnomADg_Popmax_AF','gnomADg_nhomalt']] = allvars[['gnomADg_Popmax_AF','gnomADg_nhomalt']].apply(pd.to_numeric)

for fam, dat in allvars.groupby('family_id'):
    h=dat.loc[dat['mode']!='comphet'].sort_values('highest_impact')
    com=dat.loc[dat['mode']=='comphet'].sort_values(['gene','highest_impact'])
    f=pd.concat([h,com],axis=0)
    # create directory
    os.makedirs(f'{snakemake.params.prefix}', exist_ok=True)
    f.to_csv(f'{snakemake.params.prefix}/{fam}_segregation_variants.csv',index=False)
