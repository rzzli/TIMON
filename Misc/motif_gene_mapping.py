import pandas as pd
import numpy as np
import mygene  
import pyorthomap
mg = mygene.MyGeneInfo()
ann_path="/home/zhl022/daima/projects/Motif_coocc_Mar2022/motif_lib/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv"
ann_df=pd.read_csv(ann_path,sep='\t')   # df containing different cutoff
geneout=mg.querymany(ann_df["HGNC"], scopes='hgnc',species="human",as_dataframe=True,fields=['symbol',"ensembl.gene"])
where_gene_col = list(geneout.columns).index("ensembl.gene")
for i in range(geneout.shape[0]):
    #print(i, geneout.iloc[i,:]["ensembl.gene"])
    if pd.isna(geneout.iloc[i,where_gene_col]):
        geneout.iloc[i,where_gene_col] = geneout.iloc[i,:]["ensembl"][0]['gene']
geneout.index=geneout.index.astype("int64")
df_meta1=ann_df.merge(geneout,left_on="HGNC",right_index=True,how="inner")
mouse_convert_df= pyorthomap.findOrthologsHsMm(from_filters = 'link_ensembl_gene_id', from_values = list(geneout["ensembl.gene"])).map()
mouse_convert_df=mouse_convert_df.groupby('human_ensembl_gene_id').first().reset_index()
df_merged=df_meta1.merge(mouse_convert_df,left_on="ensembl.gene",right_on="human_ensembl_gene_id")
df_merged["Model_Name"]=df_merged.apply(lambda x: x[0].split('_')[0],axis=1) # pd series containing tf names 401
df=df_merged[["Model","Model_Name","symbol","ensembl.gene","external_gene_name","mouse_ensembl_gene_id","hgnc_symbol"]]
df.columns=["Model","Model_Name","human_gene_symbol_mygene","human_ensembl_gene_id","mouse_gene_name_Biomart","mouse_ensembl_gene_id_Biomart","hgnc_symbol_Biomart"]
df.to_csv("../motif_lib/gene_mapping.csv",header=True,index=False)