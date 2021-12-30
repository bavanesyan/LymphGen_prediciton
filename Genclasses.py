def get_entrez_id(hugo_df):
    
    """hugo_df - dataframe with Hugo_Symbol column"""
    
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    
    pandas2ri.activate()

    run_get_entrez_id = ro.r(
        
        """
        library(biomaRt)

         function(hugo_df) {
            print('biomaRt is runnig and it can easily fail due to random server error, so be patient')
            
            mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
            
            tmp_ids <- getBM(attributes = c("hgnc_symbol",  'entrezgene_id'),
                     filters = "hgnc_symbol",
                     values = hugo_df$Hugo_Symbol, 
                     mart = mart)
            
            tmp_ids <- as.data.frame(tmp_ids)
            
            return(tmp_ids)
        }
        """)

    get_entrez_id_results = run_get_entrez_id(hugo_df)


    return get_entrez_id_results

def hugo_to_entrez_map(hugo_df):
    """hugo_df - dataframe with Hugo_symbol column
       Returns dataframe with all original columns + ENTREZ_ID column"""
    
    import numpy as np
    
    hugo_to_entrez_dict = np.load('/uftp/projects/GenClass_Prediction/hugo_to_entrez_dict.npy', allow_pickle=True).item()
    
    hugo_df['Hugo_Symbol'] = hugo_df['Hugo_Symbol'].map(hugo_to_entrez_dict)
    
    hugo_df.dropna(inplace=True)
    
    hugo_df.rename(columns={'Hugo_Symbol': 'ENTREZ_ID'}, inplace=True)
    
    hugo_df.ENTREZ_ID = hugo_df.ENTREZ_ID.astype(int)
            
    return hugo_df

def process_mutation_file(maf):
    
    '''Maf file should be in format of S3Cohort'''
    
    
    print('Generating Mutation file')
    
    # Get entrez_id for genes coded with hugo_symbol
    
    maf = maf.astype(str)
    
    maf = hugo_to_entrez_map(maf)
    
    # Processing
    
    if 'Start_Position' in maf.columns:
        maf = maf[['Tumor_Sample_Barcode', 'ENTREZ_ID', 
                     'Variant_Classification', 'Start_Position']]
    else:
        maf = maf[['Tumor_Sample_Barcode', 'ENTREZ_ID', 
                     'Variant_Classification']]
        maf['Start_Position'] = -1
    
    rename_dict = {
        'Tumor_Sample_Barcode': 'Sample',
        'Variant_Classification': 'Type',
        'Start_Position': 'Location'
    }
    
    maf.rename(columns = rename_dict, inplace=True)
    
    # Reannotation of mutations
    # a. TRUNC: indicating a nonsense mutation in the coding region of the gene.
    # b. MUTATION: indicating a missense or frameshift mutation in the coding region of the gene.
    # c. Synon: indicating a mutation in the 5’UTR of the gene or a synonymous mutation in the
    # coding region within 4kb of the transcription start site.
    # d. L265P: indicating a L265P mutation of the MYD88 gene. (see below)
    
    reannotation = {
        'Missense_Mutation': 'MUTATION',
        'In_Frame_Ins': 'MUTATION',
        'In_Frame_Del': 'MUTATION',
        'Frame_Shift_Ins': 'TRUNC',
        'Nonsense_Mutation': 'TRUNC',
        'Frame_Shift_Del': 'TRUNC', 
        'Nonstop_Mutation': 'TRUNC',
        'Translation_Start_Site': 'TRUNC',
        'Splice_Site': 'TRUNC', 
        'Synonymous': 'Synon',
        'Silent': 'Synon',
        "3'UTR": 'Synon', 
        'UTR5': 'Synon', 
        'UTR3': 'Synon', 
        '3’ UTR': 'Synon', 
        "5'UTR": 'Synon', 
        "3'Flank": 'NA', 
        "5'Flank": 'NA', 
        'Splice_Region': 'NA',
        'IGR': 'NA', 
        'RNA': 'NA',
    }
    
    maf.Type = maf.Type.map(reannotation)
    
    maf.dropna(inplace=True)
    
    maf.sort_values('Sample', inplace=True)
    
    print('Mutation file done')
    
    
    return maf

def get_mutation_genlist(maf_processed):
    
    
    unique_maf_genes = pd.DataFrame(maf_processed.ENTREZ_ID.unique())
    unique_maf_genes.rename(columns={0: 'ENTREZ_ID'}, inplace=True)
    unique_maf_genes.sort_values('ENTREZ_ID', inplace=True)
    
    return unique_maf_genes

def gen_arms_nt(segments):
    
    from tqdm import tqdm
    
    import warnings
    warnings.filterwarnings("ignore")
    
    from bioreactor.cna import process_segments, arm_ploidy_by_mode, normalized_total_cna, calc_ploidy_by_segments
    from gefest.utils import check_array_for_single_value
    
    segments['ploidy'] = calc_ploidy_by_segments(segments)
    
    for sample in tqdm(segments['Sample'].unique(), desc="{Generating segments_arm}"):
            
            sample_segments = segments[segments['Sample'] == sample][[
                'chrom',
                'start',
                'end',
                'n_mark',
                'minor',
                'total',
                'n_h_mark',
                'c_fraction',
                'ploidy',
                'Sample']]

            segments_processed = process_segments(sample_segments)

            arm_ploidy = arm_ploidy_by_mode(segments_processed)
            index = arm_ploidy['chrom'].astype(str).str.cat(arm_ploidy['arm'])
            index.name = 'arm'
            arm_ploidy = arm_ploidy.set_index(index)
            arm_ploidy_nt = normalized_total_cna(
                arm_ploidy.arm_total,
                check_array_for_single_value(segments['ploidy']))
            arm_ploidy_nt.name = sample


            yield arm_ploidy_nt.reindex([
                '1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q', '6p', '6q', '7p',
                '7q', '8p', '8q', '9p', '9q', '10p', '10q', '11p', '11q', '12p', '12q', '13p',
                '13q', '14p', '14q', '15p', '15q', '16p', '16q', '17p', '17q', '18p', '18q',
                '19p', '19q', '20p', '20q', '21p', '21q', '22p', '22q', '23p', '23q'
            ])

def get_segments_arm(segments):
    
    segments_arm = pd.DataFrame(gen_arms_nt(segments)).T
    
    segments_arm.fillna(0, inplace=True)
    
    segments_arm = segments_arm.astype(int)
    
    return segments_arm

def get_arm_file(segments_arm):
    
    from tqdm import tqdm
    
    """Input - arm segments. Returns a dataframe with 3 columns: Sample, Arm, Type, 
       where Type - cna event for specified arm"""
    
    cna_dict = {
        '2': 'AMP',
        '1': 'GAIN',
        '-1': 'HETLOSS',
        '-2': 'HOMDEL',
        'Gain': 'GAIN',
        'ShallowDel': 'HETLOSS',
        'Amplification': 'AMP',
        'Deletion': 'HOMDEL'
    }

    arm_df = {}

    for i in tqdm(segments_arm.columns, desc="{Generating arm_file}"):
        arm_df[i] = segments_arm[i][segments_arm[i] != 0]
    
    arm_df = pd.concat(arm_df)
    arm_df = arm_df.reset_index()
    arm_df.columns = ['Sample', 'Arm', 'Type']
    arm_df = arm_df.astype(str)
    arm_df.Type = arm_df.Type.map(cna_dict)
    arm_df.dropna(inplace=True)
    
    return(arm_df)

def get_cna_gene_file(segments_gene):
    
    from tqdm import tqdm
    
    """Input - gene segments. Genes - rows, Samples - columns. 
       Returns a dataframe with 4 columns: Sample, Hugo_Symbol, Type, Location 
       where Type - cna event for specified gene"""

    
    cna_dict = {
        '2': 'AMP',
        '1': 'GAIN',
        '-1': 'HETLOSS',
        '-2': 'HOMDEL',
        'Gain': 'GAIN',
        'ShallowDel': 'HETLOSS',
        'Amplification': 'AMP',
        'Deletion': 'HOMDEL'
    }
    
    gene_df = {}

    for i in tqdm(segments_gene.columns, desc="{Generating cna_gene_file}"):
        gene_df[i] = segments_gene[i][segments_gene[i] != 0]
    
    gene_df = pd.concat(gene_df)
    gene_df = gene_df.reset_index()
    gene_df.columns = ['Sample', 'Hugo_Symbol', 'Type']
    gene_df = gene_df.astype(str)
    gene_df.Type = gene_df.Type.map(cna_dict)
    gene_df.dropna(inplace=True)
    
    # Obtain entrez_id
    
    gene_df = hugo_to_entrez_map(gene_df)
    
    gene_df = gene_df[['Sample', 'ENTREZ_ID', 'Type']]
    
    gene_df['Location'] = -1
    
    return(gene_df)

def get_unique_cna_genes(cna_gene_file):
    
    """Input - cna gene file obtained with get_cna_gene_file. Returns dataframe with 1 column - ENTREZ_ID"""
    
    unique_cna_genes = pd.DataFrame(cna_gene_file.ENTREZ_ID.unique())
    unique_cna_genes.rename(columns={0: 'ENTREZ_ID'}, inplace=True)
    unique_cna_genes.sort_values('ENTREZ_ID', inplace=True)
    
    return unique_cna_genes

def sample_annotation(maf, cna_gene_file = None, fusion_data = None):
    
    from tqdm import tqdm
    
    """
    maf - maf file from S3Cohort
    cna_gene_file should be generated by get_cna_gene_file()
    fusion_data - raw_fusions from S3Cohort
    """
    
    sample_annotation = pd.DataFrame({
        'Sample.ID': maf.Tumor_Sample_Barcode.sort_values().unique(),
        'Copy.Number': 0,
        'BCL2.transloc': 'NA',
        'BCL6.transloc': 'NA'
    })
    
    if (cna_gene_file is None) & (fusion_data is None):
        print('Sample annotation done')
    else:
        for i in tqdm(sample_annotation['Sample.ID'], desc = '{Generating sample annotation}'):
            if (cna_gene_file is not None):
                if i in cna_gene_file.Sample.values:
                    sample_annotation.loc[sample_annotation['Sample.ID'] == i, 'Copy.Number'] = 1
            if (fusion_data is not None):
                if (i in fusion_data.loc[(fusion_data['Gene B'] == 'BCL6') & (fusion_data['Gene A'] == 'IGH@'), 'Sample'].drop_duplicates().values):
                    sample_annotation.loc[sample_annotation['Sample.ID'] == i, 'BCL6.transloc'] = 1
                if (i in fusion_data.loc[(fusion_data['Gene B'] == 'BCL2') & (fusion_data['Gene A'] == 'IGH@'), 'Sample'].drop_duplicates().values):
                    sample_annotation.loc[sample_annotation['Sample.ID'] == i, 'BCL2.transloc'] = 1
        print('Sample annotation done')

    return sample_annotation

def predict_genclasses(maf, cna_nt = None, segments = None, fusion_data = None):
    
    
    """Performs a prediction of Genclasses from maf, cna_genes, cna_segments and fusion files.
    Maf, cna_nt, segments and fusions should be obtained from S3Cohort. 
    For cna_nt genes should be in raws, samples - in columns
    Returns a dictionary with  Prediction, Compare, model.9, Fullout, Warn.set, modset. 
    Prediction key contains predicted genclasses. Compare contains statistics"""
    
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    
    #Files for prefiction
    
    print('Data processing for prediction')
    
    Mutation_flat = process_mutation_file(maf)
    
    if 'TRUNC' in Mutation_flat.Type.unique():
        trunc = True
    else:
        trunc = False
    
    mut_genelist = get_mutation_genlist(Mutation_flat)
    
    if (cna_nt is None):
        Flatfile = Mutation_flat
        print('Flatfile done')
        CGH = 1
        Copy_flat = 'NULL'
        CGH_genelist = 'NULL'
        Sample_annot = sample_annotation(maf, cna_gene_file=None, fusion_data=fusion_data)
        Sample_annot = Sample_annot.astype(str)
        Sample_annot['Copy.Number'] = Sample_annot['Copy.Number'].astype(int)
    else:   
        Copy_flat = get_cna_gene_file(cna_nt)
        Flatfile = Copy_flat.append(Mutation_flat)
        print('Flatfile done')
        CGH = 0
        CGH_genelist = get_unique_cna_genes(Copy_flat)
        Sample_annot = sample_annotation(maf, cna_gene_file=Copy_flat, fusion_data=fusion_data)
        Sample_annot = Sample_annot.astype(str)
        Sample_annot['Copy.Number'] = Sample_annot['Copy.Number'].astype(int)
    
    Flatfile = Flatfile.astype(str)
    
    if (segments is None) | (CGH == 1):
        Arm_file = 'NULL'
    else:
        segments_arm = get_segments_arm(segments)
        Arm_file = get_arm_file(segments_arm)
        Arm_file = Arm_file.astype(str)
    

    pandas2ri.activate()
    
    run_predict_genclasses = ro.r(
    """
        function(Mutation_flat, Copy_flat, Flatfile, Sample_annot,  mut_genelist, CGH_genelist, Arm_file, CGH, trunc){
    
            source('/uftp/projects/GenClass_Prediction/Model9_original.R')

            options(stringsAsFactors = FALSE)
            Warnset=NULL

            #flags that should be set before running
            CGH.class=CGH 
            Test.set=c("BN2","EZB","MCD","N1","ST2","A53")
            HasL265P=F
            Has.Trunc=trunc
            
            print('Static data')
            
            load('/uftp/projects/GenClass_Prediction/Model_weights.RData', .GlobalEnv)
            
            mut_genelist=mut_genelist[,1]
            
            if (Arm_file == 'NULL'){
                Arm_file = NULL
            }
            
            
            if (Copy_flat == 'NULL'){
                Copy_flat = NULL
            }

            if (CGH.class == 1){
                Arm_file = NULL
                CGH_genelist=unique(Full.Uber.2019[,3])
            }
            else{
                CGH_genelist = CGH_genelist[,1]
            }
            
            
            Sample_annot$BCL2.transloc = as.numeric(Sample_annot$BCL2.transloc)
            Sample_annot$BCL6.transloc = as.numeric(Sample_annot$BCL6.transloc)
            
            print('Starting prediction')
            
            result=Predict9(
                Flatfile,
                Sample_annot,
                Index.Flat,
                Genclass.lab,
                CGH.class=CGH.class,
                Full.Uber.2019,
                Fullmat.2019,
                originalpred=originalpred,
                Has.CGH=Study.samp[,2]==1,
                Arm.file=Arm_file,
                CGH.genelist=CGH_genelist,
                mut.genelist=mut_genelist,
                CalcWilcox=F,
                Test.set=Test.set,
                Has.Trunc=Has.Trunc)
            
            
            return(result)
        
        }
    """
    )
    
    predict = run_predict_genclasses(Mutation_flat, Copy_flat, Flatfile, Sample_annot,  
                                    mut_genelist, CGH_genelist, Arm_file, CGH, trunc)
    
    result = {}
    
    result['Prediction'] = predict[0]
    result['Compare'] = predict[1]
    result['model.9'] = predict[2]
    result['Fullout'] = predict[3]
    result['Warn.set'] = predict[4]
    result['modset'] = predict[5]
    
    return result






    