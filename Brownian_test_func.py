import rpy2.robjects as ro
import rpy2.robjects as robjects
import rpy2.robjects
import rpy2.robjects.numpy2ri
from scipy.sparse import issparse
import numpy as np
import pandas as pd
from tqdm import trange

def find_pathway(symbol, symbol_all):
    """
    
    Find all related pathway
    
    Input: 
    -symbol: the interest gene
    -symbol_all: all gene list before filtering 
    """
    
    rpy2.robjects.numpy2ri.activate()

    ro.r.assign("symbol_sel", symbol)
    ro.r.assign("symbol_all", symbol_all)
    #ro.r.assign("name", name)

    rpy2.robjects.numpy2ri.activate()
    robjects.r(
               '''
                library("org.Hs.eg.db")
                library(GOstats)
                library(qvalue)
                
                if (grepl("ENS", symbol_all[1])) {
                    print("The first element contains 'ENS'.")
                    ens <- symbol_sel
                    ens_all <- symbol_all
                  
                }else{
                    print("The first element does not contain 'ENS'.")

                    info <- AnnotationDbi::select(org.Hs.eg.db, keys=symbol_sel,
                        columns=c("ENSEMBL"), 
                        keytype="SYMBOL")
                    ens <- info$ENSEMBL
                    ens <- ens[!is.na(ens)]

                    info_all <- AnnotationDbi::select(org.Hs.eg.db, keys=symbol_all,
                        columns=c("ENSEMBL"), 
                        keytype="SYMBOL")
                    ens_all <- info_all$ENSEMBL
                    ens_all <- ens_all[!is.na(ens_all)]
                }
                
                
                analyze<- function(idx){
                  sel.entrez<-  ens # selected gene names (15 character)
                  all.entrez<-  ens_all # all gene names in the data matrix under study

                  # if there is more than one matched, then all are taken to the next step

                  sel.entrez<-unlist(mget(substr(sel.entrez,1,15), org.Hs.egENSEMBL2EG,ifnotfound=NA))
                  sel.entrez<-unique(sel.entrez[!is.na(sel.entrez)])
                  all.entrez<-unlist(mget(substr(all.entrez,1,15), org.Hs.egENSEMBL2EG,ifnotfound=NA))
                  all.entrez<-unique(all.entrez[!is.na(all.entrez)])

                  # this is hypergeometric testing, by GOstats

                  params <- new("GOHyperGParams", geneIds=sel.entrez[!is.na(sel.entrez)], universeGeneIds=all.entrez[!is.na(all.entrez)], ontology="BP", pvalueCutoff=1, annotation="org.Hs.eg.db")

                  over = hyperGTest(params)
                  ov<-summary(over)
                  ov<-ov[ov[,6]<=500 & ov[,6]>=12,]

                  ov<-ov[ov[,5]>10,]

                  for(i in 2:4) ov[,i]<-signif(ov[,i], 3)

                  # this is the fill in which genes caused the pathway to be significant           
                  ov<-cbind(ov, idx = rep('',nrow(ov)), name = rep('', nrow(ov)))
                  ov[,8]<-as.vector(ov[,8])
                  ov[,9]<-as.vector(ov[,9])

              for(i in 1:nrow(ov))
              {
                       this.genes<-mget(ov[i,1],org.Hs.egGO2ALLEGS)[[1]]
                       this.genes<-unique(this.genes[this.genes %in% sel.entrez])
                       that<-this.genes[1]
                       if(length(this.genes)>1)
                       {
                               for(j in 2:length(this.genes)) that<-paste(that, ", ", this.genes[j], sep="")
                       }
                       ov[i,8]<-that

                       this.genes<-unlist(mget(this.genes, org.Hs.egSYMBOL))
                       that<-this.genes[1]
                       if(length(this.genes)>1)
                       {
                               for(j in 2:length(this.genes)) that<-paste(that, ", ", this.genes[j], sep="")
                       }
                       ov[i,9]<-that

              }
                #ov <- cbind(ov, lfdr = rep(0,nrow(ov)))
                #ov$lfdr[which(is.na(ov$Pvalue)==FALSE)] <- lfdr(na.omit(ov$Pvalue))
                #ff <- lfdr(c(na.omit(ov$Pvalue)))
                return(ov)

                }
                #info <- AnnotationDbi::select(org.Hs.eg.db, keys=name,
                #                    columns=c("SYMBOL"),
                #                    keytype="ENSEMBL")
                #symbol <- info$SYMBOL
                res <- analyze(1)
               '''
               )
    # Assuming 'res' is the R dataframe
    res_r = robjects.r['res']

    # Convert R dataframe to pandas DataFrame
    res_py = pd.DataFrame(res_r)
    res_py = res_py.T

    # Extract column names from R dataframe
    colnames_r = robjects.r.colnames(res_r)
    colnames_py = list(colnames_r)

    # Set column and row names in the pandas DataFrame
    res_py.columns = colnames_py
    res_py
    pathway_df = res_py

    return(pathway_df)



def brownian_testing(pathway_df, loc, adata_data, output_csv, n_perm = 199):
    
    pd.DataFrame(columns=['ID', 'name', 'n_gene', 'pval', 'dcov', 'dcor']).to_csv(output_csv, index=False)

    
    GO_id = []
    GO_name = []
    pval = []
    est = []
    n_gene = []

    print(f"Now conduct brownian test for {len(pathway_df)} Pathways")
    
    for i in trange(1,len(pathway_df)):
        #if i%100==0:
        #    print(f'-----Round------{i}')

        this_pathway = i

        gene_str = pathway_df['name'][this_pathway]

        gene_number = pathway_df['Count'][this_pathway]

        if gene_number<500 and gene_number>10:

            gene_str = gene_str.strip("'")

            gene_str = gene_str.split(', ')

            gene_str = [item.replace('MICOS10-', '') for item in gene_str]
            gene_str = [item.replace('MMP24-AS1-', '') for item in gene_str]
            gene_str = [gene for gene in gene_str if gene in adata_data.var_names]

            gene_number = len(gene_str)

            # Assuming gene_str is defined and adata_data is your AnnData object
            data_subset = adata_data[:, gene_str].X

            # Check if the data is in a sparse matrix format
            if issparse(data_subset):
                data = data_subset.todense()
            else:
                data = data_subset


            geneName = adata_data[:,gene_str].var_names

            rpy2.robjects.numpy2ri.activate()

            ro.r.assign("data", data)
            ro.r.assign("loc", loc)
            ro.r.assign("n_perm", n_perm)

            rpy2.robjects.numpy2ri.activate()
            robjects.r(
                       '''
                        library(energy)
                        x <- data
                        y <- loc

                        set.seed(1)
                        res <- dcov.test(x, y, R=n_perm)

                        pval <- res$p.value
                        dcov <- res$estimates[1]
                        dcor <- res$estimates[2]
                       '''
                       )
            #symbol = robjects.r['symbol']
            p_val = robjects.r['pval'][0]
            dcov = robjects.r['dcov'][0]
            dcor = robjects.r['dcor'][0]
            
            # GO_id.append(pathway_df['GOBPID'][i])
            # GO_name.append(pathway_df['Term'][i])
            # pval.append(p_val)
            # #est.append(estimate)
            # DCOV.append(dcov)
            # DCOR.append(dcor)
            # #n_gene.append(gene_number)
            # n_gene.append(data_in.shape[1])

            # Construct a DataFrame for the current iteration
            current_res = pd.DataFrame({'ID': [pathway_df['GOBPID'][i]],
                                        'name': [pathway_df['Term'][i]],
                                        'n_gene': [gene_number],
                                        'pval': [p_val],
                                        'dcov': [dcov],
                                        'dcor': [dcor]})


            # Append the result of the current iteration to the CSV file
            current_res.to_csv(output_csv, mode='a', header=False, index=False)
    return(current_res)