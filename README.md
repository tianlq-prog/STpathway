# STPathway -- A Robust Statistical Approach for Finding Informative Spatially Associated Pathways
Spatial transcriptomics offers deep insights into cellular functional localization and communication by mapping gene expression to spatial locations. Traditional methods focus on selecting spatially variable genes, which often miss the complexity of biological pathways and network dynamics. Here, we introduce a novel framework, STPathway, that shifts the focus towards directly identifying functional pathways associated with spatial variability. This method adapts the Brownian distance covariance test in an innovative manner to explore the heterogeneity of biological functions over space. Unlike most other methods, this statistical testing approach is free of gene selection and parameter selection, allowing for the detection of nonlinear and complex dependencies.

![workflow](https://github.com/tianlq-prog/STpathway/blob/main/figure/frame_loc.pdf)

# Tutorial 

For the step by step tutoral, please refer to the notebook:

https://github.com/tianlq-prog/STpathway/blob/main/PDAC_example.ipynb

In the tutorial we utilize the human pancreatic ductal adenocarcinoma (PDAC) data as example and show the case of analyzing location related pathways and cancer edge realted pathways.

![PDAC data](https://github.com/tianlq-prog/STpathway/blob/main/figure/pdac_data.pdf)

We first use the function 'find_pathway' to find out all potential pathways, then use the function 'brownian_testing' to conduct testing on each pathway expression verus the spatial location. 

```pythonscript
pathway_df = find_pathway(symbol_sel.values, symbol.values)
brownian_testing(pathway_df, loc.values, adata_data, file_name, n_perm = 199)
path_res = pd.read_csv(file_name)
path_res = path_res.sort_values(by='dcor', ascending=False)
path_res = path_res.reset_index(drop=True)
```
![Location_related_table](https://github.com/tianlq-prog/STpathway/blob/main/figure/pdac_loc_table.png)

To further analyze the cancer edge-related pathways, we first identify the edge areas. The edge ratio is then adjusted according to the scale of the location axis.

```pythonscript
# Find edge nodes
edge = []

for i in range(len(info['Region'])):
    xi, yi = info['x'][i], info['y'][i]
    rule = (np.abs(info['x'] - xi) < 1.5) & (np.abs(info['y'] - yi) < 1.5)
    idx = np.where(rule)[0]
    neighbor_type = info['Region'].iloc[idx]
    if len(np.unique(neighbor_type)) == 1:
        edge.append('Inner') # 0
    else:
        edge.append('Edge')

# Add the 'edge' column to the DataFrame
info['Type'] = edge
```
For each pathway, we test the difference among the cancer inner region and edge region. 

```pythonscript
output_csv = 'PDAC_edge_example.csv'

pd.DataFrame(columns=['ID', 'name', 'n_gene', 'pval', 'dcov', 'dcor']).to_csv(output_csv, index=False)

for i in range(1,len(pathway_df)):
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

    data_in = data[cancer_in_idx,:]
    data_out = data[cancer_edge_idx,:]

    rpy2.robjects.numpy2ri.activate()

    ro.r.assign("data_in", data_in.T)
    ro.r.assign("data_out", data_out.T)

    rpy2.robjects.numpy2ri.activate()
    robjects.r(
               '''
                library(energy)
                x <- data_in
                y <- data_out

                set.seed(1)
                res <- dcov.test(x, y, R=199)

                pval <- res$p.value
                dcov <- res$estimates[1]
                dcor <- res$estimates[2]
               '''
               )
    #symbol = robjects.r['symbol']
    p_val = robjects.r['pval'][0]
    dcov = robjects.r['dcov'][0]
    dcor = robjects.r['dcor'][0]

    # Construct a DataFrame for the current iteration
    current_res = pd.DataFrame({'ID': [pathway_df['GOBPID'][i]],
                                'name': [pathway_df['Term'][i]],
                                'n_gene': [gene_number],
                                'pval': [p_val],
                                'dcov': [dcov],
                                'dcor': [dcor]})

    # Append the result of the current iteration to the CSV file
    current_res.to_csv(output_csv, mode='a', header=False, index=False)
```
To visulize the heatmap of each gene related to certain pathway, we draw the following heatmaps.

```pythonscript
fig, axes = plt.subplots(nrows=5, ncols=5, figsize=(15, 15))

axes = axes.flatten()

for i, ax in enumerate(axes):
    values_to_color = adata_data[:,gene_str[i]].X.flatten()
    ax.scatter(loc['x'].values, loc['y'].values, c=values_to_color, cmap='Reds', s=10)
    ax.set_title(f'{gene_str[i]}')
    ax.set_xlabel('X values')
    ax.set_ylabel('Y values')

plt.tight_layout()

plt.show()
```
![Gene heatmap](https://github.com/tianlq-prog/STpathway/blob/main/figure/padc_loc_gene_map.png)









