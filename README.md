# GO-term-enrichment

## Objective
To get enrichment of gene expression clusters using GO-terms or pathways

## Pathway or GO term formatting
Files downloaded from PMN (Plant Metabolic Network) or GO (gene ontology) need to be formatted correctly for the next steps.

### Pathway files
1. Obtain pathway file downloaded from PMN for your species: https://pmn.plantcyc.org/ containing the pathway and genes annotated to that pathway
2. Get the file into pathway:gene format:

        python parse_plantcyc_file_getpath-gene.py <pathway file from PMN> <index where genes are> <index where pathway name is> <index where path ID is>
        
        example:
        
        python parse_plantcyc_file_getpath-gene.py All_instances_of_Pathways_in_Zea_mays_mays.txt 3 0 5

OPTIONAL: If needed, convert gene IDs from pathway file to IDs from your expression data
NEEDED: a BLAST recipricol best match file

        python covert_geneIDs_pathfile.py <BLAST recipricol best match> <pathway file>

### GO term files
1. obtain go.obo file from http://geneontology.org/docs/download-ontology/
2. parse to get GOID: function

        python parse_GO_obo-ID-func.py <go.obo file>

3. for GO term: gene file, need gene associations

### Pfam annotation
1. obtain pfamA.tsv file from http://ftp.ebi.ac.uk/pub/databases/Pfam/Pfam-N/

### gene association
To associate genes, get gene association file from phytozome- this has gene and their pfam IDs and GO IDs
1. https://phytozome-next.jgi.doe.gov/ 
2. sign in and click on species of interest.
3. go to standard data files, select .annotation_info.txt file and download
4. use annotation file, pfamA.tsv file, and go.obo.v1.2_parsed.txt to make table with descriptions:

        python parse_phytozome_ann.py -ann_file <.annotation_info.txt> -pfam_file pfamA.tsv -go_file go.obo.v1.2_parsed.txt -pfam_ind <index with pfam IDs> -go_ind <index with goIDs> -split_by <how pfamIDs and goIDs are delimited- usually a ,>
        
        result: .annotation_info.txt.parsed.txt file with pfam and Go descriptions

## Cluster enrichment
NOTE: With pathway, GO, or cluster file, need to get rid of certain characters “, ‘, / for enrichment, specificially Test_Fisher.py to work
1. get enrichment table with a cluster file

   Need: a cluster tab-delimited file which lists gene:cluster
    
   Need: a GO term or pathway file which contains GO term: gene or pathway: gene
   
        options:
        -cl <cluster file>
        -go <pathway or GO term file>
        -genenum <integer 1 or 2>*
           * Enrichment needs a negative set to compare against. 
             -genenum 1 compares against the total number of genes in the cluster file. 
             -genenum 2 compares against the total number of genes in the pathway/GO term file.
             if there is a limited number of clusters in the cluster file (not all genes from data set) use -genenum 2
        
        example:
        python cluster_enrichment_final.py -genenum 1 -cl Maize_RPKM_nogenelen.txt_PCC.txt_clusters_0.718.txt -go All_instances_of_Pathways_in_Zea_mays_mays.txt.parsed.txt_newID.txt
        
    Output: table for enrichment file: tableforEnrichment_clusterfilename
    
2. do fisher's exact test to get p-value and/or q-value
  
  Use enrichment table to find significant clusters:
        options:
        0 = p-value only
        1 = q-value (multiple testing corrected) and p-value
        
        Notes:
        qvalue.R script must be in same folder as the Test_Fisher.py script
        NO "" in your tableforEnrichment file
  
       python Test_Fisher.py tableforEnrichment_clusterfilename.txt 1
       
  Output: tableforEnrichment_clusterfilename.txt.pqvalue
       
 OPTIONAL:
 
 3. Get only the significant under (-) or over (+) represented clusters
 
        python parse_enrichment_get_sig.py <.pqvalue file>
        

 4. Merge GO term description into results file 
 
        python merge_description.py -key [GO term key] -table [output from step 3] 
 
 5. Get only significant over (+) represented clusters for specific pathway(s)
        options:
        -dir <directory with .pqvalue files>
        -split <delimiter between cluster and pathway, usually "|">
        -path <list of pathways>
        
        python parse_enrich_get_sig_clust_for_path.py -dir ./ -split "|" -path pathA,pathB,pathC
 
 6. Obtain a binary matrix of significant clusters where genes are the row names and columns are the cluster. 1 represents the gene is present in the cluster, 0 represents the gene is absent in the cluster.
    
     options:
       
      REQUIRED:
        
      -cl = file with enrichment for significant clusters (format is .fisher.pqvalue)... if you want all clusters use -cl ALL
        
      -dir = directory with cluster files where file contains: gene \t cluster
        
      -path = file with pathway \t gene
      
      OPTIONAL:
        
      -genes = list of genes you want to extract. This option only gets a matrix that contains clusters with this list of genes
        
      -pval = p-value cutoff for cluster significance
        
      -qval = q-value cutoff for cluster significance
      
      OUTPUT:
        
      binary matrix: filename_binmatrix.txt
        
      example:
        
        python get_binmatrix_for_genes_in_sigclust.py -cl ALL -dir cluster_dir/ -path All_instances_of_Pathways_in_Zea_mays_mays.txt.parsed.txt_newID.txt 
 
## Other functions

### Get percent overlap

1. Calculate percent overlap of your co-expression clusters/modules.

        INPUT:
        
        -bin binary matrix with clusters you are interested in (from get_binmatrix_for_genes_in_sigclust.py)
        
        -mr mr output with all clusters (file is .modules.txt from MR scripts, where columns are: cluster_name | cohesive_score | genes_in_module)
        
        OUTPUT:
        
        dataframe of your cluster overlap percent (_clusteroverlap.txt)
        
        percentiles of all cluster overlap (_allcluster_percentiles.txt) - you can use this file to calculate significance of your cluster overlap. Any percentage above the cluster percentile cutoff is significant for that percentile.
        
        Example:
        
          python calc_percent_cluster_overlap.py -bin ALL_binmatrix.txt -mr Brachy_expressionmat-norm_2017nph_mod.txt_nodup_avg_MR-SP_025.modules.txt 
          
### Get feature log-ratio

1. Given an enrichment file, get percentage of GO/pathway/etc. for a given class (like Sm genes or PM genes) and the log ratio

        INPUT:
        
        enrichment table (tableforEnrichment.pqvalue) with class in second column
        
        OUTPUT:
        
        "_percent_logratio.txt" file
       
