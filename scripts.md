# This contains all the scripts for analysis of ChIP-Seq, ATAC-Seq and RNA-Seq including read alignment, peak calling, motif discovery, peak annotations and analysis, enrichment calculation, and visualizations

## WT and Nan/+ E13.5 FL RNA-Seq data from Planutis et. al. 2017 reanalysis using Salmon and DESeq2. This analysis was performed by Sree Chinta
### Obtain primary assembly and transcripts.fa.gz from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/
``` bash
grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt

salmon index -t gencode.vM25.transcripts.fa.gz -i genCodeMouseTxtomeDecoy_index --decoys decoys.txt -k 31

salmon quant -i genCodeMouseTxtome_index \
   	-l SF \
 	-r Nan31_2_AGTTCC_L008_R1_001.fastq.gz \
 	-o Nan31_2_AGTTCC_L008_R1_001.subset.salmon \
 	--useVBOpt \
 	--seqBias \
  	--validateMappings on 3 Nan files (Nan24_2, Nan31_1, Nan31_2) and 3 WT files (WT24_4, WT24_5, WT31_4) \
```
### Prepare data for DESeq2
``` bash
Unzip gencode.vM25.transcripts.fa.gz to gencode.vM25.transcripts.fa

grep -o 'ENSMUST[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]' gencode.vM25.transcripts.fa  > ensmust.txt
grep -o 'ENSMUSG[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]’ gencode.vM25.transcripts.fa  > ensmusg.txt
paste -d ',' ensmust.txt ensmusg.txt > ensmustensmusg.csv
```
### Differential gene expression analysis using DESeq2
``` R
library(DESeq2) 
library(tximport)
quantFiles <- c("<comma separated salmon quant.sf files generated from above>")
gene_map <- read.csv("ensmustensmusg.csv")
count_data <- tximport(files = quantFiles, type = “salmon”, tx2gene = gene_map, ignoreTxVersion = TRUE)
genotype = c('WT', 'Nan')
genotype = rep(genotype, each = 3)
genotype = factor(genotype)
sample_tableFinal$genotype = genotype

deseq_dataset = DESeqDataSetFromTximport(txi = count_data, colData = sample_table, design =~ genotype)
deseq_dataset = estimateSizeFactors(deseq_dataset)
normalizationFactors(deseq_dataset)
counts(deseq_dataset, normalized = TRUE)
deseq_dataset = estimateDispersions(deseq_dataset)
deseq_dataset = nbinomWaldTest(deseq_dataset)
results_table = results(deseq_dataset)
filteredResultsTable = results_table[complete.cases(results_table),]
pValFilteredResultsTable = filteredResultsTable[filteredResultsTable$padj < 0.05, ]
pValLogFilteredResultsTable = pValFilteredResultsTable[abs(pValFilteredResultsTable$log2FoldChange) > 1,]
```

## ChIP-Seq and ATAC-Seq analysis 

## Read alignment using Bowtie2 to generate alignment sam file. 

``` bash
bowtie2 --sensitive-local  \
--trim3 10 \
-x /path-to-genome-index \
-U /path-to-input-fastq \
-S /path-to-out-sam
```

## Converting the output sam files for ChIP-Seq and ATAC-Seq to an input normalized bigwig file for visualization in a genome browser such as IGV

### Step1. Convert sam file to sorted bam file using samtools
``` bash
samtools view -bu /path-to-in.sam | samtools sort -T /path-to-temp-dir  -o /path-to-outfile.bam -
```

### Step2. Index the bam file using samtools
``` bash
cd /path-to-outfile-dir
for file in ./*.bam
do
	samtools index $file
done
```

### Step3. Convert bam to input-normalized bigwig using Deeptools bamCompare
``` bash
bamCompare 	--scaleFactorsMethod None \
				--normalizeUsing RPKM \
				--operation first \
				--binSize 50 \
				--effectiveGenomeSize 2652783500 \
				--smoothLength 100 \
				--bamfile1 file1.bam \
				--bamfile2 file1-input.bam \
				--outFileFormat bigwig \
				--outFileName /path-to-outdir/file1.bw
```

### Step4. Merge multiple normalized bigwig files from replicate ChIP-Seq experiments into one bigwig file using ucsc-utils
``` bash
bigWigMerge rep1.bw rep2.bw rep3.bw file.bedGraph

sort -k1,1 -k2,2n file.bedGraph > file.bg

bedGraphToBigWig file.bg /PATH/TO/mm10/chrom.sizes file.bw
```

## Peak calling and annotations with Homer
### Step1. Create Tag Directories from aligned sam files
``` bash
makeTagDirectory /PATH-to-TD/TD /PATH-to-file/file.sam
```
### Step2. Use ChIP and Input tag directories to call differential peaks using specific options. The following was used for the EKLF ChIP-Seq data.
``` bash
getDifferentialPeaksReplicates.pl 	-t  eklf_rep1 eklf_rep2 eklf_rep3 \
								-i input_rep1 input_rep2 input_rep3 \
								-genome mm10 \
								-size 125 \
								-minDist 250 \
								-region \
								-q 0.01 \
								-f 4 \
								> eklf_peaks.txt

pos2bed.pl eklf_peaks.txt > eklf_peaks.bed ### to generate bed files for each peak file.
bed2pos.pl file.bed > file.txt ### to get Homer peak file from a bed file
```
## EKLF 7B2a ChIP-Seq motif  determination using Homer.

``` bash
findMotifsGenome.pl eklf_peaks.txt mm10 OUTPUT-Dir 
	-size given \
	-mask \
	-dumpFasta \
```

## Heatmaps and enrichment profiles for ChIP and ATAC using Deeptools
### Step1. Generate the enrichment values for specific peak regions for ChIP-Seq data e.g. EKLF-ChIP
``` bash
computeMatrix reference-point \
								--referencePoint center \
								--regionsFileName eklf_peaks.bed \
								--scoreFileName eklf_chip_input_norm.bw nan_chip_input_norm.bw \
								--outFileName eklf_peaks.gz \
								-b 2000 \
								-a 2000 \
								--smartLabels \
								--skipZeros \
								--missingDataAsZero \
								--sortUsingSamples 1
```

### Step2. Plot the heatmap
``` bash
plotHeatmap 	-m eklf_peaks.gz \
				-o eklf_peaks.png \
				--colorMap YlOrRd \
				--zMin "auto" \
				--zMax "auto" \
				--dpi 300 \
				--regionsLabel '' \
				--samplesLabel "WT EKLF" "Nan/+ EKLF" \
				--yAxisLabel "norm coverage"
```

### Step3. Plot ChIP or ATAC profile only
``` bash
plotProfile 	-m eklf_peaks.gz \
			-o eklf_peaks_profile.png \ 
			--colors firebrick salmon \
			--dpi 300 \
			--samplesLabel "WT" "Nan/+" \
			--yAxisLabel "norm coverage" \
			--perGroup
```

### For plotting ChIP/ATAC profiles of gene regions, get the bed file from UCSC table browser for the corresponding genes and then use the same methods above with the following modification
    ### For TSS use:
``` bash
computeMatrix reference-point --referencePoint TSS
```
	### For entire gene region use:
``` bash
computeMatrix scale-regions --regionBodyLength 5000
```
	### in each case the -b and -a parameters can be adjusted accordingly

## Calculating enrichment of ChIP-Seq using Homer
``` bash
	# RNAP II S5
annotatePeaks.pl tss mm10 
						-d WT-Ser5-1/ WT-Ser5-2/ Nan-Ser5-1/ Nan-Ser5-2/ \
						-list eklf_target_genes.csv \ > enrichments.txt

	# RNAP II S2
annotatePeaks.pl tss mm10 
						-d WT-Ser2-1/ WT-Ser2-2/ Nan-Ser2-1/ Nan-Ser2-2/ \
						-list eklf_target_genes.csv \ > enrichments.txt
```

## EKLF peak analysis and violin plots using Python
``` python3
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from pylab import savefig
    %matplotlib inline
    
	w_peaks = pd.read_csv('eklf_peaks.txt', sep='\t')

# Getting EKLF target genes - removing duplicate gene rows from peak list:
	def unique(peaks):
   		 for file in peaks:
        		df = pd.read_csv('{s}.txt'.format(s=peaks), sep='\t')
        		df_genes = df[['Gene Name']].drop_duplicates()
        		df_genes.to_csv('{s}_genes.csv'.format(s=peaks), index=False, header=None)
                
# Then get EKLF and Nan target genes by:
	peaks = ['eklf_peaks', 'nan_peaks']
	unique(peaks)

# Get promoter, intron, intergenic peaks from Homer peak file
	w_prom = w_peaks[w_peaks['Annotation'].str.startswith('promoter', na = False)]
	w_intron = w_peaks[w_peaks['Annotation'].str.startswith('intron', na = False)]
	w_exon = w_peaks[w_peaks['Annotation'].str.startswith('exon', na = False)]
	w_intergenic = w_peaks[w_peaks['Annotation'].str.startswith('Intergenic', na = False)]
	w_utr = w_peaks[w_peaks['Annotation'].str.contains('UTR', na = False)]

# Get numbers of each as follows:
	w_prom.shape
# Create dataframe for plotting
	wt_nan_features = pd.DataFrame([['promoter',3725,2744], ['intron',3042,8026], ['intergenic',1556,3641], ['exon',262,325], ['UTR',317,376]],columns=['feature','WT', 'Nan/+'])
# For all enhancer peaks:
	w_intron_intergenic = w_peaks[w_peaks['Annotation'].str.contains('nt', na = False)]

# EKLF targets compared with differential gene expression from RNA-Seq data
	rna_down = pd.read_csv('down_genes.txt') ## down and up genes are just the list of gene names that are down or up-regulated in Nan/+ over WT
	rna_up = pd.read_csv('up_genes.txt')
	down_w_peaks = rna_down.merge(w_peaks)
	up_n_peaks = rna_up.merge(n_peaks)
# Getting genes that are EKLF direct targets (or Nan targets with up_n_peaks object)
	down_w_peaks_uniq = down_w_peaks.loc[:,['Genename']]
	down_w_peaks_uniq.drop_duplicates(inplace=True)
# Similarly for promoter or enhancer peaks
	w_enhancer = w_intron_intergenic
	down_w_enhancer = w_enhancer.merge(rna_down, on='Gene Name')
    
# Violin plots from Chip enrichment
    ax = sns.violinplot(data=x, inner='box', linewidth=2, 
               palette=['color1','color2'], saturation=1)
    ax.set(ylabel='ChIP Enrichment')
    ax=ax.get_figure()
    ax.savefig('fig.png', dpi=300, bbox_inches='tight')

# Enhanced boxplots for Chip enrichment
    ax = sns.boxenplot(data=x, linewidth=2, 
               palette=['color1','color2'], saturation=1, k_depth = 'proportion', outlier_prop = 0.1, showfliers=False)
    ax.set(ylabel='ChIP Enrichment')
    ax=ax.get_figure()
    ax.savefig('fig.png', dpi=300, bbox_inches='tight')
    
# Quadrant plot
    ax = sns.scatterplot(data=x, x="log2 w/n s5", y="log2 w/n s2")
    ax.set(ylabel = "Difference in Gene Body pSer2", xlabel = "Difference in TSS pSer5")
    ax=ax.get_figure()
    ax.savefig('quadrant.png', dpi=300, bbox_inches='tight')
```



