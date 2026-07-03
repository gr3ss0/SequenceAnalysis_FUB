# Final Project
## Discussion of assignment
- we move from viral world to bacteria - data will be larger
- HiFi data is the best thing you can use (accuracy), otherwise LONG reads for complex regions (repetiotions) with high error rate and SHORT reads for polishing
- here is the idea LONG always and SHORT optional add-on
- **input**: FASTQ. already basecalled and demultiplexed
- annotation = prediction of protein coding regions. Optionaly prediciton for other particles (RNA, long/short)
- comparison = unlike to viruses (WGS alignment) , here would be lot missing (uncont), big variance, ...
  -  **CORE GENOME** := set of genes in all genomes (really 100%) vs **PANGENOME**:= union of all samples' genes
-  antibiotic resistance output: one file - one tsv OR multisheet excel (sheet per sample)
-  protein phylogenetics:
  - provided aminoacid seq
  - collect otrhologs from sequence
  - build a tree
  - deal with naming, ?extract from the header
- provide default config and **default sample sheet**
- balance simplicity of workflow with useability

### Project plan
- screen assignment 4 and 5 and think what makes sense to include
- if we say smth stupid because of the lack in reasearch we will be noted however points deducted
- if there are some hidden problems, we will be advised
- every tool named has to be cited. if there's no paper cite GitHub repo
  
### Progress report
- simple message in webex
### Final report
- pdf, 10 pages, depending on number and size of the Figures.
- written like an application note
- What are the steps, Why are the step, Most important parts of Config
- Write it like a Guide throuh wirkflow
- Demonstration, Description,
- For Demonstration use provided test data from Server
- get some numbers and beautiful plots from Test data. Don't do biological interpretation, just describe that this is less continuous than that. Just mention the results (so many of such genes) but don't do claims of early mutations ...
- Do not take any assumptions from Dataset (names, etc). Just .FASTQ and bacterial matter.
- 

## HPC Advice
- let snakemkae manage
  - important to mention resources and threads otherwise it will use default settings
  - snakemake can scan architecture and ask for appropriate nodes
- container tag
- container image and use all conda envs inside 
