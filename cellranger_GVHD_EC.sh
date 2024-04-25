cd '/Users/sbhise/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments/Ping1'
mkdir GVHD_EC_expt  ##CHANGE THIS
cd GVHD_EC_expt  ##CHANGE THIS
mkdir rmd res data figs cds
cd ..

ssh sbhise@rhino


/Volumes/fh/scratch/delete90/furlan_s/demux/230530_VH00738_139_AAC3LLJHV/AAC3LLJHV/outs/fastq_path/AAC3LLJHV

export fq=/fh/scratch/delete90/furlan_s/demux/230530_VH00738_139_AAC3LLJHV/AAC3LLJHV/outs/fastq_path/AAC3LLJHV    #CHANGE EVERY TIME

export id=GVHD_EC_expt  ##CHANGE THIS
export wd="/fh/scratch/delete90/furlan_s/${id}"
mkdir $wd
cd $wd
#ls -alh $fq #optionally check fastq folder                             
samps=(GVHD_EC_WT GVHD_EC_Ifngko GVHD_EC_IL6R noGVHD_EC_TCD) #edit this



#change these depending on genome and protein usage
export transcriptome=/shared/biodata/ngs/Reference/10X/refdata-cellranger-mm10-3.0.0        # check transcriptome mouse/human
protein_run=false
export featurefile=/fh/fast/furlan_s/grp/refs/totalseq/TSC_Ping/Tsc_features.csv
export PROTtag="_CSP"

#change these depending on whether vdj was performed
vdj_run=false
export VDJtag="_VDJT"
export vdjfile="/shared/biodata/ngs/Reference/10X/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0"
export vdj_type="VDJ-T"

#change this depending on the sample tag and desired merge name for souporcell
export GEXtag="_GEX"
export mergename="merge"


#rm -R !("archive") #optionally clean folder

#setup environment
export PRIORPATH=$PATH
ml SAMtools/1.14-GCC-11.2.0
ml Python/3.9.6-GCCcore-11.2.0
ml Singularity/3.5.3
export PATH=/home/sbhise/software/cellranger-7.0.0:$PATH   ##CHANGE THIS

#####CELLRANGER-MULTI######

for samp in ${samps[@]}; do
if [ "$protein_run" = true ]; then
  feature_ref_line1="[feature]"
  feature_ref_line2="ref,${featurefile}"
  samp_CSP="${samp}${PROTtag}"
  antibody_capture_line="${samp_CSP},${fq},any,${samp_CSP},Antibody Capture,"
else
  feature_ref_line1=""
  feature_ref_line2=""
  antibody_capture_line=""
fi
if [ "$vdj_run" = true ]; then
  vdj_line1="[vdj]"
  vdj_line2="ref,${vdjfile}"
  samp_VDJ="${samp}${VDJtag}"
  vdj_capture_line="${samp_VDJ},${fq},any,${samp_VDJ},${vdj_type},"
else
  vdj_line1="" 
  vdj_line2=""
  vdj_capture_line=""
fi
samp_GEX="${samp}${GEXtag}"
csv="${wd}/samp_${samp}.csv"
export samp
cat > $csv << EOL
[gene-expression]
ref,$transcriptome
expect-cells,10000
include-introns

$feature_ref_line1
$feature_ref_line2
$vdj_line1
$vdj_line2

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
$samp_GEX,$fq,any,$samp_GEX,Gene Expression,
$antibody_capture_line
$vdj_capture_line
EOL
export csv=$csv
sed -i '/^$/d' $csv
sbatch -n 1 -c 24 -p campus-new  --mem-per-cpu=21000MB --wrap='cellranger multi --id=$samp \
                   --csv=$csv --localcores=24 --localmem=480'
done

squeue -u sbhise
hitparade 


####VELOCYTO#####
#ALLCELLRANGERRUNS=$(sbatch --wrap='echo hello')
####DOWNLOAD VELOCYTO#####
ml Python/3.8.6-GCCcore-10.2.0
#ml velocyto.R/0.6-foss-2019b-R-4.0.2
source $HOME/virtualenvs/velocyto/bin/activate
ml SAMtools/1.14-GCC-11.2.0
#ml Python/3.9.6-GCCcore-11.2.0
ml Singularity/3.5.3
#pip install Cython
#pip install velocyto

export id=GVHD_EC_expt  ##CHANGE THIS
export wd="/fh/scratch/delete90/furlan_s/${id}"

mkdir $wd/velocyto
cd $wd/velocyto
samps=(GVHD_EC_WT GVHD_EC_Ifngko GVHD_EC_IL6R noGVHD_EC_TCD)
#ls -alh $fq #optionally check fastq folder

for sample in ${samps[@]}; do
export sample=$sample
mkdir ${sample}
sbatch -n 1 -c 22 -p campus-new  --mem-per-cpu=21000MB --wrap='. /home/sbhise/virtualenvs/velocyto/bin/activate && \
velocyto run -vv -b ../$sample/outs/per_sample_outs/$sample/count/sample_feature_bc_matrix/barcodes.tsv.gz -m /fh/fast/furlan_s/grp/refs/GRCm38/gencode_M24.gtf -o $sample \
-@ 22 ../$sample/outs/per_sample_outs/$sample/count/sample_alignments.bam /shared/biodata/ngs/Reference/10X/refdata-gex-mm10-2020-A/genes/genes.gtf'
done

squeue -u sbhise
hitparade

#####MERGEBAMS and SOUPORCELL
#mergebams
#install mergebams https://github.com/furlan-lab/mergebams (i installed in my home directory in subfolder called develop)
#mkdir $wd/merge
#cd $wd/merge


ml SAMtools
ml Python/3.9.6-GCCcore-11.2.0
ml Singularity/3.5.3
cd "/fh/scratch/delete90/furlan_s/${id}"
export id=Ping1
export outdir="/fh/scratch/delete90/furlan_s/${id}/merge"
mkdir $outdir
cd $outdir
export bam1=/fh/scratch/delete90/furlan_s/Ping1/PZ_1/outs/per_sample_outs/PZ_1/count/sample_alignments.bam
export bam2=/fh/scratch/delete90/furlan_s/Ping1/PZ_2/outs/per_sample_outs/PZ_2/count/sample_alignments.bam
export bam3=/fh/scratch/delete90/furlan_s/Ping1/PZ_3/outs/per_sample_outs/PZ_3/count/sample_alignments.bam
export bc1=/fh/scratch/delete90/furlan_s/Ping1/PZ_1/outs/per_sample_outs/PZ_1/count/sample_feature_bc_matrix/barcodes.tsv.gz
export bc2=/fh/scratch/delete90/furlan_s/Ping1/PZ_2/outs/per_sample_outs/PZ_2/count/sample_feature_bc_matrix/barcodes.tsv.gz
export bc3=/fh/scratch/delete90/furlan_s/Ping1/PZ_3/outs/per_sample_outs/PZ_3/count/sample_feature_bc_matrix/barcodes.tsv.gz
export labels=PZ_1,PZ_2,PZ_3
sbatch -n 1 -c 24 -p campus-new  --mem-per-cpu=21000MB --wrap='~/.local/bin/mergeBams -i $bam1,$bam2,$bam3 -l $labels -b $bc1,$bc2,$bc3 -o $outdir'
#souporcell
#install souporcell https://github.com/wheaton5/souporcell
RUNNUM=$(sbatch -n 1 -c 24 -p campus-new -M gizmo --mem-per-cpu=16000MB --wrap='samtools sort -@ 24 out.bam -o out.sorted.bam')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN="${LASTRUNTEMP//' on cluster gizmo'}"
RUNNUM=$(sbatch -n 1 -c 12 -p campus-new -M gizmo --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='samtools index -@ 12 out.sorted.bam')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN="${LASTRUNTEMP//' on cluster gizmo'}"
export REF=/fh/fast/furlan_s/grp/refs/GRCm38/refdata-gex-mm10-2020-A/fasta/genome.fa
#export VCF=/fh/fast/bleakley_m/refs/filtered_2p_1kgenomes_unchr.vcf
#you will need to change the lines below (next line and lines 5, 18, and 18 below this comment) to reflect your souporcell install
export sif=/home/sbhise/sifs/souporcell.sif
export SINGULARITY_BINDPATH="/fh/scratch,/fh/fast,/shared"
export OUTDIR=$outdir/souporcell_1
export K=1
mkdir -p $OUTDIR
RUNNUM=$(sbatch -n 1 -c 35 -p campus-new -M gizmo  --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='singularity exec $sif souporcell_pipeline.py \
                  -i out_bam.sorted.bam -b out_barcodes.tsv.gz -f $REF  --skip_remap True \
                  -t 35 -o $OUTDIR -k $K')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN="${LASTRUNTEMP//' on cluster gizmo'}"
ks=(1 2 3 4 5 6 7 8)
for k in ${ks[@]}; do
export K=$k
export OUTDIR=${outdir}/souporcell_${K}
mkdir -p $OUTDIR
RUNNUM=$(sbatch -n 1 -c 35 -p campus-new -M gizmo --dependency=afterok:$LASTRUN --mem-per-cpu=16000MB --wrap='singularity exec $sif /opt/souporcell/souporcell/target/release/souporcell -a souporcell_1/alt.mtx -r souporcell_1/ref.mtx -b out_barcodes.tsv.gz --min_alt 2 --min_ref 2 -k $K -t 35 > $OUTDIR/clusters_tmp.tsv 2> $OUTDIR/log.tsv')
LASTRUNTEMP="${RUNNUM//'Submitted batch job '}"
export LASTRUN2="${LASTRUNTEMP//' on cluster gizmo'}"
sbatch -n 1 -c 1 -p campus-new -M gizmo  --dependency=afterok:$LASTRUN2 --mem-per-cpu=16000MB --wrap='singularity exec $sif troublet -a souporcell_1/alt.mtx -r souporcell_1/ref.mtx --clusters $OUTDIR/clusters_tmp.tsv > $OUTDIR/clusters.tsv'
done


#--common_variants $VCF
##### RUN ON LOCAL COMPUTER AFTER MOUNTING SCRATCH FOLDER#######
R 
id<-"Ping1"

library(colorout)
library(openxlsx)
rdir<-file.path("/Volumes/fh/scratch/delete90/furlan_s", id)
wdir<-file.path("/Users/sbhise/Fred Hutchinson Cancer Research Center/Furlan_Lab - General/experiments", id, "data")
l<-list.files(rdir)
dirs<-c(l[grepl("^[P-S]", l)])
souporcell<-FALSE
souporcell<-FALSE
vdj<-TRUE

  for(i in dirs){
    rfiles<-list.files(file.path(rdir, i, "outs/per_sample_outs", i, "count"), full.name=T)
    if(souporcell){
      soups<-rfiles[grep("^souporcell", basename(rfiles))]
      tc<-file.path(soups, "clusters.tsv")
    }else{
      tc<-vector()
    }
    if(vdj){
      vdj_dir<-list.files(rfiles, full.name=T)[grep("^vdj", basename(list.files(rfiles, full.name=T)))]
      files<-list.files(vdj_dir, full.name=T)
      tc<-c(tc, files[grep(".csv$", basename(files))], files[grep(".fasta$", basename(files))], files[grep(".tsv$", basename(files))], files[grep(".json$", basename(files))])
    }
    rfiles<-list.files(file.path(rdir, i, "outs/per_sample_outs"), full.name=T)
    for(sampledir in rfiles){
      files<-list.files(sampledir, full.name=T)
      tc<-c(tc, files[grep(".csv$", basename(files))], files[grep(".html$", basename(files))])
      morefiles<-list.files(file.path(sampledir, "count"), full.name=T)
      tc<-c(tc, morefiles[grep(".csv$", basename(morefiles))], morefiles[grep(".matrix\\.h5$", basename(morefiles))])
      cn<-gsub(rdir, wdir, tc)
      cd<-lapply(cn, function(x) {
        message(paste0("Creating directory: ", x))
        dir.create(dirname(x), recursive=T)
      })
      ad<-lapply(1:length(tc), function(n) {
        message(paste0("Copying file: ", tc[n], " to: ", cn[n]))
        file.copy(tc[n], cn[n])
      })
    }
  }











