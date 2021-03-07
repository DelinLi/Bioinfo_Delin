#by Delin Li delin.bio@gmail.com
#20201105

while getopts m:p:a:f:g: flag
do
    case "${flag}" in
        a) gwas=${OPTARG};;
        m) maf=${OPTARG};;
        f) file=${OPTARG};;
        p) prefix=${OPTARG};;
        g) vcf=${OPTARG};;
    esac
done
#defualt vcf file [huge]
if [ -z ${vcf} ]
then
vcf="/datalus/guyongzhe/RNASeq_delin/ReSeq/Imp/Soybean.Imputed.merged.Chrs.withID.vcf"
fi

if [ ! -z ${file} ] &&  [ ! -z ${gwas} ] && [ ! -z ${maf} ] && [ ! -z ${prefix} ] && [ ! -z ${vcf} ]
then
echo "Phenotype file: $file; GWAS Type: $gwas; MAF for SNPs filtering $maf; genotype vcf file ${vcf}; output files prefix $prefix."
else 
    echo "####Please give parameters:"
    echo "####sh  gwas.run.sh -a GWAS TYPE (CMLM, FarmCPU, GEMMA or ALL) -m MAF cutoff -f phenotype file with path -p outfiles prefix -g the vcf files (uncopressed)"
    exit
fi

if [ ! -e $file ]
then
    echo "####$file not exist, please check the path"
    exit
fi
if [ ! -e $vcf ]
then
    echo "####$vcf not exist, please check the path"
    exit
fi

#:<<'END'
###The folder
mkdir ${prefix}
echo "Start at:".$("date") >${prefix}/${prefix}.log.txt

######samples in the phenotype
awk '{print $1}' $file >>${prefix}/${file}.sample.temp
vcftools --vcf  $vcf  --keep ${prefix}/${file}.sample.temp --maf $maf --out ${prefix}/${prefix} --recode 
echo "VCF extraction done:".$("date") >>${prefix}/${prefix}.log.txt
cd ${prefix} 
~/opt/tassel-5-standalone/run_pipeline.pl -Xms10g  -Xmx60g -fork1 -vcf ${prefix}.recode.vcf  -export ${prefix} -exportType Hapmap -runfork1
echo "Hapmap convertion was done:".$("date") >>${prefix}.log.txt
# the pca
echo "PCA done:".$("date") >>${prefix}.log.txt
plink --vcf ${prefix}.recode.vcf --double-id  --pca 
###Gemma
plink   --vcf ${prefix}.recode.vcf  --out ${prefix}  --double-id
#####For R 
source /datalus/guyongzhe/RNASeq_delin/Kmer/kmer_ind/plink/bashrc
Rscript  ~/opt/DL.GWAS/Fam.R $file ${prefix}
gemma -bfile ${prefix} -gk 1 -o ${prefix} -n 1
gemma -bfile ${prefix} -k output/${prefix}.cXX.txt -lmm 1 -n 1 -o ${prefix}
#END

echo "GEMMA GWAS done:".$("date") >>${prefix}.log.txt
###GAPIT & FarmCPU
Rscript ~/opt/DL.GWAS/gapit.R ${prefix} $file

###Manhattan Plots
Rscript ~/opt/DL.GWAS/manhattan.R  ${prefix}
echo "Manhattan Plots were done:".$("date") >>${prefix}.log.txt

echo "Success at:".$("date") >>${prefix}.log.txt
