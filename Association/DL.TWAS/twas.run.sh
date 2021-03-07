#by Delin Li delin.bio@gmail.com
#20201117

while getopts p:f: flag
do
    case "${flag}" in
        f) file=${OPTARG};;
        p) prefix=${OPTARG};;
    esac
done

if [ ! -z ${file} ] && [ ! -z ${prefix} ]
then
echo "Phenotype file: $file; TWAS Type: $twas; MAF for SNPs filtering $maf; genotype vcf file ${vcf}; output files prefix $prefix."
else 
    echo "####Please give parameters:"
    echo "####sh  twas.run.sh -f phenotype file with path -p outfiles prefix"
    exit
fi

if [ ! -e $file ]
then
    echo "####$file not exist, please check the path"
    exit
fi

###The folder
mkdir ${prefix}
echo "Start at:".$("date") >${prefix}/${prefix}.log.txt

######samples in the phenotype
awk '{print $1}' $file >>${prefix}/${file}.sample.temp
cd ${prefix}
mkdir tables

source /datalus/guyongzhe/RNASeq_delin/Kmer/kmer_ind/plink/bashrc

#:<<'END'
#####For R 
Rscript  ~/opt/DL.TWAS/gemma.prepare.R  ${prefix} $file

Tissues=(Leaf Seed)
for tissue in ${Tissues[*]}
do
echo $tissue
mkdir $tissue
cd $tissue
while IFS="	" read -r num trait
do
gemma -g ../tables/${tissue}.Quan5.Geno.${prefix}.txt  -p  ../tables/${tissue}.Quan5.Pheno.${prefix}.txt   -n $num -gk 1 -o ${tissue}.${trait}.Quan5
gemma -g ../tables/${tissue}.Quan5.Geno.${prefix}.txt  -p  ../tables/${tissue}.Quan5.Pheno.${prefix}.txt    -k output/${tissue}.${trait}.Quan5.cXX.txt -o ${prefix}.${tissue}.${trait}.Quan5  -eigen -n $num
gemma -g ../tables/${tissue}.Quan5.Geno.${prefix}.txt -a ../tables/${tissue}.Quan5.Map.${prefix}.txt  -k  output/${tissue}.${trait}.Quan5.cXX.txt -n $num -o ${prefix}.${tissue}.${trait}  -lmm 1  -p  ../tables/${tissue}.Quan5.Pheno.${prefix}.txt
done<../tables/Triats.txt
cd ../
done
#END

echo "GEMMA GWAS done:".$("date") >>${prefix}.log.txt
###GAPIT & FarmCPU
Rscript ~/opt/DL.TWAS/gapit.R ${prefix} $file

###Manhattan Plots
mkdir AssocTables
cd AssocTables
find ../ -name "*assoc.txt" -exec ln -s {} ./ \;
find ../ -name "*GWAS.Results.csv" -exec ln -s {} ./ \;

Rscript ~/opt/DL.TWAS/manhattan.R  ${prefix}
echo "Manhattan Plots were done:".$("date") >>${prefix}.log.txt

echo "Success at:".$("date") >>${prefix}.log.txt
