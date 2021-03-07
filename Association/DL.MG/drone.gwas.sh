#by Delin Li delin.bio@gmail.com
#20201215

while getopts p:f: flag
do
    case "${flag}" in
        f) file=${OPTARG};;
        p) prefix=${OPTARG};;
    esac
done

if [ ! -z ${file} ] && [ ! -z ${prefix} ]
then
echo "For Drone GWAS."
echo "Phenotype file: $file; output files prefix $prefix."
else 
    echo "####Please give parameters:"
    echo "####sh  drone.run.sh -f phenotype file with path -p outfiles prefix"
    exit
fi

if [ ! -e $file ]
then
    echo "####$file not exist, please check the path"
    exit
fi

#:<<'END'
###The folder
mkdir ${prefix}
echo "Start at:".$("date") >${prefix}/${prefix}.log.txt

######samples in the phenotype
awk '{print $1}' $file >>${prefix}/${file}.sample.temp
cd ${prefix} 
ln -s ../../data/Drone.bed ./${prefix}.bed
ln -s ../../data/Drone.bim ./${prefix}.bim
ln -s ../../data/plink.e* ./
cp ../../data/Drone.fam ./${prefix}.fam
ln -s ../../data/Drone.cxx.txt ./${prefix}.cXX.txt
#####For R 
source /datalus/guyongzhe/RNASeq_delin/Kmer/kmer_ind/plink/bashrc
Rscript  ~/opt/DL.MG/Fam.R $file ${prefix}
while IFS="	" read -r num trait
do
gemma -bfile ${prefix} -k  ${prefix}.cXX.txt -lmm 1 -n $num -o ${prefix}.${trait}
Rscript ~/opt/DL.MG/manhattan.R  ${trait}
done<Triats.txt
#END


echo "Success at:".$("date") >>${prefix}.log.txt
