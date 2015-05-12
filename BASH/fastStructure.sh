#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# admixture.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: Serial fastStructure analyses\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -i <plink out file> -o <output path> -s <output name> [options]\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis script runs fastStructure from k to K for n CV iterations.\e[0m"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-k <int>\tStarting K [default = 1]\e[0m"
  echo -e "\e[92m\t-K <int>\tFinishing K [default = 6]\e[0m"
  echo -e "\e[92m\t-n <int>\tNumber of CV iterations per K [default = 0]\e[0m"
  echo -e "\n"

}

IN=
OUT=
k=1
K=6
N=0

while getopts ":i:k:K:o:n:" opt; do
  case $opt in
    i) IN=${OPTARG};;
    k) k=${OPTARG};;
    K) K=${OPTARG};;
    o) OUT=${OPTARG};;
    n) N=${OPTARG};;
  esac
done

if [[ -z ${IN} ]] | [[ -z ${OUT} ]] 
then
  usage
  exit 1
fi


for((I=${k};I<=${K};I++));
do 
  echo -e "\e[91mRunning fastStructure for K = \e[94m${I}\n"
  structure.py -K ${I} --cv=${N} --input=${IN} --output=${OUT} 
done

