#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

# ==============================================================================
# admixture.sh
# ==============================================================================
echo -e "\n\e[1m\e[34m====================================\e[0m"
echo -e "\e[93m\e[1mSeqApiPop\e[0m\e[93m: Serial ADMIXTURE analyses\e[0m"
echo -e "\e[1m\e[34m====================================\e[0m"
usage()
{
  echo -e "\n\e[96mUsage\e[0m"
  echo -e "\e[92m   $0 -i <plink out file> -o <output path> -s <output name> [options]\e[0m"
  echo -e "\n\e[96mDetails\e[0m"
  echo -e "\e[92mThis script runs admixture from k to K for n iterations.\e[0m"
  echo -e "\n\e[96mOPTIONS:\e[0m"
  echo -e "\e[92m\t-k <int>\tStarting K [default = 1]\e[0m"
  echo -e "\e[92m\t-K <int>\tFinishing K [default = 6]\e[0m"
  echo -e "\e[92m\t-n <int>\tNumber of iterations per K [default = 1]\e[0m"
  echo -e "\e[92m\t-p <bool>\tSupervision flag [default = 0]\e[0m"
  echo -e "\e[92m\t\t\tSupervised analysis [-p 1] requires an additional file with the .pop suffix detailing the ancestries of the reference individuals.\e[0m"
  echo -e "\e[92m\t\t\tThe file path and prefix should match the [-i] file, and it will be detected automatically assuming this convention is followed.\e[0m"
  echo -e "\n\e[96mRequires\e[0m"
  echo -e "\e[92mAdmixture\e[0m"
  echo -e "\n"

}

ADMIX="/usr/local/bioinfo/src/ADMIXTURE/admixture_linux-1.23"
IN=
OUT=
k=1
K=6
N=1
P=0
SAMPLE=

while getopts ":i:k:K:o:n:s:p:" opt; do
  case $opt in
    i) IN=${OPTARG};;
    k) k=${OPTARG};;
    K) K=${OPTARG};;
    o) OUT=${OPTARG};;
    n) N=${OPTARG};;
    s) SAMPLE=${OPTARG};;
    p) P=${OPTARG};;
  esac
done

if [[ -z ${IN} ]] | [[ -z ${OUT} ]] | [[ -z ${SAMPLE} ]]
then
  usage
  exit 1
fi


cd ${OUT}
for((I=${k};I<=${K};I++));
do 
  echo -e "\e[91mRunning ADMIXTURE for K = \e[94m${I}\n"
  if [[ $P == 0 ]]; then ${ADMIX}/admixture -B${N} --cv ${IN}.bed ${I} | tee ${SAMPLE}-${I}.out; fi
  if [[ $P == 1 ]]; then ${ADMIX}/admixture --supervised -B${N} --cv ${IN}.bed ${I} | tee ${SAMPLE}-${I}.out; fi
done

