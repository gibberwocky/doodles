#!/bin/bash
#$ -M david.wragg@toulouse.inra.fr
#$ -m a


IN=

while getopts ":i:" opt; do
  case $opt in
    i) IN=${OPTARG};;
  esac
done

if [[ -z ${IN} ]] 
then
  exit 1
fi

#Chr	Accession.ver
#LG1	NC_007070.3
#LG2	NC_007071.3
#LG3	NC_007072.3
#LG4	NC_007073.3
#LG5	NC_007074.3
#LG6	NC_007075.3
#LG7	NC_007076.3
#LG8	NC_007077.3
#LG9	NC_007078.3
#LG10	NC_007079.3
#LG11	NC_007080.3
#LG12	NC_007081.3
#LG13	NC_007082.3
#LG14	NC_007083.3
#LG15	NC_007084.3
#LG16	NC_007085.3
#MT	NC_001566.1

sed -i 's/NC_007070.3/1/g' ${IN}
sed -i 's/NC_007071.3/2/g' ${IN}
sed -i 's/NC_007072.3/3/g' ${IN}
sed -i 's/NC_007073.3/4/g' ${IN}
sed -i 's/NC_007074.3/5/g' ${IN}
sed -i 's/NC_007075.3/6/g' ${IN}
sed -i 's/NC_007076.3/7/g' ${IN}
sed -i 's/NC_007077.3/8/g' ${IN}
sed -i 's/NC_007078.3/9/g' ${IN}
sed -i 's/NC_007079.3/10/g' ${IN}
sed -i 's/NC_007080.3/11/g' ${IN}
sed -i 's/NC_007081.3/12/g' ${IN}
sed -i 's/NC_007082.3/13/g' ${IN}
sed -i 's/NC_007083.3/14/g' ${IN}
sed -i 's/NC_007084.3/15/g' ${IN}
sed -i 's/NC_007085.3/16/g' ${IN}
sed -i 's/NC_001566.1/17/g' ${IN}



