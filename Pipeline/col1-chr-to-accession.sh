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

sed -i 's/1\t/NC_007070.3\t/' ${IN}
sed -i 's/2\t/NC_007071.3\t/' ${IN}
sed -i 's/3\t/NC_007072.3\t/' ${IN}
sed -i 's/4\t/NC_007073.3\t/' ${IN}
sed -i 's/5\t/NC_007074.3\t/' ${IN}
sed -i 's/6\t/NC_007075.3\t/' ${IN}
sed -i 's/7\t/NC_007076.3\t/' ${IN}
sed -i 's/8\t/NC_007077.3\t/' ${IN}
sed -i 's/9\t/NC_007078.3\t/' ${IN}
sed -i 's/10\t/NC_007079.3\t/' ${IN}
sed -i 's/11\t/NC_007080.3\t/' ${IN}
sed -i 's/12\t/NC_007081.3\t/' ${IN}
sed -i 's/13\t/NC_007082.3\t/' ${IN}
sed -i 's/14\t/NC_007083.3\t/' ${IN}
sed -i 's/15\t/NC_007084.3\t/' ${IN}
sed -i 's/16\t/NC_007085.3\t/' ${IN}
sed -i 's/17\t/NC_001566.1\t/' ${IN}



