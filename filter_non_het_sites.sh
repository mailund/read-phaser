#!/bin/bash

het_sites=$(mktemp ${TMPDIR}/hetsites.XXXXXX)

{
    cat $1 | awk '{print $1"_"$2}' | sort -k1,1 > ${het_sites}
}

awk '{print $1"_"$2,$3,$4,$5,$6,$7}' | sort -k1,1 | join - ${het_sites} | sed 's/_/ /' | awk '{print $1"_"$3,$2,$4,$5,$6,$7}' | sort -k1,1 | join - ${het_sites} | sed 's/_/ /'| awk -v OFS="\t" '{print $1,$3,$2,$4,$5,$6,$7}' | sort -k1,1 -k2,2n

