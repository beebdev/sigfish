#!/bin/bash

set -e

HARU_VENV=~/haru/
REF=test/rnasequin_sequences_2.4.fa
BLOW5=test/sequin_reads.blow5
THREADS=8
REF_PAF=test/rna.minimap2.paf
MY_PAF=test/test.paf

make
./sigfish dtw -g ${REF} -s ${BLOW5} -t ${THREADS} --rna --full-ref -q 500 --from-end > ${MY_PAF}

source ${HARU_VENV}/bin/activate
uncalled pafstats -r ${REF_PAF} ${MY_PAF}
