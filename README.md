# ITS_primer_trim

Input: Either forward (R1) or reverse (R2) fastq file

Output: Forward (R1) fastq with ITS2R primer sequence trimmed from input forward fastq file AND a log file that records the position of the trimmed primer sequence (when input is forward fastq (R1) file)
OR
Reverse (R2) fastq with ITS1F primer sequence trimmed from input reverse fastq file AND a log file that records the position of the trimmed primer sequence (when input is reverse fastq (R2) file)

## To trim reverse primers from R1: 
python remove_primers.py -r remove_rev_primer_from_R1 -i ${INPUT_R1_FASTQ_PATH} -o ${OUTPUT_R1_FASTQ_PATH}
(e.g. python remove_primers.py -r remove_rev_primer_from_R1 -i ./test_data/Sub10003.V1.sputum.redo_R1.fastq -o ./test_data/no_primer_Sub10003.V1.sputum.redo_R1.fastq)

## To trim forward primers from R2:
python remove_primers.py -r remove_fwd_primer_from_R2 -i ${INPUT_R2_FASTQ_PATH} -o ${OUTPUT_R2_FASTQ_PATH}
(e.g. python remove_primers.py -r remove_fwd_primer_from_R2 -i ./test_data/Sub10003.V1.sputum.redo_R2.fastq -o ./test_data/no_primer_Sub10003.V1.sputum.redo_R2.fastq)
