usage: GPCR_Selectivity_Explorer.py [-h] [-l1 LIST1 [LIST1 ...]] [-r1 REFERENCE1] [-l2 LIST2 [LIST2 ...]] [-r2 REFERENCE2] [-sg SEGMENTS [SEGMENTS ...]]
                                    [-c CONSERVATION_CUTOFF] [-sb SUBSTITUTION_MATRIX_METHOD] [-o OUTPUTFOLDER] [--custom] [-r1f REFERENCE1_RESIDUES_TABLE]
                                    [-r1A REFERENCE1_ALIGNMENT] [-r2f REFERENCE2_RESIDUES_TABLE] [-r2A REFERENCE2_ALIGNMENT]


optional arguments:
  -h, --help            show this help message and exit
  -l1 LIST1 [LIST1 ...], --list1 LIST1 [LIST1 ...]
                        Space-separated list for Reference 1 Receptor segments (e.g., adrb1_human adrb2_human adrb3_human)
  -r1 REFERENCE1, --reference1 REFERENCE1
                        Reference 1 Receptor Name in GPCRdb (e.g., adrb2_human)
  -l2 LIST2 [LIST2 ...], --list2 LIST2 [LIST2 ...]
                        Space-separated list for Reference 2 Receptor segments (e.g., drd1_human drd5_human)
  -r2 REFERENCE2, --reference2 REFERENCE2
                        Reference 2 Receptor Name in GPCRdb (e.g., drd1_human)
  -sg SEGMENTS [SEGMENTS ...], --segments SEGMENTS [SEGMENTS ...]
                        Segment information for reference receptor residues from GPCRdb (e.g., TM1 TM2 ECL2 ICL1)
  -c CONSERVATION_CUTOFF, --conservation_cutoff CONSERVATION_CUTOFF
                        Conservation Cutoff
  -sb SUBSTITUTION_MATRIX_METHOD, --substitution_matrix_method SUBSTITUTION_MATRIX_METHOD
                        Substitution matrix method. Options include: BENNER22, BENNER6, BENNER74, BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90, DAYHOFF, FENG,
                        GENETIC, GONNET1992, HOXD70, JOHNSON, JONES, LEVIN, MCLACHLAN, MDM78, NUC.4.4, PAM250, PAM30, PAM70, RAO, RISLER, SCHNEIDER, STR, TRANS
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Output folder
  --custom              use custom alginments and residue tables
  -r1f REFERENCE1_RESIDUES_TABLE, --reference1_residues_table REFERENCE1_RESIDUES_TABLE
                        Reference 1 modifed residues table excel file from GPCRdb
  -r1A REFERENCE1_ALIGNMENT, --reference1_alignment REFERENCE1_ALIGNMENT
                        Reference 1 modifed alignment excel file from GPCRdb
  -r2f REFERENCE2_RESIDUES_TABLE, --reference2_residues_table REFERENCE2_RESIDUES_TABLE
                        Reference 2 modifed residues table excel file from GPCRdb
  -r2A REFERENCE2_ALIGNMENT, --reference2_alignment REFERENCE2_ALIGNMENT
                        Reference 2 modifed alignment excel file from GPCRdb


1- Example to run:


python GPCR_Selectivity_Explorer.py -l1 adrb1_human adrb2_human adrb3_human -r1 adrb2_human -l2 drd1_human drd5_human -r2 drd1_human -sg TM1 TM2 TM3 TM4 TM5 TM6 TM7 ECL2 -c 80 -sb BLOSUM62 -o outputs_default

CUSTOM ALIGNMENTS:

To use custom alignments. First, run GPCRs Selectivity Explorer and edit (remove/add sequences) the output alignment files for List A and List B. If you add new positions (e.g., the full ECL3), you must also update the reference residue table by adding the amino acid number(s) and its corresponding GPCRdb number(s) (e.g., E193 67x10). If no GPCRdb numbers exists for the new position(s), create GPCRdb number(s), ensure each number is unique and is not already assigned to another position. Finally, use the custom alignment and residue tables as inputs for the workflow.


Example to run using custom alignments:

python GPCR_Selectivity_Explorer.py -r1 adrb2_human   -r2 drd1_human  -c 80 -r1f  ../outputs_default/adrb2_human_GPCRdb_residues_table_full_seq.xlsx -r1A ../outputs_default/adrb2_human_GPCRdb_alignment.xlsx  -r2f  ../outputs_default/drd1_human_GPCRdb_residues_table_full_seq.xlsx  -r2A  ../outputs_default/drd1_human_GPCRdb_alignment.xlsx  --custom -o outputs_custom
