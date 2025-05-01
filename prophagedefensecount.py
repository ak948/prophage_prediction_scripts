import os
import subprocess
from Bio import SeqIO
import csv
import pandas as pd
import statistics
import shutil
import numpy as np
import pylab
import matplotlib.pyplot as plt
from pylab import MaxNLocator
import gzip
from Bio import Entrez
import urllib
import time
import sys

csv.field_size_limit(sys.maxsize)

with open('bacteriaprophagedefensecount.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(["BiosampleID", "Prophage Defense System Count",
                         "Prophage Defense System Subtype ‘Abi2’ Count", "Prophage Defense System Subtype ‘AbiC’ Count",
                         "Prophage Defense System Subtype ‘AbiD’ Count", "Prophage Defense System Subtype ‘AbiE’ Count",
                         "Prophage Defense System Subtype ‘AbiG’ Count", "Prophage Defense System Subtype ‘AbiH’ Count",
                         "Prophage Defense System Subtype ‘AbiI’ Count", "Prophage Defense System Subtype ‘AbiJ’ Count",
                         "Prophage Defense System Subtype ‘AbiL’ Count", "Prophage Defense System Subtype ‘AbiN’ Count",
                         "Prophage Defense System Subtype ‘AbiP2’ Count",
                         "Prophage Defense System Subtype ‘AbiQ’ Count", "Prophage Defense System Subtype ‘AbiR’ Count",
                         "Prophage Defense System Subtype ‘AbiU’ Count",
                         "Prophage Defense System Subtype ‘Avs_II’ Count", "Prophage Defense System Subtype ‘BREX_I’ Count",
                         "Prophage Defense System Subtype ‘BREX_III’ Count",
                         "Prophage Defense System Subtype ‘CAS_Class2-Subtype-II-A’ Count",
                         "Prophage Defense System Subtype ‘CAS_Class2-Subtype-II-C’ Count",
                         "Prophage Defense System Subtype ‘CBASS_I’ Count",
                         "Prophage Defense System Subtype ‘CBASS_II’ Count",
                         "Prophage Defense System Subtype ‘DarTG’ Count", "Prophage Defense System Subtype ‘DdmDE’ Count",
                         "Prophage Defense System Subtype ‘Detocs’ Count", "Prophage Defense System Subtype ‘Dodola’ Count",
                         "Prophage Defense System Subtype ‘DRT_2’ Count", "Prophage Defense System Subtype ‘Eleos’ Count", "Prophage Defense System Subtype ‘Gabija’ Count",
                         "Prophage Defense System Subtype ‘Gao_Iet’ Count", "Prophage Defense System Subtype ‘Gao_Qat’ Count",
                         "Prophage Defense System Subtype ‘Hachiman’ Count", "Prophage Defense System Subtype ‘Hna’ Count",
                         "Prophage Defense System Subtype ‘Kiwa’ Count", "Prophage Defense System Subtype ‘Lamassu-Cap4_nuclease’ Count",
                         "Prophage Defense System Subtype ‘Lamassu-Hypothetical’ Count", "Prophage Defense System Subtype ‘Lamassu-Lipase’ Count", "Prophage Defense System Subtype ‘Lamassu-Mrr’ Count",
                         "Prophage Defense System Subtype ‘Lamassu-Protease’ Count", "Prophage Defense System Subtype ‘MazEF’ Count", "Prophage Defense System Subtype ‘Nhi’ Count",
                         "Prophage Defense System Subtype ‘NLR_like_bNACHT01’ Count", "Prophage Defense System Subtype ‘Olokun’ Count", "Prophage Defense System Subtype ‘PARIS_I’ Count",
                         "Prophage Defense System Subtype ‘PD-Lambda-1’ Count", "Prophage Defense System Subtype ‘PD-Lambda-5’ Count",
                         "Prophage Defense System Subtype ‘PD-T4-1’ Count", "Prophage Defense System Subtype ‘PD-T4-8’ Count",
                         "Prophage Defense System Subtype ‘PD-T4-9’ Count", "Prophage Defense System Subtype ‘PD-T7-2’ Count", "Prophage Defense System Subtype ‘PD-T7-3’ Count",
                         "Prophage Defense System Subtype ‘PD-T7-4’ Count", "Prophage Defense System Subtype ‘PD-T7-5’ Count", "Prophage Defense System Subtype ‘Pycsar’ Count",
                         "Prophage Defense System Subtype ‘Retron_I_C’ Count", "Prophage Defense System Subtype ‘Retron_III’ Count", "Prophage Defense System Subtype ‘Retron_XI’ Count",
                         "Prophage Defense System Subtype ‘RexAB’ Count", "Prophage Defense System Subtype ‘RloC’ Count",
                         "Prophage Defense System Subtype ‘RM_Type_I’ Count", "Prophage Defense System Subtype ‘RM_Type_II’ Count", "Prophage Defense System Subtype ‘RM_Type_IIG’ Count",
                         "Prophage Defense System Subtype ‘RM_Type_III’ Count", "Prophage Defense System Subtype ‘RM_Type_IV’ Count", "Prophage Defense System Subtype ‘RnlAB’ Count",
                         "Prophage Defense System Subtype ‘RosmerTA’ Count", "Prophage Defense System Subtype ‘Rst_HelicaseDUF2290’ Count",
                         "Prophage Defense System Subtype ‘SanaTA’ Count", "Prophage Defense System Subtype ‘SEFIR’ Count",
                         "Prophage Defense System Subtype ‘Septu’ Count", "Prophage Defense System Subtype ‘Shedu’ Count",
                         "Prophage Defense System Subtype ‘ShosTA’ Count", "Prophage Defense System Subtype ‘SoFic’ Count",
                         "Prophage Defense System Subtype ‘SpbK’ Count", "Prophage Defense System Subtype ‘Stk2’ Count",
                         "Prophage Defense System Subtype ‘Thoeris_I’ Count", "Prophage Defense System Subtype ‘Thoeris_II’ Count",
                         "Prophage Defense System Subtype ‘Tiamat’ Count", "Prophage Defense System Subtype ‘Wadjet_I’ Count",
                         "Prophage Defense System Subtype ‘Wadjet_II’ Count", "Prophage Anti-Defense System Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-CRISPR: Aca_alone’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-CBASS: acb2’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-CRISPR: acriia13’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-CRISPR: acriia16’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-CRISPR: acriia17’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-CRISPR: acriia21’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-CRISPR: acriia3’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-RM: arda_ardu’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-RM: dam’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-RM: dar_ddr_hdf_ulx’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-RecBCD: gam’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-RM: mom’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-RM: ral’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-TA: rexb’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-Thoeris: tad1’ Count",
                         "Prophage Anti-Defense System Subtype ‘Anti-Thoeris: tad2_acriia7’ Count"])

df = pd.read_csv('results.tsv', sep='\t')
genome_name = df['BiosampleID']

for name in genome_name:
    defense_system_count = 0
    anti_defense_system_count = 0
    defense_system_1 = 0
    defense_system_2 = 0
    defense_system_3 = 0
    defense_system_4 = 0
    defense_system_5 = 0
    defense_system_6 = 0
    defense_system_7 = 0
    defense_system_8 = 0
    defense_system_9 = 0
    defense_system_10 = 0
    defense_system_11 = 0
    defense_system_12 = 0
    defense_system_13 = 0
    defense_system_14 = 0
    defense_system_15 = 0
    defense_system_16 = 0
    defense_system_17 = 0
    defense_system_18 = 0
    defense_system_19 = 0
    defense_system_20 = 0
    defense_system_21 = 0
    defense_system_22 = 0
    defense_system_23 = 0
    defense_system_24 = 0
    defense_system_25 = 0
    defense_system_26 = 0
    defense_system_27 = 0
    defense_system_28 = 0
    defense_system_29 = 0
    defense_system_30 = 0
    defense_system_31 = 0
    defense_system_32 = 0
    defense_system_33 = 0
    defense_system_34 = 0
    defense_system_35 = 0
    defense_system_36 = 0
    defense_system_37 = 0
    defense_system_38 = 0
    defense_system_39 = 0
    defense_system_40 = 0
    defense_system_41 = 0
    defense_system_42 = 0
    defense_system_43 = 0
    defense_system_44 = 0
    defense_system_45 = 0
    defense_system_46 = 0
    defense_system_47 = 0
    defense_system_48 = 0
    defense_system_49 = 0
    defense_system_50 = 0
    defense_system_51 = 0
    defense_system_52 = 0
    defense_system_53 = 0
    defense_system_54 = 0
    defense_system_55 = 0
    defense_system_56 = 0
    defense_system_57 = 0
    defense_system_58 = 0
    defense_system_59 = 0
    defense_system_60 = 0
    defense_system_61 = 0
    defense_system_62 = 0
    defense_system_63 = 0
    defense_system_64 = 0
    defense_system_65 = 0
    defense_system_66 = 0
    defense_system_67 = 0
    defense_system_68 = 0
    defense_system_69 = 0
    defense_system_70 = 0
    defense_system_71 = 0
    defense_system_72 = 0
    defense_system_73 = 0
    defense_system_74 = 0
    defense_system_75 = 0
    defense_system_76 = 0
    defense_system_77 = 0
    defense_system_78 = 0
    defense_system_79 = 0
    anti_defense_system_1 = 0
    anti_defense_system_2 = 0
    anti_defense_system_3 = 0
    anti_defense_system_4 = 0
    anti_defense_system_5 = 0
    anti_defense_system_6 = 0
    anti_defense_system_7 = 0
    anti_defense_system_8 = 0
    anti_defense_system_9 = 0
    anti_defense_system_10 = 0
    anti_defense_system_11 = 0
    anti_defense_system_12 = 0
    anti_defense_system_13 = 0
    anti_defense_system_14 = 0
    anti_defense_system_15 = 0
    anti_defense_system_16 = 0

    with open('prophage_phage_defenses.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['Genome ID'] == name:
                defense_system_count +=1
                if row['Defense Subtype'] == 'Abi2':
                    defense_system_1 += 1
                elif row['Defense Subtype'] == 'AbiC':
                    defense_system_2 += 1
                elif row['Defense Subtype'] == 'AbiD':
                    defense_system_3 += 1
                elif row['Defense Subtype'] == 'AbiE':
                    defense_system_4 += 1
                elif row['Defense Subtype'] == 'AbiG':
                    defense_system_5 += 1
                elif row['Defense Subtype'] == 'AbiH':
                    defense_system_6 += 1
                elif row['Defense Subtype'] == 'AbiI':
                    defense_system_7 += 1
                elif row['Defense Subtype'] == 'AbiJ':
                    defense_system_8 += 1
                elif row['Defense Subtype'] == 'AbiL':
                    defense_system_9 += 1
                elif row['Defense Subtype'] == 'AbiN':
                    defense_system_10 += 1
                elif row['Defense Subtype'] == 'AbiP2':
                    defense_system_11 += 1
                elif row['Defense Subtype'] == 'AbiQ': defense_system_12 += 1
                elif row['Defense Subtype'] == 'AbiR': defense_system_13 += 1
                elif row['Defense Subtype'] == 'AbiU': defense_system_14 += 1
                elif row['Defense Subtype'] == 'Avz_II': defense_system_15 += 1
                elif row['Defense Subtype'] == 'BREX_I': defense_system_16 += 1
                elif row['Defense Subtype'] == 'BREX_III': defense_system_17 += 1
                elif row['Defense Subtype'] == 'CAS_Class2-Subtype-II-A': defense_system_18 += 1
                elif row['Defense Subtype'] == 'CAS_Class2-Subtype-II-C': defense_system_19 += 1
                elif row['Defense Subtype'] == 'CBASS_I': defense_system_20 += 1
                elif row['Defense Subtype'] == 'CBASS_II': defense_system_21 += 1
                elif row['Defense Subtype'] == 'DarTG': defense_system_22 += 1
                elif row['Defense Subtype'] == 'DdmDE': defense_system_23 += 1
                elif row['Defense Subtype'] == 'Detocs': defense_system_24 += 1
                elif row['Defense Subtype'] == 'Dodola': defense_system_25 += 1
                elif row['Defense Subtype'] == 'DRT_2': defense_system_26 += 1
                elif row['Defense Subtype'] == 'Eleos': defense_system_27 += 1
                elif row['Defense Subtype'] == 'Gabija': defense_system_28 += 1
                elif row['Defense Subtype'] == 'Gao_Iet': defense_system_29 += 1
                elif row['Defense Subtype'] == 'Gao_Qat': defense_system_30 += 1
                elif row['Defense Subtype'] == 'Hachiman': defense_system_31 += 1
                elif row['Defense Subtype'] == 'Hna': defense_system_32 += 1
                elif row['Defense Subtype'] == 'Kiwa': defense_system_33 += 1
                elif row['Defense Subtype'] == 'Lamassu-Cap4_nuclease': defense_system_34 += 1
                elif row['Defense Subtype'] == 'Lamassu-Hypothetical': defense_system_35 += 1
                elif row['Defense Subtype'] == 'Lamassu-Lipase': defense_system_36 += 1
                elif row['Defense Subtype'] == 'Lamassu-Mrr': defense_system_37 += 1
                elif row['Defense Subtype'] == 'Lamassu-Protease': defense_system_38 += 1
                elif row['Defense Subtype'] == 'MazEF': defense_system_39 += 1
                elif row['Defense Subtype'] == 'Nhi': defense_system_40 += 1
                elif row['Defense Subtype'] == 'NLR_like_bNACHT01': defense_system_41 += 1
                elif row['Defense Subtype'] == 'Olokun': defense_system_42 += 1
                elif row['Defense Subtype'] == 'PARIS_I': defense_system_43 += 1
                elif row['Defense Subtype'] == 'PD-Lambda-1': defense_system_44 += 1
                elif row['Defense Subtype'] == 'PD-Lambda-5': defense_system_45 += 1
                elif row['Defense Subtype'] == 'PD-T4-1': defense_system_46 += 1
                elif row['Defense Subtype'] == 'PD-T4-8': defense_system_47 += 1
                elif row['Defense Subtype'] == 'PD-T4-9': defense_system_48 += 1
                elif row['Defense Subtype'] == 'PD-T7-2': defense_system_49 += 1
                elif row['Defense Subtype'] == 'PD-T7-3': defense_system_50 += 1
                elif row['Defense Subtype'] == 'PD-T7-4': defense_system_51 += 1
                elif row['Defense Subtype'] == 'PD-T7-5': defense_system_52 += 1
                elif row['Defense Subtype'] == 'Pycsar': defense_system_53 += 1
                elif row['Defense Subtype'] == 'Retron_I_C': defense_system_54 += 1
                elif row['Defense Subtype'] == 'Retron_III': defense_system_55 += 1
                elif row['Defense Subtype'] == 'Retron_XI': defense_system_56 += 1
                elif row['Defense Subtype'] == 'RexAB': defense_system_57 += 1
                elif row['Defense Subtype'] == 'RloC': defense_system_58 += 1
                elif row['Defense Subtype'] == 'RM_Type_I': defense_system_59 += 1
                elif row['Defense Subtype'] == 'RM_Type_II': defense_system_60 += 1
                elif row['Defense Subtype'] == 'RM_Type_IIG': defense_system_61 += 1
                elif row['Defense Subtype'] == 'RM_Type_III': defense_system_62 += 1
                elif row['Defense Subtype'] == 'RM_Type_IV': defense_system_63 += 1
                elif row['Defense Subtype'] == 'RnlAB': defense_system_64 += 1
                elif row['Defense Subtype'] == 'RosmerTA': defense_system_65 += 1
                elif row['Defense Subtype'] == 'Rst_HelicaseDUF2290': defense_system_66 += 1
                elif row['Defense Subtype'] == 'SanaTA': defense_system_67 += 1
                elif row['Defense Subtype'] == 'SEFIR': defense_system_68 += 1
                elif row['Defense Subtype'] == 'Septu': defense_system_69 += 1
                elif row['Defense Subtype'] == 'Shedu': defense_system_70 += 1
                elif row['Defense Subtype'] == 'ShosTA': defense_system_71 += 1
                elif row['Defense Subtype'] == 'SoFic': defense_system_72 += 1
                elif row['Defense Subtype'] == 'SpbK': defense_system_73 += 1
                elif row['Defense Subtype'] == 'Stk2': defense_system_74 += 1
                elif row['Defense Subtype'] == 'Thoeris_I': defense_system_75 += 1
                elif row['Defense Subtype'] == 'Thoeris_II': defense_system_76 += 1
                elif row['Defense Subtype'] == 'Tiamat': defense_system_77 += 1
                elif row['Defense Subtype'] == 'Wadjet_I': defense_system_78 += 1
                elif row['Defense Subtype'] == 'Wadjet_II': defense_system_79 += 1

    with open('prophage_anti_defenses.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['Genome ID'] == name:
                anti_defense_system_count +=1
                if row['Defense Subtype'] == 'Aca_alone':
                    anti_defense_system_1 += 1
                elif row['Defense Subtype'] == 'acb2':
                    anti_defense_system_2 += 1
                elif row['Defense Subtype'] == 'acriia13':
                    anti_defense_system_3 += 1
                elif row['Defense Subtype'] == 'acriia16':
                    anti_defense_system_4 += 1
                elif row['Defense Subtype'] == 'acriia17':
                    anti_defense_system_5 += 1
                elif row['Defense Subtype'] == 'acriia21':
                    anti_defense_system_6 += 1
                elif row['Defense Subtype'] == 'acriia3':
                    anti_defense_system_7 += 1
                elif row['Defense Subtype'] == 'arda_ardu':
                    anti_defense_system_8 += 1
                elif row['Defense Subtype'] == 'dam':
                    anti_defense_system_9 += 1
                elif row['Defense Subtype'] == 'dar_ddr_hdf_ulx':
                    anti_defense_system_10 += 1
                elif row['Defense Subtype'] == 'gam':
                    anti_defense_system_11 += 1
                elif row['Defense Subtype'] == 'mom':
                    anti_defense_system_12 += 1
                elif row['Defense Subtype'] == 'ral':
                    anti_defense_system_13 += 1
                elif row['Defense Subtype'] == 'rexb':
                    anti_defense_system_14 += 1
                elif row['Defense Subtype'] == 'tad1':
                    anti_defense_system_15 += 1
                elif row['Defense Subtype'] == 'tad2_acriia7':
                    anti_defense_system_16 += 1


    with open('bacteriaprophagedefensecount.tsv', 'a') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow([f'{name}', f'{defense_system_count}',
                             f'{defense_system_1}', f'{defense_system_2}', f'{defense_system_3}', f'{defense_system_4}', f'{defense_system_5}',
                             f'{defense_system_6}', f'{defense_system_7}', f'{defense_system_8}', f'{defense_system_9}', f'{defense_system_10}',
                             f'{defense_system_11}', f'{defense_system_12}', f'{defense_system_13}', f'{defense_system_14}', f'{defense_system_15}',
                             f'{defense_system_16}', f'{defense_system_17}', f'{defense_system_18}', f'{defense_system_19}', f'{defense_system_20}',
                             f'{defense_system_21}', f'{defense_system_22}', f'{defense_system_23}', f'{defense_system_24}', f'{defense_system_25}',
                             f'{defense_system_26}', f'{defense_system_27}', f'{defense_system_28}', f'{defense_system_29}', f'{defense_system_30}',
                             f'{defense_system_31}', f'{defense_system_32}', f'{defense_system_33}', f'{defense_system_34}', f'{defense_system_35}',
                             f'{defense_system_36}', f'{defense_system_37}', f'{defense_system_38}', f'{defense_system_39}', f'{defense_system_40}',
                             f'{defense_system_41}', f'{defense_system_42}', f'{defense_system_43}', f'{defense_system_44}', f'{defense_system_45}',
                             f'{defense_system_46}', f'{defense_system_47}', f'{defense_system_48}', f'{defense_system_49}', f'{defense_system_50}',
                             f'{defense_system_51}', f'{defense_system_52}', f'{defense_system_53}', f'{defense_system_54}', f'{defense_system_55}',
                             f'{defense_system_56}', f'{defense_system_57}', f'{defense_system_58}', f'{defense_system_59}', f'{defense_system_60}',
                             f'{defense_system_61}', f'{defense_system_62}', f'{defense_system_63}', f'{defense_system_64}', f'{defense_system_65}',
                             f'{defense_system_66}', f'{defense_system_67}', f'{defense_system_68}', f'{defense_system_69}', f'{defense_system_70}',
                             f'{defense_system_71}', f'{defense_system_72}', f'{defense_system_73}',
                             f'{defense_system_74}', f'{defense_system_75}',
                             f'{defense_system_76}', f'{defense_system_77}', f'{defense_system_78}',
                             f'{defense_system_79}', f'{anti_defense_system_count}',
                             f'{anti_defense_system_1}', f'{anti_defense_system_2}', f'{anti_defense_system_3}', f'{anti_defense_system_4}', f'{anti_defense_system_5}',
                             f'{anti_defense_system_6}', f'{anti_defense_system_7}', f'{anti_defense_system_8}', f'{anti_defense_system_9}', f'{anti_defense_system_10}',
                             f'{anti_defense_system_11}', f'{anti_defense_system_12}', f'{anti_defense_system_13}', f'{anti_defense_system_14}', f'{anti_defense_system_15}',
                             f'{anti_defense_system_16}'])
