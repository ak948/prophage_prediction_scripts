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

with open('bacteriacount.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(["BiosampleID", 
                         "Virulence Gene ‘ace’ Count", "Virulence Gene ‘acm’ Count", "Virulence Gene ‘asa1’ Count",
                         "Virulence Gene ‘bopD’ Count", "Virulence Gene ‘cap8B’ Count", "Virulence Gene ‘cap8M’ Count",
                         "Virulence Gene ‘capN’ Count", "Virulence Gene ‘cpsA/uppS’ Count",
                         "Virulence Gene ‘cpsB/cdsA’ Count", "Virulence Gene ‘cpsC’ Count",
                         "Virulence Gene ‘cpsD’ Count", "Virulence Gene ‘cpsE’ Count", "Virulence Gene ‘cpsF’ Count",
                         "Virulence Gene ‘cpsG’ Count", "Virulence Gene ‘cpsH’ Count", "Virulence Gene ‘cpsI’ Count",
                         "Virulence Gene ‘cpsJ’ Count", "Virulence Gene ‘cpsK’ Count", "Virulence Gene ‘csgA’ Count",
                         "Virulence Gene ‘cylA’ Count", "Virulence Gene ‘cylB’ Count", "Virulence Gene ‘cylI’ Count",
                         "Virulence Gene ‘cylL’ Count", "Virulence Gene ‘cylM’ Count", "Virulence Gene ‘cylR1’ Count",
                         "Virulence Gene ‘cylR2’ Count", "Virulence Gene ‘cylS’ Count", "Virulence Gene ‘ebpA’ Count",
                         "Virulence Gene ‘ebpB’ Count", "Virulence Gene ‘ebpC’ Count",
                         "Virulence Gene ‘ecbA/fss3’ Count", "Virulence Gene ‘EF0149’ Count",
                         "Virulence Gene ‘EF0485’ Count", "Virulence Gene ‘EF0818’ Count",
                         "Virulence Gene ‘EF3023’ Count", "Virulence Gene ‘efaA’ Count", "Virulence Gene ‘esp’ Count",
                         "Virulence Gene ‘fimE’ Count", "Virulence Gene ‘fsrA’ Count", "Virulence Gene ‘fsrB’ Count",
                         "Virulence Gene ‘fsrC’ Count", "Virulence Gene ‘gelE’ Count", "Virulence Gene ‘gspG’ Count",
                         "Virulence Gene ‘hcp1/tssD1’ Count", "Virulence Gene ‘hlb’ Count",
                         "Virulence Gene ‘KP1_RS17280’ Count", "Virulence Gene ‘kpsT’ Count",
                         "Virulence Gene ‘neuB’ Count", "Virulence Gene ‘neuD’ Count", "Virulence Gene ‘papX’ Count",
                         "Virulence Gene ‘phzC1’ Count", "Virulence Gene ‘phzD1’ Count", "Virulence Gene ‘phzE1’ Count",
                         "Virulence Gene ‘phzF1’ Count", "Virulence Gene ‘phzG2’ Count", "Virulence Gene ‘pic’ Count",
                         "Virulence Gene ‘pilO’ Count", "Virulence Gene ‘prgB/asc10’ Count",
                         "Virulence Gene ‘rfbK1’ Count", "Virulence Gene ‘scm’ Count", "Virulence Gene ‘sgrA’ Count",
                         "Virulence Gene ‘sprE’ Count", "Virulence Gene ‘srtC’ Count", "Virulence Gene ‘tssA’ Count",
                         "Virulence Gene ‘vgrG’ Count", "Virulence Gene ‘ykgK/ecpR’ Count"])

df = pd.read_csv('filtered_results.tsv', sep='\t')
genome_name = df['BiosampleID']

for name in genome_name:
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
    defense_system_80 = 0
    defense_system_81 = 0
    defense_system_82 = 0
    defense_system_83 = 0
    defense_system_84 = 0
    defense_system_85 = 0
    defense_system_86 = 0
    defense_system_87 = 0
    defense_system_88 = 0
    defense_system_89 = 0
    defense_system_90 = 0
    defense_system_91 = 0
    defense_system_92 = 0
    defense_system_93 = 0
    defense_system_94 = 0
    defense_system_95 = 0
    defense_system_96 = 0
    defense_system_97 = 0
    defense_system_98 = 0
    defense_system_99 = 0
    defense_system_100 = 0
    defense_system_101 = 0
    defense_system_102 = 0
    defense_system_103 = 0
    defense_system_104 = 0
    defense_system_105 = 0
    defense_system_106 = 0
    defense_system_107 = 0
    defense_system_108 = 0
    defense_system_109 = 0
    defense_system_110 = 0
    defense_system_111 = 0
    defense_system_112 = 0
    defense_system_113 = 0
    defense_system_114 = 0
    defense_system_115 = 0
    defense_system_116 = 0
    defense_system_117 = 0
    defense_system_118 = 0
    defense_system_119 = 0
    defense_system_120 = 0
    defense_system_121 = 0
    defense_system_122 = 0
    defense_system_123 = 0
    defense_system_124 = 0
    defense_system_125 = 0
    defense_system_126 = 0
    defense_system_127 = 0
    defense_system_128 = 0
    defense_system_129 = 0
    defense_system_130 = 0
    defense_system_131 = 0
    defense_system_132 = 0
    defense_system_133 = 0
    defense_system_134 = 0
    defense_system_135 = 0
    defense_system_136 = 0
    vir_gene_1 = 0
    vir_gene_2 = 0
    vir_gene_3 = 0
    vir_gene_4 = 0
    vir_gene_5 = 0
    vir_gene_6 = 0
    vir_gene_7 = 0
    vir_gene_8 = 0
    vir_gene_9 = 0
    vir_gene_10 = 0
    vir_gene_11 = 0
    vir_gene_12 = 0
    vir_gene_13 = 0
    vir_gene_14 = 0
    vir_gene_15 = 0
    vir_gene_16 = 0
    vir_gene_17 = 0
    vir_gene_18 = 0
    vir_gene_19 = 0
    vir_gene_20 = 0
    vir_gene_21 = 0
    vir_gene_22 = 0
    vir_gene_23 = 0
    vir_gene_24 = 0
    vir_gene_25 = 0
    vir_gene_26 = 0
    vir_gene_27 = 0
    vir_gene_28 = 0
    vir_gene_29 = 0
    vir_gene_30 = 0
    vir_gene_31 = 0
    vir_gene_32 = 0
    vir_gene_33 = 0
    vir_gene_34 = 0
    vir_gene_35 = 0
    vir_gene_36 = 0
    vir_gene_37 = 0
    vir_gene_38 = 0
    vir_gene_39 = 0
    vir_gene_40 = 0
    vir_gene_41 = 0
    vir_gene_42 = 0
    vir_gene_43 = 0
    vir_gene_44 = 0
    vir_gene_45 = 0
    vir_gene_46 = 0
    vir_gene_47 = 0
    vir_gene_48 = 0
    vir_gene_49 = 0
    vir_gene_50 = 0
    vir_gene_51 = 0
    vir_gene_52 = 0
    vir_gene_53 = 0
    vir_gene_54 = 0
    vir_gene_55 = 0
    vir_gene_56 = 0
    vir_gene_57 = 0
    vir_gene_58 = 0
    vir_gene_59 = 0
    vir_gene_60 = 0
    vir_gene_61 = 0
    vir_gene_62 = 0
    vir_gene_63 = 0
    vir_gene_64 = 0
    vir_gene_65 = 0
    vir_gene_66 = 0
    amr_gene_1 = 0
    amr_gene_2 = 0
    amr_gene_3 = 0
    amr_gene_4 = 0
    amr_gene_5 = 0
    amr_gene_6 = 0
    amr_gene_7 = 0
    amr_gene_8 = 0
    amr_gene_9 = 0
    amr_gene_10 = 0
    amr_gene_11 = 0
    amr_gene_12 = 0
    amr_gene_13 = 0
    amr_gene_14 = 0
    amr_gene_15 = 0
    amr_gene_16 = 0
    amr_gene_17 = 0
    amr_gene_18 = 0
    amr_gene_19 = 0
    amr_gene_20 = 0
    amr_gene_21 = 0
    amr_gene_22 = 0
    amr_gene_23 = 0
    amr_gene_24 = 0
    amr_gene_25 = 0
    amr_gene_26 = 0
    amr_gene_27 = 0
    amr_gene_28 = 0
    amr_gene_29 = 0
    amr_gene_30 = 0
    amr_gene_31 = 0
    amr_gene_32 = 0
    amr_gene_33 = 0
    amr_gene_34 = 0
    amr_gene_35 = 0
    amr_gene_36 = 0
    amr_gene_37 = 0
    amr_gene_38 = 0
    amr_gene_39 = 0
    amr_gene_40 = 0
    amr_gene_41 = 0
    amr_gene_42 = 0
    amr_gene_43 = 0
    amr_gene_44 = 0
    amr_gene_45 = 0
    amr_gene_46 = 0
    amr_gene_47 = 0
    amr_gene_48 = 0
    amr_gene_49 = 0
    amr_gene_50 = 0
    amr_gene_51 = 0
    amr_gene_52 = 0
    amr_gene_53 = 0
    amr_gene_54 = 0
    amr_gene_55 = 0
    amr_gene_56 = 0
    amr_gene_57 = 0
    amr_gene_58 = 0
    amr_gene_59 = 0
    amr_gene_60 = 0
    amr_gene_61 = 0
    amr_gene_62 = 0
    amr_gene_63 = 0
    amr_gene_64 = 0
    amr_gene_65 = 0

    # with open('bacteria_phage_defenses.tsv', 'r') as file:
    #     reader = csv.DictReader(file, delimiter='\t')
    #     for row in reader:
    #         if row['Genome ID'] == name:
    #             if row['Defense Subtype'] == 'Abi2':
    #                 defense_system_1 += 1
    #             elif row['Defense Subtype'] == 'AbiA_large':
    #                 defense_system_2 += 1
    #             elif row['Defense Subtype'] == 'AbiB':
    #                 defense_system_3 += 1
    #             elif row['Defense Subtype'] == 'AbiC':
    #                 defense_system_4 += 1
    #             elif row['Defense Subtype'] == 'AbiD':
    #                 defense_system_5 += 1
    #             elif row['Defense Subtype'] == 'AbiE':
    #                 defense_system_6 += 1
    #             elif row['Defense Subtype'] == 'AbiG':
    #                 defense_system_7 += 1
    #             elif row['Defense Subtype'] == 'AbiH':
    #                 defense_system_8 += 1
    #             elif row['Defense Subtype'] == 'AbiI':
    #                 defense_system_9 += 1
    #             elif row['Defense Subtype'] == 'AbiJ':
    #                 defense_system_10 += 1
    #             elif row['Defense Subtype'] == 'AbiK':
    #                 defense_system_11 += 1
    #             elif row['Defense Subtype'] == 'AbiL': defense_system_12 += 1
    #             elif row['Defense Subtype'] == 'AbiN': defense_system_13 += 1
    #             elif row['Defense Subtype'] == 'AbiO': defense_system_14 += 1
    #             elif row['Defense Subtype'] == 'AbiP2': defense_system_15 += 1
    #             elif row['Defense Subtype'] == 'AbiQ': defense_system_16 += 1
    #             elif row['Defense Subtype'] == 'AbiR': defense_system_17 += 1
    #             elif row['Defense Subtype'] == 'AbiU': defense_system_18 += 1
    #             elif row['Defense Subtype'] == 'AbiV': defense_system_19 += 1
    #             elif row['Defense Subtype'] == 'Aditi': defense_system_20 += 1
    #             elif row['Defense Subtype'] == 'Avs_II': defense_system_21 += 1
    #             elif row['Defense Subtype'] == 'Avs_IV': defense_system_22 += 1
    #             elif row['Defense Subtype'] == 'Avs_V': defense_system_23 += 1
    #             elif row['Defense Subtype'] == 'Azaca': defense_system_24 += 1
    #             elif row['Defense Subtype'] == 'BREX_I': defense_system_25 += 1
    #             elif row['Defense Subtype'] == 'BREX_III': defense_system_26 += 1
    #             elif row['Defense Subtype'] == 'Bunzi': defense_system_27 += 1
    #             elif row['Defense Subtype'] == 'CAS_Class1-Subtype-I-B': defense_system_28 += 1  #######################################
    #             elif row['Defense Subtype'] == 'CAS_Class1-Subtype-I-C': defense_system_29 += 1 ###################################
    #             elif row['Defense Subtype'] == 'CAS_Class1-Subtype-I-E': defense_system_30 += 1 ####################################
    #             elif row['Defense Subtype'] == 'CAS_Class1-Subtype-III-A': defense_system_31 += 1 ###################################
    #             elif row['Defense Subtype'] == 'CAS_Class1-Subtype-III-B': defense_system_32 += 1 ###################################
    #             elif row['Defense Subtype'] == 'CAS_Class1-Subtype-III-D': defense_system_33 += 1 ###################################
    #             elif row['Defense Subtype'] == 'CAS_Class1-Type-III': defense_system_34 += 1 ############################################
    #             elif row['Defense Subtype'] == 'CAS_Class2-Subtype-II-A': defense_system_35 += 1 ####################################
    #             elif row['Defense Subtype'] == 'CAS_Class2-Subtype-II-C': defense_system_36 += 1 ##########################################
    #             elif row['Defense Subtype'] == 'CAS_Cluster': defense_system_37 += 1
    #             elif row['Defense Subtype'] == 'CBASS_I': defense_system_38 += 1
    #             elif row['Defense Subtype'] == 'CBASS_II': defense_system_39 += 1
    #             elif row['Defense Subtype'] == 'DarTG': defense_system_40 += 1
    #             elif row['Defense Subtype'] == 'DdmDE': defense_system_41 += 1
    #             elif row['Defense Subtype'] == 'Detocs': defense_system_42 += 1
    #             elif row['Defense Subtype'] == 'DISARM_2': defense_system_43 += 1
    #             elif row['Defense Subtype'] == 'Dnd_ABCDEFGH': defense_system_44 += 1
    #             elif row['Defense Subtype'] == 'Dodola': defense_system_45 += 1
    #             elif row['Defense Subtype'] == 'DRT_2': defense_system_46 += 1
    #             elif row['Defense Subtype'] == 'DRT_4': defense_system_47 += 1
    #             elif row['Defense Subtype'] == 'DRT6': defense_system_48 += 1
    #             elif row['Defense Subtype'] == 'DRT8': defense_system_49 += 1
    #             elif row['Defense Subtype'] == 'DRT9': defense_system_50 += 1
    #             elif row['Defense Subtype'] == 'Druantia_III': defense_system_51 += 1
    #             elif row['Defense Subtype'] == 'Dsr_II': defense_system_52 += 1
    #             elif row['Defense Subtype'] == 'Eleos': defense_system_53 += 1
    #             elif row['Defense Subtype'] == 'FS_HsdR_like': defense_system_54 += 1
    #             elif row['Defense Subtype'] == 'FS_Sma': defense_system_55 += 1
    #             elif row['Defense Subtype'] == 'Gabija': defense_system_56 += 1
    #             elif row['Defense Subtype'] == 'Gao_Her_DUF': defense_system_57 += 1
    #             elif row['Defense Subtype'] == 'Gao_Iet': defense_system_58 += 1
    #             elif row['Defense Subtype'] == 'Gao_Ppl': defense_system_59 += 1
    #             elif row['Defense Subtype'] == 'Gao_Qat': defense_system_60 += 1
    #             elif row['Defense Subtype'] == 'GAPS4': defense_system_61 += 1
    #             elif row['Defense Subtype'] == 'Hachiman': defense_system_62 += 1
    #             elif row['Defense Subtype'] == 'Hna': defense_system_63 += 1
    #             elif row['Defense Subtype'] == 'Kiwa': defense_system_64 += 1
    #             elif row['Defense Subtype'] == 'Lamassu-Cap4_nuclease': defense_system_65 += 1 #################################
    #             elif row['Defense Subtype'] == 'Lamassu-Hydrolase': defense_system_66 += 1
    #             elif row['Defense Subtype'] == 'Lamassu-Hydrolase_Protease': defense_system_67 += 1
    #             elif row['Defense Subtype'] == 'Lamassu-Hypothetical': defense_system_68 += 1
    #             elif row['Defense Subtype'] == 'Lamassu-Lipase': defense_system_69 += 1
    #             elif row['Defense Subtype'] == 'Lamassu-Mrr': defense_system_70 += 1
    #             elif row['Defense Subtype'] == 'Lamassu-Protease': defense_system_71 += 1
    #             elif row['Defense Subtype'] == 'Lamassu-Sir2': defense_system_72 += 1 ######################################
    #             elif row['Defense Subtype'] == 'MazEF': defense_system_73 += 1
    #             elif row['Defense Subtype'] == 'Menshen': defense_system_74 += 1
    #             elif row['Defense Subtype'] == 'Mok_Hok_Sok': defense_system_75 += 1
    #             elif row['Defense Subtype'] == 'Mokosh_Type_I_A': defense_system_76 += 1
    #             elif row['Defense Subtype'] == 'Mokosh_Type_I_B': defense_system_77 += 1
    #             elif row['Defense Subtype'] == 'Mokosh_TypeII': defense_system_78 += 1
    #             elif row['Defense Subtype'] == 'Nhi': defense_system_79 += 1
    #             elif row['Defense Subtype'] == 'NLR_like_bNACHT01': defense_system_80 += 1
    #             elif row['Defense Subtype'] == 'NLR_like_bNACHT09': defense_system_81 += 1
    #             elif row['Defense Subtype'] == 'Olokun': defense_system_82 += 1
    #             elif row['Defense Subtype'] == 'pAgo_S1B': defense_system_83 += 1
    #             elif row['Defense Subtype'] == 'PARIS_I': defense_system_84 += 1
    #             elif row['Defense Subtype'] == 'PARIS_II': defense_system_85 += 1
    #             elif row['Defense Subtype'] == 'PARIS_II_merge': defense_system_86 += 1
    #             elif row['Defense Subtype'] == 'PD-Lambda-1': defense_system_87 += 1 #############################################
    #             elif row['Defense Subtype'] == 'PD-Lambda-5': defense_system_88 += 1
    #             elif row['Defense Subtype'] == 'PD-T4-1': defense_system_89 += 1
    #             elif row['Defense Subtype'] == 'PD-T4-2': defense_system_90 += 1
    #             elif row['Defense Subtype'] == 'PD-T4-3': defense_system_91 += 1
    #             elif row['Defense Subtype'] == 'PD-T4-5': defense_system_92 += 1
    #             elif row['Defense Subtype'] == 'PD-T4-7': defense_system_93 += 1
    #             elif row['Defense Subtype'] == 'PD-T4-8': defense_system_94 += 1
    #             elif row['Defense Subtype'] == 'PD-T4-9': defense_system_95 += 1
    #             elif row['Defense Subtype'] == 'PD-T7-2': defense_system_96 += 1
    #             elif row['Defense Subtype'] == 'PD-T7-3': defense_system_97 += 1
    #             elif row['Defense Subtype'] == 'PD-T7-4': defense_system_98 += 1
    #             elif row['Defense Subtype'] == 'PD-T7-5': defense_system_99 += 1 #####################################################
    #             elif row['Defense Subtype'] == 'Pelif': defense_system_100 += 1
    #             elif row['Defense Subtype'] == 'PrrC': defense_system_101 += 1
    #             elif row['Defense Subtype'] == 'Pycsar': defense_system_102 += 1
    #             elif row['Defense Subtype'] == 'Retron_I_B': defense_system_103 += 1
    #             elif row['Defense Subtype'] == 'Retron_I_C': defense_system_104 += 1
    #             elif row['Defense Subtype'] == 'Retron_II': defense_system_105 += 1
    #             elif row['Defense Subtype'] == 'Retron_III': defense_system_106 += 1
    #             elif row['Defense Subtype'] == 'Retron_IV': defense_system_107 += 1
    #             elif row['Defense Subtype'] == 'Retron_VI': defense_system_108 += 1
    #             elif row['Defense Subtype'] == 'Retron_XI': defense_system_109 += 1
    #             elif row['Defense Subtype'] == 'Retron_XII': defense_system_110 += 1
    #             elif row['Defense Subtype'] == 'RexAB': defense_system_111 += 1
    #             elif row['Defense Subtype'] == 'RloC': defense_system_112 += 1
    #             elif row['Defense Subtype'] == 'RM_Type_I': defense_system_113 += 1
    #             elif row['Defense Subtype'] == 'RM_Type_II': defense_system_114 += 1
    #             elif row['Defense Subtype'] == 'RM_Type_IIG': defense_system_115 += 1
    #             elif row['Defense Subtype'] == 'RM_Type_III': defense_system_116 += 1
    #             elif row['Defense Subtype'] == 'RM_Type_IV': defense_system_117 += 1
    #             elif row['Defense Subtype'] == 'RnlAB': defense_system_118 += 1
    #             elif row['Defense Subtype'] == 'RosmerTA': defense_system_119 += 1
    #             elif row['Defense Subtype'] == 'Rst_HelicaseDUF2290': defense_system_120 += 1
    #             elif row['Defense Subtype'] == 'Rst_TIR-NLR': defense_system_121 += 1 #############################################
    #             elif row['Defense Subtype'] == 'SanaTA': defense_system_122 += 1
    #             elif row['Defense Subtype'] == 'SEFIR': defense_system_123 += 1
    #             elif row['Defense Subtype'] == 'Septu': defense_system_124 += 1
    #             elif row['Defense Subtype'] == 'Shango': defense_system_125 += 1
    #             elif row['Defense Subtype'] == 'Shedu': defense_system_126 += 1
    #             elif row['Defense Subtype'] == 'ShosTA': defense_system_127 += 1
    #             elif row['Defense Subtype'] == 'SoFic': defense_system_128 += 1
    #             elif row['Defense Subtype'] == 'SpbK': defense_system_129 += 1
    #             elif row['Defense Subtype'] == 'Stk2': defense_system_130 += 1
    #             elif row['Defense Subtype'] == 'Thoeris_I': defense_system_131 += 1
    #             elif row['Defense Subtype'] == 'Thoeris_II': defense_system_132 += 1
    #             elif row['Defense Subtype'] == 'Tiamat': defense_system_133 += 1
    #             elif row['Defense Subtype'] == 'Uzume': defense_system_134 += 1
    #             elif row['Defense Subtype'] == 'Wadjet_I': defense_system_135 += 1
    #             elif row['Defense Subtype'] == 'Wadjet_II': defense_system_136 +=1

    with open('bacteria_virulence_genes.tsv', 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['Genome ID'] == name:
                if row['Gene'] == 'ace': vir_gene_1 += 1
                elif row['Gene'] == 'acm': vir_gene_2 += 1
                elif row['Gene'] == 'asa1': vir_gene_3 += 1
                elif row['Gene'] == 'bopD': vir_gene_4 += 1
                elif row['Gene'] == 'cap8B': vir_gene_5 += 1
                elif row['Gene'] == 'cap8M': vir_gene_6 += 1
                elif row['Gene'] == 'capN': vir_gene_7 += 1
                elif row['Gene'] == 'cpsA/uppS': vir_gene_8 += 1
                elif row['Gene'] == 'cpsB/cdsA': vir_gene_9 += 1
                elif row['Gene'] == 'cpsC': vir_gene_10 += 1
                elif row['Gene'] == 'cpsD': vir_gene_11 += 1
                elif row['Gene'] == 'cpsE': vir_gene_12 += 1
                elif row['Gene'] == 'cpsF': vir_gene_13 += 1
                elif row['Gene'] == 'cpsG': vir_gene_14 += 1
                elif row['Gene'] == 'cpsH': vir_gene_15 += 1
                elif row['Gene'] == 'cpsI': vir_gene_16 += 1
                elif row['Gene'] == 'cpsJ': vir_gene_17 += 1
                elif row['Gene'] == 'cpsK': vir_gene_18 += 1
                elif row['Gene'] == 'csgA': vir_gene_19 += 1
                elif row['Gene'] == 'cylA': vir_gene_20 += 1
                elif row['Gene'] == 'cylB': vir_gene_21 += 1
                elif row['Gene'] == 'cylI': vir_gene_22 += 1
                elif row['Gene'] == 'cylL': vir_gene_23 += 1
                elif row['Gene'] == 'cylM': vir_gene_24 += 1
                elif row['Gene'] == 'cylR1': vir_gene_25 += 1
                elif row['Gene'] == 'cylR2': vir_gene_26 += 1
                elif row['Gene'] == 'cylS': vir_gene_27 += 1
                elif row['Gene'] == 'ebpA': vir_gene_28 += 1
                elif row['Gene'] == 'ebpB': vir_gene_29 += 1
                elif row['Gene'] == 'ebpC': vir_gene_30 += 1
                elif row['Gene'] == 'ecbA/fss3': vir_gene_31 += 1
                elif row['Gene'] == 'EF0149': vir_gene_32 += 1
                elif row['Gene'] == 'EF0485': vir_gene_33 += 1
                elif row['Gene'] == 'EF0818': vir_gene_34 += 1
                elif row['Gene'] == 'EF3023': vir_gene_35 += 1
                elif row['Gene'] == 'efaA': vir_gene_36 += 1
                elif row['Gene'] == 'esp': vir_gene_37 += 1
                elif row['Gene'] == 'fimE': vir_gene_38 += 1
                elif row['Gene'] == 'fsrA': vir_gene_39 += 1
                elif row['Gene'] == 'fsrB': vir_gene_40 += 1
                elif row['Gene'] == 'fsrC': vir_gene_41 += 1
                elif row['Gene'] == 'gelE': vir_gene_42 += 1
                elif row['Gene'] == 'gspG': vir_gene_43 += 1
                elif row['Gene'] == 'hcp1/tssD1': vir_gene_44 += 1
                elif row['Gene'] == 'hlb': vir_gene_45 += 1
                elif row['Gene'] == 'KP1_RS17280': vir_gene_46 += 1
                elif row['Gene'] == 'kpsT': vir_gene_47 += 1
                elif row['Gene'] == 'neuB': vir_gene_48 += 1
                elif row['Gene'] == 'neuD': vir_gene_49 += 1
                elif row['Gene'] == 'papX': vir_gene_50 += 1
                elif row['Gene'] == 'phzC1': vir_gene_51 += 1
                elif row['Gene'] == 'phzD1': vir_gene_52 += 1
                elif row['Gene'] == 'phzE1': vir_gene_53 += 1
                elif row['Gene'] == 'phzF1': vir_gene_54 += 1
                elif row['Gene'] == 'phzG2': vir_gene_55 += 1
                elif row['Gene'] == 'pic': vir_gene_56 += 1
                elif row['Gene'] == 'pilO': vir_gene_57 += 1
                elif row['Gene'] == 'prgB/asc10': vir_gene_58 += 1
                elif row['Gene'] == 'rfbK1': vir_gene_59 += 1
                elif row['Gene'] == 'scm': vir_gene_60 += 1
                elif row['Gene'] == 'sgrA': vir_gene_61 += 1
                elif row['Gene'] == 'sprE': vir_gene_62 += 1
                elif row['Gene'] == 'srtC': vir_gene_63 += 1
                elif row['Gene'] == 'tssA': vir_gene_64 += 1
                elif row['Gene'] == 'vgrG': vir_gene_65 += 1
                elif row['Gene'] == 'ykgK/ecpR': vir_gene_66 +=1

    # with open('bacteria_amr_genes.tsv', 'r') as file:
    #     reader = csv.DictReader(file, delimiter='\t')
    #     for row in reader:
    #         if f'{name}_' in row['AMR Gene ID']:
    #             if row['AMR Gene Family'] == "23S rRNA with mutation conferring resistance to macrolide antibiotics": amr_gene_1 += 1
    #             elif row["AMR Gene Family"] == "AAC(2')": amr_gene_2 += 1
    #             elif row["AMR Gene Family"] == "AAC(3)": amr_gene_3 += 1
    #             elif row["AMR Gene Family"] == "AAC(6')": amr_gene_4 += 1
    #             elif row["AMR Gene Family"] == "ADC beta-lactamases pending classification for carbapenemase activity": amr_gene_5 += 1
    #             elif row["AMR Gene Family"] == "aminoglycoside bifunctional resistance protein": amr_gene_6 += 1
    #             elif row["AMR Gene Family"] == "ANT(4')": amr_gene_7 += 1
    #             elif row["AMR Gene Family"] == "ANT(6)": amr_gene_8 += 1
    #             elif row["AMR Gene Family"] == "ANT(9)": amr_gene_9 += 1
    #             elif row["AMR Gene Family"] == "APH(2'')": amr_gene_10 += 1
    #             elif row["AMR Gene Family"] == "APH(3')": amr_gene_11 += 1
    #             elif row["AMR Gene Family"] == "ATP-binding cassette (ABC) antibiotic efflux pump": amr_gene_12 += 1
    #             elif row["AMR Gene Family"] == "ATP-binding cassette (ABC) antibiotic efflux pump; major facilitator superfamily (MFS) antibiotic efflux pump; resistance-nodulation-cell division (RND) antibiotic efflux pump; General Bacterial Porin with reduced permeability to beta-lactams": amr_gene_13 += 1
    #             elif row["AMR Gene Family"] == "BlaZ beta-lactamase": amr_gene_14 += 1
    #             elif row["AMR Gene Family"] == "Cfr 23S ribosomal RNA methyltransferase": amr_gene_15 += 1
    #             elif row["AMR Gene Family"] == "CfxA beta-lactamase": amr_gene_16 += 1
    #             elif row["AMR Gene Family"] == "chloramphenicol acetyltransferase (CAT)": amr_gene_17 += 1
    #             elif row["AMR Gene Family"] == "daptomycin resistant cls": amr_gene_18 += 1
    #             elif row["AMR Gene Family"] == "daptomycin resistant gshF": amr_gene_19 += 1
    #             elif row["AMR Gene Family"] == "daptomycin resistant liaF": amr_gene_20 += 1
    #             elif row["AMR Gene Family"] == "daptomycin resistant liaR": amr_gene_21 += 1
    #             elif row["AMR Gene Family"] == "daptomycin resistant liaS": amr_gene_22 += 1
    #             elif row["AMR Gene Family"] == "daptomycin resistant YybT": amr_gene_23 += 1
    #             elif row["AMR Gene Family"] == "EC beta-lactamase": amr_gene_24 += 1
    #             elif row["AMR Gene Family"] == "elfamycin resistant EF-Tu": amr_gene_25 += 1
    #             elif row["AMR Gene Family"] == "Erm 23S ribosomal RNA methyltransferase": amr_gene_26 += 1
    #             elif row["AMR Gene Family"] == "fosfomycin thiol transferase": amr_gene_27 += 1
    #             elif row["AMR Gene Family"] == "glycopeptide resistance gene cluster; Van ligase": amr_gene_28 += 1
    #             elif row["AMR Gene Family"] == "glycopeptide resistance gene cluster; vanR": amr_gene_29 += 1
    #             elif row["AMR Gene Family"] == "glycopeptide resistance gene cluster; vanT": amr_gene_30 += 1
    #             elif row["AMR Gene Family"] == "glycopeptide resistance gene cluster; vanU": amr_gene_31 += 1
    #             elif row["AMR Gene Family"] == "glycopeptide resistance gene cluster; vanV": amr_gene_32 += 1
    #             elif row["AMR Gene Family"] == "glycopeptide resistance gene cluster; vanXY": amr_gene_33 += 1
    #             elif row["AMR Gene Family"] == "IND beta-lactamase": amr_gene_34 += 1
    #             elif row["AMR Gene Family"] == "lincosamide nucleotidyltransferase (LNU)": amr_gene_35 += 1
    #             elif row["AMR Gene Family"] == "lsa-type ABC-F protein": amr_gene_36 += 1
    #             elif row["AMR Gene Family"] == "macrolide phosphotransferase (MPH)": amr_gene_37 += 1
    #             elif row["AMR Gene Family"] == "major facilitator superfamily (MFS) antibiotic efflux pump": amr_gene_38 += 1
    #             elif row["AMR Gene Family"] == "major facilitator superfamily (MFS) antibiotic efflux pump; resistance-nodulation-cell division (RND) antibiotic efflux pump": amr_gene_39 += 1
    #             elif row["AMR Gene Family"] == "MCR phosphoethanolamine transferase": amr_gene_40 += 1
    #             elif row["AMR Gene Family"] == "Miscellaneous ABC-F subfamily ATP-binding cassette ribosomal protection proteins": amr_gene_41 += 1
    #             elif row["AMR Gene Family"] == "msr-type ABC-F protein": amr_gene_42 += 1
    #             elif row["AMR Gene Family"] == "multidrug and toxic compound extrusion (MATE) transporter": amr_gene_43 += 1
    #             elif row["AMR Gene Family"] == "NDM beta-lactamase": amr_gene_44 += 1
    #             elif row["AMR Gene Family"] == "non-erm 23S ribosomal RNA methyltransferase (G748)": amr_gene_45 += 1
    #             elif row["AMR Gene Family"] == "OXA beta-lactamase": amr_gene_46 += 1
    #             elif row["AMR Gene Family"] == "quinolone resistance protein (qnr)": amr_gene_47 += 1
    #             elif row["AMR Gene Family"] == "resistance-nodulation-cell division (RND) antibiotic efflux pump": amr_gene_48 += 1
    #             elif row["AMR Gene Family"] == "resistance-nodulation-cell division (RND) antibiotic efflux pump; General Bacterial Porin with reduced permeability to beta-lactams": amr_gene_49 += 1
    #             elif row["AMR Gene Family"] == "SIM beta-lactamase": amr_gene_50 += 1
    #             elif row["AMR Gene Family"] == "small multidrug resistance (SMR) antibiotic efflux pump": amr_gene_51 += 1
    #             elif row["AMR Gene Family"] == "streptogramin vat acetyltransferase": amr_gene_52 += 1
    #             elif row["AMR Gene Family"] == "streptogramin vgb lyase": amr_gene_53 += 1
    #             elif row["AMR Gene Family"] == "streptothricin acetyltransferase (SAT)": amr_gene_54 += 1
    #             elif row["AMR Gene Family"] == "TEM beta-lactamase": amr_gene_55 += 1
    #             elif row["AMR Gene Family"] == "tetracycline inactivation enzyme": amr_gene_56 += 1
    #             elif row["AMR Gene Family"] == "tetracycline-resistant ribosomal protection protein": amr_gene_57 += 1
    #             elif row["AMR Gene Family"] == "trimethoprim resistant dihydrofolate reductase dfr": amr_gene_58 += 1
    #             elif row["AMR Gene Family"] == "vanH; glycopeptide resistance gene cluster": amr_gene_59 += 1
    #             elif row["AMR Gene Family"] == "vanS; glycopeptide resistance gene cluster": amr_gene_60 += 1
    #             elif row["AMR Gene Family"] == "vanW; glycopeptide resistance gene cluster": amr_gene_61 += 1
    #             elif row["AMR Gene Family"] == "vanX; glycopeptide resistance gene cluster": amr_gene_62 += 1
    #             elif row["AMR Gene Family"] == "vanY; glycopeptide resistance gene cluster": amr_gene_63 += 1
    #             elif row["AMR Gene Family"] == "vanZ; glycopeptide resistance gene cluster": amr_gene_64 += 1
    #             elif row["AMR Gene Family"] == "vga-type ABC-F protein": amr_gene_65 +=1
    #
    # with open('genomeinfo.tsv', 'r') as file:
    #     reader = csv.DictReader(file, delimiter='\t')
    #     for row in reader:
    #         if f'{name}' in row['Biosampleid']:
    #             submitter = row['Submitter']

    with open('bacteriacount.tsv', 'a') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow([f'{name}',
                             f'{vir_gene_1}', f'{vir_gene_2}', f'{vir_gene_3}', f'{vir_gene_4}', f'{vir_gene_5}',
                             f'{vir_gene_6}', f'{vir_gene_7}', f'{vir_gene_8}', f'{vir_gene_9}', f'{vir_gene_10}',
                             f'{vir_gene_11}', f'{vir_gene_12}', f'{vir_gene_13}', f'{vir_gene_14}', f'{vir_gene_15}',
                             f'{vir_gene_16}', f'{vir_gene_17}', f'{vir_gene_18}', f'{vir_gene_19}', f'{vir_gene_20}',
                             f'{vir_gene_21}', f'{vir_gene_22}', f'{vir_gene_23}', f'{vir_gene_24}', f'{vir_gene_25}',
                             f'{vir_gene_26}', f'{vir_gene_27}', f'{vir_gene_28}', f'{vir_gene_29}', f'{vir_gene_30}',
                             f'{vir_gene_31}', f'{vir_gene_32}', f'{vir_gene_33}', f'{vir_gene_34}', f'{vir_gene_35}',
                             f'{vir_gene_36}', f'{vir_gene_37}', f'{vir_gene_38}', f'{vir_gene_39}', f'{vir_gene_40}',
                             f'{vir_gene_41}', f'{vir_gene_42}', f'{vir_gene_43}', f'{vir_gene_44}', f'{vir_gene_45}',
                             f'{vir_gene_46}', f'{vir_gene_47}', f'{vir_gene_48}', f'{vir_gene_49}', f'{vir_gene_50}',
                             f'{vir_gene_51}', f'{vir_gene_52}', f'{vir_gene_53}', f'{vir_gene_54}', f'{vir_gene_55}',
                             f'{vir_gene_56}', f'{vir_gene_57}', f'{vir_gene_58}', f'{vir_gene_59}', f'{vir_gene_60}',
                             f'{vir_gene_61}', f'{vir_gene_62}', f'{vir_gene_63}', f'{vir_gene_64}', f'{vir_gene_65}',
                             f'{vir_gene_66}'])
