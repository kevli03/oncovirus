# identify short linear motifs (SLiMs)

import numpy as np
import pandas as pd
import re

motif_dict = {\
'(([DEST]|^).{0,4}[LI].C.E.{1,4}[FLMIVAWPHY].{0,8}([DEST]|$))':'LIG_Rb_LxCxE_1', \
'.P[TS]AP.':'LIG_PTAP_UEV_1', \
'[ILVMF]YP.[ILVMF]':'LIG_LYPXL_S_1', \
'[FVILMY].FG[DES]F':'LIG_G3BP_FGDF_1', \
'[LMTAFSRI][^KRG]W[DE].{3,5}[LIVMFPA]':'LIG_KLC1_WD_1', \
'[ILVM]YP...[ILVM][^P][^P][LI]':'LIG_LYPXL_L_2', \
'[DE]H.Y':'LIG_HCF-1_HBM_1', \
'([VILPF].{1,3}L.I(S))':'LIG_IRF3_LxIS_1', \
'..[VILM]..[LM][FY]D.':'LIG_Rb_pABgroove_1', \
'....[LIFVYMTE][ASGC][^P]{2}L[^P]{2}[IVMTL][GACS][D][^P][FVLMI].':'LIG_BH_BH3_1',\
'EP[IL]Y[TAG]':'LIG_CSK_EPIYA_1', \
'.KEN.':'DEG_APCC_KENbox_2', \
'([LIVMP].{0,2}(T)P..([ST]))':'DEG_SCF_FBW7_1', \
'(D(S)G.{2,3}([ST]))':'DEG_SCF_TrCP1_1', \
'.P.A.V.P[^P]':'DEG_SIAH_1', \
'.R..[PGAV][DEIP]G.':'DOC_ANK_TNKS_1', \
'R.LF':'DOC_CYCLIN_1', \
'RV.F':'DOC_PP1', \
'SP.L.LT':'DOC_PP2B_1', \
'P.E[^P].S[^P]':'DOC_USP7_MATH_2', \
'(R[^DE]{0,2}[^DEPG]([ST])(([FWYLMV].)|([^PRIKGN]P)|([^PRIKGN].{2,4}[VILMFWYP])))':'LIG_14-3-3_CanoR_1', \
'L[IVLMF].[IVLMF][DE]':'LIG_Clathr_ClatBox_1', \
'L[^P]{2,2}[HI]I[^P]{2,2}[IAV][IL]':'LIG_CORNRBOX', \
'((P[LVIPME][DENS][LM][VASTRG])|(G[LVIPME][DENS][LM][VASTRG]((K)|(.[KR]))))':'LIG_CtBP_PxDLS_1', \
'[^P].[KR].TQT':'LIG_Dynein_DLC8_1', \
'P.L.P':'LIG_MYND_1', \
'...[ST].[ACVILF]$':'LIG_PDZ_Class_1', \
'RGD':'LIG_RGD', \
'Y..V':'LIG_SH2', \
'P..P.[KR]':'LIG_SH3_2', \
'[PSAT].[QE]E':'LIG_TRAF2_1', \
'..P.E..[FYWHDE].':'LIG_TRAF6', \
'([DEN]..(Y)..[LI].{6,12}(Y)..[LI])':'LIG_TYR_ITAM', \
'((C)[^DENQ][LIVMF].$)':'MOD_CAAXbox', \
'(.(N)[^P][ST]..)':'MOD_N-GLC_1', \
'(^M{0,1}(G)[^EDRKHPFYW]..[STAGCN][^P])':'MOD_Nmyristoyl', \
'([VILMAFP](K).E)':'MOD_SUMO', \
'Y..[LMVIF]':'TRG_ENDOCYTIC_2', \
'K.{0,1}K.{2,3}$':'TRG_ER_diLys_1', \
'[KRHQSAP][DENQT]EL$':'TRG_ER_KDEL_1', \
'E...LL':'TRG_LysEnd_APsAcLL_1', \
'(([DEQ].{0,1}[LIM].{2,3}[LIVMF][^P]{2,3}[LMVF].[LMIV].{0,3}[DE])|([DE].{0,1}[LIM].{2,3}[LIVMF][^P]{2,3}[LMVF].[LMIV].{0,3}[DEQ]))':'TRG_NES_CRM1_1'}

role_dict = {\
'LIG_Rb_LxCxE_1':'binds retinoblastoma B pocket', \
'LIG_PTAP_UEV_1':'UEV domain binding PTAP motif', \
'LIG_LYPXL_S_1':'endosomal sorting of membrane proteins', \
'LIG_G3BP_FGDF_1':'binds Ras GTPase activating SH3 domain', \
'LIG_KLC1_WD_1':'binds kinesin light chain TPR region', \
'LIG_LYPXL_L_2':'endosomal sorting of membrane proteins', \
'LIG_HCF-1_HBM_1':'binds transcriptional coactivator HCF-1', \
'LIG_IRF3_LxIS_1':'interferon regulatory factor 3 binding site', \
'LIG_Rb_pABgroove_1':'binds retinoblastoma AB groove', \
'LIG_BH_BH3_1':'binds BH domains to inhibit apoptosis',\
'LIG_CSK_EPIYA_1':'binds C-terminal Src kinase SH2 domain', \
'DEG_APCC_KENbox_2':'binds to subunit, causing protein to be targeted for degradation', \
'DEG_SCF_FBW7_1':'binds box proteins of SCF complex', \
'DEG_SCF_TrCP1_1':'binds F box protein of SCF complex', \
'DEG_SIAH_1':'binds to substrate-binding domain of SIAH family', \
'DOC_ANK_TNKS_1':'interacts with ankyrin repeat domain region to facilitate PARsylation', \
'DOC_CYCLIN_1':'N/A', \
'DOC_PP1':'involved in docking', \
'DOC_PP2B_1':'N/A', \
'DOC_USP7_MATH_2':'MATH domain binding variant', \
'LIG_14-3-3_CanoR_1':'phospho-motif mediating strong interaction with 14-3-3 proteins', \
'LIG_Clathr_ClatBox_1':'interacts with beta propeller structure on Clathrin heavy chain', \
'LIG_CORNRBOX':'confers binding to nuclear receptors', \
'LIG_CtBP_PxDLS_1':'interacts with NAD-dependent repressor CtBP proteins', \
'LIG_Dynein_DLC8_1':'interacts with target-accepting grooves of Dynein light chain dimer', \
'LIG_MYND_1':'recognized by MYND domain-containing proteins', \
'LIG_PDZ_Class_1':'C-terminal class 1 PDZ-binding motif', \
'LIG_RGD':'N/A', \
'LIG_SH2':'SH2 domain binding motif', \
'LIG_SH3_2':'recognized by class II SH3 domains', \
'LIG_TRAF2_1':'major TRAF2-binding consensus motif', \
'LIG_TRAF6':'TRAF6 binding site', \
'LIG_TYR_ITAM':'immunoreceptor tyrosine-based activatory motif', \
'MOD_CAAXbox':'generic CAAX box prenylation motif', \
'MOD_N-GLC_1':'generic motif for N-glycosylation', \
'MOD_Nmyristoyl':'generic motif for N-Myristoylation site', \
'MOD_SUMO':'motif recognized for modification by SUMO-1', \
'TRG_ENDOCYTIC_2':'tyrosine-based sorting signal, interacts with subunit of AP complex', \
'TRG_ER_diLys_1':'responsible for COPI-mediated retrieval from post-ER compartments', \
'TRG_ER_KDEL_1':'Golgi-to-ER retrieving signal', \
'TRG_LysEnd_APsAcLL_1':'N/A', \
'TRG_NES_CRM1_1':'leucine-rich nuclear export signal binding to CRM1 exportin protein'}

def identify(todo):
    input_file = input('Name of input file: ')
    type_dict = {}
    for regex in motif_dict:
        type_dict[motif_dict[regex]] = 0

    motif = 0; total = 0
    data = pd.read_csv(input_file)
    extra_cols = ['has motif', 'motif seq', 'motif', 'role']
    motif_df = pd.DataFrame(columns=(list(data) + extra_cols))          #first column has data
    # motif_df = pd.DataFrame(columns=(list(data)[1:] + extra_cols))    #else

    if (todo == 'N'):
        motif_output_file = input('Name of output motif file: ')
        stats_output_file = input('Name of output statistical file: ')
        for index, row in data.iterrows():
            seq = row['protein']
            motifs = {}

            for regex in motif_dict:
                match = re.search(regex, seq)
                if match:
                    sequence = re.findall(regex, seq)
                    for i in range(len(sequence)):
                        if (type(sequence[i]) is tuple):
                            sequence[i] = sequence[i][0]
                        type_dict[motif_dict[regex]] += 1
                    motifs[motif_dict[regex]] = sequence
            
            output_dict = {}
            for col in (list(data)):            #first column has data
            # for col in (list(data)[1:]):      #else
                output_dict[col] = row[col]
            
            if (len(motifs) > 0):
                motif += 1
                for i in motifs:
                    for j in motifs[i]:
                        motif_df = motif_df.append((output_dict | {'has motif':'Yes', 'motif seq':j, 'motif':i, 'role':role_dict[i]}), ignore_index=True)
                        # replace append with concat (deprecated)
            else:
                motif_df = motif_df.append((output_dict | {'has motif':'No', 'motif seq':'none', 'motif':'none', 'role':'none'}), ignore_index=True)
            total += 1

        motif_df.index = np.arange(1, len(motif_df) + 1)
        motif_df.to_csv(motif_output_file)
        total_sum = sum(list(type_dict.values()))
        df2 = pd.DataFrame([type_dict]); df2 = df2.transpose()
        df3 = pd.DataFrame([role_dict]); df3 = df3.transpose()
        stats_df = pd.concat([df2, df3], axis = 1)
        stats_df.to_csv(stats_output_file)
        print()
        print(f'Found {total_sum} total motifs out of {motif} peptides with motifs out of {total} total peptides.')

    if (todo == 'F'):        
        prot_output_file = input('Name of protein output statistical file: ')
        pep_output_file = input('Name of peptide output statistical file: ')
        prot_count_dict = {}; pep_count_dict = {}
        for index, row in data.iterrows():
            seq = row['protein']

            for regex in motif_dict:
                match = re.search(regex, seq)
                if match:
                    sequence = re.findall(regex, seq)
                    for i in range(len(sequence)):
                        if (str(row['Virus']) + ' - ' + str(row['gene_name'])) not in prot_count_dict:
                            prot_count_dict[str(row['Virus']) + ' - ' + str(row['gene_name'])] = 1
                        else:
                            prot_count_dict[str(row['Virus']) + ' - ' + str(row['gene_name'])] += 1
                        if (str(row['Virus']) + ' - ' + str(row['gene_name']) + ' - ' + str(row['protein'])) not in pep_count_dict:
                            pep_count_dict[str(row['Virus']) + ' - ' + str(row['gene_name']) + ' - ' + str(row['protein'])] = 1
                        else:
                            pep_count_dict[str(row['Virus']) + ' - ' + str(row['gene_name']) + ' - ' + str(row['protein'])] += 1

        sum_prot = sum(list(prot_count_dict.values()))
        dfpr = pd.DataFrame([prot_count_dict]); dfpr = dfpr.transpose(); dfpr.to_csv(prot_output_file)
        dfpe = pd.DataFrame([pep_count_dict]); dfpe = dfpe.transpose(); dfpe.to_csv(pep_output_file)
        print()
        print(f'Found {sum_prot} motifs.')

    if (todo == 'A'):
        output_file = input('Name of output statistical file: ')
        for index, row in data.iterrows():
            seq = row['protein']

            for regex in motif_dict:
                match = re.search(regex, seq)
                if match:
                    type_dict[motif_dict[regex]] += 1

        total = sum(list(type_dict.values()))
        df2 = pd.DataFrame([type_dict]); df2 = df2.transpose()
        df3 = pd.DataFrame([role_dict]); df3 = df3.transpose()
        stats_df = pd.concat([df2, df3], axis = 1)
        stats_df.to_csv(output_file)
        print()
        print(f'Found {total} motifs.')

def main():
    todo = input('normal or peptide/protein full count or abbreviated count? (N/F/A): ')
    identify(todo)

main()
