# identify peptides based on "significance" or type
# *significant: direction = 2, |LFC| > 1, q_value < 0.001
# find specific or overlapping peptides

import numpy as np
import pandas as pd

def identify_type(pep_type):
    input_file = input('Name of input file: ')
    output_file = input('Name of output file: ')

    type_num = 0; total = 0
    data = pd.read_csv(input_file)
    sig_df = pd.DataFrame(columns=(list(data)[1:]))

    for index, row in data.iterrows():
        if (row['Type2'] == pep_type):
            type_num += 1
            output_dict = {}
            for col in (list(data)[1:]):
                output_dict[col] = row[col]
            sig_df = sig_df.append(output_dict, ignore_index=True)
        total += 1

    sig_df.index = np.arange(1, len(sig_df) + 1)
    sig_df.to_csv(output_file)
    print()
    print(f'Found {type_num} {pep_type} peptides out of {total} total peptides.')

def identify_sig(direction, LFC, q_value):
    input_file = input('Name of input file: ')
    output_file = input('Name of output file: ')

    sig = 0; total = 0
    data = pd.read_csv(input_file)
    sig_df = pd.DataFrame(columns=(list(data)[1:]))

    for index, row in data.iterrows():
        if (row['Direction'] == direction and abs(float(row['LFC1'])) > LFC and abs(float(row['LFC2'])) > LFC and \
            abs(float(row['LFC3'])) > LFC and float(row['Q_value_1']) < q_value and \
            float(row['Q_value_2']) < q_value and float(row['Q_value_3'] < q_value)):
            sig += 1
            output_dict = {}
            for col in (list(data)[1:]):
                output_dict[col] = row[col]
            sig_df = sig_df.append(output_dict, ignore_index=True)
        total += 1

    sig_df.index = np.arange(1, len(sig_df) + 1)
    sig_df.to_csv(output_file)
    print()
    print(f'Found {sig} significant peptides out of {total} total peptides.')

def identify_enr_dep(task):
    input_file = input('Name of input file: ')
    output_file = input('Name of output file: ')

    sig = 0; total = 0
    data = pd.read_csv(input_file)
    sig_df = pd.DataFrame(columns=(list(data)[1:]))

    if (task == 'E'):
        for index, row in data.iterrows():
            if (float(row['LFC1']) > 0 and float(row['LFC2']) > 0 and float(row['LFC3']) > 0):
                sig += 1
                output_dict = {}
                for col in (list(data)[1:]):
                    output_dict[col] = row[col]
                sig_df = sig_df.append(output_dict, ignore_index=True)
            total += 1
        print(f'Found {sig} enriched peptides out of {total} total peptides.')

    elif(task == 'D'):
        for index, row in data.iterrows():
            if (float(row['LFC1']) < 0 and float(row['LFC2']) < 0 and float(row['LFC3']) < 0):
                sig += 1
                output_dict = {}
                for col in (list(data)[1:]):
                    output_dict[col] = row[col]
                sig_df = sig_df.append(output_dict, ignore_index=True)
            total += 1
        print(f'Found {sig} depleted peptides out of {total} total peptides.')

    sig_df.to_csv(output_file)

def find_spec(file):
    spec_file = open(file, 'r')
    specs = []
    for line in spec_file:
        specs.append(line.strip())
    
    input_file = input('Name of input file: ')
    output_file = input('Name of output file: ')

    count = 0
    data = pd.read_csv(input_file)
    sig_df = pd.DataFrame(columns=(list(data)[1:]))

    for index, row in data.iterrows():
        if (row['ID'] in specs):
        #if (row['motif'] in specs):
            count += 1
            output_dict = {}
            for col in (list(data)[1:]):
                output_dict[col] = row[col]
            sig_df = sig_df.append(output_dict, ignore_index=True)

    sig_df.index = np.arange(1, len(sig_df) + 1)
    sig_df.to_csv(output_file)
    print()
    print(f'Found {count} occurences.')

def find_over():
    task = input('enriched, depleted, or overall? (input E/D/O): ')
    if (task not in 'EDO'):
        print('Invalid input')

    input_files = []
    while (True):
        next = input('Name of an input file (preferably smallest first, N if done listing): ')
        if (next == 'N'):
            break
        else:
            input_files.append(next)        
    output_file = input('Name of output file: ')

    data = []; file_num = len(input_files)
    for file in input_files:
        data.append(pd.read_csv(file))
    sig_df = pd.DataFrame(columns=(list(data[0])[1:]))

    count_dict = {}; prot_list = []
    for index, row in data[0].iterrows():
        count_dict[row['ID']] = 1

    for i in range(1, file_num):
        for index, row in data[i].iterrows():
            if (row['ID'] in count_dict):
                count_dict[row['ID']] += 1
    
    for id in count_dict.copy():
        if (count_dict[id] != file_num):
            del count_dict[id]

    LFC_dict = {}
    if (task != 'O'):
        for index, row in data[0].iterrows():
            if (row['ID'] in count_dict):
                LFC_dict[row['ID']] = row['LFC']

        for i in range(1, file_num):
            for index, row in data[i].iterrows():
                if (row['ID'] in count_dict):
                    if (LFC_dict[row['ID']] * row['LFC'] <= 0):
                        LFC_dict[row['ID']] = 0
    
        for id in count_dict.copy():
            if (LFC_dict[id] == 0):
                del count_dict[id]

        factor = 0
        if (task == 'E'):
            factor = 1
        elif (task == 'D'):
            factor = -1

        for id in count_dict.copy():
            if (LFC_dict[id] * factor <= 0):
                del count_dict[id]

    for index, row in data[0].iterrows():
        if (row['ID'] in count_dict):
            if (row['Gene'] not in prot_list):
                prot_list.append(row['Gene'])
            output_dict = {}
            for col in (list(data[0])[1:]):
                output_dict[col] = row[col]
            sig_df = sig_df.append(output_dict, ignore_index=True)

    sig_df.index = np.arange(1, len(sig_df) + 1)
    sig_df.to_csv(output_file)
    print()
    print(f'Found {len(count_dict)} peptides and {len(prot_list)} proteins.')

def main():
    todo = input('identify or find? (input I or F): ')
    if (todo == 'I'):
        task = input('type, significance, or enriched/depleted? (input T/S/E/D): ')
        if (task == 'T'):
            pep_type = input('name of type?: ')
            identify_type(pep_type)
        elif (task == 'S'):
            direction = (input('direction?: '))
            LFC = float(input('|LFC| > ?: '))
            q_value = float(input('q_value < ?: '))
            identify_sig(direction, LFC, q_value)
        elif (task in 'ED'):
            identify_enr_dep(task)
        else:
            print('Invalid input')
    elif (todo == "F"):
        task = input('specific or overlapping? (input S or O): ')
        if (task == 'S'):
            spec = input('input name of .txt file containing identifiers: ')
            find_spec(spec)
        elif (task == 'O'):
            find_over()
        else:
            print('Invalid input')
    else:
        print('Invalid input')

main()
