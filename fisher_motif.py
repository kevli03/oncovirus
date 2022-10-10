# use Fisher's exact test to determine significant changes in motif presence
# columns: specific motif count, rest motif count
# rows: result motif total, original motif total

import numpy as np
import pandas as pd
import math

def fisher():
    orig_file = input('Name of original input file: ')
    orig_total = int(input('Total motif count of original: '))
    result_file = input('Name of result input file: ')
    result_total = int(input('Total motif count of result: '))
    output_file = input('Name of output file: ')
    diff_total = orig_total - result_total
    orig_counts = {}; result_counts = {}; diff_counts = {}; roles = {}; prob_counts = {}
    orig = pd.read_csv(orig_file)
    result = pd.read_csv(result_file)
    sig_df = pd.DataFrame(columns=['Motif', 'p-value', 'Role'])

    for index, row in orig.iterrows():
        orig_counts[row['Motif']] = row['Count']
        roles[row['Motif']] = row['Role']

    for index, row in result.iterrows():
        result_counts[row['Motif']] = row['Count']

    for id in orig_counts:
        diff_counts[id] = orig_counts[id] - result_counts[id]

    for id in orig_counts:
        r = math.comb(result_total, result_counts[id])
        d = math.comb(diff_total, diff_counts[id])
        o = math.comb(orig_total, orig_counts[id])
        prob_counts[id] = (r * d) / o

    for id in prob_counts.copy():
        if (prob_counts[id] > 0.05):
            del prob_counts[id]

    for id in prob_counts:
        sig_df = sig_df.append({'Motif':id, 'p-value':prob_counts[id], 'Role':roles[id]}, ignore_index = True)

    sig_df.to_csv(output_file)
    print('# significant:', len(prob_counts))

def main():
    fisher()

main()