"""
    app.py - Data analysis project. This is the first version of a project that
    determines the top 30 gene expression mean values in a set of research
    study subjects with a presumptive diagnosis of Alzheimer's disease (AD).
    Provides option to compute top 30 values with study subjects without AD.

    Rexford Cristal
    Last Modified: 04/21/2019
"""

import flask
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

app = flask.Flask(__name__)

DONORS = pd.read_csv('data/DonorInformation.csv')
SAMPLES = pd.read_csv('data/columns-samples.csv')
GENE_LIST = pd.read_csv('data/rows-genes.csv')
FPKM = pd.read_csv('data/fpkm_table_normalized.csv')
DS = pd.merge(DONORS, SAMPLES, on='donor_id')

def compute(ad):
    """ Computes top 30 gene expression mean values in a set of research study subjects with presumptive AD. """

    if ad == 'True':
        ds_ad = DS[ (DS['nincds_arda_diagnosis'] == "Probable Alzheimer'S Disease") | (DS['nincds_arda_diagnosis'] == "Possible Alzheimer'S Disease") ]
    else:
        ds_ad = DS[ (DS['nincds_arda_diagnosis'] != "Probable Alzheimer'S Disease") & (DS['nincds_arda_diagnosis'] != "Possible Alzheimer'S Disease") ]

    rnaseq = list(ds_ad['rnaseq_profile_id'].values)
    rnaseq = [ str(e) for e in rnaseq ]
    columns = ['gene_id \ rnaseq_profile_id'] + rnaseq

    gene_expression = FPKM.loc[:, columns]
    gene_expression2 = gene_expression.set_index('gene_id \ rnaseq_profile_id')
    gene_expression2['mean'] = gene_expression2.mean(axis=1)
    gene_expression3 = gene_expression2.sort_values('mean', ascending=False)

    gene_results = pd.merge(gene_expression3, GENE_LIST, left_on='gene_id \ rnaseq_profile_id', right_on='gene_id')
    gene_results_top_30 = gene_results.loc[0:29, ('mean', 'gene_symbol')]

    gene_symbols = gene_results_top_30['gene_symbol'].values
    gene_mean_values = gene_results_top_30['mean'].values
    ind = np.arange(len(gene_symbols))
    plt.bar(ind, gene_mean_values)
    plt.xticks(ind, gene_symbols, rotation='vertical')
    plt.xlabel('Gene Symbols')
    plt.ylabel('Mean Gene Expression (FPKM)')
    plt.title('ad={}'.format(ad))
    plt.margins(0.2)
    plt.subplots_adjust(left=0.3, bottom=0.3)
    plt.savefig('static/images/result.png')

    # Print results for debugging.
    result_string = 'Gene Symbol: Mean\n'
    index = 0
    for gene in gene_symbols:
        result_string = result_string + '{}: {}\n'.format(gene, gene_mean_values[index])
        index = index + 1
    print(result_string)

@app.route('/')
def display_input_form():
    return flask.render_template('input_form.html')

@app.route('/result')
def display_result():
    compute(flask.request.args.get('ad'))
    return flask.render_template('result.html', image_path='images/result.png')

@app.after_request
def set_response_headers(response):
    """ Prevent browser from caching image. """

    response.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
    return response

if __name__ == '__main__':
    app.run(debug=True, port=8000)
