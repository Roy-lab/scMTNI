import os
import pandas as pd
import argparse

def main(args):
    indir = args.indir + '/'
    outdir = args.outdir + '/'
    regfile = args.regfile
    splitgene = args.splitgene

    files = pd.read_table(args.filelist, header=None)
    species_names = files[0].values
    expression_files = files[1].values

    if len(regfile) > 0:
        with open(regfile) as file:
            regulators = {gene.strip() for gene in file}
    else:
        regulators = None

    # Collect the set of all genes from all expression files
    all_genes = set()
    for expression_file in expression_files:
        data = pd.read_table(expression_file, index_col=0)
        species_genes = set(data.index.astype(str))
        all_genes.update(species_genes)

    # all_regulators is the intersection of regulators and all genes actually in the expression data
    if regulators is not None:
        all_regulators = regulators & all_genes
    else:
        all_regulators = all_genes

    # Write global AllGenes.txt file
    all_genes = sorted(all_genes)
    with open(indir + 'AllGenes.txt', "w") as file:
        file.write("\n".join(all_genes))

    # Write global TFs.txt file.
    all_regulators = sorted(all_regulators)
    with open(indir + 'TFs.txt', "w") as file:
        file.write("\n".join(all_regulators))

    ## split AllGenes.txt into 50 genes per run
    if splitgene:
        splits_dir = indir + 'gene_splits/'
        os.makedirs(splits_dir, exist_ok=True)

        j = 0
        for i in range(0, len(all_genes), splitgene):
            genes = all_genes[i: i + splitgene]
            filename = splits_dir + 'AllGenes' + str(j) + '.txt'
            with open(filename, "w") as file:
                file.write("\n".join(genes))
            j += 1

        print(f"Created {j} gene splits.")

    ## generate config file: testdata_config.txt

    cols = {
        0: species_names,
        1: expression_files,
        2: outdir + species_names,
    }

    if args.motifs:
        cols[3] = indir + species_names+'_network.txt'
        config_filename = indir + 'testdata_config.txt'
    else:
        cols[3] = indir + 'motif_empty.txt'
        config_filename = indir + 'testdata_config_noprior.txt'

    config = pd.DataFrame(cols)
    config.to_csv(config_filename, sep='\t', mode='w', index=None, header=False, float_format='%g')

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--filelist', help='file list', type=str, default='')
    parser.add_argument('--regfile', help='regulator file', type=str, default='')
    parser.add_argument('--indir', help='directory of input data', type=str, default='')
    parser.add_argument('--outdir', help='output directory', type=str, default='Results/')
    parser.add_argument('--motifs', help='add motifs as prior or not', type=bool, default=0)
    parser.add_argument('--splitgene', help='split all genes into X genes per run(0: no split, >0: X)', type=int, default=0)
    args = parser.parse_args()
    main(args)
