import sys
import argparse
from dee2conn import DEE2
from dee2converter import SummarizedExperiment, ConvertedMatrix, ConvertedListVector

# TODO create a more comprehensive help section for CLI functions.

parser = argparse.ArgumentParser(
    description='getDEE2 R package wrapper. Can run getDEE2 functions provided the right arguments')
parser.add_argument('-f', '--function', metavar='', help='name of function to run.')
parser.add_argument('-s', '--species', metavar='', help='Organism of interest to search for.')
parser.add_argument('-d', '--dataset', metavar='', help='Accession numbers or keyword.')
parser.add_argument('-o', '--outfile', metavar='', help='Output file name.')
parser.add_argument('-col', '--column', metavar='', help='Column to query.')
parser.add_argument('-b', '--bundle', metavar='', help='Bundles file to use.')
parser.add_argument('-z', '--dee2-zip', metavar='', help='DEE2 zip file to use.')
parser.add_argument('-m', '--metadata', metavar='', help='Metadata file to use.')
parser.add_argument('-c', '--counts', metavar='', help='Counts to filter by.')
parser.add_argument('--enable-debugging', metavar='',
                    help='Enables a limited debugger to see the R console results')

parser.add_argument_group()
args = parser.parse_args()


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def main():
    function = args.function
    species = args.species
    dataset = args.dataset
    outfile = args.outfile

    dee2 = DEE2()

    dee2.data_set = dataset
    dee2.species = species

    if args.dee2_zip != 'None':
        dee2.dee2_zip = args.dee2_zip

    if args.metadata != 'None':
        dee2.metadata = args.metadata

    if args.column != 'None':
        dee2.col = args.column

    if args.bundle != 'None':
        dee2.bundles = args.bundle

    if args.counts != 'None':
        dee2.counts = args.counts

    if function == 'getDEE2_bundle' and len(dataset.split(',')) > 1:
        raise ValueError(f'Function {function} only accepts one Accession Number.'
                         f' Value {",".join(dataset.split(","))} is not valid')

    func_call = getattr(dee2, function)()

    if isinstance(func_call, SummarizedExperiment):
        # this is identical to the R: assays(x)$counts
        # This needs to be implemented somehow
        # results = func_call.assays.data.listData.counts

        results = func_call.colData.to_pd()
        results = dee2.convert_to_tsv(results, outfile)
    elif isinstance(func_call, ConvertedMatrix):
        results = func_call.df_matrix
        results = dee2.convert_to_tsv(results, outfile)
    elif isinstance(func_call, ConvertedListVector):
        results = func_call.MetadataFull.to_pd()
        results = dee2.convert_to_tsv(results, outfile)
    else:
        results = dee2.convert_to_tsv(func_call, outfile)

    return results


if __name__ == '__main__':
    main()
