import sys
import argparse
from dee2conn import DEE2
from dee2converter import SummarizedExperiment, ConvertedMatrix, ConvertedListVector

# TODO create a more comprehensive help section for CLI functions.
# TODO Learn more about Galaxy ObjectStore: https://galaxyproject.org/admin/objectstore/
# TODO Redo xml commands and structure as it has been passed through planemo linter
# TODO Test with planemo
# TODO add a function that parses arguments and sets them as dee2 class attributes.
# TODO add error handling. What happens when there is a NullType because of wrong accessions?
#   what happens when the error that some accessions are wrong is risen? (only happens in R)

parser = argparse.ArgumentParser(
    description='getDEE2 R package wrapper. Can run getDEE2 functions provided the right arguments')
parser.add_argument('-f', '--function', metavar='', help='name of function to run.')
parser.add_argument('-of', '--other-functions', nargs='+', metavar='', help='name of function to run.')
parser.add_argument('-s', '--species', metavar='', help='Organism of interest to search for.')
parser.add_argument('-d', '--dataset', metavar='', help='Accession numbers or keyword.')
parser.add_argument('-o', '--outfile', metavar='', help='Output file name.')

# leaving the infile args as possible future implementation outside of galaxy. It is redundant at the moment.
# parser.add_argument('-i', '--infile', metavar='', help='Infile file name.')

parser.add_argument('-col', '--column', metavar='', help='Column to query.')
parser.add_argument('-b', '--bundle', metavar='', help='Bundles file to use.')
parser.add_argument('-m', '--metadata', metavar='', help='Metadata file to use.')
parser.add_argument('-c', '--counts', metavar='', help='Counts to filter by.')
parser.add_argument('--legacy', metavar='', help='Use legacy mode.')
parser.add_argument('--set-base-url', metavar='', help='Changes the base URL to fetch data from.')
parser.add_argument('--enable-debugging', metavar='',
                    help='Enables a limited debugger to see the R console results')

parser.add_argument_group()
args = parser.parse_args()


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def main():
    results = None
    function = args.function
    other_functions = args.other_functions
    species = args.species
    dataset = args.dataset
    outfile = args.outfile
    column = args.column
    bundle = args.bundle
    metadata = args.metadata
    counts = args.counts
    legacy = args.legacy
    url = args.set_base_url
    debugging = args.enable_debugging

    print(f'''
        All Args:
        function: {function}
        other_functions: {other_functions}
        species: {species}
        dataset: {dataset}
        outfile: {outfile}
        column: {column}
        bundle: {bundle}
        metadata: {metadata}
        counts: {counts}
        legacy: {legacy}
        url: {url}
        debugging: {debugging}
        ''')

    dee2 = DEE2(supress_r_warnings=True)

    dee2.data_set = dataset
    dee2.species = species

    if function == 'getDEE2_bundle' and len(dataset.split(',')) > 1:
        raise ValueError(f'Function {function} only accepts one Accession Number.'
                         f' Value {",".join(dataset.split(","))} is not valid')

    if 'load' in function:
        func_call = getattr(dee2, function)(metadata)
    else:
        func_call = getattr(dee2, function)()

    if isinstance(func_call, SummarizedExperiment):
        # this is identical to the R: assays(x)$counts
        # This needs to be implemented somehow
        # results = func_call.assays.data.listData.counts

        results = func_call.colData.to_pd()
        results = dee2.convert_to_csv(results, outfile)
        print(results)
    elif isinstance(func_call, ConvertedMatrix):
        results = func_call.df_matrix
        results = dee2.convert_to_csv(results, outfile)
    elif isinstance(func_call, ConvertedListVector):
        results = func_call.MetadataFull.to_pd()
        results = dee2.convert_to_csv(results, outfile)
    else:
        results = dee2.convert_to_csv(func_call, outfile)

    return results


if __name__ == '__main__':
    main()
