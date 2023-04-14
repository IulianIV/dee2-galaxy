import sys
import argparse
from dee2conn import DEE2

# TODO for each selectable function in dee2.xml add a help attribute to mention what attributes should be set for it
#   to properly work
# TODO create a more comprehensive help section for CLI functions.
# TODO Learn more about Galaxy ObjectStore: https://galaxyproject.org/admin/objectstore/
# TODO Redo xml commands and structure as it has been passed through planemo linter
# TODO Test with planemo

parser = argparse.ArgumentParser(
    description='getDEE2 R package wrapper. Can run getDEE2 functions provided the right arguments')
parser.add_argument('-f', '--function', metavar='', help='name of functions to run.')
parser.add_argument('-t', '--type', metavar='', help='Type of function. Makes processing easier')
parser.add_argument('-s', '--species', metavar='', help='Organism of interest to search for.')
parser.add_argument('-d', '--dataset', metavar='', help='Accession numbers or keyword.')
parser.add_argument('-o', '--outfile', metavar='', help='Output file name.')
parser.add_argument('-i', '--infile', metavar='', help='Infile file name.')
parser.add_argument('-col', '--column', metavar='', help='Column to query.')
parser.add_argument('-b', '--bundle', metavar='', help='Bundles file to use.')
parser.add_argument('-m', '--metadata', metavar='', help='Metadata file to use.')
parser.add_argument('-c', '--counts', metavar='', help='Counts to filter by.')
parser.add_argument('--legacy', action='store_true', help='Use legacy mode.')
parser.add_argument('--set-base-url', metavar='', help='Changes the base URL to fetch data from.')
parser.add_argument('--enable-debugging', action='store_true',
                    help='Enables a limited debugger to see the R console results')

parser.add_argument_group()
args = parser.parse_args()


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def main():
    function = args.function
    type_ = args.type
    species = args.species
    dataset = args.dataset
    outfile = args.outfile
    infile = args.infile
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
        species: {species}
        dataset: {dataset}
        outfile: {outfile}
        infile: {infile}
        column: {column}
        bundle: {bundle}
        metadata: {metadata}
        counts: {counts}
        legacy: {legacy}
        url: {url}
        debugging: {debugging}
        f_type: {type_}
        ''')

    dee2 = DEE2(supress_r_warnings=True)

    dee2.data_set = dataset
    dee2.species = species

    if function == 'getDEE2':

        results = dee2.getDEE2().colData.to_pd()

        csv_file = dee2.convert_to_csv(results, outfile, mode='w', sep='\t')

        return csv_file

    # if type_ == 'loader' and infile is not None:
    #     try:
    #         load_func = getattr(dee2, function)(infile)
    #         return print(load_func)
    #     except AttributeError:
    #         return print("Called load function does not exist.")


if __name__ == '__main__':
    main()
