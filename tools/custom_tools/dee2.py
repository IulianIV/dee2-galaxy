import sys
from dee2_conn import DEE2


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def main():
    selector, data_set, species = '', '', ''

    try:
        selector = sys.argv[1]
        data_set = sys.argv[2]
        species = sys.argv[3]
        # use this variable to grab a getDEE2 function selected from galaxy
        # function = sys.argv[4]

    except Exception:
        stop_err("Usage: python dee2.py data_set_type species")

    if selector == 'accession':

        dee2 = DEE2(supress_r_warnings=False)

        dee2.data_set = data_set
        dee2.species = species

        # query = dee2.queryDEE2()
        get_dee2 = dee2.getDEE2()
        print(f'final results: {get_dee2}')


if __name__ == '__main__':
    main()

# library('getDEE2')
# lsf.str("package:getDEE2")
# ls("package:getDEE2")
# bundles <- list_bundles("celegans")
# bundles
# browseVignettes("getDEE2")
# # if a species name is wrong, it throws error and examples of species
# mdat <- getDEE2Metadata("celegans")
# head(mdat)
#
# mdat1 <- mdat[which(mdat$SRP_accession %in% "SRP009256"),]
# type(mdat)
# SRRvec <- as.vector(mdat1$SRR_accession)
# SRRvec
#
# suppressPackageStartupMessages(library("SummarizedExperiment"))
# x <- getDEE2("celegans",SRRvec,metadata=mdat,counts="GeneCounts")

