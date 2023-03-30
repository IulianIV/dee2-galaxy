import rpy2.rinterface
from rpy2.robjects.packages import importr
from rpy2 import robjects
from rpy2.robjects.vectors import ListVector, BoolVector, DataFrame
from rpy2.robjects.methods import RS4
from rpy2.rinterface import NULL


# TODO All "loadXYZ" functions rely on zipfile processing.
#   could all these functions be combined in a single larger functions for easier access?
#   How should the zip file be processed if there is a choice to not save it locally?
# TODO create a decorator for the species validation?
# TODO Add type hints too all functions and arguments if necessary
# TODO check this situation that occurs on multiple functions
#  @param ... Additional parameters to be passed to download.file.
# TODO It seems to be specific on situations where it requires SRRVector
#   as far as I can tell, when passing a SRR Vector it cant be passed with a single value.
#   Nonetheless, it seems that there are places where "query" is required, which is similar to SRRVector
#   but I think single values can be placed.
# TODO test all `load` functions against a downloaded zip file
# TODO when the users passes a path to to output file check if that path exists, otherwise create it. Add this info
#   as tip at `.xml` file level
# TODO adding out files from the input does not seem to save any files... might be because of Galaxy itself
#   try and see if the file can be saved in history then downloaded

# Set of valid species recognized by dee2.io
valid_species = {'athaliana', 'celegans', 'dmelanogaster', 'drerio', 'ecoli',
                 'hsapiens', 'mmusculus', 'rnorvegicus', 'scerevisiae'}

# Set of valid column names recognized by dee2.io
valid_cols = {'file_name', 'date_added', 'time_added', 'file_size',
              'SRP_accession', 'GSE_accession'}

# Set of valid counts names recognized by dee2.io
valid_counts = {'GeneCounts', 'TxCounts', 'Tx2Gene'}

FALSE = BoolVector((False,))
TRUE = BoolVector((True,))


class DEE2:

    def __init__(self, supress_r_warnings: bool = True):
        """
        Initializes the connection to the getDEE2 R package

        :param supress_r_warnings: Completely blocks warning prompts to the terminal.
        """
        if supress_r_warnings:
            self._supress_r_warnings()

        # R library imports
        self.r_base = importr('base')
        self.dee2_conn = importr('getDEE2')

        self.species = None
        self.data_set = None

        self.metadata = NULL
        self.outfile = NULL

    @staticmethod
    def _supress_r_warnings():
        """
        This basically completely suppresses warning outputs to R console.
        In some cases the R installation, therefore the R's C API does not add packages in all folders, which
        results in empty package folders. Those folders are consistently reported as warnings by R
        and pushed to console.
        """
        import rpy2.rinterface_lib.callbacks

        def no_callback(*args):
            return None

        rpy2.rinterface_lib.callbacks.consolewrite_warnerror = no_callback

    def _check_species(self, species: str) -> bool:
        return self._check_values(species, valid_species)

    def _check_cols(self, col: str) -> bool:
        return self._check_values(col, valid_cols)

    def _check_counts(self, count: str) -> bool:
        return self._check_values(count, valid_counts)

    @staticmethod
    def _check_values(data: str, data_set: set) -> bool:
        """
        Provides a way to validate function arguments against known values. Helps keep DRY.
        :param data: value to check against 'data_set'
        :param data_set: Set datatype which contains valid values
        :return: True or False depending on existence of data in data_set
        """
        if data not in data_set:
            print(f'''
            Provided value {data} is not found in the list. Check spelling and try again.
            Valid values are {', '.join(data_set)}.
            ''')

            return False
        else:
            return True

    def convert_to_srr_vector(self, data: [str, None] = None) -> ListVector:
        splitter = robjects.r['strsplit']
        combine = robjects.r['c']

        data = self.data_set or data

        srr_split = splitter(data, split=',')
        srr_vector = combine(srr_split)

        # because we pass a list, the vector itself becomes part of a list, such as [[1]]
        # to grab it correctly the first index must be accessed
        return srr_vector[0]

    def getDEE2(self, species: str = None, srr_vector: [str, ListVector] = None, counts: str = 'GeneCounts',
                metadata=NULL, outfile: [str, NULL] = NULL,
                legacy=FALSE, base_url='http://dee2.io/cgi-bin/request.sh?') -> [None, RS4]:
        """
        Runs the getDEE2() R function from getDEE2 R package.

        The getDEE2 function fetches gene expression data from the DEE2 database of RNA sequencing data and returns it
        as a SummarizedExperiment object.

        :param counts: A string, either 'GeneCounts', 'TxCounts' or 'Tx2Gene'.
            When 'GeneCounts' is specified, STAR gene level counts are returned.
            When 'TxCounts' is specified, kallisto transcript counts are returned.
            When 'Tx2Gene' is specified, kallisto counts aggregated (by sum) on gene
            are returned. If left blank, "GeneCounts" will be fetched.
        :param species: A character string matching a species of interest.
        :param srr_vector: A character string or vector thereof of SRA run accession numbers
        :param metadata: metadata optional R object of DEE2 metadata to query.
        :param outfile: An optional file name for the downloaded dataset.
        :param legacy: Whether data should be returned in the legacy (list) format. Default is FALSE. Leave this FALSE
            if you want to receive data as Summarized experiment.
        :param base_url: The base URL of the service. Leave this as the default URL unless you want to download from
            a 3rd party mirror.
        :return: a SummarizedExperiment object.
        """

        species = self.species or species
        metadata = self.metadata or metadata
        outfile = self.outfile or outfile

        # TODO add instance level counts
        if not self._check_counts(counts):
            return None

        if not self._check_species(species):
            return None

        get_dee2 = robjects.r['getDEE2']
        srr_vector = self.convert_to_srr_vector() or srr_vector
        data = get_dee2(species, srr_vector, counts, metadata, outfile, legacy, base_url)

        return data

    # TODO bundles is expected to be a table how to pass this to R from Python if previously generated?
    # TODO Check which counts can be added
    # TODO col needs to be selected from Galaxy and passed on. Adding it as class attribute does not make sense
    #   in this case
    # TODO col will have to be replaced with None after the previous TODO is done
    # TODO the 'query' (or in our case the srr_vector) field IS NOT a vector field. It is a standard string
    #   representing a dataset name. Therefore it only accepts A SINGLE VALUE.
    #   Converting it to srr_vector is redundant.
    # TODO edit out the srr_vector, add and validate it as single string and in the front-end, when this function
    #   is chosen show a text saying that from the data_set input ONLY THE FIRST value will be kept and that it only
    #   accepts single values. At the moment tests will be done as is.
    def getDEE2_bundle(self, species: str = None, srr_vector: [str, ListVector] = None, col: str = 'SRP_accession',
                       bundles=NULL, counts: str = 'GeneCounts', legacy=FALSE,
                       base_url: str = "http://dee2.io/huge/") -> [None, RS4]:

        """
        Runs the getDEE2_bundle() R function from getDEE2 R package.

        Get a DEE2 project bundle.

        The getDEE2_bundle function fetches gene expression data from DEE2.

        This function will only work if all SRA runs have been successfully
        processed for an SRA project. This function returns a
        SummarizedExperiment object.

        :param species: A character string matching the species of interest.
        :param srr_vector: A character string, such as the SRA project accession number
            or the GEO series accession number.
        :param col: the column name to be queried, usually "SRP_accession" for SRA
            project accession or "GSE_accession" for GEO series accession.
        :param counts: A string, either 'GeneCounts', 'TxCounts' or 'Tx2Gene'.
            When 'GeneCounts' is specified, STAR gene level counts are returned.
            When 'TxCounts' is specified, kallisto transcript counts are returned.
            When 'Tx2Gene' is specified, kallisto counts aggregated (by sum) on gene
            are returned. If left blank, "GeneCounts" will be fetched.
        :param bundles: optional table of previously downloaded bundles.
            providing this will speed up performance if multiple queries are made in a
            session. If left blank, the bundle list will be fetched again.
        :param legacy: Whether data should be returned in the legacy (list) format.
            Default is FALSE. Leave this FALSE if you want to receive data as Summarized experiment.
        :param base_url: The base URL of the service. Leave this as the default URL
            unless you want to download from a 3rd party mirror.
        :returns: a SummarizedExperiment object.
        """

        species = self.species or species

        if not self._check_species(species):
            return None

        if not self._check_cols(col):
            return None

        # TODO add instance level counts
        if not self._check_counts(counts):
            return None

        get_dee2_bundles = robjects.r['getDEE2_bundle']

        srr_vector = self.convert_to_srr_vector() or srr_vector

        data = get_dee2_bundles(species, srr_vector, col, counts, bundles, legacy, base_url)

        return data

    def getDEE2Metadata(self, species: str = None, outfile: [str, NULL] = NULL) -> [None, DataFrame]:
        """
        Runs the getDEE2Metadata() R function for getDEE2 R package.

        Get DEE2 Metadata.

        This function fetches the short metadata for the species of interest.

        :param species: A character string matching a species of interest.
        :param outfile: Optional filename
        :returns: a table of metadata.
        """
        species = self.species or species
        outfile = self.outfile or outfile

        if not self._check_species(species):
            return None

        get_dee2_metadata = robjects.r['getDEE2Metadata']

        data = get_dee2_metadata(species, outfile)

        return data

    # TODO what to do here front-end wise when no file is provided? just print it?
    def list_bundles(self, species: str = None):
        """
        Runs the list_bundles() R function for getDEE2 R package.

        Get a table of all completed projects at DEE2.

        :param species: A character string matching a species of interest.
        :returns: a table of project bundles available at DEE2.io/huge
        """
        species = self.species or species

        if not self._check_species(species):
            return None

        bundles = robjects.r['list_bundles']

        data = bundles(species)

        return data

    @staticmethod
    def _call_load_functions(func_name: str, zip_name: str = None) -> [str, DataFrame]:
        """
        Makes sure the code respects the DRY principle.

        All getDEE2 R functions which contain a leading 'load' require the same argument 'zipname'
        It makes sense to create a function which handles the calls, instead of duplicating code.
        """

        if func_name is None:
            return print("You have not set a load function to call. Please set the 'load_func' attribute.")

        grab_func = robjects.r[func_name]

        data = grab_func(zip_name)

        return data

    def _call_get_dee2_dependent(self, func_name: str, get_dee2, counts: str = None):
        """
        Makes sure the code respects the DRY principle.

        Since the getDEE2 R functions: se, Tx2Gene and srx_agg all have the same arguments
        it makes sense to use a single function to call any of the 3 above.
        """

        # TODO add instance level counts
        if not self._check_counts(counts):
            return None

        get_dee2 = self.getDEE2(legacy=TRUE) or get_dee2

        grab_func = robjects.r[func_name]

        if counts is None:
            data = grab_func(get_dee2)
        else:
            data = grab_func(get_dee2, counts)

        print(f'type of this function: {type(data)}')

        return data

    def loadFullMeta(self, zip_name: str) -> [str, DataFrame]:
        """
        Runs the loadFullMeta() R function from getDEE2 R package.

        Load Full Metadata

        This function loads the full metadata, which contains many fields.
        :param zip_name: Path to the zipfile.
        :returns: a dataframe of full metadata.
        """

        return self._call_load_functions('loadFullMeta', zip_name)

    def loadGeneCounts(self, zip_name: str) -> [str, DataFrame]:
        """
        Runs the loadGeneCounts() R function from getDEE2 R package.

        Load Gene Counts

        This function loads STAR gene level counts from a downloaded zip file.
        :param zip_name: Path to the zipfile.
        :returns: a dataframe of gene expression counts.
        """

        return self._call_load_functions('loadGeneCounts', zip_name)

    def loadGeneInfo(self, zip_name: str) -> [str, DataFrame]:
        """
        Runs the loadGeneInfo() R function from getDEE2 R package.

        Load Gene Info

        This function loads gene information. This information includes gene names
        and lengths which is useful for downstream analysis.
        :param zip_name: Path to the zipfile.
        :returns: a dataframe of gene information.
        """

        return self._call_load_functions('loadGeneInfo', zip_name)

    def loadQcMx(self, zip_name: str) -> [str, DataFrame]:
        """
        Runs the loadQcMx() R function from getDEE2 R package.

        Load Quality Control Info.

        This function loads quality control data. More information about the QC
        metrics is available from the project GitHub page:
        https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md
        :param zip_name: Path to the zipfile.
        :returns: a dataframe of quality control metrics.
        """

        return self._call_load_functions('loadQcMx', zip_name)

    def loadSummaryMeta(self, zip_name: str) -> [str, DataFrame]:
        """
        Runs the loadSummaryMeta() R function from getDEE2 R package.

        Load Summary Metadata

        This function loads the summary metadata, which are the most relevant SRA
        accession numbers.
        :param zip_name: Path to the zipfile.
        :returns: a dataframe of summary metadata.
        """

        return self._call_load_functions('loadSummaryMeta', zip_name)

    def loadTxCounts(self, zip_name: str) -> [str, DataFrame]:
        """
        Runs the loadTxCounts() R function from getDEE2 R package.

        Load Transcript Counts

        This function loads Kallisto transcript level counts from a downloaded zip file.

        :param zip_name: Path to the zipfile.
        :returns: a dataframe of transcript expression counts.
        """

        return self._call_load_functions('loadTxCounts', zip_name)

    def loadTxInfo(self, zip_name: str) -> [str, DataFrame]:
        """
        Runs the loadTxInfo() R function from getDEE2 R package.

        Load Transcript Info

        This function loads transcript information. This information includes
        transcript lengths, corresponding parent gene accession and gene symbol
        that might be useful for downstream analysis.

        :param zip_name: Path to the zipfile.
        :returns: a dataframe of transcript info.
        """

        return self._call_load_functions('loadTxInfo', zip_name)

    # TODO bundles is expected to be a table how to pass this to R from Python if previously generated?
    # TODO col needs to be selected from Galaxy and passed on. Adding it as class attribute does not make sense
    #   in this case
    # TODO the 'query' (or in our case the srr_vector) field IS NOT a vector field. It is a standard string
    #   representing a dataset name. Therefore it only accepts A SINGLE VALUE.
    #   Converting it to srr_vector is redundant.
    # TODO edit out the srr_vector, add and validate it as single string and in the front-end, when this function
    #   is chosen show a text saying that from the data_set input ONLY THE FIRST value will be kept and that it only
    #   accepts single values. At the moment tests will be done as is.
    def query_bundles(self, species: str = None, srr_vector: [str, ListVector] = None, col: str = 'SRP_accession',
                      bundles=NULL) -> [None, ListVector]:
        """
        Runs the query_bundles() R function from getDEE2 R package.

        Query whether a project bundle is available from DEE2.

        This function sends a query to check whether a dataset is available or not.

        :param species: A character string matching a species of interest.
        :param srr_vector: ('query' in getDEE2) A character string, such as  the SRA project accession
                            number or the GEO series accession number
        :param col: The column name to be queried, usually "SRP_accession" for
                    SRA project accession or "GSE_accession" for GEO series accession.
        :param bundles: optional table of previously downloaded bundles.
        :returns: a list of datasets that are present and absent.
        """
        species = self.species or species

        if not self._check_species(species):
            return None

        if not self._check_cols(col):
            return None

        query_bundles = robjects.r['query_bundles']

        srr_vector = self.convert_to_srr_vector() or srr_vector

        data = query_bundles(species, srr_vector, col, bundles)

        return data

    def queryDEE2(self, species: str = None, srr_vector: [str, ListVector] = None, metadata=NULL) -> [None, ListVector]:
        """
        Runs the queryDEE2() R function from getDEE2 R package.

        Query Whether a DEE2 Dataset is Available.

        This function sends a query to check whether a dataset is available or not.

        :param species: A character string matching a species of interest.
        :param srr_vector: A character string or vector thereof of SRA run accession numbers
        :param metadata: optional R object of DEE2 metadata to query.
        :returns: a list of datasets that are present and absent
        """
        species = self.species or species
        srr_vector = self.convert_to_srr_vector() or srr_vector
        metadata = self.metadata or metadata

        if not self._check_species(species):
            return None

        query_dee2 = robjects.r['queryDEE2']

        data = query_dee2(species, srr_vector, metadata)

        print(f'type of this function: {type(data)}')

        return data

    # Depends on getDEE2
    def se(self, get_dee2=None, counts: str = "GeneCounts"):
        """
        Runs the se() R function from getDEE2 R package.

        Create summarizedExperiment object.

        This function creates a SummarizedExperiment object from a legacy getDEE2
        dataset.

        :param get_dee2: a getDEE2 object. Not setting 'legacy' to True breaks the function.
        :param counts: select "GeneCounts" for STAR based gene counts,
            "TxCounts" for kallisto transcript level counts or
            "Tx2Gene" for transcript counts aggregated to gene level.
            Default is "GeneCounts"
        :returns: a SummarizedExperiment object
        """

        return self._call_get_dee2_dependent('se', get_dee2, counts)

    # TODO type hints for getDEE2 - what type is it?
    def srx_agg(self, get_dee2: getDEE2 = None, counts: str = 'GeneCounts'):
        """
        Runs the srx_agg() R function from getDEE2 R package.

        Summarized run data to experiments

        Sometimes, each SRA experiment data is represented in two or more runs, and
        they need to be aggregated.

        :param get_dee2: a getDEE2 object. Not setting 'legacy' to True breaks the function.
        :param counts: select "GeneCounts" for STAR based gene counts,
            "TxCounts" for kallisto transcript level counts or
            "Tx2Gene" for transcript counts aggregated to gene level.
            Default is "GeneCounts"
        :returns: a dataframe with gene expression data summarized to SRA experiment
            accession numbers rather than run accession numbers.
        """

        return self._call_get_dee2_dependent('srx_agg', get_dee2, counts)

    def Tx2Gene(self, get_dee2=None):
        """
        Runs the Tx2Gene() R function from getDEE2 R package.

        Aggregate Transcript Counts to Gene-Level Counts.

        This function converts Kallisto transcript-level expression estimates to
        gene-level estimates. Counts for each transcript are summed to get an
        aggregated gene level score.
        :param get_dee2: a getDEE2 object. Not setting 'legacy' to True breaks the function.
        :returns: a dataframe of gene expression counts.
        """
        return self._call_get_dee2_dependent('Tx2Gene', get_dee2)
