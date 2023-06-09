import pandas.core.frame
import rpy2.robjects.vectors
from rpy2.robjects.vectors import Matrix, ListVector, DataFrame
from rpy2.robjects.methods import RS4
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

# TODO see if RS4_AutoType method can be used here. More info:
#   https://rpy2.github.io/doc/v3.5.x/html/robjects_oop.html#automated-mapping-of-user-defined-classes
# TODO rpy2 DataFrame objects have a 'to_csvfile' method. Maybe this could be used to make
#   processing faster

RPY2_MATRICES = (rpy2.robjects.vectors.Matrix,
                 rpy2.robjects.vectors.FloatMatrix,
                 rpy2.robjects.vectors.IntMatrix,
                 rpy2.robjects.vectors.StrMatrix,
                 rpy2.robjects.vectors.BoolMatrix,
                 rpy2.robjects.vectors.ByteMatrix,
                 rpy2.robjects.vectors.ComplexMatrix)


class ConvertedDEE2Object:
    """
    Class that handles conversion from R's S4 OOP Objects to Python customized Classes.
    At the moment it does not offer much functionality besides:
    - accessing a S4 objects properties, like in R;
    - conversion to pandas DataFrame
    - list of S4 Object Slotnames
    """
    def __getattr__(self, item):
        try:
            return ro2pyo(self.slots[item])
        except LookupError:
            raise AttributeError(f'{type(self).__name__} object has no attribute {item}')

    def list_slots(self):
        return f'{type(self).__name__} object slots: {",".join(tuple(self.slotnames()))}'

    def to_pd(self) -> pandas.core.frame.DataFrame:
        return pd_from_r_df(self)


class SummarizedExperiment(RS4, ConvertedDEE2Object):
    pass


class DFrame(RS4, ConvertedDEE2Object):
    pass


class SimpleAssays(RS4, ConvertedDEE2Object):
    pass


class SimpleList(RS4, ConvertedDEE2Object):
    pass


class ConvertedListVector(ListVector, ConvertedDEE2Object):
    def __getattr__(self, item):
        try:
            return ro2pyo(self.rx2(item))
        except LookupError:
            raise AttributeError(f'{type(self).__name__} object has no attribute {item}')


class ConvertedMatrix(Matrix, pandas.DataFrame):
    """
    Class that handles conversion from R's Matrix Objects to Python customized Classes.
    It handles automatic conversion from any type of R Matrix (FloatMatrix, ByteMatrix etc.)
    All regular properties & attributes of rpy2 Matrix class are inherited and can be accessed
        through these classes 'r_matrix' attribute.

    Additionally, the class covers conversion from R Matrix to Pandas DataFrame
    """
    def __init__(self, matrix: RPY2_MATRICES):
        super().__init__()
        self.r_matrix = matrix
        self.indices = self._indices
        self.df_matrix = self._df_matrix

    @property
    def _df_matrix(self):
        n_cols = self.r_matrix.ncol
        col_names = self.r_matrix.colnames
        matrix_dict = dict()
        matrix_df = pandas.DataFrame

        for col in range(1, n_cols + 1):
            matrix_dict[col_names[col - 1]] = [row for row in self.r_matrix.rx(True, col)]

        matrix_dict['indices'] = self.indices

        matrix_df = matrix_df.from_dict(matrix_dict)
        matrix_df.set_index('indices', drop=True, inplace=True)

        return matrix_df

    @property
    def _indices(self):
        indices = [row_name for row_name in self.r_matrix.rownames]
        return indices


def ro2pyo(obj: [RS4, RPY2_MATRICES]):
    """
    Function that determines which subclass of ConvertedDEE2Object or other ConvertedObject to call
    """
    if 'SummarizedExperiment' in obj.rclass:
        res = SummarizedExperiment(obj)
    elif 'DFrame' in obj.rclass:
        res = DFrame(obj)
    elif 'SimpleAssays' in obj.rclass:
        res = SimpleAssays(obj)
    elif 'SimpleList' in obj.rclass:
        res = SimpleList(obj)
    elif isinstance(obj, RPY2_MATRICES):
        res = ConvertedMatrix(obj)
    elif isinstance(obj, ListVector):
        res = ConvertedListVector(obj)
    else:
        res = obj
    return res


def pd_from_r_df(r_df: [DataFrame, ConvertedDEE2Object]) -> pandas.core.frame.DataFrame:
    """
    Function that handles conversion from R DataFrame to Pandas DataFrame
    """
    try:
        r_df = ro.DataFrame(r_df)
    except ValueError:
        r_df = ro.DataFrame({'results': r_df})

    with (ro.default_converter + pandas2ri.converter).context():
        pd = ro.conversion.get_conversion().rpy2py(r_df)

    return pd


def convert_rdf_to_pd(func):
    """
    Handy wrapper for 'pd_from_r_df' that helps to automatically convert R DataFrame results from getDEE2 functions
    to Python Pandas DataFrame
    """

    def wrapper(*args, **kwargs):
        results = func(*args, **kwargs)

        converted = pd_from_r_df(results)

        return converted

    return wrapper


def convert_ro2pyo(func):
    """
    Wrapper that helps with automatic conversion of R S4 object to Python Objects
    """
    def wrapper(*args, **kwargs):
        obj = func(*args, **kwargs)
        return ro2pyo(obj)

    return wrapper


# TODO If I am not mistaken this is a ListVector too. Therefore, this should be included, if possible, in the
#   Matrix & List Class above.
#   Although, the ConvertedListVector class specifically calls for a Lists' member. This converts it as a whole.
def convert_query(func):
    """
    Wrapper that helps with the conversion of R ListVector to Pandas DataFrame
    """

    def wrapper(*args, **kwargs):
        obj = func(*args, **kwargs)

        bundle = dict()
        for k, v in obj.items():
            bundle[k] = [x for x in v]

        bundle_df = pandas.DataFrame.from_dict(bundle, orient='index')
        bundle_df = bundle_df.transpose()

        return bundle_df

    return wrapper
