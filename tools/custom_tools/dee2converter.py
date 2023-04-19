import pandas.core.frame
import rpy2.robjects.vectors
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.methods import RS4
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from abc import ABC

# TODO see if RS4_AutoType method can be used here. More info:
#   https://rpy2.github.io/doc/v3.5.x/html/robjects_oop.html#automated-mapping-of-user-defined-classes


class ConvertedDEE2Object(ABC):
    def __getattr__(self, item):
        try:
            return rse2pyse(self.slots[item])
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


def rse2pyse(obj: RS4):
    if 'SummarizedExperiment' in obj.rclass:
        res = SummarizedExperiment(obj)
    elif 'DFrame' in obj.rclass:
        res = DFrame(obj)
    elif 'SimpleAssays' in obj.rclass:
        res = SimpleAssays(obj)
    elif 'SimpleList' in obj.rclass:
        res = SimpleList(obj)
    else:
        res = obj
    return res


def pd_from_r_df(r_df: [rpy2.robjects.vectors.DataFrame, ConvertedDEE2Object]) -> pandas.core.frame.DataFrame:
    try:
        r_df = ro.DataFrame(r_df)
    except ValueError:
        r_df = ro.DataFrame({'results': r_df})

    with (ro.default_converter + pandas2ri.converter).context():
        pd = ro.conversion.get_conversion().rpy2py(r_df)

    return pd


def convert_rdf_to_pd(func):

    def wrapper(*args, **kwargs):
        results = func(*args, **kwargs)

        # TODO Use data from here: https://rpy2.github.io/doc/v3.5.x/html/vector.html#rpy2.robjects.vectors.Array
        #   to process and parse Matrix and FloatMatrix types.
        #   They are necessary for srx_agg() to work.
        # for row in results.rx(True, 1):
        #     print(row)

        if 'outfile' in kwargs and kwargs['outfile'] is not None:
            return results

        converted = pd_from_r_df(results)

        return converted

    return wrapper


def convert_rse2pyse(func):

    def wrapper(*args, **kwargs):
        obj = func(*args, **kwargs)
        return rse2pyse(obj)

    return wrapper


def convert_query(func):

    def wrapper(*args, **kwargs):
        obj = func(*args, **kwargs)

        bundle = dict()

        for k, v in obj.items():
            bundle[k] = [x for x in v]

        bundle_df = pandas.DataFrame.from_dict(bundle)

        return bundle_df

    return wrapper

