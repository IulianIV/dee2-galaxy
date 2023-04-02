import pandas.core.frame
import rpy2.robjects.vectors
from rpy2.robjects.methods import RS4
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from abc import ABC

# TODO Auto-convert all results to Pandas DF and add that as typing hints
# TODO In se objects, should we go deep enough to also convert assays.data.listData from list to
#   something more pythonic, to make other attributes, such as 'elementType' accessible?
# TODO maybe add a choice to not auto-convert to pandas df?


class ConvertedDEE2Object(ABC):
    def __getattr__(self, item):
        try:
            return rse2pyse(self.slots[item])
        except LookupError:
            raise AttributeError(f'{type(self).__name__} object has no attribute {item}')

    def list_slots(self):
        return f'{type(self).__name__} object slots: {",".join(tuple(self.slotnames()))}'

    @property
    def to_pd(self):
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

        converted = pd_from_r_df(results)

        return converted

    return wrapper


def convert_rse2pyse(func):

    def wrapper(*args, **kwargs):
        obj = func(*args, **kwargs)
        return rse2pyse(obj)

    return wrapper
