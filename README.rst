R's getDEE2 Python Wrapper
****

Installation
====

At the moment the tool is not available in the Galaxy Toolshed.
Nonetheless, it can still be used in a local Galaxy instance.

1. You need to have R installed. `How to install R <https://cran.r-project.org/doc/FAQ/R-FAQ.html#How-can-R-be-installed_003f>`_.
2. Once R is installed, you need to install getDEE2. `How to install getDEE2 <https://www.bioconductor.org/packages/release/bioc/html/getDEE2.html>`_.

The easiest way to make use of the toll is to download the files
found under :code:`tools/custom_tools` code and add them in a sub-folder from
the Galaxy Tool folder (i.e. the same way it is named in the repository)

After they are added you just need to add a reference to the tool in
an existing section or a new one. Similarly to :code:`config/tool_conf.xml.sample`

After all these steps are done, the tool should be usable in your local Galaxy instance.

Requirements
====

The tool is dependent on **rpy2** and **pandas**.

* `rpy2 <https://rpy2.github.io/>`_ is used to create a connection between python and the getDEE2 library
* `pandas <https://pandas.pydata.org/>`_ is needed because some results from the R DataFrame
    objects are converted to pandas DataFrames

Usage
====

Besides the obvious usage with Galaxy, this tool is fully usable with a command line.

Example usage

:code:`dee2.py -f getDEE2 -s ecoli -d "SRR1613487,SRR1613488"`

This produces an R SummarizedExperiment object which in turn gets converted
to Python classes and can be queried and manipulated.

The way the module is constructed makes it possible to use it as a package too. Just import :code:`dee2conn.py`
and :code:`dee2converter.py` and you can leverage the full extent of the module.

For example, in R:

.. code-block:: R

    x <- getDEE2('ecoli', c("SRR1613487","SRR1613488")
    x@colData

Would produce a table of values

The Python equivalent:

.. code-block:: python

    from dee2conn import DEE2

    dee2 = DEE2()
    dee2.data_set = 'SRR1613487,SRR1613488'
    dee2.species = 'ecoli'

    # this runs the getDEE2 function
    getdee2 = dee2.getDEE2()

    # this accesses the same values as the previous R command
    colData = getdee2.colData

    # this would convert the results to pandas DataFrame for further processing
    pd_df = colData.to_pd


All :code:`getDEE2` R results are converted to :code:`rpy2` objects and then to custom Python classes. The custom classes inherit
from base rpy2 classes, thus they ar basically the same. What they offer instead is a more user-friendly usage.

All converted objects can be transformed to Pandas Dataframes for further manipulation.

The module also offers full usability of the R module through Python. By default, R messages are turned off.
To see the R shell in progress set :code:`supress_r_warnings` to :code:`False` on your DEE2 instance.

Besides the possibility to manipulate data through pandas, the module also offers a :code:`convert_to_tsv` method which
attempts to convert the previously created pandas DataFrames to a CSV file.

Info about getDEE2 functions can be found on the packages' website.
