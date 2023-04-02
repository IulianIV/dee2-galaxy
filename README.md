# R's getDEE2 Python & Galaxy Implementation

## Installation

At the moment the tool is not available in the Galaxy Toolshed.
Nonetheless, it can still be used in a local Galaxy instance.

1. You need to have R installed. How to install R.
2. Once R is installed, you need to install getDEE2. How to install getDEE2.

The easiest way to make use of the toll is to download the files
found under `tools/custom_tools` and add them in a sub-folder from
the Galaxy Tool folder (i.e. the same way it is named in the repository)

After they are added you just need to add a reference to the tool in
an existing section or a new one. Similarly to `config/tool_conf.xml.sample`

After all these steps are done, the tool should be usable in your local Galaxy instance.

## Requirements

The tool is dependent on **rpy2** and **pandas**.

* rpy2 is used to create a connection between python and the getDEE2 library
* pandas is needed because some results from the R DataFrame objects are
    converted to pandas DataFrames

## Usage

Besides the obvious usage with Galaxy, this tool is fully usable with a command line.

Example usage

`dee2.py -f getDEE2 -s ecoli -d "SRR1613487,SRR1613488"`

This produces an R SummarizedExperiment object which in turn gets converted
to Python classes and can be queried and manipulated.

For example, in R:

```
x <- getDEE2('ecoli', c("SRR1613487","SRR1613488")
x@colData
```
Would produce a table of values

The Python equivalent:

```python
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
```

Info about getDEE2 functions can be found on the packages' website.

R functions which return SummarizedExperiments object need to explicitly
be converted to pandas DataFrames since some attributes of those classes
are not convertible.

Functions which return R DataFrame objects get automatically converted.