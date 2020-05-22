import sys
from RDAT_utils import *
import pandas as pd

# bare bones. Usage: first argument is input json, second is output json.
# Note: write to zipped json! These are huge for chemical mapping
#
df = pd.read_json(sys.argv[1])

cdf = write_concatenated_dataframe(df)
cdf = filter_data(cdf)

cdf.to_json(sys.argv[2])