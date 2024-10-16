# Description
The tool helps to count transition frequency matrix from a sequence of events. 
This set of functions is written for a project, dedicated to a study of spider behavior sequences. An approach was inspired by PAST software (Hammer et al., 2001): the tool, implemented in PAST was not enough versatile and supported encoding of events by digits only, that seemed uncomfortable for a quick work.
# Disclaimer
Testing of an algorythm was conducted via comparison with results of an analysis of the same data in PAST. 
# Functions
## Main functions
```txt_seq_read(file, start='1', stop='2', include=True, printres=True):```

Reads text file with events in a column. Different observations are separated by an empty line. Start and end of an observation can be specified additionally (in our data it was '1' and '2' respectively, thus we put these values as default values).
Returns 2D array, each member contains a sequence of string values.
<strong>Arguments:</strong>

<i>file (string) - </i> positional - the name of a txt file to analyze in a working folder.\
<i> start (string) </i> keyword, default "1". Means the beginning of an observation.\
<i> stop (string) </i> keyword, default "2". Means the finishing of an observation. \
<i> include (boolean) </i> keyword, default True. Setting of a false value allow to drop start and stop events from analysis. \
<i> printres (boolean) </i> keyword, default True. Can be switched of to avoid from printing of a result to a console. 
\
\
```blanc_table(data, start='1', stop='2', order=True, include=False):```\
Gets 2D array from ```txt_seq_read``` and convert it to a blanc pandas DataFrame. Can Try to sort events by their usual position. Made as a separate function in order to make possible to add events, which absent in original data to compare different behavior patterns. \
<strong>Arguments:</strong>\

<i>data (list) - </i> positional - 2D array of events.\
<i> start (string) </i> keyword, default "1". Means the beginning of an observation.\
<i> stop (string) </i> keyword, default "2". Means the finishing of an observation. \
<i> order (boolean) </i> keyword, default True. Finds a positions of events in each observation and tries to sort table attributes according to position. Switching off leads to albhabetical sorting of events in a table.\
<i> include (boolean) </i> keyword, default True. Setting of a false value allow to drop start and stop events from analysis. 

```freq_table(blanc, data, flatten=False, saveres=True)```\
Gets empty table and 2D array of data. In experimental needs array can be flattened. Returns pandas DataFrame with a chart of transitions from row value to column value. \
<strong>Arguments:</strong>\
<i> saveres (boolean) </i> keyword, default True. Records a result as xlsx file.  

```prob_table(blanc, trans, rd=2)```\
Convert frequencies into probabilities. Gets blanc table and freq_table

```permutate(data, iters=9999)```\
Gets 2D array of prepared data and resamples it. Resampling method is to copy each observation for pre-defined number of times, and each time randomly spare 2 events. Keyword 'iters' means a number of copies\

```get_seq_from_lens(file, sep = ' - ', flatten=False)```
Actually may work not very clear still because it is written for only one instance.
Gets values from txt, where events and their start times are listed in a column. separator can b 

# References
An algorithm of an analysis is described in Davis, J.C. 1986. Statistics and Data Analysis in Geology. John Wiley & Sons.
# ToDo list:
We are planning to add graph building, simple filters and chi square module. In order to make faster, we will not use pandas engine to count transitions in next versions. We also plan to build a lib, now you can just fetch "markovish.py" to your working directory.\
Typical pipeline is in "example.ipynb", sample data in 'sample data.txt'