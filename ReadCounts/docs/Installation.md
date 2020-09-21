
# readCounts Installation #

readCounts is currently available only as a Python source.


## Python 2.7 readCounts ##

1. Unpack the downloaded release:

    ```
    tar xzf readCounts-*.tgz
    ln -s readCounts-* readCounts
    ```

2. Locate your Python binary and ensure it is version 2.7:

    ```
    python --version

    /path/to/python2.7 --version 
    ``` 

   We refer to the Python binary as `python` below, please substitute
   whatever path and version numbers are required to run Python 2.7 on
   your system. We recommend the Enthought Python Distribution (EPD) which
   pre-installs all but the pysam third-party dependencies needed by readCounts.

3. Ensure the necessary third-party Python modules are installed. pysam version >= 0.8.1 is required. 

    ```
    pysam, numpy, scipy
    ```

   For the configuration and execution GUI (optional):
    
    ```
    wx
    ```

   For Excel format SNV input files (optional):

    ```
    xlrd, openpyxl
    ```

   The existence of required modules can be tested as follows (demonstrated here for `scipy`):

    ```
    python2.7 -c "from scipy import __version__; print __version__"
    ```

4. The readCounts program is located in the src subdirectory:

    ```
    python ./readCounts/src/readCounts.py -h
    
    python ./readCounts/src/readCounts.py
    ```

5. Test the installation using the provided example data:

    ```
    cd readCounts
    python ./src/readCounts.py -r "data/example-*.bam" -s "data/example-SNV.tsv"
    -o output.txt
    ```
