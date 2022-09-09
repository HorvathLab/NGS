
# SCReadCounts Installation #

SCReadCounts is available as a self-contained binary package for 64-bit Linux systems, as Python source, and MacOS (Darwin). The self-contained binary package is appropriate for most Linux and MacOS users. The pythonic version requires pysam, numpy and scipy along with other packages (See the install instructions for more details). 


## Binary Package Installation ##

### Linux ###

1. Unpack the download.
    ```
    tar xzf SCReadCounts-*.Linux-x86_64.tgz
    cd SCReadCounts-*.Linux-x86_64
    ```
2. See the graphical user interface.
    ```
    bin/scReadCounts
    ```
3. Command-line help.
    ```
    bin/scReadCounts -h
    ```
4. Run the examples.
    ```
    cd data
    ./example.sh
    ```

### MacOS ###
1. Unpack the download.
    ```
    % tar xzf SCReadCounts-*.Darwin-x86_64.tgz
    % cd SCReadCounts-*.Darwin-x86_64.tgz
    ```
2. See the graphical user interface.
    ```
    % bin/scReadCounts
    ```
3. Command-line help.
    ```
    % bin/scReadCounts -h
    ```
4. Run the examples.
    ```
    % cd data
    % ./example.sh
    ```

## Conda-based Installation ##

1. Create a conda environment for **HorvathLab** tools (if not done previously)
    ```
    % conda create -n HorvathLab
    ```
2. Install in the conda environment
    ```
    % conda install -n HorvathLab -c bioconda screadcounts
    ```
3. Run scReadCounts:
    ```
    % conda run -n HorvathLab scReadCounts
    ```
4. Update scReadCounts to latest version:
    ```
    % conda update -n HorvathLab scReadCounts
    ```

## Python Source Installation ##

1. Unpack the downloaded:
    ```
    tar xzf SCReadCounts-*.Python-3.7.tgz
    cd SCReadCounts-*.Python-3.7.tgz
    ```
2. Install the necessary, and optional (if desired), Python 3 packages:
    ```
    pip3 install -r src/requirements.txt 
    pip3 install -r src/optional_requirements.txt
    ```
3. If `python3` is not on your path or is called something else
    ```
    export PYTHON3=<path to python3>
    ```
    or
    ```
    setenv PYTHON3 <path to python3>
    ```
4. See the graphical user interface (if wxPython is installed).
    ```
    bin/scReadCounts
    ```
5. Command-line help.
    ```
    bin/scReadCounts -h
    ```
6. Run the examples.
    ```
    cd data
    ./example.sh
    ```

## See Also

[SCReadCounts Home](..), [Usage](Usage.md), [Input Files](InputFiles.md), [Output Files](OutputFiles.md)
