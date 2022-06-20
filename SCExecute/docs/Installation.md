
# SCExecute Installation #

SCExecute is available as a self-contained binary package for 64-bit Linux and MacOS systems and as Python source. The self-contained binary package is appropriate for most Linux and MacOS users. The pysam package, plus a variety of common third-party python packages including numpy and scipy must be installed to use SCExecute in Python source form. See the install instructions below for more details. 

## Binary Package Installation ##

### Linux ###
1. Unpack the download.
    ```
    tar xzf SCExecute-*.Linux-x86_64.tgz
    cd SCExecute-*.Linux-x86_64
    ```
2. See the graphical user interface.
    ```
    bin/scExecute
    ```
3. Command-line help.
    ```
    bin/scExecute -h
    ```
4. Run the examples.
    ```
    cd data
    ./example.sh
    ```
### MacOS ###
1. Unpack the download.
    ```
    % tar xzf SCExecute-*.Darwin-x86_64.tgz
    % cd SCExecute-*.Darwin-x86_64
    ```
2. See the graphical user interface.
    ```
    % bin/scExecute
    ```
3. Command-line help.
    ```
    % bin/scExecute -h
    ```
4. Run the examples.
    ```
    % cd data
    % ./example.sh
    ```

## Python Source Installation ##

1. Unpack the downloaded:
    ```
    tar xzf SCExecute-*.Python-3.7.tgz
    cd SCExecute-*.Python-3.7
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
    bin/scExecute
    ```
5. Command-line help.
    ```
    bin/scExecute -h
    ```
6. Run the examples.
    ```
    cd data
    ./example.sh
    ```

## See Also

[SCExecute Home](..), [Usage](Usage.md), [Input Files](InputFiles.md), [Cell Barcodes](Barcodes.md), [Command/Template Substitution](CommandSubst.md), [Examples](Examples.md)
