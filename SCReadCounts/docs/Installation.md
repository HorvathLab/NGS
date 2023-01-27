
# SCReadCounts Installation #

SCReadCounts is available as a self-contained binary package for 64-bit Intel-based Linux and MacOS (Darwin), and as Python Source. The self-contained binary packages are appropriate for most users. The pythonic version requires pysam, numpy and scipy along with other packages (See the install instructions for more details). Conda-based installation of Python source provides a simple, platform-independent installation and update procedure. 

* [Binary Package Installation](#binary-package-installation)
  * [64-bit Linux](#64-bit-linux)
  * [Intel-based MacOS](#intel-based-macos)
* [Conda-based Installation](#conda-based-installation)
* [Python Source Installation](#python-source-installation)


## Binary Package Installation ##

### 64-bit Linux ###

1. Unpack the download.
    ```
    % tar xzf SCReadCounts-*.Linux-x86_64.tgz
    % cd SCReadCounts-*.Linux-x86_64
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

### Intel-based MacOS ###
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
    % conda run -n HorvathLab --live-stream scReadCounts
    ```
4. Command-line help:
    ```
    % conda run -n HorvathLab --live-stream scReadCounts -h
    ```
5. Update scReadCounts to latest version:
    ```
    % conda update -n HorvathLab -c bioconda screadcounts
    ```

## Python Source Installation ##

1. Unpack the downloaded:
    ```
    % tar xzf SCReadCounts-*.Python-3.7.tgz
    % cd SCReadCounts-*.Python-3.7.tgz
    ```
2. Install the necessary, and optional (if desired), Python 3 packages:
    ```
    % pip3 install -r src/requirements.txt 
    % pip3 install -r src/optional_requirements.txt
    ```
3. If `python3` is not on your path or is called something else
    ```
    % export PYTHON3=<path to python3>
    ```
    or
    ```
    % setenv PYTHON3 <path to python3>
    ```
4. See the graphical user interface (if wxPython is installed).
    ```
    % bin/scReadCounts
    ```
5. Command-line help.
    ```
    % bin/scReadCounts -h
    ```
6. Run the examples.
    ```
    % cd data
    % ./example.sh
    ```

## See Also

[SCReadCounts Home](..), [Usage](Usage.md), [Input Files](InputFiles.md), [Output Files](OutputFiles.md)
