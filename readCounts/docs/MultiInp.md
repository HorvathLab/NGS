# readCounts Multi Inputs

## Description

readCounts enable the users to run multiple inputs which can be done through
writing a shell script as illustrated below:

### Running Multiple Inputs:

    #!/bin/sh
    while read line
    do
    {
    python src/readCounts.py \
     -s "variants.txt" \
     -r $line".sorted.bam" \
     -o $line"_10rnaed.csv" \
     -m 10 \
     -F False \
     -f False \
     -U False \
     -t 20
    }
    done < list

## Shell Script Description:

The shell script is an example of how one can run multiple inputs (samples) to
do so, firstly, you need to have a list which contains all of the inputs file
names without the file extensions. Then, the list will go through a loop that
iterates one
sample(input) at a time. The readCounts will run the iterated input and once
done it will go to the next input of list until it reaches the end of the list.

