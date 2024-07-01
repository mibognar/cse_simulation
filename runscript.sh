#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <filename_to_run> <parameter1> <parameter2>"
    exit 1
fi

# Assign the parameters to variables
FILENAME_TO_RUN=$1
PARAM1=$2
PARAM2=$3

# Run the R script in the background with the provided parameters
nohup Rscript "$FILENAME_TO_RUN" 100 "$PARAM1" "$PARAM2" > output.log 2>&1 &

# Notify the user
echo "R script $FILENAME_TO_RUN is running in the background with parameters $PARAM1 and $PARAM2."
echo "Output is being logged to output.log"
