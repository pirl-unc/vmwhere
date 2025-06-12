#!/bin/bash

# Example usage of vmwhere with test data
echo "Running vmwhere visualize with example data..."


# Check if vmwhere is installed
if ! command -v vmwhere &> /dev/null; then
    echo "Error: vmwhere is not installed. Please run 'pip install -e .' from the project root."
    exit 1
fi

# Check if input CSV file exists
if [[ ! -f "output/example_sample_allele_calls_results.csv" ]]; then
    echo "Error: Input CSV file 'example_sample_allele_calls_results.csv' not found!"
    echo "Make sure you've run the profile command first to generate the CSV file for the example sample."
    exit 1
fi


# kick off profiling function
vmwhere visualize \
    --sample_csv output/example_sample_allele_calls_results.csv \
    --chr chr6 \
    --start 6706603 \
    --output_pdf output/chr6_example_region_visualization.pdf 

# Check if the command succeeded
if [ $? -eq 0 ]; then
    echo "Example visualization complete! Check output/ directory for example results."
    echo "Generated file: output/chr6_example_region_visualization.pdf"
else
    echo "Error: vmwhere visualize command failed."
    exit 1
fi

