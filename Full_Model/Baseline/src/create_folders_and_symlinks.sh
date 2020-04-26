#!/bin/bash
# Create folders
mkdir -p ../input
mkdir -p ../output

# Create symlinks
ln -s ../../input_raw/formatlab_all_part1_21cty.csv ../input/part1_21cty.csv
ln -s ../../input_raw/formatlab_all_part2_21cty.csv ../input/part2_21cty.csv
