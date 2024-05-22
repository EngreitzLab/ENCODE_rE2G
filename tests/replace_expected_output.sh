#!/bin/bash
# Used when we want to change the expected output to match the latest test output run
set -e

rm -rf expected_output/generic/K562_chr22/*
cp -R test_output/generic/K562_chr22/dhs_intact_hic expected_output/generic/K562_chr22/