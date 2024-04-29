#!/bin/bash
for i in {1..100}
do
  echo Verifying $i
  pace2024verifier $EXAMPLES_DIR/${i}.gr $OUTPUT_DIR/${i}.sol > $OUTPUT_DIR/${i}.ver
done
