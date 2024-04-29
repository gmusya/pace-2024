#!/bin/bash
for i in {1..100}
do
  echo Solving $i
  $HEURISTIC_EXECUTABLE < $EXAMPLES_DIR/${i}.gr > $OUTPUT_DIR/${i}.sol 2>$OUTPUT_DIR/${i}.log
done
