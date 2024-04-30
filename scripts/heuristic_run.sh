#!/bin/bash
mkdir -p $OUTPUT_DIR

for i in {1..100}
do
  echo Solving $i
  time $HEURISTIC_EXECUTABLE < $EXAMPLES_DIR/${i}.gr > $OUTPUT_DIR/${i}.sol 2>$OUTPUT_DIR/${i}.log
done

cat $OUTPUT_DIR/*.log | grep final_score > $OUTPUT_DIR/result.log
