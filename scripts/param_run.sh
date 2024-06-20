#!/bin/bash
mkdir -p $OUTPUT_DIR

for i in {97..100}
do
  echo Solving $i
  time $PARAM_EXECUTABLE < $EXAMPLES_DIR/${i}.gr > $OUTPUT_DIR/${i}.sol 2>$OUTPUT_DIR/${i}.log
done

cat $OUTPUT_DIR/{1..100}.log | grep final_score > $OUTPUT_DIR/result.log
