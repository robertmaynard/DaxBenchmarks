#bin/bash

./BenchmarkCuda --file=/temp/test.nhdr --pipeline=1 > cuda_threshold.csv
./BenchmarkCuda --file=/temp/test.nhdr --pipeline=2 > cuda_mc.csv

./BenchmarkTBB --file=/temp/test.nhdr --pipeline=1 > tbb_threshold.csv
./BenchmarkTBB --file=/temp/test.nhdr --pipeline=2 > tbb_mc.csv

./BenchmarkSerial --file=/temp/test.nhdr --pipeline=1 > serial_threshold.csv
./BenchmarkSerial --file=/temp/test.nhdr --pipeline=2 > serial_mc.csv
