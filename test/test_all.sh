# Compiles and test the following settings and writes results
#
# Tobias Lindstrøm Jensen, Aalborg University, Nov. 2012
# tlj@es.aau.dk

cd ../src
make benchmark

cd ../cvxgen5010
make

cd ../data
rm sol_problem20040*
rm sol_problem260100*
rm sol_problem5010*

cd ../test
make benchmark_large P=problem20040 NF=121 M=200 N=40
make benchmark_large P=problem260100 NF=121 M=260 N=100
make benchmark_small P=problem5010 NF=487 M=50 N=10