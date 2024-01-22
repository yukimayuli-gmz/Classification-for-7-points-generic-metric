# Classification-for-7-points-generic-metric

## Description

A R function to search information which be used to classify 7 points generic metric by our classification. Users can search the information of `the Base Case`, `Gromov Structure`, `Extended Graph`, `the Quartets be checked`, `the` $K_{2,3}$ `be checked`, `the number of generic metric case in this base case` and `a permutation for input generic metric to correspond to this information`.

In `rdatalist` directory, there is information for our classification which be saved as some Rdata files. And some files will be loaded when the search function executed. Our classification is executed in R, and all information we used in classification is in R memory and environment. This is the reason why we provide a R function for users to search classification information of generic metrics.

In `main.R` file, there is the main function `search_case_number_by_generic_metric` of the search function. In `functions.R` file, there are some sub-functions called by the main function. Users need to run these two scripts before using the search function.

In `test.R` file, there is a section of test code with an example of a generic metric for testing purposes. Users can use this code to check whether the search function is working properly.

This function uses sample Linear Programming, user need to install `gurobi` or `Rglpk` package in R. If you have any questions about the code, please feel free to contact us. [`guomengzhen@picb.ac.cn`](guomengzhen@picb.ac.cn)

## Usage

First, download all files, run `functions.R` and `main.R` in R.<br>
Then run the search function.<br>

    search_case_number_by_generic_metric(M_metric, data_path = {the path for storing the folder rdatalist})

## Arguments

When use the search function in R, there are 2 parameters can be set.<br>
* `M_metric`<br>
The distance matrix of a generic metric (the matrix should be a suitable generic metric and a symmetric matrix).<br>
* `data_path`<br>
The path for storing the folder `rdatalist`, default is `getwd()` which means the current path used by R.

## Value

The output of the search function is a list in R with 6 objects.<br>
* `case_number`<br>
The number of generic metric case in this base case. In each base case, all feasible generic metrics are arranged in order when classifying, and this number can be used to search the inequalities and extreme rays corresponding to this general metric (which have not been uploaded yet).<br>
* `permutation`<br>
A permutation of 7 taxa for input generic metric. The permutated new metric with new taxa ordering is corresponds to other output information. The permutation $(p_1, p_2, ..., p_7)$ is means put column $p_1$ of input matrix to column 1 of new matrix.<br>
* `Base_case`<br>
The distance matrix of a generic metric (the matrix should be a suitable generic metric and a symmetric matrix).<br>
* `Gromov_structure_and_extended_graph`<br>
The path for storing the folder `rdatalist`, default is `getwd()` which means the current path used by R.<br>
* `Quartets`<br>
The distance matrix of a generic metric (the matrix should be a suitable generic metric and a symmetric matrix).<br>
* `K23`<br>
The path for storing the folder `rdatalist`, default is `getwd()` which means the current path used by R.
