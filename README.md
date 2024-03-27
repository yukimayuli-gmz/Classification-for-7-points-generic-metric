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
The number of generic metric case in this base case. In each base case, all feasible generic metrics are arranged in order when classifying, and this number can be used to search the inequalities and extreme rays corresponding to this general metric (which have been uploaded in [Zenodo](https://doi.org/10.5281/zenodo.10551370)).<br>
* `permutation`<br>
A permutation of 7 taxa for input generic metric. The permutated new metric with new taxa ordering is corresponds to other output information. The permutation $(p_1, p_2, ..., p_7)$ is means put column $p_1$ of input matrix to column 1 of new matrix.<br>
* `Base_case`<br>
The number of the `Base Case` to which input generic metric belongs.<br>
* `Gromov_structure_and_extended_graph`<br>
The `Gromov Structure` and `Extended Graph` to which input generic metric belongs. This is an integer matrix with 6 columns. Column 1 represents the taxa 1 to 7. Column 1, 2 and 3 represent the three taxa which form the minimal Gromov product at the taxa in colomn 1 in Gromov Structure. The number in column 4 is 8 or 9, the number is 8 means using `Choice 1` to get a extended $(3,x)$-soluiton graph, the number is 9 means using `Choice 2` or `Choice 3` to get a extended $(3,x)$-soluiton graph. Column 5 and 6 represent the two taxa be selected when obtaining the extended graph.<br>
* `Quartets`<br>
Quartets that are consistent with the general metric. We only record those quartets that remain undetermined in the early steps of the classification process and need to be checked. The `Quartets` object is 0 means all 35 quartets are determined and there are no quartets need to be checked. Otherwise, this is a matrix or vector (if only one row) with 5 columns. A row is $(x, y_1, y_2, y_3, z)$ represents a quartets with $d(x,y_{z}) > d(y_{z'}, y_{z''})$, $z \in$ {1,2,3} and $z \neq z', z''$.<br>
* `K23`<br>
$K_{2,3}$ that are consistent with the general metric. We only record those $K_{2,3}$ that remain undetermined in the early steps of the classification process and need to be checked. The $K_{2,3}$ object is 0 means all $K_{2,3}$ are determined and there are no $K_{2,3}$ need to be checked. Otherwise, this is a matrix or vector (if only one row) with 6 or 7 columns. In a row, the first 6 column is $(x_1, x_2, y_1, y_2, z, w)$ represents a $K_{2,3}$ with $x_{1}x_{2}, y_{1}y_{2}$ are the two edges and $z$ is the `central point`. $w = 1$ means that only one edge of $x_{1}x_{2}$ and $y_{1}y_{2}$ can form a triangle with $z$ and get a feasible $K_{2,3}$ graph, so, during the classification process, we do not need to check them. If there are only one row with $w = 2$, we also did not check them, because there are only two feasible generic metrics, each satisfying the two distinct $K_{2,3}$ graphs. If there are more than two rows with $w=2$, we will add column 7 to the matrix to indicate which $K_{2,3}$ graph the generic metric be consistent with. If the column 7 is 1, it indicates that the edge $x_{1}x_{2}$ forms a triangle with $z$. Conversely, if the column 7 is 2, it indicates that the edge $y_{1}y_{2}$ forms a triangle with $z$.

## Inequalities and Extreme Rays

We have uploaded all the `inequalities` and `extreme rays` data corresponding to every general metrics to [Zenodo](https://doi.org/10.5281/zenodo.10551370). They are quite extensive, so we segmented them based on 37 base cases. User can use `case_number` to find the inequalities and extreme rays corresponding to a general metric.

When saving the extreme ray data, we saved the data for every 10,000 cases in a text file `rays_case_{i}_{j}.txt`, so the $i$ is the number of base case and $j$ is the `case_number` divided by 10000 and rounded up. And used

        ###
        case_number
        ###
        
to indicate the separation at the end of extreme rays for each case.

We have classified all prime metrics of 7 points into 13,182 equivalence classes (the prime metrics of the `Metric Cone` are the extreme rays of every generic metrics). And we saved them in the file `isomorphic_prime_rays_list_with_7_points.txt`. Within the file `isomorphic_prime_rays_list_with_7_points.txt`, we further divided these equivalence classes into 7,783 different categories based on the frequency of each number in the prime metrics. In the `isomorphic_prime_rays_list_with_7_points.txt` file, the first line of each category represents a vector indicating the frequency of appearance of different numbers starting from 0 in that category of prime metrics. Each subsequent line represents an equivalence class of prime metrics. Finally, the 7,783 different categories are separated by

        ###
        category number
        ###.

We also saved the label of every extreme rays in the 13,182 equivalence classes, for every 10,000 cases in a text file `rays_label_case_{i}_{j}.txt`, so the $i$ is the number of base case and $j$ is the `case_number` divided by 10000 and rounded up. In the file storing labels, each line represents a vector with 4 numbers (a, b, c, d), where $a$ represents the `case_number`, $b$ represents the ray number in this case, $c$ indicates which of the 7,783 categories this ray belongs to, and $d$ indicates which equivalence class within this category the ray belongs to.

## Gromov Structure

We also upload the `Gromov Product Structure` data with 8 and 9 points to [Zenodo](https://doi.org/10.5281/zenodo.10551370). We devide the classes of 8 and 9 points into three parts: <br>
* `8_points_gromov_structure_with_8_degree_0.txt` and `9_points_gromov_structure_with_9_degree_0.txt`<br>
We suppose 8 (or 9) appears only once in 8 (or 9) triples of a *Gromov product structure*.<br>
* `8_points_gromov_structure_with_8_degree_1.txt` and `9_points_gromov_structure_with_9_degree_1.txt`<br>
We suppose 8 (or 9) appears twice in 8 (or 9) triples of a *Gromov product structure*.<br>
* `8_points_gromov_structure_2_regular_multigraph_*.txt` and `9_points_gromov_structure_2_regular_multigraph_*.txt`<br>
We suppose all points appear three times in 8 (or 9) triples of a *Gromov product structure*. Then we consider 8 (or 9) edges by deleting the vertexes from all triples, the the edges form a 2-regular multigraph. The asterisk (*) in the file name represents the shape of the 2-regular multigraph. For example, in "4_2_2", it indicates that the 2-regular multigraph includes one square and two edges.<br>

We use $(x;y,z)$ to represent a triple, which is consistent with the symbols used in our article. Another commonly used representation is $\Delta_{xyz}$.

## Extreme metric

We also upload the list of 3,918 equivalence classes of `extreme metrics` with 8 points to [Zenodo](https://doi.org/10.5281/zenodo.10551370). The file name is `isomorphic_extreme_metrics_list_with_8_points.txt`. The format of this file is the same as the file `isomorphic_prime_rays_list_with_7_points.txt` we described ealier. Note that we manually added the trivial split at the end of the list, which is the item 1,213.

## Reference
