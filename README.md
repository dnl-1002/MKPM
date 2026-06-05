# Modularity community detection based on maximal $k$-plex

The parallel version is in the file "Parallel".

## Datasets
We provide some small datasets in the file "Datasets". 

The `.txt` files are the original datasets, while files with the suffix `.bin` are preprocessed files we use.

To ensure reproducibility of dataset preprocessing, espically the vertex ID using in break tie of degeneracy order, we use `toBin` to convert orginal dataset to binary format, which is the same as [ListPlex](https://github.com/joey001/ListPlex).


`toBin` is also included in “Sequential” and "Parallel" .
```
 ./toBin <input> <output>
```
**Example**
```bash
 ./toBin ./Dataset/jazz.txt ./Dataset/jazz.bin
```
All the datasets can be found on SNAP and LAW
## Experiment
### Tested Environment
* RedHatEnterpriseServer 7.9
* C++ 14
* g++ - 7.2.0
### Parallel Algorithm
#### Compile The Code
```
cd Parallel
make
```
#### Usage Instruction
**Note**: Before testing the code, we bind each thread to one CPU core.

```
export OMP_PROC_BIND=true OMP_PLACES=cores
./PlexEnum <dataset> -k <k> -q <lb> [-tau <timeout threshold>(default 0.1)] [-t <thread number>]
```
k is the maximum number of non-neighbors of a k-plex

q is the lower bound of k-plex.

tau is the user-defined task timeout threshold

t is the number of threads used

**Run The Example**
```
./PlexEnum ../Dataset/jazz.bin -k 4 -q 12 -tau 0.1 -t 16

