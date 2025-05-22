## [BEAST 2](http://beast2.org) package for inferring transmission trees

Bayesian Reconstruction and Evolutionary Analysis of Transmission Histories (BREATH) is a BEAST 2 package for inferring transmission trees.


## Using BREATH

See the [tutorial](https://github.com/rbouckaert/transmission/blob/main/doc/tutorial/README.md) on how to use BREATH.


## Installing by hand

Install [BEAST 2](http://beast2.org).

Download [BREATH.v0.0.1.zip](https://github.com/rbouckaert/transmission/releases/download/v0.0.1/BREATH.v0.0.1.zip). 
Then, create a BREATH subdirectory:

```
for Windows in Users\<YourName>\BEAST\2.7\BREATH
for Mac in /Users/<YourName>\/Library/Application Support/BEAST/2.7/BREATH
for Linux /home/<YourName>/.beast/2.7/BREATH
```

Here `<YourName>` is the username you use.
Unzip the file [BREATH.v0.0.1.zip](https://github.com/rbouckaert/transmission/releases/download/v0.0.1/BREATH.v0.0.1.zip) in the newly created directory.


## Build from code

An alternative way to install is to build from the source code. 
First, get code for beast2, BeastFX and this repository. Then run

```bash
ant install
```

to install the package.

## Reference

Caroline Colijn, Matthew David Hall, Remco Bouckaert.
Taking a BREATH (Bayesian Reconstruction and Evolutionary Analysis of Transmission Histories) to simultaneously infer phylogenetic and transmission trees for partially sampled outbreaks.
[BioRxiv](https://www.biorxiv.org/content/10.1101/2024.07.11.603095), 2024
[doi:10.1101/2024.07.11.603095](https://doi.org/10.1101/2024.07.11.603095)
