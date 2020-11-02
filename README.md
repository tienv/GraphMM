# GraphMM
Package GraphMM performs multiple hypothesis testing using graph-based mixture model, which take advantage of clustering pattern to improve testing power. 

Read about GraphMM https://github.com/tienv/GraphMMproject/tree/master/paper_biostat/GraphMM.pdf

## Requirements
graphmm.R can be installed on any system, and requires the following key elements:

* R (≥ 3.3.0)
* Open MP support
* igraph library
* RcppArmadillo library (== 0.9.900.3.0)

## Installation Instructions

`library(devtools)` 

`install_github('igraph')`

```
url <- "https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.9.900.3.0.tar.gz"
pkgFile <- "RcppArmadillo_0.9.900.3.0.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
```

`install_github(repo = "tienv/GraphMM", build_opts = c("--no-resave-data"))`

## Load package and get started by looking at help files
`library(GraphMM)`

`help(package = GraphMM)`

## Uninstall package
`remove.packages('GraphMM')`

## Errors 
### 1: clang error (Mac OS)
If you see this error on your MacOS:

`clang: error: unsupported option ‘-fopenmp’`

You can install another clang and point R to use that clang, which supports the -fopenmp parameter.

Install llvm on your mac

`brew install llvm`

create a Makevars file

`touch ~/.R/Makevars`

Add these lines to the Makevars file

```
# comment out first line 'CC= ... if there are errors with compiling a package
CC=/usr/local/opt/llvm/bin/clang -fopenmp
CXX=/usr/local/opt/llvm/bin/clang++

# Also potentially CXX11 (for C++11 compiler)
CXX11=/usr/local/opt/llvm/bin/clang++

# -O3 should be faster than -O2 (default) level optimisation ..
CFLAGS=-g -O3 -Wall -pedantic -std=gnu99 -mtune=native -pipe
CXXFLAGS=-g -O3 -Wall -pedantic -std=c++11 -mtune=native -pipe
LDFLAGS=-L/usr/local/opt/gettext/lib -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib
CPPFLAGS=-I/usr/local/opt/gettext/include -I/usr/local/opt/llvm/include
```
This should resolve the '-fopenmp' issue on mac

## 2: igraph installation

If you are facing issues while installing devtools or igraph library, you can try the following to fix the issue. 

Following the instructions of Dean Attali (https://www.digitalocean.com/community/tutorials/how-to-set-up-r-on-ubuntu-14-04)

```
$ sudo apt-get -y install libcurl4-gnutls-dev libxml2-dev libssl-dev
$ sudo su  
$ R
```
```
install.packages('devtools', repos='http://cran.rstudio.com/')
install.packages("igraph", repos='http://cran.rstudio.com/')
```
Since the package is installed by root, it can be used by all users of the system.
