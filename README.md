# Data Unfolding using RooUnfold Package

This code is to unfold data present in the `root` file. The name of the plots
accepted by the code is hard-coded inside. There are two cpp files:
```
fillResponse.C                 -> fill response histograms
anal_muon_unfold_basic.cc      -> to unfold full range of the plot (could be a little backdated)
anal_muon_unfold_basic_half.cc -> to unfold half of the plot (i.e. -ve and +ve seperately)
```

Install `RooUnfold` using the following commands.
```
svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk RooUnfold
cd RooUnfold
make
```
Then source it using the following by adding in `.bashrc`.
```
# RooUnfold
export ROOUNFOLD=/home/surya/products/RooUnfold/
export LD_LIBRARY_PATH=$ROOUNFOLD:$LD_LIBRARY_PATH
```
More info about installation can be found [here](https://github.com/suryamondal/various_commands/tree/main/package_installation).

**Note**: Source the same ``CERN ROOT`` used during compilation of `RooUnfold`.

### How to run:
Check the file `execute`.
