# Data Unfolding using RooUnfold Package

This code is to unfold data present in the `root` file. The name of the plots
accepted by the code is hard-coded inside.

Install `RooUnfold` using the following commands.
```
svn co https://svnsrv.desy.de/public/unfolding/RooUnfold/trunk RooUnfold
cd RooUnfold
make
```
The source it using the following.
```
# RooUnfold
export ROOUNFOLD=/home/surya/products/RooUnfold/
export LD_LIBRARY_PATH=$ROOUNFOLD:$LD_LIBRARY_PATH
```
More info about installation can be found [here](https://github.com/suryamondal/various_commands/tree/main/package_installation).

**Note**: Source the same ``ROOT`` used during compilation of `RooUnfold`.

### How to run:
Check the file `execute`.