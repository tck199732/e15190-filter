# e15190-filter
Transport models are used to simulate the dynamics of the system in heavy ion collisions at intermediate energies. However, the results obtained assume no experimental condition and thus poses ambiguity when compared to data. This repository offers a light-weighted program to filter out the "visible" content from the raw output. These codes are usually written in fortran and thus their raw output would be large-sized text file. For experimentalist, it is advantageous to use a [ROOT](https://root.cern/) file to store millions of events. The main program will output a new ROOT file where particles not seen by detectors are discarded. This repository only offers filter for the NSCL E15190 experiment. 

### 1. Installation options
---------------------------------
Simply clone the repository
```
git clone https://github.com/tck199732/e15190-filter.git
```
The only dependency is [ROOT](https://root.cern/). To use this repository, either use conda or docker (local installation is not recommended but should work without problem):
- install [`miniconda3`](https://docs.conda.io/projects/miniconda/en/latest/) and then 
```
conda env create -f environment.yml --prefix ./env
```
- use the official [docker image](https://hub.docker.com/u/rootproject) or [singularity](https://apptainer.org/)

### 2. Get started
----------------------------------------------
Simply `. activate.sh` to activate the environment. The user should only modified [`tree.hpp`](./tree.hpp) to change the branch settings according to their ROOT file. After compiling the program,
```
./main.exe -r Ca40Ni58E56 -i ${input.root} -o ${output.root} -n {chain-name} -f cms
```
To see all options, `./main.exe -h`,
```
Optional arguments:
  -h, --help                                                     shows help message and exits 
  -v, --version                                                  prints version information and exits 
  -r, --reaction, reaction name                                  [default: "Ca40Ni58E56"]
  -i, --input-files, list of input paths, separated by comma     [nargs: 1 or more] [required]
  -o, --output-file, output file path                            [default: "output.root"]
  -n, --chain-name, name of the TChain                           [default: "AMD"]
  -f, --reference-frame, reference frame of the input particles  [default: "cms"]
  -d, --debug, debug mode                                        
  --progress-bar, show progress bar
```

### 3. Remarks
-------------------------------------------------------
- the filters applied are simple and only include
    - detector signal correction
    - $\theta_{\mathrm{lab}}$, $\phi$ cut
    - threshold energy

- the detector systems in E15190 mainly include 
    - Microball
    - HiRA10
    - veto wall
    - neutron wall A/B





