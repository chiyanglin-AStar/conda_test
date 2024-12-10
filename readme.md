## this is conda environment in gitpod 

ref from [URL](https://anaconda.org/conda-forge/fenics-dolfin)

         [url](https://anaconda.org/conda-forge/fenics-dolfinxy)
### run fenics environment 

    conda create -n fenics-env -c conda-forge fenics matplotlib

    conda activate fenics-env

### fenics environment 2

    conda install conda-forge::fenics-dolfinx

### fenics install 

    !sudo apt-get install software-properties-common 
    !sudo apt-get install -y -qq software-properties-common python-software-properties module-init-tools
    !sudo add-apt-repository -y ppa:fenics-packages/fenics
    !sudo apt-get update -qq
    !sudo apt install -y --no-install-recommends fenics

*or* 

    sudo apt install python3-dolfinx-complex


*or*

    sudo apt install -y --no-install-recommends fenics

