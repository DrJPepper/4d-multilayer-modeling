# 4D Multilayer Modeling

<p float="left" align="center">
  <img src="https://drjpepper.github.io/images/4d_modeler/diagram.png" width="70%" /> 
</p>

This repo contains the code from my thesis project developing a deformation
modeling framework for epithelium-inspired, responsive 4D materials.

## Dependencies
- Build tools: CMake, Make and GCC
- Libraries:
    - Qt6
    - VTK
    - OpenMP
    - OpenCL
    - FMT
    - Eigen3

## Installation

```
git clone git@github.com:DrJPepper/4d-multilayer-modeling.git
git clone git@github.com:DrJPepper/4d-modeling-example-inputs.git
cd 4d-multilayer-modeling
ln -sr ../4d-modeling-example-inputs ./input
mkdir build
cd build
cmake ..
make
```
