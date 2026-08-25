# 4D Multilayer Modeling

<p float="left" align="center">
  <img src="https://drjpepper.github.io/images/4d_modeler/diagram.png" width="70%" /> 
</p>

This repo contains the code from my
[thesis project](https://researchdiscovery.drexel.edu/esploro/outputs/doctoral/An-epithelium-inspired-deformation-modeling-framework-for/991022193194504721)
developing a deformation
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

Installation instructions assume you are using Linux, but may also work on Mac.

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

## Running

Grids are generated using TOML configuration files, examples of which can be
found in the
[4d-modeling-example-inputs repo](https://github.com/DrJPepper/4d-modeling-example-inputs).
All config files also have an associated model file consisting of either a
triangle mesh or the control points for a bicubic Bezier patch. From the
`build` directory, a config file can be passed by using the `-f` option:

    ./4d_multilayer_modeler -f ../input/configs/4dp/egg_chair.toml

This will initiate the surface analysis process, after which the GUI will open
and show the results of the analysis, followed by the generated grid, then the
physics simulation interface.

Details on the customizable parameters in the config files can be found in the
[example input repo](https://github.com/DrJPepper/4d-modeling-example-inputs).

## Attributions

I made modifications to [tqdm.h](https://github.com/tqdm/tqdm.cpp) and
[ShapeOp](https://www.shapeop.org/index.php), and thus have included them in
this repository. `tqdm.h` is licensed under the MIT licence, and ShapeOp is
licensed under the Mozilla Public License Version 2.0.
