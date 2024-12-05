# HydroModels.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://chooron.github.io/HydroModels.jl/)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Overview

HydroModels.jl is a modern hydrological modeling framework that extends and enhances SUPERFLEX's design philosophy. Built on Julia language and the SciML (Scientific Machine Learning) ecosystem, it combines flexible model construction with computational efficiency, particularly supporting deep learning integration in hydrological modeling.

## Key Features

- **Flexible Model Construction**: Supports development of lumped, semi-distributed, and distributed hydrological models
- **Deep Learning Integration**: Enables neural network integration for enhanced flux calculations and dynamic parameter estimation
- **Computational Efficiency**: Leverages Julia's high-performance capabilities and the SciML ecosystem
- **Gradient-Based Optimization**: Supports advanced parameter optimization techniques
- **Comprehensive Framework**: Provides tools for both traditional hydrological modeling and modern machine learning approaches

## Framework Capabilities

HydroModels.jl offers:
- Easy implementation and customization of hydrological models
- Integration with Julia's scientific computing ecosystem
- Support for various modeling approaches:
  - Traditional conceptual models
  - Neural network enhanced models
  - Distributed hydrological systems

## Case Studies

The framework has been validated through various applications, including:
1. M50 Model Implementation
2. GR4J-based Distributed Hydrological Model

## Getting Started

Check out our tutorials to start using HydroModels.jl:
- [Run a Bucket Model](tutorials/run_a_bucket.md)
- [Run ExpHydro Model](tutorials/run_a_exphydro_model.md)

## Installation

```julia
using Pkg
Pkg.add("HydroModels")
```

## Contributing

HydroModels.jl is an open-source project available on Github. We welcome contributions from the community to help advance deep learning applications in hydrology.

## Documentation

For detailed information about using HydroModels.jl, please refer to our comprehensive documentation and tutorials in this site.
