# CFDeez
CFDeez is a 2D finite-volume computational fluid dynamics (CFD) application and library,
developed entirely from scratch. While primarily created as a learning tool, it can also serve as a helpful resource
for those interested in CFD.

> [!WARNING]
> Under active development, bugs to be expected.

<p align="center">
  <img src="./demo.png" alt="Demo Image" width="400"/>
  <p align="center"><i>This shows a passive scalar, like dye or smoke, carried by fluid around an obstacle</i></p>
</p>

## Motivations
CFD has a reputation of being difficult. So the goal was to see if 2 univeristy students with big egos but no real understanding
of anything could pull it off. Development is still underway so our egos are in a superposition of being massive and also
heavily bruised.

## Quickstart
To get started with CFDeez, you’ll need the [Odin](https://odin-lang.org/) compiler installed, as no pre-built binaries will be provided for early versions.

```bash
# Clone the repository
git clone https://github.com/Rwn-A/cfdeez
cd cfdeez

# Build the CLI tool
odin build src/cmd/cli -out:cfdeez -o:aggressive

# Run a basic example
./cfdeez ./examples/00_basic.fml
```

After this results will be saved to a `.out` directory located beside your executable.

## Project Structure
CFDeez is comprised of 3 components.
1. CFD Library `/src/cfd`
    - The core cfd code is its own library, it includes code for the fields, discritizations, linear systems and more
    - Includes optional sub-packages for outputting fields, reading `.msh` files, and pre-built solver algorithims.
2. Configuration Language `/src/fml`
    - CFDeez uses a custom JSON esque configuration language called FML, FML can be used standalone
    although I can't think of many use cases outside of this application.
3. Application Interface `/src/cmd`
    - The configuration language and CFD library are wrapped together into a single
      executable with a predefined configuration schema.

> [!TIP]
>This documentation covers using cfdeez as a complete application.
>If you're interested in using the CFD library or the FML configuration system on their own,
>see the individual `README.md` files in the `src/cfd/` and `src/fml/` directories.

## Features and Limitations
CFDeez is not designed as production CFD software. It is CPU only and single threaded.
While we might add GPU processing and/or threading, it is not a guarantee and outside the scope of this project for now.
The primary goal of this project is to understand how CFD works, and write readable code that represents the physics. Specific
discritization schemes are currently not configurable.

Features:
- 🔃 2D incompressible flow (mostly works but is unstable so not part of the repo yet.)
- ✅ 2D transport
- ✅ Unstructured meshes
- ✅ Transient and steady-state
- ✅ Output to VTU and CSV
- ✅ Configurable

Currently, we have no 1st party post-processing tools, best bet is to output to vtu and use something like
[ParaView]((https://www.paraview.org/))

### Potential Features Under Consideration
In no particular order, the below features are things we are interested in but have no current plans to attempt to implement.
- Turbelence Modelling
- 3D
- 1st Party Mesh Generation
- Compressible Flow
- Configurable schemes
- Native GUI application

## Accuracy
This is a engineering simulation in spirit so we have tried to use accurate methods. That being said, most discritizations are only first-order accurate, but implicit for stability.
Higher order methods may be implemented in the future.

## Setup & Configuration
The best way to learn the configuration options is to review the examples directory.
Each example is numbered in order of complexity. There are also some prebuilt `.msh` files available.

### Writing your own config
A configuration is 2 parts.
1. `.fml` file - Most of the configuration options are here, this is the file you pass to the executable.
2. `.msh` file - this defines the mesh for your simulation, the path to this file is referenced in the main config file.

To build a `.msh` file [gmsh](https://gmsh.info/) is required. Or you can convert from another format using something like
[meshio](https://github.com/nschloe/meshio).

**Note:** only version 2.2 of the `.msh` format is supported for now.

### Config Options
Below is a reference of all options available in the configuration file.

**Note:** Defining some of the below options will require that other options also be defined. The application will tell you
if your missing an option and you can refer here for what it means and how to define it.
> [!NOTE]
> Section is blank for now until project matures.


## License
The CFD library and application code are licensed under the **GPL** license.
The FML configuration language is licensed separately under the **MIT** license.

For specific licensing information, please refer to the LICENSE file in the root directory
and the LICENSE file located in the `src/fml` directory.

## Contributing
Any and all contributions welcome, that does **not** mean any and all contributions will be accepted.