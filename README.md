# CAFLOOD FLUIDIT

This is an improved version of the latest published open source version of University of Exeter's CAFLOOD cellular automata
based 2D flood analysis engine. The original source is: https://git.exeter.ac.uk/caddies/caddies-caflood.

The software is licensed under GPL2 as the original.

## Compiling

You first need to checkout the `caddies-api` project (from https://git.exeter.ac.uk/caddies/caddies-api). This project
must then be checked out into that projects under `apps/caddies-caflood`.
Compilation is done by invoking `cmake` in the api project root directory.

The different solvers, e.g. OpenMP and OpenCL, can be enabled by supplying the implmentation directory for `cmake`:

```bash
cd <caddies source root directory>
cmake -DCAAPI_OCL_EVENTS=disable -DCAAPI_SPECIFIC_IMPL_DIR=$(pwd)/impls/square-cell/vn-neighbours/1-levels/opencl -B build-opencl
cmake -DCAAPI_SPECIFIC_IMPL_DIR=$(pwd)/impls/square-cell/vn-neighbours/1-levels/openmp -B build-openmp
```

The original building manual is available at https://engineering.exeter.ac.uk/media/universityofexeter/emps/research/cws/downloads/caddies/CADDIES_Build_instructions.docx

## Running the simulations

In order to run the simulation, you must have all the input files available in a directory and configured the `sim.csv`
description of the problem, as per the original manual:
https://engineering.exeter.ac.uk/media/universityofexeter/emps/research/cws/downloads/caddies/CADDIES-manual-caflood-110.zip

## New features

This repository contains several features not found in the original open source version of the CAFLOOD simulator.
 * Support for full coupling between other engines, such as 1D models like EPASWMM - done via socket based communication
 * Support for culverts connecting two points in the model (i.e. pipe like structures) via the coupling engine
 * Support for per cell Manning and infiltration coeffiecients
 * Support for per cell permeability coefficient
 * Support for per cell initial water level
 * Support for per cell, time-varying rainfall intensities
 * Various fixes for compiling the OpenCL and OpenMP versions on different platforms (Linux, Windows, MacOS X)

Upcoming features include
 * Better mass balance reporting
 * More optimizations

All the features and their usage is fully described by the source code.

## Coupled simulations

Fluidit has implemented a proprietary engine for handling the actual coupling of 1D-2D simulations between CAFLOOD and
EPASWMM engines. However, it communicates with CAFLOOD using an open, simple protocol over a TCP/IP socket. Besides
the description of the protocol below and the code present in `CouplingManager.cpp`, no further support in implementing
a coupled simulations will be provided.

One must have a coupling server process listening for a port, and then pass the port number to CAFLOOD with
`--port <portnumber>` command line argument. CAFLOOD will then connect to the port.

Every time step CAFLOOD will output a `\n` terminated string to the socket:
`HEAD,<time in seconds as double>,<depth>,<head>[...]` where head and depth (both in meters) are repeated for every coupling point
described in the `coupling.csv` (or what ever name) was given in the `sim.csv` configuration file. Likewise CAFLOOD will read
lines, that must conform to format `FLOW,<time in seconds as double>,<signed flow through the coupling point in m³/s>[...]`,
where flow must be present for every coupling point.

Either side may also send lines with commands `END` to mark end of values, `STOP` to ask halting the process, and
`WAITING,<time as decimal seconds>` to let the other side know, that results will be expected up to the given time.

The `coupling.csv` format is simple, the file lists all the coupling points, one per line: `<name>,<x coordinate>,<y coordinate>`.
The linking to the coupling file is done in the `sim.csv` by adding a line: `Coupling CSV,<name of the file>`.
