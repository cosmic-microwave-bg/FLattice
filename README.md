# FLattice

Updated on January 28th, 2020 ver. 0.0.2

Add the "initial_velocity" parameter to parameter.json so that we can set initial velocity of fields easily.

---

FLattice is the Fast Lattice code to caluclate the scalar feild evolution ( still in the experimental stage ).

The following items are necessary to run the code ( or fix the relevant parts ).

- cmake
- Intel compiler
- fftw

I assume that FFTW files are installed in  `/usr/local/lib` and `/usr/local/include`. If you installed FFTW in other directory, change `include_directories` and `link_directories` in the cmake file.

## How to use

Here, I introduce a simple procedure to use FLattice.

1. Create make file by `cmake`

   Move to the download directory and create `Makefile` in a `build` directory.

   ```bash
   cd /PATH/to/the/FLattice-master
   mkdir build
   cd build
   cmake ..
   ```

   `cmake` create the `Makefile` automatically.

2. Run

   Typing

   ```bash
   make
   ```

   in the `build` directory creates the executable file as `FLattice`. Then, run the code by

   ```bash
   ./FLattice
   ```

3. Check the result

   The data such as time, field averages, field variances, etc.  are written in `status.txt`. The field values in each time step are stored in `.vtk` files in `data` directory when `dim` is 2 or 3. I recommend you to use `Paraview.app` to visualize them.

## Basic Concept

 We have four directories in `src` directory.

- `src/parameter/`

  In FLattice we globally define the simulatioin parameters. Simulation parameters are set from `parameter.json`  by `setParameter()` function in `parameter.cpp` at runtime. When you want to add another paraemters, make sure that `parameter.hpp` and `paramete.cpp` files are properly modified.

- `src/evolutioin_scheme/`

  All Classes to calculate the time evolutioin are included in the `evolution` directory. I still only implement the 2nd-order or 4th-order simplectic integration scheme ( leap-frog method ), but I will add other integration scheme such as 4th-order Runge-Kutta method.
  
- `src/physical_quantity/`

  By default, the enegy density and the Q-ball charge density can be calculated in `Energy` class and `Charge` class. You can add another classes by inheriting the `PhysicalQuanitity` class if necessary.

The other directory, `src/include/` includes  two important files,  `model.hpp` and `simulator.hpp`.

- `src/include/model.hpp`

  In `model.hpp` file, we globally define the poteintial `V` and its derivative  `dV`  of fields. Harmonic oscillator potential is set initially, so set  whether you simulate under the expanding universe or not.

- `src/include/simulator.hpp`

  Class `Simulator` defined in `simlator.hpp` manages all classes necessary to simulation, such as `Field`, `evolution_shceme`, `physical_quantity` and so on. This class is pretty messy and it will be refactored in the future release.
