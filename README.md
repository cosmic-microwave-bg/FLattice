# FLattice

---

Updated on June 13th, 2019

We totally refactored the codes. Please read the **Basic Concept** below.

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

- `parameter`

  In FLattice we globally define the simulatioin parameters. Simulation parameters are set from `parameter.json`  by `setParameter()` function in `parameter.cpp` at runtime. When you want to add another paraemters, make sure that 'parameter.hpp' and 'paramete.cpp' files are properly modified.

Ther other three directories are structured . Each class is derived by its base class and its instance is created in the corresponding `create-` functions.

- `model`

  The  `model` directory includes the files of various model, i.e. poteintials. By default, `harmonic_osciilator.hpp` and 'qball.hpo' are included, so add an another model or rewrite the existing files if necessary. However, please be careful when you inhrerits the `Model` class  because the member function of `Model` class  `dV` is frequently called in solving the equaitons of motions and it slows down the calculation speed (about 1.2 - 2 times).

- `evolutioin_scheme`

  All Classes to calculate the time evolutioin are included in the `evolution` directory. I still only implement the 2nd-order or 4th-order simplectic integration scheme ( leap-frog method ), but I will add other integration scheme such as 4th-order Runge-Kutta method.
  
- 'physical_quantity'

  By default, the enegy density and the Q-ball charge density can be calculated in `Energy` class and `Charge` class. You can add an another class by inheriting the `PhysicalQuanitity` class if necessary.

