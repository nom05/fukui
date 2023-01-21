# fukui v220331
A Program to Compute Fukui Indices and Atomic Overlap Matrices

## Authors
  - Nicolás Otero Martínez     -     nom05 (at) uvigo.es - University of Vigo
  - Marcos Mandado Alonso      -   mandado (at) uvigo.es - University of Vigo
  - Ricardo A. Mosquera Castro -  mosquera (at) uvigo.es - University of Vigo

## How to use it
Type: ./fukui.x wfnfile stk\_file/cube\_file [# atom]
 - Extension for stk and wfn files are optional, extension for cube file is mandatory.
 - stk\_file can be either a FORTRAN binary file or a ordinary text file.
   They have to include, following this order: x, y, z, Gauss quadrature weight and electron density values.
   If the ordinary text format is considered, then the format to read must be specified. the label free can be used to use free format specification.
 - [# atom] is optional and enables gauss. func. skipping.                                                                                                                          
   If wfnfile.skd exists, it reads cut-off1 inside.
 - Sigma/pi separation: create .sg or .som (if both files exist, the first file sought is the .sg file).                                                                            
 - Only pi MOs: .sg/.som + .pi (command> touch wfnfile.pi).                                                                                                                         
 - Skip points with gauss. quad. weight < cut-off2: touch wfnfile.skg                                                                                                               
 - Skip points with ref. dens < cut-off3: touch wfnfile.skr.                                                                                                                        
 - cut-off2 and cut-off3 can be changed by editing the values in the previously created files.
 - Enable parallelization: Max         -> touch wfnfile.proc                                                                                                                          
                           Set # procs -> edit  wfnfile.proc                                                                                                                          
 - Atomic Overlap Matrix (AOM) calculation: touch wfnfile.aom                                                                                                                       
 - Change max dim. of array allocation: edit  wfnfile.mxal                                                                                                                          
   Include in this file (all lines are mandatory):                                                                                                                                  
   - Line 1: Max # atoms.                                                                                                                                                      
   - Line 2: Max #    mol. orbs.                                                                                                                                               
   - Line 3: Max # points.                                                                                                                                                     
   - Line 4: Max # prim. funcs.

## MIT license
See corresponding file called [`LICENSE`](LICENSE) for more details.

## Compile the code
See manual for detailed information.

Dependencies: A Fortran compiler

After cloning the code of the repository:

 $ mkdir build && cd build/
 $ cmake ..
 $ make

To force CMake to compile the code with your favourite compiler, prepend the cmake line 
   with e.g. FC=gfortran:

 $ cmake -DCMAKE\_Fortran\_COMPILER=ifort ..

