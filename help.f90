! ** ./help.f90 >> Help file of the fukui program
!
!  Copyright (c) 2022  Nicolás Otero Martínez - Marcos Mandado Alonso - Ricardo A. Mosquera Castro
!  This file is part of the fukui program available in:
!      https://github.com/nom05/fukui
!
!  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
!  by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!  
!  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License along with this code.  If not, see 
!  <http://www.gnu.org/licenses/>.


subroutine help

  character (len = 100) :: file1
  call getarg (0,file1)
  write (*,'( &
    & "Type: ",A," wfnfile stk_file/cube_file [# atom]"                                                            ,/,10X, &
    &   " * Extension for stk and wfn files are optional, extension for cube file is mandatory."                   ,/,10X, &
    &   " * stk_file can be either a FORTRAN binary file or a ordinary text file."                                 ,/,13X, &
    &      "They have to include, following this order: x, y, z, Gauss quadrature"                                 ,/,13X, &
    &      "weight and electron density values."                                                                   ,/,13X, &
    &      "If the ordinary text format is considered, then the format to read must"                               ,/,13X, &
    &      "be specified. the label free can be used to use free format specification."                            ,/,10X, &
    &    " * [# atom] is optional and enables gauss. func. skipping."                                              ,/,13X, &
    &      "If wfnfile.skd exists, it reads cut-off1 inside."                                                      ,/,10X, &
    &    " * Sigma/pi separation: create .sg or .som (if both files exist, the first file sought is the .sg file).",/,10X, &
    &    " * Only pi MOs: .sg/.som + .pi (command> touch wfnfile.pi)."                                             ,/,10X, &
    &    " * Skip points with gauss. quad. weight < cut-off2: touch wfnfile.skg"                                   ,/,10X, &
    &    " * Skip points with ref. dens < cut-off3: touch wfnfile.skr."                                            ,/,10X, &
    &    " * cut-off2 and cut-off3 can be changed by editing the values in the"                                    ,/,13X, &
    &      "previously created files."                                                                             ,/,10X, &
    &    " * Enable parallelization: Max        ->touch wfnfile.proc"                                              ,/,37X, &
    &                                "Set # procs->edit  wfnfile.proc"                                             ,/,10X, &
    &    " * Atomic Overlap Matrix (AOM) calculation: touch wfnfile.aom"                                           ,/,10X, &
    &    " * Change max dim. of array allocation: edit  wfnfile.mxal"                                              ,/,13X, &
    &      "Include in this file (all lines are mandatory):"                                                       ,/,17X, &
    &         " o Line 1: Max # atoms."                                                                            ,/,17X, &
    &         " o Line 2: Max #    mol. orbs."                                                                     ,/,17X, &
    &         " o Line 3: Max # points."                                                                           ,/,17X, &
    &         " o Line 4: Max # prim. funcs."                                                                      ,       &
    &                        /)') trim(file1)
end subroutine help
