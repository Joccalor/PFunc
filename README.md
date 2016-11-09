PFunc
===

PFunc is a set of tools for fitting and analyzing function-valued traits. It was designed for quantifying preference functions (hence the name) in biological studies of sexual selection and evolution of mate preferences. PFunc accepts a bivariate dataset composed of stimuli and responses to those stimuli as its input. It then fits cubic splines to these data using the gam function in R, and it displays the curves, along with several useful measurements of the shape of the curve. Users have several options for adjusting the curves and the tools used to measure the curves, and for outputting results for further analysis.

Table of Contents
---

1. Files
2. Setup
 * Linux or Mac OS X
 * Windows
3. Usage
  * Running the PFunc GUI
  * Running PFunc from the R command line
4. Acknowledgements
5. License


Files
---
`PFunc_GUI.py` - the main file to run for the full PFunc GUI experience  
`PFunc_Rcode.r` - the supporting R code that fits the splines and extracts the useful metrics. This code *can* be run on its own in R without the GUI.  
`ExampleData_Horizontal.csv` - example data in one of two possible layouts  
`ExampleData_Vertical.csv` - example data in one of two possible layouts  
`icon.png` - icon file  
`PFuncPath.txt` - may help users who encounter path-related errors  
`COPYING` - the full GPLv3 license  
`README.md` and `README.txt` - README files with the same content in two different file formats.

Setup
---

Usage
---

Acknowledgements
---
Many thanks to Rafael Rodriguez, Kasey Fowler-Finn, Gerlinde Hoebel, David Gray, Darren Rebar, and Michael Reichert for their testing and feedback during development.

License
---
Copyright 2016 Joseph Kilmer

PFunc is distributed under the GNU General Public License v3. A full copy of the license is available in the accompanying file called COPYING.txt.

PFunc is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PFunc is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
