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
#### Linux or Mac OS X
1. Install and set up R  
  1. Visit <https://cran.r-project.org/> and follow the links to download the latest version of R for your operating system. Install R as you would a normal program. Note for Linux users, you may instead install R from the command line.  
  2. Install the mgcv package in R by opening R and entering the following command: `install.packages("mgcv")`. You will be prompted to select a mirror from which to download the package. Any option will work, but it is typically best to select the option closest to your location. Follow any instructions R gives you for installing the package.  
2. Install and set up Python  
  1. Visit <https://www.python.org/> and follow the links to download the latest version of Python 3.x for your operating system. Install python as you would a normal program. Note for Linux users, you may instead install Python from the command line.  
  2. Install the necessary Python libraries. Open up the terminal (if you are not sure how, run a search for "terminal" on your computer. Alternatively, Mac users will find a link to the terminal in their Applications > Utilities folder). Enter the following commands:  
    `pip3 install rpy2`  
    `pip3 install matplotlib`  
  Follow any instructions that either one gives you to complete the setup.
3. Running PFunc


#### All users
1. Install and set up R  
  1. Visit <https://cran.r-project.org/> and follow the links to download the latest version of R for your operating system. Install R as you would a normal program.  
  **Note for Windows users:** The setup wizard will give you the option of installing the 32-bit files or the 64-bit files. Install both.  
  **Note for Linux users:** You may instead install R from the command line.  
  2. Install the mgcv package in R by opening R and entering the following command: `install.packages("mgcv")`. You will be prompted to select a mirror from which to download the package. Any option will work, but it is typically best to select the option closest to your location. Follow any instructions R gives you for installing the package.  
2. Install and set up Python  
  1. Visit <https://www.python.org/> and follow the links to download the latest version of Python 3.x for your operating system. Install python as you would a normal program.  
  **Note for Windows users:** the first screen of the installation wizard will ask if you want to add Python to PATH. Make sure to select this option.  
  **Note for Linux users:** you may instead install Python from the command line.  
  2. Install the matplotlib Python library.  
  Open up the terminal (if you are not sure how, run a search for "terminal" or "command prompt" on your computer). Enter the following command:  
    `pip3 install matplotlib`  
  Follow any instructions to complete the setup.
  3. Install the rpy2 Python library.  
  **For Linux and Mac users:** Enter the following command in the terminal, just like with the previous step:  
  `pip3 install rpy2`  
  **For Windows users:** Go to <http://www.lfd.uci.edu/~gohlke/pythonlibs/#rpy2> and download one of the rpy2.whl files. Start by trying one of the ones at the end of the list (`rpy2-2.8.3-cp35-cp35m-win32.whl` worked for me). If the following steps don't work with the one you tried, then try another.  
  * Download the rpy2.whl file, and take note of the directory path where you save it.
3. Running PFunc


#### Windows


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
