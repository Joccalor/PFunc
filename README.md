PFunc
===

PFunc is a set of tools for fitting and analyzing function-valued traits. It was designed for quantifying preference functions (hence the name) in biological studies of sexual selection and evolution of mate preferences. PFunc accepts a bivariate dataset composed of stimuli and responses to those stimuli as its input. It then fits cubic splines to these data using the gam function in R, and it displays the curves, along with several useful measurements of the shape of the curve. Users have several options for adjusting the curves and the tools used to measure the curves, and for outputting results for further analysis.

Table of Contents
---

1. Files
2. Setup
3. Usage
4. Troubleshooting
5. Acknowledgements
6. Contact
7. License

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

#### 1. Install and set up R  
  1. Visit <https://cran.r-project.org/> and follow the links to download the latest version of R for your operating system. Install R as you would a normal program.  
  **Note for Windows users:** The setup wizard will give you the option of installing the 32-bit files or the 64-bit files. Install both.  
  **Note for Linux users:** You may instead install R from the command line.  

  2. Install the mgcv package in R by opening R and entering the following command: `install.packages("mgcv")`. You will be prompted to select a mirror from which to download the package. Select an option that is relatively close to your location. Follow any instructions R gives you for installing the package.  

#### 2. Install and set up Python  
  1. Visit <https://www.python.org/> and follow the links to download the latest version of Python 3.x for your operating system. Install python as you would a normal program.  
  **Note for Windows users:** the first screen of the installation wizard will ask if you want to add Python to PATH. Make sure to select this option.  
  **Note for Linux users:** you may instead install Python from the command line.  

  2. Install the matplotlib Python library. To do this, open up the terminal (if you are not sure how, run a search for "terminal" or "command prompt" on your computer).  
  **For Linux and Mac users:** enter this command `pip3 install matplotlib`  
  **For Windows users:** enter this command `pip install matplotlib`  
  Follow any instructions to complete the setup.

  3. Install the rpy2 Python library.  
  **For Linux and Mac users:** Enter the following command in the terminal, just like with the previous step: `pip3 install rpy2`  
  **For Windows users:** Go to <http://www.lfd.uci.edu/~gohlke/pythonlibs/#rpy2> and download one of the rpy2.whl files. Start by trying one of the ones at the end of the list (`rpy2-2.8.3-cp35-cp35m-win32.whl` worked for me). If the following steps don't work with the one you tried, then try another.

    * Download the rpy2.whl file, and take note of the directory path where you save it.  
    * Enter the following command in the command prompt (like in the previous step), substituting the directory path below for the one where you saved the rpy2.whl file:  
  `pip install “C:\Users\Joey\Downloads\rpy2-2.8.3-cp35-cp35m-win32.whl”`


Usage
---
#### Data input  
PFunc is expecting data with two main components: a set of **stimuli** (which are plotted along the *x*-axis) and a set of **responses** to those stimulus values (which are plotted on the *y*-axis).

You have the option of formatting your data either in a **horizontal** layout (as in the `ExampleData_Horizontal.csv` file) or a **vertical** layout (as in the `ExampleData_Vertical.csv` file).

Regardless of which layout option you choose, you must save your data as a **.csv** file. PFunc does not accept other file formats.

#### Running the PFunc GUI

#### Running PFunc from the R command line

Troubleshooting
---


Acknowledgements
---
Many thanks to Rafael Rodriguez, Kasey Fowler-Finn, Gerlinde Hoebel, David Gray, Darren Rebar, and Michael Reichert for their testing and feedback during development.

Contact
---
For comments or questions, contact Joey Kilmer at jtkilmer@uwm.edu  
Get the latest version of PFunc at <https://github.com/Joccalor/PFunc>

License
---
Copyright (C) 2016 Joseph Kilmer

PFunc is distributed under the GNU General Public License v3. A full copy of the license is available in the accompanying file called COPYING.txt.

PFunc is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PFunc is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
