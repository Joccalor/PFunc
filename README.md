PFunc
===
PFunc is a set of tools for fitting and analyzing function-valued traits. It was designed for quantifying preference functions (hence the name) in biological studies of sexual selection and evolution of mate preferences. PFunc accepts a bivariate dataset composed of stimuli and responses to those stimuli as its input. It then fits cubic splines to these data using the gam function in R, and it displays the curves, along with several useful measurements of the shape of the curve. Users have several options for adjusting the curves and the tools used to measure the curves, and for outputting results for further analysis.

Table of Contents
---
1. Files
2. Setup
  * Install and Set Up R
  * Install and Set Up Python
3. Usage
  * Data Input
  * Running the PFunc GUI
    * Startup
    * The Interface
    * Settings
    * Group-Level Splines
    * Output
  * Running PFunc from the R Command Line
    * Startup
    * Arguments of PFunc
    * Examples
4. Acknowledgements
5. Contact
6. License

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
Follow these steps to setup PFunc on Windows, Mac OS X or Linux systems.

#### 1. Install and Set Up R  
  1. Visit <https://cran.r-project.org/> and follow the links to download the latest version of R for your operating system. Install R as you would a normal program.  
  **Note for Windows users:** The setup wizard will give you the option of installing the 32-bit files or the 64-bit files. Install both.  
  **Note for Linux users:** You may instead install R from the command line.  

  2. Install the mgcv package in R by opening R and entering the following command: `install.packages("mgcv")`. You will be prompted to select a mirror from which to download the package. Select an option that is relatively close to your location. Follow any instructions R gives you for installing the package.  

#### 2. Install and Set Up Python  
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
  `pip install "C:\Users\joey\Downloads\rpy2-2.8.3-cp35-cp35m-win32.whl"`

Usage
---
Below is a bare-bones overview of using PFunc. For more details on usage as well as general theory underlying the study of preference functions, check out our paper [*add citation when paper is published*].

### Data Input  
PFunc is expecting data with two main components: a set of **stimuli** (which are plotted along the *x*-axis) and a set of **responses** to those stimulus values (which are plotted on the *y*-axis).

You have the option of formatting your data in one of two ways:

* In the **horizontal** layout, the first column of data contains the stimulus values, and each subsequent column contains the responses of an individual to those stimulus values. See the `ExampleData_Horizontal.csv` file.
* In the **vertical** layout, data are organized into a minimum of three columns: one with individual identities, one with stimulus values presented to those individuals, and one with the responses of those individuals to those stimuli. See the `ExampleData_Vertical.csv` file.

Regardless of which layout option you choose, you must save your data as a **.csv** file.

### Running the PFunc GUI  
#### Startup  
You have several different options for running PFunc. Make sure that `PFunc_GUI.py` and `PFunc_Rcode.r` are together in the same directory. If you encounter an error, see the Troubleshooting section below.
 * **Option 1: Run PFunc through IDLE**  
 The most straightforward way to run PFunc is through IDLE. Right-click (or the equivalent on a Mac) on the `PFunc_GUI.py` file, and choose to "Open With..." IDLE (or "Edit With..." on Windows). This will start up IDLE with the code for PFunc displayed in a window. Click on the window of colorful code to make sure it is the current active window, and then run the code either by pressing F5 or by clicking Run > Run Module in the menu at the top bar.
 * **Option 2: Run PFunc from the command line**  
   * Take note of the directory where PFunc is saved.
   * Open up an instance of your terminal (command prompt in Windows) and navigate to the PFunc directory using the `cd` command. Examples: `cd "C:\Users\joey\PFunc"` in Windows; `cd "/home/joey/PFunc"` in Linux; `cd "/Users/joey/PFunc"` in Mac OS X.
   * Then enter this command for Mac or Linux: `python3 PFunc_GUI.py`, or this command for Windows: `python PFunc_GUI.py`.
 * **Troubleshooting**  
 You might encounter an error while trying to run PFunc if rpy2 has trouble finding R on your computer. One sign of this error is if the error message mentions R_HOME or PATH. To solve this, you'll have to find where R is installed on your computer and follow the instructions in the `PFuncPath.txt` file that accompanies this program.

#### The Interface  
Once you open your data file, PFunc will display a page of graphs, each one showing data from a single individual with a spline fit through the points. A page displays up to nine individuals, and you can navigate to different pages with the controls at the bottom of the screen.

You can select a graph by clicking on it, and you can enlarge a graph by double-clicking on it.

When a graph is selected, its smoothing parameter is displayed in the Smoothing box on the control panel on the right-hand side of the window, and several important metrics are displayed in the Summary box right below the Smoothing box.

##### *Smoothing Parameter*
The smoothing parameter controls how stiff or how wobbly a spline is. A low smoothing parameter yields a wobbly spline, and at extremely low values, the spline runs through every data point. A high smoothing parameter yields a spline that is a gentle curve, and at extremely high values, the spline is a straight line. The gam function in R (which fits the splines) chooses an optimal smoothing parameter for you. In most cases, it is not necessary to alter the smoothing parameter. However, should you wish to change it, you have two main options:

* **Specify your own smoothing parameter.** You can use the "-" and "+" buttons next to the smoothing parameter to increase and decrease the value. You may also type in your own smoothing parameter and set it by pressing the Enter key. When a smoothing parameter has been altered, the small magenta dots in the bottom left corner of the graphs change to cyan. If at any point you want to revert to the default smoothing parameter, simply press the Reset button in the Smoothing box (this will also change the cyan dot back to magenta).  
* **Adjust the limits imposed upon smoothing parameters.** PFunc was built with the philosophy that preferences function splines should be neither too wobbly nor too stiff. Therefore, by default, PFunc restricts all smoothing parameter values to be between 0.05 and 5. We found this to be a reasonable range for the types of data we were working with, but you may find a different range to be more suitable for your data. You can change this range in the SP Limits setting (described below). Alternatively, you may decide that you want no limits on your smoothing parameters, and instead want to see what default the program gives you, no matter how big or how small. You have the option of turning off smoothing parameter limits altogether (also described in the SP Limits setting below). When you make changes to the SP Limits setting, only the graphs with magenta (and not cyan) dots will change. This prevents you from accidentally wiping away any intentional changes to specific smoothing parameters.

##### *Summary*
The summary window displays several useful measurements of the spline in the selected graph:

* **Peak Stimulus:** This is the x-axis value that corresponds to the highest point of the curve, as illustrated by the red vertical lines on the graphs.
* **Peak Height:** This is the y-axis value that corresponds to the highest point of the curve.
* **Tolerance:** This is the width of the curve at a certain height. By default it is measured 1/3 of the way down from the peak of the curve. See the Tolerance settings below to change the proportion of the drop or to measure tolerance at a set y-axis value for all curves.
* **HD Strength:** This is height-dependent strength. Strength is a measure of how far the curve drops away from the peak. This measure is height-dependent because the same absolute drop in a curve will be a larger proportion of a curve that is close to the baseline, compared to a curve that is high up. It is calculated as the coefficient of variation, which is the standard deviation divided by the mean.
* **HI Strength:** This is height-independent strength. Strength is a measure of how far the curve drops away from the peak. This measure is height-independent because the same absolute drop in a curve will be the same, regardless of how high up the curve is. It is calculated as the standard deviation divided by the range.
* **Responsiveness:** This is the average height of the curve.
* **Smoothing:** This displays the smoothing parameter used to generate the curve with the above values. Most of the time it is the same value as is in the Smoothing box.

##### *Messages*
This is a box for PFunc to display messages to you. The messages might explicitly confirm an action (for example, when you open a new file, PFunc wipes away the previous smoothing values that you set, and the Messages box tells you this so you know you're starting with a clean slate), or they might give the user a warning. Messages are accompanied by a timestamp.

#### Settings  
You have access to a variety of settings that give you control over how your data are analyzed. Below are descriptions of the available settings. PFunc will always start up with the default settings, but you can save any changes you make to the settings through the File menu ("Save Current Settings"), and then load them again later ("Load Previous Settings"). You can also revert to the default settings ("Restore Default Settings").

* **View** Here you can toggle the visibility of 5 different things on or off, depending on how you want to view your graphs.
* **SP Limits** Here you can adjust the limits to the smoothing parameters that are chosen for your splines. By default, PFunc limits them to valuse between 0.05 and 5, but you have full control over the minimum and maximum values here, and you also have the option to turn these limits off entirely with the checkbox.
* **Tolerance** This gives you fine control over how tolerance (the width of the curve at a given height) is measured.
  * First you have the option of measuring tolerance relative to the height of the peak (e.g. 1/3 of the way down from the top), or measuring it at a set y-axis value. To switch between these two options, click one of the two radio buttons next to the appropriate settings: the radio button next to "Drop from peak" is for the relative option, and the radio button next to "At set value" is for the absolute option.
  * Once you've selected *how* tolerance should be measured, you can control *where* it is measured. If you've chosen the relative option, you can control the drop from the peak. This takes the distance from the peak of the curve to the floor (which you can also set), and measures tolerance at a particular proportion of that distance down. If, on the other hand, you've chosen to measure tolerance at an absolute height, you then can specify that exact height with the "At set value" setting.
  * Finally, you can decide whether you want a broad measure of tolerance or a strict one. The difference comes into play when there are secondary peaks in your curves. If you choose the broad option, then tolerance will be measured as the total width of the curve under the primary peak plus any relevant secondary peaks. If you choose the strict option, then tolerance will be measured as the width of the curve *only* under the primary peak.
* **Peak** You may encounter some curves that have a peak in the middle and then a high point at a value at the edge of your x-axis range. If the height of this edge portion is higher than the central peak, then PFunc will naturally call it the peak. But if you want information on that central peak instead, you can adjust this Peak Within value. What it does is it searches for a peak within a certain central proportion of the curve, and if it finds one, it uses that. For example, if you set Within to equal 0.7, PFunc will first search for a peak in the middle 70% of the curve.

#### Group-Level Splines
You can combine individual splines to form group-level splines. This can come in handy if you want to generate splines at the replicate-, treatment-, family-, population-, or species-level. Group splines can even be combined to form higher-order group splines.

To do this, first go to Advanced > Construct Group-Level Spline..., then choose which individuals belong in the group (use the Shift and Ctrl keys to select multiple individuals). You can decide whether you want the mean or median of these individual splines, and you can choose to name the new group-level spline. Once you are finished, press Okay, and your new group-level spline will be added alongside your other splines.

Note that this does not affect your input data file; if you want to retain these values, you'll need to output them (see below). Also note that group-level splines may be best fit with lower smoothing parameters than individual-level splines.

#### Output  
Here are descriptions of the various output options available under the File menu.

* Output Spline Figures: Creates a pdf with all the graphs using the current smoothing parameters.
* Output Spline Summaries: Creates a spreadsheet that contains all of the information in the Summary window for every individual.
* Output Spline Points: PFunc extracts the y-values at 200 evenly-spaced points along the curve, and it saves these as a spreadsheet. This is useful if you want to plot your curves in a different program.
* Output Tolerance Points: Creates a spreadsheet containing all of the x-axis values that correspond to the upper and lower limits of tolerance--that is, the start and stop points of the horizontal blue lines in the graphs.

### Running PFunc from the R Command Line
If you are comfortable working in the R command line environment, you may use Pfunc without the GUI. Note that when PFunc is used this way, data **must** be set up in the horizontal format (as in `ExampleData_Horizontal.csv`), never in the vertical format (as in `ExampleData_Vertical.csv`).

#### Startup
Open R. Set the working directory to wherever the `PFunc_Rcode.r` file is saved on your computer.  
* One way to do this is by finding the directory path and entering the command `setwd(directory/path)`. You can check your current working directory with the command `getwd()`.
* There is another way to do this if you are not running R directly from your terminal. Depending on your version of R, you will find the option either under "File > Change dir..." or "Misc > Change Working Directory..."
* Load the script into R with this command `source(PFunc_Rcode.r)`
* Load your data into R with a command like `mydata <- read.csv(datafile.csv)`
* Now you can run the default analysis on your data with the command `PFunc(mydata)`. See below for a full set of options.  

#### Arguments of PFunc
Here is a list of all the arguments that you can specify when calling PFunc directly from R.  

* `input.data` - the name of your dataset in the R environment (no default)

* `diagnose.col` - an optional argument used only when assessing individual functions. To see the output for a single individual without creating a new file, set this number equal to the corresponding column in your datasheet (default = 0).

* `diagnose.sp` - an optional argument used only when assessing individual functions. By setting this to a positive number, you can specify the smoothing parameter of individual functions that you are examining (default = -1).

* `auto.sort` - an option to sort the data by increasing values of the stimulus variable (default = TRUE).

* `blacklist` - if there are columns of your data that you wish to exclude from analysis, place their numbers in this argument with the concatenate `c()` function. For example, to exclude the third and fourth columns, set `blacklist = c(3, 4)` (default = NULL).

* `whitelist` - if there are only a few specific columns of your data that you wish to analyze and exclude the rest, you can place the included columns' numbers in the whitelist with the concatenate `c()` function. For example, to include only the third and fourth columns in analysis, set `whitelist = c(3, 4)` (default = NULL).

* `sp.binding` - when true, it constrains the auto-generated smoothing parameters to be between `min.sp` and `max.sp` (see below) (default = TRUE).

* `sp.assign` - allows you to specify your own smoothing parameters for individual columns rather than letting the script generate optimized values. Must be in the form of a two-column matrix in which the first column contains the column numbers for which you want to specify smoothing parameters. The second column contains the corresponding desired smoothing parameters. For example, `new.sp <- matrix(c(2, 1, 3, .5, 4, 0), ncol=2, byrow=T)`. Here, the individual in input.data column 2 gets an sp of 1, the individual in input.data column 3 gets an sp of .5, and the one in column 4 gets an sp of 0. To run with these values, set `sp.assign = new.sp` (default = 0).

* `max.sp` - if you are constraining auto-generated smoothing parameters with the `sp.binding` argument (above), `max.sp` is the upper limit (default = 5).

* `min.sp` - if you are constraining auto-generated smoothing parameters with the `sp.binding` argument (above), `min.sp` is the lower limit (default = 0.05).

* `peak.within` - sets the window in which the script will first look for the peak of a closed function. If no peak is found, it will then expand its search outward. The default value is 0.7, meaning that it will first look inside of the middle 70% of the function.

* `summary.out` - filename for the output summary that includes peak, tolerance, and other measures of each function. If FALSE, no file will be made (default = "spline_summaries.csv").

* `graph.out` - filename for the output graphs for all analyzed functions. If FALSE, no file will be made (default = "spline_graphs.pdf").

* `points.out` - filename for the output points that make up each function. Useful for plotting functions in external programs. If FALSE, no file will be made. You can specify the number of points you want to output for each function in the `predictions` argument below (default = FALSE).

* `max.y` - for the output graphs in `graph.out` (above), this argument specifies the height of the vertical axis. Leave at default value of 0 for automatic.

* `pdf.row` - for the output graphs in `graph.out` (above), this argument specifies the number of rows of graphs each page of the pdf will have (default = 4).

* `pdf.col` - for the output graphs in `graph.out` (above), this argument specifies the number of columns of graphs each page of the pdf will have (default = 3).

* `graph.points` - an option to turn data points on or off in the output graphs (default = TRUE).

* `graph.peak` - an option to turn on or off a red vertical line at the peak of the function (default = TRUE).

* `graph.tol` - an option to turn on or off a blue horizontal line showing the width calculated for tolerance (default = TRUE).

* `graph.sp` - an option to print the smoothing parameter of each function above its corresponding graph (default = FALSE).

* `graph.se` - an option to plot the standard error around the spline (default = FALSE).

* `drop` - specify how far down from the peak of the function that tolerance (the width) should be measured (default = 1/3).

* `tol.mode` - how would you like the script to calculate tolerance? The default value "norm" looks at the width of the function under the main peak plus any other areas at that same y value that may be under secondary peaks. The alternative value "strict" only counts the width of the part of the function that contains the main peak and ignores other peaks, even if they rise above the chosen y value.

* `tol.floor` - this specifies the baseline of the data. It is used in conjunction with `drop` to determine the height where tolerance should be measured. For example, if `drop` equals 1/3, and `floor` equals 0, then tolerance will be measured at 1/3 of the distance from the peak to 0 (default = 0).

* `predictions` - if you want to extract points along your splines using `points.out`, this argument tells the script how many points you would like. The default value 0 outputs points only at the input stimulus values. Anything greater than 0 will output that many points evenly spaced along the x-axis.

* `ghost` - if you are specifying your own smoothing parameters with the `assign.sp` argument, and if you are creating an output graph file with the `graph.out` argument, then setting ghost to TRUE will plot splines with the original default smoothing parameters in gray alongside your new splines. This is useful for comparing the effects of your altered smoothing parameters (default = FALSE).

* `allfromsplines` - an option to specify how you would like strength and responsiveness to be calculated. With the default value of TRUE, they will be calculated from the splines. Change this to FALSE only if you want to calculate strength and responsiveness from your input data points.

#### Examples
The following examples assume that your data file is called "mydata" in the R environment.

* To output spline summaries and a pdf of the graphs using automatic smoothing parameters bound between 0.05 and 5, use this command `PFunc(mydata)`

* To look at an individual spline (say, the one in the third column of your spreadsheet) without outputting any files, use this command `PFunc(mydata, 3)`

* To see what that individual spline would look like with a new smoothing parameter (say, 0.7), use this command `PFunc(mydata, 3, 0.7)`

* To keep this change for the final output, start by constructing a matrix as described in `sp.assign` above. Then set `sp.assign` to the name of the matrix, like this:  
```
new.sp <- matrix(c(3, 0.7), ncol=2, byrow=T)
PFunc(mydata, sp.assign=new.sp)
```

* By default, the height at which tolerance is measured is a proportion of the distance between the peak of a curve and the value of the `floor` argument. If you prefer to measure tolerance at the same y-axis height across all individuals (say, 2.5), regardless of the height of the peak, set `floor` equal to your desired height, and set `drop` equal to 1, like this `PFunc(mydata, floor=2.5, drop=1)`.

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
