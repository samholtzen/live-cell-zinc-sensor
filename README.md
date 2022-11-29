# Live-Cell Zinc Sensor Pipeline
##### SCRIPTS USED IN PUBLICATION “HUMAN CELLS EXPERIENCE A ZN2+ IN EARLY G1” BY RAKSHIT & HOLTZEN ET AL
## Introduction
### Long-term live-cell free zinc measurements
Many groups have done several population-level analyses of zinc throughout the cell cycle, but very few if any have published on long-term live-cell zinc measurements in asynchronously cycling cells. This pipeline will aid in analyzing long-term live-cell imaging of cells expressing the Zap suite of genetically encoded FRET sensors for measuring free zinc.
### EllipTrack
This pipeline is formatted for use with the EllipTrack output from a computing cluster, but you can also use it with EllipTrack run on your local computer. The outputs are slightly different between each, so make sure you familiarize yourself with which you should expect and how to treat it.
#### From a local computer:
EllipTrack run on a local computer will have one folder for all the segmentation parameters and track visualization movies (if you chose to save these). Additionally, in the results folder, it will have one of each of the following files:
-	Jitter Correction – This file describes the movement in the x and y directions of an image between frames.
-	Segmentation – This file stores the properties of each segmented ellipse, including parameters about its size and shape.
-	Probabilities – This file stores the predictions made by the machine learning algorithm to classify each ellipse.
-	Tracks – This file stores the information about linking tracks together between frames in order to follow a cell as it moves and divides.
-	Signals – This file stores the extracted fluorescent signals and physical properties of the cell, such as size and shape. We will use this file in our analyses.

Each .mat file contains a cell of data structures in the shape of the 96-well plate you processed. For example, if you processed movies from wells A1 to D12, the outputs will all be of shape 4x12. As such, it’s very important to remember which columns or rows you treated with which treatments, so you can back-calculate which elements of the signals file correspond to what treatments. We suggest keeping a .csv file of your treatment conditions in each well of the plate so you can read it into the pipeline programmatically.
#### From the computing cluster:
EllipTrack run on the computing cluster will have one folder that contains all the output files for one movie. Putting these individual files in a single file is hard, so there is a script that takes all the native file structure and places it into a single storage structure.
## The Pipeline
The actual pipeline is set up in discrete steps that can be run modularly or together using the main_file.m script. A short description of each script and how the scripts interact is listed below. Refer to the figure to track the flow of information from one script to another.
### The Main File
The main_file.m is the wrapper file that controls the pipeline. This includes all of the parameters for the analysis. The file will read in a struct_cell file and extract information from the file name, including the date of the experiment that you put in the filename and the cell type.

There are several switch/case expressions that will change the parameters of analysis based on the name of your input. For example, if you marked the struct_cell file with the cell type of ‘scr’, the switch/case expression will check if the cell type matches any known cases, and if it does, it will assign parameters for color, treatment conditions, and what conditions you want to plot. If your cell type does not show up here, you may add more cases as needed.

The last parameter that is set before running the scripts is the FRET sensor filter. The dynamic range of the ZapCV2 sensor lies somewhere between 3 and 8 units, which means that anything beyond that is most likely aberrantly expressing cells; however, if you have a validated sensor that has a different dynamic range, change the parameters to the minimum and maximum value allowed by your sensor.
### Running the Pipeline
#### Step 1: Combining individual folders
To combine the subfolders and subfiles into one structure, you should use the combine_folders.m file. If you used EllipTrack locally instead of on the computing cluster, you can skip this step. Below are the instructions for processing an entire directory copied from the computing cluster.
1.	Open the combine_folders.m file in MATLAB.
2.	Click “Run” or type the following into the command window

`combine_folders`

3.	A user interface will open. Navigate to the directory where the individual folders are listed. The folder should be labeled results/ and will contain all folders in the format ‘1_1_1’ etc.
4.	Press OK.
5.	Wait. A full data structure corresponding to a 96-well plate could be upwards of 5 GB or more, so it will take MATLAB some time to load and unload the individual files.
6.	Once it is done, check the results/ folder and make sure there is a file labeled ‘signals.mat’. This is your combined signals file that you’ll be using for the next steps.
##### Step 2: Extracting signals
Once you have your ‘signals.mat’ file, you will need to extract and condense the signals. The ‘signals.mat’ file is absurdly large, and it contains far more information than we need to do our analysis. As such, there is a script that processes the file and extracts the desired information.
### Subfunctions
#### link_mother_daughter.m
This file takes in a cell array containing a data structure from the ‘signals.mat’ file. It works on an individual cell array, or a cell array with multiple elements. It works backwards from the end towards the beginning and links mother and daughter cells together to create a lineage and link all the signals to each other for all of the properties you are interested in. For example, it will link the nuclear intensity, FRET, and CFP intensities of an entire lineage and store it in an array.
#### get_proliferation.m
This file takes in a cell array containing a data structure from the ‘signals.mat’ file. It works on an individual cell array, or a cell array with multiple elements. It works by cycling through all cell tracks in a well and adds up the number of ellipses at every frame. It then divides the entire vector by the first element to get a normalized cell count.

Below are the steps to take to extract the information:
1.	Open the export_lineage_struct.m file in MATLAB.
2.	Press “Run” or type the following into the command window

`combine_folders`

3.	A user interface will open. Navigate to the ‘signals.mat’ file that you just created in the previous step and click on it.
4.	Press OK.
5.	In the command window, the following prompt will be printed:

`What cell type is this? Input a string.`

6.	Type in your cell type surrounded by single quotes. For example, if this structure corresponds to a scrambled control, type in ‘scr’ and press Enter.
7.	In the command window the following prompt will be printed:

`What date was this experiment conducted on? Input a number in the format YYYMMDD.`

8.	Enter a number that corresponds to the date you conducted this experiment and press enter.
9.	Wait. Line 10 will take several minutes, as the ‘signals.mat’ file is very large.
10.	Once it is done, open the output directory and ensure that your file saved in the correct location, with the correct name.
#### Step 3: Analyzing Data
##### Changing Run Parameters
1.	Open main_file.m.
2.	Go to lines 26-45. Scan through the cases in the switch/case method and see if your cell type is present. If it is, continue to step 5. If not, do the following.
3.	Add another line after the last case block and type out the word “case”. Add a space, then type in your cell type in single quotes as it matches what you put in the export_lineage_struct.m file.
4.	Fill in three variables. The first two should be the same length as the number of columns in your experiment, and the third should contain only integers between 1 and the number of columns in your experiment:
5.	The colors you want to use on the plots (colors_cell). Should be a cell array of strings the same length as the number of columns in your experiment. This will take a MATLAB color name, a MATLAB short color name, an RGB triplet, or a Hexadecimal Code.
6.	The treatment conditions you used in each column of the experiment (condition_cell). Should be a cell array of strings the same length as the number of columns in your experiment.
7.	The column numbers of the conditions you want to plot (conditions_to_plot). Should be an array of integers between 1 and the number of columns in your experiment.
8.	Go to lines 60-61. Change these numbers to the minimum and maximum FRET values you expect for your zinc sensor. It will change depend on the sensor, but a rule of thumb for ZapCV2 is a FRET value between 3 and 8.
##### Running The Pipeline
1.	Click the Run button at the top. MATLAB will open a user interface for you to select one or more struct_cell files to process.
2.	Use Ctrl-Click or Command-Click to highlight the ones you want to process and Press OK.
3.	As each step is finished, you will briefly see the plots show up on the screen, then close. The plots and FRET features CSV files are saved in the output folder set on lines 51-54.
