# infn-magnetometer
INFN g-2 magnetometer analysis scripts

## Instructions

### gm2ita setup
When logging on gm2ita, please run these commands:
```
source setup_USER.sh
source /opt/cernroot/bin/thisroot.sh
```
In order to get the ROOT envinronment to run scripts

### Cloning repository
To get the code, please do:

Option 1) You'll need an assigned ssh key on github and you'll be able to commit/push your edits:
```
git clone git@github.com:giro94/infn-magnetometer.git
```
Option 2) Read-only mode:
```
git clone https://github.com/giro94/infn-magnetometer.git
```

### **Before analyzing any file, fix infinity symbols**
```
./fix_files_symbols.sh /path/to/input/folder/
```

### Inspect the traces in a file
Interactively:
```
root -l display_traces\(\"/path/to/input/file.csv\"\)
```
Save to output:
```
root -l -b -q display_traces\(\"/path/to/input/file.csv\"\,\"output.root\"\)
```

### HWP scan
```
cd HWP_scan_analysis
root -l -b -q plot_HWPscan\(\"/path/to/input/folder/\"\,\"\{-5,0,5,...,25\}\)
```
* The first argument is the path to the data folder.
* The second argument is the vector of HWP angles **chronologically** used during the scan.
* The code will create an output file named "folder_output.root" in the same folder you're running the script from.
* After the scan is complete, please write an "angles.txt" file containing one column listing the HWP angles. Therefore you'll be able to do the following instead:
```
root -l -b -q plot_HWPscan\(\"/path/to/input/folder/\"\)
```
and the code will read from the "angles.txt" file.


### Eddy currents analysis
```
cd Eddy_analysis
root -l -b -q analyze_eddycurrents.C\(\"/path/to/input/folder/\"\,\"output.root\"\)
```
This will average all the traces in the folder, align the kick in time, perform blumlein, SNR, baseline trends, and produce many useful histograms.

Optional arguments:
```
root -l -b -q analyze_eddycurrents.C\(\"/path/to/input/folder/\"\,\"output.root\"\,MAXFILES\)
```
To limit reading maximum MAXFILES files in the input folder

