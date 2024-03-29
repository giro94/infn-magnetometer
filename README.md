# infn-magnetometer
INFN g-2 magnetometer analysis scripts

#### Table of Contents  
[General instructions](#general-instructionss)  
* [gm2ita](#gm2ita-setup)
* [Cloning repository](#cloning-repository)
* [Fix symbols](#before-analyzing-any-file-fix-infinity-symbols)
* [Note on running scripts](#note-on-running-scripts)
  
[Script instructions](#script-instructions)  
* [Inspect a file](#inspect-the-traces-in-a-file)
* [Average many files](#inspect-the-average-traces-in-a-folder)
* [HWP scan](#hwp-scan)
* [Magnet ramp](#magnet-ramp-analysis)
* [Eddy currents](#eddy-currents-analysis)

[Build your own code](#build-your-own-code)

## General instructions

### gm2ita setup
When logging on gm2ita, please run these commands in order to get the ROOT envinronment to run scripts:
```
source setup_USER.sh
source /opt/cernroot/bin/thisroot.sh
```

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

### Note on running scripts
Every time you want to run a ROOT macro, you have two ways of giving arguments:  
1. escape each special character: `root -l -b -q macro.C\(\"string_argument\"\,\{array,of,values\}\,1000\)`
2. put single quotes: `root -l -b -q 'macro.C("string_argument",{2,5,4},1000)'`

In the following examples, I will always use option 1.

## Script instructions

### Inspect the traces in a file
Interactively:
```
root -l display_traces.C\(\"/path/to/input/file.csv\"\)
```
or, save to output:
```
root -l -b -q display_traces.C\(\"/path/to/input/file.csv\"\,\"output.root\"\)
```


### Inspect the average traces in a folder
Interactively:
```
root -l average_traces.C\(\"/path/to/input/folder/\"\)
```
or, save to output:
```
root -l -b -q average_traces.C\(\"/path/to/input/folder/\"\,\"output.root\"\)
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

### Magnet ramp analysis
```
cd Ramp_up_analysis
root -l -b -q plot_rampup.C\(\"/path/to/input/folder/\"\,\"output.root\"\,\"magnet_current.csv\"\)
```
The magnet current is a csv that can be retrieved from SEEQ (need special authorization, ask Matt).
  The downloaded file is actually a .xlsx, and you'll need to export the 2nd tab in csv format.
You don't need to provide the third argument. If missing, the script will plot the measured B-A vs timestamps.


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


## Build your own code
The general purpose header `analysis_tools.C` contains useful functions to help you read all the data.  
A minimal working code could be:
```
#include "analysis_tools.C"

vector<TString> files = getListOfFiles(folder);
int Nfiles = files.size();

//Reading first file to get trace info
pair<int,map<TString,int>> headers = getFileLengthAndHeaders(Form("%s/%s",folder.Data(),files[0].Data()));
int Nlines = headers.first;
map<TString,int> map_varnames = headers.second;
int Nvars = map_varnames.size();

for (int fi=0; fi<Nfiles; fi++){
		
  TString fname = files[fi];
  TDatime datetime = getFileTime(fname);
  TString filepath = Form("%s/%s",folder.Data(),fname.Data());

  vector<vector<double>> traces = readFileTraces(filepath,Nvars);
  vector<double> trace_time = traces[map_varnames["Time"]];
  vector<double> trace_A = traces[map_varnames["Channel A"]];
  vector<double> trace_B = traces[map_varnames["Channel B"]];
  vector<double> trace_C = traces[map_varnames["Channel C"]];
  vector<double> trace_avgC = traces[map_varnames["average(C)"]];
        
  /* .... */
}
```



