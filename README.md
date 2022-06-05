<div align="center">
<img src="/image/MSRIDD_ICON.png" width="80px">
</div>

# MS-RIDD: Mass Spectrometry Radical Induced Dissociation Decipherer - A software for tandem mass spectra analysis of radical induced dissociation
The current program supports tandem MS/MS spectra analysis for oxygen attachment dissociation (OAD), which is one of the radical induced dissociation. The assembly is licensed under the CC-BY 4.0. Please send feedback, bug reports, and questions to Issues page.
***

## Requirements
Python version 3.7.3  
Packages included in requirements_pip.yaml are also required to be installed.

## Usage
MS-RIDD supports the export format of peak list and alignment result from MS-DIAL software.
For preparaing an input data, you need mzML file of OAD-MS/MS and text library that reflects the result of CID-MS/MS data.
To prepare a text library, please see the tutorial of MS-DIAL: https://mtbinfo-team.github.io/mtbinfo.github.io/

### Single analysis
1. Select "File > New project > Single analysis"
2. Select the project file folder
3. Select an input data (file must be text format (.txt) from MS-DIAL export)
4. Select data format (If a file is from peak list export, please select PeakList.)
5. Click start button

### Batch analysis
For a batch analysis. you need (a. Alignment result data, (b. PeakID data, (c. Peak list data of each sample.
All data can be exported from MS-DIAL.

1. Select "File > New project > Batch analysis"
2. Select the project file folder
3. Select an alignment data
4. Select PeakID data
5. Select peak list data of each sample
6. Click start button
***

## About data files in MS-RIDD project
After a data analysis, MS-RIDD generates several data files 
namely "analysis_table.pkl", "cid_result.pickle", "extracted_msms.pickle", "graph_info.pickle", "oad_result.pickle", and "structure_info.pickle".  
All these files should be in a same directory to open a MS-RIDD project again.
***

## Developers
Lead developer: Haruki Uchino (Keio University & RIKEN)  
Current main developer: Haruki Uchino (Keio University & RIKEN)  
Supervisor: Hiroshi Tsugawa (Tokyo University of Agriculture and Technology & RIKEN)
***

### Source code license
The source code is licensed under GNU LESSER GENERAL PUBLIC LICENSE (LGPL) version 3.  
See LGPL.txt for full text of the license.
***
