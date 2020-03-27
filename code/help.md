# Help

Please do not hesitate to contact the team with questions (see contact information listed in the main tab)

## Contact Rates data
The model could be adjusted for different schools and institutions by loading a new contact matrix

A square matrix in CSV format for the cohorts.
* all entries in [0,1]
* Diagonals should normally be largest
* Columns and rows indicating the cohort names

For example, a school with 7 grades may have:
  
 |G1|G2|G3|G4|G5|G6|G7
---|---|---|---|---|---|---|---
G1| 1|0.05|0.05|0.05|0.05|0.05|0.05
G2| 0.05|1|0.3|0.1|0.1|0.1|0.1
G3| 0.05|0.3|1|0.1|0.1|0.1|0.1
G4| 0.05|0.1|0.1|1|0.3|0.3|0.2
G5| 0.05|0.1|0.1|0.3|1|0.3|0.3
G6| 0.05|0.1|0.1|0.3|0.3|1|0.3
G7| 0.05|0.1|0.1|0.2|0.3|0.3|1

## Symptoms data
A table in CSV format for the symptom scores
* CSV format with Header: Time,ViralShedding,TotalSymptomScore,SystemicSymptoms
* Each row -> day since infection. 

Here is the data from Carrat et al. using in modeling influenza (day from infection):
 
 Time   | ViralShedding | TotalSymptomScore | SystemicSymptoms | RespiratorySymptoms | NasalSymptoms
|---|---|---|---|---|---|
1|1.88648930|0.251207729|0.115168539326|0.145505617978|0.180898876404
2|2.99897529|0.666666666|0.934269662921|0.661235955056|0.792696629213
3|2.63071806|0.850241545|0.747191011236|0.954494382022|0.939325842697
4|2.16136266|0.666666666|0.600561797753|0.908988764045|0.908988764045
5|1.54065323|0.468599033|0.302247191011|0.555056179775|0.701685393258
6|1.07134294|0.178743961|0.170786516854|0.560112359551|0.469101123596
7|0.73664978|0.0628019323|0.0848314606742|0.494382022472|0.170786516854
8|0.30103888|0.05797|0.0696629213|0.362921348315|0.0241573033708
9|0.35336993|0.0|0.0848314606|0.135393258427|0.0

## Notes
"New ILI cases" plot is a simplified model of newly-detected symptomatic cases. 
It is defined here as 0.25 of infected in day 1, 0.5 of infected in day 2 and 0.25 of cases in day 3.  The sum is multiplied by the symptom rate (0.88 is typical for fever and influenza).  User can select a different symptom rate to model another pathogen or symptom type (e.g. respiratory or nasal).

"Peak duration" is the interval between the first and last day with at least two new ILI cases. 
Two days is chosen to make it robust to noise at the start and end of the outbreak.

The Calibration Dataset should have exactly two columns: first column "date"" in the format "YYYY-MM-DD" and a second column reporting cases named either "New_ILI" or "New_absent"

