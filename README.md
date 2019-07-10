# MATLAB_GET_EER_Binary_Search

Performs a Binary Search approach to estimate the Equal Error Rate (EER) for very large biometric databases. It was created by Lee Friedman, Ph.D. following a similar Python program yet to be uploaded to GitHub by Mr. Vladyslav Prokopenko.  It is not an exact copy and it has completely different stopping conditions.

It is posted and available in support of this submitted article:

Lee Friedman, Hal S. Stern, Vladyslav Prokopenko and Oleg V. Komogorthsev (Submitted, July 2019).
Relationship between Number of Subjects and Biometric Authentication Equal Error Rates
IEEE TRANSACTIONS ON BIOMETRICS, BEHAVIOR, AND IDENTITY SCIENCE, Submitted

Description: This is a MATLAB script and associated MATLAB functions. It is designed to compute the Equal Error Rate (EER) for very large data sets, with numbers of subjects up to 100,000 or beyond.  Standard ROC analyses run in memory requires that all the similarity scores be held in memory.  With 100,000 (10^5) subjects, there are 10,000,000,000 (10^10) similarity scores.  These will not fit in the memory of many standard desktop computers.   This program solves that problem.


Table of Contents: 

(1) MATLAB_GET_EER_Binary_Search.m - Main Program

(2) QuickDistance.m - Function to compute similarity scores

(3) fastAUC.m - Function to perform a tradition ROC analysis in memory for smaller datasets (N<= 10,000).
see:  

https://www.mathworks.com/matlabcentral/fileexchange/42860-fast-auc-calculator-and-roc-curve-plotter

(4) autoArrangeFigures.m - Function to tile figures on a display.  See: 

https://www.mathworks.com/matlabcentral/fileexchange/48480-automatically-arrange-figure-windows

(5) ..\InputData\ (Practice Data):

'Band_6_NFeat_010_NumberOfSubjects_00001000.csv';% 1,000 Subjects, 10 Features
'Band_6_NFeat_010_NumberOfSubjects_00010000.csv';%10,000 Subjects, 10 Features
'Band_6_NFeat_010_NumberOfSubjects_00050000.csv';%50,000 Subjects, 10 Features

(6) ..\OutputFiles\ contains some program output files.

Installation: Just download the zip and unpack it anywhere.  It uses relative paths. Open the main file in MATLAB. It was deveoped under MATLAB Version: 9.6.0.1099231 (R2019a) Update 1

Usage: The Main Program is self-documented. Just read the comments, modify it accordingly and run it.  It should run with the 10,000 subject database without modification.

License: MIT License
