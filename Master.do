********************************
* PIPELINE ANALYSIS MASTER DOC
********************************
* MASTER.DO

* project: staph pipeline
* author: Amy Mason
* start date: July 2016

*********************************
* This project looked at 3 pipeline analysis of staph aureus results and compared their ability to identify genotypes for antibiotic resistance; 
* all three methods are compared site by site; then all the methods are compared to the lab results
* The three methods considered are: Mykrobe, Typewriter and Genefinder. 

* raw data files in E:\users\amy.mason\Pipeline_27_07_2016\Inputs
* stata databases in E:\users\amy.mason\Pipeline_27_07_2016\Datasets
* Table and list output in E:\users\amy.mason\Pipeline_27_07_2016\Excel_output
* Image output in E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs

sysdir set PERSONAL "E:\users\amy.mason\ADO"
sysdir set STBPLUS  "E:\users\amy.mason\ADO\stbplus"
sysdir set OLDPLACE "E:\users\amy.mason\ADO\oldplace"
set logtype text

*******************
cd E:\users\amy.mason\Pipeline_27_07_2016\

set more off
* data input: this takes the raw input files and turns them into stata datasets, with some minimal renaming of variables but no cleaning
run inputs.do


* clean data : these two programs clean up the data in the stata datasets; first the phenotype data and then the three methods data
cd E:\users\amy.mason\Pipeline_27_07_2016\
run clean_pheno.do
cd E:\users\amy.mason\Pipeline_27_07_2016\
run clean_predict.do

* prediction of phenotypes; uses the cleaned site predictions to make antibiotic predictions
cd E:\users\amy.mason\Pipeline_27_07_2016\
run create_predict_anti.do
cd E:\users\amy.mason\Pipeline_27_07_2016\
run create_predict_viru.do

* DOWN TO HERE

* analyse data
cd E:\users\amy.mason\Pipeline_27_07_2016\
*compares A/P site predictions; makes agreement graph 1
run Analysis_Pipeline.do 

*compares phenotype predictions; makes agreement graph 2
*compares phenotype predictions to goldstandard; makes agreement graph 3; makes sensitivity and specifitity tables/graphs

cd E:\users\amy.mason\Pipeline_27_07_2016\
run Analysis_Pipeline2.do 

* compares virulence predictions
cd E:\users\amy.mason\Pipeline_27_07_2016\
run Analysis_Pipeline3.do 

exit

cd E:\users\amy.mason\Pipeline_27_07_2016\
*discrepancies check; Sarah asked me to check the discrepancies aren't due to threshold problems
*TO DO: ROC transfer
run site.do
* makes ROC graphs and finds optimal cutoffs
run ROC.do

