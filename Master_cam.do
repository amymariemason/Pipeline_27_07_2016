********************************
* PIPELINE ANALYSIS MASTER DOC
********************************
* MASTER.DO
* VERSION: CAMBRIDGE

* project: staph pipeline
* author: Amy Mason
* start date: July 2016

*********************************
* This project looked at 3 pipeline analysis of staph aureus results and compared their ability to identify genotypes for antibiotic resistance; 
* all three methods are compared site by site; then all the methods are compared to the lab results
* The three methods considered are: Mykrobe, Typewriter and Genefinder. 

* raw data files in \\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs
* stata databases in \\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Datasets
* Table and list output in \\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Excel_output
* Image output in \\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Graphs_Outputs

sysdir set PERSONAL "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\ADO"
sysdir set STBPLUS  "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\ADO\stbplus"
sysdir set OLDPLACE "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\ADO\oldplace"
set logtype text

*******************
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\"

set more off
* data input: this takes the raw input files and turns them into stata datasets, with some minimal renaming of variables but no cleaning
run inputs.do


* clean data : these two programs clean up the data in the stata datasets; first the phenotype data and then the three methods data
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\"
run clean_pheno.do
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\"
run clean_predict.do

* prediction of phenotypes; uses the cleaned site predictions to make antibiotic predictions
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\"
run create_predict_anti.do
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\"
run create_predict_viru.do

* DOWN TO HERE

* analyse data
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\"
*compares A/P site predictions; makes agreement graph 1
run Analysis_Pipeline.do 

*compares phenotype predictions; makes agreement graph 2
*compares phenotype predictions to goldstandard; makes agreement graph 3; makes sensitivity and specifitity tables/graphs

cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\"
run Analysis_Pipeline2.do 

* compares virulence predictions
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\"
run Analysis_Pipeline3.do 

exit

cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\"
*discrepancies check; Sarah asked me to check the discrepancies aren't due to threshold problems
*TO DO: ROC transfer
run site.do
* makes ROC graphs and finds optimal cutoffs
run ROC.do

