********************************
* PIPELINE ANALYSIS MASTER DOC
********************************
* MASTER.DO

* project: staph pipeline
* author: Amy Mason
* start date: July 2016
* raw data files in E:\users\amy.mason\Pipeline_27_07_2016\Inputs

sysdir set PERSONAL "E:\users\amy.mason\ADO"
sysdir set STBPLUS  "E:\users\amy.mason\ADO\stbplus"
sysdir set OLDPLACE "E:\users\amy.mason\ADO\oldplace"
set logtype text

*******************
cd E:\users\amy.mason\Pipeline_27_07_2016\

set more off
* data input
run inputs016.do


* clean data
cd E:\users\amy.mason\Pipeline_27_07_2016\
run clean_pheno.do
cd E:\users\amy.mason\Pipeline_27_07_2016\
run clean_predict16.do

* prediction of phenotypes
cd E:\users\amy.mason\Pipeline_27_07_2016\
run create_predict_anti.do
cd E:\users\amy.mason\Pipeline_27_07_2016\
run create_predict_viru.do

* analyse data
cd E:\users\amy.mason\Pipeline_27_07_2016\
run Analysis_Pipeline.do *compares A/P site predictions; makes agreement graph 1

cd E:\users\amy.mason\Pipeline_27_07_2016\
run Analysis_Pipeline2.do *compares phenotype predictions; makes agreement graph 2