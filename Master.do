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
*run inputs.do


* clean data
cd E:\users\amy.mason\Pipeline_27_07_2016\
run clean_pheno.do
cd E:\users\amy.mason\Pipeline_27_07_2016\
run clean_predict.do
cd E:\users\amy.mason\Pipeline_27_07_2016\
run create_predict_anti.do
cd E:\users\amy.mason\Pipeline_27_07_2016\
run create_predict_viru.do


* create prediction panals


* analyse data
