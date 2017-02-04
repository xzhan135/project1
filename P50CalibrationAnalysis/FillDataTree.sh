#!/bin/bash

export ANALYSIS_DIR=/home/xzhan135/PROSPECT_CODE/project1/P50CalibrationAnalysis
export DATATREE=/home/xzhan135/data/P50Data

$ANALYSIS_DIR/FillDataTree $DATATREE/Data/BKGD-P50D-Data-S029.list
$ANALYSIS_DIR/FillDataTree $DATATREE/Data/Cf252-P50D-Data-S028-F016.list
$ANALYSIS_DIR/FillDataTree $DATATREE/Data/Na22-P50D-Data-S028-F014.list
$ANALYSIS_DIR/FillDataTree $DATATREE/Data/Na22-P50D-Data-S028-F015.list
$ANALYSIS_DIR/FillDataTree $DATATREE/Data/Cs137-P50D-Data-S028-F003.list
$ANALYSIS_DIR/FillDataTree $DATATREE/Data/Cs137-P50D-Data-S028-F010.list
