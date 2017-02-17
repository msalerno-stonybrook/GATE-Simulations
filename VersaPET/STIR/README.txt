--------------------------------Sorting Code----------------------------------
g++ -o VersaPET_processing VersaPET_processing.c `root-config --cflags --libs`
./VersaPET_processing root_output_linesource VP_linesource

./VersaPET_processing YourRootFile VP_linesource
-------------------------------------------------------------------------------
g++ -o VersaPET_processing_gap VersaPET_processing_gap.c `root-config --cflags --libs`
./VersaPET_processing_gap root_output_linesource VP_linesource_gap
-----------------------------Sinogram Dimensions------------------------------
***VersaPET***
- Without gaps
data type:
sx = 191
sy = 96
sz = 256

- With gaps
data type:
sx = 215
sy = 108
sz = 256
------------------------------------------------------------------------------

-------------------------------Image Dimensions-------------------------------
***VersaPET***
data type:
sx = 127
sy = 127
sz = 31
-------------------------------------------------------------------------------
