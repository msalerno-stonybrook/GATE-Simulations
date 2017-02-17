--------------------------------GATE------------------------------------------
To adjust the Gate macros for different scanners you only need to adjust
- Geometry.mac
- Digitzer.mac
- Phantom mac files
- Source mac files

simu_voxel.mac is the main mac file that you need to run -> it calls all of the other required mac files
- you can run this file with gate by going to the directory and typing Gate simu_voxel.mac
-------------------------------------------------------------------------------

--------------------------------Sorting Code----------------------------------
g++ -o VersaPET_processing VersaPET_processing.c `root-config --cflags --libs`
./VersaPET_processing root_output_linesource VP_linesource

./VersaPET_processing YourRootFile VP_linesource
-------------------------------------------------------------------------------
g++ -o VersaPET_processing_gap VersaPET_processing_gap.c `root-config --cflags --libs`
./VersaPET_processing_gap root_output_linesource VP_linesource_gap
-------------------------------------------------------------------------------
g++ -o mMR_processing mMR_processing.c `root-config --cflags --libs`
./mMR_processing root_mMR_linesource mMR_linesource
-------------------------------------------------------------------------------

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

***mMR***
- Without gaps
data type:
sx = 
sy = 
sz = 

- With gaps
data type:
sx = 
sy = 
sz = 
-------------------------------------------------------------------------------

-------------------------------Image Dimensions-------------------------------
***VersaPET***
data type:
sx = 127
sy = 127
sz = 31


***mMR***
data type:
sx = 344
sy = 344
sz = 127
-------------------------------------------------------------------------------

------------------------------------ROOT--------------------------------------
- to view root output file (tree), while in root use TBrowser with the following command

new TBroswer()
-------------------------------------------------------------------------------
