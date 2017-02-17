compile sinogram processing code

g++ -o mMR_processing_adjust_add_gaps mMR_processing_adjust_add_gaps.c `root-config --cflags --libs`

./mMR_processing_adjust_add_gaps root_output_linesource.root GATE_sino




g++ -o mMR_processing_adjust mMR_processing_adjust.c `root-config --cflags --libs`

./mMR_processing_adjust root_output_linesource.root GATE_sino_2





g++ -o VersaPET_processing_gap VersaPET_processing_gap.c `root-config --cflags --libs`

./VersaPET_processing_gap root_output_linesource.root GATE_sino_VP