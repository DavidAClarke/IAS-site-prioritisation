#!/bin/sh
zig4 -r species_wgt_CAZ/species_wgt_CAZ.dat species_wgt_CAZ/species_wgt_CAZ.spp species_wgt_CAZ/species_wgt_CAZ_out/species_wgt_CAZ.txt 0 0 1 0 --grid-output-formats=compressed-tif --image-output-formats=png --use-threads=4
