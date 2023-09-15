#!/bin/sh
zig4 -r species_wgt_CAZ_KBA/species_wgt_CAZ_KBA.dat species_wgt_CAZ_KBA/species_wgt_CAZ_KBA.spp species_wgt_CAZ_KBA/species_wgt_CAZ_KBA_out/species_wgt_CAZ_KBA.txt 0 0 1 0 --grid-output-formats=compressed-tif --image-output-formats=png --use-threads=4
