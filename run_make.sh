# make --makefile ~/collabs/HeitmanLab/calcineurin_reg/calcineurin_reg.mk --directory /nfs/gems_sata/heitman/calcineurin_reg/rnaseq
# make all  -l 16 -j # for parallel make
make -l 8 -j 8 --makefile ~/collabs/HeitmanLab/calcineurin_reg/calcineurin_reg.mk --directory /nfs/gems_sata/heitman/calcineurin_reg/rnaseq
