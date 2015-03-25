#!/bin/bash

pmeta=$( ls -t pickup.ckpt*.meta |head -1 )
psuff=${pmeta:7:5}
# contains whitespace, should be ok though
iter=$( grep timeStepNumber $pmeta | sed -e 's/.*\[ \(.*\) \];/\1/' )

# strip from data
sed -i.bak -n -e '/pickupSuff/I!p' data
sed -i.bak -e "s/nIter0=\(.*\)/nIter0=$iter,\n pickupSuff='$psuff',/I" data

