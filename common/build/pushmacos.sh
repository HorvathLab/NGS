#!/bin/sh
rsync -avz --progress -e ssh ./*.macOS-*.tgz nedwards@edwardslab.bmcb.georgetown.edu:projects/Horvath-Lab/NGS/dist
