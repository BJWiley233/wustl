#! /bin/bash

# added myself as user in install_and_add_user.sh so can run without sudp
docker build -t my_image:version1.0 .
docker run -i -t my_image:version1.0
bedtools | head -n4
