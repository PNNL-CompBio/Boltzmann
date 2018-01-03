#!/bin/csh
/home/dbaxter/cmake-3.9.3/build/bin/cmake -DCMAKE_INSTALL_PREFIX=/home/dbaxter/sundials/sundials-2.7.0/install -DEXAMPLES_INSTALL_PATH=/home/dbaxter/sundials/sundials-2.7.0/install/examples -DFCMIX_ENABLE=ON -DCMAKE_BUILD_PREFIX=/home/dbaxter/sundials/sundials-2.7.0/build -DCMAKE_ROOT=/home/dbaxter/cmake-3.9.3/ /home/dbaxter/sundials/sundials-2.7.0/src/cvode

