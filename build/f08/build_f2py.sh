#!/bin/bash

python -m numpy.f2py efermi.f90 -m efermif08 -h efermif08.pyf --overwrite-signature
python -m numpy.f2py \
    --f90exec="gfortran" \
    --opt="-O3" \
    --f90flags="--fast-math --std=f2008" \
    -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION \
    -c efermif08.pyf efermi.f90
