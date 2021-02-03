#!/bin/bash

python -m numpy.f2py efermi.f -m efermif77 -h efermif77.pyf --overwrite-signature
python -m numpy.f2py \
    --f77exec="gfortran" \
    --opt="-O3" \
    --f77flags="--fast-math" \
    -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION \
    -c efermif77.pyf efermi.f
