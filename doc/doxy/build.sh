#!/bin/bash

VERSION=$(grep -rnE "define CRY_(MAJOR|MINOR|PATCH)" ../../include/cry/version.h | awk '{ printf("%d", $3); if (NR!=3) printf("."); }')
echo $VERSION

sed -i "s/\\\$VERSION\\\$/v$VERSION/g" Doxyfile.in

doxygen Doxyfile.in
