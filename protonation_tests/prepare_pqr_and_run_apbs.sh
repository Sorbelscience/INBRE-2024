#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <pdb_file> <pH>"
    exit 1
fi

PDB_FILE=$1
PH=$2

# Convert the PDB file to PQR format using PDB2PQR
echo "Converting $PDB_FILE to ${PDB_FILE%.pdb}_ph${PH}.pqr at pH $PH..."
pdb2pqr --ff=AMBER --with-ph=$PH $PDB_FILE ${PDB_FILE%.pdb}_ph${PH}.pqr

if [ $? -ne 0 ]; then
    echo "pdb2pqr failed. Exiting."
    exit 1
fi

# Prepare input file for APBS
echo "Preparing input file for APBS..."
cat > ${PDB_FILE%.pdb}_ph${PH}.in <<EOF
read
    mol pqr ${PDB_FILE%.pdb}_ph${PH}.pqr
end
elec
    mg-auto
    dime 65 65 65
    cglen 80 80 80
    fglen 60 60 60
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl sdh
    ion charge -1 conc 0.15 radius 1.8
    ion charge 1 conc 0.15 radius 1.8
    pdie 2.0
    sdie 78.54
    srfm smol
    chgm spl2
    srad 1.4
    swin 0.3
    sdens 10.0
    temp 298.15
    calcenergy no
    calcforce no
    write pot dx ${PDB_FILE%.pdb}_ph${PH}.pot
end
quit
EOF

# Run APBS
echo "Running APBS..."
apbs ${PDB_FILE%.pdb}_ph${PH}.in

if [ $? -ne 0 ]; then
    echo "APBS failed. Exiting."
    exit 1
fi

echo "APBS calculation completed successfully."
