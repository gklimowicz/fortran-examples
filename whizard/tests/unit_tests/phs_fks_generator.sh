#!/bin/sh
### Check WHIZARD phase space generator for FKS regions
echo "Running script $0"
exec ./run_whizard_ut.sh --check phs_fks_generator
