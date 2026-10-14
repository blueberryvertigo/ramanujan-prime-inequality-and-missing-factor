#!/usr/bin/env bash
# Regenerate every output used by the paper with the final scripts (after the referee fixes).
export PYTHONIOENCODING=utf-8
run() { python -B -u "$@"; }
getLC() { grep -oE "verified contiguously for L in \[[0-9.]+" "$1" | grep -oE "[0-9.]+$"; }
( run fks_fine_R.py 5.558691  2900 3160 PUB     KLN618 > final_C_PUB.txt 2>&1
  LC=$(getLC final_C_PUB.txt)
  ( run hybrid_D.py 5.558691 PUB "$LC" KLN618 > final_D_PUB.txt 2>&1 ) &
  ( run zfr_B.py "$LC" KLN618 > final_Bprime.txt 2>&1 ) &
  wait ) &
( run fks_fine_R.py 4.896     2400 3160 BTYBEL3 BEL618 > final_C_BTYBEL.txt 2>&1
  LC=$(getLC final_C_BTYBEL.txt)
  run hybrid_D.py 4.896 BTYBEL3 "$LC" BEL618 > final_D_BTYBEL.txt 2>&1 ) &
run fks_fine_R.py 5.558691  2850 3160 PUBBEL  BEL618 > final_C_PUBBEL.txt 2>&1 &
run fks_fine_R.py 5.558691  2900 3160 MTY77   KLN    > final_C_MTY77.txt 2>&1 &
run fks_fine_R.py 5.5666305 2900 3160 FKSR    KLN    > final_C_FKSR.txt 2>&1 &
run fks_fine_R.py 4.896     2300 3160 BTYKLN  KLN    > final_C_BTYKLN.txt 2>&1 &
wait
echo ALL_DONE > final_done.txt
