#!/bin/bash

mkdir RMSD-OpenFF
cp ../computermsd.out .
cp ../computematch.out .

./computermsd.out
./computematch.out > OpenFF-MM.report
# for specific subset ---> awk '$2 ~ /GNT-/ {print $0}' OpenFF-MM.report |grep MD | tee fingerprint.report
awk '/ MD  /' OpenFF-MM.report > fingerprint.report
awk '/NO Num QM orphans/' OpenFF-MM.report > orphans
awk '{s+=$(NF)} END{print s}' orphans
# MD  0 NC  6 NM  6 NO  0 EO  0 PC  0 LQ  0 LF  0 BF  5 MO  6 BE  2 QS  0 FS  0 SC 66 NL 60
awk '{SUM+=$4}END{print "MD "SUM}' fingerprint.report
awk '{SUM+=$6}END{print "NC "SUM}' fingerprint.report
awk '{SUM+=$8}END{print "NM "SUM}' fingerprint.report
awk '{SUM+=$10}END{print "NO "SUM}' fingerprint.report
awk '{SUM+=$12}END{print "EO "SUM}' fingerprint.report
awk '{SUM+=$14}END{print "PC "SUM}' fingerprint.report
awk '{SUM+=$16}END{print "LQ "SUM}' fingerprint.report
awk '{SUM+=$18}END{print "LF "SUM}' fingerprint.report
awk '{SUM+=$20}END{print "BF "SUM}' fingerprint.report
awk '{SUM+=$22}END{print "MO "SUM}' fingerprint.report
awk '{SUM+=$24}END{print "BE "SUM}' fingerprint.report
awk '{SUM+=$26}END{print "QS "SUM}' fingerprint.report
awk '{SUM+=$28}END{print "FS "SUM}' fingerprint.report
awk '{SUM+=$30}END{print "SC "SUM}' fingerprint.report
awk '{SUM+=$32}END{print "NL "SUM}' fingerprint.report

