#!/bin/bash

EIDIR=/ensembl/easy-import
CONFDIR=/import/conf
DATADIR=/import/data
SETUPDB=0
IMPORTVCF=0

DEFAULTINI="$CONFDIR/default.ini"
OVERINI="$CONFDIR/overwrite.ini"

while getopts "sid:" OPTION
do
  case $OPTION in
    s)  SETUPDB=1;;        # setup_database.pl
    i)  IMPORTVCF=1;;      # import_vcf.pl
    v)  VARIANTS=$OPTARG;; # core database name
  esac
done

# check VARIANTS has been specified
if [ -z ${VARIANTS+x} ]; then
  echo "ERROR: VARIANTS variable (-e VARIANTS=variants) has not been set"
  exit
fi

if ! [ -d $VARIANTS ]; then
  mkdir -p $DATADIR/$VARIANTS
fi

cd $DATADIR/$VARIANTS

if ! [ -d log ]; then
  mkdir -p log
fi

# check if $DEFAULTINI file exists
if ! [ -s $DEFAULTINI ]; then
  DEFAULTINI=
fi

# check if $OVERINI file exists
if ! [ -s $OVERINI ]; then
  OVERINI=
fi

if ! [ -z $INI ]; then
  OVERINI="$CONFDIR/$INI $OVERINI"
fi

# check main ini file exists
if ! [ -s $CONFDIR/$VARIANTS.ini ]; then
  echo "ERROR: no VARIANTS $VARIANTS.ini exists in conf dir"
  exit
fi

DBINI=$CONFDIR/$VARIANTS.ini

if ! [ $SETUPDB -eq 0 ]; then
  echo "setting up variation database"
  perl $EIDIR/variation/setup_database.pl $DEFAULTINI $DBINI $OVERINI &> >(tee log/setup_database.err)
fi

if ! [ $IMPORTVCF -eq 0 ]; then
  echo "importing variants from vcf file"
  perl $EIDIR/variation/import_vcf.pl $DEFAULTINI $DBINI $OVERINI &> >(tee log/import_vcf.err)
fi

cd ../
