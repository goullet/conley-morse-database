#!/bin/bash

modeldir=$1
resultsdir=$2

subdirs=`find "$resultsdir" -mindepth 1 -maxdepth 1 -type d -name '*' | sed 's:^\./::'`

for i in $subdirs
do
  # name of the run
  name=`echo $i | rev | cut -d'/' -f1 | rev`

  #cp -r $i/Zoo $modeldir/$name/

  cp scriptMerge.sh scriptMerge_$name.sh
  sed -i.bak "s|CASE|$resultsdir/$name|g" scriptMerge_$name.sh
  sed -i.bak "s|MODEL|$modeldir/$name|g" scriptMerge_$name.sh 

  qsub scriptMerge_$name.sh

done

rm *.bak

cp $resultsdir/database.db $modeldir/
cp $resultsdir/indexSQL.html $modeldir/
cp $resultsdir/checkradio-form.php $modeldir/
cp $resultsdir/index.html $modeldir/
cp $resultsdir/getListDatabases.php $modeldir/
