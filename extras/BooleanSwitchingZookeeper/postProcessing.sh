#!/bin/bash 

# Destination directory where the different cases are
destdir=$1

# code path
codepath=$2

# Merge all the individual databases into one database for the model 
./mergeMultipleDatabases.sh $destdir $destdir/database.db 

# optimize the database
./optimizeSQLDatabase.sh $destdir/

# put together the morsegraphs files in Zoo/graphs
subdirs=`find "$destdir" -mindepth 1 -maxdepth 1 -type d -name "*" | sed 's:^\./::'`

for i in $subdirs
do
  rm -rf $i/Zoo

  mkdir $i/Zoo
  mkdir $i/Zoo/graphs 

  mv $i/MGCC* $i/Zoo/graphs/
  mv $i/HASSE* $i/Zoo/graphs/

  # Add the html/php files for the webpage (Shaun's stuff)
  cp -r template/*  $i/Zoo/
  ( cd $i/Zoo ; ./COMPUTE graphs )

  # Add the SQL html/php files for the parameter webpage
  MakeWebpage/QueryDatabase/templateForParameter/createWebpage.sh $codepath $i/Zoo/
  cp MakeWebpage/QueryDatabase/templateForParameter/checkradio-form.php $i/Zoo/
  cp MakeWebpage/QueryDatabase/templateForParameter/extractMorseSetMGCC.php $i/Zoo/
  cp MakeWebpage/QueryDatabase/templateForParameter/extractMorseSetMGCC.sh $i/Zoo/
  cp MakeWebpage/QueryDatabase/templateForParameter/extractParameterGraphMGCC.php $i/Zoo/
  cp MakeWebpage/QueryDatabase/templateForParameter/extractParameterGraphMGCC.sh $i/Zoo/
  cp MakeWebpage/QueryDatabase/templateForParameter/seeParameterGraphMGCC.php $i/Zoo/
  cp MakeWebpage/QueryDatabase/templateForParameter/seeParameterGraphMGCC.sh $i/Zoo/
  cp MakeWebpage/QueryDatabase/templateForParameter/indexParameterGraph.html $i/Zoo/
  cp MakeWebpage/QueryDatabase/templateForParameter/readGraphviz.php $i/Zoo/
done

# Add the SQL html/php files for the model webpage
  MakeWebpage/QueryDatabase/templateForModel/createWebpage.sh $codepath $destdir/
  cp MakeWebpage/QueryDatabase/templateForModel/checkradio-form.php $destdir/
  cp MakeWebpage/ListAllDatabases/index.html $destdir/
  cp MakeWebpage/ListAllDatabases/getListDatabases.php $destdir/
