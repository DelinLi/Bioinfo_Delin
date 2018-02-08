#!/bin/sh
###Download from http://blog.csdn.net/yangjun2/article/details/8080807
###Compare files' md5 from two files, i.e. 1)file from fileI exist in fileII? 2) If yes, do they share same md5 
###Modified by Delin LI delin.bio@gmail.com
###Date Mar 31, 2015

usage(){
  echo "usage:'MD5_Compare.sh fileI fileII' compare files' md5 code from fileI with the ones from fileII"
}


if [ $# -ne 2 ];
then
	usage
	exit 1
fi

cat $1 | while read myline
do
	e0=`echo $myline |awk '{print $1}'`
	e1=`echo $myline |awk '{print $2}'`
	count=`grep $e1 $2|wc -l`

	if [ $count -ne 1 ] ;
	then
		echo "$e1 more than one time in $2"
		exit 1
	fi
 
	te1=`grep $e1 $2|awk '{print $2}' `
	te0=`grep $e1 $2|awk '{print $1}' `
	if [ "$e1"x = "$te1"x ];
	then 
		if [ "$te0"x = "$e0"x ] ; 
		then
			echo "file:"$e1" equals!"
		else  
			echo "file:"$e1" not equals!"
		fi
   else
   	echo "file:"$e1" is not exist in $2"
   fi
done
