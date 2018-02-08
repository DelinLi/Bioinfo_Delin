return=`ps -u Data2Bio | grep "SCREEN -d -m ssh -C2qTnN -[D] 8080 delin@129.186.85.56"`
#echo $return
if [ ${#return} -eq 0 ] 
then
#       echo "No"
        screen -d -m ssh  -C2qTnN -D 8080  delin@129.186.85.56
fi
