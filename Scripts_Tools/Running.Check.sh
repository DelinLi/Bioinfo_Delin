#ssh -C2qTnN -D 8080 delin@129.186.85.56
Status=$(ps -A | grep -w "C2qTnN" | grep "ssh")
echo $Status
if ! [ -n "$Status" ]
then
	ssh -C2qTnN -D 8080 delin@129.186.85.56
fi
