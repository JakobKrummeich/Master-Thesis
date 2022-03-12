startIndex=$1
endIndex=$2

for (( index=startIndex; index<=${endIndex}; index++ )); do

	echo "Canceling job ${index}"
	scancel ${index}

done

