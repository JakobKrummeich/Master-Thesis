parentDirectory=$1
fileToRemove=$2

for dir in ${parentDirectory}*/ ; do

	echo "${dir}/${fileToRemove}"
	rm "${dir}/${fileToRemove}"

done

