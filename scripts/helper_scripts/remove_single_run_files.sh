parentDirectory=$1
fileToRemove=$2

for dir in ${parentDirectory}*/ ; do

	rm "${dir}/${fileToRemove}"

done

