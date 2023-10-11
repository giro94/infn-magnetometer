#usage:
#./fix_files_symbols.sh <folder> [extension=csv]

if [[ $# -lt 1 ]] ; then
	echo "Usage:"
	echo "./fix_files_symbols.sh <folder> [extension=csv]"
fi

folder=${1}
ext=${2:-csv}

if [ -d "${folder}" ] ; then
	sed -i -e 's/\xE2\x88\x9E/10000/g' -e 's/NaN/10000/g' -e 's/inf/10000/g' ${folder}/*.${ext}
else
	echo "The argument you provided is not a folder! (${folder})"
fi

