#!/bin/bash

# Checks that are missing in uncrustify
mycheck()
{
    file=$1

    # Detect multibytes
    if [[ ! -z $(file $file | grep UTF-8) ]]
    then
        echo "> Found multibyte character"
        grep --color='auto' -P -n "[^\x00-\x7F]" $file
    fi

    # Add whitespace before and after comments
    sed -i 's|\([^[:space:]]\)\*\/|\1 \*\/|g' $file
    sed -i 's|\/\*\*\([^[:space:]<]\)|\/\*\* \1|g' $file
    sed -i 's|\/\*\([^[:space:]*]\)|\/\* \1|g' $file

    # Remove extra whitespaces before pointers
    sed -i 's|\([a-zA-Z0-9]\)   *\(\*[\^s]\)|\1 \2|g' $file
}


dirs=$1

if [[ $dirs == "" ]]
then
    dirs="../include ../src ../test"
fi
echo "Checking $dirs"

files=$(find $dirs -type f | grep -E ".(h|c)$")
for file in $files
do
    mycheck $file
    uncrustify -c uncrustify.cfg --no-backup -f $file -o $file
    #clang-format -style=file -i -verbose
done


