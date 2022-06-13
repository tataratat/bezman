#!/bin/sh

# script to copy the headers to all the source files and header files
# original script from 
# https://unix.stackexchange.com/questions/20641/how-to-prepend-a-license-header-recursively-for-all-h-and-cpp-files-in-a-direc
# slightly modified
#!/bin/bash
for dir in ./src ./example ./tests
do
  for i in `find $dir \( -name '*.hpp' -o -name '*.cpp' -o -name '*.inc' \)`
  do
    echo $i
    if ! grep -q Copyright $i
    then
      echo "Add Copyright to " $i
      printf "/*\n" >> $i.new 
      cat LICENSE.md >> $i.new 
      printf "*/\n\n" >> $i.new
      cat $i >> $i.new
      mv $i.new $i
    fi
  done
done
