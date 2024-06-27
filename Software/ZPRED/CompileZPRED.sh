#!/bin/bash

clear
Compiler=/usr/bin/g++
Fldr=$(pwd)/; cFldr=$Fldr"C++/"; xFldr=$Fldr"bin/";

C_File=$cFldr"ZPRED.cpp"
X_File=$xFldr"ZPRED.exe"
Fldr2=$cFldr"alglib/src/"
other_C_Files=""
cd $Fldr2
for file in *; do
	len=${#file}	
	ext=${file:(( $len - 3 )):3}
	if [ $ext == "cpp" ]
		then
			other_C_Files=$other_C_Files$Fldr2$file" "
	fi
done
cd $Fldr
other_C_Files=$other_C_Files$cFldr"zpredUtils.cpp "$cFldr"mobility.cpp "$cFldr"solution.cpp "
# Remove Previous Executable
rm -f $X_File
# Compile Executable
$Compiler -std=gnu++11 -o $X_File -I$cFldr -I$Fldr2 $C_File $other_C_Files -lm -w -lrt -lboost_system -lboost_filesystem -lcrypto -lpthread
# Run Executable
#$X_File -i zpredInput5.txt

