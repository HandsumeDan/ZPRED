#!/bin/bash

clear
Compiler=/usr/bin/g++
Fldr=$(pwd)/; cFldr=$Fldr"C++/"; xFldr=$Fldr"bin/";

C_File=$cFldr"ZPRED.cpp"
X_File=$xFldr"ZPRED.exe"
soFile=$Fldr"libZPRED.so"
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
rm -f $soFile
# Compile Executable
#$Compiler -std=gnu++11 -o $X_File -I$cFldr -I$Fldr2 $C_File $other_C_Files -lm -w -lrt -lboost_system -lboost_filesystem -lcrypto -lpthread
# Run Executable
#$X_File -i zpredInput5.txt

# Make Shared Object Library
oFldr=$Fldr"objects/"
mkdir $oFldr
# Generate Object Files
cd $oFldr
$Compiler --std=gnu++11 -c -fPIC -I$cFldr -I$Fldr2 $other_C_Files -w

# Gather Objects Generated
O_Files=""
for file in *; do
	O_Files=$O_Files$oFldr$file" "
done

# Compile Shared Object Library
cd $Fldr
#g++ -fPIC -Xlinker --add-needed -Xlinker --no-as-needed -shared -o "libTestFunctions.so" $O_Files
$Compiler -shared -o $soFile $O_Files
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$Fldr
rm -f $O_Files
rmdir $oFldr
# Test Shared Object Library Function
$Compiler --std=gnu++11 -o $X_File -I$cFldr -I$Fldr2 $C_File -L$Fldr -lZPRED -lm -w -lrt -lboost_system -lboost_filesystem -lcrypto -lpthread
#$X_File -i zpredInput.txt
