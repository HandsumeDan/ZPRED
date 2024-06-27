#!/bin/bash

clear
Fldr=$(pwd)/
Compiler=/usr/bin/g++
# Update OS
#sudo apt-get update;
# Install the Build Essential package
#sudo apt-get install build-essential;
# Install Qt Creator
#sudo apt-get install qtcreator;
# Set Qt Version 5 to default
#sudo apt install qt5-default;
# Acquire Additional Stuff
#sudo apt-get install qt5-doc qtbase5-examples qtbase5-doc-html;
# For Audio\
#sudo apt-get install sox;
# Boost Libraries
#sudo apt-get install libboost-all-dev;
# Adaptive Poisson-Boltzmann Solver
#sudo apt-get install apbs;
# PyMOL Viewer
#sudo apt-get install pymol;
# INSTALL OTHER SOFTWARES
# Compile ZPRED
#cd $Fldr"Software/ZPRED/"
#./CompileZPRED.sh
# Move into User Interface Folder
cd $Fldr"UI/"
# Make Qt Project File
proFile=$Fldr"UI/UI.pro"
xFile="UI"
qmake -project -o $proFile;

mocMACRO="-DQT_NO_DEBUG -DQT_WIDGETS_LIB -DQT_GUI_LIB -DQT_CORE_LIB"
preINCLUDE="-std=gnu++11 -c -m64 -pipe -O2 -Wall -W -D_REENTRANT -fPIC "$mocMACRO
INCLUDE="-I. -isystem /usr/include/x86_64-linux-gnu/qt5 -isystem /usr/include/x86_64-linux-gnu/qt5/QtWidgets -isystem /usr/include/x86_64-linux-gnu/qt5/QtGui -isystem /usr/include/x86_64-linux-gnu/qt5/QtCore -I/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64"
moc=/usr/lib/x86_64-linux-gnu/qt5/bin/moc
mocINCLUDE="-I/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64 -I"$Fldr"UI/ -I/usr/include/x86_64-linux-gnu/qt5 -I/usr/include/x86_64-linux-gnu/qt5/QtWidgets -I/usr/include/x86_64-linux-gnu/qt5/QtGui -I/usr/include/x86_64-linux-gnu/qt5/QtCore -I/usr/include/c++/5 -I/usr/include/x86_64-linux-gnu/c++/5 -I/usr/include/c++/5/backward -I/usr/lib/gcc/x86_64-linux-gnu/5/include -I/usr/local/include -I/usr/lib/gcc/x86_64-linux-gnu/5/include-fixed -I/usr/include/x86_64-linux-gnu -I/usr/include ui.h"
#
$Compiler $preINCLUDE $INCLUDE -o main.o main.cpp -w -lboost_system -lboost_filesystem
$moc $mocMACRO $mocINCLUDE -o moc_ui.cpp
$Compiler $preINCLUDE $INCLUDE -o moc_ui.o moc_ui.cpp -w -lboost_system -lboost_filesystem
$Compiler $preINCLUDE $INCLUDE -o ui.o ui.cpp -w -lboost_system -lboost_filesystem
# Open (.pro) file and add "QT += core gui widgets"
echo "QT += core gui widgets" >> $proFile;
#echo "QT += multimedia multimediawidgets core gui widgets" >> $proFile;
# Platform Specific Project
qmake $proFile;
$Compiler -m64 -Wl,-O1 -o $xFile main.o ui.o moc_ui.o -L/usr/X11R6/lib64 -lQt5Widgets -lQt5Gui -lQt5Core -lGL -lpthread -w -lboost_system -lboost_filesystem
# Run Generated MakeFile
make;
cd ..
# Run UI
./"UI/"$xFile;


#rm -f *.o
#rm -f $Fldr"moc_button.cpp"
#rm -f $Fldr"moc_gui.cpp"
#rm -f $Fldr"moc_predefs.h"
#rm -f $Fldr"Makefile"
#rm -f $Fldr"ZPRED_UI"
#qmake -project
#qmake $Fldr"ZPRED_UI.pro"
#make
#./ZPRED_UI
