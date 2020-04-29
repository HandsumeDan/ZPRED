clear
Fldr=$(pwd)/
# Update OS
sudo apt-get update;
# Install the Build Essential package
sudo apt-get install build-essential;
# Install Qt Creator
sudo apt-get install qtcreator;
# Set Qt Version 5 to default
sudo apt install qt5-default;
# Acquire Additional Stuff
#sudo apt-get install qt5-doc qtbase5-examples qtbase5-doc-html;
# For Audio\
sudo apt-get install sox;
# INSTALL OTHER SOFTWARES
# Compile ZPRED
#. $Fldr"Software/ZPRED/C++/CompileZPRED.sh"
# Move into User Interface Folder
cd $Fldr"UI/"
# Make Qt Project File
proFile=$Fldr"UI/UI.pro"
xFile="UI"
qmake -project -o $proFile;
#
g++ -std=gnu++11 -c -m64 -pipe -O2 -Wall -W -D_REENTRANT -fPIC -DQT_NO_DEBUG -DQT_WIDGETS_LIB -DQT_GUI_LIB -DQT_CORE_LIB -I. -I. -isystem /usr/include/x86_64-linux-gnu/qt5 -isystem /usr/include/x86_64-linux-gnu/qt5/QtWidgets -isystem /usr/include/x86_64-linux-gnu/qt5/QtGui -isystem /usr/include/x86_64-linux-gnu/qt5/QtCore -I. -I/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64 -o main.o main.cpp -w -lboost_system -lboost_filesystem
#
/usr/lib/x86_64-linux-gnu/qt5/bin/moc -DQT_NO_DEBUG -DQT_WIDGETS_LIB -DQT_GUI_LIB -DQT_CORE_LIB -I/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64 -I/home/handsomedan/Desktop/forDoug/UI -I$Fldr"UI/" -I/usr/include/x86_64-linux-gnu/qt5 -I/usr/include/x86_64-linux-gnu/qt5/QtWidgets -I/usr/include/x86_64-linux-gnu/qt5/QtGui -I/usr/include/x86_64-linux-gnu/qt5/QtCore -I/usr/include/c++/5 -I/usr/include/x86_64-linux-gnu/c++/5 -I/usr/include/c++/5/backward -I/usr/lib/gcc/x86_64-linux-gnu/5/include -I/usr/local/include -I/usr/lib/gcc/x86_64-linux-gnu/5/include-fixed -I/usr/include/x86_64-linux-gnu -I/usr/include ui.h -o moc_ui.cpp
#
g++ -std=gnu++11 -c -m64 -pipe -O2 -Wall -W -D_REENTRANT -fPIC -DQT_NO_DEBUG -DQT_WIDGETS_LIB -DQT_GUI_LIB -DQT_CORE_LIB -I. -I. -isystem /usr/include/x86_64-linux-gnu/qt5 -isystem /usr/include/x86_64-linux-gnu/qt5/QtWidgets -isystem /usr/include/x86_64-linux-gnu/qt5/QtGui -isystem /usr/include/x86_64-linux-gnu/qt5/QtCore -I. -I/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64 -o moc_ui.o moc_ui.cpp -w -lboost_system -lboost_filesystem
#
g++ -std=gnu++11 -c -m64 -pipe -O2 -Wall -W -D_REENTRANT -fPIC -DQT_NO_DEBUG -DQT_WIDGETS_LIB -DQT_GUI_LIB -DQT_CORE_LIB -I. -I. -isystem /usr/include/x86_64-linux-gnu/qt5 -isystem /usr/include/x86_64-linux-gnu/qt5/QtWidgets -isystem /usr/include/x86_64-linux-gnu/qt5/QtGui -isystem /usr/include/x86_64-linux-gnu/qt5/QtCore -I. -I/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64 -o ui.o ui.cpp -w -lboost_system -lboost_filesystem
# Open (.pro) file and add "QT += core gui widgets"
echo "QT += core gui widgets" >> $proFile;
#echo "QT += multimedia multimediawidgets core gui widgets" >> $proFile;
# Platform Specific Project
qmake $proFile;
g++ -m64 -Wl,-O1 -o $xFile main.o ui.o moc_ui.o -L/usr/X11R6/lib64 -lQt5Widgets -lQt5Gui -lQt5Core -lGL -lpthread -w -lboost_system -lboost_filesystem
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
