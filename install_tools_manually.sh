#!/usr/bin/sh

#### install CIRI
echo "Downloading CIRI..."
wget http://sourceforge.net/projects/ciri/files/CIRI-full/CIRI-full_v2.0.zip 
echo "Done!"
echo "Installing CIRI..."
unzip CIRI-full_v2.0.zip
rm CIRI-full_v2.0.zip
mv CIRI-full_v2.0 CIRI
sed -i "1i\\#\!/usr/bin/perl" ./CIRI/bin/CIRI_v2.0.6/CIRI2.pl
chmod a+x ./CIRI/bin/CIRI_v2.0.6/*
echo "export PATH=\"$(pwd)/CIRI/bin/CIRI_v2.0.6:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
echo "Finish CIRI installation!"

#### install Find_circ
echo "Downloading Find_circ..."
git clone http://github.com/marvin-jens/find_circ.git 
echo "Done!"
echo "Installing Find_circ..."
sed -i '1d' ./find_circ/unmapped2anchors.py
sed -i "1i\\#\!/usr/bin/python" ./find_circ/unmapped2anchors.py
sed -i '1d' ./find_circ/find_circ.py
sed -i "1i\\#\!/usr/bin/python" ./find_circ/find_circ.py
chmod a+x ./find_circ/*
echo "export PATH=\"$(pwd)/find_circ:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
echo "Finish Find_circ installation!"


