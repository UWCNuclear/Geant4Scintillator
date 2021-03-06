*The Geant4 meeting recordings are on the group Google Drive.*

***Aliyu and Julien review the code in the first two videos*** **:-)**

*Instructions to install Geant4 on Ubuntu:* https://github.com/UWCNuclear/UbuntuSetUp

To move files from Windows to your home directory on Ubuntu, paste in your Ubuntu terminal:

    mv /mnt/c/Users/PHYISICS/WHERE-IT-IS-ON-WINDOWS/filename ~

# Downloading York's OpPhoton code

 Paste in your terminal:
 
    cd ~/G4
    git clone https://github.com/UWCNuclear/Geant4Scintillator.git
    mv OpPhoton ..

# Using York's OpPhoton code

The directory lives in the Ubuntu subsystem at:

    ~/G4/OpPhoton

The detector construction is in:

    ~/G4/OpPhoton/src/OpPhotonDetectorConstruction_EJ200.cc
    
To edit the detector geometry:
    cd OpPhoton/src/
    gedit OpPhotonDetectorConstruction_EJ200.cc &

To create a build directory (the first time only):

    cd ~/G4/OpPhoton
    mkdir build
    cd build
      
To compile the simulation:

    cd OpPhoton/build/
    cmake ..
    make -j
      

For a visualisation of the setup:

    ./OpPhoton
      
To simulate an event, click on the green arrow at the top or type "/run/beamOn 1" in the Session box at the bottom.
 
To open the macro and change the number of events to run:

    gedit annihilGamma.mac &

***Before running the code,*** make sure that you edit the name of the previous output file or it will get overwritten when the simulation runs again:

    mv test.root ../output/test_SHAPE_LENGTH_REFLECTIVITY.root (e.g. test_Rod_10cm_975.root)

To run the code:

    ./OpPhoton -m annihilGamma.mac
        
To adapt the analysis code for the parameters that you need (stick length, number of files to analyse, their names, the names of the output images, and histogram legends):  

    mv test.root ../output/test_SHAPE_LENGTH_REFLECTIVITY.root
    cd ../output
    gedit OpPhoton.C &

Make sure to save different copies of OpPhoton.C to compare different parameters on the same plots.

To analyse the results:

    root -l OpPhoton.C
      
***More details:*** OpPhoton/README


# Geant4 tips!

- Geant4 assumes dimensions are in cm if the units are not specified. To be safe, always specify units using "star"m, "star"cm, "star"mm, etc.
- To edit a line, first make a copy of the line and comment out the original line.
