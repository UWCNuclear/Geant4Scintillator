***The Geant4 meeting recordings are on the group Google Drive.***

*Aliyu and Julien review the code in the first two videos* :-)

# Using York's OpPhoton code

In the MANDELA desktop at UniZulu, the directory lives in the Ubuntu subsystem at:

    /G4/G4WORK/OpPhoton

The detector construction is in:

    /OpPhoton/src/OpPhotonDetectorConstruction_EJ200.cc
    
To edit the detector geometry:
    cd OpPhoton/src/
    gedit OpPhotonDetectorConstruction_EJ200.cc &

To create a build directory (the first time only):

    cd OpPhoton
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

***Before running the code,*** make sure that you edit the name of the previous output file or it will get overwritten when the simulation runs again.

To run the code:

    ./OpPhoton -m annihilGamma.mac
        
To analyse the results:

    cd ../output
    root -l OpPhoton.C
        
To adapt the analysis code for the parameters that you need (stick length, number of files to analyse, their names,  and the names of the output images):  

    gedit OpPhoton.C &
    root -l OpPhoton.C
  

***More details:*** OpPhoton/README

# Geant4 tips!

- Geant4 assumes dimensions are in cm if the units are not specified. To be safe, always specify units using "star"m, "star"cm, "star"mm, etc.
- To edit a line, first make a copy of the line and comment out the original line.
- Save different copies of OpPhoton.C to compare different parameters on the same plots.
