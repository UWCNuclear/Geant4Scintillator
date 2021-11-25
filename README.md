***The Geant4 meeting recordings are on the group Google Drive.***

*Aliyu and Julien review the code in the first two videos :-)*

# Using York's Optical Photon code OpPhoton

In the MANDELA desktop at UniZulu, The directory lives in the Ubuntu subsystem at:

    /G4/G4WORK/OpPhoton

The detector construction is in:

    /OpPhoton/src/OpPhotonDetectorConstruction_EJ200.cc

To compile the simulation the first time:

    cd OpPhoton
    mkdir build
    cd build
    cmake ..
    make -j
      
For a visualisation of the setup:

    ./OpPhoton
      
To simulate an vent, you can click on the green arrow at the top or type "/run/beamOn 1" in the Session box at the bottom.
 
To open the macro and change the number of events to run:

    gedit annihilGamma.mac
        
To run the code:

    ./OpPhoton -m annihilGamma.mac
        
To analyse the results:

    cd ../output
    root -l OpPhoton.C
        
To adapt the analysis code for the parameters that you need (stick length, number of files to analyse, their names,  and the names of the output images):  

    gedit OpPhoton.C
    root -l OpPhoton.C
  

***More details:*** OpPhoton/README

# Geant4 tips

-- Geant4 assumes dimensions are in cm if the units are not specified. To be safe, always specify units using "star"m, "star"cm, "star"mm, etc.
