using Geant4
using Geant4.SystemOfUnits
using Geant4.SystemOfUnits: m, mm, cm3

const numberOfCrystal_DOI = 1
const numberOfCrystal_tangential = 16
const numberOfCrystal_axial = 16

# size of the crystals in the DOI (x), tangential (y) and axial (z) direction. 
const sizeOfCrystal_DOI = 20.5 * mm
const sizeOfCrystal_tangential = 2.8 * mm
const sizeOfCrystal_axial = 2.8 * mm

#Inter crystal gap between crystals. 
const crystalGap_DOI = 0.0
const crystalGap_tangential = 0.08 * mm
const crystalGap_axial = 0.08 * mm

#Define Aluminum thickness cover for the detector
const AluminumCoverThickness = 0.3 * mm

#Number of PET detectors per ring
const numberOfDetector_perRing = 40

#Number of rings in the PET system 
const numberOfRings = 4

#Radius of the PET system. 
#Note that the radius is defined between to opposing sensitive detectors 
# (It does not include detector block cover).

const scannerRadius = 330 * mm

#Gap between two adjacent rings. There are three gaps in four ring PET system
const ringGap = 10.4 * mm

                              
function petDetectorConstruction(fCheckOverlaps=true)

	#size of the world
	worldSizeX = 2 * m
	worldSizeY = 2 * m
	worldSizeZ = 4 * m


	#Define World
    air  = FindOrBuildMaterial(nist, "G4_AIR")

    solid_world = G4Box("world",  worldSizeX/2.0, worldSizeY/2.0, worldSizeZ/2.0)   
    world_logicalV = G4LogicalVolume(move!(solid_world), air, "world_logicalV")
    world_physicalV  = G4PVPlacement(nothing,     # no rotation
                              G4ThreeVector(),    # at (0,0,0)
                              world_logicalV,     # its logical volume
                              "world_physicalV",  # its name
                              nothing,            # its mother volume
                              false,              # no boolean operation
                              0,                  # copy number
                              fCheckOverlaps)

    SetVisAttributes(logicworld, G4VisAttributes!GetInvisible())
	#G4Box* solid_world = new G4Box("world", worldSizeX/2., worldSizeY/2., worldSizeZ/2.);

	# Define a logical volume for the world volume
	#world_logicalV = new G4LogicalVolume(solid_world, air, "world_logicalV", 0,0,0);


	# Define the physical world volume
	#world_physicalV = new G4PVPlacement(0,G4ThreeVector(),world_logicalV, "world_physicalV", 0, false, 0, fCheckOverlaps);
	#world_logicalV->SetVisAttributes (G4VisAttributes::Invisible);

	#NOTE!!!
	#The scanner specification (like size and number of crystals in each detctor) are given in the "doiPETGlobalParameters.hh" header file.

	#Each block detector is identified with its unique number, blockIndex.
	blockIndex = 0

	#Each crystal is identified with its unique number, crystalIndex
	crystalIndex = 0

	#This is to place the phantom. 

    #ConstructPhantom(world_logicalV);

	# Define air volume (box) to fill the detector block. 
    # Crystal elements (scintillators) is then placed.
	sizeOfAirBox_DOI = (numberOfCrystal_DOI * sizeOfCrystal_DOI
    ) + (numberOfCrystal_DOI - 1) * crystalGap_DOI
	sizeOfAirBox_axial = (numberOfCrystal_axial * sizeOfCrystal_axial
    ) + (numberOfCrystal_axial - 1) * crystalGap_axial
	sizeOfAirBox_tangential = (numberOfCrystal_tangential * sizeOfCrystal_tangential
    ) + (numberOfCrystal_tangential - 1) * crystalGap_tangential

	#This will pass the size of the detector block so that the position of the PMT will be calculated with respect to the axis of the detector block.
	#see doiPETAnalysis.cc
	#pAnalysis = doiPETAnalysis::GetInstance();
	#pAnalysis->GetSizeOfDetector(sizeOfAirBox_DOI,sizeOfAirBox_tangential, sizeOfAirBox_axial);

	println("size of crytal element: =", sizeOfCrystal_tangential,
            " ", sizeOfCrystal_axial, " ", sizeOfCrystal_DOI)

	println("Size of detector block (without Al cover): ",sizeOfAirBox_tangential,
            " ", sizeOfAirBox_axial, " ", sizeOfAirBox_DOI)


	#Define the size of the detector block. 
	sizeOfBlockDetector_DOI = sizeOfAirBox_DOI + AluminumCoverThickness
	sizeOfBlockDetector_axial = sizeOfAirBox_axial + AluminumCoverThickness
	sizeOfBlockDetector_tangential = sizeOfAirBox_tangential + AluminumCoverThickness

	#Define solid shape for the detector block
    
	blockDetector = G4Box("blockDetector",sizeOfBlockDetector_DOI/2,
                          sizeOfBlockDetector_tangential/2,
                          sizeOfBlockDetector_axial/2)

	#Define the logical volume for the detector block
    Aluminum  = FindOrBuildMaterial(nist, "G4_Al")
    world_logicalV = G4LogicalVolume(move!(solid_world), air, "world_logicalV")
	blockDetector_logicalV = G4LogicalVolume(move!(blockDetector),Aluminum,
                             "blockDetector_logicalV")


	#Define air (box) inside the detector block. Crystal elements will be placed in it.
	airBox = G4Box("airBox", sizeOfAirBox_DOI/2, sizeOfAirBox_tangential/2,sizeOfAirBox_axial/2)

	#Define the logical volume
	airBox_logicalV = G4LogicalVolume(move!(airBox),air,"airBox_logicalV")

	#Define its physical volume and place it inside the detector block
	#airBox_physicalV = new G4PVPlacement (0,G4ThreeVector(0,0,0),airBox_logicalV,"airBox_physicalV", blockDetector_logicalV,false,0,fCheckOverlaps);

    airBox_physicalV  = G4PVPlacement(nothing,         # no rotation
                              G4ThreeVector(),         # at (0,0,0)
                              airBox_logicalV,         # its logical volume
                              "airBox_physicalV",      # its name
                              blockDetector_logicalV,  # its mother volume
                              false,                   # no boolean operation
                              0,                       # copy number
                              fCheckOverlaps)

	#Arrange the PET ring and place the PET detectors in the ring(s) 

	for rng in 0:numberOfRings -1 
	
		#place the detectors in a ring along the axial direction. 
        #Note that the ring gap between two adjcent rings is measured 
        #from scintillator to scintillator.  It does not include the Aluminum thickness cover.

		detectorPositionZ = (rng - numberOfRings/2.0 
        + 0.5) * (sizeOfBlockDetector_axial + ringGap - AluminumCoverThickness)

		for i in 0:numberOfDetector_perRing -1
		
			#The azimuthal angle to arrange the detectors in a ring
			thetaDetector = (i*2*π/numberOfDetector_perRing)

			#The radius of the scanner is measured from opposing crystal (scintillator) faces. 
            #It does not include the Aluminum thickness cover.

			detectorPositionX = (scannerRadius + 
                                 sizeOfBlockDetector_DOI/2 - 
                                 AluminumCoverThickness/2) * cos(thetaDetector)

			detectorPositionY = (scannerRadius + sizeOfBlockDetector_DOI/2 - 
                                 AluminumCoverThickness/2) * sin(thetaDetector)

			#Define the rotation matrix for correct placement of detetors
			rotm_PET = G4RotationMatrix()
			rotateZ(rotm_PET, thetaDetector)

			uz_PET = G4ThreeVector(detectorPositionX,detectorPositionY,detectorPositionZ)
			#transform = G4Transform3D(rotm_PET,uz_PET)

			#Define the physical volume of the detectors.
            blockDetector_physicalV  = G4PVPlacement(rotm_PET, 
                                                     uz_PET ,      
                                                     blockDetector_logicalV,
                                                     "blockDetector_physicalV",
                                                     world_logicalV,
                                                     false, 
                                                     blockIndex,
                                                     fCheckOverlaps)

	
			blockIndex++
			prinln("Ring =", rng, 
                   " ", detectorPositionX - (
                  (sizeOfBlockDetector_DOI - AluminumCoverThickness)/2)*cos(thetaDetector),
                  " ", detectorPositionY- (
                  (sizeOfBlockDetector_DOI- AluminumCoverThickness)/2)*sin(thetaDetector),
                  " ", detectorPositionZ)

        end
	end

	#Define the solid crystal
	CrystalSolid = G4Box("Crystal", sizeOfCrystal_DOI/2.0, 
                         sizeOfCrystal_tangential/2.0, sizeOfCrystal_axial/2.0)

	#Define the local volume of the crystal
	crystal_logicalV = G4LogicalVolume(move!(CrystalSolid),
                                       crystalMaterial,"Crystal_logicalV")

	# Place the crystals inside the detectors and give them a unique number with crystalIndex.

	for i_DOI in 0:numberOfCrystal_DOI-1

		crystalPositionX=(i_DOI-numberOfCrystal_DOI/2.0 + 0.5)*(sizeOfCrystal_DOI + crystalGap_DOI)

		for i_axial in 0:numberOfCrystal_axial-1 
			crystalPositionZ = (i_axial-numberOfCrystal_axial/2.0 
            + 0.5)*(sizeOfCrystal_axial + crystalGap_axial)

			for i_tan in 0:numberOfCrystal_tangential - 1
				crystalPositionY=(i_tan- numberOfCrystal_tangential/2.0 
                + 0.5)*(sizeOfCrystal_tangential + crystalGap_tangential)

				println("crystalIndex = ", crystalIndex,
                        " crystalPositionX = ", crystalPositionX,
                        " crystalPositionY = ", crystalPositionY,
                        " crystalPositionZ =", crystalPositionZ)

				#place the crystal inside the block detector. 
				crystal_physicalV = G4PVPlacement (0, 
                G4ThreeVector (crystalPositionX,crystalPositionY,crystalPositionZ), 
                crystal_logicalV, "Crystal_physicalV", airBox_logicalV,false,
                crystalIndex,false)

				crystalIndex++;
            end
		end
	end
