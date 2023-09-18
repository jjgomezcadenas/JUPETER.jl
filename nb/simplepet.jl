### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 10c0c25c-51c8-4df9-a173-5219b0066b0d
begin
	using Geant4
	using Geant4.SystemOfUnits
	using Geant4.PhysicalConstants
	using Geant4.SystemOfUnits: mm, m, cm3, ns
	import CairoMakie
	using GeometryBasics, Rotations, IGLWrap_jll
end

# ╔═╡ 2981407b-fe86-49f4-a56e-eb8750ea4803
(1966 - 3.34 * 15 * 36) / 35

# ╔═╡ 9b0557ec-2eaa-4179-9e6d-4d29eed78585
6 * 16 * 8

# ╔═╡ f9c9986c-c7c6-4b21-85a6-7f37d119a7ec
radius(sc, nc, nd)  = sc * nc * nd / 2π

# ╔═╡ 8abdeb8e-5124-465c-87cf-e6a7d05b4d1b
radius(6mm, 16, 51)

# ╔═╡ 6b19af15-4cdc-4bfa-8d86-02e742df8ba1
 (1600π - 2π*radius(6mm, 16, 51))/50

# ╔═╡ 5683139e-0409-4c2c-9636-3df16d394410

#
#draw(pv[], wireframe=true)

# ╔═╡ b43e765f-c4ae-4e73-8fb8-d11d5a80a6db
mutable struct PetDetector <: G4JLDetector

	# World
	const worldSizeXY::Float64
	const worldSizeZ::Float64

	# crystal type 

	const crstName::String
	
    # size of the crystals in the DOI (x), tangential (y) and axial (z) direction.
	const sizeOfCrystal_DOI::Float64 
	const sizeOfCrystal_tangential::Float64 
	const sizeOfCrystal_axial::Float64 
	
	# number of crystals
	const numberOfCrystal_DOI::Int             
	const numberOfCrystal_tangential::Int      
	const numberOfCrystal_axial::Int 
	
	#Number of PET detectors per ring (40)
	const numberOfDetector_perRing::Int

	#Number of rings in the PET system (4)
	const numberOfRings::Int
	
	#Inter crystal gap between crystals. (0, 0.08, 0.08)
	const crystalGap_DOI::Float64
	const crystalGap_tangential::Float64
	const crystalGap_axial::Float64

	#Define Aluminum thickness cover for the detector 0.3 * mm
	const aluminumCoverThickness::Float64

	#Radius of the PET system (330 mm) 
	#Note that the radius is defined between to opposing sensitive detectors 
	# (It does not include detector block cover).
	const scannerRadius::Float64
	
	#Gap between two adjacent rings. (10.4 mm)
	const ringGap::Float64	 

	# overlaps
    const checkOverlaps::Bool

	# Crystal material
	crystalMaterial::CxxPtr{G4Material}

	function PetDetector(; worldSizeXY = 1.0m,
						   worldSizeZ = 1.0m,
						   crstName = "LYSO",
	                       sizeOfCrystal_DOI = 40.0mm,
						   sizeOfCrystal_tangential = 3.0mm,
						   sizeOfCrystal_axial = 3.0mm,
                           numberOfCrystal_DOI=1,             
						   numberOfCrystal_tangential=16,      
						   numberOfCrystal_axial=16, 
	                       numberOfDetector_perRing = 40,
	                       numberOfRings=4,
						   crystalGap_DOI=0.08mm,
						   crystalGap_tangential=0.08mm,
                           crystalGap_axial=0.08mm,
	                       aluminumCoverThickness=0.3mm,
                           scannerRadius=330.0mm, 
	                       ringGap=10.4mm,
						   checkOverlaps = true)

		self = new(worldSizeXY, 
			       worldSizeZ,
			       crstName,
				   sizeOfCrystal_DOI,
	     		   sizeOfCrystal_tangential,
	     		   sizeOfCrystal_axial,
	     		   numberOfCrystal_DOI,             
		 		   numberOfCrystal_tangential,      
				   numberOfCrystal_axial, 
				   numberOfDetector_perRing,
				   numberOfRings,  
				   crystalGap_DOI,
				   crystalGap_tangential,
				   crystalGap_axial,
				   aluminumCoverThickness,
				   scannerRadius, 
				   ringGap,
				   checkOverlaps)

		nist     = G4NistManager!Instance()
		if crstName == "LYSO"
            lyso = G4Material("lyso", density= 7.4*g/cm3, ncomponents=4)
            Lu = FindOrBuildElement(nist,"Lu")
            Si = FindOrBuildElement(nist,"Si")
            O  = FindOrBuildElement(nist,"O")
            Y  = FindOrBuildElement(nist,"Y")
            AddElement(lyso, Lu, fractionmass=.71)
            AddElement(lyso, Si, fractionmass=.07)
            AddElement(lyso, O, fractionmass=.18)
            AddElement(lyso, Y, fractionmass=.04)
            self.crystalMaterial = CxxPtr(lyso)
        elseif crstName == "BGO"
            self.crystalMaterial  = FindOrBuildMaterial(nist, "G4_BGO")
        else
            self.crystalMaterial  = FindOrBuildMaterial(nist, "G4_CESIUM_IODIDE")
        end

		return self
    end
end

    
    

# ╔═╡ d3b7d9f0-c1e3-4b7c-8fdb-60c1f87c5288
pdet= PetDetector(worldSizeXY = 2.0m,
				  worldSizeZ = 2.0m, 
				  sizeOfCrystal_tangential = 6.0mm,
				  sizeOfCrystal_axial = 6.0mm,
				  numberOfCrystal_tangential=16,      
				  numberOfCrystal_axial=16,
                  numberOfDetector_perRing = 51,
                  scannerRadius=800.0mm, 
				  ringGap=5.0mm,
                  numberOfRings=16)


# ╔═╡ fc732bc2-cc9e-4416-bd6a-71438b158402
function worldConstruction(det::PetDetector)::CxxPtr{G4VPhysicalVolume}
	
	nist     = G4NistManager!Instance()
	air  = FindOrBuildMaterial(nist, "G4_AIR")
	
    solid_world      = G4Box("world",  det.worldSizeXY/2.0, 
		                               det.worldSizeXY/2.0, 
		                               det.worldSizeZ/2.0)   
    world_logicalV   = G4LogicalVolume(move!(solid_world), air, "world_logicalV")
    world_physicalV  = G4PVPlacement(nothing,     # no rotation
                              G4ThreeVector(),    # at (0,0,0)
                              world_logicalV,     # its logical volume
                              "world_physicalV",  # its name
                              nothing,            # its mother volume
                              false,              # no boolean operation
                              0,                  # copy number
                              det.checkOverlaps)

    SetVisAttributes(world_logicalV, G4VisAttributes!GetInvisible())
	return world_physicalV 
end

# ╔═╡ 9e436e02-ccf7-4ee6-b104-1e0c34991131
function sizeBlock(det::PetDetector)
	(; worldSizeXY, worldSizeZ, crstName,
	   sizeOfCrystal_DOI,
	   sizeOfCrystal_tangential,
	   sizeOfCrystal_axial,
       numberOfCrystal_DOI,             
	   numberOfCrystal_tangential,      
	   numberOfCrystal_axial, 
       numberOfDetector_perRing,
       numberOfRings,
       crystalGap_DOI,
       crystalGap_tangential,
       crystalGap_axial,
       aluminumCoverThickness,
       scannerRadius, 
       ringGap,
       checkOverlaps, crystalMaterial)  = det

	
	sizeOfAirBox_DOI = (numberOfCrystal_DOI * sizeOfCrystal_DOI
    ) + (numberOfCrystal_DOI - 1) * crystalGap_DOI
	sizeOfAirBox_axial = (numberOfCrystal_axial * sizeOfCrystal_axial
    ) + (numberOfCrystal_axial - 1) * crystalGap_axial
	sizeOfAirBox_tangential = (numberOfCrystal_tangential * sizeOfCrystal_tangential
    ) + (numberOfCrystal_tangential - 1) * crystalGap_tangential

	#Define the size of the detector block. 
	sizeOfBlockDetector_DOI = sizeOfAirBox_DOI + aluminumCoverThickness
	sizeOfBlockDetector_axial = sizeOfAirBox_axial + aluminumCoverThickness
	sizeOfBlockDetector_tangential = sizeOfAirBox_tangential + aluminumCoverThickness
	
	(sizeOfAirBox_DOI, sizeOfAirBox_axial, sizeOfAirBox_tangential,
	sizeOfBlockDetector_DOI, sizeOfBlockDetector_axial, sizeOfBlockDetector_tangential)
end

# ╔═╡ 3f4f386c-11e1-44b1-89f8-abed4e719775
"""
- Define air volume (box) to fill the detector block. 
Crystal elements (scintillators) are then placed inside.
"""
function airboxConstruction(det::PetDetector)::Tuple{CxxPtr{G4VPhysicalVolume},
	                                                 CxxPtr{G4LogicalVolume}}
	
	nist     = G4NistManager!Instance()
	aluminum  = FindOrBuildMaterial(nist, "G4_Al")
	air  = FindOrBuildMaterial(nist, "G4_AIR")
	
	(sizeOfAirBox_DOI, sizeOfAirBox_axial, sizeOfAirBox_tangential,
	sizeOfBlockDetector_DOI, sizeOfBlockDetector_axial, sizeOfBlockDetector_tangential) = sizeBlock(det)

	blockDetector = G4Box("blockDetector",sizeOfBlockDetector_DOI/2,
                          sizeOfBlockDetector_tangential/2,
                          sizeOfBlockDetector_axial/2)

	#Define the logical volume for the detector block
	blockDetector_logicalV = G4LogicalVolume(move!(blockDetector), aluminum,
                             "blockDetector_logicalV")

	#Define air (box) inside the detector block. Crystal elements will be placed in it.
	airBox = G4Box("airBox", sizeOfAirBox_DOI/2, 
		           sizeOfAirBox_tangential/2,sizeOfAirBox_axial/2)
	#Define the logical volume
	airBox_logicalV = G4LogicalVolume(move!(airBox),air,"airBox_logicalV")

	#Define its physical volume and place it inside the detector block
    airBox_physicalV  = G4PVPlacement(nothing,         # no rotation
                              G4ThreeVector(),         # at (0,0,0)
                              airBox_logicalV,         # its logical volume
                              "airBox_physicalV",      # its name
                              blockDetector_logicalV,  # its mother volume
                              false,                   # no boolean operation
                              0,                       # copy number
                              det.checkOverlaps) 
	return airBox_physicalV, CxxPtr(blockDetector_logicalV) 

end

# ╔═╡ f8387935-8726-4db9-8679-56dff32719e7
function blockConstruction(det::PetDetector)::CxxPtr{G4VPhysicalVolume}
    (; worldSizeXY, worldSizeZ, crstName,
	   sizeOfCrystal_DOI,
	   sizeOfCrystal_tangential,
	   sizeOfCrystal_axial,
       numberOfCrystal_DOI,             
	   numberOfCrystal_tangential,      
	   numberOfCrystal_axial, 
       numberOfDetector_perRing,
       numberOfRings,
       crystalGap_DOI,
       crystalGap_tangential,
       crystalGap_axial,
       aluminumCoverThickness,
       scannerRadius, 
       ringGap,
       checkOverlaps, crystalMaterial)  = det

	nist     = G4NistManager!Instance()
	aluminum  = FindOrBuildMaterial(nist, "G4_Al")
	air  = FindOrBuildMaterial(nist, "G4_AIR")
	
    #Define World
    world_physicalV = worldConstruction(det)
	world_logicalV  = GetLogicalVolume(world_physicalV)
	airBox_physicalV, blockDetector_logicalV = airboxConstruction(det)

	(sizeOfAirBox_DOI, sizeOfAirBox_axial, sizeOfAirBox_tangential,
	sizeOfBlockDetector_DOI, sizeOfBlockDetector_axial, sizeOfBlockDetector_tangential) = sizeBlock(det)

	#Each block detector is identified with its unique number, blockIndex.
	blockIndex = 0
	
	for rng in 0:numberOfRings -1
	
		#place the detectors in a ring along the axial direction. 
        #Note that the ring gap between two adjcent rings is measured 
        #from scintillator to scintillator.  It does not include the Aluminum 				# thickness cover.

		println(" ring number =", rng)
		
		detectorPositionZ = (rng - numberOfRings/2.0 
        + 0.5) * (sizeOfBlockDetector_axial + ringGap -
		          aluminumCoverThickness)

		for i in 0:numberOfDetector_perRing -1
		#for i in 0:0

			#println(" detector number =", i)
		
			#The azimuthal angle to arrange the detectors in a ring
			thetaDetector = (i*2*π/numberOfDetector_perRing)

			#println(" thetaDetector =", thetaDetector)

			#The radius of the scanner is measured from opposing crystal (scintillator) faces. 
            #It does not include the Aluminum thickness cover.

			detectorPositionX = (scannerRadius + 
                                 sizeOfBlockDetector_DOI/2 - 
                                 aluminumCoverThickness/2) * cos(thetaDetector)

			#detectorPositionX = 0.0

			detectorPositionY = (scannerRadius + sizeOfBlockDetector_DOI/2 - 
                                 aluminumCoverThickness/2) * sin(thetaDetector)

			#println(" detectorPositionX =", detectorPositionX)
			#println(" detectorPositionY =", detectorPositionY)
			
			#Define the transformation matrix for correct placement of detetors

			rotm_PET = G4RotationMatrix()
			rotateZ(rotm_PET, thetaDetector)
			trans = G4Transform3D(rotm_PET,G4ThreeVector(detectorPositionX, detectorPositionY,detectorPositionZ))
			
			# Place block detectors

			println(" placing copy: blockIndex =", blockIndex)
			
            blockDetector_physicalV  = G4PVPlacement(trans ,      
                                                     blockDetector_logicalV[],
                                                     "blockDetector_physicalV",
                                                     world_logicalV[],
                                                     false, 
                                                     blockIndex,
                                                     checkOverlaps)

	
			blockIndex+=1
			

        end
	end
	return world_physicalV 
end

# ╔═╡ ba184fb2-1b74-4728-ab9e-eb213a178808
#pv = petDetectorConstruction(pdet)
pv = blockConstruction(pdet)

# ╔═╡ 12e279bb-3a7d-45bd-9c46-526807c365ea
function crystalConstruction(det::PetDetector)::CxxPtr{G4VPhysicalVolume}
    (; worldSizeXY, worldSizeZ, crstName,
	   sizeOfCrystal_DOI,
	   sizeOfCrystal_tangential,
	   sizeOfCrystal_axial,
       numberOfCrystal_DOI,             
	   numberOfCrystal_tangential,      
	   numberOfCrystal_axial, 
       numberOfDetector_perRing,
       numberOfRings,
       crystalGap_DOI,
       crystalGap_tangential,
       crystalGap_axial,
       aluminumCoverThickness,
       scannerRadius, 
       ringGap,
       checkOverlaps, crystalMaterial)  = det


    nist     = G4NistManager!Instance()
	aluminum  = FindOrBuildMaterial(nist, "G4_Al")
	air  = FindOrBuildMaterial(nist, "G4_AIR")
	
    #Define World
    world_physicalV = worldConstruction(det)
	world_logicalV  = GetLogicalVolume(world_physicalV)
	airBox_physicalV, blockDetector_logicalV = airboxConstruction(det)
	airBox_logicalV  = GetLogicalVolume(airBox_physicalV)

	(sizeOfAirBox_DOI, sizeOfAirBox_axial, sizeOfAirBox_tangential,
	sizeOfBlockDetector_DOI, sizeOfBlockDetector_axial, sizeOfBlockDetector_tangential) = sizeBlock(det)

	#Define the solid crystal
	crystalSolid = G4Box("Crystal", sizeOfCrystal_DOI/2.0, 
                         sizeOfCrystal_tangential/2.0, sizeOfCrystal_axial/2.0)

	#Define the local volume of the crystal
	crystal_logicalV = G4LogicalVolume(move!(crystalSolid),
                                       crystalMaterial,"crystal_logicalV")

	# Place the crystals inside the detectors and give them a unique number with crystalIndex.

	crystalIndex = 0
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
				crystal_physicalV = G4PVPlacement(nothing, 
                G4ThreeVector(crystalPositionX,crystalPositionY,crystalPositionZ), 
                crystal_logicalV, "crystal_physicalV", airBox_logicalV[],false,
                crystalIndex,false)

				crystalIndex+=1
            end
		end
	end
	
    return airBox_physicalV              # return a pointer to the G4PhysicalVolume
end


# ╔═╡ f2364278-3b4e-4cb7-bf1c-eac0eefd08ff
pv2 = crystalConstruction(pdet)

# ╔═╡ 5517fd83-f6cf-4f3a-adfd-6e61e44944c0
draw(pv2[], wireframe=true)

# ╔═╡ 09220519-0785-4819-a1c3-a5e1e92aeead
function petDetectorConstruction(det::PetDetector)::CxxPtr{G4VPhysicalVolume}
    (; worldSizeXY, worldSizeZ, crstName,
	   sizeOfCrystal_DOI,
	   sizeOfCrystal_tangential,
	   sizeOfCrystal_axial,
       numberOfCrystal_DOI,             
	   numberOfCrystal_tangential,      
	   numberOfCrystal_axial, 
       numberOfDetector_perRing,
       numberOfRings,
       crystalGap_DOI,
       crystalGap_tangential,
       crystalGap_axial,
       aluminumCoverThickness,
       scannerRadius, 
       ringGap,
       checkOverlaps, crystalMaterial)  = det

	println(" ++ pet Detector construction\n, worldSizeXY =", worldSizeXY,
	        " worldSizeZ = ", worldSizeZ,
		    " crstName   =",  crstName,
		    " sizeOfCrystal_DOI = ", sizeOfCrystal_DOI,
	        " sizeOfCrystal_tangential =", sizeOfCrystal_tangential,
	        " sizeOfCrystal_axial =", sizeOfCrystal_axial,
            " numberOfCrystal_DOI =", numberOfCrystal_DOI,          
	        " numberOfCrystal_tangential =", numberOfCrystal_tangential,
		    " numberOfCrystal_axial =", numberOfCrystal_axial,
            " numberOfDetector_perRing =",numberOfDetector_perRing,
            " numberOfRings =", numberOfRings,
            " crystalGap_DOI =", crystalGap_DOI,
            " crystalGap_tangential =",crystalGap_tangential,
            " crystalGap_axial =", crystalGap_axial,
            " aluminumCoverThickness =", aluminumCoverThickness,
            " scannerRadius =", scannerRadius, 
            " ringGap =", ringGap,
            " checkOverlaps =", checkOverlaps)
	
	nist     = G4NistManager!Instance()
	aluminum  = FindOrBuildMaterial(nist, "G4_Al")
	air  = FindOrBuildMaterial(nist, "G4_AIR")
	
    #Define World
    

    solid_world      = G4Box("world",  worldSizeXY/2.0, 
		                               worldSizeXY/2.0, 
		                               worldSizeZ/2.0)   
    world_logicalV   = G4LogicalVolume(move!(solid_world), air, "world_logicalV")
    world_physicalV  = G4PVPlacement(nothing,     # no rotation
                              G4ThreeVector(),    # at (0,0,0)
                              world_logicalV,     # its logical volume
                              "world_physicalV",  # its name
                              nothing,            # its mother volume
                              false,              # no boolean operation
                              0,                  # copy number
                              checkOverlaps)

    SetVisAttributes(world_logicalV, G4VisAttributes!GetInvisible())

	#Each block detector is identified with its unique number, blockIndex.
	blockIndex = 0

	#Each crystal is identified with its unique number, crystalIndex
	crystalIndex = 0
	
	# Define air volume (box) to fill the detector block. 
    # Crystal elements (scintillators) are then placed inside.
	
	sizeOfAirBox_DOI = (numberOfCrystal_DOI * sizeOfCrystal_DOI
    ) + (numberOfCrystal_DOI - 1) * crystalGap_DOI
	sizeOfAirBox_axial = (numberOfCrystal_axial * sizeOfCrystal_axial
    ) + (numberOfCrystal_axial - 1) * crystalGap_axial
	sizeOfAirBox_tangential = (numberOfCrystal_tangential * sizeOfCrystal_tangential
    ) + (numberOfCrystal_tangential - 1) * crystalGap_tangential

	#Define the size of the detector block. 
	sizeOfBlockDetector_DOI = sizeOfAirBox_DOI + aluminumCoverThickness
	sizeOfBlockDetector_axial = sizeOfAirBox_axial + aluminumCoverThickness
	sizeOfBlockDetector_tangential = sizeOfAirBox_tangential + aluminumCoverThickness


	println(" sizeOfBlockDetector_DOI =", sizeOfBlockDetector_DOI, 
	        " sizeOfBlockDetector_axial =", sizeOfBlockDetector_axial,
		    " sizeOfBlockDetector_tangential =",sizeOfBlockDetector_tangential)

	#Define solid shape for the detector block
	blockDetector = G4Box("blockDetector",sizeOfBlockDetector_DOI/2,
                          sizeOfBlockDetector_tangential/2,
                          sizeOfBlockDetector_axial/2)

	#Define the logical volume for the detector block
	blockDetector_logicalV = G4LogicalVolume(move!(blockDetector), aluminum,
                             "blockDetector_logicalV")

	#Define air (box) inside the detector block. Crystal elements will be placed in it.
	airBox = G4Box("airBox", sizeOfAirBox_DOI/2, 
		           sizeOfAirBox_tangential/2,sizeOfAirBox_axial/2)
	#Define the logical volume
	airBox_logicalV = G4LogicalVolume(move!(airBox),air,"airBox_logicalV")

	#Define its physical volume and place it inside the detector block
    airBox_physicalV  = G4PVPlacement(nothing,         # no rotation
                              G4ThreeVector(),         # at (0,0,0)
                              airBox_logicalV,         # its logical volume
                              "airBox_physicalV",      # its name
                              blockDetector_logicalV,  # its mother volume
                              false,                   # no boolean operation
                              0,                       # copy number
                              checkOverlaps) 

	#for rng in 0:numberOfRings -1
	for rng in 0:0
	
		#place the detectors in a ring along the axial direction. 
        #Note that the ring gap between two adjcent rings is measured 
        #from scintillator to scintillator.  It does not include the Aluminum 				# thickness cover.

		println(" ring number =", rng)
		
		detectorPositionZ = (rng - numberOfRings/2.0 
        + 0.5) * (sizeOfBlockDetector_axial + ringGap - aluminumCoverThickness)

		detectorPositionZ = 0

		println(" detectorPositionZ =", detectorPositionZ)


		#for i in 0:numberOfDetector_perRing -1
		for i in 0:0

			println(" detector number =", i)
		
			#The azimuthal angle to arrange the detectors in a ring
			thetaDetector = (i*2*π/numberOfDetector_perRing)

			println(" thetaDetector =", thetaDetector)

			#The radius of the scanner is measured from opposing crystal (scintillator) faces. 
            #It does not include the Aluminum thickness cover.

			detectorPositionX = (scannerRadius + 
                                 sizeOfBlockDetector_DOI/2 - 
                                 aluminumCoverThickness/2) * cos(thetaDetector)

			detectorPositionX = 0.0

			detectorPositionY = (scannerRadius + sizeOfBlockDetector_DOI/2 - 
                                 aluminumCoverThickness/2) * sin(thetaDetector)

			println(" detectorPositionX =", detectorPositionX)
			println(" detectorPositionY =", detectorPositionY)
			
			#Define the transformation matrix for correct placement of detetors

			rotm_PET = G4RotationMatrix()
			rotateZ(rotm_PET, thetaDetector)
			trans = G4Transform3D(rotm_PET,G4ThreeVector(detectorPositionX, detectorPositionY,detectorPositionZ))
			
			# Place block detectors

			println(" placing copy: blockIndex =", blockIndex)
			
            blockDetector_physicalV  = G4PVPlacement(trans ,      
                                                     blockDetector_logicalV,
                                                     "blockDetector_physicalV",
                                                     world_logicalV,
                                                     false, 
                                                     blockIndex,
                                                     checkOverlaps)

	
			blockIndex+=1
			

        end
	end

	#Define the solid crystal
	crystalSolid = G4Box("Crystal", sizeOfCrystal_DOI/2.0, 
                         sizeOfCrystal_tangential/2.0, sizeOfCrystal_axial/2.0)

	#Define the local volume of the crystal
	crystal_logicalV = G4LogicalVolume(move!(crystalSolid),
                                       crystalMaterial,"crystal_logicalV")

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
				crystal_physicalV = G4PVPlacement(nothing, 
                G4ThreeVector(crystalPositionX,crystalPositionY,crystalPositionZ), 
                crystal_logicalV, "crystal_physicalV", airBox_logicalV,false,
                crystalIndex,false)

				crystalIndex+=1
            end
		end
	end
	
    return world_physicalV              # return a pointer to the G4PhysicalVolume
end



# ╔═╡ 23a8ae05-b356-408c-be2d-0f73b07977bf
Geant4.getConstructor(::PetDetector)::Function = petDetectorConstruction

# ╔═╡ 27339a42-5579-11ee-25a2-cda93eb9e6f2
"""
Returns the root dir on top of the notebook (nb) dir
"""
function findrdir()
	nbdir = split(@__DIR__,"/")
	reduce(string, [string(x,"/") for x in nbdir[1:end-1]])
end

# ╔═╡ daacda17-5c47-4e09-add1-657525eeefb9
using Pkg; Pkg.activate(findrdir())

# ╔═╡ Cell order:
# ╠═daacda17-5c47-4e09-add1-657525eeefb9
# ╠═10c0c25c-51c8-4df9-a173-5219b0066b0d
# ╠═d3b7d9f0-c1e3-4b7c-8fdb-60c1f87c5288
# ╠═2981407b-fe86-49f4-a56e-eb8750ea4803
# ╠═9b0557ec-2eaa-4179-9e6d-4d29eed78585
# ╠═f9c9986c-c7c6-4b21-85a6-7f37d119a7ec
# ╠═8abdeb8e-5124-465c-87cf-e6a7d05b4d1b
# ╠═6b19af15-4cdc-4bfa-8d86-02e742df8ba1
# ╠═ba184fb2-1b74-4728-ab9e-eb213a178808
# ╠═5683139e-0409-4c2c-9636-3df16d394410
# ╠═f2364278-3b4e-4cb7-bf1c-eac0eefd08ff
# ╠═5517fd83-f6cf-4f3a-adfd-6e61e44944c0
# ╠═b43e765f-c4ae-4e73-8fb8-d11d5a80a6db
# ╠═fc732bc2-cc9e-4416-bd6a-71438b158402
# ╠═9e436e02-ccf7-4ee6-b104-1e0c34991131
# ╠═3f4f386c-11e1-44b1-89f8-abed4e719775
# ╠═f8387935-8726-4db9-8679-56dff32719e7
# ╠═12e279bb-3a7d-45bd-9c46-526807c365ea
# ╠═09220519-0785-4819-a1c3-a5e1e92aeead
# ╠═23a8ae05-b356-408c-be2d-0f73b07977bf
# ╠═27339a42-5579-11ee-25a2-cda93eb9e6f2
