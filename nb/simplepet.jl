### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 10c0c25c-51c8-4df9-a173-5219b0066b0d
begin
	using Geant4
	using Geant4.SystemOfUnits
	using Geant4.PhysicalConstants
	using Geant4.SystemOfUnits: mm, m, micrometer, cm3, ns
	import CairoMakie
	using GeometryBasics, Rotations, IGLWrap_jll
end

# ╔═╡ d96ff12a-912b-42f1-aa5a-e2e3e094013b
begin
	using PlutoUI
end

# ╔═╡ d3b7d9f0-c1e3-4b7c-8fdb-60c1f87c5288



# ╔═╡ d611396b-42c5-45f1-a7d2-a9ca35e201b4
md"""
## Total Body PET
### Define the geometry of the setup
"""

# ╔═╡ 7d57da90-c068-49be-86a0-44b9bdb65650
md""" Select crystal : $(@bind xname PlutoUI.Select(["LYSO", "BGO", "CSI"]))"""

# ╔═╡ a6ec6b55-c3a4-422f-87e4-85bdc4f0653f
md""" Select size of crystal along x (mm) : $(@bind xt PlutoUI.NumberField(1.0:0.1:100.0; default=40.0))"""


# ╔═╡ 16cc8b24-b929-4db5-9600-088ecd47200d
md""" Select size of crystal along y (mm) : $(@bind yt PlutoUI.NumberField(1.0:0.1:100.0; default=6.0))"""

# ╔═╡ 7f2f5fe2-07b7-42df-92b5-a399093a13e4
md""" Select size of crystal along z (mm) : $(@bind zt PlutoUI.NumberField(1.0:0.1:100.0; default=6.0))"""

# ╔═╡ e6a841fb-cdbc-47fc-acef-ad93b9ec0970
begin
	sizeOfCrystal_x = xt * mm
	sizeOfCrystal_tangential = yt * mm 
	sizeOfCrystal_axial = zt * mm 
	crho = Dict("LYSO"=>11.0mm, "BGO"=>11.0mm, "CsI"=>11.0mm)
	cdecayTime=Dict("LYSO"=>40*ns, "BGO"=>300*ns, "CsI"=>1000*ns) 
	clightYieldPhKeV = Dict("LYSO"=>25, "BGO"=>8, "CsI"=>100)
	xz = 2*crho[xname]
	nothing
end

# ╔═╡ cf75476c-44f6-4a81-b14d-28d43ebc60f4
md"""
#### Crystal
- crystal = $xname
- radiation length = $(crho[xname]) mm
- decay time = $(cdecayTime[xname]) ns
- Light Yield = $(clightYieldPhKeV[xname]) photons/keV
- Crystal size: x = $xt mm; y = $yt mm; z = $zt mm;
"""

# ╔═╡ 47823b34-fa27-4964-9ffa-7731999718f0
md"""
#### Module (aka Block)
"""

# ╔═╡ d8f64915-fa87-43cf-8071-118f3e4024cc
md""" Select number of crystals per module along x : $(@bind numberOfCrystal_DOI PlutoUI.NumberField(1:10; default=1))"""

# ╔═╡ b8fd8ee9-458b-4a31-a7cb-3ade41937385
md""" Select number of crystals per module tangential (along y) : $(@bind numberOfCrystal_tangential PlutoUI.NumberField(1:20; default=8))"""

# ╔═╡ 2de84055-45d5-494f-aaf6-1d1c541143c2
md""" Select number of crystals per module axial (along z) : $(@bind numberOfCrystal_axial PlutoUI.NumberField(1:20; default=8))"""

# ╔═╡ 22474b14-f74b-44d5-b42e-2eee640d3ae8
md""" Select crystal gap (thickness of reflector in μm) : $(@bind crystalGap PlutoUI.NumberField(10.0:100.0; default=80.0))""" 

# ╔═╡ d443a4a2-652d-4a44-bd80-8e8516c122f6
begin
	
end

# ╔═╡ ec06700d-867e-4dab-88de-0873cba1e7d9
md""" Select al cover thickness (μm): $(@bind  actmu PlutoUI.NumberField(10.0:500.0; default=300.0))""" 

# ╔═╡ 61dfd680-0d2a-40ff-91a5-7b167e2e753d


# ╔═╡ 32b7f676-3db7-4208-b3f4-61ac2839bb3d
md"""
#### Wheel (aka RING)
"""

# ╔═╡ 03f8baef-c40c-4abd-b0ae-2c4dd3c0538f
md""" Select number of modules per Wheel  : $(@bind numberOfDetector_perRing PlutoUI.NumberField(1:100; default=51))"""

# ╔═╡ a1a32ef5-5d2a-4e21-ad57-0c0513d9c5ef
md""" Select number of Wheels in scanner  : $(@bind numberOfRings PlutoUI.NumberField(1:100; default=20))"""

# ╔═╡ fb17313d-6d76-4732-acca-058d0b8209f7
md""" Select ring gap in mm  : $(@bind rngap PlutoUI.NumberField(1.0:20.0; default=5.0))"""

# ╔═╡ 094936fc-019f-4ea2-be29-39ec888db81a
md"""
#### Scanner radius
"""

# ╔═╡ 735a817c-3aa8-492c-bace-03436f51220e
md""" Select scanner radius in mm  : $(@bind scnr PlutoUI.NumberField(100.0:1000.0; default=400.0))"""

# ╔═╡ 95a18519-34e1-4fcb-ab97-c12d37a7ab89
begin
	scannerRadius           = scnr * mm
	aluminumCoverThickness  = actmu * micrometer
	crystalGap_DOI          = crystalGap * micrometer
	crystalGap_tangential   = crystalGap * micrometer
	crystalGap_axial        = crystalGap * micrometer
	ringGap                 = rngap * mm
	nothing
end

# ╔═╡ 12270d91-2a57-4ac2-a277-7e3c806ddb3a
md"""
#### Pet detector
"""

# ╔═╡ 5683139e-0409-4c2c-9636-3df16d394410

#
md" Draw block? $(@bind drawblk PlutoUI.CheckBox(default=false))" 


# ╔═╡ 554b8b2e-5c5a-4ba1-8eb2-e97c104c3327
md" Draw scanner? $(@bind drawscanner PlutoUI.CheckBox(default=false))" 


# ╔═╡ b43e765f-c4ae-4e73-8fb8-d11d5a80a6db
mutable struct PetDetector <: G4JLDetector

	# World
	#const worldSizeXY::Float64
	#const worldSizeZ::Float64

	# crystal type 

	const crstName::String
	
    # size of the crystals in the DOI (x), tangential (y) and axial (z) direction.
	const sizeOfCrystal_x::Float64 
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

	# Computed variables
	# Crystal material
	crystalMaterial::CxxPtr{G4Material}

	# air box takes into account the space between crystals 
	# e.g., the reflectors. It could be filled with Teflon rather than aire
	# but for simple simulations does not make a difference. 
	
	sizeOfAirBox_DOI::Float64
	sizeOfAirBox_axial::Float64 
	sizeOfAirBox_tangential::Float64 

	# the size of the block is obtained adding the Al cover.
	sizeOfBlockDetector_DOI::Float64
	sizeOfBlockDetector_axial::Float64
	sizeOfBlockDetector_tangential::Float64

	# axial length 
	axialLength::Float64

	# World size.
	worldSizeXY::Float64
	worldSizeZ::Float64

	function PetDetector(; crstName = "LYSO",
							sizeOfCrystal_x = 40.0mm,
							sizeOfCrystal_tangential = 6.0mm,
							sizeOfCrystal_axial = 6.0mm,
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

		self = new(crstName,
				   sizeOfCrystal_x,
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
	
		self.sizeOfAirBox_DOI = (numberOfCrystal_DOI * sizeOfCrystal_x) + (
			                numberOfCrystal_DOI - 1) * crystalGap_DOI
		
		self.sizeOfAirBox_axial = (numberOfCrystal_axial * sizeOfCrystal_axial) + (
		                    numberOfCrystal_axial - 1) * crystalGap_axial
		
		self.sizeOfAirBox_tangential = (numberOfCrystal_tangential * 
		sizeOfCrystal_tangential) + (
		numberOfCrystal_tangential - 1) * crystalGap_tangential

		self.sizeOfBlockDetector_DOI = self.sizeOfAirBox_DOI + aluminumCoverThickness
		self.sizeOfBlockDetector_axial = self.sizeOfAirBox_axial + aluminumCoverThickness
		self.sizeOfBlockDetector_tangential = self.sizeOfAirBox_tangential + aluminumCoverThickness

		self.axialLength = numberOfRings * self.sizeOfBlockDetector_axial +(numberOfRings -1) * ringGap
		
		self.worldSizeXY = scannerRadius * 2 + self.sizeOfAirBox_DOI + 50mm
		self.worldSizeZ = self.axialLength + 50mm

		
	
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

    
    

# ╔═╡ c464eff2-733a-40ce-9b5e-03eae2c63913
pdet= PetDetector(crstName = xname,
	              sizeOfCrystal_x = sizeOfCrystal_x,
				  sizeOfCrystal_tangential = sizeOfCrystal_tangential,
	              sizeOfCrystal_axial = sizeOfCrystal_axial,
				  numberOfCrystal_DOI = numberOfCrystal_DOI,
				  numberOfCrystal_tangential = numberOfCrystal_tangential,      
				  numberOfCrystal_axial = numberOfCrystal_axial,
                  numberOfDetector_perRing = numberOfDetector_perRing,
				  numberOfRings = numberOfRings,
				  crystalGap_DOI = crystalGap_DOI,
			      crystalGap_tangential = crystalGap_tangential,
				  crystalGap_axial = crystalGap_axial,
				  aluminumCoverThickness = aluminumCoverThickness, 
                  scannerRadius = scannerRadius, 
				  ringGap=ringGap,
                  checkOverlaps=true)


# ╔═╡ fb0098b4-fe81-4a72-bcc6-bc6c7d0eedf9
md"""
#### Blocks
- size BlockDetector DOI (x) mm = $(pdet.sizeOfBlockDetector_DOI)
- size BlockDetector axial (z) mm = $(pdet.sizeOfBlockDetector_axial)
- size BlockDetector tangential (y) mm $(pdet.sizeOfBlockDetector_tangential)
"""

# ╔═╡ 47110cdf-c511-4f95-b1cb-81236e436659
function num_crystals(det::PetDetector) 
	crystal_module = (det.numberOfCrystal_DOI *det.numberOfCrystal_tangential *  
	                  det.numberOfCrystal_axial)
	
	crystals_ring = det.numberOfDetector_perRing * crystal_module
	
	det.numberOfRings * crystals_ring
end

# ╔═╡ 88f2904f-b1bf-406b-a0e4-69b533fe8281
md"""
#### Scanner
- number of crystals in scanner = $(num_crystals(pdet))
- Scanner diameter (mm) =$(pdet.scannerRadius * 2) (WorldXY = $(pdet.worldSizeXY)  )
- Scanner axial length (mm) =$(pdet.axialLength) (WorldZ = $(pdet.worldSizeZ)  )
"""

# ╔═╡ e61962d4-cb13-43d1-ac1d-550d0c69ff8e
num_crystals(pdet) 

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
	   sizeOfCrystal_x,
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

	
	sizeOfAirBox_DOI = (numberOfCrystal_DOI * sizeOfCrystal_x
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
	   sizeOfCrystal_x,
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
if drawscanner
	let
		pv = blockConstruction(pdet)
		draw(pv[], wireframe=true)
	end
end

# ╔═╡ 12e279bb-3a7d-45bd-9c46-526807c365ea
function crystalConstruction(det::PetDetector)::CxxPtr{G4VPhysicalVolume}
    (; worldSizeXY, worldSizeZ, crstName,
	   sizeOfCrystal_x,
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
	crystalSolid = G4Box("Crystal", sizeOfCrystal_x/2.0, 
                         sizeOfCrystal_tangential/2.0, sizeOfCrystal_axial/2.0)

	#Define the local volume of the crystal
	crystal_logicalV = G4LogicalVolume(move!(crystalSolid),
                                       crystalMaterial,"crystal_logicalV")

	# Place the crystals inside the detectors and give them a unique number with crystalIndex.

	crystalIndex = 0
	for i_DOI in 0:numberOfCrystal_DOI-1

		crystalPositionX=(i_DOI-numberOfCrystal_DOI/2.0 + 0.5)*(sizeOfCrystal_x + crystalGap_DOI)

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
if drawblk
	let
		pv = crystalConstruction(pdet)
		draw(pv[], wireframe=true)
	end
end

# ╔═╡ 09220519-0785-4819-a1c3-a5e1e92aeead
function petDetectorConstruction(det::PetDetector)::CxxPtr{G4VPhysicalVolume}
    (; worldSizeXY, worldSizeZ, crstName,
	   sizeOfCrystal_x,
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
		    " sizeOfCrystal_x = ", sizeOfCrystal_x,
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
	
	sizeOfAirBox_DOI = (numberOfCrystal_DOI * sizeOfCrystal_x
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
# ╠═d96ff12a-912b-42f1-aa5a-e2e3e094013b
# ╠═d3b7d9f0-c1e3-4b7c-8fdb-60c1f87c5288
# ╠═d611396b-42c5-45f1-a7d2-a9ca35e201b4
# ╠═7d57da90-c068-49be-86a0-44b9bdb65650
# ╠═a6ec6b55-c3a4-422f-87e4-85bdc4f0653f
# ╠═16cc8b24-b929-4db5-9600-088ecd47200d
# ╠═7f2f5fe2-07b7-42df-92b5-a399093a13e4
# ╠═e6a841fb-cdbc-47fc-acef-ad93b9ec0970
# ╠═cf75476c-44f6-4a81-b14d-28d43ebc60f4
# ╠═47823b34-fa27-4964-9ffa-7731999718f0
# ╠═d8f64915-fa87-43cf-8071-118f3e4024cc
# ╠═b8fd8ee9-458b-4a31-a7cb-3ade41937385
# ╠═2de84055-45d5-494f-aaf6-1d1c541143c2
# ╠═22474b14-f74b-44d5-b42e-2eee640d3ae8
# ╠═d443a4a2-652d-4a44-bd80-8e8516c122f6
# ╠═ec06700d-867e-4dab-88de-0873cba1e7d9
# ╠═61dfd680-0d2a-40ff-91a5-7b167e2e753d
# ╠═32b7f676-3db7-4208-b3f4-61ac2839bb3d
# ╠═03f8baef-c40c-4abd-b0ae-2c4dd3c0538f
# ╠═a1a32ef5-5d2a-4e21-ad57-0c0513d9c5ef
# ╠═fb17313d-6d76-4732-acca-058d0b8209f7
# ╠═094936fc-019f-4ea2-be29-39ec888db81a
# ╠═735a817c-3aa8-492c-bace-03436f51220e
# ╠═95a18519-34e1-4fcb-ab97-c12d37a7ab89
# ╠═12270d91-2a57-4ac2-a277-7e3c806ddb3a
# ╠═c464eff2-733a-40ce-9b5e-03eae2c63913
# ╠═fb0098b4-fe81-4a72-bcc6-bc6c7d0eedf9
# ╠═88f2904f-b1bf-406b-a0e4-69b533fe8281
# ╠═5683139e-0409-4c2c-9636-3df16d394410
# ╠═f2364278-3b4e-4cb7-bf1c-eac0eefd08ff
# ╠═554b8b2e-5c5a-4ba1-8eb2-e97c104c3327
# ╠═ba184fb2-1b74-4728-ab9e-eb213a178808
# ╠═47110cdf-c511-4f95-b1cb-81236e436659
# ╠═e61962d4-cb13-43d1-ac1d-550d0c69ff8e
# ╠═b43e765f-c4ae-4e73-8fb8-d11d5a80a6db
# ╠═fc732bc2-cc9e-4416-bd6a-71438b158402
# ╠═9e436e02-ccf7-4ee6-b104-1e0c34991131
# ╠═3f4f386c-11e1-44b1-89f8-abed4e719775
# ╠═f8387935-8726-4db9-8679-56dff32719e7
# ╠═12e279bb-3a7d-45bd-9c46-526807c365ea
# ╠═09220519-0785-4819-a1c3-a5e1e92aeead
# ╠═23a8ae05-b356-408c-be2d-0f73b07977bf
# ╠═27339a42-5579-11ee-25a2-cda93eb9e6f2
