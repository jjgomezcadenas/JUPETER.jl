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

# ╔═╡ da44a5a8-63d4-4f51-9ae5-4fd7ec1eeb9a
begin
	using Geant4
	using Geant4.SystemOfUnits
	using Geant4.SystemOfUnits: mm, cm3, ns
	using CairoMakie
	using GeometryBasics, Rotations, IGLWrap_jll
end

# ╔═╡ 318da72f-2518-42a2-a24d-ac4c0bae6b1b
using PlutoUI

# ╔═╡ b7ae5979-1af1-4a33-9b2f-7a5094876ba3
md""" Select crystal size : $(@bind xt NumberField(1.0:0.1:10.0; default=3.5))"""

# ╔═╡ ed1f83d7-2f71-47fb-9492-2402600fdb9f
md"""
#### Define crystals of LYSO, BGO and CsI (cryogenic)
"""

# ╔═╡ dd9b973e-11d1-4ae6-a37a-131e67aadb64


# ╔═╡ 2ff98bd5-55ff-4d89-a8f4-29c7631ba775
md"""
#### Define Particle gun
"""

# ╔═╡ 958fdf5a-baa6-4214-ab04-a7c0b3b43d69
md"""
#### Define application
"""

# ╔═╡ d6700feb-63f4-4d90-8a2c-25aaa9f57289
solidcrst = G4Box("CRYS", 0.5 , 0.5 , 0.5)

# ╔═╡ 1178e753-b33a-4ded-a71d-86aa107311e8
md"""
## Functions
"""

# ╔═╡ ff07c1de-5ac7-4076-8fd9-acd265d917d7
md"""
### Materials & Elements
"""

# ╔═╡ 934d809c-2233-4201-b01b-85feddc41f31
"""
	Defines a Material
"""
struct Material
	name::String # ej LYSO
	density::Float64 # ej: 7.4*g/cm3
	nelms::Int  # ej 4
	elnames::Vector{String} 
	elcomp::Vector{<:Real} 
	g4material::CxxPtr{G4Material}
	
	function Material(name::String; density::Float64, nelms::Int,
	                  elnames::Vector{String},
	                  elcomp::Vector{<:Real})

		# NB, rho must include g4 units
    	g4m = G4Material(name, density= density, ncomponents=nelms)
		nist     = G4NistManager!Instance()
		
		for (i,el) in enumerate(elnames)
			xel = FindOrBuildElement(nist,el)
			if typeof(elcomp[i]) == Float64
				AddElement(g4m, xel, fractionmass=elcomp[i])
			else
				AddElement(g4m, xel, natoms=elcomp[i])
			end
		end
            
		new(name, density, nelms, elnames, elcomp, CxxPtr(g4m))
	end

	function Material(g4name::String)
		nist = G4NistManager!Instance()
		g4m  = FindOrBuildMaterial(nist, g4name)
		new(g4name, 0.0, 0, [""], [0], g4m)
	end
end
	

# ╔═╡ c7e2a059-1322-427f-ad41-4ce00fbe9d7c
"""
	Defines a Scintillating Material
"""
struct ScintMaterial
	material::Material 
	decayTime::Float64 # use G4, eg, 40 ns LYSO
	lightYield::Float64 # e.g LYSO photons/keV	25
	x0::Float64 # e.g., LYSO 11 mm

	function ScintMaterial(material; decayTime::Float64, lightYield::Float64,
	                       x0::Float64)
		new(material, decayTime, lightYield, x0)
	end
end

# ╔═╡ 7d4ce84f-8694-45fe-b435-0819554695ce
md"""
### Detector construction
"""

# ╔═╡ 6019ebf8-a12b-446c-9e14-66ab340918a7
"""
Construct a box detector
"""
struct BoxDetector <: G4JLDetector
	name::String
	material::Material
    x::Float64 # use units
	y::Float64
    z::Float64
	lcrst::CxxPtr{G4LogicalVolume}
	
    function BoxDetector(name::String; material::Material, 
	                     x::Float64, y::Float64, z::Float64)
		
    	scrst  = G4Box(name, 0.5 * x, 0.5 * y, 0.5 * z)
    	logvol = G4LogicalVolume(move!(scrst),  material.g4material, name)
		lcrst  = CxxPtr(logvol)
    	new(name, material, x, y, z, lcrst)
    end
end



# ╔═╡ 220ae589-e928-43b6-860d-5624416d8b7b
begin
	g4lyso = Material("LYSO", density= 7.36*g/cm3, nelms=4,
	               elnames = ["Lu","Si","O","Y"],
                   elcomp = [0.714467891, 0.063714272, 0.181479787, 0.04033805])
	lyso =ScintMaterial(g4lyso, decayTime=40*ns, lightYield=25.0, x0=11*mm)
	
	crstlyso = BoxDetector("Crystal_LYSO"; material=g4lyso, x = xt*mm, y = xt*mm, z = lyso.x0)
end

# ╔═╡ e5e475ae-8dbc-4659-930b-f45d11d78dc9
begin
	g4bgo = Material("BGO", density= 7.13*g/cm3, nelms=3,
                 elnames = ["Bi", "Ge", "O"], 
	             elcomp  = [4, 3, 12])
	bgo =ScintMaterial(g4bgo, decayTime=300*ns, lightYield=8.0, x0=11*mm)
	
	crstbgo = BoxDetector("Crystal_BGO"; material=g4bgo, x = xt*mm, y = xt*mm, z = bgo.x0)
end

# ╔═╡ 4bb28472-64fa-4b1b-beff-d457b7fffcf5
begin
	g4csi = Material("CsI", density= 4.51*g/cm3, nelms=2,
                 elnames = ["Cs", "I"], 
	             elcomp  = [1, 1])
	csi =ScintMaterial(g4csi, decayTime=1000*ns, lightYield=100.0, x0=18.6*mm)
	crstcsi = BoxDetector("Crystal_CsI"; material=g4csi, x = xt*mm, y = xt*mm, z = csi.x0)
end

# ╔═╡ b3e0c7ad-0ddf-4786-99d5-825af569dd47
particlegun = G4JLGunGenerator(particle = "gamma", 
                               energy = 511keV, 
                               direction = G4ThreeVector(0,0,1), 
                               position = G4ThreeVector(0,0,-crstcsi.z/2))


# ╔═╡ 6fc012ec-9502-4b66-9acd-61f1a72a4d31
logiccrst = G4LogicalVolume(move!(solidcrst),  crstcsi.material.g4material, "CRYS")

# ╔═╡ 66045be2-05e4-4215-a179-a45cd9a0febc
crstcsi.lcrst

# ╔═╡ db502188-df2d-4ae2-be31-480601a873c7
"""
	Construct the detector. A simple crystal inside an air box
"""
function crstDetectorConstruction(det::BoxDetector)::CxxPtr{G4VPhysicalVolume}

	air = Material("G4_AIR")
	world = BoxDetector("WORLD_AIR"; material=air, 
	                     x = 2 * det.x, y = 2*det.y, z = 20.0*mm) 
	checkOverlaps=true
	
	#logicworld = world.lcrst
    #physWorld  = G4PVPlacement(nothing,          # no rotation
    #                          G4ThreeVector(),  # at (0,0,0)
    #                          logicworld,       # its logical volume
    #                          "World",          # its name
     #                         nothing,          # its mother volume
     #                         false,            # no boolean operation
     #                         0,                # copy number
     #                         checkOverlaps)                # overlaps checking 

	world_sizeXY   = 2 * det.x
    world_sizeZ    = 20.0*mm
    nist           = G4NistManager!Instance()
    worldMaterial  = FindOrBuildMaterial(nist, "G4_AIR")
    
    solidworld = G4Box("World",  0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ)   
    logicworld = G4LogicalVolume(move!(solidworld), worldMaterial, "World")
    physWorld  = G4PVPlacement(nothing,          # no rotation
                              G4ThreeVector(),  # at (0,0,0)
                              logicworld,       # its logical volume
                              "World",          # its name
                              nothing,          # its mother volume
                              false,            # no boolean operation
                              0,                # copy number
                              checkOverlaps)
    
	#logiccrst = det.lcrst
    #G4PVPlacement(nothing,           # no rotation
    #              G4ThreeVector(),   # at (0,0,0)
    #              logiccrst,          # its logical volume
    #              "CRYSTAL",            # its name
    #              logicworld,        # its mother  volume
    #              false,             # no boolean operation
    #              0,                 # copy number
    #              checkOverlaps)     # overlaps checking

	solidcrst = G4Box("CRYSTAL", 0.5 * det.x, 0.5 * det.y, 0.5 * det.z)
    logiccrst = G4LogicalVolume(move!(solidcrst),  det.material.g4material, "CRYSTAL")
    G4PVPlacement(nothing,           # no rotation
                  G4ThreeVector(),   # at (0,0,0)
                  logiccrst,          # its logical volume
                  "CRYSTAL",            # its name
                  logicworld,        # its mother  volume
                  false,             # no boolean operation
                  0,                 # copy number
                  checkOverlaps)
    SetVisAttributes(logicworld, G4VisAttributes!GetInvisible())
    SetVisAttributes(logiccrst, G4VisAttributes(G4Colour(1.0, 1.0, 0.0)))
    
    return physWorld              # return a pointer to the G4PhysicalVolume
end



# ╔═╡ 95232662-11ad-4d34-9f66-beaac6139504
"""
Pass the BoxDetector object to G4
"""
Geant4.getConstructor(::BoxDetector)::Function = crstDetectorConstruction

# ╔═╡ 24ef73f6-06fe-46be-a845-77a0ec02d93f
md"""
### Data structures
"""

# ╔═╡ 4e20c222-fdea-4955-b83d-f6a23b4310df
"""
	Define a Track 
"""
struct Track
    particle::String
    charge::Int
    energy::Float64
    trkid::Int
    edep::Vector{Float64}
    points::Vector{Point3{Float64}}
end



# ╔═╡ 77290633-897d-4786-9f36-62b8383dac8a
"""
	Define simulation data structure
"""
mutable struct CrstSimData <: G4JLSimulationData
    #---Run data
    fParticle::String
    fEkin::Float64
    fEdep::Float64
    #---trigger
    trigger::Bool
    #---tracks
    tracks::Vector{Track}
    CrstSimData() = new("", 0.0, 0.0, false, [])
end

# ╔═╡ 0f710cfb-5162-4199-a994-6c6e7b7c373a
md"""
### User actions
"""

# ╔═╡ 9206489b-9363-46ba-b1cc-14fe27114898
"""
Actions before starting the run
- stores the particle type and initial kinetic energy 
  of the generated primary particle

"""
function beginrun!(run::G4Run, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    gun = app.generator.data.gun
    data.fParticle = gun |> GetParticleDefinition |> GetParticleName |> String
    data.fEkin = gun |> GetParticleEnergy
    return
end


# ╔═╡ 8ffe91c9-8e4f-4513-b3c4-c557801d2145
"""
Actions before starting the event
- Clear the trigger and the list of tracks for the current event

"""
function beginevent!(::G4Event, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    data.trigger = false
    empty!(data.tracks)
    data.fEdep = 0.0
    return
end

# ╔═╡ 185e8d8d-e0d7-451a-b119-1e8be2849424
"""
Actions before tracking
- pushes a new `Track` with the particle name, 
  charge, intial energy and initial point

"""
function pretrackaction!(track::G4Track, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    tracks = data.tracks
    p = GetPosition(track)[]
    particle = track |> GetParticleDefinition
    name = particle |> GetParticleName |> String
    charge = particle |> GetPDGCharge |> Int
    energy = track |> GetKineticEnergy
    trkid = track |> GetTrackID
    push!(tracks, Track(name, charge, energy, trkid, [], 
         [Point3{Float64}(x(p),y(p),z(p))]))
    return
end

# ╔═╡ dda576f9-9d27-4ef6-b8ed-4ab8302c5758
"""
Actions at the end of tracking
- Used to set the trigger if the initial particle 
  exists the world without a sizeable interaction

"""
function posttrackaction!(track::G4Track, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    id = track |> GetTrackID
    energy = track |> GetKineticEnergy
    if id == 1 && energy > 0.99  * data.fEkin # primary particle did not losse any energy
        data.trigger = true
    end
    return
end

# ╔═╡ 7160dda2-bece-4f60-89b3-90be5942e77e
"""
Step actions
- Pushes points to the latest `Track` in the track list
- Accumulated edep

"""
function stepaction!(step::G4Step, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    tracks = data.tracks
    solid_name = step |> GetPreStepPoint |> GetTouchable |> GetVolume |> GetName |> String 
    #println("solid name = ", solid_name)

    #process_name = step |> GetPostStepPoint |> GetProcessDefinedStep |> GetProcessName

    #println("process name = ", process_name)

    if solid_name == "CRYSTAL"
        p = GetPosition(GetPostStepPoint(step))
        edep = GetTotalEnergyDeposit(step)
        data.fEdep += edep
        push!(tracks[end].points, Point3{Float64}(x(p),y(p),z(p)))
        push!(tracks[end].edep, edep)
    end
    return
end

# ╔═╡ da3595d6-ab96-4b84-9702-91297c8f66e3
app = G4JLApplication(; detector = crstlyso,             # detector with parameters
                        simdata = CrstSimData(),      # simulation data structure
                        generator = particlegun,      # primary particle generator
                        nthreads = 0,                 # # of threads (0 = no MT)
                        physics_type = QBBC,          # what physics list 
                        stepaction_method = stepaction!,           
                        begineventaction_method = beginevent!,    
                        pretrackaction_method = pretrackaction!,       
                        posttrackaction_method = posttrackaction!,     
                        beginrunaction_method = beginrun!              
                      );

# ╔═╡ 276b6734-220c-4f5f-9b2d-eb6db73bb711
configure(app)

# ╔═╡ 757d0bb2-0299-43ed-bd11-41996c190513
initialize(app)

# ╔═╡ 19e2b51d-4f4d-4443-b221-441c46ac2f37
md"""
### Utility functions
"""

# ╔═╡ f651b5b2-5091-11ee-0fd9-fda723590e78
function findrdir()
	nbdir = split(@__DIR__,"/")
	reduce(string, [string(x,"/") for x in nbdir[1:end-1]])
end

# ╔═╡ 87195243-2626-4509-b58d-bb86b6128f91
using Pkg; Pkg.activate(findrdir())

# ╔═╡ ede57068-01de-4aa4-a182-be4dfa1928f4
begin
	srcdir = string(findrdir(), "src")
	#include(joinpath(srcdir, "crstDet.jl"))
	#include(joinpath(srcdir, "crstSimData.jl"))
	#include(joinpath(srcdir, "crstUserActions.jl"))
end

# ╔═╡ Cell order:
# ╠═87195243-2626-4509-b58d-bb86b6128f91
# ╠═da44a5a8-63d4-4f51-9ae5-4fd7ec1eeb9a
# ╠═318da72f-2518-42a2-a24d-ac4c0bae6b1b
# ╠═ede57068-01de-4aa4-a182-be4dfa1928f4
# ╠═b7ae5979-1af1-4a33-9b2f-7a5094876ba3
# ╠═ed1f83d7-2f71-47fb-9492-2402600fdb9f
# ╠═220ae589-e928-43b6-860d-5624416d8b7b
# ╠═dd9b973e-11d1-4ae6-a37a-131e67aadb64
# ╠═e5e475ae-8dbc-4659-930b-f45d11d78dc9
# ╠═4bb28472-64fa-4b1b-beff-d457b7fffcf5
# ╠═2ff98bd5-55ff-4d89-a8f4-29c7631ba775
# ╠═b3e0c7ad-0ddf-4786-99d5-825af569dd47
# ╠═958fdf5a-baa6-4214-ab04-a7c0b3b43d69
# ╠═da3595d6-ab96-4b84-9702-91297c8f66e3
# ╠═276b6734-220c-4f5f-9b2d-eb6db73bb711
# ╠═757d0bb2-0299-43ed-bd11-41996c190513
# ╠═d6700feb-63f4-4d90-8a2c-25aaa9f57289
# ╠═6fc012ec-9502-4b66-9acd-61f1a72a4d31
# ╠═66045be2-05e4-4215-a179-a45cd9a0febc
# ╠═1178e753-b33a-4ded-a71d-86aa107311e8
# ╠═ff07c1de-5ac7-4076-8fd9-acd265d917d7
# ╠═934d809c-2233-4201-b01b-85feddc41f31
# ╠═c7e2a059-1322-427f-ad41-4ce00fbe9d7c
# ╠═7d4ce84f-8694-45fe-b435-0819554695ce
# ╠═6019ebf8-a12b-446c-9e14-66ab340918a7
# ╠═db502188-df2d-4ae2-be31-480601a873c7
# ╠═95232662-11ad-4d34-9f66-beaac6139504
# ╠═24ef73f6-06fe-46be-a845-77a0ec02d93f
# ╠═4e20c222-fdea-4955-b83d-f6a23b4310df
# ╠═77290633-897d-4786-9f36-62b8383dac8a
# ╠═0f710cfb-5162-4199-a994-6c6e7b7c373a
# ╠═9206489b-9363-46ba-b1cc-14fe27114898
# ╠═8ffe91c9-8e4f-4513-b3c4-c557801d2145
# ╠═185e8d8d-e0d7-451a-b119-1e8be2849424
# ╠═dda576f9-9d27-4ef6-b8ed-4ab8302c5758
# ╠═7160dda2-bece-4f60-89b3-90be5942e77e
# ╠═19e2b51d-4f4d-4443-b221-441c46ac2f37
# ╠═f651b5b2-5091-11ee-0fd9-fda723590e78
