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
import PlutoUI

# ╔═╡ b7ae5979-1af1-4a33-9b2f-7a5094876ba3
md""" Select crystal size : $(@bind xt PlutoUI.NumberField(1.0:0.1:10.0; default=3.5))"""

# ╔═╡ 709c0542-69fd-4709-a170-51e90dc019cc
md""" Select World xy : $(@bind wxt PlutoUI.NumberField(1.0:0.1:10.0; default=5.0))"""

# ╔═╡ 02508a2b-a483-46a7-aa5c-31f20d259ad6
md""" Select World z : $(@bind wz PlutoUI.NumberField(1.0:0.1:30.0; default=20.0))"""

# ╔═╡ 82ec0e6e-68db-4f51-8537-7176142455fe
md""" Select zgun : $(@bind zgun PlutoUI.NumberField(5.5:0.1:10.0; default=10.0))"""

# ╔═╡ ed1f83d7-2f71-47fb-9492-2402600fdb9f
md"""
#### Define crystals of LYSO, BGO and CsI (cryogenic)
"""

# ╔═╡ 04b01b62-ec6b-47c1-83ef-a179fec1b141
md"""
#### Define World (a box with air)
"""

# ╔═╡ 7aa3defb-3d9d-4032-9f8c-736aeb63aaa6
md"""
#### Define crystal detector for LYSO
"""

# ╔═╡ 2ff98bd5-55ff-4d89-a8f4-29c7631ba775
md"""
#### Define Particle gun
"""

# ╔═╡ b3e0c7ad-0ddf-4786-99d5-825af569dd47
particlegun = G4JLGunGenerator(particle = "gamma", 
                               energy = 511keV, 
                               direction = G4ThreeVector(0,0,1), 
                               position = G4ThreeVector(0,0,-zgun))


# ╔═╡ 958fdf5a-baa6-4214-ab04-a7c0b3b43d69
md"""
#### Define application
"""

# ╔═╡ 2584759f-42af-45d9-a5d1-1857798d94e2
md"Run app? $(@bind runapp PlutoUI.CheckBox(default=false))" 

# ╔═╡ 760d76c1-5cd7-46ab-8bd7-294b3a808ac3
md"""
#### Trigger and draw individual events
"""

# ╔═╡ 6da5fd9d-ddcf-467d-9511-6fd69e573890
md"trigger and draw events? $(@bind nxt PlutoUI.CheckBox(default=false))" 

# ╔═╡ 61a849b0-98d4-44e5-a930-c1156427a266
@bind go PlutoUI.Button("Trigger event")

# ╔═╡ 5fd011d6-fa2b-4ef4-ba61-06b3e90a41c3
md"Run analysis for a single crystal? $(@bind anx PlutoUI.CheckBox(default=false))" 

# ╔═╡ 548267cf-2d21-4508-b590-4b4ae54797b5
md""" number of events to run : $(@bind nrun PlutoUI.NumberField(100:100000; default=10000))"""

# ╔═╡ 3f674191-d5db-487e-aba5-941e531cb38b
@bind rgp PlutoUI.Button("Run!")

# ╔═╡ ba11b408-e4f0-4d69-9874-3b977914134b
md"""
##### Run with BGO 
"""

# ╔═╡ ca7335db-11ca-40cc-baf4-ee48d9fd136d
md"""
##### Run with CsI 
"""

# ╔═╡ 3619681e-6c74-471a-bd0e-5128a728d0d4
md"""
##### Plot and compare results
"""

# ╔═╡ 1178e753-b33a-4ded-a71d-86aa107311e8
md"""
## Functions
"""

# ╔═╡ 282e1e9c-c3db-4e38-bfd7-65e3a30051ca
md"""
### Geometry
"""

# ╔═╡ 934d809c-2233-4201-b01b-85feddc41f31
"""
	Defines a Geant 4 Material 
"""
struct g4Material
	name::String # ej LYSO
	density::Float64 # ej: 7.4*g/cm3
	nelms::Int  # ej 4
	elnames::Vector{String} 
	elcomp::Vector{<:Real} 
	g4mat::CxxPtr{G4Material}
	
	function g4Material(name::String; density::Float64, nelms::Int,
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

	function g4Material(g4name::String)
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
	material::g4Material 
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
#### Solid construction
"""

# ╔═╡ 6019ebf8-a12b-446c-9e14-66ab340918a7
"""
Construct a G4box 

"""
struct g4Box
	material::g4Material
    x::Float64 # use units
	y::Float64
    z::Float64
	lvol::CxxPtr{G4LogicalVolume}
	
    function g4Box(; material::g4Material, 
	               x::Float64, y::Float64, z::Float64)
		
    	scrst  = G4Box(material.name, 0.5 * x, 0.5 * y, 0.5 * z)
    	logvol = G4LogicalVolume(move!(scrst),  material.g4mat, material.name)
		lvol  = CxxPtr(logvol)
    	new(material, x, y, z, lvol)
    end
end



# ╔═╡ 220ae589-e928-43b6-860d-5624416d8b7b
begin
	g4lyso   = g4Material("LYSO", density= 7.36*g/cm3, nelms=4,
	                    elnames = ["Lu","Si","O","Y"],
                        elcomp = [0.714467891, 0.063714272, 0.181479787, 0.04033805])
	lyso     = ScintMaterial(g4lyso, decayTime=40*ns, lightYield=25.0, x0=11*mm)

	crstlyso = g4Box(material=g4lyso, x = xt*mm, y = xt*mm, z = lyso.x0)
end

# ╔═╡ e5e475ae-8dbc-4659-930b-f45d11d78dc9
begin
	g4bgo   = g4Material("BGO", density= 7.13*g/cm3, nelms=3,
                 elnames = ["Bi", "Ge", "O"], 
	             elcomp  = [4, 3, 12])
	bgo     = ScintMaterial(g4bgo, decayTime=300*ns, lightYield=8.0, x0=11*mm)
	
	crstbgo = g4Box(material=g4bgo, x = xt*mm, y = xt*mm, z = bgo.x0)
end

# ╔═╡ 4bb28472-64fa-4b1b-beff-d457b7fffcf5
begin
	g4csi   = g4Material("CsI", density= 4.51*g/cm3, nelms=2,
                 elnames = ["Cs", "I"], 
	             elcomp  = [1, 1])
	csi     = ScintMaterial(g4csi, decayTime=1000*ns, lightYield=100.0, x0=18.6*mm)
	crstcsi = g4Box(material=g4csi, x = xt*mm, y = xt*mm, z = csi.x0)
end

# ╔═╡ 0f3c5d45-0c38-40c0-9855-f257a0d846e5
begin
	air = g4Material("G4_AIR")
	world = g4Box(material=air, x = wxt * mm, y = wxt * mm, z = wz * mm) 
end

# ╔═╡ 689f4c8b-3530-422d-bc8e-65e407c45b86
"""
Construct a scintillating crystal
"""
struct crstDetector <: G4JLDetector
	world::g4Box
	crystal::g4Box
	checkOverlaps::Bool
    
    function crstDetector(; world::g4Box, crystal::g4Box, checkOverlaps=true)
    	new(world, crystal, checkOverlaps)
    end
end

# ╔═╡ f7fc9254-9d90-4882-b954-0116b5112ee2
lysodet = crstDetector(world=world, crystal=crstlyso, checkOverlaps=true)

# ╔═╡ 36d46822-3cc7-4714-a8c7-9da03ac60263
typeof(lysodet)

# ╔═╡ cb43ae23-eeca-4d6a-b282-51c54b579871
# NB, checkOverlaps set to false so that we can "overwrite" the previous crystal
bgodet = crstDetector(world=world, crystal=crstbgo, checkOverlaps=false)

# ╔═╡ 24ef73f6-06fe-46be-a845-77a0ec02d93f
md"""
### Data structures to extract information from event
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
if runapp
	app = G4JLApplication(; detector = lysodet,          # detector with parameters
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
	configure(app)
	initialize(app)
end

# ╔═╡ 3eaab173-16b8-4ceb-8dba-da6d3c22d694
reinitialize(app, bgodet)

# ╔═╡ 85a24190-8d20-4234-82a1-dd50a1891ab2
md"""
### Analysis functions
"""

# ╔═╡ e640afcd-92c3-4b3e-a95c-640677092f80
"""
Shoots the gun until we get an interaction
"""
function nexttrigger(app)
    data = app.simdata[1]
    beamOn(app,1)
    n = 1
    while data.trigger
        beamOn(app,1)
        n += 1
    end
    println("Got a trigger after $n generated particles")
end

# ╔═╡ 8948aeb7-5b7b-49c6-978a-4d6ec91aa1bc
"""
Run for nrun event and collect the deposited energy edep
"""
function run1(app; nrun=100)
    EDEP = zeros(nrun)
    nout = 0
    data = app.simdata[1]
    
    for i in 1:nrun
        beamOn(app,1)
        if data.trigger 
            nout+=1
        else
            data = app.simdata[1]
            EDEP[i] = data.fEdep
        end
    end
    EDEP
end

# ╔═╡ ee08a598-b0a1-4f64-a4e7-dc5d28598d20
if anx
	rgp
	edlyso = run1(app; nrun=nrun)
	f = Figure(resolution=(400,400))
	hist(f[1,1],edlyso, bins = 10, color = :red, strokewidth = 1, strokecolor = :black)
	hist(f[1,2],filter(x->x>0 && x < 0.510, edlyso), bins = 20, color = :red, strokewidth = 1, strokecolor = :black)
	f
end

# ╔═╡ 5bc0af4f-337e-4485-bbf2-84d38b9b09f0
edbgo = run1(app; nrun=nrun)

# ╔═╡ a8c3e6a2-c691-4b2e-a6fd-b872215fa053
begin
	csidet = crstDetector(world=world, crystal=crstcsi, checkOverlaps=false)
	reinitialize(app, csidet)
	edcsi = run1(app; nrun=nrun)
end

# ╔═╡ 5279e090-38ba-4c97-994d-c0611e19b5dc
md"""
- Crystal = LYSO
- Dimensions: xy = $(xt) mm, z = $(lyso.x0/mm) mm
- Fraction of events that interact in crystal = $(length(filter(x->x > 0, edlyso)) /length(edlyso))
- Fraction of events in photopeak = , $(length(filter(x->x > 0.509, edlyso)) /length(edlyso))

- Crystal = BGO
- Dimensions: xy = $(xt) mm, z = $(bgo.x0/mm) mm
- Fraction of events that interact in crystal = $(length(filter(x->x > 0, edbgo)) /length(edbgo))
- Fraction of events in photopeak = , $(length(filter(x->x > 0.509, edbgo)) /length(edbgo))

- Crystal = CsI
- Dimensions: xy = $(xt) mm, z = $(csi.x0/mm) mm
- Fraction of events that interact in crystal = $(length(filter(x->x > 0, edcsi)) /length(edcsi))
- Fraction of events in photopeak = , $(length(filter(x->x > 0.509, edcsi)) /length(edcsi))
"""

# ╔═╡ 665614bf-3f7c-45d5-b316-bbe9bd058603
"""
Draw detector

"""
function drawdetector!(s, app)
    world = GetWorldVolume()
    Geant4.draw!(s, world, wireframe=true)
    return s
end

# ╔═╡ ef947bdd-393f-45cb-95ab-aedb2cbd76e8
"""
Draw event
"""
function drawevent!(s, app, verbose=2)
    data = app.simdata[1]
    if verbose > 0
        println("Run info: Particle =", data.fParticle, 
                " Ekin = ", data.fEkin,
                " Edep = ", data.fEdep)
    end
    for t in data.tracks
        if t.particle == "gamma"
            if verbose > 1
                println(" gamma: energy = ", t.energy, " trkid =", t.trkid, 
					    " edep =", t.edep, 
                        " nof points = ", length(t.points), " points: = ", t.points)
                scatter!(s, t.points, markersize=7, color=:blue)
            end
        elseif t.particle == "e-"
            if verbose > 1
                println(" e-: energy = ", t.energy, " trkid =", t.trkid, 
					    " edep =", t.edep,
                        " nof points = ", length(t.points), " points: = ", t.points)
            end
            scatter!(s, t.points, markersize=7, color=:red)
        end
        
    end
end


# ╔═╡ 643549e6-0e11-48a7-911f-b351363cc6b1
if nxt
	go
	nexttrigger(app)

	fig = Figure(resolution=(500,500))
	s = LScene(fig[1,1])
	drawdetector!(s, app)
	drawevent!(s, app)
	display(fig)
end

# ╔═╡ 105e2342-7135-49a1-9458-0327be7b3ddb
function plot_edep(edep1, edep2, edep3)
	f = Figure(resolution=(500,500))
	hist(f[1,1],edep1, bins = 200, color = :red, strokewidth = 1, 
		 strokecolor = :black)
	hist(f[1,2],filter(x->x>0 && x < 0.510, edep1), bins = 20, color = :red, strokewidth = 1, strokecolor = :black)
	hist(f[2,1],edep2, bins = 200, color = :red, strokewidth = 1, 
		 strokecolor = :black)
	hist(f[2,2],filter(x->x>0 && x < 0.510, edep2), bins = 20, color = :red, strokewidth = 1, strokecolor = :black)
	hist(f[3,1],edep3, bins = 200, color = :red, strokewidth = 1, 
		 strokecolor = :black)
	hist(f[3,2],filter(x->x>0 && x < 0.510, edep3), bins = 20, color = :red, strokewidth = 1, strokecolor = :black)
	f
end

# ╔═╡ 92538392-a68c-4ef1-a302-a2dafbe07ba2
plot_edep(edlyso, edbgo, edcsi)

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

# ╔═╡ f49001d9-c1cb-41d9-8659-aca18e5b89c7
md"""
### Geant4-Julia Interface
"""

# ╔═╡ 92cfd1aa-ec72-4137-aaf5-b15c051925eb
G4JPVPlacement(r::Nothing, d::G4ThreeVector, 
	          l::CxxPtr{G4LogicalVolume}, s::String, 
              p::Nothing, b1::Bool, n::Int, b2::Bool=false) = G4PVPlacement(CxxPtr{G4RotationMatrix}(C_NULL), d, l, s, 
	          CxxPtr{G4LogicalVolume}(C_NULL), b1, n, b2)

# ╔═╡ 72da82b1-8fb9-45a8-9130-66d7f24b668f
G4JPVPlacement(r::Nothing, d::G4ThreeVector, 
	          l::CxxPtr{G4LogicalVolume}, s::String, 
              p::CxxPtr{G4LogicalVolume}, b1::Bool, n::Int, b2::Bool=false) = G4PVPlacement(CxxPtr{G4RotationMatrix}(C_NULL), d, l, s, 
	          p, b1, n, b2)

# ╔═╡ db502188-df2d-4ae2-be31-480601a873c7
"""
	Construct the detector. A simple crystal inside an air box
"""
function crstDetectorConstruction(det::crstDetector)::CxxPtr{G4VPhysicalVolume}
	
	logicworld = det.world.lvol
	println("Construct World: logicworld type ->", typeof(logicworld))
	
    physWorld  = G4JPVPlacement(nothing,          # no rotation
                                G4ThreeVector(),  # at (0,0,0)
                                logicworld,       # its logical volume
                                "World",          # its name
                                nothing,          # its mother volume
                                false,            # no boolean operation
                                0,                # copy number
                                det.checkOverlaps)                # overlaps checking 

	logiccrst = det.crystal.lvol  
    G4JPVPlacement(nothing, G4ThreeVector(),   logiccrst,    
                  "CRYSTAL",  logicworld,  false, 0, det.checkOverlaps)     

    SetVisAttributes(logicworld, G4VisAttributes!GetInvisible())
    SetVisAttributes(logiccrst, G4VisAttributes(G4Colour(1.0, 1.0, 0.0)))
    
    return physWorld              # return a pointer to the G4PhysicalVolume
end



# ╔═╡ 95232662-11ad-4d34-9f66-beaac6139504
"""
Pass the BoxDetector object to G4
"""
Geant4.getConstructor(::crstDetector)::Function = crstDetectorConstruction

# ╔═╡ Cell order:
# ╠═87195243-2626-4509-b58d-bb86b6128f91
# ╠═da44a5a8-63d4-4f51-9ae5-4fd7ec1eeb9a
# ╠═318da72f-2518-42a2-a24d-ac4c0bae6b1b
# ╠═ede57068-01de-4aa4-a182-be4dfa1928f4
# ╠═b7ae5979-1af1-4a33-9b2f-7a5094876ba3
# ╠═709c0542-69fd-4709-a170-51e90dc019cc
# ╠═02508a2b-a483-46a7-aa5c-31f20d259ad6
# ╠═82ec0e6e-68db-4f51-8537-7176142455fe
# ╠═ed1f83d7-2f71-47fb-9492-2402600fdb9f
# ╠═220ae589-e928-43b6-860d-5624416d8b7b
# ╠═e5e475ae-8dbc-4659-930b-f45d11d78dc9
# ╠═4bb28472-64fa-4b1b-beff-d457b7fffcf5
# ╠═04b01b62-ec6b-47c1-83ef-a179fec1b141
# ╠═0f3c5d45-0c38-40c0-9855-f257a0d846e5
# ╠═7aa3defb-3d9d-4032-9f8c-736aeb63aaa6
# ╠═f7fc9254-9d90-4882-b954-0116b5112ee2
# ╠═36d46822-3cc7-4714-a8c7-9da03ac60263
# ╠═2ff98bd5-55ff-4d89-a8f4-29c7631ba775
# ╠═b3e0c7ad-0ddf-4786-99d5-825af569dd47
# ╠═958fdf5a-baa6-4214-ab04-a7c0b3b43d69
# ╠═2584759f-42af-45d9-a5d1-1857798d94e2
# ╠═da3595d6-ab96-4b84-9702-91297c8f66e3
# ╠═760d76c1-5cd7-46ab-8bd7-294b3a808ac3
# ╠═6da5fd9d-ddcf-467d-9511-6fd69e573890
# ╠═61a849b0-98d4-44e5-a930-c1156427a266
# ╠═643549e6-0e11-48a7-911f-b351363cc6b1
# ╠═5fd011d6-fa2b-4ef4-ba61-06b3e90a41c3
# ╠═548267cf-2d21-4508-b590-4b4ae54797b5
# ╠═3f674191-d5db-487e-aba5-941e531cb38b
# ╠═ee08a598-b0a1-4f64-a4e7-dc5d28598d20
# ╠═ba11b408-e4f0-4d69-9874-3b977914134b
# ╠═cb43ae23-eeca-4d6a-b282-51c54b579871
# ╠═3eaab173-16b8-4ceb-8dba-da6d3c22d694
# ╠═5bc0af4f-337e-4485-bbf2-84d38b9b09f0
# ╠═ca7335db-11ca-40cc-baf4-ee48d9fd136d
# ╠═a8c3e6a2-c691-4b2e-a6fd-b872215fa053
# ╠═3619681e-6c74-471a-bd0e-5128a728d0d4
# ╠═92538392-a68c-4ef1-a302-a2dafbe07ba2
# ╠═5279e090-38ba-4c97-994d-c0611e19b5dc
# ╠═1178e753-b33a-4ded-a71d-86aa107311e8
# ╠═282e1e9c-c3db-4e38-bfd7-65e3a30051ca
# ╠═934d809c-2233-4201-b01b-85feddc41f31
# ╠═c7e2a059-1322-427f-ad41-4ce00fbe9d7c
# ╠═7d4ce84f-8694-45fe-b435-0819554695ce
# ╠═6019ebf8-a12b-446c-9e14-66ab340918a7
# ╠═689f4c8b-3530-422d-bc8e-65e407c45b86
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
# ╠═85a24190-8d20-4234-82a1-dd50a1891ab2
# ╠═e640afcd-92c3-4b3e-a95c-640677092f80
# ╠═8948aeb7-5b7b-49c6-978a-4d6ec91aa1bc
# ╠═665614bf-3f7c-45d5-b316-bbe9bd058603
# ╠═ef947bdd-393f-45cb-95ab-aedb2cbd76e8
# ╠═105e2342-7135-49a1-9458-0327be7b3ddb
# ╠═19e2b51d-4f4d-4443-b221-441c46ac2f37
# ╠═f651b5b2-5091-11ee-0fd9-fda723590e78
# ╠═f49001d9-c1cb-41d9-8659-aca18e5b89c7
# ╠═92cfd1aa-ec72-4137-aaf5-b15c051925eb
# ╠═72da82b1-8fb9-45a8-9130-66d7f24b668f
