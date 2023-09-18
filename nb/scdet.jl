### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 4c346534-bded-48ed-91d7-17dea3bd6933
begin
	using Geant4
	using Geant4.SystemOfUnits
	using Geant4.PhysicalConstants
	using Geant4.SystemOfUnits: mm, cm3, ns
	import CairoMakie
	using GeometryBasics, Rotations, IGLWrap_jll
end

# ╔═╡ e3411132-8907-4c6d-85e6-02460a31f331
begin
	import PlutoUI
	using FHist
	using Printf
	using Plots
end

# ╔═╡ 05b9f480-5390-11ee-2101-4b947107b1e9
md"""
# scdet

A complete example of Geant4 application in Julia. Original code at https://github.com/settwi/g4-basic-scintillation and example also available as parte of the Geant4.jl package.

It defines a scintillation setup, including the propagation of optical photons
"""

# ╔═╡ fa607f1f-b43b-4e12-bb4b-ebf18cc651c9
md"""
- Include the Geant4 specific packages
"""

# ╔═╡ 57cc79b3-7a00-4cd1-a935-3c7380870b8c
md"""
- Include other packages 
"""

# ╔═╡ 5ffed418-69ab-40ee-be61-5582f20df6e8
PlutoUI.TableOfContents(title="Scdet example (PM, JJGC)", indent=true)

# ╔═╡ 0fdc6a21-97bc-47e1-bf25-9f7f27defa85
#include(joinpath(@__DIR__, "scindet.jl"))

# ╔═╡ b73b55f8-883d-466a-bfa5-c2ff42cb1167
md"""
# Application code
"""

# ╔═╡ 4273bde0-b5fe-4e7c-8c2a-65b38a96ca2e
md"""
### Build detector
"""

# ╔═╡ 9555f147-771e-4aa5-8fe7-c784a290cebd
md"""
### Define particle gun
"""

# ╔═╡ 5a09e358-5955-4628-86e5-44f571933876
particlegun = G4JLGunGenerator(particle = "gamma", 
                               energy = 30keV, 
                               direction = G4ThreeVector(0, 0, -1), 
                               position = G4ThreeVector(0, 0, 2cm))


# ╔═╡ 973eac13-e50d-4da8-b284-31475d9853c5
md"""
### Create the application
"""

# ╔═╡ dce6c5ed-66f2-447f-ba12-5155cb460948
md"""
### Configure, Initialize and Run
"""

# ╔═╡ fae57a76-f596-4fd6-b447-818e9bb6c446
md"""
# Functions
"""

# ╔═╡ 874dc33d-bb41-43f2-8b43-b25811eb08d1
md"""
## G4 application
"""

# ╔═╡ ef5d20a8-d572-4309-9380-e8b17b9029ae
md"""
### Analysis functions
"""

# ╔═╡ bf408fe1-9995-4699-84b6-2868c16e69e8
md"""
### Simulation Data (normally filled by actions)
"""

# ╔═╡ 7787a486-eb86-450a-a98e-b332d39dcbda
const Hist1D64 = Hist1D{Float64, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}


# ╔═╡ a2a1eddc-1d75-4026-b69f-09093b3baf63
mutable struct ScintSimData <: G4JLSimulationData
    scintPhotonsPerEvent::Int64
    scintPhotonsHisto::Hist1D64
    siHitsHisto::Hist1D64
    crysEdepHisto::Hist1D64
    ScintSimData() = new()
end

# ╔═╡ c84b2410-e33a-4d69-a96c-813bf931a523
"""
Adds ScintSimData

"""
function add!(x::ScintSimData, y::ScintSimData)
    x.scintPhotonsHisto += y.scintPhotonsHisto
    x.siHitsHisto += y.siHitsHisto
    x.crysEdepHisto += y.crysEdepHisto
end


# ╔═╡ 920e3820-b4c2-4132-9b93-6625514aec87
"""
Plot Simulation data
"""
function do_plot(data::ScintSimData)
    (;scintPhotonsHisto, siHitsHisto, crysEdepHisto) = data
    lay = @layout [° °; °]
    plot(layout=lay, show=true, size=(1400,1000))
    plot!(subplot=1, scintPhotonsHisto, title="# scintillating photons/event", 
		  xlabel="# photons", show=true)
    plot!(subplot=2, siHitsHisto, title="# hits in Silicon/event", 
		  xlabel="# hits")
    plot!(subplot=3, crysEdepHisto, title="energy deposited in crystal", 
		  xlabel="keV")
end

# ╔═╡ 670a8f07-cd3d-44c8-aa76-577814f02495
md"""
### Physics list
"""

# ╔═╡ ab3b41d2-e043-4afa-bc44-55694f6cc403
"""
Define physics list
"""
struct ScintPhysicsList <: G4VUserPhysicsList
    function ScintPhysicsList(verbose)
        pl = FTFP_BERT(verbose)
        ReplacePhysics(pl, move!(G4EmStandardPhysics_option4(verbose)))
        RegisterPhysics(pl, move!(G4OpticalPhysics(verbose)))
        # need to enable scintillation
        optpar = G4OpticalParameters!Instance()
        SetProcessActivation(optpar, "Scintillation", true);
        # I have found Cherenkov radiation to be error-prone
        SetProcessActivation(optpar, "Cerenkov", true);
        return pl
    end 
end

# ╔═╡ 1371985e-2392-4500-bc67-d9ed60ac6117
md"""
### Detector definition
"""

# ╔═╡ 40af1a84-df87-461c-a656-25b986641328
"""
Constructs the struct defining the detector. 
"""
mutable struct ScintDetector <: G4JLDetector
    # main input parameters
    const worldSize::Float64
    const doTeflon::Bool
    const doAluminum::Bool
    const crys_size::Float64
    const si_side::Float64
    const si_thick::Float64
    const tef_thick::Float64
    const al_thick::Float64
    const air_gap::Float64
    const checkOverlaps::Bool
    # calculated variables
    vacMat::CxxPtr{G4Material}
    cebr3Mat::CxxPtr{G4Material}
    qzMat::CxxPtr{G4Material}
    alMat::CxxPtr{G4Material}
    tefMat::CxxPtr{G4Material}
    siMat::CxxPtr{G4Material}
    beMat::CxxPtr{G4Material}
    # constructor with defaults values for parameters
    function ScintDetector(; worldSize = 5cm,
                             doTeflon = true,
                             doAluminum = true,
                             crys_size = 2cm,
                             si_side = 1cm,
                             si_thick = 2mm,
                             tef_thick = 1mm,
                             al_thick = 1mm,
                             air_gap = 10mm,
                             checkOverlaps = true )
        new(worldSize, doTeflon, doAluminum, crys_size, si_side, 
			       si_thick, tef_thick, al_thick, air_gap, checkOverlaps)
    end
end


# ╔═╡ 3e47feeb-9ac4-4e19-86ad-42c944614a84
crysdet = ScintDetector()

# ╔═╡ fb8ce755-a055-468f-bf07-214d852b65d6
"""
Pass the function to the constructor
"""
#Geant4.getConstructor(::ScintDetector)::Function = scintConstruct

# ╔═╡ 433d8de4-cf0e-4e3c-aa4c-99fcdd6562b6
md"""
### User Actions
"""

# ╔═╡ 2c5ade51-5d9c-4a27-ba40-b23a24ddd74a
md"""
- Enumerates Hit types. 
- Define two type of hits corresponding to the scintillation cristal and the silicon
    detector
"""

# ╔═╡ 8f27a4aa-634f-4185-a2be-72d9e025eb20

@enum HitType Si ScintCryst


# ╔═╡ 2cb2f53b-adda-4b6c-9059-12c3b36a9a60
"""
Defines a hit:

- The fields are arrival time, deposited energy, hit type (si or cristal) and position
"""
struct Hit
    arrivalTime::Float64
    depositedEnergy::Float64
    hittype::HitType
    position::G4ThreeVector
    function Hit(time, pos, edep, typ)
        hit = new(time, edep, typ, G4ThreeVector())
        assign(hit.position, pos)
        return hit
    end
end

# ╔═╡ 736418a2-2a6d-42f6-95fd-eb3a41a4b372
md"""
- Define Crystal Sensitive Detector & Silicon Sensitive Detector 
"""

# ╔═╡ a8d80ee8-15c3-47a4-b94a-2e3d8d013b16
"""
Crystal Sensitive Detector SD collected data
"""
struct CrystalSDData <: G4JLSDData
    hitcollection::Vector{Hit}
    CrystalSDData() = new(Hit[])
end

# ╔═╡ c289d0a6-50b4-4645-9d8a-61fd783e6491
"""
Silicon Sensitive Detector Sd collected data
"""
struct SiSDData <: G4JLSDData
    hitcollection::Vector{Hit}
    SiSDData() = new(Hit[])
end

# ╔═╡ ef5ec36c-e282-485a-8fc3-433f4626f35a
"""
Init method to process hits in crytal (SD)
- Empty the hit collection 
"""
function crystal_initialize(::G4HCofThisEvent, data::CrystalSDData)::Nothing
    empty!(data.hitcollection)
    return
end

# ╔═╡ 52afa901-e0d3-4299-8011-cf96d2155b0a
md"""
**crystal_processHits**
- Process hits in crystal 
- Ignore optical photons. For all other particles store a hit
"""

# ╔═╡ 567fca33-3bf8-472c-b422-cd1e0a08562f

let G4OpticalPhoton, first=true
    global function crystal_processHits(step::G4Step, ::G4TouchableHistory,
		                                data::CrystalSDData)::Bool
        if first
            G4OpticalPhoton = FindParticle("opticalphoton")
            first = false
        end
        part = step |> GetTrack |> GetParticleDefinition
        part == G4OpticalPhoton && return false 
        edep = step |> GetTotalEnergyDeposit
        edep <  0. && return false
        pos = step |> GetPostStepPoint |> GetPosition
        push!(data.hitcollection, Hit(0., pos, edep, ScintCryst))
        return true
    end
end

# ╔═╡ d6b7ef9d-8d26-439c-af4f-c1a5abe8633a
md"""
- Create SD instance for crystal 
"""

# ╔═╡ e5531168-e245-4d32-9e47-fc26b6d78aab
crystal_SD = G4JLSensitiveDetector("Crystal_SD", 
	         CrystalSDData();      # SD name an associated data are mandatory
             processhits_method=crystal_processHits,     # method also mandatory
             initialize_method=crystal_initialize)       # intialize method


# ╔═╡ cbc8085c-e0ff-471f-9fed-a86eb88efa55
"""
Init method to process hits in silicon (SD)
"""
function si_initialize(::G4HCofThisEvent, data::SiSDData)::Nothing
    empty!(data.hitcollection)
    return
end

# ╔═╡ 7d31ec0c-6bf5-4da3-97db-4419da2d12ba
"""
Process hits in silicon (do nothing)
"""
function si_processHits(step::G4Step, ::G4TouchableHistory, data::SiSDData)::Bool
    return false
end

# ╔═╡ 3b63f6c8-cf82-4e3d-9714-0b9e7fb9e32b
md"""
- Create SD instance for silicon
"""

# ╔═╡ 59c1a1e1-f6fc-481a-b7e5-b9fab665f9f3
silicon_SD = G4JLSensitiveDetector("Silicon_SD", 
	                                SiSDData();               
                                    processhits_method=si_processHits,     
                                    initialize_method=si_initialize)       

# ╔═╡ b84f3ccc-6cb8-40fa-9d6d-d93090d1e620
"""
Process optical photons (hit also includes time)
"""
function processOptical(step::G4Step, data::SiSDData)
    edep = step |> GetTotalEnergyDeposit
    edep == 0 && return
    pos =  step |> GetPostStepPoint |> GetPosition
    t = step |> GetPostStepPoint |> GetGlobalTime
    push!(data.hitcollection, Hit(t, pos, edep, Si))
    return
end

# ╔═╡ c27b3002-4711-404e-9b2c-7c76d66364a5
md"""
#### User actions
"""

# ╔═╡ ba9c2dde-5e31-4778-8c14-60d6d1669d30
"""
Begin Run Action (book histograms, old HBOOK style reborn in package FHist)
"""
function beginrun(run::G4Run, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    # init histos
    data.scintPhotonsHisto = Hist1D(;bins=2.0:0.5:50.0)
    data.siHitsHisto = Hist1D(;bins=20.0:20.0:2000.0)
    data.crysEdepHisto = Hist1D(;bins=28.0:0.05:32.0)
    nothing
end

# ╔═╡ 8a77e357-cd2d-4183-b0ab-4fef8e6282eb
"""
End Run Action (handle threads if relevant)
"""
function endrun(run::G4Run, app::G4JLApplication)::Nothing
    #---end run action is called for each workwer thread and the master one
    if G4Threading!G4GetThreadId() == -1   
        data = app.simdata[1]
    #---This is the master thread, so we need to add all the simuation results
		
        for d in app.simdata[2:end]
            add!(data, d)
        end
        nEvt = GetNumberOfEvent(run)
        G4JL_println("Run ended with a total of $nEvt events") 
    end
end

# ╔═╡ 34f1d993-450e-4cea-803a-b733fff9dcdc
"""
Begin Event Action (clear number of scint photons per event)
"""
function beginevent(evt::G4Event, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    data.scintPhotonsPerEvent = 0
    nothing
end

# ╔═╡ 857cdecb-27a2-4dec-a02d-d6deef70c869
"""
End Event Action
- Store sinctillation photons
- get and store hit collection for crystal and silicon 
"""
function endevent(evt::G4Event, app::G4JLApplication)
    data = getSIMdata(app)
    push!(data.scintPhotonsHisto, data.scintPhotonsPerEvent)
    
	# get the data from the SDs
    cryshits = getSDdata(app, "Crystal_SD").hitcollection
    sihits = getSDdata(app, "Silicon_SD").hitcollection
    push!(data.siHitsHisto, length(sihits))
    push!(data.scintPhotonsHisto, length(cryshits))
    push!(data.crysEdepHisto, isempty(cryshits) ? 0.0 : sum(
		  (h.depositedEnergy for h in cryshits))/keV)
    nothing
end

# ╔═╡ 483f9846-3568-4817-86de-36ee556666d6
md"""
- User Stepping action
- Count the number of scintillation photons per event
- For optical photons:
- Find if step in geometrical boundary and the status is "Detection"
- In this case, call processOptical and get the hit for the optical photon
"""

# ╔═╡ 1b55751d-4aff-41f6-b022-6b7f1e874f5f
"""
Finds optical boundaries
"""
function findOpticalBoundary(step::G4Step)::CxxPtr{G4OpBoundaryProcess}
    pv  = step |> GetTrack |> GetDefinition |> GetProcessManager |> GetProcessList
    for p in pv
        if p |> GetProcessName |> String == "OpBoundary"
            return CxxPtr{G4OpBoundaryProcess}(p)
        end
    end
    # G4JL_println("ImpSteppingAction/findOpticalBoundary: issue finding optical boundary\n")
    return CxxPtr{G4OpBoundaryProcess}(C_NULL)
end


# ╔═╡ 54b88049-5f28-4029-8897-a9aaed82f954
let G4OpticalPhoton, first=true
    global function stepping(step::G4Step, app::G4JLApplication)::Nothing
        if first
            G4OpticalPhoton = FindParticle("opticalphoton")
            first = false
        end

        boundary = findOpticalBoundary(step)
        boundary == C_NULL && return

        data = getSIMdata(app)
        #---trackScintillation
        # secs = step |> GetSecondaryInCurrentStep # not wrapped (missing definition of type const std::vector<const G4Track *> *)
        nsec = step |> GetNumberOfSecondariesInCurrentStep
        if nsec > 0
            allsecs =  step |> GetSecondary
            lsec = length(allsecs)
            secs = allsecs[][lsec-nsec+1:lsec]
            # store the number of optical photons generated 
			# to be saved at end of event
            numOps = count(t -> t |> GetParticleDefinition == G4OpticalPhoton, secs)
            data.scintPhotonsPerEvent += numOps
        end
        #---process optical
        track = step |> GetTrack
        if track |> GetParticleDefinition == G4OpticalPhoton
            prePV = step |> GetPreStepPoint |> GetPhysicalVolume
            postPV = step |> GetPostStepPoint |> GetPhysicalVolume
            preName = prePV == C_NULL ? "" : prePV |> GetName |> String
            postName = postPV == C_NULL ? "" : postPV |> GetName |> String

            yesSilicon = occursin("si", preName) || occursin("si", postName)

            if step |> GetPostStepPoint |> GetStepStatus == fGeomBoundary
                stat = boundary |> GetStatus
                if stat == Detection
                    sddata = getSDdata(app, "Silicon_SD")
                    processOptical(step, sddata)
                    return
                elseif stat in (Absorption, TotalInternalReflection, 
				                StepTooSmall, Transmission, 
				                FresnelRefraction, FresnelReflection)
                    return
                elseif stat == NoRINDEX
                    G4JL_println("NoRINDEX: \n pre vol = $preName \n post vol = $postName")
                else
                    if yesSilicon
                        G4JL_println("something weird: $stat \n pre volume $preName \n post volume $postName")
                        return
                    end
                end
            end
        end
    end
end

# ╔═╡ 3953460d-d0ae-4bb6-95d0-4d57daf8100d
app = G4JLApplication(; detector = crysdet,          # detector with parameters
                        simdata = ScintSimData(),    # simulation data structure
                        generator = particlegun,     # primary particle generator
                        nthreads = 0,                # # of threads (0 = no MT)
                        physics_type = ScintPhysicsList,   # what physics list 
                        stepaction_method = stepping,      # step action
                        begineventaction_method = beginevent,  # init per-event data
						endeventaction_method = endevent,  
                        beginrunaction_method = beginrun,   # begin run action
                        endrunaction_method = endrun,       # end run action
                        sdetectors = ["si_log" => silicon_SD,
                                      "cebr3_log" => crystal_SD]    # mapping of LVs to SDs (+ means multiple LVs with same name)
                      );

# ╔═╡ b92515e3-fc13-4776-9d2c-e2a8c843a24f
begin
	configure(app)
	
end

# ╔═╡ 2853e73b-3e20-4c53-ab0a-1e467bb5cbd6
initialize(app)

# ╔═╡ 06e633da-dfee-4a84-b11e-1ce376c076eb
beamOn(app,1000)

# ╔═╡ 79d7aab9-6f9f-44b1-9fd9-f57617a17672
md"""
## Utility functions
"""

# ╔═╡ aa6c726f-f36b-4cd5-a892-f36b160bafe8
"""
Returns the root dir on top of the notebook (nb) dir
"""
function findrdir()
	nbdir = split(@__DIR__,"/")
	reduce(string, [string(x,"/") for x in nbdir[1:end-1]])
end

# ╔═╡ ed6d6e22-93a7-4825-a018-e26b30f08539
using Pkg; Pkg.activate(findrdir())

# ╔═╡ 0d98b72e-2eb2-44c5-a474-75a04c8656bf
md"""
## G4 application parameters
"""

# ╔═╡ dbc332f0-ec5b-4747-8dd1-5ec0c2710bff
begin
	const useSpline = true
	const VACUUM_TEMPERATURE = 2.73kelvin
	const VACUUM_PRESSURE = 1.322e-11pascal
	const SATELLITE_TEMP = 283kelvin
	
	const BR_MASS_FRAC = 0.631108
	const CE_MASS_FRAC = 1 - BR_MASS_FRAC
	
	const AL_REFR_IDX_ENERGIES = StdVector(Float64[1.9997*eV, 2.5200*eV, 3.0020*eV, 3.5024*eV, 4.0255*eV])
	const AL_REFR_IDX_REAL = StdVector(Float64[0.798, 0.465, 0.333, 0.255, 0.202])
	const AL_REFR_IDX_IMAG = StdVector(Float64[5.495, 4.71, 3.933, 3.341, 2.865])
	
	const CEBR3_DENSITY = 5.1 * g/cm3
	const CEBR3_SCINT_RESLN_SCALE = 1
	# Quarati et al (not encapsulated)
	const CEBR3_SCINT_YIELD = 60000 / MeV
	const CEBR3_DECAY_TIME_CONSTANT = 20 * ns
	
	# Quarati et al, 2.5mm thickness
	# energies of cebr3 emission spectrum
	const CEBR3_SCINT_OPTICAL_ENERGIES = StdVector(Float64[
	    2.58558e+00*eV, 2.60375e+00*eV, 2.62735e+00*eV, 2.65348e+00*eV, 2.68122e+00*eV, 2.71175e+00*eV, 2.74185e+00*eV, 2.77033e+00*eV,
	    2.79940e+00*eV, 2.83029e+00*eV, 2.86188e+00*eV, 2.88540e+00*eV, 2.90805e+00*eV, 2.93235e+00*eV, 2.94140e+00*eV, 2.95443e+00*eV,
	    2.96758e+00*eV, 2.98485e+00*eV, 2.99962e+00*eV, 3.01317e+00*eV, 3.03098e+00*eV, 3.04899e+00*eV, 3.06159e+00*eV, 3.07005e+00*eV,
	    3.08139e+00*eV, 3.09426e+00*eV, 3.10724e+00*eV, 3.12032e+00*eV, 3.12911e+00*eV, 3.13942e+00*eV, 3.14683e+00*eV, 3.15726e+00*eV,
	    3.16475e+00*eV, 3.17228e+00*eV, 3.17984e+00*eV, 3.19049e+00*eV, 3.20429e+00*eV, 3.21976e+00*eV, 3.23224e+00*eV, 3.24482e+00*eV,
	    3.25274e+00*eV, 3.25910e+00*eV, 3.27189e+00*eV, 3.28155e+00*eV, 3.28802e+00*eV, 3.29941e+00*eV, 3.31417e+00*eV, 3.32574e+00*eV,
	    3.33572e+00*eV, 3.34576e+00*eV, 3.35249e+00*eV, 3.36094e+00*eV, 3.36943e+00*eV, 3.37967e+00*eV, 3.38998e+00*eV, 3.40035e+00*eV,
	    3.41079e+00*eV, 3.42129e+00*eV, 3.43185e+00*eV, 3.43893e+00*eV, 3.44782e+00*eV, 3.45676e+00*eV, 3.46574e+00*eV, 3.47477e+00*eV,
	    3.48385e+00*eV, 3.49480e+00*eV, 3.50767e+00*eV, 3.52063e+00*eV, 3.53556e+00*eV, 3.55820e+00*eV, 3.57729e+00*eV, 3.59658e+00*eV,
	    3.62001e+00*eV, 3.64375e+00*eV, 3.67588e+00*eV, 3.70036e+00*eV, 3.75876e+00*eV, 3.81466e+00*eV, 3.86776e+00*eV
	])
	
	# relative intensities of emission spectrum
	const CEBR3_SCINT_OPTICAL_INTENSITIES = StdVector(Float64[
	    3.95471e-04, 4.54059e-04, 4.54059e-04, 4.54059e-04, 4.54059e-04, 5.71236e-04, 6.88413e-04, 7.47001e-04,
	    8.64178e-04, 9.81354e-04, 1.15712e-03, 1.39147e-03, 1.74300e-03, 2.09453e-03, 2.32889e-03, 2.79759e-03,
	    3.26630e-03, 4.02795e-03, 4.84818e-03, 5.66842e-03, 7.13313e-03, 8.77360e-03, 1.00625e-02, 1.11757e-02,
	    1.24647e-02, 1.40465e-02, 1.55113e-02, 1.68002e-02, 1.77376e-02, 1.86750e-02, 1.93781e-02, 2.01397e-02,
	    2.06670e-02, 2.11357e-02, 2.15459e-02, 2.17802e-02, 2.20731e-02, 2.23661e-02, 2.27176e-02, 2.31277e-02,
	    2.34207e-02, 2.35964e-02, 2.38894e-02, 2.40652e-02, 2.43581e-02, 2.50026e-02, 2.59400e-02, 2.67602e-02,
	    2.74633e-02, 2.78734e-02, 2.81663e-02, 2.83421e-02, 2.82835e-02, 2.81663e-02, 2.76390e-02, 2.67016e-02,
	    2.54127e-02, 2.40066e-02, 2.18974e-02, 2.03155e-02, 1.83235e-02, 1.63315e-02, 1.43981e-02, 1.25233e-02,
	    1.08242e-02, 9.06654e-03, 7.19172e-03, 5.66842e-03, 4.32089e-03, 2.97336e-03, 1.86018e-03, 1.15712e-03,
	    6.88413e-04, 3.36883e-04, 1.61118e-04, 4.39412e-05, 4.39412e-05, 4.39412e-05, 4.39412e-05
	])


	# 2022 February 04: refractive index experimentially is about 2.09 for the emisison range, and relatively constant.
	# update to reflect that here.
	const CEBR3_REFR_IDX_ENERGIES = StdVector(Float64[ 1e-3*eV, 5*eV ])
	const CEBR3_REFR_IDXS = StdVector(Float64[ 2.09, 2.09 ])
	
	const CEBR3_ABS_LEN_ENERGIES = StdVector(Float64[
	    1.00000e+00*eV, 1.18329e+00*eV, 2.06497e+00*eV, 2.38979e+00*eV, 2.71462e+00*eV, 2.94664e+00*eV, 3.17865e+00*eV, 3.68910e+00*eV,
	    4.10673e+00*eV, 4.38515e+00*eV, 4.66357e+00*eV, 4.98840e+00*eV, 5.17401e+00*eV, 5.31323e+00*eV, 5.45244e+00*eV, 5.68445e+00*eV,
	    5.96288e+00*eV, 6.14849e+00*eV, 6.28770e+00*eV, 6.47332e+00*eV, 6.61253e+00*eV, 6.84455e+00*eV, 7.03016e+00*eV, 7.12297e+00*eV,
	    7.30858e+00*eV, 7.49420e+00*eV, 7.58701e+00*eV, 7.67981e+00*eV, 7.77262e+00*eV, 7.86543e+00*eV, 8.00464e+00*eV, 8.05104e+00*eV,
	    8.19026e+00*eV, 8.32947e+00*eV, 8.46868e+00*eV, 8.56148e+00*eV, 8.70070e+00*eV, 8.79350e+00*eV, 8.93271e+00*eV, 9.11833e+00*eV,
	    9.35035e+00*eV, 9.48956e+00*eV, 9.67517e+00*eV, 9.76798e+00*eV, 9.86079e+00*eV, 1.00000e+01*eV, 1.01856e+01*eV, 1.03248e+01*eV,
	    1.04640e+01*eV, 1.06497e+01*eV, 1.08817e+01*eV, 1.12529e+01*eV, 1.14849e+01*eV, 1.16705e+01*eV, 1.19026e+01*eV, 1.22274e+01*eV,
	    1.25986e+01*eV, 1.27378e+01*eV, 1.28770e+01*eV, 1.30626e+01*eV, 1.32019e+01*eV, 1.32947e+01*eV, 1.34339e+01*eV, 1.35731e+01*eV,
	    1.38051e+01*eV, 1.40139e+01*eV, 1.42691e+01*eV, 1.45940e+01*eV, 1.47332e+01*eV, 1.49188e+01*eV, 1.51972e+01*eV, 1.53364e+01*eV,
	    1.54756e+01*eV, 1.56613e+01*eV, 1.58005e+01*eV, 1.60325e+01*eV, 1.62181e+01*eV, 1.64037e+01*eV, 1.65893e+01*eV, 1.68213e+01*eV,
	    1.69606e+01*eV, 1.71694e+01*eV, 1.74710e+01*eV, 1.77494e+01*eV, 1.83063e+01*eV, 1.86311e+01*eV, 1.90487e+01*eV, 1.92343e+01*eV,
	    1.94664e+01*eV, 1.97912e+01*eV, 2.01856e+01*eV, 2.03480e+01*eV, 2.05336e+01*eV, 2.06729e+01*eV, 2.08121e+01*eV, 2.09977e+01*eV,
	    2.10441e+01*eV, 2.11833e+01*eV, 2.13689e+01*eV, 2.16009e+01*eV, 2.17865e+01*eV, 2.20650e+01*eV, 2.22970e+01*eV, 2.24362e+01*eV,
	    2.25290e+01*eV, 2.26682e+01*eV, 2.28074e+01*eV
	])
	
	# absorption lengths are super long so just make them 1 km
	const CEBR3_ABS_LEN = StdVector(fill(1.0km, length(CEBR3_ABS_LEN_ENERGIES)))
	
	const QZ_REFR_IDX_ENERGIES = StdVector(Float64[
	        1.85051e-01*eV, 1.91570e-01*eV, 1.98311e-01*eV, 2.05306e-01*eV, 2.12556e-01*eV, 2.20025e-01*eV, 2.27786e-01*eV, 2.35801e-01*eV,
	        2.44111e-01*eV, 2.52720e-01*eV, 2.61625e-01*eV, 2.70826e-01*eV, 2.80380e-01*eV, 2.90293e-01*eV, 3.00495e-01*eV, 3.11049e-01*eV,
	        3.22037e-01*eV, 3.33380e-01*eV, 3.45168e-01*eV, 3.57303e-01*eV, 3.69881e-01*eV, 3.82904e-01*eV, 3.96369e-01*eV, 4.10408e-01*eV,
	        4.24894e-01*eV, 4.39816e-01*eV, 4.55322e-01*eV, 4.71423e-01*eV, 4.87935e-01*eV, 5.05233e-01*eV, 5.22919e-01*eV, 5.41416e-01*eV,
	        5.60507e-01*eV, 5.80179e-01*eV, 6.00699e-01*eV, 6.21786e-01*eV, 6.43739e-01*eV, 6.66223e-01*eV, 6.89951e-01*eV, 7.14195e-01*eV,
	        7.39321e-01*eV, 7.65334e-01*eV, 7.92231e-01*eV, 8.20001e-01*eV, 8.49207e-01*eV, 8.79320e-01*eV, 9.10310e-01*eV, 9.42129e-01*eV,
	        9.75485e-01*eV, 1.00964e+00*eV, 1.04540e+00*eV, 1.08189e+00*eV, 1.12000e+00*eV, 1.15981e+00*eV, 1.20023e+00*eV, 1.24282e+00*eV,
	        1.28668e+00*eV, 1.33202e+00*eV, 1.37883e+00*eV, 1.42740e+00*eV, 1.47776e+00*eV, 1.52991e+00*eV, 1.58365e+00*eV, 1.63957e+00*eV,
	        1.69725e+00*eV, 1.75715e+00*eV, 1.81902e+00*eV, 1.88311e+00*eV, 1.94944e+00*eV, 2.01830e+00*eV, 2.08939e+00*eV, 2.16302e+00*eV,
	        2.23919e+00*eV, 2.31789e+00*eV, 2.39954e+00*eV, 2.48416e+00*eV, 2.57175e+00*eV, 2.66232e+00*eV, 2.75643e+00*eV, 2.85349e+00*eV,
	        2.95411e+00*eV, 3.05756e+00*eV, 3.16528e+00*eV, 3.27740e+00*eV, 3.39218e+00*eV, 3.51230e+00*eV, 3.63590e+00*eV, 3.76394e+00*eV,
	        3.89642e+00*eV, 4.03332e+00*eV, 4.17596e+00*eV, 4.32302e+00*eV, 4.47596e+00*eV, 4.63319e+00*eV, 4.79629e+00*eV, 4.96533e+00*eV,
	        5.14031e+00*eV, 5.32121e+00*eV, 5.50796e+00*eV, 5.70304e+00*eV, 5.90401e+00*eV
	])
	
	const QZ_REFR_IDXS = StdVector(Float64[
	        1.15965e+00, 1.19733e+00, 1.22809e+00, 1.25373e+00, 1.27537e+00, 1.29375e+00, 1.30964e+00, 1.32341e+00,
	        1.33548e+00, 1.34612e+00, 1.35553e+00, 1.36388e+00, 1.37137e+00, 1.37810e+00, 1.38412e+00, 1.38956e+00,
	        1.39450e+00, 1.39898e+00, 1.40307e+00, 1.40678e+00, 1.41017e+00, 1.41327e+00, 1.41612e+00, 1.41875e+00,
	        1.42115e+00, 1.42336e+00, 1.42540e+00, 1.42730e+00, 1.42903e+00, 1.43065e+00, 1.43214e+00, 1.43353e+00,
	        1.43482e+00, 1.43601e+00, 1.43713e+00, 1.43817e+00, 1.43915e+00, 1.44005e+00, 1.44091e+00, 1.44171e+00,
	        1.44246e+00, 1.44317e+00, 1.44384e+00, 1.44448e+00, 1.44509e+00, 1.44566e+00, 1.44621e+00, 1.44674e+00,
	        1.44725e+00, 1.44773e+00, 1.44821e+00, 1.44867e+00, 1.44912e+00, 1.44957e+00, 1.45001e+00, 1.45045e+00,
	        1.45089e+00, 1.45132e+00, 1.45177e+00, 1.45221e+00, 1.45267e+00, 1.45314e+00, 1.45362e+00, 1.45412e+00,
	        1.45463e+00, 1.45517e+00, 1.45572e+00, 1.45631e+00, 1.45693e+00, 1.45758e+00, 1.45826e+00, 1.45899e+00,
	        1.45976e+00, 1.46057e+00, 1.46144e+00, 1.46238e+00, 1.46337e+00, 1.46444e+00, 1.46558e+00, 1.46680e+00,
	        1.46812e+00, 1.46953e+00, 1.47105e+00, 1.47270e+00, 1.47447e+00, 1.47640e+00, 1.47847e+00, 1.48071e+00,
	        1.48315e+00, 1.48579e+00, 1.48868e+00, 1.49182e+00, 1.49526e+00, 1.49900e+00, 1.50310e+00, 1.50761e+00,
	        1.51257e+00, 1.51804e+00, 1.52408e+00, 1.53085e+00, 1.53836e+00
	])
	
	const TEFLON_REFR_IDX_ENERGIES = StdVector(Float64[ 0.1*eV, 100*eV ])
	const TEFLON_REFR_IDXS = StdVector(Float64[ 1.38, 1.38 ])
	# 23 February 2022
	# Teflon reflects ~97% of UV light so just make it 1 for simplicity.
	const TEFLON_REFLECTIVITY = StdVector(Float64[ 1., 1. ])
	const TEFLON_TRANSMITTANCE = StdVector(Float64[  0., 0. ])
	
	
	const SI_DENSITY = 2.33 * g/cm3
	
	const SI_DET_EFF_ENERGIES = StdVector(reverse(Float64[                 # !!!! Note that G4 expects increaing evenrgy numbers
	    3.98830e+00*eV, 3.90793e+00*eV, 3.86124e+00*eV, 3.78585e+00*eV, 3.71335e+00*eV, 3.64358e+00*eV, 3.57638e+00*eV, 3.51161e+00*eV,
	    3.44915e+00*eV, 3.42478e+00*eV, 3.40076e+00*eV, 3.36535e+00*eV, 3.33066e+00*eV, 3.29669e+00*eV, 3.25245e+00*eV, 3.19880e+00*eV,
	    3.15714e+00*eV, 3.13671e+00*eV, 3.10656e+00*eV, 3.05757e+00*eV, 3.01011e+00*eV, 2.97319e+00*eV, 2.92829e+00*eV, 2.86766e+00*eV,
	    2.82587e+00*eV, 2.78528e+00*eV, 2.75364e+00*eV, 2.71508e+00*eV, 2.68501e+00*eV, 2.63395e+00*eV, 2.59170e+00*eV, 2.55079e+00*eV,
	    2.51767e+00*eV, 2.47273e+00*eV, 2.42329e+00*eV, 2.35844e+00*eV, 2.31343e+00*eV, 2.27009e+00*eV, 2.20805e+00*eV, 2.15889e+00*eV,
	    2.09815e+00*eV, 2.04937e+00*eV, 2.01111e+00*eV, 1.97024e+00*eV, 1.92716e+00*eV, 1.87862e+00*eV, 1.82213e+00*eV, 1.78523e+00*eV,
	    1.74349e+00*eV, 1.68883e+00*eV, 1.64583e+00*eV, 1.59179e+00*eV, 1.54610e+00*eV, 1.51233e+00*eV, 1.47103e+00*eV, 1.43193e+00*eV,
	    1.40292e+00*eV, 1.38094e+00*eV
	]))
	
	const SI_DET_EFF = StdVector(fill(1., length(SI_DET_EFF_ENERGIES)))
	# set to zero so that we either detect or reflect
	const SI_TRANSMITTANCE = StdVector(fill(0., length(SI_DET_EFF_ENERGIES)))
	const SI_BROADCOM_NUMBERS = StdVector(reverse(Float64[
	    2.62014e-01, 2.70252e-01, 2.79863e-01, 3.00458e-01, 3.23799e-01, 3.49886e-01, 3.67735e-01, 3.78719e-01,
	    3.85584e-01, 3.86957e-01, 3.89703e-01, 3.93822e-01, 3.99314e-01, 4.04805e-01, 4.13043e-01, 4.24027e-01,
	    4.35011e-01, 4.41876e-01, 4.48741e-01, 4.51487e-01, 4.45995e-01, 4.36384e-01, 4.24027e-01, 4.06178e-01,
	    3.93822e-01, 3.81465e-01, 3.70481e-01, 3.59497e-01, 3.48513e-01, 3.32037e-01, 3.18307e-01, 3.04577e-01,
	    2.92220e-01, 2.78490e-01, 2.63387e-01, 2.42792e-01, 2.30435e-01, 2.19451e-01, 2.04348e-01, 1.93364e-01,
	    1.79634e-01, 1.67277e-01, 1.56293e-01, 1.46682e-01, 1.38444e-01, 1.30206e-01, 1.19222e-01, 1.09611e-01,
	    1.01373e-01, 8.76430e-02, 7.94050e-02, 6.70481e-02, 5.88101e-02, 5.33181e-02, 4.37071e-02, 3.54691e-02,
	    2.86041e-02, 2.44851e-02
	]))
	
	const SI_REFR_IDX_ENERGY = StdVector(reverse(Float64[
	    4.69993e+00*eV, 4.67512e+00*eV, 4.65057e+00*eV, 4.62455e+00*eV, 4.60053e+00*eV, 4.57506e+00*eV, 4.54988e+00*eV, 4.52497e+00*eV,
	    4.50033e+00*eV, 4.47435e+00*eV, 4.45026e+00*eV, 4.42485e+00*eV, 4.39972e+00*eV, 4.37488e+00*eV, 4.35032e+00*eV, 4.32453e+00*eV,
	    4.30053e+00*eV, 4.27532e+00*eV, 4.25040e+00*eV, 4.22433e+00*eV, 4.20001e+00*eV, 4.17455e+00*eV, 4.14940e+00*eV, 4.12456e+00*eV,
	    4.10001e+00*eV, 4.07441e+00*eV, 4.05045e+00*eV, 4.02546e+00*eV, 3.99949e+00*eV, 3.97513e+00*eV, 3.94980e+00*eV, 3.92479e+00*eV,
	    3.90010e+00*eV, 3.87451e+00*eV, 3.85044e+00*eV, 3.82549e+00*eV, 3.79970e+00*eV, 3.77540e+00*eV, 3.75028e+00*eV, 3.72549e+00*eV,
	    3.69992e+00*eV, 3.67469e+00*eV, 3.64981e+00*eV, 3.62527e+00*eV, 3.60001e+00*eV, 3.57509e+00*eV, 3.54950e+00*eV, 3.52528e+00*eV,
	    3.50040e+00*eV, 3.47489e+00*eV, 3.44976e+00*eV, 3.42498e+00*eV, 3.39962e+00*eV, 3.37464e+00*eV, 3.35002e+00*eV, 3.32486e+00*eV,
	    3.30009e+00*eV, 3.27481e+00*eV, 3.24991e+00*eV, 3.22540e+00*eV, 3.19959e+00*eV, 3.17501e+00*eV, 3.15001e+00*eV, 3.12539e+00*eV,
	    3.10038e+00*eV, 3.07500e+00*eV, 3.05004e+00*eV, 3.02474e+00*eV, 2.99986e+00*eV, 2.97467e+00*eV, 2.94990e+00*eV, 2.92485e+00*eV,
	    2.90022e+00*eV, 2.87533e+00*eV, 2.85021e+00*eV, 2.82488e+00*eV, 2.80000e+00*eV, 2.77494e+00*eV, 2.74970e+00*eV, 2.72493e+00*eV,
	    2.70000e+00*eV, 2.67496e+00*eV, 2.64980e+00*eV, 2.62512e+00*eV, 2.59979e+00*eV, 2.57496e+00*eV, 2.55007e+00*eV, 2.52514e+00*eV,
	    2.50019e+00*eV, 2.47523e+00*eV, 2.44980e+00*eV, 2.42488e+00*eV, 2.40000e+00*eV, 2.37518e+00*eV, 2.34997e+00*eV, 2.32485e+00*eV,
	    2.29984e+00*eV, 2.27494e+00*eV, 2.25017e+00*eV, 2.22513e+00*eV, 2.19986e+00*eV, 2.17516e+00*eV, 2.14989e+00*eV, 2.12484e+00*eV,
	    2.10000e+00*eV, 2.07505e+00*eV, 2.05000e+00*eV, 2.02489e+00*eV, 2.00007e+00*eV, 1.97490e+00*eV, 1.95005e+00*eV, 1.92492e+00*eV,
	    1.90014e+00*eV, 1.87514e+00*eV, 1.84996e+00*eV, 1.82491e+00*eV, 1.80000e+00*eV, 1.77501e+00*eV, 1.74995e+00*eV, 1.72512e+00*eV,
	    1.70004e+00*eV, 1.67501e+00*eV, 1.65004e+00*eV, 1.62496e+00*eV, 1.60000e+00*eV, 1.57500e+00*eV, 1.55000e+00*eV, 1.52502e+00*eV,
	    1.49993e+00*eV
	]))
	
	const SI_REFR_IDX_REAL = StdVector(reverse(Float64[
	    1.82500e+00, 1.85400e+00, 1.89900e+00, 1.95400e+00, 2.01400e+00, 2.09900e+00, 2.18800e+00, 2.29300e+00,
	    2.42600e+00, 2.57400e+00, 2.72300e+00, 2.87700e+00, 3.03700e+00, 3.23200e+00, 3.40600e+00, 3.66800e+00,
	    3.92600e+00, 4.20500e+00, 4.45400e+00, 4.65900e+00, 4.79100e+00, 4.88500e+00, 4.92800e+00, 4.98500e+00,
	    4.98600e+00, 5.00900e+00, 5.02500e+00, 5.02400e+00, 5.02900e+00, 5.03400e+00, 5.03900e+00, 5.03500e+00,
	    5.03200e+00, 5.05100e+00, 5.06800e+00, 5.06600e+00, 5.09100e+00, 5.12000e+00, 5.13500e+00, 5.15900e+00,
	    5.17700e+00, 5.21800e+00, 5.24700e+00, 5.27700e+00, 5.31000e+00, 5.36300e+00, 5.42900e+00, 5.48700e+00,
	    5.60100e+00, 5.75200e+00, 5.96700e+00, 6.25100e+00, 6.56300e+00, 6.82100e+00, 6.98900e+00, 7.00400e+00,
	    6.89900e+00, 6.73500e+00, 6.54800e+00, 6.36900e+00, 6.19100e+00, 6.04100e+00, 5.89700e+00, 5.76400e+00,
	    5.64900e+00, 5.54300e+00, 5.45000e+00, 5.35600e+00, 5.27900e+00, 5.19700e+00, 5.12600e+00, 5.05600e+00,
	    4.99400e+00, 4.93300e+00, 4.87200e+00, 4.82500e+00, 4.77000e+00, 4.72200e+00, 4.67700e+00, 4.63100e+00,
	    4.59100e+00, 4.54900e+00, 4.51000e+00, 4.47100e+00, 4.43500e+00, 4.40200e+00, 4.36800e+00, 4.33700e+00,
	    4.30900e+00, 4.28100e+00, 4.25200e+00, 4.22400e+00, 4.20000e+00, 4.17700e+00, 4.15300e+00, 4.13100e+00,
	    4.10900e+00, 4.08900e+00, 4.06800e+00, 4.04900e+00, 4.02900e+00, 4.01200e+00, 3.99500e+00, 3.97700e+00,
	    3.96000e+00, 3.94410e+00, 3.92910e+00, 3.91440e+00, 3.90000e+00, 3.88600e+00, 3.87230e+00, 3.85890e+00,
	    3.84500e+00, 3.83300e+00, 3.82050e+00, 3.80840e+00, 3.79600e+00, 3.78510e+00, 3.77390e+00, 3.76310e+00,
	    3.75300e+00, 3.74230e+00, 3.73240e+00, 3.72280e+00, 3.71600e+00, 3.70460e+00, 3.69600e+00, 3.68770e+00,
	    3.67800e+00
	]))
	
	const SI_REFR_IDX_IMAG = StdVector(reverse(Float64[
	    4.17500e+00, 4.26200e+00, 4.34100e+00, 4.43900e+00, 4.54300e+00, 4.65900e+00, 4.73500e+00, 4.85500e+00,
	    4.94100e+00, 5.03100e+00, 5.11100e+00, 5.15800e+00, 5.23300e+00, 5.28000e+00, 5.33200e+00, 5.34700e+00,
	    5.35300e+00, 5.27200e+00, 5.14200e+00, 4.94900e+00, 4.74000e+00, 4.54800e+00, 4.36400e+00, 4.20600e+00,
	    4.06700e+00, 3.94100e+00, 3.82700e+00, 3.74400e+00, 3.64500e+00, 3.57500e+00, 3.49200e+00, 3.44400e+00,
	    3.39200e+00, 3.34700e+00, 3.28600e+00, 3.25800e+00, 3.21300e+00, 3.17400e+00, 3.14600e+00, 3.11200e+00,
	    3.08300e+00, 3.05000e+00, 3.03200e+00, 3.00900e+00, 3.00500e+00, 3.00800e+00, 3.01800e+00, 3.02800e+00,
	    3.05600e+00, 3.09200e+00, 3.10500e+00, 3.03700e+00, 2.86900e+00, 2.56200e+00, 2.16300e+00, 1.75000e+00,
	    1.40400e+00, 1.06800e+00, 8.71000e-01, 6.95000e-01, 5.87000e-01, 4.86000e-01, 4.40000e-01, 3.85000e-01,
	    3.26000e-01, 2.97000e-01, 2.59000e-01, 2.32000e-01, 2.03000e-01, 1.87000e-01, 1.71000e-01, 1.54000e-01,
	    1.42000e-01, 1.30000e-01, 1.21000e-01, 1.08000e-01, 1.01000e-01, 9.50000e-02, 8.60000e-02, 8.00000e-02,
	    6.90000e-02, 6.26000e-02, 5.68000e-02, 4.82000e-02, 4.27000e-02, 3.93000e-02, 3.62000e-02, 3.29000e-02,
	    2.80000e-02, 2.40000e-02, 2.28000e-02, 1.89000e-02, 1.44000e-02, 1.22000e-02, 1.09000e-02, 1.00000e-02,
	    8.00000e-03, 7.50000e-03, 6.00000e-03, 6.40000e-03, 5.20000e-03, 4.70000e-03, 4.50000e-03, 3.80000e-03,
	    3.50000e-03, 2.60000e-03, 2.00000e-03, 2.40000e-03, 1.70000e-03, 1.30000e-03, 8.00000e-04, 1.30000e-03,
	    1.70000e-03, 2.10000e-03, 1.70000e-03, 2.40000e-03, 2.90000e-03, 2.60000e-03, 2.40000e-03, 2.40000e-03,
	    2.90000e-03, 2.60000e-03, 4.00000e-03, 4.40000e-03, 1.20000e-03, 2.90000e-03, 4.70000e-03, 5.00000e-03,
	    4.60000e-03
	]))
	
	const BE_DENSITY = 1.85 * g / cm3
	const BE_NUM_COMPONENTS = 1
	const BE_REFR_IDX_ENERGIES = StdVector(Float64[
	    2.00000e-02*eV, 2.05626e-02*eV, 2.11411e-02*eV, 2.17360e-02*eV, 2.23475e-02*eV, 2.29762e-02*eV, 2.36228e-02*eV, 2.42873e-02*eV,
	    2.49706e-02*eV, 2.56733e-02*eV, 2.63953e-02*eV, 2.71383e-02*eV, 2.79017e-02*eV, 2.86868e-02*eV, 2.94934e-02*eV, 3.03236e-02*eV,
	    3.11769e-02*eV, 3.20538e-02*eV, 3.29552e-02*eV, 3.38829e-02*eV, 3.48358e-02*eV, 3.58160e-02*eV, 3.68233e-02*eV, 3.78601e-02*eV,
	    3.89251e-02*eV, 4.00207e-02*eV, 4.11457e-02*eV, 4.23039e-02*eV, 4.34941e-02*eV, 4.47177e-02*eV, 4.59763e-02*eV, 4.72699e-02*eV,
	    4.85984e-02*eV, 4.99674e-02*eV, 5.13732e-02*eV, 5.28177e-02*eV, 5.43028e-02*eV, 5.58311e-02*eV, 5.74027e-02*eV, 5.90176e-02*eV,
	    6.06784e-02*eV, 6.23851e-02*eV, 6.41408e-02*eV, 6.59455e-02*eV, 6.77991e-02*eV, 6.97089e-02*eV, 7.16672e-02*eV, 7.36861e-02*eV,
	    7.57572e-02*eV, 7.78893e-02*eV, 8.00828e-02*eV, 8.23323e-02*eV, 8.46482e-02*eV, 8.70309e-02*eV, 8.94805e-02*eV, 9.19969e-02*eV,
	    9.45867e-02*eV, 9.72501e-02*eV, 9.99792e-02*eV, 1.02798e-01*eV, 1.05689e-01*eV, 1.08663e-01*eV, 1.11718e-01*eV, 1.14864e-01*eV,
	    1.18091e-01*eV, 1.21410e-01*eV, 1.24832e-01*eV, 1.28344e-01*eV, 1.31956e-01*eV, 1.35668e-01*eV, 1.39485e-01*eV, 1.43409e-01*eV,
	    1.47444e-01*eV, 1.51592e-01*eV, 1.55857e-01*eV, 1.60242e-01*eV, 1.64750e-01*eV, 1.69387e-01*eV, 1.74152e-01*eV, 1.79051e-01*eV,
	    1.84089e-01*eV, 1.89269e-01*eV, 1.94592e-01*eV, 2.00068e-01*eV, 2.05698e-01*eV, 2.11483e-01*eV, 2.17432e-01*eV, 2.23552e-01*eV,
	    2.29839e-01*eV, 2.36304e-01*eV, 2.42954e-01*eV, 2.49792e-01*eV, 2.56818e-01*eV, 2.64043e-01*eV, 2.71472e-01*eV, 2.79112e-01*eV,
	    2.86961e-01*eV, 2.95039e-01*eV, 3.03340e-01*eV, 3.11871e-01*eV, 3.20646e-01*eV, 3.29666e-01*eV, 3.38940e-01*eV, 3.48476e-01*eV,
	    3.58284e-01*eV, 3.68365e-01*eV, 3.78728e-01*eV, 3.89385e-01*eV, 4.00336e-01*eV, 4.11593e-01*eV, 4.23183e-01*eV, 4.35078e-01*eV,
	    4.47322e-01*eV, 4.59916e-01*eV, 4.72843e-01*eV, 4.86155e-01*eV, 4.99835e-01*eV, 5.13903e-01*eV, 5.28357e-01*eV, 5.43218e-01*eV,
	    5.58512e-01*eV, 5.74214e-01*eV, 5.90373e-01*eV, 6.06992e-01*eV, 6.24071e-01*eV, 6.41607e-01*eV, 6.59666e-01*eV, 6.78213e-01*eV,
	    6.97324e-01*eV, 7.16920e-01*eV, 7.37080e-01*eV, 7.57850e-01*eV, 7.79138e-01*eV, 8.01087e-01*eV, 8.23596e-01*eV, 8.46771e-01*eV,
	    8.70614e-01*eV, 8.95128e-01*eV, 9.20310e-01*eV, 9.46155e-01*eV, 9.72806e-01*eV, 1.00020e+00*eV, 1.02832e+00*eV, 1.05725e+00*eV,
	    1.08701e+00*eV, 1.11758e+00*eV, 1.14896e+00*eV, 1.18136e+00*eV, 1.21458e+00*eV, 1.24875e+00*eV, 1.28388e+00*eV, 1.31999e+00*eV,
	    1.35714e+00*eV, 1.39532e+00*eV, 1.43457e+00*eV, 1.47493e+00*eV, 1.51642e+00*eV, 1.55910e+00*eV, 1.60296e+00*eV, 1.64805e+00*eV,
	    1.69442e+00*eV, 1.74208e+00*eV, 1.79111e+00*eV, 1.84150e+00*eV, 1.89332e+00*eV, 1.94659e+00*eV, 2.00136e+00*eV, 2.05766e+00*eV,
	    2.11555e+00*eV, 2.17505e+00*eV, 2.23625e+00*eV, 2.29915e+00*eV, 2.36385e+00*eV, 2.43035e+00*eV, 2.49872e+00*eV, 2.56903e+00*eV,
	    2.64133e+00*eV, 2.71562e+00*eV, 2.79206e+00*eV, 2.87060e+00*eV, 2.95137e+00*eV, 3.03437e+00*eV, 3.11973e+00*eV, 3.20754e+00*eV,
	    3.29780e+00*eV, 3.39051e+00*eV, 3.48593e+00*eV, 3.58398e+00*eV, 3.68485e+00*eV, 3.78855e+00*eV, 3.89508e+00*eV, 4.00466e+00*eV,
	    4.11743e+00*eV, 4.23328e+00*eV, 4.35231e+00*eV, 4.47483e+00*eV, 4.60070e+00*eV, 4.73005e+00*eV, 4.86327e+00*eV, 4.99997e+00*eV
	])
	
	const BE_REFR_IDX_REAL = StdVector(Float64[
	    9.40750e+01, 9.15210e+01, 8.89980e+01, 8.65080e+01, 8.40500e+01, 8.16260e+01, 7.92350e+01, 7.68790e+01,
	    7.45580e+01, 7.22730e+01, 7.00250e+01, 6.78140e+01, 6.56410e+01, 6.35060e+01, 6.14110e+01, 5.93560e+01,
	    5.73420e+01, 5.53690e+01, 5.34370e+01, 5.15480e+01, 4.97020e+01, 4.78990e+01, 4.61390e+01, 4.44230e+01,
	    4.27510e+01, 4.11230e+01, 3.95400e+01, 3.80000e+01, 3.65050e+01, 3.50540e+01, 3.36480e+01, 3.22840e+01,
	    3.09650e+01, 2.96880e+01, 2.84540e+01, 2.72610e+01, 2.61110e+01, 2.50010e+01, 2.39320e+01, 2.29020e+01,
	    2.19100e+01, 2.09570e+01, 2.00410e+01, 1.91610e+01, 1.83160e+01, 1.75060e+01, 1.67300e+01, 1.59860e+01,
	    1.52740e+01, 1.45920e+01, 1.39400e+01, 1.33170e+01, 1.27220e+01, 1.21540e+01, 1.16110e+01, 1.10940e+01,
	    1.06010e+01, 1.01310e+01, 9.68310e+00, 9.25680e+00, 8.85100e+00, 8.46510e+00, 8.09800e+00, 7.74920e+00,
	    7.41770e+00, 7.10290e+00, 6.80410e+00, 6.52050e+00, 6.25150e+00, 5.99640e+00, 5.75460e+00, 5.52560e+00,
	    5.30870e+00, 5.10340e+00, 4.90910e+00, 4.72540e+00, 4.55170e+00, 4.38770e+00, 4.23280e+00, 4.08660e+00,
	    3.94870e+00, 3.81870e+00, 3.69630e+00, 3.58110e+00, 3.47280e+00, 3.37100e+00, 3.27540e+00, 3.18580e+00,
	    3.10190e+00, 3.02340e+00, 2.95010e+00, 2.88170e+00, 2.81810e+00, 2.75890e+00, 2.70410e+00, 2.65350e+00,
	    2.60680e+00, 2.56380e+00, 2.52460e+00, 2.48880e+00, 2.45640e+00, 2.42730e+00, 2.40130e+00, 2.37830e+00,
	    2.35820e+00, 2.34090e+00, 2.32640e+00, 2.31450e+00, 2.30520e+00, 2.29840e+00, 2.29400e+00, 2.29210e+00,
	    2.29240e+00, 2.29510e+00, 2.30000e+00, 2.30700e+00, 2.31630e+00, 2.32760e+00, 2.34100e+00, 2.35650e+00,
	    2.37400e+00, 2.39340e+00, 2.41470e+00, 2.43790e+00, 2.46290e+00, 2.48960e+00, 2.51800e+00, 2.54810e+00,
	    2.57960e+00, 2.61250e+00, 2.64660e+00, 2.68200e+00, 2.71830e+00, 2.75550e+00, 2.79330e+00, 2.83160e+00,
	    2.87020e+00, 2.90890e+00, 2.94740e+00, 2.98560e+00, 3.02330e+00, 3.06020e+00, 3.09610e+00, 3.13090e+00,
	    3.16440e+00, 3.19640e+00, 3.22680e+00, 3.25540e+00, 3.28220e+00, 3.30700e+00, 3.32980e+00, 3.35050e+00,
	    3.36900e+00, 3.38530e+00, 3.39940e+00, 3.41130e+00, 3.42090e+00, 3.42820e+00, 3.43320e+00, 3.43600e+00,
	    3.43660e+00, 3.43480e+00, 3.43080e+00, 3.42460e+00, 3.41610e+00, 3.40540e+00, 3.39240e+00, 3.37720e+00,
	    3.35970e+00, 3.34000e+00, 3.31800e+00, 3.29370e+00, 3.26700e+00, 3.23810e+00, 3.20680e+00, 3.17310e+00,
	    3.13700e+00, 3.09840e+00, 3.05740e+00, 3.01380e+00, 2.96770e+00, 2.91900e+00, 2.86780e+00, 2.81390e+00,
	    2.75740e+00, 2.69830e+00, 2.63660e+00, 2.57230e+00, 2.50550e+00, 2.43620e+00, 2.36460e+00, 2.29070e+00,
	    2.21470e+00, 2.13680e+00, 2.05730e+00, 1.97630e+00, 1.89420e+00, 1.81120e+00, 1.72780e+00, 1.64440e+00'
	])
	
	const BE_REFR_IDX_IMAG = StdVector(Float64[
	    1.60420e+02, 1.58200e+02, 1.55990e+02, 1.53780e+02, 1.51590e+02, 1.49390e+02, 1.47210e+02, 1.45030e+02,
	    1.42850e+02, 1.40680e+02, 1.38510e+02, 1.36350e+02, 1.34190e+02, 1.32040e+02, 1.29890e+02, 1.27750e+02,
	    1.25610e+02, 1.23480e+02, 1.21360e+02, 1.19240e+02, 1.17130e+02, 1.15030e+02, 1.12940e+02, 1.10860e+02,
	    1.08790e+02, 1.06730e+02, 1.04680e+02, 1.02650e+02, 1.00630e+02, 9.86220e+01, 9.66310e+01, 9.46580e+01,
	    9.27010e+01, 9.07630e+01, 8.88440e+01, 8.69440e+01, 8.50660e+01, 8.32080e+01, 8.13730e+01, 7.95610e+01,
	    7.77720e+01, 7.60070e+01, 7.42670e+01, 7.25520e+01, 7.08630e+01, 6.91990e+01, 6.75620e+01, 6.59520e+01,
	    6.43680e+01, 6.28120e+01, 6.12830e+01, 5.97820e+01, 5.83080e+01, 5.68620e+01, 5.54440e+01, 5.40530e+01,
	    5.26900e+01, 5.13550e+01, 5.00460e+01, 4.87650e+01, 4.75110e+01, 4.62830e+01, 4.50830e+01, 4.39080e+01,
	    4.27590e+01, 4.16360e+01, 4.05390e+01, 3.94660e+01, 3.84180e+01, 3.73940e+01, 3.63950e+01, 3.54180e+01,
	    3.44650e+01, 3.35350e+01, 3.26270e+01, 3.17400e+01, 3.08760e+01, 3.00320e+01, 2.92090e+01, 2.84060e+01,
	    2.76240e+01, 2.68600e+01, 2.61160e+01, 2.53900e+01, 2.46820e+01, 2.39920e+01, 2.33200e+01, 2.26650e+01,
	    2.20260e+01, 2.14030e+01, 2.07960e+01, 2.02050e+01, 1.96290e+01, 1.90680e+01, 1.85210e+01, 1.79880e+01,
	    1.74680e+01, 1.69630e+01, 1.64700e+01, 1.59900e+01, 1.55220e+01, 1.50670e+01, 1.46230e+01, 1.41910e+01,
	    1.37700e+01, 1.33610e+01, 1.29620e+01, 1.25730e+01, 1.21950e+01, 1.18270e+01, 1.14690e+01, 1.11200e+01,
	    1.07810e+01, 1.04510e+01, 1.01290e+01, 9.81720e+00, 9.51360e+00, 9.21860e+00, 8.93200e+00, 8.65370e+00,
	    8.38350e+00, 8.12130e+00, 7.86710e+00, 7.62080e+00, 7.38230e+00, 7.15150e+00, 6.92830e+00, 6.71280e+00,
	    6.50490e+00, 6.30450e+00, 6.11160e+00, 5.92620e+00, 5.74820e+00, 5.57770e+00, 5.41450e+00, 5.25860e+00,
	    5.10990e+00, 4.96830e+00, 4.83380e+00, 4.70620e+00, 4.58540e+00, 4.47120e+00, 4.36350e+00, 4.26200e+00,
	    4.16660e+00, 4.07700e+00, 3.99310e+00, 3.91450e+00, 3.84120e+00, 3.77280e+00, 3.70910e+00, 3.64990e+00,
	    3.59510e+00, 3.54420e+00, 3.49730e+00, 3.45400e+00, 3.41420e+00, 3.37770e+00, 3.34430e+00, 3.31390e+00,
	    3.28630e+00, 3.26140e+00, 3.23900e+00, 3.21910e+00, 3.20130e+00, 3.18580e+00, 3.17230e+00, 3.16060e+00,
	    3.15080e+00, 3.14270e+00, 3.13620e+00, 3.13120e+00, 3.12750e+00, 3.12500e+00, 3.12380e+00, 3.12350e+00,
	    3.12410e+00, 3.12550e+00, 3.12750e+00, 3.13000e+00, 3.13280e+00, 3.13580e+00, 3.13880e+00, 3.14150e+00,
	    3.14390e+00, 3.14560e+00, 3.14660e+00, 3.14640e+00, 3.14500e+00, 3.14200e+00, 3.13720e+00, 3.13020e+00,
	    3.12100e+00, 3.10910e+00, 3.09440e+00, 3.07660e+00, 3.05540e+00, 3.03080e+00, 3.00240e+00, 2.97020e+00
	])
	
	nothing

end

# ╔═╡ b5450e23-7848-4f4e-a278-14c83cdd78c7
"""
Define the materials needed for the application
- Modifies the ScintDetector struct adding the needed materials

"""
function defineMaterials!(det::ScintDetector)
    # to avoid warnings of redefinitions
    G4Material!GetMaterial("vacumm", false) != C_NULL && return

    nist = G4NistManager!Instance()

    #---vacuum-------------------------------------------------------------- 
    vacMat = G4Material("vacumm", z=1.0, a=1.01g/mole, 
		                density=universe_mean_density, state=kStateGas, 
                        temperature=VACUUM_TEMPERATURE, pressure=VACUUM_PRESSURE)
	
    vacPt = G4MaterialPropertiesTable()
    AddProperty(vacPt, "RINDEX", CEBR3_REFR_IDX_ENERGIES, 
		        StdVector(fill(1., length(CEBR3_REFR_IDX_ENERGIES))), useSpline)
                SetMaterialPropertiesTable(vacMat, move!(vacPt))
    det.vacMat = CxxPtr(vacMat)

    #---cebr3----------------------------------------------------------------
	# Define material 
    ce = FindOrBuildElement(nist, "Ce")
    br = FindOrBuildElement(nist, "Br")
    cebr3 = G4Material("cebr3", density=CEBR3_DENSITY, ncomponents=2, 
		               state=kStateSolid, temperature=SATELLITE_TEMP)
    AddElement(cebr3, br, BR_MASS_FRAC)
    AddElement(cebr3, ce, CE_MASS_FRAC)

    #  optical parameters
    opp = G4OpticalParameters!Instance()
    SetScintFiniteRiseTime(opp, false)

    scintPt = G4MaterialPropertiesTable()

    AddConstProperty(scintPt, "SCINTILLATIONYIELD", CEBR3_SCINT_YIELD)
    AddProperty(scintPt, "RINDEX", 
		        CEBR3_REFR_IDX_ENERGIES, 
		        CEBR3_REFR_IDXS, 
		        useSpline)
	
    AddProperty(scintPt, "SCINTILLATIONCOMPONENT1", 
		        CEBR3_SCINT_OPTICAL_ENERGIES, 
		        CEBR3_SCINT_OPTICAL_INTENSITIES, 
		        useSpline)
	
    AddConstProperty(scintPt, "SCINTILLATIONTIMECONSTANT1", 
		             CEBR3_DECAY_TIME_CONSTANT)

    # # of photons emitted = RESOLUTION_SCALE * sqrt(mean # of photons)
    AddConstProperty(scintPt, "RESOLUTIONSCALE", CEBR3_SCINT_RESLN_SCALE)
	
    AddProperty(scintPt, "ABSLENGTH", 
		        CEBR3_ABS_LEN_ENERGIES, 
		        CEBR3_ABS_LEN, 
		        useSpline)
	
    # skip optical Rayleigh scattering (not important)
    # skip Mie scattering (doesn't apply)
    SetMaterialPropertiesTable(cebr3, move!(scintPt))
	
    det.cebr3Mat = CxxPtr(cebr3)

    #---quartz-----------------------------------------------------------------
    qz = FindOrBuildMaterial(nist, "G4_SILICON_DIOXIDE")
    qzPt = G4MaterialPropertiesTable()
    AddProperty(qzPt, "RINDEX", QZ_REFR_IDX_ENERGIES, QZ_REFR_IDXS )
    # do not add real/imag refractive indices to materials. they areo nly for surface indices.
    # qzPt->AddProperty(kREFR_IDX_REAL, QZ_REFR_IDX_ENERGIES, QZ_REFR_IDXS);
    # qzPt->AddProperty(kREFR_IDX_IMAG, QZ_REFR_IDX_ENERGIES, std::vector<G4double>(QZ_REFR_IDXS.size(), 0.0));
    SetMaterialPropertiesTable(qz, move!(qzPt))
    det.qzMat = qz

    #---teflon-------------------------------------------------------------------
    teflon = FindOrBuildMaterial(nist, "G4_TEFLON")
    tefPt = G4MaterialPropertiesTable()
    AddProperty(tefPt, "RINDEX", 
		        TEFLON_REFR_IDX_ENERGIES, 
		        TEFLON_REFR_IDXS, 
		        useSpline)
    AddProperty(tefPt, "REFLECTIVITY", 
		        TEFLON_REFR_IDX_ENERGIES, 
		        TEFLON_REFLECTIVITY, 
		        useSpline)
    SetMaterialPropertiesTable(teflon, move!(tefPt))
    det.tefMat = CxxPtr(teflon)

    #---Silicon-------------------------------------------------------------------
    si = FindOrBuildMaterial(nist, "G4_Si")
    simPt = G4MaterialPropertiesTable()
    # set to one and apply detection in post-processing (?)
    AddProperty(simPt, "EFFICIENCY", 
		        SI_DET_EFF_ENERGIES, 
		        SI_DET_EFF, 
		        useSpline)
	
    refl = StdVector(fill(0., length(SI_BROADCOM_NUMBERS)))
	
    AddProperty(simPt, "REFLECTIVITY", 
		        SI_DET_EFF_ENERGIES, 
		        refl, 
		        useSpline)
	
    AddProperty(simPt, "TRANSMITTANCE", 
		        SI_DET_EFF_ENERGIES, 
		        SI_TRANSMITTANCE, useSpline)
	
    AddProperty(simPt, "RINDEX", 
		        SI_REFR_IDX_ENERGY, 
		        SI_REFR_IDX_REAL, 
		        useSpline)
	
    AddProperty(simPt, 
		        "REALRINDEX", 
		        SI_REFR_IDX_ENERGY, 
		        SI_REFR_IDX_REAL, 
		        useSpline)
    AddProperty(simPt, "IMAGINARYRINDEX", 
		        SI_REFR_IDX_ENERGY, 
		        SI_REFR_IDX_IMAG, 
		        useSpline)
    
	SetMaterialPropertiesTable(si, move!(simPt))
    det.siMat = si

    #---aluminum--------------------------------------------------------------------
    alElt = FindOrBuildElement(nist, "Al")
    alElt == C_NULL &&  error("Nist manager gave nullptr for aluminum")

    al = G4Material("aluminum", density=2.712g/cm3, ncomponents=1, 
		            state=kStateSolid, temperature=SATELLITE_TEMP)
    AddElement(al, alElt, 1.0)
    alPt = G4MaterialPropertiesTable()
	
    # DO NOT ADD JUST RINDEX!!!!!!!!!!!!!!!!!!!1
    # IT WILL CAUSE PROGRAM TO CRASH IF YOU PUT A SKIN AROUND THE MATERIAL
    # AddProperty(alPt, "RINDEX", AL_REFR_IDX_ENERGIES, AL_REFR_IDX_REAL)
    
	AddProperty(alPt, "REALRINDEX", AL_REFR_IDX_ENERGIES, AL_REFR_IDX_REAL)
    AddProperty(alPt, "IMAGINARYRINDEX", AL_REFR_IDX_ENERGIES, AL_REFR_IDX_IMAG)
    SetMaterialPropertiesTable(al, move!(alPt))
    det.alMat = CxxPtr(al)

    #---beryllium-------------------------------------------------------------------
    beElt = FindOrBuildElement(nist, "Be")
    be = G4Material("beryllium", density=BE_DENSITY, ncomponents=1, 
		            state=kStateSolid, temperature=SATELLITE_TEMP)
    AddElement(be, beElt, 1)
    bePt = G4MaterialPropertiesTable()
    AddProperty(bePt, "RINDEX", BE_REFR_IDX_ENERGIES, BE_REFR_IDX_REAL)
    AddProperty(bePt, "REALRINDEX", BE_REFR_IDX_ENERGIES, BE_REFR_IDX_REAL)
    AddProperty(bePt, "IMAGINARYRINDEX", BE_REFR_IDX_ENERGIES, BE_REFR_IDX_IMAG)
    SetMaterialPropertiesTable(be, move!(bePt))
    det.beMat = CxxPtr(be)
    
    return
end


# ╔═╡ 2c1dd74c-a6f3-4c11-a7e5-daf0df670d79
"""
Construct the detector geometry
"""
function scintConstruct(det::ScintDetector)::CxxPtr{G4VPhysicalVolume}
    
	(; worldSize, doTeflon, doAluminum, crys_size, si_side, si_thick, tef_thick, al_thick, air_gap, checkOverlaps) = det
 
    #---Materials---------------------------------------------------------------
    defineMaterials!(det)

    #---World-------------------------------------------------------------------
    worldBox = G4Box("World", 0.5*worldSize, 0.5*worldSize, 0.5*worldSize)        
    worldLogVol = G4LogicalVolume(worldBox, det.vacMat, "World")
    worldVis = G4VisAttributes(true, G4Colour(1, 1, 1, 0.05))
    SetVisAttributes(worldLogVol, worldVis)
    worldPlacement = G4PVPlacement(nothing, G4ThreeVector(), worldLogVol, 
		                           "World", nothing, false, 0, 
		                           checkOverlaps)

    #---cryst--------------------------------------------------------------------
    crystBox = G4Box("cebr3_box", crys_size/2, crys_size/2, crys_size/2)
    crystLogVol = G4LogicalVolume(crystBox, det.cebr3Mat, "cebr3_log")
    crystVis = G4VisAttributes(true, G4Colour(0.35, 0.5, 0.92, 0.8))
    SetVisAttributes(crystLogVol, crystVis)
    crystPlacement = G4PVPlacement(nothing, G4ThreeVector(), crystLogVol, 
		                           "cebr3_phys", worldLogVol, false, 0, 
		                           checkOverlaps)

    surf = G4OpticalSurface("cebr3_optical_edge")
    # See G4 documentation on UNIFIED model
    SetType(surf, dielectric_LUTDAVIS)
    SetModel(surf, DAVIS)
    SetFinish(surf, Rough_LUT);
    G4LogicalBorderSurface("cebr3_logical_border_surf", 
		                   CxxPtr(crystPlacement), 
		                   CxxPtr(worldPlacement), 
		                   CxxPtr(surf))

    #---teflon--------------------------------------------------------------------
    if doTeflon
        precutHalfThick = (air_gap + tef_thick + crys_size)/2
        preCutPtfeBox = G4Box("precut_ptfe", 
			                   precutHalfThick, 
			                   precutHalfThick, 
			                   precutHalfThick)
        addAirGapHalfXyz = 0.5 * (crys_size + air_gap)

        crystPlusAirGap = G4Box("", 
			                    addAirGapHalfXyz, 
			                    addAirGapHalfXyz, 
			                    addAirGapHalfXyz)
		
        translateSubtract = G4ThreeVector(0, 0, 0)
		
        slicedPtfe = G4SubtractionSolid("sliced_ptfe", 
			                            move!(preCutPtfeBox), 
			                            move!(crystPlusAirGap), 
			                            CxxPtr{G4RotationMatrix}(C_NULL), 
			                            translateSubtract)
		
        translateSubtract = G4ThreeVector(0, 0, 1.8*cm)
		
        cutCap = G4Box("", 
			           precutHalfThick + 5*mm, 
			           precutHalfThick + 5*mm, 
			           precutHalfThick)
		
        slicedPtfeOpen = G4SubtractionSolid("sliced_ptfe_open", 
			                                move!(slicedPtfe), 
			                                move!(cutCap), 
			                                CxxPtr{G4RotationMatrix}(C_NULL),
			                                translateSubtract)
		
        ptfeLogVol = G4LogicalVolume(slicedPtfeOpen, det.tefMat, "ptfe_log")
        tefVis = G4VisAttributes(true, G4Colour(0, 1, 0, 0.1))
        SetVisAttributes(ptfeLogVol, tefVis)

        # attachTeflonOpticalSurf
        tefSuf = G4OpticalSurface("teflon_optical_surf")
        SetType(tefSuf, dielectric_dielectric)
        SetModel(tefSuf, unified)
        SetFinish(tefSuf, groundfrontpainted)
        SetSigmaAlpha(tefSuf, 0.)
        SetMaterialPropertiesTable(tefSuf, GetMaterialPropertiesTable(det.tefMat))
        G4LogicalSkinSurface("teflon_skin_surf", 
			                 CxxPtr(ptfeLogVol), 
			                 CxxPtr(tefSuf))
        
        G4PVPlacement(nothing, G4ThreeVector(), ptfeLogVol, 
			          "ptfe_phys", worldLogVol, false, 0, checkOverlaps)
    end

    #---aluminum-------------------------------------------------------------------
    if doAluminum
        sliceCapThick = 0.5 * (tef_thick + crys_size + al_thick + 2*air_gap) + 5*mm
        precutHalfThick = sliceCapThick - 5*mm
        sliceOutHalfThick = precutHalfThick - 0.5 * al_thick
      
        precut = G4Box("", 
			           precutHalfThick, 
			           precutHalfThick, 
			           precutHalfThick)
		
        cutThisOut = G4Box("", 
			               sliceOutHalfThick, 
			               sliceOutHalfThick, 
			               sliceOutHalfThick)
		
        cutOffEnd = G4Box("", sliceCapThick, sliceCapThick, 6*mm)
		
        hollow = G4SubtractionSolid("hollow_al", 
			                        move!(precut), 
			                        move!(cutThisOut), 
			                        CxxPtr{G4RotationMatrix}(C_NULL), 
			                        G4ThreeVector())
		
        translateSubtract = G4ThreeVector(0, 0, sliceCapThick)
        open = G4SubtractionSolid("hollow_open_al", 
			                      move!(hollow), 
			                      move!(cutOffEnd), 
			                      CxxPtr{G4RotationMatrix}(C_NULL), 
			                      translateSubtract)
        
        alLogVol = G4LogicalVolume(open, det.alMat, "al_log")
        alVis = G4VisAttributes(true, G4Colour(1, 0., 0., 0.1))
        SetVisAttributes(alLogVol, alVis)

        # attachAlOpticalSurface(alLogVol)
        alSurf = G4OpticalSurface("al_optical_surf")
        SetModel(alSurf, unified)
        SetType(alSurf, dielectric_metal)
        SetFinish(alSurf, polished)
        SetMaterialPropertiesTable(alSurf, GetMaterialPropertiesTable(det.alMat))
        G4LogicalSkinSurface("al_skin_surf", CxxPtr(alLogVol), CxxPtr(alSurf))
        G4PVPlacement(nothing, G4ThreeVector(), alLogVol, 
			          "al_phys", worldLogVol, false, 0, checkOverlaps)
    end

    #---sillicon----------------------------------------------------------------
    siBox = G4Box("si_box", si_side/2, si_side/2, si_thick/2)
    siLogVol = G4LogicalVolume(siBox, det.siMat, "si_log")
    siVis = G4VisAttributes(true, G4Colour(0, 0, 1, 0.5))
    SetVisAttributes(siLogVol, siVis)
  
    translate = G4ThreeVector(0, 0, crys_size/2 + si_thick/2)
    G4PVPlacement(nothing, translate, siLogVol, "si_phys", 
		          worldLogVol, false, 0, checkOverlaps)

    siSurf = G4OpticalSurface("si_surf")
    SetModel(siSurf, unified)
    SetFinish(siSurf, polished)
    SetType(siSurf, dielectric_dielectric)
    SetMaterialPropertiesTable(siSurf, GetMaterialPropertiesTable(det.siMat))
    G4LogicalSkinSurface("si_skin", CxxPtr(siLogVol), CxxPtr(siSurf))
  
    #  return the world
    return worldPlacement
end



# ╔═╡ 5bd4dc5f-28c8-4bf5-b521-e2946e938483
Geant4.getConstructor(::ScintDetector)::Function = scintConstruct

# ╔═╡ Cell order:
# ╠═05b9f480-5390-11ee-2101-4b947107b1e9
# ╠═ed6d6e22-93a7-4825-a018-e26b30f08539
# ╠═fa607f1f-b43b-4e12-bb4b-ebf18cc651c9
# ╠═4c346534-bded-48ed-91d7-17dea3bd6933
# ╠═57cc79b3-7a00-4cd1-a935-3c7380870b8c
# ╠═e3411132-8907-4c6d-85e6-02460a31f331
# ╠═5ffed418-69ab-40ee-be61-5582f20df6e8
# ╠═0fdc6a21-97bc-47e1-bf25-9f7f27defa85
# ╠═b73b55f8-883d-466a-bfa5-c2ff42cb1167
# ╠═4273bde0-b5fe-4e7c-8c2a-65b38a96ca2e
# ╠═3e47feeb-9ac4-4e19-86ad-42c944614a84
# ╠═9555f147-771e-4aa5-8fe7-c784a290cebd
# ╠═5a09e358-5955-4628-86e5-44f571933876
# ╠═973eac13-e50d-4da8-b284-31475d9853c5
# ╠═5bd4dc5f-28c8-4bf5-b521-e2946e938483
# ╠═3953460d-d0ae-4bb6-95d0-4d57daf8100d
# ╠═dce6c5ed-66f2-447f-ba12-5155cb460948
# ╠═b92515e3-fc13-4776-9d2c-e2a8c843a24f
# ╠═2853e73b-3e20-4c53-ab0a-1e467bb5cbd6
# ╠═06e633da-dfee-4a84-b11e-1ce376c076eb
# ╠═fae57a76-f596-4fd6-b447-818e9bb6c446
# ╠═874dc33d-bb41-43f2-8b43-b25811eb08d1
# ╠═ef5d20a8-d572-4309-9380-e8b17b9029ae
# ╠═c84b2410-e33a-4d69-a96c-813bf931a523
# ╠═920e3820-b4c2-4132-9b93-6625514aec87
# ╠═bf408fe1-9995-4699-84b6-2868c16e69e8
# ╠═7787a486-eb86-450a-a98e-b332d39dcbda
# ╠═a2a1eddc-1d75-4026-b69f-09093b3baf63
# ╠═670a8f07-cd3d-44c8-aa76-577814f02495
# ╠═ab3b41d2-e043-4afa-bc44-55694f6cc403
# ╠═1371985e-2392-4500-bc67-d9ed60ac6117
# ╠═40af1a84-df87-461c-a656-25b986641328
# ╠═b5450e23-7848-4f4e-a278-14c83cdd78c7
# ╠═2c1dd74c-a6f3-4c11-a7e5-daf0df670d79
# ╠═fb8ce755-a055-468f-bf07-214d852b65d6
# ╠═433d8de4-cf0e-4e3c-aa4c-99fcdd6562b6
# ╠═2c5ade51-5d9c-4a27-ba40-b23a24ddd74a
# ╠═8f27a4aa-634f-4185-a2be-72d9e025eb20
# ╠═2cb2f53b-adda-4b6c-9059-12c3b36a9a60
# ╠═736418a2-2a6d-42f6-95fd-eb3a41a4b372
# ╠═a8d80ee8-15c3-47a4-b94a-2e3d8d013b16
# ╠═c289d0a6-50b4-4645-9d8a-61fd783e6491
# ╠═ef5ec36c-e282-485a-8fc3-433f4626f35a
# ╠═52afa901-e0d3-4299-8011-cf96d2155b0a
# ╠═567fca33-3bf8-472c-b422-cd1e0a08562f
# ╠═d6b7ef9d-8d26-439c-af4f-c1a5abe8633a
# ╠═e5531168-e245-4d32-9e47-fc26b6d78aab
# ╠═cbc8085c-e0ff-471f-9fed-a86eb88efa55
# ╠═7d31ec0c-6bf5-4da3-97db-4419da2d12ba
# ╠═3b63f6c8-cf82-4e3d-9714-0b9e7fb9e32b
# ╠═59c1a1e1-f6fc-481a-b7e5-b9fab665f9f3
# ╠═b84f3ccc-6cb8-40fa-9d6d-d93090d1e620
# ╠═c27b3002-4711-404e-9b2c-7c76d66364a5
# ╠═ba9c2dde-5e31-4778-8c14-60d6d1669d30
# ╠═8a77e357-cd2d-4183-b0ab-4fef8e6282eb
# ╠═34f1d993-450e-4cea-803a-b733fff9dcdc
# ╠═857cdecb-27a2-4dec-a02d-d6deef70c869
# ╠═483f9846-3568-4817-86de-36ee556666d6
# ╠═54b88049-5f28-4029-8897-a9aaed82f954
# ╠═1b55751d-4aff-41f6-b022-6b7f1e874f5f
# ╠═79d7aab9-6f9f-44b1-9fd9-f57617a17672
# ╠═aa6c726f-f36b-4cd5-a892-f36b160bafe8
# ╠═0d98b72e-2eb2-44c5-a474-75a04c8656bf
# ╠═dbc332f0-ec5b-4747-8dd1-5ec0c2710bff
