using GLMakie
using Oceananigans
using Oceananigans.Units
using Printf

# Parameters
Δz = 4          # Vertical resolution
Lz = 128        # Extent of vertical domain
Nz = Int(Lz/Δz) # Vertical resolution
f₀ = 0.0        # Coriolis parameter (s⁻¹)
N² = 1e-5       # Buoyancy gradient (s⁻²)
arch = CPU()

# Jᵇ = α * g * Q / (ρ₀ cₚ) where Q is heat flux, α is the thermal expansion coefficient
Jᵇ = 0.0        # Surface buoyancy flux (m² s⁻³)
τˣ = 1e-4       # Surface kinematic momentum flux (m s⁻¹)
stop_time = 2days

grid = RectilinearGrid(arch,
                       size = (Nx, Ny, Nz),
                       x = (0, 2Lz),
                       y = (0, 2Lz),
                       z = (-Lz, 0),
                       topology = (Periodic, Periodic, Bounded))

# Coriolis force
coriolis = FPlane(f=f₀)

# Boundary conditions
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵇ))
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τˣ))

# Turbulence closure
#closure = convective_adjustment
closure = CATKEVerticalDiffusivity()

model = NonhydrostaticModel(; grid, coriolis,
                            advection = WENO(order=9),
                            timestepper = :RungeKutta3,
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            boundary_conditions = (; b=b_bcs, u=u_bcs))
                                    
bᵢ(z) = N² * z
set!(model, b=bᵢ)

simulation = Simulation(model; Δt=10.0, stop_time)
conjure_time_step_wizard!(simulation, IterationInterval(10), cfl=0.5)

outputs = merge(model.velocities, model.tracers)

filename = "wind_driven_deepening.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, outputs; filename,
                                                      schedule = TimeInterval(10minutes),
                                                      overwrite_existing = true)

progress(sim) = @info string("Iter: ", iteration(sim), " t: ", prettytime(sim),
                             ", max(b): ", maximum(model.tracers.b))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

@info "Running a simulation of $model..."

run!(simulation)

