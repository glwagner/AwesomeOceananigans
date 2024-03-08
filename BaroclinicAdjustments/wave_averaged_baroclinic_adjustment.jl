using Oceananigans
using Oceananigans.Units
using Printf

arch = GPU()

# Domain parameters
Lx = 10kilometers # east-west extent [m]
Ly = 10kilometers # north-south extent [m]
Lz = 160           # depth [m]
f₀ = 1e-4          # [s⁻¹] Coriolis parameter 
Nx = 512
Ny = 512
Nz = 128

# Front / initial condition parameters
hᵢ = 50            # Initial mixed layer depth (m)
N² = 1e-8          # [s⁻²] sub-mixed-layer buoyancy frequency / stratification
M² = 2e-8          # [s⁻²] horizontal buoyancy gradient
Δy = 1kilometer    # [m] width of the region of the front
Δb = Δy * M²       # buoyancy jump associated with the front
ϵb = 1e-2 * Δb     # noise amplitude

# Surface forcing parameters
u★ = 5e-3  # [m s⁻¹] ocean-side friction velocity
θ★ = 30    # [ᵒ] wind direction (relative to along-front direction)
Jᵇ = 1e-9  # [m² s⁻³] surface buoyancy flux 

# Gravity wave parameters
a = 0.8          # [m] wave amplitude
λ = 60           # [m] wavelength
g = 9.81         # [m s⁻²] gravitational acceleration
k = 2π / λ       # [m⁻¹] wavenumber
σ = sqrt(g * k)  # [s⁻¹] wave frequency
Uˢ = a^2 * k * σ # [m s⁻¹] surface Stokes drift
θˢ = 30          # [ᵒ] Stokes drift direction (relative to along-front direction)

# Output parameters
save_fields_interval = 0.5day
output_directory = "nobackup/users/glwagner/awesome-oceananigans" #"."

grid = RectilinearGrid(arch,
                       size = (Nx, Ny, Nz),
                       halo = (4, 4, 4),
                       x = (-Lx/2, Lx/2),
                       y = (-Ly/2, Ly/2),
                       z = (-Lz, 0),
                       topology = (Periodic, Bounded, Bounded))

# Possibly add bathymetry
# δx = δy = 4kilometer
# h = 50 # meters
# @inline bump(x, y) = -Lz + h * exp(-x^2 / 2δx^2 - y^2 / 2δy^2)
# grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bump))

τˣ = -u★^2 * cosd(θ★)
τʸ = -u★^2 * sind(θ★)

u_top_bc = FluxBoundaryCondition(τˣ)
v_top_bc = FluxBoundaryCondition(τʸ)
u_bcs = FieldBoundaryConditions(top=u_top_bc)
v_bcs = FieldBoundaryConditions(top=v_top_bc)

b_top_bc = FluxBoundaryCondition(Jᵇ)
b_bcs = FieldBoundaryConditions(top=b_top_bc)

return 2k * Uˢ * exp(2k * z) * cos_θˢ

@inline function ∂z_Uˢ(z, t, parameters)
    k = parameters.k
    Uˢ = parameters.Uˢ
    return 2k * Uˢ * exp(2k * z)
end

@inline ∂z_uˢ(z, t, p) = p.cos_θˢ * ∂z_Uˢ(z, t, p)
@inline ∂z_vˢ(z, t, p) = p.sin_θˢ * ∂z_Uˢ(z, t, p)

stokes_drift_parameters = (; k, Uˢ, cos_θˢ = cosd(θˢ), sin_θˢ = sind(θˢ))
stokes_drift = UniformStokesDrift(; ∂z_uˢ, ∂z_vˢ, parameters=stokes_drift_parameters)

model = NonhydrostaticModel(; grid, stokes_drift,
                            boundary_conditions = (; u=u_bcs, b=b_bcs),
                            timestepper = :RungeKutta3,
                            coriolis = FPlane(f=f₀),
                            buoyancy = BuoyancyTracer(),
                            tracers = :b,
                            advection = WENO())

ramp(y, Δy) = min(max(0, y/Δy + 1/2), 1)

function bᵢ(x, y, z)
    z̃ = min(z, -hᵢ)
    return N² * z̃ + Δb * ramp(y, Δy) + ϵb * randn()
end

set!(model, b=bᵢ)

simulation = Simulation(model, Δt=1minute, stop_time=20days)
conjure_time_step_wizard!(simulation, IterationInterval(20), cfl=0.5, max_Δt=20minutes)
b = model.tracers.b
u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)

wall_clock = Ref(time_ns())

function print_progress(sim)
    u, v, w = model.velocities
    progress = 100 * (time(sim) / sim.stop_time)
    elapsed = (time_ns() - wall_clock[]) / 1e9

    @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s\n",
            progress, iteration(sim), prettytime(sim), prettytime(elapsed),
            maximum(abs, interior(u)),
            maximum(abs, interior(v)),
            maximum(abs, interior(w)),
            prettytime(sim.Δt))

    wall_clock[] = time_ns()
    
    return nothing
end

add_callback!(simulation, print_progress, IterationInterval(100))

filename = "wave_averaged_baroclinic_adjustment"

slicers = (east = (grid.Nx, :, :),
           north = (:, grid.Ny, :),
           bottom = (:, :, 1),
           top = (:, :, grid.Nz))

for side in keys(slicers)
    indices = slicers[side]

    simulation.output_writers[side] = JLD2OutputWriter(model, (; b, ζ);
                                                       filename = filename * "_$(side)_slice",
                                                       dir = output_directory,
                                                       schedule = TimeInterval(save_fields_interval),
                                                       overwrite_existing = true,
                                                       indices)
end

run!(simulation)

