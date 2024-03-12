# # Submesoscale baroclinic adjustment
#
# Roughly follows
#
# Hamlington et al.,
# "Langmuir-Submesoscale Interactions: Descriptive Analysis of Multiscale Frontal Spindown Simulations",
# Journal of Physical Oceananography, 2014

using Oceananigans
using Oceananigans.Units

arch = GPU()
Nx = 2048
Ny = 2048
Nz = 96
save_fields_interval = 1hour
output_directory = "." #"nobackup/users/glwagner/awesome-oceananigans" #"."
filename = "submesoscale_baroclinic_adjustment_large_eddy_simulation"

include("submesoscale_baroclinic_adjustment_setup.jl")

grid = RectilinearGrid(arch,
                       size = (Nx, Ny, Nz),
                       halo = (4, 4, 4),
                       x = (-Lx/2, Lx/2),
                       y = (-Ly/2, Ly/2),
                       z = (-Lz, 0),
                       topology = (Periodic, Bounded, Bounded))

model = NonhydrostaticModel(; grid, stokes_drift,
                            boundary_conditions = (; u=u_bcs, b=b_bcs),
                            timestepper = :RungeKutta3,
                            coriolis = FPlane(f=f₀),
                            buoyancy = BuoyancyTracer(),
                            tracers = :b,
                            advection = WENO(order=9))

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

slicers = (east = (grid.Nx, :, :),
           north = (:, grid.Ny, :),
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

