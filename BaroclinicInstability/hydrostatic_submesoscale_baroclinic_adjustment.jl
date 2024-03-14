# # Hydrostatic submesoscale baroclinic adjustment
#
# The setup roughly follows
#
# Hamlington et al.,
# "Langmuir-Submesoscale Interactions: Descriptive Analysis of Multiscale Frontal Spindown Simulations",
# Journal of Physical Oceananography, 2014

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity, RiBasedVerticalDiffusivity
using Oceananigans.TurbulenceClosures: FivePointHorizontalFilter

arch = GPU()
Nx = Ny = 256 
Nz = 80
weno_order = 5
save_fields_interval = 1hour
output_directory = "." #"nobackup/users/glwagner/awesome-oceananigans" #"."

include("submesoscale_baroclinic_adjustment_setup.jl")

setup = langmuir_density_front_setup() #u★=0.0)
Lx, Ly, Lz = setup.domain

grid = RectilinearGrid(arch,
                       size = (Nx, Ny, Nz),
                       halo = (7, 7, 7),
                       x = (-Lx/2, Lx/2),
                       y = (-Ly/2, Ly/2),
                       z = (-Lz, 0),
                       topology = (Periodic, Bounded, Bounded))

closures = [
    CATKEVerticalDiffusivity()
    RiBasedVerticalDiffusivity()
#    RiBasedVerticalDiffusivity(horizontal_Ri_filter=FivePointHorizontalFilter())
    ConvectiveAdjustmentVerticalDiffusivity(convective_κz=1e-1, convective_νz=1e-2)
]

for closure in closures

    prefix = if closure isa CATKEVerticalDiffusivity
                "catke"
             elseif closure isa RiBasedVerticalDiffusivity
                "ri_based"
             elseif closure isa ConvectiveAdjustmentVerticalDiffusivity
                "convective_adjustment"
             end

    filename = "$(prefix)_modified_hydrostatic_submesoscale_baroclinic_adjustment_Nx$(Nx)_Nz$(Nz)"

    @show closure

    model = HydrostaticFreeSurfaceModel(; grid, closure,
                                        boundary_conditions = setup.boundary_conditions,
                                        coriolis = setup.coriolis,
                                        buoyancy = setup.buoyancy,
                                        tracers = (setup.tracers..., :e),
                                        momentum_advection = WENO(order=weno_order),
                                        tracer_advection = WENO(order=weno_order))

    set!(model; e=1e-6, setup.initial_condition...)

    simulation = Simulation(model, Δt=1minute, stop_time=30days)
    conjure_time_step_wizard!(simulation, IterationInterval(20), cfl=0.3, max_Δt=20minutes)
    b = model.tracers.b
    e = model.tracers.e
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

        simulation.output_writers[side] = JLD2OutputWriter(model, (; b, ζ, e);
                                                           filename = filename * "_$(side)_slice",
                                                           dir = output_directory,
                                                           schedule = TimeInterval(save_fields_interval),
                                                           overwrite_existing = true,
                                                           indices)
    end

    run!(simulation)

end
