using Oceananigans
using Oceananigans.Units
using Printf

arch = CPU()
Lx = 100kilometers # east-west extent [m]
Ly = 100kilometers # north-south extent [m]
Lz = 400           # depth [m]
δx = 5kilometers   # Width of the bump in x
δy = 20kilometers  # Width of the bump in y
h = 100            # bump height (m)
N² = 1e-5 # [s⁻²] buoyancy frequency / stratification
M² = 1e-6 # [s⁻²] horizontal buoyancy gradient
Δy = 10kilometers  # width of the region of the front
Δb = Δy * M²       # buoyancy jump associated with the front
ϵb = 1e-2 * Δb     # noise amplitude
save_fields_interval = 0.5day

grid = RectilinearGrid(arch,
                       size = (48, 48, 24),
                       halo = (4, 4, 4),
                       x = (-Lx/2, Lx/2),
                       y = (-Ly/2, Ly/2),
                       z = (-Lz, 0),
                       topology = (Periodic, Bounded, Bounded))

@inline bump(x, y) = -Lz + h * exp(-x^2 / 2δx^2 - y^2 / 2δy^2)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bump))

model = NonhydrostaticModel(; grid,
                            timestepper = :RungeKutta3,
                            coriolis = BetaPlane(latitude = -45),
                            buoyancy = BuoyancyTracer(),
                            tracers = :b,
                            advection = WENO())

ramp(y, Δy) = min(max(0, y/Δy + 1/2), 1)
bᵢ(x, y, z) = N² * z + Δb * ramp(y, Δy) + ϵb * randn()
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
            maximum(abs, u), maximum(abs, v), maximum(abs, w), prettytime(sim.Δt))

    wall_clock[] = time_ns()
    
    return nothing
end

add_callback!(simulation, print_progress, IterationInterval(100))

filename = "bumpy_baroclinic_adjustment"

slicers = (east = (grid.Nx, :, :),
           north = (:, grid.Ny, :),
           bottom = (:, :, 1),
           top = (:, :, grid.Nz))

for side in keys(slicers)
    indices = slicers[side]

    simulation.output_writers[side] = JLD2OutputWriter(model, (; b, ζ);
                                                       filename = filename * "_$(side)_slice",
                                                       schedule = TimeInterval(save_fields_interval),
                                                       overwrite_existing = true,
                                                       indices)
end

run!(simulation)

#=
using GLMakie

sides = keys(slicers)
slice_filenames = NamedTuple(side => filename * "_$(side)_slice.jld2" for side in sides)

b_timeserieses = (east   = FieldTimeSeries(slice_filenames.east, "b"),
                  north  = FieldTimeSeries(slice_filenames.north, "b"),
                  bottom = FieldTimeSeries(slice_filenames.bottom, "b"),
                  top    = FieldTimeSeries(slice_filenames.top, "b"))

xb, yb, zb = nodes(b_timeserieses.east)

xb = xb ./ 1e3 # convert m -> km
yb = yb ./ 1e3 # convert m -> km

Nx, Ny, Nz = size(grid)
x_xz = repeat(x, 1, Nz)
y_xz_north = y[end] * ones(Nx, Nz)
z_xz = repeat(reshape(z, 1, Nz), Nx, 1)

x_yz_east = x[end] * ones(Ny, Nz)
y_yz = repeat(y, 1, Nz)
z_yz = repeat(reshape(z, 1, Nz), grid.Ny, 1)

x_xy = x
y_xy = y
z_xy_top = z[end] * ones(grid.Nx, grid.Ny)
z_xy_bottom = z[1] * ones(grid.Nx, grid.Ny)
nothing #hide

# Then we create a 3D axis. We use `zonal_slice_displacement` to control where the plot of the instantaneous
# zonal average flow is located.

fig = Figure(size = (1600, 800))

zonal_slice_displacement = 1.2

ax = Axis3(fig[2, 1],
           aspect=(1, 1, 1/5),
           xlabel = "x (km)",
           ylabel = "y (km)",
           zlabel = "z (m)",
           xlabeloffset = 100,
           ylabeloffset = 100,
           zlabeloffset = 100, 
           limits = ((x[1], zonal_slice_displacement * x[end]), (y[1], y[end]), (z[1], z[end])),
           elevation = 0.45,
           azimuth = 6.8,
           xspinesvisible = false,
           zgridvisible = false,
           protrusions = 40,
           perspectiveness = 0.7)

n = length(times)

b_slices = (east   = interior(b_timeserieses.east[n], 1, :, :),
            north  = interior(b_timeserieses.north[n], :, 1, :),
            bottom = interior(b_timeserieses.bottom[n], :, :, 1),
            top    = interior(b_timeserieses.top[n], :, :, 1))

clims = 1.1 .* extrema(b_timeserieses.top[n][:])

kwargs = (colorrange=clims, colormap=:deep)
surface!(ax, x_yz_east, y_yz, z_yz;    color = b_slices.east, kwargs...)
surface!(ax, x_xz, y_xz_north, z_xz;   color = b_slices.north, kwargs...)
surface!(ax, x_xy, y_xy, z_xy_bottom ; color = b_slices.bottom, kwargs...)
surface!(ax, x_xy, y_xy, z_xy_top;     color = b_slices.top, kwargs...)

display(fig)
=#
