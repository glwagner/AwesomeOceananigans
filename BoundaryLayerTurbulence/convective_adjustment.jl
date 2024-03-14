using GLMakie
using Oceananigans
using Oceananigans.Units
using Printf

using Oceananigans.TurbulenceClosures: ConvectiveAdjustmentVerticalDiffusivity
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

# Parameters
Δz = 4          # Vertical resolution
Lz = 256        # Extent of vertical domain
Nz = Int(Lz/Δz) # Vertical resolution
f₀ = 1e-4       # Coriolis parameter (s⁻¹)
N² = 1e-5       # Buoyancy gradient (s⁻²)

# Jᵇ = α * g * Q / (ρ₀ cₚ) where Q is heat flux, α is the thermal expansion coefficient
Jᵇ = +1e-7      # Surface buoyancy flux (m² s⁻³)
τˣ = 0.0        # Surface kinematic momentum flux (m s⁻¹)
stop_time = 2days

κ_ca = 1.0 # m² s⁻¹
convective_adjustment = ConvectiveAdjustmentVerticalDiffusivity(convective_κz=κ_ca, 
                                                                convective_νz=κ_ca)

grid = RectilinearGrid(size=Nz, z=(-Lz, 0), topology=(Flat, Flat, Bounded))

# Coriolis force
coriolis = FPlane(f=f₀)

# Boundary conditions
b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵇ))
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τˣ))

# Turbulence closure
#closure = convective_adjustment
closure = CATKEVerticalDiffusivity()

model = HydrostaticFreeSurfaceModel(; grid, closure, coriolis,
                                    tracers = (:b, :e),
                                    buoyancy = BuoyancyTracer(),
                                    boundary_conditions = (; b=b_bcs, u=u_bcs))
                                    
bᵢ(z) = N² * z
set!(model, b=bᵢ, e=1e-6)

simulation = Simulation(model; Δt=2minutes, stop_time)

diffusivities = (κᵘ = model.diffusivity_fields.κᵘ,
                 κᶜ = model.diffusivity_fields.κᶜ)

outputs = merge(model.velocities, model.tracers, diffusivities)

filename = "convective_adjustment.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, outputs; filename,
                                                      schedule = TimeInterval(20minutes),
                                                      overwrite_existing = true)

progress(sim) = @info string("Iter: ", iteration(sim), " t: ", prettytime(sim),
                             ", max(b): ", maximum(model.tracers.b))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

@info "Running a simulation of $model..."

run!(simulation)

bt = FieldTimeSeries(filename, "b")
ut = FieldTimeSeries(filename, "u")
vt = FieldTimeSeries(filename, "v")
κᶜt = FieldTimeSeries(filename, "κᶜ")
κᵘt = FieldTimeSeries(filename, "κᵘ")

zc = znodes(bt)
zf = znodes(κᶜt)
Nt = length(bt.times)

fig = Figure(size=(1500, 600))

slider = Slider(fig[2, 1:3], range=1:Nt, startvalue=1)
n = slider.value

buoyancy_label = @lift "Buoyancy at t = " * prettytime(bt.times[$n])
velocities_label = @lift "Velocities at t = " * prettytime(bt.times[$n])
diffusivities_label = @lift "Eddy diffusivities at t = " * prettytime(bt.times[$n])

axb = Axis(fig[1, 1], xlabel=buoyancy_label, ylabel="z (m)")
axu = Axis(fig[1, 2], xlabel=velocities_label, ylabel="z (m)")
axκ = Axis(fig[1, 3], xlabel=diffusivities_label, ylabel="z (m)")

xlims!(axb, -grid.Lz * N², 0)
xlims!(axu, -0.1, 0.1)
xlims!(axκ, -1e-1, 5e-1)

colors = [:black, :blue, :red, :orange]

bn  = @lift interior(bt[$n], 1, 1, :)
un  = @lift interior(ut[$n], 1, 1, :)
vn  = @lift interior(vt[$n], 1, 1, :)
κᶜn = @lift interior(κᶜt[$n], 1, 1, :)
κᵘn = @lift interior(κᵘt[$n], 1, 1, :)

lines!(axb, bn,  zc)
lines!(axu, un,  zc, label="u")
lines!(axu, vn,  zc, label="v", linestyle=:dash)
lines!(axκ, κᶜn, zf, label="κᶜ")
lines!(axκ, κᵘn, zf, label="κᵘ", linestyle=:dash)

axislegend(axu, position=:rb)

display(fig)

#=
record(fig, "convective_adjustment.mp4", 1:Nt, framerate=24) do nn
    n[] = nn
end
=#

