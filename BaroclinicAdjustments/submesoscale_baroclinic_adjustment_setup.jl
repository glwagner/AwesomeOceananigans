# # Submesoscale baroclinic adjustment
#
# Roughly follows
#
# Hamlington et al.,
# "Langmuir-Submesoscale Interactions: Descriptive Analysis of Multiscale Frontal Spindown Simulations",
# Journal of Physical Oceananography, 2014

using Oceananigans
using Oceananigans.Units
using Printf

# Domain parameters
Lx = 10kilometers # east-west extent [m]
Ly = 10kilometers # north-south extent [m]
Lz = 160           # depth [m]

# Front / initial condition parameters
f₀ = 1e-4          # [s⁻¹] Coriolis parameter 
hᵢ = 50            # Initial mixed layer depth (m)
N² = 1e-8          # [s⁻²] sub-mixed-layer buoyancy frequency / stratification
M² = 2e-8          # [s⁻²] horizontal buoyancy gradient
Δy = 500           # [m] width of the region of the front
Δb = Δy * M²       # buoyancy jump associated with the front
ϵb = 1e-2 * Δb     # noise amplitude

# Note that via thermal wind, f ∂z u ~ ∇ₕ b ~ M².
# Therefore the Richardson number,
#
# Ri ≡ N² / | ∂z u |²
#
# can be written

Ri = N² * f² / M²^2

# Eddy length scales and growth rates:

L = 2π * M² * hᵢ / f^2 * sqrt(2/5) * sqrt(1 + Ri)
τ = sqrt(54/5) * sqrt(1 + Ri) / f

@info """ Some simulation parameters:

         Richardson number: $Ri
         Eddy length scale: $L
    Eddy growth time-scale: $(prettytime(τ))
"""

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

ramp(y, Δy) = min(max(0, y/Δy + 1/2), 1)

function bᵢ(x, y, z)
    z̃ = min(z, -hᵢ)
    return N² * z̃ + Δb * ramp(y, Δy) + ϵb * randn()
end


