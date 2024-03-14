using GLMakie
using Oceananigans
using Statistics

Nx = 128
Nz = 40
name = "modified_hydrostatic_submesoscale_baroclinic_adjustment"

catke_filename    =    "catke_$(name)_Nx$(Nx)_Nz$(Nz)_top_slice.jld2"
ri_based_filename = "ri_based_$(name)_Nx$(Nx)_Nz$(Nz)_top_slice.jld2"

bt_catke = FieldTimeSeries(catke_filename, "b")
ζt_catke = FieldTimeSeries(catke_filename, "ζ")
#et_catke = FieldTimeSeries(catke_filename, "e")

bt_ri_based = FieldTimeSeries(ri_based_filename, "b")
ζt_ri_based = FieldTimeSeries(ri_based_filename, "ζ")

xb, yb, zb = nodes(bt_catke)
xζ, yζ, zζ = nodes(ζt_catke)

xb = xb ./ 1e3
yb = yb ./ 1e3
xζ = xζ ./ 1e3
yζ = yζ ./ 1e3

t = bt_catke.times
Nt = length(bt_catke)
f = 1e-4

fig = Figure(size=(1200, 1000))

axb_catke    = Axis(fig[1, 1], aspect=1, xlabel="x (km)", ylabel="y (km)")
#axe_catke    = Axis(fig[1, 2], aspect=1, xlabel="x (km)", ylabel="y (km)")
axζ_catke    = Axis(fig[1, 2], aspect=1, xlabel="x (km)", ylabel="y (km)")
axb_ri_based = Axis(fig[2, 1], aspect=1, xlabel="x (km)", ylabel="y (km)")
axζ_ri_based = Axis(fig[2, 2], aspect=1, xlabel="x (km)", ylabel="y (km)")

slider = Slider(fig[3, 1:2], range=1:Nt, startvalue=1)
n = slider.value

bn_catke = @lift interior(bt_catke[$n], :, :, 1)
en_catke = @lift interior(et_catke[$n], :, :, 1)
ζn_catke = @lift interior(ζt_catke[$n], :, :, 1) ./ f

bn_ri_based = @lift interior(bt_ri_based[$n], :, :, 1)
ζn_ri_based = @lift interior(ζt_ri_based[$n], :, :, 1) ./ f

δb = 5e-6

bcolorrange = @lift begin
    #bmin = minimum(bt_ri_based[$n])
    #bmax = maximum(bt_ri_based[$n])
    #(bmin, bmax)
    bavg = mean(bt_ri_based[$n])
    (bavg - δb, bavg + δb)
end

     heatmap!(axb_catke,    xb, yb, bn_catke,    colormap=:thermal; colorrange=bcolorrange)
     #heatmap!(axe_catke,    xb, yb, en_catke,    colormap=:solar, colorrange=(0, 1e-5))
hm = heatmap!(axb_ri_based, xb, yb, bn_ri_based, colormap=:thermal; colorrange=bcolorrange)
Colorbar(fig[1:2, 0], hm, flipaxis=false, label="Buoyancy")

     heatmap!(axζ_catke,    xζ, yζ, ζn_catke,    colormap=:balance, colorrange=(-2, 2))
hm = heatmap!(axζ_ri_based, xζ, yζ, ζn_ri_based, colormap=:balance, colorrange=(-2, 2))
Colorbar(fig[1:2, 3], hm, label="Rossby number")

display(fig)

#=
record(fig, "$name.mp4") do io
    for nn = 1:Nt
        n[] = nn
        @info "Drawing frame $nn of $Nt..."
        recordframe!(io)  # record a new frame
    end
end
=#
