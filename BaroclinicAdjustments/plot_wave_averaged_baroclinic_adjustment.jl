using Oceananigans
using GLMakie

prefix = "wave_averaged_baroclinic_adjustment"

sides = (:east, :north, :bottom, :top)
slice_filenames = NamedTuple(side => prefix * "_$(side)_slice.jld2" for side in sides)

bts = (east   = FieldTimeSeries(slice_filenames.east, "b"),
       north  = FieldTimeSeries(slice_filenames.north, "b"),
       bottom = FieldTimeSeries(slice_filenames.bottom, "b"),
       top    = FieldTimeSeries(slice_filenames.top, "b"))

xb, yb, zb = nodes(bts.east)

xb = xb ./ 1e3 # convert m -> km
yb = yb ./ 1e3 # convert m -> km

grid = bts.east.grid
Nx, Ny, Nz = size(grid)
x = xnodes(bts.top)
y = ynodes(bts.top)
z = znodes(bts.top)

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

times = bts.top.times
n = length(times)

b_slices = (east   = interior(bts.east[n], 1, :, :),
            north  = interior(bts.north[n], :, 1, :),
            bottom = interior(bts.bottom[n], :, :, 1),
            top    = interior(bts.top[n], :, :, 1))

clims = 1.1 .* extrema(bts.top[n][:])

kwargs = (colorrange=clims, colormap=:deep)
surface!(ax, x_yz_east, y_yz, z_yz;    color = b_slices.east, kwargs...)
surface!(ax, x_xz, y_xz_north, z_xz;   color = b_slices.north, kwargs...)
surface!(ax, x_xy, y_xy, z_xy_bottom ; color = b_slices.bottom, kwargs...)
surface!(ax, x_xy, y_xy, z_xy_top;     color = b_slices.top, kwargs...)

display(fig)

