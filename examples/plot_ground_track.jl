#!/usr/bin/env julia
"""
Ground track plotting example

Shows where a satellite passes over Earth

Run: julia --project=. examples/plot_ground_track.jl
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using Plots
using Printf
include("../src/OrbitalMechanics.jl")
using .OrbitalMechanics

println("Computing ground track...")

# ISS orbit
orbit = KeplerianElements(
    6780.0,
    0.0006,
    deg2rad(51.6),
    deg2rad(120.0),
    deg2rad(90.0),
    deg2rad(0.0)
)

# Compute ground track for 2 orbits with high resolution
T = orbital_period(orbit.a)
duration = 2 * T
dt = 10.0  # 10 second steps for smooth track

println("Period: $(T/60) minutes, computing $(duration/60) minutes...")

ground_track = compute_ground_track(orbit, duration, dt)

# Extract lat/lon and wrap longitude to [-180, 180]
lats_all = [rad2deg(gt.lat) for gt in ground_track]
lons_raw = [rad2deg(gt.lon) for gt in ground_track]

# Build arrays with NaN breaks for dateline crossings
lons = Float64[]
lats = Float64[]
push!(lons, mod(lons_raw[1] + 180, 360) - 180)
push!(lats, lats_all[1])

for i in 2:length(lons_raw)
    lon_wrapped = mod(lons_raw[i] + 180, 360) - 180
    # Check for dateline crossing (jump > 180 degrees)
    if abs(lon_wrapped - lons[end]) > 180
        # Insert NaN to break the line at dateline
        push!(lons, NaN)
        push!(lats, NaN)
    end
    push!(lons, lon_wrapped)
    push!(lats, lats_all[i])
end

println("Computed $(length(ground_track)) points")

# Filter out NaN values for stats
valid_lats = filter(!isnan, lats)
valid_lons = filter(!isnan, lons)
println("Valid points: $(length(valid_lats)) latitudes, $(length(valid_lons)) longitudes")
if length(valid_lats) > 0
    println("Latitude range: [$(minimum(valid_lats)), $(maximum(valid_lats))]")
    println("Longitude range: [$(minimum(valid_lons)), $(maximum(valid_lons))]")
else
    println("WARNING: All latitudes are NaN!")
end

# Create plot
plot(lons, lats,
    xlabel="Longitude [deg]",
    ylabel="Latitude [deg]",
    title="ISS Ground Track (2 Orbits, ~$(round(Int, duration/60)) min)",
    label="Ground Track",
    linewidth=2,
    color=:red,
    xlims=(-180, 180),
    ylims=(-60, 60),
    size=(1200, 600),
    grid=true,
    legend=:topright
)

# Add start marker
scatter!([lons[1]], [lats[1]],
    marker=:circle,
    markersize=8,
    color=:green,
    label="Start"
)

# Add gridlines for reference
hline!([-51.6, 51.6], 
    linestyle=:dash, 
    color=:gray, 
    alpha=0.5, 
    label="Max Latitude (±51.6°)"
)

hline!([0], 
    linestyle=:dash, 
    color=:black, 
    alpha=0.3, 
    label="Equator"
)

vline!([0], 
    linestyle=:dash, 
    color=:gray, 
    alpha=0.5, 
    label="Prime Meridian"
)

savefig("ground_track.png")
println("\nSaved as 'ground_track.png'")
