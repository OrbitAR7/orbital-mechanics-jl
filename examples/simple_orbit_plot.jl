#!/usr/bin/env julia
"""
Simple orbit plotting example

Creates a basic 3D visualization of a satellite orbit

Run: julia --project=. examples/simple_orbit_plot.jl
"""

# Load the orbital mechanics code
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using Plots
using LinearAlgebra
include("../src/OrbitalMechanics.jl")
using .OrbitalMechanics

println("Creating orbit visualization...")

# Define an ISS-like orbit
orbit = KeplerianElements(
    6780.0,           # Semi-major axis [km]
    0.0006,           # Eccentricity (almost circular)
    deg2rad(51.6),    # Inclination (ISS tilt)
    deg2rad(120.0),   # RAAN
    deg2rad(90.0),    # Argument of perigee
    deg2rad(0.0)      # Starting true anomaly
)

# Propagate orbit for one full period
T = orbital_period(orbit.a)
num_points = 500  # More points for smoother orbit
times = range(0, T, length=num_points)

println("Orbital period: $(T/60) minutes")
println("Altitude: $(orbit.a - R_EARTH) km")

# Store positions
positions = zeros(3, num_points)

for (i, t) in enumerate(times)
    # Propagate to this time
    oe_t = propagate_kepler(orbit, t)
    # Convert to position/velocity
    sv = keplerian_to_cartesian(oe_t)
    positions[:, i] = sv.r
end

# Create a 3D plot with equal aspect ratio
plot3d(
    positions[1, :], positions[2, :], positions[3, :],
    xlabel="X [km]", ylabel="Y [km]", zlabel="Z [km]",
    title="ISS Orbit (~400 km altitude)",
    label="Orbit",
    linewidth=2,
    color=:blue,
    legend=:topright,
    size=(800, 800),
    aspect_ratio=:equal,  # This is key for correct visualization!
    camera=(30, 30)
)

# Add Earth sphere at origin (radius = 6378.137 km)
θ = range(0, 2π, length=60)
φ = range(0, π, length=60)
R = R_EARTH  # Should be ~6378 km

x = [R * sin(p) * cos(t) for p in φ, t in θ]
y = [R * sin(p) * sin(t) for p in φ, t in θ]
z = [R * cos(p) for p in φ, t in θ]

surface!(x, y, z, 
    color=:lightblue, 
    alpha=0.5, 
    label="Earth (R≈$(round(Int, R)) km)",
    colorbar=false
)

# Mark the starting point
scatter3d!([positions[1, 1]], [positions[2, 1]], [positions[3, 1]],
    marker=:circle,
    markersize=6,
    color=:red,
    label="Start"
)

println("Plot created! Close the window when done viewing.")
savefig("orbit_3d.png")
println("Saved as 'orbit_3d.png'")
