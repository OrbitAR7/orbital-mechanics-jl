# Learning Orbital Mechanics with Julia

This is my project for learning Julia by implementing basic orbital mechanics calculations. It includes satellite orbit propagation, orbital transfers, and trajectory planning.

## What it does

- Convert between different orbit representations (Keplerian elements ↔ position/velocity)
- Predict where satellites will be over time
- Calculate fuel costs for orbital maneuvers (like going from LEO to GEO)
- Solve the Lambert problem (finding trajectories between two points)
- Visualize orbits in 3D and plot ground tracks showing where satellites pass over Earth

## Running the examples

```bash
# First time setup
cd orbital-mechanics-jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run the examples (text output)
julia --project=. examples/orbit_visualization.jl
julia --project=. examples/orbital_transfers.jl
julia --project=. examples/lambert_problem.jl

# Run plotting examples (creates PNG files)
julia --project=. examples/simple_orbit_plot.jl       # Creates orbit_3d.png
julia --project=. examples/plot_ground_track.jl       # Creates ground_track.png
```

## Simple example

```julia
using OrbitalMechanics

# Define an ISS-like orbit
orbit = KeplerianElements(
    6780.0,           # altitude + Earth radius [km]
    0.0006,           # how elliptical (0 = circle, <1 = ellipse)
    deg2rad(51.6),    # tilt relative to equator
    deg2rad(120.0),   # where orbit crosses equator going north
    deg2rad(90.0),    # where lowest point is
    deg2rad(0.0)      # current position in orbit
)

# Get position and velocity
sv = keplerian_to_cartesian(orbit)
println("Position: ", sv.r, " km")
println("Velocity: ", sv.v, " km/s")

# How much fuel to get to GEO?
transfer = hohmann_transfer(6778.0, 42164.0)
println("Δv needed: ", transfer.Δv_total, " km/s")
```

## What's in each example

**orbit_visualization.jl** - Demonstrates orbit propagation with J2 perturbations, shows how orbital elements evolve over time

**orbital_transfers.jl** - Compares Hohmann vs bi-elliptic transfers, analyzes plane change costs, shows constellation deployment strategies

**lambert_problem.jl** - Solves trajectory design problems, includes Earth-Mars transfer examples and rendezvous scenarios

**simple_orbit_plot.jl** - Creates a 3D visualization of ISS orbit around Earth (generates orbit_3d.png)

**plot_ground_track.jl** - Plots satellite ground track showing two complete orbits with Earth rotation effects (generates ground_track.png)
<img width="1200" height="600" alt="ground_track" src="https://github.com/user-attachments/assets/ca861525-5144-4840-8bcf-332297ee7efc" />

## Learning resources

If you're new to orbital mechanics:
- [Orbital Mechanics for Engineering Students](https://www.elsevier.com/books/orbital-mechanics-for-engineering-students/curtis/978-0-08-097747-8) by Howard Curtis
- [Fundamentals of Astrodynamics](https://www.dover.com/books/0486600610/) by Bate, Mueller, and White

For learning Julia:
- [Julia Documentation](https://docs.julialang.org/)
- [Think Julia](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html)

## Basic concepts

**Keplerian Elements** describe an orbit with 6 numbers:
- `a` - size of the orbit
- `e` - how circular it is
- `i` - tilt angle
- `Ω` - where it crosses the equator
- `ω` - orientation of the ellipse
- `ν` - current position

**State Vector** is position + velocity at a specific time

## Why I'm using Julia

- Fast (like C++) but easier to write
- Can use Greek letters in code (ν, ω, Ω) which matches the math
- Good for scientific computing
- Fun to learn!

## Notes

This is a learning project - there might be bugs or things that could be done better. If you spot something, feel free to suggest improvements!
