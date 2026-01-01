#!/usr/bin/env julia
"""
Orbit visualization example

Playing around with orbital mechanics - seeing how satellites move
and how Earth's shape (J2) affects their orbits.

Run: julia --project=. examples/orbit_visualization.jl
"""

# Load the orbital mechanics code
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using Printf
using LinearAlgebra
include("../src/OrbitalMechanics.jl")
using .OrbitalMechanics

# Helper functions to make output look nice
function print_header(title::String)
    println("\n", "="^60)
    println("  ", title)
    println("="^60)
end

function print_section(title::String)
    println("\n  ", title)
    println("  ", "-"^50)
end

# Main code
function main()
    print_header("Orbital Mechanics in Julia")
    
    # Let's model the ISS orbit
    print_section("1. Initial Orbit Definition (ISS-like)")
    
    # ISS orbital elements
    oe_iss = KeplerianElements(
        6780.0,           # a: Semi-major axis [km]
        0.0006,           # e: Eccentricity
        deg2rad(51.6),    # i: Inclination [rad]
        deg2rad(120.0),   # Ω: RAAN [rad]
        deg2rad(90.0),    # ω: Argument of perigee [rad]
        deg2rad(0.0)      # ν: True anomaly [rad]
    )
    
    println("\n", oe_iss)
    
    # Calculate some basic orbital properties
    T = orbital_period(oe_iss.a)
    h = specific_angular_momentum(oe_iss)
    ε = orbital_energy(oe_iss.a)
    
    @printf("\n  Orbital Period:    %.2f min\n", T / 60)
    @printf("  Angular Momentum:  %.2f km²/s\n", h)
    @printf("  Specific Energy:   %.4f km²/s²\n", ε)
    @printf("  Apogee Altitude:   %.2f km\n", oe_iss.a * (1 + oe_iss.e) - R_EARTH)
    @printf("  Perigee Altitude:  %.2f km\n", oe_iss.a * (1 - oe_iss.e) - R_EARTH)
    
    # Converting to position and velocity vectors
    print_section("2. State Vector at Epoch")
    
    sv = keplerian_to_cartesian(oe_iss)
    println("\n", sv)
    
    # Check if conversion works both ways (it should!)
    oe_recovered = cartesian_to_keplerian(sv)
    @printf("\n  Round-trip error in a: %.6e km\n", abs(oe_iss.a - oe_recovered.a))
    
    # Propagate the orbit over time
    print_section("3. Two-Body Propagation (1 orbit)")
    
    println("\n  Time [min]    r [km]      v [km/s]    ν [deg]")
    println("  ----------    ------      --------    -------")
    
    for fraction in 0.0:0.125:1.0
        t = fraction * T
        oe_prop = propagate_kepler(oe_iss, t)
        sv_prop = keplerian_to_cartesian(oe_prop)
        
        @printf("  %10.2f    %.2f    %.4f    %7.2f\n",
                t/60, norm(sv_prop.r), norm(sv_prop.v), rad2deg(oe_prop.ν))
    end
    
    # =========================================================================
    # J2 Perturbation Effects
    # =========================================================================
    print_section("4. J2 Perturbation Effects (24 hours)")
    
    duration = 24 * 3600.0  # 24 hours
    
    # Compute J2 drift rates
    p = oe_iss.a * (1 - oe_iss.e^2)
    n = sqrt(GM_EARTH / oe_iss.a^3)
    J2_factor = 1.5 * J2_EARTH * (R_EARTH / p)^2 * n
    
    Ω_dot = -J2_factor * cos(oe_iss.i)  # rad/s
    ω_dot = J2_factor * (2 - 2.5 * sin(oe_iss.i)^2)  # rad/s
    
    @printf("\n  RAAN precession rate:      %.4f deg/day\n", rad2deg(Ω_dot) * 86400)
    @printf("  Arg. perigee drift rate:   %.4f deg/day\n", rad2deg(ω_dot) * 86400)
    
    # Propagate and show element evolution
    println("\n  Time [hr]    Ω [deg]     ω [deg]     Change")
    println("  ---------    -------     -------     ------")
    
    Ω_init = rad2deg(oe_iss.Ω)
    ω_init = rad2deg(oe_iss.ω)
    
    for hours in [0, 4, 8, 12, 16, 20, 24]
        t = hours * 3600.0
        oe_j2 = propagate_j2(oe_iss, t)
        
        ΔΩ = rad2deg(oe_j2.Ω) - Ω_init
        Δω = rad2deg(oe_j2.ω) - ω_init
        
        @printf("  %9d    %7.2f     %7.2f     ΔΩ=%.2f°, Δω=%.2f°\n",
                hours, rad2deg(oe_j2.Ω), rad2deg(oe_j2.ω), ΔΩ, Δω)
    end
    
    # =========================================================================
    # Ground Track
    # =========================================================================
    print_section("5. Ground Track (1 orbit)")
    
    track = compute_ground_track(oe_iss, T, T/20; use_j2=true)
    
    println("\n  Time [min]    Lat [deg]    Lon [deg]    Alt [km]")
    println("  ----------    ---------    ---------    --------")
    
    for (i, gt) in enumerate(track)
        if i % 2 == 1 || i == length(track)
            @printf("  %10.2f    %9.2f    %9.2f    %8.2f\n",
                    gt.t/60, rad2deg(gt.lat), rad2deg(gt.lon), gt.alt)
        end
    end
    
    # =========================================================================
    # Numerical vs Analytical Comparison
    # =========================================================================
    print_section("6. Numerical vs Analytical Propagation")
    
    # Numerical propagation with Cowell's method
    times, states = propagate_cowell(sv, T, 10.0; include_j2=true)
    
    println("\n  Comparison after 1 orbit:")
    
    # Final states
    sv_analytical = keplerian_to_cartesian(propagate_j2(oe_iss, T))
    sv_numerical = states[end]
    
    Δr = norm(sv_analytical.r - sv_numerical.r)
    Δv = norm(sv_analytical.v - sv_numerical.v)
    
    @printf("  Position difference: %.6f km\n", Δr)
    @printf("  Velocity difference: %.6e km/s\n", Δv)
    
    # =========================================================================
    # Different Orbit Types
    # =========================================================================
    print_section("7. Orbit Type Comparison")
    
    orbits = [
        ("LEO (ISS)",        6780.0,  0.0006, 51.6),
        ("GPS (MEO)",       26560.0,  0.01,   55.0),
        ("GEO",             42164.0,  0.0,     0.0),
        ("Molniya (HEO)",   26554.0,  0.74,   63.4),
        ("Sun-Sync",         7080.0,  0.001,  98.2),
    ]
    
    println("\n  Orbit Type        a [km]       Period       v_peri [km/s]")
    println("  ----------        ------       ------       -------------")
    
    for (name, a, e, inc) in orbits
        T_orb = orbital_period(a)
        r_peri = a * (1 - e)
        v_peri = vis_viva(r_peri, a)
        
        if T_orb < 3600
            period_str = @sprintf("%.1f min", T_orb / 60)
        elseif T_orb < 86400
            period_str = @sprintf("%.1f hr", T_orb / 3600)
        else
            period_str = @sprintf("%.1f days", T_orb / 86400)
        end
        
        @printf("  %-16s  %8.1f  %10s       %.3f\n", name, a, period_str, v_peri)
    end
    
    # =========================================================================
    # Summary
    # =========================================================================
    print_header("Analysis Complete")
    
    println("""
    
  Key Observations:
  
  1. J2 causes RAAN precession (westward for prograde, ~$(round(rad2deg(Ω_dot) * 86400, digits=2))°/day for ISS)
  2. Argument of perigee also drifts due to J2
  3. Numerical and analytical methods agree well for short propagation
  4. Different orbit types serve different mission needs
  
  Julia Benefits Demonstrated:
  - Clean mathematical syntax (≈, ², Greek letters)
  - Fast numerical computation
  - Type system for safety
  - Easy interoperability
    """)
end

# Run main
main()
