#!/usr/bin/env julia
"""
Example: Orbital Transfer Maneuvers

Demonstrates:
- Hohmann transfer
- Bi-elliptic transfer
- Plane change maneuvers
- Combined transfers
- Delta-v budget analysis

Author: Hamoud Alshammari
Date: 2025
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using Printf
include("../src/OrbitalMechanics.jl")
using .OrbitalMechanics

function print_header(title::String)
    println("\n", "="^60)
    println("  ", title)
    println("="^60)
end

function main()
    print_header("Orbital Transfer Maneuvers")
    
    # =========================================================================
    # LEO to GEO Transfer
    # =========================================================================
    println("\n  Scenario: LEO (400 km) → GEO Transfer")
    println("  ", "-"^50)
    
    r_leo = R_EARTH + 400.0    # LEO radius [km]
    r_geo = 42164.0            # GEO radius [km]
    
    @printf("\n  Initial orbit (LEO):  %.0f km altitude (r = %.0f km)\n", 
            r_leo - R_EARTH, r_leo)
    @printf("  Final orbit (GEO):    %.0f km altitude (r = %.0f km)\n", 
            r_geo - R_EARTH, r_geo)
    
    # --- Hohmann Transfer ---
    println("\n  [1] Hohmann Transfer:")
    hohmann = hohmann_transfer(r_leo, r_geo)
    
    @printf("      Δv₁ (LEO departure):     %.4f km/s\n", hohmann.Δv1)
    @printf("      Δv₂ (GEO insertion):     %.4f km/s\n", hohmann.Δv2)
    @printf("      Total Δv:                %.4f km/s\n", hohmann.Δv_total)
    @printf("      Transfer time:           %.2f hours\n", hohmann.tof / 3600)
    
    # --- Bi-elliptic Transfer ---
    # For LEO to GEO, bi-elliptic is NOT more efficient, but let's show it
    println("\n  [2] Bi-elliptic Transfer (rb = 100,000 km):")
    rb = 100000.0
    bielliptic = bielliptic_transfer(r_leo, r_geo, rb)
    
    @printf("      Δv₁ (LEO departure):     %.4f km/s\n", bielliptic.Δv1)
    @printf("      Δv₂ (at apogee):         %.4f km/s\n", bielliptic.Δv2)
    @printf("      Δv₃ (GEO insertion):     %.4f km/s\n", bielliptic.Δv3)
    @printf("      Total Δv:                %.4f km/s\n", bielliptic.Δv_total)
    @printf("      Transfer time:           %.2f hours\n", bielliptic.tof / 3600)
    
    # Comparison
    println("\n  Comparison:")
    @printf("      Hohmann Δv:      %.4f km/s (baseline)\n", hohmann.Δv_total)
    @printf("      Bi-elliptic Δv:  %.4f km/s (%+.1f%%)\n", 
            bielliptic.Δv_total, 
            100 * (bielliptic.Δv_total - hohmann.Δv_total) / hohmann.Δv_total)
    
    # =========================================================================
    # When is Bi-elliptic Better?
    # =========================================================================
    print_header("Bi-elliptic vs Hohmann Analysis")
    
    println("\n  Finding the crossover point...")
    println("\n  For bi-elliptic to beat Hohmann, r₂/r₁ must be > 11.94")
    println("\n  r₂/r₁ ratio    Hohmann Δv    Bi-elliptic Δv    Better")
    println("  -----------    ----------    --------------    ------")
    
    r1_fixed = R_EARTH + 300.0
    
    for ratio in [2, 5, 10, 11.94, 15, 20, 50]
        r2 = r1_fixed * ratio
        rb_opt = r2 * 2  # Not truly optimal, but reasonable
        
        h = hohmann_transfer(r1_fixed, r2)
        b = bielliptic_transfer(r1_fixed, r2, rb_opt)
        
        better = h.Δv_total < b.Δv_total ? "Hohmann" : "Bi-elliptic"
        
        @printf("  %11.2f    %10.4f    %14.4f    %s\n",
                ratio, h.Δv_total, b.Δv_total, better)
    end
    
    # =========================================================================
    # Plane Change Maneuver
    # =========================================================================
    print_header("Plane Change Analysis")
    
    println("\n  Cost of inclination change at different altitudes:")
    println("\n  Altitude [km]    v [km/s]    Δv for 10°    Δv for 28°")
    println("  -------------    --------    ----------    ----------")
    
    for alt in [400, 1000, 5000, 20000, 35786]
        r = R_EARTH + alt
        a = r  # Assume circular orbit
        v = orbital_velocity(r, a)
        
        Δv_10 = plane_change(v, deg2rad(10))
        Δv_28 = plane_change(v, deg2rad(28))
        
        @printf("  %13d    %8.3f    %10.3f    %10.3f\n",
                alt, v, Δv_10, Δv_28)
    end
    
    println("\n  Key insight: Plane changes are cheaper at higher altitudes!")
    println("  This is why GTO inclination changes are done at apogee.")
    
    # =========================================================================
    # Combined Maneuver: LEO to GEO with Inclination Change
    # =========================================================================
    print_header("LEO to GEO with Inclination Change")
    
    println("\n  Scenario: 400 km LEO @ 28.5° → GEO @ 0° (Cape Canaveral)")
    
    Δi = deg2rad(28.5)
    
    # Option 1: Plane change at LEO, then Hohmann
    v_leo = orbital_velocity(r_leo, r_leo)
    Δv_plane_leo = plane_change(v_leo, Δi)
    
    println("\n  [Option 1] Plane change at LEO, then Hohmann:")
    @printf("      Plane change Δv (LEO):   %.4f km/s\n", Δv_plane_leo)
    @printf("      Hohmann Δv:              %.4f km/s\n", hohmann.Δv_total)
    @printf("      Total:                   %.4f km/s\n", Δv_plane_leo + hohmann.Δv_total)
    
    # Option 2: Hohmann, then plane change at GEO
    v_geo = orbital_velocity(r_geo, r_geo)
    Δv_plane_geo = plane_change(v_geo, Δi)
    
    println("\n  [Option 2] Hohmann, then plane change at GEO:")
    @printf("      Hohmann Δv:              %.4f km/s\n", hohmann.Δv_total)
    @printf("      Plane change Δv (GEO):   %.4f km/s\n", Δv_plane_geo)
    @printf("      Total:                   %.4f km/s\n", hohmann.Δv_total + Δv_plane_geo)
    
    # Option 3: Combined maneuver at GEO insertion
    # This is the most efficient - combine the circularization and plane change
    a_transfer = (r_leo + r_geo) / 2
    v_transfer_apo = vis_viva(r_geo, a_transfer)  # Velocity on transfer orbit at apogee
    
    # Combined Δv using cosine rule
    Δv_combined = sqrt(v_geo^2 + v_transfer_apo^2 - 2*v_geo*v_transfer_apo*cos(Δi))
    
    println("\n  [Option 3] Combined maneuver at apogee:")
    @printf("      Δv₁ (LEO departure):     %.4f km/s\n", hohmann.Δv1)
    @printf("      Δv₂ (combined @ apo):    %.4f km/s\n", Δv_combined)
    @printf("      Total:                   %.4f km/s\n", hohmann.Δv1 + Δv_combined)
    
    println("\n  Combined maneuver saves Δv by doing plane change at lowest velocity!")
    
    # =========================================================================
    # LEO Constellation Deployment
    # =========================================================================
    print_header("LEO Constellation Deployment")
    
    println("\n  Deploying 24 satellites into 4 orbital planes (6 per plane)")
    println("  Orbit: 550 km altitude, 53° inclination")
    
    r_constellation = R_EARTH + 550.0
    v_const = orbital_velocity(r_constellation, r_constellation)
    T_const = orbital_period(r_constellation)
    
    @printf("\n  Orbital velocity:  %.3f km/s\n", v_const)
    @printf("  Orbital period:    %.1f min\n", T_const / 60)
    
    # RAAN spacing: 90° between planes
    # Phase spacing: 60° between satellites in same plane
    
    println("\n  RAAN drift from J2 for plane spacing:")
    
    # Different inclinations for different RAAN drift rates
    for inc in [53.0, 70.0, 87.0, 97.0]
        a = r_constellation
        e = 0.0
        p = a * (1 - e^2)
        n = sqrt(GM_EARTH / a^3)
        
        Ω_dot = -1.5 * J2_EARTH * (R_EARTH / p)^2 * n * cosd(inc)
        drift_per_day = rad2deg(Ω_dot) * 86400
        
        days_for_90deg = abs(90 / drift_per_day)
        
        @printf("      i = %5.1f°: RAAN drift = %+6.2f°/day (90° in %.1f days)\n",
                inc, drift_per_day, days_for_90deg)
    end
    
    # =========================================================================
    # Summary
    # =========================================================================
    print_header("Transfer Summary")
    
    println("""
    
  Key Results:
  
  1. LEO → GEO Hohmann transfer requires ~3.94 km/s total Δv
  2. Bi-elliptic beats Hohmann only when r₂/r₁ > 11.94
  3. Plane changes are much cheaper at higher altitudes
  4. Combined maneuvers (circularize + plane change) save significant Δv
  5. J2 precession can be used for constellation RAAN spacing
  
  Practical Application:
  - GTO launches from Cape Canaveral (28.5°) to GEO require ~4.2 km/s total
  - This is why GEO satellites need large propulsion systems
  - Electric propulsion can reduce mass but increases transfer time
    """)
end

main()
