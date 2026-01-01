#!/usr/bin/env julia
"""
Example: Lambert Problem - Interplanetary Trajectory Design

Demonstrates:
- Lambert problem solution
- Transfer trajectory computation
- Porkchop plot concept
- Earth-Mars transfer analysis

Author: Hamoud Alshammari
Date: 2025
"""

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using Printf
include("../src/OrbitalMechanics.jl")
using .OrbitalMechanics
using LinearAlgebra

function print_header(title::String)
    println("\n", "="^60)
    println("  ", title)
    println("="^60)
end

"""
Simple circular orbit position at given true anomaly
"""
function circular_orbit_position(r::Float64, θ::Float64)
    return [r * cos(θ), r * sin(θ), 0.0]
end

function main()
    print_header("Lambert Problem - Trajectory Design")
    
    # =========================================================================
    # Simple Example: Earth Orbit Transfer
    # =========================================================================
    println("\n  Example 1: Transfer between two points in Earth orbit")
    println("  ", "-"^50)
    
    # Two points on a circular orbit
    r_orbit = R_EARTH + 500.0  # 500 km altitude
    
    # Position 1: θ = 0°
    r1 = circular_orbit_position(r_orbit, 0.0)
    
    # Position 2: θ = 120°
    r2 = circular_orbit_position(r_orbit, deg2rad(120))
    
    @printf("\n  r₁ = [%.2f, %.2f, %.2f] km\n", r1...)
    @printf("  r₂ = [%.2f, %.2f, %.2f] km\n", r2...)
    @printf("  Transfer angle: 120°\n")
    
    # Solve for different times of flight
    println("\n  TOF [min]    Δv₁ [km/s]    Δv₂ [km/s]    Total [km/s]")
    println("  ---------    ----------    ----------    ------------")
    
    # Circular orbit velocity
    v_circ = sqrt(GM_EARTH / r_orbit)
    
    for tof_min in [30, 45, 60, 90, 120]
        tof = tof_min * 60.0
        
        try
            v1, v2 = lambert_solve(r1, r2, tof; μ=GM_EARTH)
            
            # Delta-v from circular orbit
            v1_circ = [-v_circ * sin(0.0), v_circ * cos(0.0), 0.0]
            v2_circ = [-v_circ * sin(deg2rad(120)), v_circ * cos(deg2rad(120)), 0.0]
            
            Δv1 = norm(v1 - v1_circ)
            Δv2 = norm(v2_circ - v2)
            
            @printf("  %9d    %10.4f    %10.4f    %12.4f\n",
                    tof_min, Δv1, Δv2, Δv1 + Δv2)
        catch e
            @printf("  %9d    No solution (transfer time too short)\n", tof_min)
        end
    end
    
    # =========================================================================
    # Interplanetary: Earth to Mars
    # =========================================================================
    print_header("Earth-Mars Transfer (Simplified)")
    
    # Simplified model: coplanar circular orbits
    r_earth = 1.0 * AU  # Earth orbit radius
    r_mars = 1.524 * AU  # Mars orbit radius
    
    @printf("\n  Earth orbit: %.3f AU (%.0f km)\n", r_earth/AU, r_earth)
    @printf("  Mars orbit:  %.3f AU (%.0f km)\n", r_mars/AU, r_mars)
    
    # Hohmann reference
    hohmann = hohmann_transfer(r_earth, r_mars; μ=GM_SUN)
    
    println("\n  Hohmann Transfer (reference):")
    @printf("      Δv at Earth: %.3f km/s\n", hohmann.Δv1)
    @printf("      Δv at Mars:  %.3f km/s\n", hohmann.Δv2)
    @printf("      Total Δv:    %.3f km/s\n", hohmann.Δv_total)
    @printf("      Transfer time: %.1f days\n", hohmann.tof / 86400)
    
    # Lambert solutions for different transfer times
    println("\n  Lambert Solutions for Different Transfer Times:")
    println("\n  TOF [days]    Δv_depart    Δv_arrive    Total [km/s]    vs Hohmann")
    println("  ----------    ---------    ---------    ------------    ----------")
    
    # Earth at departure: θ = 0°
    r1_sun = circular_orbit_position(r_earth, 0.0)
    v_earth = sqrt(GM_SUN / r_earth)
    
    for tof_days in [150, 200, 250, 259, 300, 350, 400]
        tof = tof_days * 86400.0
        
        # Mars position at arrival (simplified: assume Mars moves ~0.524°/day)
        mars_motion = 0.524 * tof_days  # degrees
        θ_mars = deg2rad(mars_motion)
        r2_sun = circular_orbit_position(r_mars, θ_mars)
        v_mars = sqrt(GM_SUN / r_mars)
        
        try
            v1, v2 = lambert_solve(r1_sun, r2_sun, tof; μ=GM_SUN)
            
            # Departure Δv (from Earth's circular orbit)
            v_earth_vec = [0.0, v_earth, 0.0]  # Earth velocity at θ=0
            Δv_depart = norm(v1 - v_earth_vec)
            
            # Arrival Δv (to Mars's circular orbit)
            v_mars_vec = [-v_mars * sin(θ_mars), v_mars * cos(θ_mars), 0.0]
            Δv_arrive = norm(v_mars_vec - v2)
            
            total = Δv_depart + Δv_arrive
            diff = total - hohmann.Δv_total
            
            @printf("  %10d    %9.3f    %9.3f    %12.3f    %+.3f\n",
                    tof_days, Δv_depart, Δv_arrive, total, diff)
        catch e
            @printf("  %10d    No solution\n", tof_days)
        end
    end
    
    # =========================================================================
    # Porkchop Plot Concept
    # =========================================================================
    print_header("Porkchop Plot Concept")
    
    println("""
    
  A porkchop plot shows total Δv as a function of:
  - Departure date (x-axis)
  - Arrival date or TOF (y-axis)
  
  The contours form a shape resembling a porkchop.
  
  Optimal transfers appear as local minima in the plot.
  
  For Earth-Mars:
  - Launch windows occur approximately every 26 months
  - Optimal TOF is around 250-280 days
  - Minimum energy transfer (Hohmann) takes ~259 days
    """)
    
    # Simplified porkchop: Δv vs TOF
    println("  Simplified Porkchop (Δv vs TOF):")
    println("\n  TOF [days]    Total Δv [km/s]")
    println("  ----------    ---------------")
    
    min_dv = Inf
    optimal_tof = 0
    
    for tof_days in 150:25:450
        tof = tof_days * 86400.0
        
        mars_motion = 0.524 * tof_days
        θ_mars = deg2rad(mars_motion)
        r2_sun = circular_orbit_position(r_mars, θ_mars)
        
        try
            v1, v2 = lambert_solve(r1_sun, r2_sun, tof; μ=GM_SUN)
            
            v_earth_vec = [0.0, v_earth, 0.0]
            v_mars_vec = [-sqrt(GM_SUN/r_mars) * sin(θ_mars), 
                          sqrt(GM_SUN/r_mars) * cos(θ_mars), 0.0]
            
            Δv_total = norm(v1 - v_earth_vec) + norm(v_mars_vec - v2)
            
            # ASCII visualization
            bar_len = min(50, max(1, round(Int, (Δv_total - 5) * 5)))
            bar = "█"^bar_len
            
            @printf("  %10d    %6.3f  %s\n", tof_days, Δv_total, bar)
            
            if Δv_total < min_dv
                min_dv = Δv_total
                optimal_tof = tof_days
            end
        catch
            continue
        end
    end
    
    @printf("\n  Optimal TOF: %d days with Δv = %.3f km/s\n", optimal_tof, min_dv)
    @printf("  Hohmann TOF: %d days with Δv = %.3f km/s\n", 
            round(Int, hohmann.tof/86400), hohmann.Δv_total)
    
    # =========================================================================
    # Rendezvous Problem
    # =========================================================================
    print_header("Rendezvous Problem")
    
    println("\n  Scenario: Chase satellite catching target satellite")
    
    # Target in circular orbit
    r_target = R_EARTH + 400.0
    T_target = orbital_period(r_target)
    
    # Initial phase angle (chaser behind target)
    phase_angle = deg2rad(10)  # 10° behind
    
    @printf("\n  Target orbit: %.0f km altitude\n", r_target - R_EARTH)
    @printf("  Orbital period: %.2f min\n", T_target / 60)
    @printf("  Initial phase: %.1f° behind target\n", rad2deg(phase_angle))
    
    # For rendezvous, chaser must arrive when target reaches the same position
    # This requires solving Lambert with constraint
    
    println("\n  Transfer options (to same point, different # of target orbits):")
    println("\n  Target Orbits    TOF [min]    Δv₁ [km/s]    Δv₂ [km/s]")
    println("  -------------    ---------    ----------    ----------")
    
    r1_chase = circular_orbit_position(r_target, -phase_angle)
    r2_chase = circular_orbit_position(r_target, 0.0)  # Target position
    
    v_target = sqrt(GM_EARTH / r_target)
    v1_circ = [v_target * sin(phase_angle), v_target * cos(phase_angle), 0.0]
    v2_circ = [0.0, v_target, 0.0]
    
    for n_orbits in 1:5
        # Time for target to complete n orbits minus phase angle
        tof = T_target * n_orbits - phase_angle / (2π) * T_target
        
        try
            v1, v2 = lambert_solve(r1_chase, r2_chase, tof; μ=GM_EARTH)
            
            Δv1 = norm(v1 - v1_circ)
            Δv2 = norm(v2_circ - v2)
            
            @printf("  %13d    %9.2f    %10.4f    %10.4f\n",
                    n_orbits, tof/60, Δv1, Δv2)
        catch
            @printf("  %13d    No solution\n", n_orbits)
        end
    end
    
    # =========================================================================
    # Summary
    # =========================================================================
    print_header("Summary")
    
    println("""
    
  Lambert's Problem Applications:
  
  1. Interplanetary trajectory design
     - Find optimal launch windows
     - Minimize total Δv for transfers
  
  2. Rendezvous operations
     - Space station resupply
     - Satellite servicing
  
  3. Mission planning
     - Porkchop plots for launch window analysis
     - Multi-flyby trajectory optimization
  
  Key Insights:
  - Hohmann is optimal only for specific conditions
  - Shorter transfers cost more Δv
  - Launch window selection critical for interplanetary
  - Rendezvous timing affects fuel cost
    """)
end

main()
