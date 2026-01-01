"""
    OrbitalMechanics

A Julia package for orbital mechanics computations including:
- Keplerian orbital elements and conversions
- Orbit propagation (two-body, J2 perturbation)
- Lambert problem solver
- Hohmann and bi-elliptic transfers
- Ground track computation
- Visibility analysis

Author: Hamoud Alshammari
Date: 2025
"""
module OrbitalMechanics

using LinearAlgebra
using Printf

export 
    # Constants
    GM_EARTH, R_EARTH, J2_EARTH, OMEGA_EARTH,
    GM_SUN, AU,
    # Types
    KeplerianElements, StateVector, GroundTrack,
    # Conversions
    keplerian_to_cartesian, cartesian_to_keplerian,
    ecef_to_lla, lla_to_ecef,
    # Propagation
    propagate_kepler, propagate_j2, propagate_cowell,
    # Maneuvers
    hohmann_transfer, bielliptic_transfer, plane_change,
    # Lambert
    lambert_solve,
    # Analysis
    orbital_period, orbital_energy, orbital_velocity,
    specific_angular_momentum, vis_viva,
    compute_ground_track,
    # Anomalies
    mean_to_eccentric, eccentric_to_true, true_to_eccentric,
    mean_to_true, true_to_mean

# =============================================================================
# Constants
# =============================================================================

"Earth gravitational parameter [km³/s²]"
const GM_EARTH = 398600.4418

"Earth equatorial radius [km]"
const R_EARTH = 6378.137

"Earth J2 coefficient"
const J2_EARTH = 0.00108263

"Earth rotation rate [rad/s]"
const OMEGA_EARTH = 7.2921159e-5

"Sun gravitational parameter [km³/s²]"
const GM_SUN = 1.32712440018e11

"Astronomical Unit [km]"
const AU = 1.495978707e8

# =============================================================================
# Types
# =============================================================================

"""
    KeplerianElements

Classical orbital elements.

# Fields
- `a`: Semi-major axis [km]
- `e`: Eccentricity [-]
- `i`: Inclination [rad]
- `Ω`: Right ascension of ascending node [rad]
- `ω`: Argument of periapsis [rad]
- `ν`: True anomaly [rad]
"""
struct KeplerianElements
    a::Float64  # Semi-major axis [km]
    e::Float64  # Eccentricity
    i::Float64  # Inclination [rad]
    Ω::Float64  # RAAN [rad]
    ω::Float64  # Argument of periapsis [rad]
    ν::Float64  # True anomaly [rad]
end

"""
    StateVector

Position and velocity in Cartesian coordinates.

# Fields
- `r`: Position vector [km]
- `v`: Velocity vector [km/s]
"""
struct StateVector
    r::Vector{Float64}  # Position [km]
    v::Vector{Float64}  # Velocity [km/s]
end

"""
    GroundTrack

Ground track point with time.

# Fields
- `t`: Time [s]
- `lat`: Latitude [rad]
- `lon`: Longitude [rad]
- `alt`: Altitude [km]
"""
struct GroundTrack
    t::Float64
    lat::Float64
    lon::Float64
    alt::Float64
end

# =============================================================================
# Anomaly Conversions
# =============================================================================

"""
    mean_to_eccentric(M, e; tol=1e-12, max_iter=50)

Convert mean anomaly to eccentric anomaly using Newton-Raphson.

# Arguments
- `M`: Mean anomaly [rad]
- `e`: Eccentricity
- `tol`: Convergence tolerance
- `max_iter`: Maximum iterations

# Returns
- `E`: Eccentric anomaly [rad]
"""
function mean_to_eccentric(M::Float64, e::Float64; tol::Float64=1e-12, max_iter::Int=50)
    # Initial guess
    E = M > π ? M - e : M + e
    
    for _ in 1:max_iter
        f = E - e * sin(E) - M
        f_prime = 1 - e * cos(E)
        dE = -f / f_prime
        E += dE
        
        if abs(dE) < tol
            break
        end
    end
    
    return E
end

"""
    eccentric_to_true(E, e)

Convert eccentric anomaly to true anomaly.
"""
function eccentric_to_true(E::Float64, e::Float64)
    ν = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2))
    # Normalize to [0, 2π]
    return mod(ν, 2π)
end

"""
    true_to_eccentric(ν, e)

Convert true anomaly to eccentric anomaly.
"""
function true_to_eccentric(ν::Float64, e::Float64)
    return 2 * atan(sqrt((1 - e) / (1 + e)) * tan(ν / 2))
end

"""
    mean_to_true(M, e)

Convert mean anomaly to true anomaly.
"""
function mean_to_true(M::Float64, e::Float64)
    E = mean_to_eccentric(M, e)
    return eccentric_to_true(E, e)
end

"""
    true_to_mean(ν, e)

Convert true anomaly to mean anomaly.
"""
function true_to_mean(ν::Float64, e::Float64)
    E = true_to_eccentric(ν, e)
    M = E - e * sin(E)
    # Normalize to [0, 2π]
    return mod(M, 2π)
end

# =============================================================================
# Orbital Element Conversions
# =============================================================================

"""
    keplerian_to_cartesian(oe; μ=GM_EARTH)

Convert Keplerian orbital elements to Cartesian state vector.

# Arguments
- `oe`: KeplerianElements
- `μ`: Gravitational parameter [km³/s²]

# Returns
- `StateVector`: Position and velocity vectors
"""
function keplerian_to_cartesian(oe::KeplerianElements; μ::Float64=GM_EARTH)
    a, e, i, Ω, ω, ν = oe.a, oe.e, oe.i, oe.Ω, oe.ω, oe.ν
    
    # Semi-latus rectum
    p = a * (1 - e^2)
    
    # Position in perifocal frame
    r_pqw = p / (1 + e * cos(ν)) * [cos(ν), sin(ν), 0.0]
    
    # Velocity in perifocal frame
    v_pqw = sqrt(μ / p) * [-sin(ν), e + cos(ν), 0.0]
    
    # Rotation matrix: perifocal -> inertial
    R = rotation_matrix(Ω, i, ω)
    
    r = R * r_pqw
    v = R * v_pqw
    
    return StateVector(r, v)
end

"""
    cartesian_to_keplerian(sv; μ=GM_EARTH)

Convert Cartesian state vector to Keplerian orbital elements.
"""
function cartesian_to_keplerian(sv::StateVector; μ::Float64=GM_EARTH)
    r, v = sv.r, sv.v
    r_mag = norm(r)
    v_mag = norm(v)
    
    # Specific angular momentum
    h = cross(r, v)
    h_mag = norm(h)
    
    # Node vector
    n = cross([0.0, 0.0, 1.0], h)
    n_mag = norm(n)
    
    # Eccentricity vector
    e_vec = ((v_mag^2 - μ / r_mag) * r - dot(r, v) * v) / μ
    e = norm(e_vec)
    
    # Specific orbital energy
    ε = v_mag^2 / 2 - μ / r_mag
    
    # Semi-major axis
    a = -μ / (2 * ε)
    
    # Inclination
    i = acos(clamp(h[3] / h_mag, -1.0, 1.0))
    
    # RAAN
    if n_mag > 1e-10
        Ω = acos(clamp(n[1] / n_mag, -1.0, 1.0))
        if n[2] < 0
            Ω = 2π - Ω
        end
    else
        Ω = 0.0
    end
    
    # Argument of periapsis
    if n_mag > 1e-10 && e > 1e-10
        ω = acos(clamp(dot(n, e_vec) / (n_mag * e), -1.0, 1.0))
        if e_vec[3] < 0
            ω = 2π - ω
        end
    else
        ω = 0.0
    end
    
    # True anomaly
    if e > 1e-10
        ν = acos(clamp(dot(e_vec, r) / (e * r_mag), -1.0, 1.0))
        if dot(r, v) < 0
            ν = 2π - ν
        end
    else
        ν = 0.0
    end
    
    return KeplerianElements(a, e, i, Ω, ω, ν)
end

"""
    rotation_matrix(Ω, i, ω)

Compute rotation matrix from perifocal to inertial frame.
"""
function rotation_matrix(Ω::Float64, i::Float64, ω::Float64)
    cΩ, sΩ = cos(Ω), sin(Ω)
    ci, si = cos(i), sin(i)
    cω, sω = cos(ω), sin(ω)
    
    return [
        cΩ*cω - sΩ*sω*ci  -cΩ*sω - sΩ*cω*ci   sΩ*si;
        sΩ*cω + cΩ*sω*ci  -sΩ*sω + cΩ*cω*ci  -cΩ*si;
        sω*si              cω*si               ci
    ]
end

# =============================================================================
# Coordinate Transformations
# =============================================================================

"""
    ecef_to_lla(r_ecef)

Convert ECEF coordinates to geodetic (WGS84).

# Returns
- `(lat, lon, alt)`: Latitude [rad], Longitude [rad], Altitude [km]
"""
function ecef_to_lla(r_ecef::Vector{Float64})
    x, y, z = r_ecef
    
    # WGS84 parameters
    a = R_EARTH
    f = 1 / 298.257223563
    e² = 2f - f^2
    
    lon = atan(y, x)
    
    # Iterative solution for latitude
    p = sqrt(x^2 + y^2)
    lat = atan(z, p * (1 - e²))  # Initial guess
    
    for _ in 1:10
        N = a / sqrt(1 - e² * sin(lat)^2)
        lat = atan(z + e² * N * sin(lat), p)
    end
    
    N = a / sqrt(1 - e² * sin(lat)^2)
    alt = p / cos(lat) - N
    
    return (lat, lon, alt)
end

"""
    lla_to_ecef(lat, lon, alt)

Convert geodetic coordinates to ECEF.
"""
function lla_to_ecef(lat::Float64, lon::Float64, alt::Float64)
    a = R_EARTH
    f = 1 / 298.257223563
    e² = 2f - f^2
    
    N = a / sqrt(1 - e² * sin(lat)^2)
    
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = (N * (1 - e²) + alt) * sin(lat)
    
    return [x, y, z]
end

# =============================================================================
# Orbital Mechanics Functions
# =============================================================================

"""
    orbital_period(a; μ=GM_EARTH)

Compute orbital period from semi-major axis.
"""
orbital_period(a::Float64; μ::Float64=GM_EARTH) = 2π * sqrt(a^3 / μ)

"""
    orbital_energy(a; μ=GM_EARTH)

Compute specific orbital energy.
"""
orbital_energy(a::Float64; μ::Float64=GM_EARTH) = -μ / (2a)

"""
    orbital_velocity(r, a; μ=GM_EARTH)

Compute orbital velocity using vis-viva equation.
"""
orbital_velocity(r::Float64, a::Float64; μ::Float64=GM_EARTH) = sqrt(μ * (2/r - 1/a))

"""
    vis_viva(r, a; μ=GM_EARTH)

Vis-viva equation: v² = μ(2/r - 1/a)
"""
vis_viva(r::Float64, a::Float64; μ::Float64=GM_EARTH) = sqrt(μ * (2/r - 1/a))

"""
    specific_angular_momentum(oe; μ=GM_EARTH)

Compute specific angular momentum magnitude.
"""
function specific_angular_momentum(oe::KeplerianElements; μ::Float64=GM_EARTH)
    p = oe.a * (1 - oe.e^2)
    return sqrt(μ * p)
end

# =============================================================================
# Orbit Propagation
# =============================================================================

"""
    propagate_kepler(oe, Δt; μ=GM_EARTH)

Propagate orbit using Kepler's equation (two-body only).

# Arguments
- `oe`: Initial orbital elements
- `Δt`: Time to propagate [s]
- `μ`: Gravitational parameter

# Returns
- New `KeplerianElements` at time t + Δt
"""
function propagate_kepler(oe::KeplerianElements, Δt::Float64; μ::Float64=GM_EARTH)
    # Mean motion
    n = sqrt(μ / oe.a^3)
    
    # Current mean anomaly
    M0 = true_to_mean(oe.ν, oe.e)
    
    # Propagate mean anomaly
    M = mod(M0 + n * Δt, 2π)
    
    # Convert back to true anomaly
    ν_new = mean_to_true(M, oe.e)
    
    return KeplerianElements(oe.a, oe.e, oe.i, oe.Ω, oe.ω, ν_new)
end

"""
    propagate_j2(oe, Δt; μ=GM_EARTH)

Propagate orbit including J2 secular effects.

Includes:
- RAAN precession: Ω̇ = -3/2 * J2 * (R/p)² * n * cos(i)
- Argument of periapsis drift: ω̇ = 3/4 * J2 * (R/p)² * n * (5cos²(i) - 1)
- Mean anomaly rate correction

# Arguments
- `oe`: Initial orbital elements
- `Δt`: Time to propagate [s]

# Returns
- New `KeplerianElements` with J2 secular effects
"""
function propagate_j2(oe::KeplerianElements, Δt::Float64; μ::Float64=GM_EARTH)
    a, e, i, Ω, ω, ν = oe.a, oe.e, oe.i, oe.Ω, oe.ω, oe.ν
    
    # Semi-latus rectum
    p = a * (1 - e^2)
    
    # Mean motion
    n = sqrt(μ / a^3)
    
    # J2 perturbation rates
    J2_factor = 1.5 * J2_EARTH * (R_EARTH / p)^2 * n
    
    # RAAN precession (negative for prograde orbits)
    Ω_dot = -J2_factor * cos(i)
    
    # Argument of periapsis drift
    ω_dot = J2_factor * (2 - 2.5 * sin(i)^2)
    
    # Mean motion correction
    n_corrected = n * (1 + J2_factor * sqrt(1 - e^2) * (1 - 1.5 * sin(i)^2))
    
    # Propagate mean anomaly
    M0 = true_to_mean(ν, e)
    M = mod(M0 + n_corrected * Δt, 2π)
    
    # Update elements
    Ω_new = mod(Ω + Ω_dot * Δt, 2π)
    ω_new = mod(ω + ω_dot * Δt, 2π)
    ν_new = mean_to_true(M, e)
    
    return KeplerianElements(a, e, i, Ω_new, ω_new, ν_new)
end

"""
    propagate_cowell(sv, Δt, dt; μ=GM_EARTH, include_j2=true)

Propagate orbit using Cowell's method (numerical integration).

Uses RK4 integrator with optional J2 perturbation.

# Arguments
- `sv`: Initial state vector
- `Δt`: Total propagation time [s]
- `dt`: Integration time step [s]
- `include_j2`: Include J2 perturbation

# Returns
- `(times, states)`: Arrays of times and state vectors
"""
function propagate_cowell(sv::StateVector, Δt::Float64, dt::Float64; 
                          μ::Float64=GM_EARTH, include_j2::Bool=true)
    
    function acceleration(r::Vector{Float64})
        r_mag = norm(r)
        
        # Two-body acceleration
        a = -μ * r / r_mag^3
        
        # J2 perturbation
        if include_j2
            z² = r[3]^2
            r² = r_mag^2
            factor = 1.5 * J2_EARTH * μ * R_EARTH^2 / r_mag^5
            
            a[1] += factor * r[1] * (5 * z² / r² - 1)
            a[2] += factor * r[2] * (5 * z² / r² - 1)
            a[3] += factor * r[3] * (5 * z² / r² - 3)
        end
        
        return a
    end
    
    # RK4 step function
    function rk4_step(r::Vector{Float64}, v::Vector{Float64}, h::Float64)
        k1_v = acceleration(r)
        k1_r = v
        
        k2_v = acceleration(r + 0.5h * k1_r)
        k2_r = v + 0.5h * k1_v
        
        k3_v = acceleration(r + 0.5h * k2_r)
        k3_r = v + 0.5h * k2_v
        
        k4_v = acceleration(r + h * k3_r)
        k4_r = v + h * k3_v
        
        r_new = r + h/6 * (k1_r + 2k2_r + 2k3_r + k4_r)
        v_new = v + h/6 * (k1_v + 2k2_v + 2k3_v + k4_v)
        
        return r_new, v_new
    end
    
    # Integration
    n_steps = ceil(Int, Δt / dt)
    times = zeros(n_steps + 1)
    states = Vector{StateVector}(undef, n_steps + 1)
    
    r, v = copy(sv.r), copy(sv.v)
    times[1] = 0.0
    states[1] = StateVector(copy(r), copy(v))
    
    for i in 1:n_steps
        r, v = rk4_step(r, v, dt)
        times[i+1] = i * dt
        states[i+1] = StateVector(copy(r), copy(v))
    end
    
    return times, states
end

# =============================================================================
# Orbital Maneuvers
# =============================================================================

"""
    hohmann_transfer(r1, r2; μ=GM_EARTH)

Compute Hohmann transfer parameters between circular orbits.

# Arguments
- `r1`: Initial orbit radius [km]
- `r2`: Final orbit radius [km]

# Returns
- `(Δv1, Δv2, Δv_total, tof)`: Delta-v at each point [km/s] and time of flight [s]
"""
function hohmann_transfer(r1::Float64, r2::Float64; μ::Float64=GM_EARTH)
    # Transfer orbit semi-major axis
    a_t = (r1 + r2) / 2
    
    # Velocities on circular orbits
    v1 = sqrt(μ / r1)
    v2 = sqrt(μ / r2)
    
    # Velocities on transfer orbit
    v_t1 = sqrt(μ * (2/r1 - 1/a_t))
    v_t2 = sqrt(μ * (2/r2 - 1/a_t))
    
    # Delta-v
    Δv1 = abs(v_t1 - v1)
    Δv2 = abs(v2 - v_t2)
    
    # Time of flight (half period of transfer orbit)
    tof = π * sqrt(a_t^3 / μ)
    
    return (Δv1=Δv1, Δv2=Δv2, Δv_total=Δv1+Δv2, tof=tof)
end

"""
    bielliptic_transfer(r1, r2, rb; μ=GM_EARTH)

Compute bi-elliptic transfer parameters.

# Arguments
- `r1`: Initial orbit radius [km]
- `r2`: Final orbit radius [km]  
- `rb`: Intermediate apoapsis radius [km]

# Returns
- Named tuple with delta-v values and time of flight
"""
function bielliptic_transfer(r1::Float64, r2::Float64, rb::Float64; μ::Float64=GM_EARTH)
    # Semi-major axes
    a1 = (r1 + rb) / 2
    a2 = (rb + r2) / 2
    
    # Velocities
    v1 = sqrt(μ / r1)
    v2 = sqrt(μ / r2)
    
    v_t1_peri = sqrt(μ * (2/r1 - 1/a1))
    v_t1_apo = sqrt(μ * (2/rb - 1/a1))
    v_t2_apo = sqrt(μ * (2/rb - 1/a2))
    v_t2_peri = sqrt(μ * (2/r2 - 1/a2))
    
    # Delta-v
    Δv1 = abs(v_t1_peri - v1)
    Δv2 = abs(v_t2_apo - v_t1_apo)
    Δv3 = abs(v2 - v_t2_peri)
    
    # Time of flight
    tof = π * (sqrt(a1^3 / μ) + sqrt(a2^3 / μ))
    
    return (Δv1=Δv1, Δv2=Δv2, Δv3=Δv3, Δv_total=Δv1+Δv2+Δv3, tof=tof)
end

"""
    plane_change(v, Δi)

Compute delta-v for a plane change maneuver.

# Arguments
- `v`: Orbital velocity at maneuver point [km/s]
- `Δi`: Inclination change [rad]

# Returns
- `Δv`: Required delta-v [km/s]
"""
plane_change(v::Float64, Δi::Float64) = 2 * v * sin(Δi / 2)

# =============================================================================
# Lambert Problem
# =============================================================================

"""
    lambert_solve(r1, r2, tof; μ=GM_EARTH, prograde=true)

Solve Lambert's problem using universal variable formulation.

# Arguments
- `r1`: Initial position vector [km]
- `r2`: Final position vector [km]
- `tof`: Time of flight [s]
- `prograde`: True for prograde transfer

# Returns
- `(v1, v2)`: Velocity vectors at r1 and r2 [km/s]
"""
function lambert_solve(r1::Vector{Float64}, r2::Vector{Float64}, tof::Float64;
                       μ::Float64=GM_EARTH, prograde::Bool=true)
    
    r1_mag = norm(r1)
    r2_mag = norm(r2)
    
    # Cross product for direction
    cross_r = cross(r1, r2)
    
    # Transfer angle
    cos_Δν = dot(r1, r2) / (r1_mag * r2_mag)
    
    if prograde
        Δν = cross_r[3] >= 0 ? acos(cos_Δν) : 2π - acos(cos_Δν)
    else
        Δν = cross_r[3] < 0 ? acos(cos_Δν) : 2π - acos(cos_Δν)
    end
    
    # Geometric parameters
    # Handle edge case when cos_Δν is close to 1
    if abs(1 - cos_Δν) < 1e-10
        return ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    end
    A = sin(Δν) * sqrt(r1_mag * r2_mag / (1 - cos_Δν))
    
    # Stumpff functions
    function stumpff_C(z)
        if z > 1e-6
            return (1 - cos(sqrt(z))) / z
        elseif z < -1e-6
            return (1 - cosh(sqrt(-z))) / z
        else
            return 1/2 - z/24 + z^2/720
        end
    end
    
    function stumpff_S(z)
        if z > 1e-6
            sz = sqrt(z)
            return (sz - sin(sz)) / sz^3
        elseif z < -1e-6
            sz = sqrt(-z)
            return (sinh(sz) - sz) / sz^3
        else
            return 1/6 - z/120 + z^2/5040
        end
    end
    
    # Newton iteration for z
    z = 0.0
    for iter in 1:100
        C = stumpff_C(z)
        S = stumpff_S(z)
        
        # Check if C is valid
        if C <= 0
            z += 0.1
            continue
        end
        
        y = r1_mag + r2_mag + A * (z * S - 1) / sqrt(C)
        
        if y < 0
            z += 0.1
            continue
        end
        
        χ = sqrt(y / C)
        
        F = χ^3 * S + A * sqrt(y) - sqrt(μ) * tof
        
        # Compute derivative more carefully
        if abs(C) > 1e-10
            dFdz = χ^3 * (stumpff_S(z + 1e-8) - S) / 1e-8 + A/2 * (C * y + A) / (2 * C * sqrt(y * C))
        else
            dFdz = sqrt(2) / 40 * y^1.5 + A/8 * (sqrt(y) + A * sqrt(1 / (2y)))
        end
        
        # Prevent divergence
        if abs(dFdz) < 1e-15
            break
        end
        
        dz = -F / dFdz
        
        # Limit step size
        if abs(dz) > 1.0
            dz = sign(dz) * 1.0
        end
        
        z += dz
        
        if abs(dz) < 1e-10
            break
        end
    end
    
    # Final computation with safety checks
    C = stumpff_C(z)
    S = stumpff_S(z)
    
    if C <= 0
        # Failed to converge, return zero velocities
        return ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    end
    
    y = r1_mag + r2_mag + A * (z * S - 1) / sqrt(C)
    
    if y < 0
        # Invalid solution
        return ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    end
    
    f = 1 - y / r1_mag
    g = A * sqrt(y / μ)
    g_dot = 1 - y / r2_mag
    
    v1 = (r2 - f * r1) / g
    v2 = (g_dot * r2 - r1) / g
    
    return (v1, v2)
end

# =============================================================================
# Ground Track
# =============================================================================

"""
    compute_ground_track(oe, duration, dt; μ=GM_EARTH, use_j2=true)

Compute ground track for an orbit.

# Arguments
- `oe`: Initial orbital elements
- `duration`: Total time [s]
- `dt`: Time step [s]
- `use_j2`: Include J2 perturbation

# Returns
- Vector of `GroundTrack` points
"""
function compute_ground_track(oe::KeplerianElements, duration::Float64, dt::Float64;
                               μ::Float64=GM_EARTH, use_j2::Bool=true)
    
    n_points = ceil(Int, duration / dt) + 1
    track = Vector{GroundTrack}(undef, n_points)
    
    current_oe = oe
    
    for i in 1:n_points
        t = (i - 1) * dt
        
        # Get ECI position
        sv = keplerian_to_cartesian(current_oe; μ=μ)
        
        # Rotate to ECEF (simplified - just longitude rotation)
        θ = OMEGA_EARTH * t
        R_eci_ecef = [
            cos(θ)  sin(θ)  0;
            -sin(θ) cos(θ)  0;
            0       0       1
        ]
        
        r_ecef = R_eci_ecef * sv.r
        
        # Convert to LLA
        lat, lon, alt = ecef_to_lla(r_ecef)
        
        track[i] = GroundTrack(t, lat, lon, alt)
        
        # Propagate
        if use_j2
            current_oe = propagate_j2(current_oe, dt; μ=μ)
        else
            current_oe = propagate_kepler(current_oe, dt; μ=μ)
        end
    end
    
    return track
end

# =============================================================================
# Pretty Printing
# =============================================================================

function Base.show(io::IO, oe::KeplerianElements)
    println(io, "KeplerianElements:")
    @printf(io, "  a = %.3f km\n", oe.a)
    @printf(io, "  e = %.6f\n", oe.e)
    @printf(io, "  i = %.4f° (%.6f rad)\n", rad2deg(oe.i), oe.i)
    @printf(io, "  Ω = %.4f° (%.6f rad)\n", rad2deg(oe.Ω), oe.Ω)
    @printf(io, "  ω = %.4f° (%.6f rad)\n", rad2deg(oe.ω), oe.ω)
    @printf(io, "  ν = %.4f° (%.6f rad)", rad2deg(oe.ν), oe.ν)
end

function Base.show(io::IO, sv::StateVector)
    println(io, "StateVector:")
    @printf(io, "  r = [%.6f, %.6f, %.6f] km\n", sv.r...)
    @printf(io, "  v = [%.6f, %.6f, %.6f] km/s\n", sv.v...)
    @printf(io, "  |r| = %.6f km\n", norm(sv.r))
    @printf(io, "  |v| = %.6f km/s", norm(sv.v))
end

end # module
