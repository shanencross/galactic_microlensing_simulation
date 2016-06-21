"""
rate_calculation_testing_2.py
@author Shanen Cross
"""

import sys
import os

import numpy
import astropy
from scipy.optimize import fmin

import reading_in_star_population
import time_advancing

STAR_POP_DIR = os.path.join(sys.path[0], "star_population_tables")
#STAR_POP_FILENAME = "1466028123.767236.resu"
#STAR_POP_FILENAME = "1466028463.709599.resu"
#STAR_POP_FILENAME = "1466032757.632040.resu"
STAR_POP_FILENAME = "1466196244.700497.resu"

STAR_POP_FILEPATH = os.path.join(STAR_POP_DIR, STAR_POP_FILENAME)

def main():
    pass
    """
    pseudocode:

    -read in star pop data

    -time advance star pop data for a given period of time with a given time step size
    -from this, get set of advanced star populations for each step from start-time to end-time

    -calculate whether any source stars passed within a corresponding len's einstein ring
    -for each time step:
        -for each pair of lens and source:
            -Get Einstein ring radius
            -get angular separation between source and lens: magnitude of ((l_source, b_source) - (l_lens, b_lens))
            -if angular separation is less than the einstein ring radius, this counts as an event...
            -...but somehow don't count the same event occuring over multiple steps?

    -Or, alternatively:
        -calculate only the start and end points for time advancement of a start_pop: inital and final star_pop
        -have a function that determines both the angular separation and Einstein ring radius between star pops 
         at any given time between start and end point
        -Also have the function deliver the minimum angular separation - einstein ring difference
        -for each pair of lens and source:
            -call angular separation einstein radius minimum difference function and check if it is less than the 0 (if so it is a lensing event)

    -Say for lens star AB and source star XY we have initial positions ab_i = (a_i, b_i) and xy_i = (x_i, y_i) and final positions 
     ab_f = (a_f, b_f) and xy_f = (x_f, y_f)
        -final position calculated by p_ab(u_ab, ab_i, t_f) = (a_f, b_f) = (v_ab, w_ab) * t_f + (a_i, b_i) where (v_ab, w_ab) is AB velocity 
         u_ab and t_f is user-specified time period
        -equivalent done for XY 
    -angular separation at time t is ang_separation(u_ab, ab_i, u_xy, xy_i, t) = p_ab(u_ab, ab_i, t) - p_xy(u_xy, xy_i, t)
    -einstein ring radius at time t is get_einstein_rad(dist_ab(v_r_ab, ab_i, t), dist_xy(v_r_xy, xy_i, t), m_ab)
        -dist_ab(t) ~= dist_ab_i (basically constant?)
        -OR actually we have radial velocity v_r, so we can time advance distance as well:
         dist_ab(v_r_ab, ab_i, t) = v_r * t + ab_i
        -get_einstein_rad(d_ab, d_xy, m_ab) = sqrt( (4*g*m_ab/c**2) * (1/d_ab - 1/d_xy) )
    -difference_function(u_ab, v_r_ab, ab_i, u_xy, v_r_xy, xy_i, t) \
        = ang_separation(u_ab, ab_i, u_xy, xy_i, t) - get_einstein_rad(dist_ab(v_r_ab, ab_i, t), dist_xy(v_r_xy, xy_i, t), m_ab)

    -if min(difference_function(~) <= 0:
        -event_counter += 1
    ...

    Ok, so actually getting the minimization function by solving d(difference_function)/dt = 0 for t and plugging the t value(s) back
    into the original difference function (and picking whichever result is smaller if their are multiple t solutions to avoid getting max)...
    is going to be a pain mathematically, might need mathematica, then we can compare this to computational method later.

    
    

    """

    star_pop = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH)
    print
    #advanced_star_pop = time_advancing.time_advance(star_pop)
    lens_index = 5200
    source_index = 5470
    lens = star_pop[lens_index]
    source = star_pop[source_index]

    microlensing_happened = do_they_microlens(lens, source)
    
    print microlensing_happened    

    #print advanced_star_pop

def get_angular_Einstein_radius(dist_initial_lens, dist_initial_source, v_r_lens, v_r_source, time, mass_lens):
    dist_lens  = get_position_or_dist(v_r_lens, dist_initial_lens, time)
    dist_source = get_position_or_dist(v_r_source, dist_initial_source, time)

    theta_E = get_angular_Einstein_Radius_with_distances(dist_lens, dist_source, mass_lens)
    return theta_E

def get_angular_Einstein_radius_with_distances(dist_lens, dist_source, mass_lens):
    #D_lens = lens["Dist"] * units.kpc
    #mass = lens["Mass"] * units.solMass
    #D_source = source["Dist"] * units.kpc
    theta_E = np.sqrt( ( (4*G*mass_lens)/c**2 ) * ( 1/dist_lens - 1/dist_source ) )
   
    theta_E = (theta_E.decompose() * units.rad).to(units.mas)
    return theta_E

def get_angular_separation(position_initial_lens, position_initial_source, u_lens, u_source, time)
    position_lens = get_position_or_dist(u_lens, position_initial_lens, time)
    position_source = get_position_or_dist(u_source, position_initial_source, time)

    separation = position_source - position_lens
    pass

def get_position_or_dist(velocity, position_or_dist_initial, time):
    position_or_dist = velocity * time + position_or_dist_initial)
    return position_or_dist
    

def do_they_microlens(lens, source):

    dist_initial_lens = lens["Dist"] * units.kpc
    dist_inital_source = source["Dist"] * units.kpc
    mass_lens = lens["Mass"] * units.solMass
    u_lens = (lens["mul"] * units.mas, lens["mub" * units.mas])
    u_source = (source["mul"] * units.mas, source["mub" * units.mas])
    v_r_lens = lens["Vr"] * units.kpc # May not be correct units?
    v_r_source = source["Vr"] * units.kpc # May not be correct units?
    position_initial_lens = (lens["l"] * units.mas, lens["b"] * units.mas)
    position_initial_source = (source["l"] * units.mas, source["b"] * units.mas)

    # Technically not always accurate since radial distance changes with time,
    # But how likely is it practically that one a source will begin in front of a lens,
    # then move behind it, and then behind in a lensing event? At the time-scales we're looking at?
    if dist_initial_source > dist_initial_lense:
        return False

    min_Einstein_ring_separation = \
        fmin(get_angular_separation(position_initial_lens, position_initial_source, u_lens, u_source, time) - \
        get_angular_Einstein_radius(dist_initial_lens, dist_initial_source, v_r_lens, v_r_source, time, mass_lens))
    

    return False



if __name__ == "__main__":
    main()
