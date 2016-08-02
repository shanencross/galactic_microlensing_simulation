"""
true_observation_time.py
@ author Shanen Cross
"""
from astropy import units

"""
def get_true_observation_time(obs_number, period, night_duration, day_night_duration=24*units.h):
    observation_time = int(obs_number * period / night_duration) * day_night_duration \
                       + ( obs_number * period ) % night_duration

    return observation_time
"""

def get_true_observation_time(time, start_time_offset=0*units.h,
                              night_duration=10*units.h, day_night_duration=24*units.h):
    observation_time = int(time / night_duration) * day_night_duration \
                       + time % night_duration

    return observation_time

def get_true_observation_times(duration=154*units.h, period=17.7*units.h,
                          night_duration=10*units.h, day_night_duration=24*units.h):

    if night_duration > day_night_duration:
        print "WARNING: Night duration is longer than day+night parameter."
        print "Returning empty list."
        return []

    time_list = []
    obs_number = 0

    while True:
        time_observing = obs_number * period
        time = get_true_observation_time(time_observing,
                                         night_duration=night_duration,
                                         day_night_duration=day_night_duration)

        print (obs_number * period), time
        if time >= duration:
            break
        else:
            time_list.append(time)
            obs_number += 1

        #print time, duration
    if isinstance(duration, units.Quantity):
        time_list = units.Quantity(time_list)

    return time_list

def main():
    #time_list = get_true_observation_times(220, 17.7, 10, 24)
    #time_list = get_true_observation_times(duration=220*units.h, night_duration = 24*units.h)
    time_list = get_true_observation_times(duration=220*units.h)
    print(time_list)

if __name__ == "__main__":
    main()
