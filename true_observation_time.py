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

def get_true_observation_time(time, start_time=0*units.h,
                              night_duration=10*units.h, day_night_duration=24*units.h):
    observation_time = int((time + start_time)/ night_duration) * day_night_duration \
                       + (time + start_time) % night_duration

    return observation_time

def get_true_observation_times(start_time = 0*units.h, duration=154*units.h, period=17.7*units.h,
                          night_duration=10*units.h, day_night_duration=24*units.h):
    # start_time is not working properly
    # actually: start_time is the number of nighttime hours that pass before
    # observing begins, with time 0 being fixed as the start of a night

    if night_duration > day_night_duration:
        print "WARNING: Night duration is longer than day+night parameter."
        print "Returning empty list."
        return []

    time_list = []
    obs_number = 0

    while True:
        time_observing = obs_number * period
        time = get_true_observation_time(time_observing, start_time=start_time,
                                         night_duration=night_duration,
                                         day_night_duration=day_night_duration)

        print time_observing, time
        #if time >= start_time + duration:
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
    period = 17.7*units.h
    time_list = get_true_observation_times(period=period*5, duration=550*units.h, start_time=period*1)
    print(time_list)

if __name__ == "__main__":
    main()
