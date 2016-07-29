from astropy import units

def get_observation_time(obs_number, obs_period, night_duration, day_night_duration = 24):
    observation_time = int(obs_number * obs_period / night_duration) * day_night_duration \
                       + ( obs_number * obs_period ) % night_duration

    return observation_time

def get_observation_times(duration=154*units.h, obs_period=17.7*units.h,
                          night_duration=10*units.h, day_night_duration=24*units.h):

    time_list = []
    obs_number = 0
    while True:
        time = get_observation_time(obs_number, obs_period,
                                    night_duration, day_night_duration)

        if time >= duration:
            break
        else:
            time_list.append(time)
            obs_number += 1

        print time, duration

    if isinstance(duration, units.Quantity):
        time_list = units.Quantity(time_list)

    return time_list

def main():
    #time_list = get_observation_times(220, 17.7, 10, 24)
    time_list = get_observation_times(duration=220*units.h)
    print(time_list)

if __name__ == "__main__":
    main()
