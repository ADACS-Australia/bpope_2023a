import numpy as np
import starry

starry.config.lazy = False

def create_star():
    star = starry.Primary(starry.Map(deg=0, udeg=2, amp=1.0), m=1.0, r=1.0, prot=1.0)
    star.map[1] = 0.40
    star.map[2] = 0.26
    return star

def create_planet():
    planet = starry.Secondary(starry.Map(ydeg=5, amp=5e-3), m=0, 
                        r=0.1,
                        porb=1.0,
                        prot=1.0,
                        Omega=30,
                        ecc=0.3,
                        w=30,
                        t0=0,
                         )
    return planet

def create_system(star, planet):
    system = starry.System(star, planet)
    return system

def create_time_array():
    # time = np.linspace(-0.25, 3.25, 10000)
    time = np.linspace(-0.25, 1.00, 4)
    return time

def calculate_system_light_curve(system, time_array):
    flux_system = system.flux(time_array)
    return flux_system

def calculate_sp_light_curve(system, time_array):
    flux_star, flux_planet = system.flux(time_array, total=False)
    return flux_star, flux_planet

def main():
    star = create_star()
    planet = create_planet()
    system = create_system(star, planet)
    times = create_time_array()
    flux_s = calculate_system_light_curve(system, times)
    # fs, fp = calculate_sp_light_curve(system, times)
    # return fs, fp
    return flux_s

if __name__ == "__main__":
    main()

    