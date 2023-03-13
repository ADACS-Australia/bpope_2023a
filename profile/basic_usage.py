import starry
import numpy as np
starry.config.lazy = False

def bc_01(ydeg):
    """"default map: calculate flux"""
    map = starry.Map(ydeg=ydeg)
    mf = map.flux()
    return mf

def bc_02(udeg):
    """limb-darkened map: calculate flux"""
    map = starry.Map(udeg=udeg)
    mf = map.flux()
    return mf

if __name__ == "__main__":
    
    # mf = bc_01(5)
    # print(mf)
    mf = bc_02(2)
    print(mf)

