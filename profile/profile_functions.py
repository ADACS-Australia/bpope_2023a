from starry._core.core import OpsLD

def get_X_value(theta, xo,yo,zo,ro,inc,obl,u,f):
    return OpsLD.eval_X(theta,xo,yo,zo,ro,inc,obl,u,f)