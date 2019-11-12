import numpy as np
from . import icecore 


def alt_sensor(time_model, tas_ann, psl_ann, pr_ann, d18Opr):

    # sensor model
    d18O_ice = icecore.ice_sensor(time_model, d18Opr, pr_ann)
    # diffuse model
    ice_diffused = icecore.ice_archive(d18O_ice, pr_ann, tas_ann, psl_ann)

    pseudo_value = ice_diffused[::-1]
    pseudo_time = time_model

    res = {
        'pseudo_value': pseudo_value,
        'd18O_ice': d18O_ice,
        'ice_diffused': ice_diffused,
    }

    return res
