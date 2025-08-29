import requests


# XXX: Is this the lowest energy structures of
# binary compounds relative to phase separation into pure elements."?
# It seems to match up
def get_enthalpy_of_formation_per_atom(element1, element2):
    response = requests.get(
        "https://aflowlib.org/API/aflux/?enthalpy_formation_atom(*),species({},{}),nspecies(2),$paging(1,1)".format(
            element1, element2
        )
    )
    return response.json()[0]["enthalpy_formation_atom"]
