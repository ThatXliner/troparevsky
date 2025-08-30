# TODO: Cache data requests to disk
# And also parallelize them
import json
import math
from functools import lru_cache
from itertools import combinations
from linecache import cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import requests

# from optimade.client import OptimadeClient

# TODO: we're re-implementing but using optimade
# OQMD OPTIMADE endpoint

# client = OptimadeClient(include_providers="oqmd")

# client = OptimadeClient(include_providers="oqmd")


def generate_enthalpy_matrix(elements):
    n = len(elements)
    matrix = np.full((n, n), np.nan)  # fill with NaN for masking

    for i in range(n):
        for j in range(i):
            el1 = elements[i]
            el2 = elements[j]
            try:
                value = Hf(el1, el2)
                if value is not None:
                    matrix[i, j] = value
                    # Symmetry
                    matrix[j, i] = value
            except Exception as e:
                print(f"Failed for {el1}-{el2}: {e}")

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 8))
    cmap = plt.cm.get_cmap("coolwarm")  # diverging colormap
    cax = ax.matshow(matrix, cmap=cmap, vmin=-1000, vmax=0)  # MeV range

    # Axis labels
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(elements, rotation=90)
    ax.set_yticklabels(elements)

    # Annotate values
    for (i, j), val in np.ndenumerate(matrix):
        if not np.isnan(val):
            ax.text(
                j,
                i,
                f"{val:.0f}",
                ha="center",
                va="center",
                fontsize=8,
                color="black",
            )

    plt.title("Binary Compound Formation Enthalpy (meV/atom)", pad=20)
    fig.colorbar(cax, label="ΔHf (meV/atom)")
    plt.tight_layout()
    plt.show()


def main():
    # Original cantor alloy
    print(
        can_form_single_phase_HEA(
            [
                ("Co", 0.2),
                ("Cr", 0.2),
                ("Fe", 0.2),
                ("Mn", 0.2),
                ("Ni", 0.2),
            ],
            1000,
        )
    )
    # elements = [
    #     "Ti",
    #     "V",
    #     "Cr",
    #     "Mn",
    #     "Fe",
    #     "Co",
    #     "Ni",
    #     "Cu",
    #     "Zn",
    #     "Nb",
    #     "Mo",
    #     "Ru",
    #     "Rh",
    #     "Pd",
    #     "Ta",
    #     "W",
    #     "Re",
    #     "Os",
    # ]  # example list
    # generate_enthalpy_matrix(elements)


Alloy = list[tuple[str, float]]


def can_form_single_phase_HEA(alloy, annealing_temperature):
    """Annealing temperature is in Kelvin"""
    binary_enthalpies = [
        Hf(element1, element2)
        for (element1, _concentrations) in alloy
        for (element2, _concentrations) in alloy
        if element1 != element2
    ]
    return all(
        get_lower_bound(alloy, annealing_temperature)
        < binary_enthalpy
        < get_upper_bound(alloy)
        for binary_enthalpy in binary_enthalpies
    )


def get_upper_bound(_alloy):
    """
    > The upper limit corresponds to the largest value of H_f for which the
    > alloy does not phase separate due to the immiscibility of
    > any pair of elements.

    In the paper, this was set to 37 meV/atom "to include all known singe-phase alloys
    and it is consistent with the thermodynamic model pre-
    sented by Poletti and Battezzati"
    """
    return 37  # mEv/atom


def get_lower_bound(alloy, annealing_temperature):
    """As defined in the paper"""
    # Calculate the lower bound for the alloy
    return -annealing_temperature * mixing_entropy(alloy)


R = 8.31446261815324  # J/(mol*K)
k_b = 0.08617333262  # meV/K


# NOTE: surely the total number of moles doesn't actually matter here...
def mixing_entropy(alloy):
    # total_moles = sum(amount for _, amount in alloy)
    # https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Thermodynamics/Ideal_Systems/Entropy_of_Mixing
    return (
        #
        # -total_moles *
        # -R
        -k_b
        * math.fsum(
            # mole fraction = atomic fraction
            # because the only species are pure elements
            (concentration) * math.log(concentration)
            for (_, concentration) in alloy
            # (amount / total_moles) * math.log(amount / total_moles)
            # for (_, amount) in alloy
        )
    )


# def lowest_enthalpy_of_formation(alloy):
#     enthalpies = [
#         Hf(element1, element2)
#         for (element1, _concentration) in alloy
#         for (element2, _concentration) in alloy
#         if element1 != element2
#     ]
#     return min(enthalpies)


# XXX: Is this the lowest energy structures of
# binary compounds relative to phase separation into pure elements."?
# It seems to match up
cache_file = Path(__file__).parent.parent / "cache.json"

if not cache_file.exists():
    cache_file.touch()
    cache_file.write_text(json.dumps({}))


@lru_cache(maxsize=None)
def Hf(element1, element2):
    cache = json.loads(cache_file.read_text())
    key = str(tuple(sorted((element1, element2))))
    if key in cache:
        return cache[key]
    # Gets the enthalpy of formation
    # of the lowest energy structure of the given binary compound
    response = requests.get(
        "https://aflowlib.org/API/aflux/?enthalpy_formation_atom(*),species({},{}),nspecies(2),$paging(1,1)".format(
            element1, element2
        )
    )
    print(response, element1, element2)
    print(
        "Hf of",
        element1,
        element2,
        "is",
        response.json()[0]["enthalpy_formation_atom"] * 1000,
    )
    cache_file.write_text(
        json.dumps(
            {
                **cache,
                key: response.json()[0]["enthalpy_formation_atom"] * 1000,
            }
        )
    )
    # turn eV/atom to mEv/atom
    return response.json()[0]["enthalpy_formation_atom"] * 1000


# def Hf(element1, element2):
#     """Calculate the enthalpy of formation of a binary compound.

#         This is defined as "the lowest energy structures of binary compounds relative to phase
#     separation into pure elements."

#     However, there is a remark in the conclusion of the paper which says that a "better approach would be to
#     reference the ΔHf of the binary compounds relative to that of the corresponding binary solid solution"
#     since there are a few exceptions in which the solid solution is itself very stable. This remark is not
#     implemented in this function yet.

#         Args:
#             element1 (str): The first element in the binary compound.
#             element2 (str): The second element in the binary compound.

#         Returns:
#             float: The enthalpy of formation of the binary compound.
#     """
#     client.get('elements HAS "Ag"')
#     entries = mpr.get_entries_in_chemsys([element1, element2])
#     pd = PhaseDiagram(entries)  # I'm going to ignore the type errors here

#     # Filter entries to only those that are *binary* (contain only element1 and element2)
#     binary_entries = [
#         e
#         for e in entries
#         if set(e.composition.elements) == {Element(element1), Element(element2)}
#     ]

#     if not binary_entries:
#         return None

#     # Find entry with lowest formation energy per atom
#     lowest_entry = min(binary_entries, key=lambda e: pd.get_form_energy_per_atom(e))
#     energy_mev = pd.get_form_energy_per_atom(lowest_entry) * 1000
#     return {
#         "entry": lowest_entry,
#         "composition": lowest_entry.composition.formula,
#         "formation_energy_mev_per_atom": energy_mev,
#         "energy_above_hull_mev_per_atom": pd.get_e_above_hull(lowest_entry) * 1000,
#     }


if __name__ == "__main__":
    main()
