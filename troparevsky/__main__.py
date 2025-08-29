import os
from functools import lru_cache

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
    elements = [
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Nb",
        "Mo",
        "Ru",
        "Rh",
        "Pd",
        "Ta",
        "W",
        "Re",
        "Os",
    ]  # example list
    generate_enthalpy_matrix(elements)


def mixing_entropy(alloy): ...


# XXX: Is this the lowest energy structures of
# binary compounds relative to phase separation into pure elements."?
# It seems to match up
@lru_cache(maxsize=None)
def Hf(element1, element2):
    response = requests.get(
        "https://aflowlib.org/API/aflux/?enthalpy_formation_atom(*),species({},{}),nspecies(2),$paging(1,1)".format(
            element1, element2
        )
    )
    print(
        "Hf of",
        element1,
        element2,
        "is",
        response.json()[0]["enthalpy_formation_atom"] * 1000,
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
