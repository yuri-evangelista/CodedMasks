import numpy as np

def filter_source(data, ra, dec, verbose=False):
    mask = np.ones(len(data), dtype=bool)
    mask &= (np.isclose(data["RA"], ra) & np.isclose(data["DEC"], dec))
    filtered = data[mask]
    if verbose:
        print("Selected", len(filtered), "out of", len(data), "events")

    return filtered