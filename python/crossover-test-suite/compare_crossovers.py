"""Tools for determining the difference between two sets of crossovers."""
from itertools import groupby, product, zip_longest
from json import load
import math


def distance(cx_1, cx_2):
    """Find the distance between two crossovers."""
    if cx_1 is None or cx_2 is None:
        return math.inf
    return cx_1["cx_geom"].distance(cx_2["cx_geom"])


def find_closest(cx_segment_pair_1, cx_segment_pair_2):
    """Match the crossovers between the pairs to the closest crossover in the other pair."""

    # Add id's to each crossover
    cx_1 = []
    for cx_i, cx in enumerate(cx_segment_pair_1):
        cx_n = cx.copy()
        cx_n["id"] = cx_i
        cx_1.append(cx_n)
    cx_2 = []
    for cx_i, cx in enumerate(cx_segment_pair_2):
        cx_n = cx.copy()
        cx_n["id"] = cx_i
        cx_2.append(cx_n)

    # Calculate distance between every pair of crossovers
    distances = []
    for pair in product(cx_1, cx_2):
        distances.append((distance(*pair), pair))

    # Iterate distances by lowest and keep first match of each crossover
    matched_1 = set()
    matched_2 = set()
    matches = []
    for pair in sorted(distances, key=lambda d: d[0]):
        
        # Skip pairings involving cx's we've already matched
        if pair[1][0]["id"] in matched_1:
            continue
        if pair[1][1]["id"] in matched_2:
            continue
        matches.append(pair)
        matched_1.add(pair[1][0]["id"])
        matched_2.add(pair[1][1]["id"])

    # Add any crossovers without matches
    if len(matched_1) < len(cx_1):
        for cx in cx_1:
            if cx["id"] in matched_1:
                continue
            matches.append((math.inf, cx, None))
    if len(matched_2) < len(cx_2):
        for cx in cx_2:
            if cx["id"] in matched_2:
                continue
            matches.append((math.inf, None, cx))

    return matches
        

def label_pair(seg1_name, seg2_name):
    """
    Produce a label given a pair of segments 
    such that `label_pair(seg1_name, seg2_name) == label_pair(seg2_name, seg1_name)`.
    """
    return " ".join(sorted((seg1_name, seg2_name)))


def group_segment_pairs(cx_set):
    """Group data from crossover sets into pairs of segments."""
    groups = groupby(cx_set, lambda row: label_pair(row["seg1_name"], row["seg2_name"]))
    return {k: list(v) for k, v in groups}


def find_differences(cx_set_1, cx_set_2):
    """Compare two sets of crossovers and find the biggest differences."""
    # TODO[Reece]: Within a pair of segments, order the cx's by distance and find missing ones.
    
    cx_set_1 = group_segment_pairs(cx_set_1)
    cx_set_2 = group_segment_pairs(cx_set_2)

    all_segment_pairs = set(cx_set_1.keys()) | set(cx_set_2.keys())

    for segment_pair in all_segment_pairs:
        if segment_pair in cx_set_1:
            if segment_pair in cx_set_2:
                find_closest(cx_set_1[segment_pair], cx_set_2[segment_pair])

        # TODO[Reece]: Report missing pairs


if __name__ == "__main__":
    from plot_crossovers import load_data
    data = load_data()
    find_differences(data["0m crossovers"], data["1m crossovers"])


# TODO[REECE]: Find crossovers without matches and allow interface to filter for these crossovers
