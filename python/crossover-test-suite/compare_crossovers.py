"""Tools for determining the difference between two sets of crossovers."""
from itertools import groupby, zip_longest
from json import load
import math


MAX_DISTANCE = 2  # Do not match crossovers more than 2 meters apart


def distance(cx_1, cx_2):
    """Find the distance between two crossovers."""
    if cx_1 is None or cx_2 is None:
        return math.inf
    return cx_1["cx_geom"].distance(cx_2["cx_geom"])


def closest_to_point(cx_1, other_crossovers):
    """Retrieve the closest point to cx_1 from other_crossovers"""
    distances = [(distance(cx_1, cx_2), cx_2) for cx_2 in other_crossovers]
    closest = min(distances, key=lambda d: d[0])
    if closest[0] > MAX_DISTANCE:
        return None
    other_crossovers.remove(closest[1])
    return closest


def find_closest(cx_segment_pair_1, cx_segment_pair_2):
    """Find the closest cx in the other pair to each cx in the first pair."""
    matches = []
    other_crossovers = cx_segment_pair_2[:]

    no_matches = []
    last_cx = 0
    for cx_num, cx_1 in enumerate(cx_segment_pair_1):
        last_cx = cx_num
        if not other_crossovers:
            # More points in pair_1 than pair_2
            break
        
        closest_point = closest_to_point(cx_1, other_crossovers)
        if closest_point is None:
            no_matches.append(cx_1)
            continue

        matches.append((closest_point[0], cx_1, closest_point[1]))

    if other_crossovers:
        # More points in pair_2 than pair_1
        for cx_1, cx_2 in zip_longest(no_matches, other_crossovers, fillvalue=None):
            if cx_1 is None or cx_2 is None:
                matches.append((distance(cx_1, cx_2), cx_1, cx_2))
            else:
                # TODO[Reece]: Handle extras
    elif len(matches) < len(cx_segment_pair_1):
        # More points in pair_1 than pair_2
        for cx_1 in cx_segment_pair_1[last_cx:]:
            matches.append((math.inf, cx_1, None))

    print([match[0] for match in matches], len(matches), len(cx_segment_pair_1), len(cx_segment_pair_2))
    print()

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
