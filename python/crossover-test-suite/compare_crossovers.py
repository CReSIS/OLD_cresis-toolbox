"""Tools for determining the difference between two sets of crossovers."""
from collections import namedtuple
from itertools import groupby, product, chain
from json import load
import math


# The distance between crossovers past which will be reported
ACCEPTABLE_DISTANCE = 2


def distance(cx_1, cx_2):
    """Find the distance between two crossovers."""
    if cx_1 is None or cx_2 is None:
        return math.inf
    return cx_1["cx_geom"].distance(cx_2["cx_geom"])


Match = namedtuple("Match", "distance cx_pair segment_pair")


def find_closest(cx_segment_pair_1, cx_segment_pair_2, segment_pair):
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
        distances.append(Match(distance(*pair), pair, segment_pair))

    # Iterate distances by lowest and keep first match of each crossover
    matched_1 = set()
    matched_2 = set()
    matches = []
    for match in sorted(distances, key=lambda d: d.distance):
        
        # Skip pairings involving cx's we've already matched
        if match.cx_pair[0]["id"] in matched_1:
            continue
        if match.cx_pair[1]["id"] in matched_2:
            continue
        matches.append(match)
        matched_1.add(match.cx_pair[0]["id"])
        matched_2.add(match.cx_pair[1]["id"])

    # Add any crossovers without matches
    if len(matched_1) < len(cx_1):
        for cx in cx_1:
            if cx["id"] in matched_1:
                continue
            matches.append(Match(math.inf, (cx, None), segment_pair))
    if len(matched_2) < len(cx_2):
        for cx in cx_2:
            if cx["id"] in matched_2:
                continue
            matches.append(Match(math.inf, (None, cx), segment_pair))

    return matches
        

def label_pair(seg1_name, seg2_name):
    """
    Produce a label given a pair of segments 
    such that `label_pair(seg1_name, seg2_name) == label_pair(seg2_name, seg1_name)`.
    """
    return " ".join(sorted((seg1_name, seg2_name)))


def group_segment_pairs(cx_set):
    """Group data from crossover sets into pairs of segments."""
    key = lambda row: label_pair(row["seg1_name"], row["seg2_name"])
    # Iters passed to groupby must be sorted by the key
    groups = groupby(sorted(cx_set, key=key), key=key)
    return {k: list(v) for k, v in groups}


def find_differences(cx_set_1, cx_set_2):
    """Compare two sets of crossovers and find the distances between crossovers."""

    cx_set_1 = group_segment_pairs(cx_set_1)
    cx_set_2 = group_segment_pairs(cx_set_2)
    all_segment_pairs = set(cx_set_1.keys()) | set(cx_set_2.keys())

    segment_pair_comparisons = {}

    for segment_pair in all_segment_pairs:
        if segment_pair in cx_set_1 and segment_pair in cx_set_2:
            segment_pair_comparisons[segment_pair] = find_closest(cx_set_1[segment_pair], cx_set_2[segment_pair], segment_pair)
        elif segment_pair in cx_set_1:
            segment_pair_comparisons[segment_pair] = [Match(math.inf, (cx, None), segment_pair) for cx in cx_set_1[segment_pair]]
        elif segment_pair in cx_set_2:
            segment_pair_comparisons[segment_pair] = [Match(math.inf, (None, cx), segment_pair) for cx in cx_set_2[segment_pair]]

    ordered_distances = sorted(chain(*segment_pair_comparisons.values()), key=lambda m: m.distance, reverse=True)

    return ordered_distances


if __name__ == "__main__":
    from plot_crossovers import load_data, SIMP_RES
    data = load_data()
    # NOTE[REECE]: Always putting 0m (full res) and then simplified for consistency. This order is assumed elsewhere
    cx_distances = find_differences(data["0m crossovers"], data[f"{SIMP_RES} crossovers"])

    print("Total crossovers:", len([cx.distance for cx in cx_distances]))
    print("Total Distance (no inf):", sum(cx.distance for cx in cx_distances if cx.distance is not math.inf))
    print("# Missing:", len([cx for cx in cx_distances if cx.distance is math.inf]))
    print(f"Total past {ACCEPTABLE_DISTANCE}m difference:", len([cx.distance for cx in cx_distances if cx.distance > ACCEPTABLE_DISTANCE]))
    print(f"Total distance past {ACCEPTABLE_DISTANCE}m difference:", sum([cx.distance for cx in cx_distances if cx.distance > ACCEPTABLE_DISTANCE and cx.distance is not math.inf]))
    print("Unnacceptable differences:")
    print("{:^10} {:^20} {:^23}".format("Distance", "Angles", "Segments"), sep="|")
    for cx in cx_distances:
        if cx.distance <= ACCEPTABLE_DISTANCE:
            break
        angles = (f"{float(cx.cx_pair[0]['angle']) if cx.cx_pair[0] is not None else 'None':<6.4}",
                  f"{float(cx.cx_pair[1]['angle']) if cx.cx_pair[1] is not None else 'None':<6.4}")
        print(f"{cx.distance:<10.4}", f"{str(angles):<20}", cx.segment_pair, sep="|")