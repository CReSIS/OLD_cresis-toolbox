"""Tools for determining the difference between two sets of crossovers."""
from collections import namedtuple
from itertools import groupby, product, chain
from json import load
import math


# The distance between crossovers past which will be reported
ACCEPTABLE_DISTANCE = 2


def distance(cx_A, cx_B):
    """Find the distance between two crossovers."""
    if cx_A is None or cx_B is None:
        return math.inf
    return cx_A["cx_geom"].distance(cx_B["cx_geom"])


Match = namedtuple("Match", "distance cx_pair segment_pair")


def find_closest(cx_segment_pair_A, cx_segment_pair_B, segment_pair):
    """Match the crossovers between the pairs to the closest crossover in the other pair."""

    # Add id's to each crossover
    cx_A = []
    for cx_i, cx in enumerate(cx_segment_pair_A):
        cx_n = cx.copy()
        cx_n["id"] = cx_i
        cx_A.append(cx_n)
    cx_B = []
    for cx_i, cx in enumerate(cx_segment_pair_B):
        cx_n = cx.copy()
        cx_n["id"] = cx_i
        cx_B.append(cx_n)

    # Calculate distance between every pair of crossovers
    distances = []
    for pair in product(cx_A, cx_B):
        distances.append(Match(distance(*pair), pair, segment_pair))

    # Iterate distances by lowest and keep first match of each crossover
    matched_A = set()
    matched_B = set()
    matches = []
    for match in sorted(distances, key=lambda d: d.distance):
        
        # Skip pairings involving cx's we've already matched
        if match.cx_pair[0]["id"] in matched_A:
            continue
        if match.cx_pair[1]["id"] in matched_B:
            continue
        matches.append(match)
        matched_A.add(match.cx_pair[0]["id"])
        matched_B.add(match.cx_pair[1]["id"])

    # Add any crossovers without matches
    if len(matched_A) < len(cx_A):
        for cx in cx_A:
            if cx["id"] in matched_A:
                continue
            matches.append(Match(math.inf, (cx, None), segment_pair))
    if len(matched_B) < len(cx_B):
        for cx in cx_B:
            if cx["id"] in matched_B:
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


def find_differences(cx_set_A, cx_set_B):
    """Compare two sets of crossovers and find the distances between crossovers."""

    cx_set_A = group_segment_pairs(cx_set_A)
    cx_set_B = group_segment_pairs(cx_set_B)
    all_segment_pairs = set(cx_set_A.keys()) | set(cx_set_B.keys())

    segment_pair_comparisons = {}

    for segment_pair in all_segment_pairs:
        if segment_pair in cx_set_A and segment_pair in cx_set_B:
            segment_pair_comparisons[segment_pair] = find_closest(cx_set_A[segment_pair], cx_set_B[segment_pair], segment_pair)
        elif segment_pair in cx_set_A:
            segment_pair_comparisons[segment_pair] = [Match(math.inf, (cx, None), segment_pair) for cx in cx_set_A[segment_pair]]
        elif segment_pair in cx_set_B:
            segment_pair_comparisons[segment_pair] = [Match(math.inf, (None, cx), segment_pair) for cx in cx_set_B[segment_pair]]

    ordered_distances = sorted(chain(*segment_pair_comparisons.values()), key=lambda m: m.distance, reverse=True)

    return ordered_distances


if __name__ == "__main__":
    from plot_crossovers import load_data, TARGETA, TARGETB
    data = load_data()
    cx_distances = find_differences(data[f"{TARGETA} crossovers"], data[f"{TARGETB} crossovers"])

    print("Total crossovers:", len([cx.distance for cx in cx_distances]))
    print("Total Distance (no inf):", sum(cx.distance for cx in cx_distances if cx.distance is not math.inf))
    print("# Missing:", len([cx for cx in cx_distances if cx.distance is math.inf]))
    print(f"Total past {ACCEPTABLE_DISTANCE}m difference:", len([cx.distance for cx in cx_distances if cx.distance > ACCEPTABLE_DISTANCE]))
    print(f"Total distance past {ACCEPTABLE_DISTANCE}m difference:", sum([cx.distance for cx in cx_distances if cx.distance > ACCEPTABLE_DISTANCE and cx.distance is not math.inf]))
    
    total_A_points = 0
    for segment in data[f"{TARGETA} segments"]:
        total_A_points += int(segment["num_points"])
    total_B_points = 0
    for segment in data[f"{TARGETB} segments"]:
        total_B_points += int(segment["num_points"])
    
    print("Number of points reduced to", total_B_points / total_A_points)

    print("Unnacceptable differences:")
    print("{:^10} {:^20} {:^23}".format("Distance", "Angles", "Segments"), sep="|")
    for cx in cx_distances:
        if cx.distance <= ACCEPTABLE_DISTANCE:
            break
        angles = (f"{float(cx.cx_pair[0]['angle']) if cx.cx_pair[0] is not None else 'None':<6.4}",
                  f"{float(cx.cx_pair[1]['angle']) if cx.cx_pair[1] is not None else 'None':<6.4}")
        print(f"{cx.distance:<10.4}", f"{str(angles):<20}", cx.segment_pair, sep="|")
