"""Queries to obtain data for plotting."""

import csv
import json
from typing import List, Tuple

from plot_crossovers import DATA_DIR

import psycopg2


def query(query_str, *args, **kwargs):
    """Run a query against the postgres database identified in postgres-conn.json"""

    with open("postgres-conn.json") as f:
        with psycopg2.connect(**json.load(f)[DB]) as conn:
            with conn.cursor() as cursor:
                cursor.execute(query_str, *args, **kwargs)

                q = query_str.strip().upper()
                if q.startswith("UPDATE"):
                    return cursor.rowcount
                else:
                    return cursor.fetchall(), cursor.description


def run_query_all_segments(current_resolution: float, system: str = "rds"):
    """
    Retrieve the segments information for all segments from the database and write
    to a CSV in the DATA_DIR.

    current_resolution : The resolution of the current database for marking the CSV
    system: the system from which to query segments.
    """

    segments, desc = query(f"select name, st_npoints(geom) as num_points, geom from {system}_segments")  # type: ignore
    
    with open(DATA_DIR / f"{current_resolution}m segments.csv", "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow([header[0] for header in desc])
        writer.writerows(segments)


def run_query_all_crossovers(current_resolution: float, system: str = "rds"):
    """
    Retrieve the crossover information for all crossovers from the database and write
    to a CSV in the DATA_DIR.

    current_resolution : The resolution of the current database for marking the CSV
    system: the system from which to query crossovers.
    """

    crossovers, desc = query(f"""select angle, co.geom as cx_geom, seg1.name as seg1_name, 
                           seg2.name as seg2_name, pp1.id as pp1_id, pp2.id as pp2_id, 
                           pp1.geom as pp1_geom, pp2.geom as pp2_geom from {system}_crossovers co
                           join {system}_point_paths pp1 on co.point_path_1_id=pp1.id
                           join {system}_point_paths pp2 on co.point_path_2_id=pp2.id
                           join {system}_segments seg1 on seg1.id=pp1.segment_id
                           join {system}_segments seg2 on seg2.id=pp2.segment_id;""")  # type: ignore
    
    with open(DATA_DIR / f"{current_resolution}m crossovers.csv", "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow([header[0] for header in desc])
        writer.writerows(crossovers)


def run_query_segments(current_resolution: float, system: str = "rds", 
                       segments: Tuple[str, str] = ('20190512_01', '20190516_01')):
    """
    Retrieve the segments information for two provided segments from the database 
    and write to a CSV in the DATA_DIR.

    current_resolution : The resolution of the current database for marking the CSV
    system: the system from which to query segments.
    segments: a tuple of two segments to retrieve information for.
    """

    segments, desc = query(f"""with geom1 as (select geom from {system}_segments where name=%s),
                                geom2 as (select geom from {system}_segments where name=%s)

                              select distinct seg.name, st_npoints(geom) as num_points, ss.name, ss.location_id, geom 
                              from {system}_segments seg join {system}_seasons ss on ss.id=seg.season_id where ss.location_id=1
                              and (seg.name in %s or seg.name in 
                              (select seg.name from geom1, geom2, {system}_segments seg
                              where st_intersects(seg.geom, geom1.geom) or st_intersects(seg.geom, geom2.geom)));""",
                              [segments[0], segments[1], segments])  # type: ignore
    
    with open(DATA_DIR / f"{current_resolution}m segments.csv", "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow([header[0] for header in desc])
        writer.writerows(segments)


def run_query_crossovers(current_resolution: float, system: str = "rds",
                         segments: Tuple[str, str] = ('20190512_01', '20190516_01')):
    """
    Retrieve the crossover information for two provided segments from the database
    and write to a CSV in the DATA_DIR.

    current_resolution : The resolution of the current database for marking the CSV
    system: the system from which to query crossovers.
    segments: a tuple of two segments to retrieve information for.
    """

    crossovers, desc = query(f"""select angle, co.geom as cx_geom, seg1.name as seg1_name, 
                                 seg2.name as seg2_name, pp1.id as pp1_id, pp2.id as pp2_id, 
                                 pp1.geom as pp1_geom, pp2.geom as pp2_geom from {system}_crossovers co
                                 join {system}_point_paths pp1 on co.point_path_1_id=pp1.id
                                 join {system}_point_paths pp2 on co.point_path_2_id=pp2.id
                                 join {system}_segments seg1 on seg1.id=pp1.segment_id
                                 join {system}_segments seg2 on seg2.id=pp2.segment_id
                                 where seg1.name in %s or seg2.name in %s;""",
                                 [segments, segments])  # type: ignore
    
    with open(DATA_DIR / f"{current_resolution}m crossovers.csv", "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow([header[0] for header in desc])
        writer.writerows(crossovers)


def run_delete_crossovers(system: str = "rds"):
    """Remove all crossovers from the database.
    
    system: the system from which to delete crossovers.
    """

    print("Deleting crossovers")
    query(f"delete from {system}_crossovers;")
    print("Setting crossover_calc false")
    query(f"update {system}_segments set crossover_calc=false;")


if __name__ == "__main__":
    RESOLUTION = .1
    DB = "ops0"
    # run_delete_crossovers()
    run_query_segments(RESOLUTION)
    run_query_crossovers(RESOLUTION)
