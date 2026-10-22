'''
flightplan.py – Drone Flight Plan GUI
Generates KMZ (DJI waypoints with turns) + PDF flight plan
Requirement: fpdf2, pyproj (see environment.yml)

Methane measurement – DJI M400 – EU A1/A3 Open Category
Lawnmower pattern with 90-degree turns, configurable
extensions (e.g. default 5 m), configurable hover per line end.

USAGE:
    python flightplan.py

ARCHITECTURE:
    Geometry           – lat/lon ↔ meters (local tangent plane)
    Polygon validation – simple polygon check (no self-intersection)
    Rectangle finding  – bounding box in rotated coords
    Route generation   – lawnmower + 90-degree turns + configurable extensions
    Convex Hull        – Andrew's monotone chain
    Safety check       – A3 building distances (Overpass API)
    Airspace check     – no-fly zones via OSM (Overpass API)
    Weather check      – Open-Meteo API (forecast + archive)
    KMZ export         – DJI waypoint format (turn + hover + speed)
    PDF Generator      – 3-page flight plan (fpdf2)
    GUI                – tkinter (canvas + scrollable content)

COORDINATE INPUT (polygon field, one per line):
    Decimal:   56.493029 9.560055
    DMS (min): 56°29.6055'N 9°33.6107'E
    DMS (sec): 56°29'34.19"N 9°33'45.45"E
    (All formats can be mixed in the same file.)

ROUTE STRUCTURE:
    START (extension m before area)
    [Line 0] → [90° x2 turn] → [Line 1] → ... → [Line N]
      → END (extension m after last line)

RETURN VALUE OF generate_lawnmower:
    (waypoints, wp_actions, line_length, total_distance,
     n_lines, turn_warning)

APIs USED:
    Overpass API      – buildings + airspace (OSM data)
    Open-Meteo        – weather forecast (7 days) + archive

OUTPUT FILES:
    *.kmz             – KMZ with DJI waypoints (DJI Pilot 2)
    *.pdf             – 3-page flight plan (PDF)
'''

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import zipfile, re, math, os, json, time
import urllib.request, urllib.parse
from datetime import datetime
from fpdf import FPDF
from pyproj import Transformer

# ════════════════════════════════════════════
#  Geometry
# ════════════════════════════════════════════

def latlon_to_meters(lat, lon, ref_lat):
    '''
    Converts latitude/longitude to local tangent plane coordinates (meters).

    Uses a simple equirectangular approximation:
      x = lon * 111320 * cos(ref_lat)  (meters east)
      y = lat * 111320                 (meters north)

    Args:
        lat     – latitude in decimal degrees
        lon     – longitude in decimal degrees
        ref_lat – reference latitude (degrees) for the cos correction

    Returns:
        (x, y) in meters relative to the equator/prime meridian.
        In practice only differences matter (used around ref_lat).
    '''
    x = lon * 111320 * math.cos(math.radians(ref_lat))
    y = lat * 111320
    return x, y

def meters_to_latlon(x, y, ref_lat, cx=0, cy=0):
    '''
    Converts local tangent plane coordinates (meters) back to lat/lon.

    Uses the inverse of the equirectangular approximation:
      lon = (x + cx) / (111320 * cos(ref_lat))
      lat = (y + cy) / 111320

    Args:
        x       – meters east (local tangent plane)
        y       – meters north (local tangent plane)
        ref_lat – reference latitude (degrees) for the cos correction
        cx      – x-offset in meters (default 0)
        cy      – y-offset in meters (default 0)

    Returns:
        (lat, lon) in decimal degrees
    '''
    lon = (x + cx) / (111320 * math.cos(math.radians(ref_lat)))
    lat = (y + cy) / 111320
    return lat, lon

# ════════════════════════════════════════════
#  Polygon validation
# ════════════════════════════════════════════

def _segments_intersect(p1, p2, p3, p4):
    '''
    Checks if line segment p1-p2 intersects line segment p3-p4.

    Uses the cross-product (orientation) test: two segments intersect
    if and only if each segment's endpoints are on opposite sides of
    the other segment's line.

    Args:
        p1, p2 – endpoints of segment 1, (x, y) tuples
        p3, p4 – endpoints of segment 2, (x, y) tuples

    Returns:
        True if the segments properly intersect, False otherwise.
        (Collinear/touching cases are not considered intersecting.)
    '''
    def cross(o, a, b):
        '''Cross product of vectors OA and OB.
        Positive = counter-clockwise, negative = clockwise,
        zero = collinear.'''
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    d1 = cross(p3, p4, p1)
    d2 = cross(p3, p4, p2)
    d3 = cross(p1, p2, p3)
    d4 = cross(p1, p2, p4)
    if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and \
       ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)):
        return True
    return False

def is_polygon_simple(polygon_ll):
    '''
    Checks if a polygon is simple (not self-intersecting).

    Converts the polygon to local meters, then tests all non-adjacent
    edge pairs for intersection using _segments_intersect.

    Args:
        polygon_ll – list of (lat, lon) tuples

    Returns:
        True if the polygon is simple (valid, no self-intersection),
        False if any two non-adjacent edges cross.
    '''
    ref_lat = sum(p[0] for p in polygon_ll) / len(polygon_ll)
    poly_m = [latlon_to_meters(lat, lon, ref_lat)
              for lat, lon in polygon_ll]
    n = len(poly_m)
    for i in range(n):
        p1 = poly_m[i]
        p2 = poly_m[(i + 1) % n]
        for j in range(i + 2, n):
            if j == (n - 1) and i == 0:
                continue
            p3 = poly_m[j]
            p4 = poly_m[(j + 1) % n]
            if _segments_intersect(p1, p2, p3, p4):
                return False
    return True

# ════════════════════════════════════════════
#  Largest rectangle INSIDE polygon
# ════════════════════════════════════════════

def _point_in_poly(px, py, poly):
    '''
    Determines if a point is inside a polygon (ray-casting algorithm).

    Casts a horizontal ray from the point to the right and counts how
    many polygon edges it crosses. Odd count = inside, even = outside.

    Args:
        px, py – point coordinates in local meters
        poly   – list of (x, y) tuples in local meters

    Returns:
        True if the point is inside the polygon, False otherwise.
    '''
    inside = False
    n = len(poly)
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if (yi > py) != (yj > py):
            x_int = (xj - xi) * (py - yi) / (yj - yi) + xi
            if px < x_int:
                inside = not inside
        j = i
    return inside

def find_inscribed_rect(poly_m, line_bearing):
    '''
    Finds the bounding box of the rotated polygon.
    Rotates the polygon so the flight bearing becomes the x-axis,
    then returns the bounding box (x_min, y_min, x_max, y_max).
    '''
    theta = math.radians(line_bearing - 90)
    cos_t, sin_t = math.cos(theta), math.sin(theta)
    rotated = [(x * cos_t - y * sin_t,
                x * sin_t + y * cos_t) for x, y in poly_m]
    x_min = min(p[0] for p in rotated)
    x_max = max(p[0] for p in rotated)
    y_min = min(p[1] for p in rotated)
    y_max = max(p[1] for p in rotated)
    return (x_min, y_min, x_max, y_max)

# ════════════════════════════════════════════
#  Route generation (with 90-degree turns)
# ════════════════════════════════════════════

def generate_lawnmower(polygon_ll, line_spacing, line_bearing, altitude,
                       hover_s=0, extension=5):
    """
    Generates a lawnmower (snake) flight route with 90-degree turns.
    Returns:
        (waypoints, wp_actions, line_length, total_distance,
         n_lines, turn_warning)
        wp_actions – list of (turn_ccw, hover_s, wp_type) tuples
    """
    if len(polygon_ll) < 3:
        return [], [], 0, 0, 0, None

    ref_lat = sum(lat for lat, _ in polygon_ll) / len(polygon_ll)
    poly_m = [latlon_to_meters(lat, lon, ref_lat)
              for lat, lon in polygon_ll]
    cx_avg = sum(p[0] for p in poly_m) / len(poly_m)
    cy_avg = sum(p[1] for p in poly_m) / len(poly_m)

    if not is_polygon_simple(polygon_ll):
        hull_m = convex_hull(poly_m)
        poly_c = [(x - cx_avg, y - cy_avg) for x, y in hull_m]
    else:
        poly_c = [(x - cx_avg, y - cy_avg) for x, y in poly_m]

    rect = find_inscribed_rect(poly_c, line_bearing)
    rx1, ry1, rx2, ry2 = rect
    x_start, x_end = min(rx1, rx2), max(rx1, rx2)
    y_start, y_end = min(ry1, ry2), max(ry1, ry2)
    line_length = x_end - x_start

    if line_length < 1e-6:
        return [], [], 0, 0, 0, None

    theta = math.radians(line_bearing - 90)
    cos_t, sin_t = math.cos(theta), math.sin(theta)

    def unrotate(rx, ry):
        x = rx * cos_t + ry * sin_t
        y = -rx * sin_t + ry * cos_t
        return meters_to_latlon(x + cx_avg, y + cy_avg, ref_lat)

    # --- Check if extensions fit inside polygon ---
    turn_warning = None
    n_lines_est = int((y_end - y_start) / line_spacing) + 1
    for li in range(n_lines_est):
        y_line = y_start + li * line_spacing
        if li % 2 == 0:
            extreme_x = x_end + extension
        else:
            extreme_x = x_start - extension
        if not _point_in_poly(extreme_x, y_line, poly_c):
            turn_warning = (
                f"Turn {li+1}: {extension} m extension extends beyond "
                f"the polygon boundary "
                f"(spacing={line_spacing:.1f} m). "
                f"Consider a larger polygon or smaller spacing.")
            break

    # --- Generate waypoints ---
    waypoints = []
    wp_actions = []
    total_distance = 0.0
    n_lines = 0
    y = y_start

    # START: extension m before (x-direction)
    lat_s0, lon_s0 = unrotate(x_start - extension, y_start)
    waypoints.append((lat_s0, lon_s0, altitude))
    wp_actions.append((0, 0, "start"))
    total_distance += extension

    while y <= y_end + line_spacing * 0.5:
        if n_lines % 2 == 0:
            line_start_m = (x_start, y)
            line_end_m = (x_end, y)
        else:
            line_start_m = (x_end, y)
            line_end_m = (x_start, y)

        # --- 90° turns between lines ---
        if n_lines > 0:
            if n_lines % 2 == 1:
                turn_x = x_end + extension
            else:
                turn_x = x_start - extension
            y_prev = y - line_spacing
            turn_ccw = 90 if n_lines % 2 == 1 else 270

            lat_c1, lon_c1 = unrotate(turn_x, y_prev)
            waypoints.append((lat_c1, lon_c1, altitude))
            wp_actions.append((turn_ccw, 0, "turn"))

            lat_c2, lon_c2 = unrotate(turn_x, y)
            waypoints.append((lat_c2, lon_c2, altitude))
            wp_actions.append((turn_ccw, 0, "turn"))
            total_distance += extension + line_spacing + extension

        # --- Line ---
        lat1, lon1 = unrotate(line_start_m[0], line_start_m[1])
        lat2, lon2 = unrotate(line_end_m[0], line_end_m[1])
        waypoints.append((lat1, lon1, altitude))
        wp_actions.append((0, hover_s, "line_start"))
        waypoints.append((lat2, lon2, altitude))
        wp_actions.append((0, hover_s, "line_end"))
        total_distance += line_length
        n_lines += 1
        y += line_spacing

    if n_lines == 0:
        return [], [], 0, 0, 0, None

    # END
    y_last = y_start + (n_lines - 1) * line_spacing
    if n_lines % 2 == 1:
        end_x = x_end + extension
    else:
        end_x = x_start - extension
    lat_e1, lon_e1 = unrotate(end_x, y_last)
    waypoints.append((lat_e1, lon_e1, altitude))
    wp_actions.append((0, hover_s, "end"))
    total_distance += extension

    return (waypoints, wp_actions, line_length, total_distance,
            n_lines, turn_warning)

# ════════════════════════════════════════════
#  Convex Hull
# ════════════════════════════════════════════

def convex_hull(points):
    '''
    Computes the convex hull of a set of 2D points using Andrew's
    monotone chain algorithm.

    Sorts the points, then builds the lower and upper hulls by
    removing points that make non-left turns (cross <= 0).

    Args:
        points – iterable of (x, y) tuples

    Returns:
        List of (x, y) tuples forming the convex hull in
        counter-clockwise order. No duplicate closing point.
    '''
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts
    
    def cross(o, a, b):
        '''Cross product of vectors OA and OB.
        Positive = left turn, negative = right turn, zero = collinear.'''
        return ((a[0] - o[0]) * (b[1] - o[1]) -
                (a[1] - o[1]) * (b[0] - o[0]))
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]

def convex_hull_latlon(polygon_ll):
    '''
    Computes the convex hull of a polygon given in lat/lon.

    Converts to local meters, computes the convex hull, then converts
    back to lat/lon (no offset needed – coordinates are already absolute).

    Args:
        polygon_ll – list of (lat, lon) tuples

    Returns:
        List of (lat, lon) tuples forming the convex hull.
        Returns the input unchanged if fewer than 3 points.
    '''
    if len(polygon_ll) < 3:
        return list(polygon_ll)
    ref_lat = sum(lat for lat, _ in polygon_ll) / len(polygon_ll)
    poly_m = [latlon_to_meters(lat, lon, ref_lat)
              for lat, lon in polygon_ll]
    hull_m = convex_hull(poly_m)
    return [meters_to_latlon(x, y, ref_lat) for x, y in hull_m]

# ════════════════════════════════════════════
#  Rectangle from center
# ════════════════════════════════════════════

def rect_from_center(center_lat, center_lon, bearing_deg, dist1_m, dist2_m):
    '''
    Generates a rectangle defined by its center, direction, and half-lengths.

    The rectangle has two perpendicular directions:
      - Direction 1 = bearing_deg (from north, clockwise)
      - Direction 2 = bearing_deg + 90

    The four corners are generated by combining the two directions
    with +/- signs:
      (1, 1)  – forward along both
      (-1, 1) – back along dir 1, forward along dir 2
      (-1,-1) – back along both
      (1,-1)  – forward along dir 1, back along dir 2

    Args:
        center_lat  – center latitude in decimal degrees
        center_lon  – center longitude in decimal degrees
        bearing_deg – direction 1 in degrees from north (0 = north,
                      90 = east, clockwise)
        dist1_m     – half length along direction 1 in meters
        dist2_m     – half length along direction 2 in meters

    Returns:
        List of 4 (lat, lon) tuples – the rectangle corners in
        counter-clockwise order.
    '''
    b1 = math.radians(bearing_deg)
    b2 = math.radians(bearing_deg + 90)
    u1 = (math.sin(b1), math.cos(b1))
    u2 = (math.sin(b2), math.cos(b2))
    mplon = 111320 * math.cos(math.radians(center_lat))
    mpdat = 111320
    corners = []
    for s1, s2 in [(1,1), (-1,1), (-1,-1), (1,-1)]:
        dlat = (s1 * dist1_m * u1[1] + s2 * dist2_m * u2[1]) / mpdat
        dlon = (s1 * dist1_m * u1[0] + s2 * dist2_m * u2[0]) / mplon
        corners.append((center_lat + dlat, center_lon + dlon))
    return corners

# ════════════════════════════════════════════
#  Safety check (A3 distances)
# ════════════════════════════════════════════

def _point_to_polygon_distance(px, py, polygon):
    '''
    Computes the minimum distance from a point to a polygon's boundary.

    For each edge of the polygon, projects the point onto the edge's
    line (clamped to the segment) and computes the distance to the
    closest point on that edge. Returns the minimum over all edges.

    Args:
        px, py  – point coordinates in local meters
        polygon – list of (x, y) tuples in local meters

    Returns:
        Minimum distance in meters from the point to the polygon
        boundary. (Positive value regardless of inside/outside.)
    '''
    min_dist = float('inf')
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        dx, dy = x2 - x1, y2 - y1
        len2 = dx * dx + dy * dy
        if len2 < 1e-10:
            dist = math.sqrt((px - x1) ** 2 + (py - y1) ** 2)
        else:
            t = max(0, min(1, ((px - x1) * dx + (py - y1) * dy) / len2))
            cx, cy = x1 + t * dx, y1 + t * dy
            dist = math.sqrt((px - cx) ** 2 + (py - cy) ** 2)
        min_dist = min(min_dist, dist)
    return min_dist

def check_buildings_a3(polygon_ll, threshold=150):
    '''
    Checks building distances around a polygon for A3 drone rules.

    Queries the Overpass API (OpenStreetMap) for all buildings and
    building nodes within a 250 m buffer of the polygon. Computes the
    minimum distance from each building's centroid to the polygon
    boundary and compares it against the threshold (A3 requires
    150 m from buildings, 50 m from people).

    Args:
        polygon_ll – list of (lat, lon) tuples defining the flight area
        threshold  – minimum required distance in meters (default 150)

    Returns:
        dict with keys:
            min_distance    – closest building distance in meters,
                              or None if no buildings found
            ok              – True if min_distance >= threshold
                              (or no buildings found)
            building_name   – name of the nearest building,
                              or None, or error message on network failure
            buildings_found – total number of buildings/nodes found
    '''
    if len(polygon_ll) < 3:
        return {"min_distance": None, "ok": False,
                "building_name": None, "buildings_found": 0}
    lats = [lat for lat, _ in polygon_ll]
    lons = [lon for _, lon in polygon_ll]
    lat_center = sum(lats) / len(lats)
    lat_buf = 250 / 111320
    lon_buf = 250 / (111320 * math.cos(math.radians(lat_center)))
    south = min(lats) - lat_buf
    north = max(lats) + lat_buf
    west = min(lons) - lon_buf
    east = max(lons) + lon_buf
    query = f"""
    [out:json][timeout:10];
    way({south},{west},{north},{east})["building"];
    out center tags;
    node({south},{west},{north},{east})["building"];
    out center tags;
    """
    url = "https://overpass-api.de/api/interpreter"
    data = urllib.parse.urlencode({"data": query}).encode("utf-8")
    try:
        req = urllib.request.Request(
            url, data=data,
            headers={"User-Agent": "FlightPlan-SafetyCheck/1.0"})
        with urllib.request.urlopen(req, timeout=20) as resp:
            result = json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        return {"min_distance": None, "ok": False,
                "building_name": f"Network error: {e}",
                "buildings_found": 0}
    elements = result.get("elements", [])
    if not elements:
        return {"min_distance": None, "ok": True,
                "building_name": None, "buildings_found": 0}
    poly_m = [latlon_to_meters(lat, lon, lat_center)
              for lat, lon in polygon_ll]
    min_dist = float('inf')
    nearest_name = None
    for el in elements:
        lat = el.get("lat") or el.get("center", {}).get("lat")
        lon = el.get("lon") or el.get("center", {}).get("lon")
        if lat is None or lon is None:
            continue
        bx, by = latlon_to_meters(lat, lon, lat_center)
        dist = _point_to_polygon_distance(bx, by, poly_m)
        if dist < min_dist:
            min_dist = dist
            tags = el.get("tags", {})
            nearest_name = (tags.get("name") or
                            tags.get("building") or "Unknown building")
    if min_dist == float('inf'):
        return {"min_distance": None, "ok": True,
                "building_name": None, "buildings_found": 0}
    return {"min_distance": min_dist,
            "ok": min_dist >= threshold,
            "building_name": nearest_name,
            "buildings_found": len(elements)}

# ════════════════════════════════════════════
#  Airspace check (OSM/Overpass)
# ════════════════════════════════════════════

def check_airspace_osm(polygon_ll, buffer_m=200):
    '''
    Checks airspace restrictions around a polygon using OSM/Overpass API.

    Queries the Overpass API for restricted zones (no-fly, military,
    aerodromes, heliports, prisons, hospitals) within a buffer of the
    polygon. Classifies each zone by severity and computes its distance
    to the polygon boundary.

    Severity levels:
      PROHIBITED – no_fly drone/yes: do not fly
      PERMISSION – military, aerodrome, heliport, prison:
                   requires explicit permission
      CAUTION    – hospital, other: verify manually

    Distance thresholds:
      < 200 m – PROHIBITED/PERMISSION flagged as problem
      < 500 m – CAUTION flagged as "check manually"

    Args:
        polygon_ll – list of (lat, lon) tuples defining the flight area
        buffer_m   – search buffer in meters (default 200)

    Returns:
        dict with keys:
            ok          – True if no PROHIBITED problems found
            error       – None, or error message on network failure
            problems    – list of problem description strings
            zones_found – list of dicts with name, type, distance_m
    '''
    if len(polygon_ll) < 3:
        return {"ok": False, "error": "Polygon too small",
                "problems": []}
    lats = [p[0] for p in polygon_ll]
    lons = [p[1] for p in polygon_ll]
    ref_lat = sum(lats) / len(lats)
    lat_buf = buffer_m / 111320
    lon_buf = buffer_m / (111320 * math.cos(math.radians(ref_lat)))
    south = min(lats) - lat_buf
    north = max(lats) + lat_buf
    west = min(lons) - lon_buf
    east = max(lons) + lon_buf
    query = f"""
    [out:json][timeout:15];
    (
    nwr({south},{west},{north},{east})["no_fly"="drone"];
    nwr({south},{west},{north},{east})["no_fly"="yes"];
    nwr({south},{west},{north},{east})["military"="yes"];
    nwr({south},{west},{north},{east})["boundary"="military"];
    nwr({south},{west},{north},{east})["aeroway"="aerodrome"];
    nwr({south},{west},{north},{east})["aeroway"="heliport"];
    nwr({south},{west},{north},{east})["building"="prison"];
      node({south},{west},{north},{east})["building"="hospital"];
      node({south},{west},{north},{east})["amenity"="prison"];
    );
    out center tags;
    """
    url = "https://overpass-api.de/api/interpreter"
    data = urllib.parse.urlencode({"data": query}).encode("utf-8")
    try:
        req = urllib.request.Request(
            url, data=data,
            headers={"User-Agent": "FlightPlan-AirspaceCheck/1.0"})
        with urllib.request.urlopen(req, timeout=20) as resp:
            result = json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        return {"ok": False, "error": f"Overpass API error: {e}",
                "problems": []}
    elements = result.get("elements", [])
    if not elements:
        return {"ok": True, "problems": [], "zones_found": []}
    poly_m = [latlon_to_meters(lat, lon, ref_lat)
              for lat, lon in polygon_ll]
    px_min = min(p[0] for p in poly_m)
    px_max = max(p[0] for p in poly_m)
    py_min = min(p[1] for p in poly_m)
    py_max = max(p[1] for p in poly_m)
    problems = []
    zones_found = []
    seen = set()
    for el in elements:
        lat = (el.get("lat") or el.get("center", {}).get("lat"))
        lon = (el.get("lon") or el.get("center", {}).get("lon"))
        if lat is None or lon is None:
            continue
        bm = latlon_to_meters(lat, lon, ref_lat)
        if (bm[0] < px_min - 500 or bm[0] > px_max + 500 or
                bm[1] < py_min - 500 or bm[1] > py_max + 500):
            continue
        tags = el.get("tags", {})
        zname = tags.get("name", "Unknown")

        # --- Determine zone type first ---
        if "no_fly" in tags:
            ztype = "Drone prohibition"
            severity = "PROHIBITED"
        elif "military" in tags or tags.get("boundary") == "military":
            ztype = "Military zone"
            severity = "PERMISSION"
        elif tags.get("aeroway") == "aerodrome":
            ztype = "Airport/CTR"
            severity = "PERMISSION"
        elif tags.get("aeroway") == "heliport":
            ztype = "Heliport"
            severity = "PERMISSION"
        elif tags.get("building") == "prison" or \
                tags.get("amenity") == "prison":
            ztype = "Prison"
            severity = "PERMISSION"
        elif tags.get("building") == "hospital":
            ztype = "Hospital"
            severity = "CAUTION"
        else:
            ztype = "Airspace zone"
            severity = "CAUTION"

        # --- Dedup: skip if same (type, name) already seen ---
        key = f"{ztype}:{zname}"
        if key in seen:
            continue
        seen.add(key)

        # --- Calculate distance and record ---
        dist = _point_to_polygon_distance(bm[0], bm[1], poly_m)
        zones_found.append({"name": zname, "type": ztype,
                            "distance_m": dist})
        if severity == "PROHIBITED" and dist < 200:
            problems.append(
                f"PROHIBITED: {zname} ({ztype}) - "
                f"{dist:.0f} m from polygon. DO NOT FLY")
        elif severity == "PERMISSION" and dist < 200:
            problems.append(
                f"{zname} ({ztype}) - {dist:.0f} m. "
                f"PERMISSION required")
        elif dist < 500:
            problems.append(
                f"{zname} ({ztype}) - {dist:.0f} m. Check manually")
    return {
        "ok": not any(p.startswith("PROHIBITED") for p in problems),
        "problems": problems,
        "zones_found": zones_found,
        "error": None,
    }

# ════════════════════════════════════════════
#  Weather (Open-Meteo)
# ════════════════════════════════════════════

WEATHER_CODES_FOG = {45, 48}
WEATHER_CODES_RAIN = {51, 53, 55, 56, 57, 61, 63, 65, 66, 67, 80, 81, 82}
WEATHER_CODES_DESC = {
    0: "Clear sky", 1: "Mainly clear", 2: "Partly cloudy", 3: "Overcast",
    45: "FOG", 48: "FOG (rime)",
    51: "Light drizzle", 53: "Drizzle", 55: "Heavy drizzle",
    61: "Light rain", 63: "Rain", 65: "Heavy rain",
    66: "Freezing rain", 67: "Heavy freezing rain",
    80: "Rain showers", 81: "Rain showers", 82: "Violent rain showers",
    95: "Thunderstorm", 96: "Thunderstorm + hail",
    99: "Thunderstorm + hail (heavy)",
}

def check_weather_range(lat, lon, date_str=None, start_time="12:00",
                        end_time=None, interval_min=10):
    '''
    Checks weather conditions for a flight time window using the
    Open-Meteo API.

    For each checkpoint in the time window, retrieves hourly weather
    data and checks against drone flight requirements:
      - Wind speed: 2-5 m/s (outside range = problem)
      - Rain: < 0.05 mm/h (above = problem)
      - Fog: any fog code = problem

    Uses the archive API for past dates, forecast API for today
    and future dates (up to 7 days ahead).

    Date formats accepted: DD.MM.YYYY, DD/MM/YYYY, YYYY-MM-DD
    Time format: HH:MM (24-hour)

    Args:
        lat          – latitude in decimal degrees
        lon          – longitude in decimal degrees
        date_str     – date string (default: today)
        start_time   – start time "HH:MM" (default: "12:00")
        end_time     – end time "HH:MM" (default: start + 60 min)
        interval_min – checkpoint interval in minutes (default: 10)

    Returns:
        dict with keys:
            ok         – True if no problems found
            points     – list of weather data dicts (one per checkpoint)
            problems   – list of problem description strings
            date       – effective date (YYYY-MM-DD)
            start_time – start time string
            end_time   – end time string
            error      – None, or error message on failure

        Each dict in "points" contains:
            time            – "HH:MM"
            wind_m_s        – wind speed in m/s
            wind_dir_deg    – wind direction in degrees
            wind_gusts      – gusts in m/s
            rain_mm_h       – precipitation in mm/h
            fog             – True if fog detected
            rain            – True if rain detected
            temp_c          – temperature in °C
            rh_pct          – relative humidity in %
            cloud_cover_pct – cloud cover in %
            description     – human-readable weather code description
    '''
    today = datetime.now().strftime("%Y-%m-%d")

    # Parse user date (accepts DD.MM.YYYY, DD/MM/YYYY, YYYY-MM-DD)
    if date_str:
        dt_parsed = None
        for fmt in ("%d.%m.%Y", "%d/%m/%Y", "%Y-%m-%d"):
            try:
                dt_parsed = datetime.strptime(date_str.strip(), fmt)
                break
            except ValueError:
                continue
        if dt_parsed is None:
            return {"ok": False,
                    "error": "Invalid date. Use format: DD.MM.YYYY "
                             "(e.g. 04.10.2026)",
                    "points": [], "problems": [],
                    "date": date_str,
                    "start_time": start_time, "end_time": ""}
        effective_date = dt_parsed.strftime("%Y-%m-%d")
    else:
        effective_date = today

    try:
        sh, sm = map(int, start_time.split(":"))
        start_min = sh * 60 + sm
    except Exception:
        return {"ok": False, "error": "Invalid start time",
                "points": [], "problems": [],
                "date": effective_date,
                "start_time": start_time, "end_time": ""}
    if end_time:
        try:
            eh, em = map(int, end_time.split(":"))
            end_min = eh * 60 + em
        except Exception:
            return {"ok": False, "error": "Invalid end time",
                    "points": [], "problems": [],
                    "date": effective_date,
                    "start_time": start_time, "end_time": ""}
    else:
        end_min = start_min + 60
    if end_min < start_min:
        return {"ok": False, "error": "End time before start time",
                "points": [], "problems": [],
                "date": effective_date,
                "start_time": start_time, "end_time": ""}

    # Check if the date is too far in the future
    dt_today = datetime.strptime(today, "%Y-%m-%d")
    dt_eff = datetime.strptime(effective_date, "%Y-%m-%d")
    days_ahead = (dt_eff - dt_today).days
    if days_ahead > 7:
        return {"ok": False,
                "error": f"The date {effective_date} is {days_ahead} days in the future. "
                         f"Forecast is only available for the next 7 days. "
                         f"Please use a date within the week.",
                "points": [], "problems": [],
                "date": effective_date,
                "start_time": start_time, "end_time": ""}

    use_archive = (effective_date < today)
    if use_archive:
        url = (
            "https://archive-api.open-meteo.com/v1/archive"
            f"?latitude={lat}&longitude={lon}"
            f"&start_date={effective_date}&end_date={effective_date}"
            f"&hourly=temperature_2m,precipitation,rain,"
            f"weather_code,cloud_cover,wind_speed_10m,"
            f"wind_direction_10m,wind_gusts_10m,"
            f"relative_humidity_2m"
            f"&wind_speed_unit=ms&timezone=auto")
    else:
        url = (
            "https://api.open-meteo.com/v1/forecast"
            f"?latitude={lat}&longitude={lon}"
            f"&hourly=temperature_2m,precipitation,rain,"
            f"weather_code,cloud_cover,wind_speed_10m,"
            f"wind_direction_10m,wind_gusts_10m,"
            f"relative_humidity_2m"
            f"&wind_speed_unit=ms&forecast_days=7&timezone=auto")
    try:
        req = urllib.request.Request(
            url, headers={"User-Agent": "FlightPlan-Weather/1.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        return {"ok": False, "error": str(e),
                "points": [], "problems": [],
                "date": effective_date,
                "start_time": start_time,
                "end_time": f"{end_min // 60:02d}:{end_min % 60:02d}"}
    hourly = data.get("hourly", {})
    times = hourly.get("time", [])
    if not times:
        return {"ok": False, "error": "No weather data found",
                "points": [], "problems": [],
                "date": effective_date,
                "start_time": start_time,
                "end_time": f"{end_min // 60:02d}:{end_min % 60:02d}"}
    n_times = len(times)
    checkpoints = list(range(start_min, end_min + 1,
                              max(1, interval_min)))
    if not checkpoints:
        checkpoints = [start_min]

    def _parse_api_time(t):
        for fmt in ("%Y-%m-%dT%H:%M:%S",
                    "%Y-%m-%dT%H:%M",
                    "%Y-%m-%d %H:%M:%S",
                    "%Y-%m-%d %H:%M"):
            try:
                return datetime.strptime(t[:19], fmt)
            except ValueError:
                continue
        return None

    points = []
    problems = []
    for cp_min in checkpoints:
        cp_h, cp_m = divmod(cp_min, 60)
        target = f"{effective_date}T{cp_h:02d}:{cp_m:02d}"
        dt_target = datetime.strptime(target, "%Y-%m-%dT%H:%M")
        best_idx = 0
        best_diff = float('inf')
        for i, t in enumerate(times):
            dt_t = _parse_api_time(t)
            if dt_t is None:
                continue
            diff = abs((dt_t - dt_target).total_seconds())
            if diff < best_diff:
                best_diff = diff
                best_idx = i
        code = hourly.get("weather_code", [0] * n_times)[best_idx]
        p = {
            "time": f"{cp_h:02d}:{cp_m:02d}",
            "wind_m_s": hourly.get("wind_speed_10m",
                                   [0] * n_times)[best_idx],
            "wind_dir_deg": hourly.get("wind_direction_10m",
                                       [0] * n_times)[best_idx],
            "wind_gusts": hourly.get("wind_gusts_10m",
                                     [0] * n_times)[best_idx],
            "rain_mm_h": (hourly.get("rain", [0] * n_times)[best_idx]
                          or hourly.get("precipitation",
                                        [0] * n_times)[best_idx]),
            "fog": code in WEATHER_CODES_FOG,
            "rain": code in WEATHER_CODES_RAIN,
            "temp_c": hourly.get("temperature_2m",
                                 [0] * n_times)[best_idx],
            "rh_pct": hourly.get("relative_humidity_2m",
                                 [0] * n_times)[best_idx],
            "cloud_cover_pct": hourly.get("cloud_cover",
                                          [0] * n_times)[best_idx],
            "description": WEATHER_CODES_DESC.get(
                code, f"Code {code}"),
        }
        points.append(p)
        if p["wind_m_s"] < 2:
            problems.append(
                f"{p['time']}: Wind too low "
                f"({p['wind_m_s']:.1f} m/s)")
        elif p["wind_m_s"] > 5:
            problems.append(
                f"{p['time']}: Wind too high "
                f"({p['wind_m_s']:.1f} m/s)")
        if p["rain_mm_h"] > 0.05:
            problems.append(
                f"{p['time']}: Rain ({p['rain_mm_h']:.1f} mm/h)")
        if p["fog"]:
            problems.append(f"{p['time']}: FOG")
    return {
        "ok": len(problems) == 0,
        "points": points,
        "problems": problems,
        "date": effective_date,
        "start_time": start_time,
        "end_time": f"{end_min // 60:02d}:{end_min % 60:02d}",
        "error": None,
    }

# ════════════════════════════════════════════
#  PDF Generator
# ════════════════════════════════════════════

def _find_font(basename):
    '''
    Finds a TTF font file on the system across Windows, macOS, and Linux.

    Searches common font directories in order of likelihood.
    Used by FlightPlanPDF to locate Arial for the PDF output.

    Args:
        basename – font filename (e.g. "arial.ttf")

    Returns:
        Full path to the font file, or None if not found.
        (If None, the PDF falls back to the built-in Helvetica.)
    '''
    candidates = [
        os.path.join(os.environ.get("WINDIR", r"C:\Windows"),
                     "Fonts", basename),
        os.path.join(os.path.expanduser("~"), "Fonts", basename),
        # macOS
        "/System/Library/Fonts/" + basename,
        "/Library/Fonts/" + basename,
        os.path.join(os.path.expanduser("~"), "Library",
                     "Fonts", basename),
        # Linux
        "/usr/share/fonts/truetype/msttcorefonts/" + basename,
        "/usr/share/fonts/truetype/" + basename,
        "/usr/share/fonts/" + basename,
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None

class FlightPlanPDF(FPDF):
    '''
    PDF generator for the flight plan document (3 pages).
    Uses fpdf2 to produce a structured flight plan with:
      - Page 1: Flight details, mission, weather, battery, background
      - Page 2: Pre-flight, in-flight, post-flight, abort criteria
      - Page 3: Turn time, key figures
    Font: Arial if found on system (via _find_font),
          otherwise falls back to built-in Helvetica.
    '''
    def __init__(self):
        '''
        Initializes the PDF: sets up font, first page,
        and auto page-break (margin 15 mm).
        '''
        super().__init__()
        font_dir = _find_font("arial.ttf")
        if font_dir:
            d = os.path.dirname(font_dir)
            self.add_font("Arial", "",
                          os.path.join(d, "arial.ttf"))
            self.add_font("Arial", "B",
                          os.path.join(d, "arialbd.ttf"))
            self.add_font("Arial", "I",
                          os.path.join(d, "ariali.ttf"))
            self.add_font("Arial", "BI",
                          os.path.join(d, "arialbi.ttf"))
            self.FONT = "Arial"
        else:
            self.FONT = "Helvetica"
        self.add_page()
        self.set_auto_page_break(auto=True, margin=15)

    def header(self):
        '''
        Draws the page header on pages 2+ (skips page 1).
        Shows title, page number, and a horizontal separator line.
        '''
        if self.page_no() > 1:
            self.set_font(self.FONT, "I", 8)
            self.set_text_color(100, 100, 100)
            self.cell(0, 5,
                      "Flight Plan - Methane - DJI M400",
                      align="L")
            self.cell(0, 5, f"Page {self.page_no()}",
                      align="R", new_x="LMARGIN",
                      new_y="NEXT")
            self.set_draw_color(100, 100, 100)
            self.set_line_width(0.1)
            self.line(self.l_margin, 15,
                      self.w - self.r_margin, 15)
            self.ln(4)

    def _section_box(self, title, color):
        '''
        Draws a colored section header box.

        Args:
            title – section title text
            color – (R, G, B) tuple for fill and border color
        '''
        self.set_fill_color(*color)
        self.set_draw_color(*color)
        self.set_line_width(0.6)
        self.set_font(self.FONT, "B", 11)
        self.set_text_color(*color)
        h = 8
        self.rect(self.l_margin, self.get_y(),
                  self.w - self.l_margin - self.r_margin, h)
        self.set_xy(self.l_margin + 2, self.get_y() + 1)
        self.cell(0, 6, title)
        self.set_y(self.get_y() + h + 3)

    def _bullet(self, text, bold_prefix=""):
        '''
        Draws a bullet-list item with optional bold prefix.

        The text is indented 6 mm from the left margin.
        If bold_prefix is given, it is printed in bold before
        the rest of the text (e.g. "Source: " + value).
        Text wraps automatically to fit the page width.

        Args:
            text        – bullet text (wraps)
            bold_prefix – optional bold text before the bullet
        '''
        block_x = self.l_margin + 6
        self.set_x(block_x)
        if bold_prefix:
            self.set_font(self.FONT, "B", 10)
            self.cell(self.get_string_width(bold_prefix),
                      5, bold_prefix)
            self.set_font(self.FONT, "", 10)
        else:
            self.set_font(self.FONT, "", 10)
        text_x = self.get_x()
        avail = self.w - text_x - self.r_margin
        self.multi_cell(avail, 5, text,
                        new_x="LMARGIN", new_y="NEXT")
        self.ln(1)

    def _item_row(self, label, value):
        '''
        Draws a label-value row (used in Flight Details section).

        The label is bold and left-aligned in a 60 mm column.
        The value wraps in the remaining width.

        Args:
            label – field name (e.g. "Date:")
            value – field value (e.g. "04.10.2026")
        '''
        self.set_font(self.FONT, "B", 10)
        self.cell(60, 5, "  " + label)
        self.set_font(self.FONT, "", 10)
        self.multi_cell(0, 5, value,
                        new_x="LMARGIN", new_y="NEXT")
        self.ln(0.5)

def _weather_table(pdf, fields):
    '''
    Draws the weather requirements table on the PDF.
    Shows rows comparing flight requirements against
    measured values, with an OK/NO status column:
      - Wind: 2-5 m/s, stable
      - Rain/fog: no
      - Temperature (info only)
      - Relative humidity (info only)
    Args:
        pdf    – FlightPlanPDF instance
        fields – dict with weather measurement values (strings)
    '''
    x = pdf.l_margin + 2
    w = pdf.w - pdf.l_margin - pdf.r_margin - 4
    col1, col2, col3 = w * 0.52, w * 0.28, w * 0.20
    row_h = 7
    pdf.set_fill_color(240, 240, 240)
    pdf.set_draw_color(0, 0, 0)
    pdf.set_line_width(0.3)
    y = pdf.get_y()
    pdf.set_font(pdf.FONT, "B", 10)
    pdf.rect(x, y, w, row_h, style="DF")
    pdf.set_xy(x, y)
    pdf.cell(col1, row_h, "  Requirement")
    pdf.cell(col2, row_h, "  Measured value")
    pdf.cell(col3, row_h, "  OK?",
             new_x="LMARGIN", new_y="NEXT")
    pdf.set_font(pdf.FONT, "", 10)
    for req, val, ok in [
        ("Wind 2-5 m/s, stable",
         fields["wind_meas"], fields["wind_ok"]),
        ("Rain/fog: no",
         fields["rain"], fields["rain_ok"]),
        ("Temp:",
         f"{fields['temp_c']}\u00b0C", "-"),
        ("RH:",
         f"{fields['rh_pct']}%", "-"),
    ]:
        y = pdf.get_y()
        pdf.rect(x, y, w, row_h, style="D")
        pdf.set_xy(x, y)
        pdf.cell(col1, row_h, "  " + req)
        pdf.cell(col2, row_h, "  " + val)
        pdf.cell(col3, row_h, "  " + ok,
                 new_x="LMARGIN", new_y="NEXT")
    pdf.line(x, pdf.get_y(), x + w, pdf.get_y())

def _battery_box(pdf):
    '''
    Draws the battery management protocol box on the PDF.

    A red-highlighted box showing the critical battery levels
    and the required action at each level:
      30% – RTH (Return to Home) or safe landing
      20% – Land immediately
      10% – Automatic landing (failsafe)

    Displayed with red text (220, 53, 69) on a light red
    background (255, 240, 240) to emphasize urgency.

    Args:
        pdf – FlightPlanPDF instance
    '''
    x = pdf.l_margin + 2
    w = pdf.w - pdf.l_margin - pdf.r_margin - 4
    h = 20
    y = pdf.get_y()
    pdf.set_fill_color(255, 240, 240)
    pdf.set_draw_color(220, 53, 69)
    pdf.set_line_width(0.6)
    pdf.rect(x, y, w, h, style="DF")
    pdf.set_text_color(220, 53, 69)
    pdf.set_font(pdf.FONT, "B", 10)
    pdf.set_xy(x + 3, y + 2)
    pdf.cell(0, 5, "30%:  RTH or safe landing")
    pdf.set_xy(x + 3, y + 7)
    pdf.cell(0, 5, "20%:  LAND IMMEDIATELY")
    pdf.set_xy(x + 3, y + 12)
    pdf.cell(0, 5, "10%:  Automatic landing")
    pdf.set_y(y + h + 3)
    pdf.set_text_color(0, 0, 0)

def _abort_box(pdf):
    '''
    Draws the abort criteria box on the PDF.

    A red-highlighted box listing 5 conditions that require the
    pilot to immediately abort the flight mission:
      * Unauthorised person/animal in flight area
      * A3 distances cannot be maintained (150 m / 50 m)
      * Wind direction shift > 15 deg during flight
      * Sensor shows unreasonably high concentration
      * Loss of link or GPS (failsafe takes over)

    Displayed with red text (220, 53, 69) on a light red
    background (255, 240, 240) to emphasize urgency.

    Args:
        pdf – FlightPlanPDF instance
    '''
    x = pdf.l_margin + 2
    w = pdf.w - pdf.l_margin - pdf.r_margin - 4
    items = [
        "Unauthorised person/animal in flight area",
        "A3 distances cannot be maintained (150m/50m)",
        "Wind direction shift > 15 deg. during flight",
        "Sensor shows unreasonably high concentration",
        "Loss of link or GPS - failsafe takes over",
    ]
    h = len(items) * 5.5 + 4
    y = pdf.get_y()
    pdf.set_fill_color(255, 240, 240)
    pdf.set_draw_color(220, 53, 69)
    pdf.set_line_width(0.6)
    pdf.rect(x, y, w, h, style="DF")
    pdf.set_text_color(220, 53, 69)
    pdf.set_font(pdf.FONT, "B", 9)
    for i, item in enumerate(items):
        pdf.set_xy(x + 3, y + 2 + i * 5.5)
        pdf.cell(0, 5, "*  " + item)
    pdf.set_y(y + h + 3)
    pdf.set_text_color(0, 0, 0)

def _key_figures_table(pdf, fields):
    '''
    Draws the key figures summary table on the PDF.
    
    A two-column table (Parameter | Value) listing all critical
    flight parameters for quick reference:
      - Flight altitude, speed, wind
      - Line spacing, turn type, hover per line end
      - Background, battery landing, MTOM
      - Min. distance (A3), certificate
    The "fields" dict must contain:
      altitude, speed, wind_meas, line_spacing, hover_s

    Args:
        pdf    – FlightPlanPDF instance
        fields – dict with flight parameter values
    '''
    x = pdf.l_margin + 2
    w = pdf.w - pdf.l_margin - pdf.r_margin - 4
    col1, col2 = w * 0.5, w * 0.5
    row_h = 7
    rows = [
        ("Flight altitude", f"{fields['altitude']} m"),
        ("Speed", f"{fields['speed']} m/s"),
        ("Wind", fields["wind_meas"]),
        ("Line spacing", f"{fields['line_spacing']} m"),
        ("Turn type", f"90 deg x2 ({fields['extension']} m extension)"),
        ("Hover per line end", f"{fields['hover_s']} s"),
        ("Background", "30-60 s, before + after"),
        ("Battery landing", "20%"),
        ("MTOM", "< 25 kg"),
        ("Min. distance (A3)", "150 m / 50 m"),
        ("Certificate", "A1/A3"),
    ]
    pdf.set_fill_color(240, 240, 240)
    pdf.set_draw_color(0, 0, 0)
    pdf.set_line_width(0.3)
    y = pdf.get_y()
    pdf.set_font(pdf.FONT, "B", 10)
    pdf.rect(x, y, w, row_h, style="DF")
    pdf.set_xy(x, y)
    pdf.cell(col1, row_h, "  Parameter")
    pdf.cell(col2, row_h, "  Value",
             new_x="LMARGIN", new_y="NEXT")
    pdf.set_font(pdf.FONT, "", 10)
    for param, val in rows:
        y = pdf.get_y()
        pdf.rect(x, y, w, row_h, style="D")
        pdf.set_xy(x, y)
        pdf.cell(col1, row_h, "  " + param)
        pdf.cell(col2, row_h, "  " + val,
                 new_x="LMARGIN", new_y="NEXT")
    pdf.line(x, pdf.get_y(), x + w, pdf.get_y())

def generate_pdf(path, fields):
    """
    Generates the complete 3-page flight plan PDF.
    Page 1 – Flight details, mission, weather, battery, background:
      1. FLIGHT DETAILS   – date, kmz_name, pilot, spotter
      2. MISSION          – altitude, speed, line spacing, n_lines,
                            total distance, flight area
      3. WEATHER          – requirements vs. measured values (table)
      4. BATTERY PROTOCOL – 30%/20%/10% thresholds
      5. BACKGROUND       – measurement protocol (before/after)
    Page 2 – Pre-flight, in-flight, post-flight, abort criteria:
      6. PRE-FLIGHT       – 13 checklist items
      7. IN-FLIGHT        – 8 checklist items
      8. POST-FLIGHT      – 8 checklist items
      9. ABORT CRITERIA   – 5 conditions for immediate abort
    Page 3 – Turn time, key figures:
      10. TURN TIME       – geometry, formula, filtering
      11. KEY FIGURES     – 12-row parameter summary table
    Args:
        path   – output file path (e.g. "mission.pdf")
        fields – dict with all flight plan values:
                  date, kml_name, pilot, spotter,
                  altitude, speed, line_spacing, n_lines,
                  total_dist, flight_time,
                  rect_w, rect_h,
                  wind_meas, wind_ok,
                  rain, rain_ok,
                  temp_c, rh_pct,
                  hover_s
    """
    pdf = FlightPlanPDF()
    pdf.set_font(pdf.FONT, "B", 24)
    pdf.set_text_color(0, 0, 0)
    pdf.cell(0, 12, "Flight Plan", align="C",
             new_x="LMARGIN", new_y="NEXT")
    pdf.set_font(pdf.FONT, "", 13)
    pdf.cell(0, 8,
             "Methane Measurement - DJI M400 - A1/A3 "
             "- EU Open Category",
             align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(2)
    pdf.set_draw_color(0, 102, 204)
    pdf.set_line_width(0.5)
    pdf.line(pdf.l_margin, pdf.get_y(),
             pdf.w - pdf.r_margin, pdf.get_y())
    pdf.ln(8)

    # Flight Details section
    pdf._section_box("1. FLIGHT DETAILS", (0, 102, 204))
    pdf._item_row("Date:", fields["date"])
    pdf._item_row("KMZ file:", fields["kml_name"])
    pdf._item_row("Pilot:", fields["pilot"])
    pdf._item_row("Spotter:", fields["spotter"])
    pdf.ln(4)

    # Mission section
    pdf._section_box("2. MISSION", (0, 153, 76))
    pdf._bullet(
        f"Lawnmower pattern with 90-degree turns, "
        f"{fields['extension']} m extensions, "
        f"{fields['hover_s']} s hover per line end",
        "Pattern: ")
    pdf._bullet(
        f"{fields['altitude']} m  |  "
        f"{fields['speed']} m/s",
        "Altitude/Speed: ")
    pdf._bullet(
        f"{fields['line_spacing']} m  |  "
        f"{fields['n_lines']} lines",
        "Spacing/Lines: ")
    pdf._bullet(
        f"{fields['total_dist']} m  |  "
        f"{fields['flight_time']} s",
        "Total distance/Time: ")
    pdf._bullet(
        f"{fields['rect_w']} m x {fields['rect_h']} m  ",
        "Flight area: ")
    pdf.ln(4)

    # Weather section
    pdf._section_box(
        "3. WEATHER (FLY ONLY IF ALL ARE MET)",
        (218, 165, 32))
    _weather_table(pdf, fields)
    pdf.ln(4)

    # Battery protocol section
    pdf._section_box("4. BATTERY PROTOCOL", (220, 53, 69))
    _battery_box(pdf)
    pdf.ln(4)

    # Background measurement section
    pdf._section_box("5. BACKGROUND (MEASURED EVERY TIME)",
                     (218, 165, 32))
    pdf._bullet(
        "Before mission: 30-60 s in clean air "
        "(upwind/crosswind, same altitude)")
    pdf._bullet("After mission: same point, "
                "same altitude")
    pdf._bullet("Background = average "
                "of before + after")
    pdf._bullet("In analysis: CH4_delta "
                "= CH4 - background")
    pdf.add_page()

    # Pre-flight checklist section
    pdf._section_box("6. PRE-FLIGHT", (0, 102, 204))
    for item in [
        "KMZ file: opened + verified in "
        "Google Earth / DJI Pilot 2",
        "A3 distances confirmed: "
        "150 m buildings, 50 m people",
        "Safety at release: "
        "no smoking, no ignition sources",
        "A1/A3 certificate valid + "
        "registration + insurance",
        "NOTAMs, drone zones, "
        "geo-zone checks completed",
        "Visual inspection: frame, "
        "propellers, motors",
        "Batteries > 80%, no damage",
        "MTOM < 25 kg confirmed",
        "Methane sensor calibrated + documented",
        "GPS >= 12 satellites, RTK "
        "(if relevant)",
        "RTH altitude above obstacles",
        "Anemometer ready",
        "Spare equipment: propellers, "
        "batteries, charger, first aid",
    ]:
        pdf._bullet(item)
    pdf.ln(4)

    # In-flight checklist section
    pdf._section_box("7. IN-FLIGHT", (204, 102, 0))
    for item in [
        "Background measurement flown (30-60 s)",
        "Lines parallel with wind, "
        "consistent altitude",
        "90-degree turns (2x) at each "
        f"line end ({fields['extension']} m extension)",
        "Log: time, position, altitude "
        "at each pass",
        "Turns logged "
        "(CURVED CONTINUES)",
        "Battery monitored continuously",
        "VLOS maintained (spotter watches "
        "for people/obstacles)",
        "After mission: background "
        "measurement at same point",
    ]:
        pdf._bullet(item)
    pdf.ln(4)

    # Post-flight checklist section
    pdf._section_box("8. POST-FLIGHT", (0, 153, 76))
    for item in [
        "Visual inspection after landing",
        "Data downloaded + backup to 2 locations",
        "SD card ID / file names noted",
        "Photo documentation: setup, "
        "wind measurement, drone (before/after)",
        "Flight log book completed",
        "Background (before/after) calculated - "
        "do they match?",
        "Turn data (hover periods) "
        "filtered out",
        "CH4 profile plotted per line",
    ]:
        pdf._bullet(item)
    pdf.ln(4)

    # Abort criteria section
    pdf._section_box("9. ABORT CRITERIA", (220, 53, 69))
    _abort_box(pdf)
    pdf.add_page()

    # Turn time section
    pdf._section_box("10. TURN TIME", (218, 165, 32))
    pdf._bullet(
        f"2x 90 deg ({fields['extension']} m extension, "
        f"cross, {fields['extension']} m in)",
        "Turn geometry: ")
    pdf._bullet(
        f"t_turn = (spacing + {2*fields['extension']}) / v "
        f"+ t_acc + t_dec", "Formula: ")
    pdf._bullet(
        f"Hover periods ({fields['hover_s']} s) "
        "= filtered out",
        "Non-representative: ")
    pdf.ln(4)

    # Key figures section
    pdf._section_box("11. KEY FIGURES", (220, 53, 69))
    _key_figures_table(pdf, fields)

    # Save PDF to file
    pdf.output(path)

# ════════════════════════════════════════════
#  GUI
# ════════════════════════════════════════════

class FlightPlanApp:
    '''
    Main GUI class for the Drone Flight Plan application.

    Builds a scrollable tkinter interface with all flight plan
    configuration fields, safety/airspace/weather checks, route
    generation, and KMZ+PDF export.

    GUI Sections (top to bottom):
      1.  Flight Details        – date, pilot, spotter
      2.  Polygon               – text field (decimal or DMS) + KMZ import
      3.  Rectangle from center – alternative polygon definition
      4.  Wind & Weather        – wind direction, speed, rain, temp, RH
      5.  Parameters            – spacing, altitude, speed, repeats,
                                  flight direction (parallel/cross)
      6.  Convex hull           – checkbox for hull-based route
      7.  Safety check          – A3 building distances (Overpass API)
      8.  Airspace check        – no-fly zones (OSM/Overpass)
      9.  Weather               – Open-Meteo forecast/archive
      10. Export                – output base name, Generate route, Save (KMZ+PDF)

    Coordinate input accepts:
      Decimal: 56.493029 9.560055
      DMS (min): 56°29.6055'N 9°33.6107'E
      DMS (sec): 56°29'34.19"N 9°33'45.45"E
    (All formats can be mixed in the same polygon field.)

    Args:
        root – tkinter.Tk() root window
    '''

    def __init__(self, root):
        '''
        Initializes the GUI: sets title, geometry (640x1350),
        all instance variables for route data and check results,
        then builds the layout.
        '''
        self.root = root
        self.root.title("Drone Flight Plan")
        self.root.geometry("640x1350")
        self.root.resizable(True, True)
        self.waypoints = []
        self.wp_actions = []
        self.line_length = 0
        self.total_distance = 0
        self.n_lines = 0
        self.rect_w = 0
        self.rect_h = 0
        self.hull_ll = None
        self.safety_result = None
        self.weather_result = None
        self.airspace_result = None
        self.turn_warning = None
        self.poly_simple = True
        self.line_bearing = 0
        self._build_layout()

    def _build_layout(self):
        '''
        Builds the main window layout:
          - Header bar (title label)
          - Scrollable canvas area (content frame inside)
          - Vertical scrollbar
          - Status bar (bottom)

        Mouse wheel scrolls the canvas. Text widgets handle their
        own scroll (excluded via isinstance check).
        '''
        header = ttk.Frame(
            self.root, padding=(10, 10, 10, 5))
        header.pack(fill="x")
        ttk.Label(header, text="Drone Flight Plan",
                  font=("Arial", 16, "bold")
                  ).pack(side="left")
        canvas = tk.Canvas(self.root)
        scrollbar = ttk.Scrollbar(
            self.root, orient="vertical",
            command=canvas.yview)
        self.content = ttk.Frame(canvas)
        self.content.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")))
        canvas.create_window(
            (0, 0), window=self.content,
            anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both",
                    expand=True, padx=(10, 0),
                    pady=5)
        scrollbar.pack(side="left", fill="y", pady=5)

        def _on_mousewheel(event):
            if isinstance(event.widget, tk.Text):
                return
            canvas.yview_scroll(
                -1 if event.delta > 0 else 1,
                "units")
        canvas.bind_all("<MouseWheel>", _on_mousewheel)

        self._build_content()
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(
            self.root,
            textvariable=self.status_var,
            foreground="blue",
            padding=(10, 5)
        ).pack(fill="x", side="bottom")

    def _update_rect_labels(self, *args):
        '''
        Updates the centre coordinate labels based on the EPSG code.
        EPSG 4326 (empty) → "Center lat:" / "lon:"
        Other EPSG (e.g. 25832) → "East (m):" / "North (m):"
        '''
        epsg = self.rect_epsg_var.get().strip()
        if not epsg or epsg == "4326":
            self.rect_lat_label.config(text="Center lat:")
            self.rect_lon_label.config(text="lon:")
        else:
            self.rect_lat_label.config(text="East (m):")
            self.rect_lon_label.config(text="North (m):")

    def _build_content(self):
        '''
        Builds all GUI sections inside the scrollable content frame.
        '''
        c = self.content

        # -- 1. Flight details --
        s1 = ttk.LabelFrame(
            c, text="1. Flight Details", padding=10)
        s1.pack(fill="x", pady=5)
        self.s1_vars = {}
        for label, key, default in [
            ("Date:", "date", datetime.now().strftime("%d.%m.%Y")),
            ("Pilot:", "pilot", ""),
            ("Spotter:", "spotter", ""),
        ]:
            row = ttk.Frame(s1)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=label, width=20,
                      anchor="w").pack(side="left")
            var = tk.StringVar(value=default)
            ttk.Entry(row, textvariable=var,
                      width=20).pack(side="left", padx=5)
            self.s1_vars[key] = var

        # -- 2. Polygon --
        poly_frame = ttk.LabelFrame(
            c, text="2. Polygon (lat lon per line)",
            padding=10)
        poly_frame.pack(fill="x", pady=5)
        self.poly_text = tk.Text(
            poly_frame, height=6,
            font=("Consolas", 10))
        self.poly_text.pack(fill="x")
        self.poly_text.insert("1.0",
            "55.6750 12.5680\n"
            "55.6755 12.5685\n"
            "55.6760 12.5682\n"
            "55.6765 12.5690\n"
            "55.6770 12.5695\n"
            "55.6775 12.5700\n"
            "55.6773 12.5705\n"
            "55.6768 12.5702\n"
            "55.6765 12.5708\n"
            "55.6760 12.5703\n"
            "55.6755 12.5708\n"
            "55.6750 12.5700")
        btn_row = ttk.Frame(poly_frame)
        btn_row.pack(fill="x", pady=(5, 0))
        ttk.Button(
            btn_row, text="Load polygon...",
            command=self.import_polygon
        ).pack(side="left", padx=2)
        self.poly_count = tk.StringVar(value="12 vertices")
        ttk.Label(btn_row, textvariable=self.poly_count,
                  foreground="gray").pack(side="left", padx=10)

        # -- 3. Rectangle from center --
        rect_input_frame = ttk.LabelFrame(
            c, text="3. Rectangle from center (alternative method)",
            padding=10)
        rect_input_frame.pack(fill="x", pady=5)
        row = ttk.Frame(rect_input_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="EPSG code:",
                  width=12, anchor="w").pack(side="left")
        self.rect_epsg_var = tk.StringVar(value="25832")
        ttk.Entry(row, textvariable=self.rect_epsg_var,
                  width=10).pack(side="left", padx=4)
        ttk.Label(row, text="(empty = 4326)",
                  foreground="gray").pack(side="left", padx=4)
        row = ttk.Frame(rect_input_frame)
        row.pack(fill="x", pady=2)
        self.rect_lat_label = ttk.Label(row, text="East (m):",
                  width=12, anchor="w")
        self.rect_lat_label.pack(side="left")
        self.rect_center_lat_var = tk.StringVar(value="534446.9")
        ttk.Entry(row, textvariable=self.rect_center_lat_var,
                  width=14).pack(side="left", padx=4)
        self.rect_lon_label = ttk.Label(row, text="North (m):", width=10)
        self.rect_lon_label.pack(side="left")
        self.rect_center_lon_var = tk.StringVar(value="6261083.2")
        ttk.Entry(row, textvariable=self.rect_center_lon_var,
                  width=14).pack(side="left", padx=4)
        row = ttk.Frame(rect_input_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Direction 1 (deg from N):",
                  width=24, anchor="w").pack(side="left")
        self.rect_bearing1_var = tk.StringVar(value="0")
        ttk.Entry(row, textvariable=self.rect_bearing1_var,
                  width=10).pack(side="left", padx=4)
        row = ttk.Frame(rect_input_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Half length 1 (m):",
                  width=24, anchor="w").pack(side="left")
        self.rect_dist1_var = tk.StringVar(value="50")
        ttk.Entry(row, textvariable=self.rect_dist1_var,
                  width=10).pack(side="left", padx=4)
        ttk.Label(row, text="(center to edge)",
                  foreground="gray").pack(side="left", padx=4)
        row = ttk.Frame(rect_input_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Half length 2 (m):",
                  width=24, anchor="w").pack(side="left")
        self.rect_dist2_var = tk.StringVar(value="7.5")
        ttk.Entry(row, textvariable=self.rect_dist2_var,
                  width=10).pack(side="left", padx=4)
        ttk.Label(row, text="(direction 2 = direction 1 + 90)",
                  foreground="gray").pack(side="left", padx=4)
        ttk.Button(rect_input_frame,
                   text="Generate rectangle -> polygon",
                   command=self.generate_rect_polygon
        ).pack(anchor="w", pady=(5, 0))
        self.rect_result_var = tk.StringVar(value="")
        ttk.Label(rect_input_frame,
                  textvariable=self.rect_result_var,
                  foreground="gray").pack(anchor="w", pady=(3, 0))
        self.rect_epsg_var.trace_add(
            "write", self._update_rect_labels)

        # -- 4. Wind & Weather --
        s3 = ttk.LabelFrame(c, text="4. Wind & Weather", padding=10)
        s3.pack(fill="x", pady=5)
        self.wind_var = tk.StringVar(value="270")
        self.wind_speed_var = tk.StringVar(value="3")
        self.wx_rain_var = tk.StringVar(value="no")
        self.wx_temp_var = tk.StringVar(value="")
        self.wx_rh_var = tk.StringVar(value="")
        for label, var in [
            ("Wind direction (deg):", self.wind_var),
            ("Wind speed (m/s):", self.wind_speed_var),
            ("Rain/fog:", self.wx_rain_var),
            ("Temperature (°C):", self.wx_temp_var),
            ("RH (%):", self.wx_rh_var),
        ]:
            row = ttk.Frame(s3)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=label,
                      width=26, anchor="w").pack(side="left")
            ttk.Entry(row, textvariable=var,
                      width=10).pack(side="left", padx=5)

        # -- 5. Parameters --
        param_frame = ttk.LabelFrame(
            c, text="5. Parameters", padding=10)
        param_frame.pack(fill="x", pady=5)
        self.param_vars = {}
        for label, key, default in [
            ("Line spacing (m):", "line_spacing", 5),
            ("Altitude (m):", "altitude", 30),
            ("Speed (m/s):", "speed", 1),
        ]:
            row = ttk.Frame(param_frame)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=label,
                      width=20, anchor="w").pack(side="left")
            var = tk.StringVar(value=str(default))
            ttk.Entry(row, textvariable=var,
                      width=8).pack(side="left", padx=5)
            self.param_vars[key] = var

        # -- Hover time --
        hover_row = ttk.Frame(param_frame)
        hover_row.pack(fill="x", pady=2)
        ttk.Label(hover_row, text="Hover per line end (s):",
                  width=20, anchor="w").pack(side="left")
        self.hover_var = tk.IntVar(value=0)
        ttk.Spinbox(hover_row, from_=0, to=120,
                    textvariable=self.hover_var,
                    width=6).pack(side="left", padx=5)
        ttk.Label(hover_row, text="(0 = no hover)",
                  foreground="gray").pack(side="left", padx=4)

        # -- Extension --
        ext_row = ttk.Frame(param_frame)
        ext_row.pack(fill="x", pady=2)
        ttk.Label(ext_row, text="Extension (m):",
                  width=20, anchor="w").pack(side="left")
        self.extension_var = tk.IntVar(value=5)
        ttk.Spinbox(ext_row, from_=0, to=50,
                    textvariable=self.extension_var,
                    width=6).pack(side="left", padx=5)
        ttk.Label(ext_row, text="(0 = no extension)",
                  foreground="gray").pack(side="left", padx=4)

        # -- Flight direction (parallel / orthogonal) --
        dir_row = ttk.Frame(param_frame)
        dir_row.pack(fill="x", pady=2)
        ttk.Label(dir_row, text="Flight direction:",
                  width=20, anchor="w").pack(side="left")
        self.flight_dir_var = tk.StringVar(value="parallel")
        ttk.Radiobutton(dir_row, text="Parallel",
                        variable=self.flight_dir_var,
                        value="parallel"
        ).pack(side="left", padx=4)
        ttk.Radiobutton(dir_row, text="Cross (+90°)",
                        variable=self.flight_dir_var,
                        value="orthogonal"
        ).pack(side="left", padx=4)

        # -- Turn warning --
        self.turn_warn_var = tk.StringVar(value="")
        self.turn_warn_label = ttk.Label(
            param_frame, textvariable=self.turn_warn_var,
            foreground="red", wraplength=480, justify="left")
        self.turn_warn_label.pack(anchor="w", pady=(4, 0))

        # -- Rectangle info --
        self.rect_info_var = tk.StringVar(value="")
        ttk.Label(param_frame,
                  textvariable=self.rect_info_var,
                  foreground="gray").pack(anchor="w", pady=(4, 0))

        # -- 6. Convex hull --
        hull_frame = ttk.LabelFrame(c, text="6. Convex hull", padding=10)
        hull_frame.pack(fill="x", pady=5)
        self.hull_use_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(hull_frame,
                        text="Generate route in convex "
                             "hull (instead of polygon)",
                        variable=self.hull_use_var
        ).pack(anchor="w", pady=2)
        self.hull_info_var = tk.StringVar(value="")
        ttk.Label(hull_frame,
                  textvariable=self.hull_info_var,
                  foreground="gray").pack(anchor="w", pady=(2, 0))

        # -- 7. Safety check --
        safe_frame = ttk.LabelFrame(
            c, text="7. Safety check (A3 distances)", padding=10)
        safe_frame.pack(fill="x", pady=5)
        ttk.Button(safe_frame,
                   text="Check building distances...",
                   command=self.check_safety).pack(side="left", pady=2)
        self.safety_var = tk.StringVar(value="")
        self.safety_label = ttk.Label(
            safe_frame, textvariable=self.safety_var,
            foreground="gray", wraplength=480, justify="left")
        self.safety_label.pack(fill="x", pady=(5, 0))

        # -- 8. Airspace check --
        airspace_frame = ttk.LabelFrame(
            c, text="8. Airspace check (automated via OSM)",
            padding=10)
        airspace_frame.pack(fill="x", pady=5)
        ar_row2 = ttk.Frame(airspace_frame)
        ar_row2.pack(fill="x", pady=2)
        ttk.Button(ar_row2, text="Check airspace...",
                   command=self.check_airspace_gui
        ).pack(side="left", padx=2)
        self.airspace_var = tk.StringVar(value="")
        self.airspace_label = ttk.Label(
            airspace_frame, textvariable=self.airspace_var,
            foreground="gray", wraplength=480, justify="left")
        self.airspace_label.pack(fill="x", pady=(5, 0))

        # -- 9. Weather --
        weather_frame = ttk.LabelFrame(
            c, text="9. Weather (Open-Meteo) - flight time",
            padding=10)
        weather_frame.pack(fill="x", pady=5)
        wrow = ttk.Frame(weather_frame)
        wrow.pack(fill="x", pady=2)
        ttk.Label(wrow, text="Date (DD.MM.YYYY):",
                  width=20, anchor="w").pack(side="left")
        self.weather_date_var = tk.StringVar(value="")
        ttk.Entry(wrow, textvariable=self.weather_date_var,
                  width=12).pack(side="left", padx=4)
        ttk.Label(wrow, text="(DD.MM.YYYY, empty = today)",
                  foreground="gray").pack(side="left", padx=4)
        wrow = ttk.Frame(weather_frame)
        wrow.pack(fill="x", pady=2)
        ttk.Label(wrow, text="Start time:",
                  width=10, anchor="w").pack(side="left")
        self.weather_start_var = tk.StringVar(value="12:00")
        ttk.Entry(wrow, textvariable=self.weather_start_var,
                  width=6).pack(side="left", padx=4)
        ttk.Label(wrow, text="End time:").pack(
            side="left", padx=(10, 0))
        self.weather_end_var = tk.StringVar(value="")
        ttk.Entry(wrow, textvariable=self.weather_end_var,
                  width=6).pack(side="left", padx=4)
        ttk.Label(wrow, text="(empty = auto from flight time)",
                  foreground="gray").pack(side="left", padx=4)
        wrow = ttk.Frame(weather_frame)
        wrow.pack(fill="x", pady=2)
        ttk.Label(wrow, text="Interval (min):",
                  width=14, anchor="w").pack(side="left")
        self.weather_interval_var = tk.StringVar(value="10")
        ttk.Entry(wrow, textvariable=self.weather_interval_var,
                  width=6).pack(side="left", padx=4)
        wrow2 = ttk.Frame(weather_frame)
        wrow2.pack(fill="x", pady=2)
        ttk.Button(wrow2,
                   text="Check weather for flight time...",
                   command=self.check_weather_gui
        ).pack(side="left")
        self.weather_auto_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(wrow2, text="Auto-fill wind",
                        variable=self.weather_auto_var
        ).pack(side="left", padx=10)
        self.weather_var = tk.StringVar(value="")
        self.weather_label = ttk.Label(
            weather_frame, textvariable=self.weather_var,
            foreground="gray", wraplength=520,
            justify="left", font=("Consolas", 8))
        self.weather_label.pack(fill="x", pady=(5, 0))

        # -- 10. Export --
        fname_frame = ttk.LabelFrame(
            c, text="10. Export", padding=10)
        fname_frame.pack(fill="x", pady=5)

        row = ttk.Frame(fname_frame)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Filename:",
                  width=11, anchor="w").pack(side="left")
        self.filename_var = tk.StringVar(value="mission")
        ttk.Entry(row, textvariable=self.filename_var,
                  width=24).pack(side="left", padx=4)
        ttk.Label(row, text="(.kmz + .pdf)",
                  foreground="gray").pack(side="left")

        row = ttk.Frame(fname_frame)
        row.pack(fill="x", pady=2)
        ttk.Button(row, text="Generate route",
                   command=self.generate_route).pack(
            side="left", padx=5)
        self.export_btn = ttk.Button(
            row, text="Save (KMZ + PDF)",
            command=self.export_all, state="disabled")
        self.export_btn.pack(side="left", padx=5)
        self.wp_count = tk.StringVar(value="")
        ttk.Label(row, textvariable=self.wp_count,
                  foreground="gray").pack(side="left", padx=10)

    # ─── Methods ───────────────────────────────────
    def generate_rect_polygon(self):
        '''
        Generates a rectangle polygon from center coordinates,
        direction, and half-lengths. Writes the 4 corners to the
        polygon text field in DMS format.

        Coordinate input:
          - EPSG empty or "4326": coordinates are WGS84 decimal
            degrees (lat, lon).
          - EPSG other (e.g. "25832"): coordinates are projected
            (East, North in metres) and are transformed to
            WGS84 via pyproj before use.

        If flight direction is set to "orthogonal" (Cross +90),
        the two half-lengths are swapped so that the longer
        axis is perpendicular to the wind.

        Returns:
            None (writes to GUI fields).

        Raises:
            ValueError:
                "Invalid coordinate. Use e.g. 56.493029 or
                 56°29.6055'N" – if parsing fails
                "Lat must be -90 to 90" – if latitude out of range
                "Lon must be -180 to 180" – if longitude out of range
                "Half length 1 must be > 0 m" – if d1 <= 0
                "Half length 2 must be > 0 m" – if d2 <= 0
        '''
        try:
            center_lat = self._parse_single_coord(
                self.rect_center_lat_var.get())
            center_lon = self._parse_single_coord(
                self.rect_center_lon_var.get())
            if center_lat is None or center_lon is None:
                raise ValueError(
                    "Invalid coordinate. Use e.g. "
                    "56.493029 or 56°29.6055'N")
            epsg_str = self.rect_epsg_var.get().strip()
            if epsg_str:
                epsg = "EPSG:" + epsg_str
                t = Transformer.from_crs(epsg, "EPSG:4326",
                                         always_xy=True)
                center_lon, center_lat = t.transform(
                    center_lat, center_lon)
            bearing = float(self.rect_bearing1_var.get()) % 360
            d1 = float(self.rect_dist1_var.get())
            d2 = float(self.rect_dist2_var.get())
            if self.flight_dir_var.get() == "orthogonal":
                d1, d2 = d2, d1
            if not (-90 <= center_lat <= 90):
                raise ValueError("Lat must be -90 to 90")
            if not (-180 <= center_lon <= 180):
                raise ValueError("Lon must be -180 to 180")
            if d1 <= 0:
                raise ValueError("Half length 1 must be > 0 m")
            if d2 <= 0:
                raise ValueError("Half length 2 must be > 0 m")
            corners = rect_from_center(
                center_lat, center_lon, bearing, d1, d2)
            self.poly_text.delete("1.0", "end")
            self.poly_text.insert(
                "1.0", "\n".join(
                    self._to_dms(lat, lon)
                    for lat, lon in corners))
            self.poly_count.set("4 vertices (rectangle from center)")
            swap_note = " (swapped for Cross mode)" \
                if self.flight_dir_var.get() == "orthogonal" else ""
            self.rect_result_var.set(
                f"Rectangle: {2*d1:.0f} m x {2*d2:.0f} m | "
                f"direction {bearing:.0f} / {(bearing+90)%360:.0f}"
                f"{swap_note}")
        except ValueError as e:
            messagebox.showerror("Error", str(e))
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def check_airspace_gui(self):
        '''
        Runs the OSM airspace check for the current polygon.
        '''
        try:
            polygon = self._get_polygon()
        except ValueError as e:
            messagebox.showerror("Error", str(e))
            return
        self.airspace_var.set("Fetching airspace data...")
        self.airspace_label.config(foreground="gray")
        self.root.update_idletasks()
        result = check_airspace_osm(polygon)
        self.airspace_result = result
        if result["error"]:
            self.airspace_var.set(f"Error: {result['error']}")
            self.airspace_label.config(foreground="red")
            return
        if result["ok"] and not result["problems"]:
            self.airspace_var.set(
                "OK: No airspace restrictions found in area\n"
                "(OSM data - always check the National "
                "Airspace drone map as well)")
            self.airspace_label.config(foreground="green")
        elif result["ok"]:
            self.airspace_var.set(
                "WARNING - check manually:\n"
                + "\n".join(result["problems"])
                + "\n\n(OSM data - verify in national "
                  "airspace authority)")
            self.airspace_label.config(foreground="orange")
        else:
            self.airspace_var.set(
                "PROHIBITED / PERMISSION:\n"
                + "\n".join(result["problems"])
                + "\n\n(OSM data - verify in national "
                  "airspace authority)")
            self.airspace_label.config(foreground="red")

    def check_weather_gui(self):
        '''
        Runs the weather check for the configured flight time.
        Auto-fills wind, temp, and RH if checkbox is enabled.
        '''
        try:
            polygon = self._get_polygon()
            lat = sum(point[0] for point in polygon) / len(polygon)
            lon = sum(point[1] for point in polygon) / len(polygon)
            date_str = self.weather_date_var.get().strip() or None
            start_time = self.weather_start_var.get().strip() or "12:00"
            end_time = self.weather_end_var.get().strip() or None
            interval_min = int(self.weather_interval_var.get())

            if end_time is None and self.waypoints:
                speed = float(self.param_vars["speed"].get())
                if speed <= 0:
                    raise ValueError("Speed must be greater than zero")
                start_hour, start_minute = map(int, start_time.split(":"))
                duration_min = max(
                    1, math.ceil(self.total_distance / speed / 60))
                end_minute = start_hour * 60 + start_minute + duration_min
                if end_minute >= 24 * 60:
                    raise ValueError(
                        "Weather checks currently support a single-day "
                        "time range only")
                end_time = f"{end_minute // 60:02d}:{end_minute % 60:02d}"
        except (ValueError, TypeError) as e:
            messagebox.showerror("Weather check", str(e))
            return

        self.weather_var.set("Fetching weather data...")
        self.weather_label.config(foreground="gray")
        self.root.update_idletasks()
        result = check_weather_range(
            lat, lon, date_str=date_str, start_time=start_time,
            end_time=end_time, interval_min=interval_min)
        self.weather_result = result
        if result["error"]:
            self.weather_var.set(f"Weather error: {result['error']}")
            self.weather_label.config(foreground="red")
            return

        lines = [
            f"{point['time']}  wind {point['wind_m_s']:.1f} m/s "
            f"from {point['wind_dir_deg']:.0f} deg, "
            f"gusts {point['wind_gusts']:.1f} m/s, "
            f"rain {point['rain_mm_h']:.1f} mm/h, "
            f"{point['temp_c']:.0f}\u00b0C, "
            f"RH {point['rh_pct']:.0f}%, "
            f"{point['description']}"
            for point in result["points"]]
        status = "OK" if result["ok"] else "WARNING"
        summary = (
            f"{status}: {result['date']} {result['start_time']} - "
            f"{result['end_time']}\n" + "\n".join(lines))
        if result["problems"]:
            summary += "\n" + "\n".join(result["problems"])
        self.weather_var.set(summary)
        self.weather_label.config(
            foreground="green" if result["ok"] else "orange")

        if self.weather_auto_var.get() and result["points"]:
            points = result["points"]
            self.wind_speed_var.set(
                f"{sum(point['wind_m_s'] for point in points) / len(points):.1f}")
            sin_sum = sum(
                math.sin(math.radians(point["wind_dir_deg"]))
                for point in points)
            cos_sum = sum(
                math.cos(math.radians(point["wind_dir_deg"]))
                for point in points)
            direction = math.degrees(math.atan2(sin_sum, cos_sum)) % 360
            self.wind_var.set(f"{direction:.0f}")
            avg_temp = sum(p["temp_c"] for p in points) / len(points)
            avg_rh = sum(p["rh_pct"] for p in points) / len(points)
            self.wx_temp_var.set(f"{avg_temp:.1f}")
            self.wx_rh_var.set(f"{avg_rh:.0f}")
            if any(p["rain"] or p["fog"] for p in points):
                self.wx_rain_var.set("yes")
            else:
                self.wx_rain_var.set("no")

    def import_polygon(self):
        '''
        Opens a file dialog to import a polygon from a KML or KMZ file.
        '''
        path = filedialog.askopenfilename(
            title="Select KML/KMZ file",
            filetypes=[("KML/KMZ", "*.kml *.kmz"),
                       ("All files", "*.*")])
        if not path:
            return
        try:
            if path.lower().endswith(".kmz"):
                with zipfile.ZipFile(path) as zf:
                    name = next(
                        n for n in zf.namelist()
                        if n.endswith(".kml"))
                    content = zf.read(name).decode("utf-8")
            else:
                with open(path, encoding="utf-8") as f:
                    content = f.read()
            match = re.search(
                r"<coordinates>\s*(.*?)\s*</coordinates>",
                content, re.DOTALL)
            if not match:
                messagebox.showerror(
                    "Error", "No <coordinates> found in file")
                return
            points = []
            for token in match.group(1).split():
                parts = token.split(",")
                if len(parts) >= 2:
                    points.append(
                        (float(parts[1]), float(parts[0])))
            if len(points) > 1 and points[0] == points[-1]:
                points.pop()
            if len(points) < 3:
                messagebox.showerror(
                    "Error", "Polygon must have at least 3 vertices")
                return
            self.poly_text.delete("1.0", "end")
            self.poly_text.insert(
                "1.0", "\n".join(
                    self._to_dms(lat, lon) for lat, lon in points))
            self.poly_count.set(f"{len(points)} vertices")
            self.status_var.set(f"OK: Loaded {len(points)} vertices")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def check_safety(self):
        '''
        Runs the A3 building distance check for the current polygon.
        '''
        try:
            polygon = self._get_polygon()
        except ValueError as e:
            messagebox.showerror("Error", str(e))
            return
        self.safety_var.set("Checking building distances...")
        self.safety_label.config(foreground="gray")
        self.root.update_idletasks()
        result = check_buildings_a3(polygon, threshold=150)
        self.safety_result = result
        if result["min_distance"] is None:
            if result["buildings_found"] == 0:
                self.safety_var.set(
                    "OK: No buildings found within 250 m\n"
                    "    A3 distance met (150 m / 50 m)")
                self.safety_label.config(foreground="green")
            else:
                self.safety_var.set(
                    f"Warning: {result['building_name']}")
                self.safety_label.config(foreground="orange")
            return
        dist = result["min_distance"]
        name = result["building_name"]
        n = result["buildings_found"]
        if result["ok"]:
            self.safety_var.set(
                f"OK: Nearest building {dist:.0f} m ({name})\n"
                f"    {n} buildings | Requirement: 150 m | "
                f"Margin: {dist-150:.0f} m")
            self.safety_label.config(foreground="green")
        else:
            self.safety_var.set(
                f"FAIL: Nearest building {dist:.0f} m ({name})\n"
                f"    TOO CLOSE! Requirement: 150 m | "
                f"Short by: {150-dist:.0f} m")
            self.safety_label.config(foreground="red")

    def generate_route(self):
        '''
        Generates the flight route from the current polygon and parameters.
        Uses 90-degree turns (configurable extensions) with bounding
        box rectangle.
        '''
        try:
            polygon = self._get_polygon()
            params = self._get_params()
            wind = float(self.wind_var.get()) % 360
            line_bearing = (wind + 180) % 360
            if self.flight_dir_var.get() == "orthogonal":
                line_bearing = (line_bearing + 90) % 360
            self.line_bearing = line_bearing
            v = self.hover_var.get()
            hover_s = v if v is not None else 0
            extension = int(self.extension_var.get() or 5)
            self.poly_simple = is_polygon_simple(polygon)
            self.hull_ll = convex_hull_latlon(polygon)
            if not self.poly_simple:
                self.hull_info_var.set(
                    f"\u26a0 Polygon is self-intersecting \u2192 "
                    f"convex hull used automatically "
                    f"({len(self.hull_ll)} vertices)")
            elif self.hull_use_var.get():
                self.hull_info_var.set(
                    f"Convex hull used (user-selected) "
                    f"({len(self.hull_ll)} hull vs. "
                    f"{len(polygon)} polygon vertices)")
            else:
                self.hull_info_var.set(
                    f"{len(self.hull_ll)} hull vertices "
                    f"(polygon: {len(polygon)} vertices)")
            if not self.poly_simple or self.hull_use_var.get():
                route_area = self.hull_ll
            else:
                route_area = polygon
            (self.waypoints,
             self.wp_actions,
             self.line_length,
             self.total_distance,
             self.n_lines,
             self.turn_warning
            ) = generate_lawnmower(
                route_area,
                params["line_spacing"],
                line_bearing,
                params["altitude"],
                hover_s=hover_s,
                extension=extension
            )
            if self.turn_warning:
                self.turn_warn_var.set("WARNING: " + self.turn_warning)
                self.turn_warn_label.config(foreground="red")
            else:
                self.turn_warn_var.set("")
            if self.waypoints:
                self.rect_w = self.line_length
                self.rect_h = (
                    self.n_lines * params["line_spacing"])
                self.rect_info_var.set(
                    f"Rectangle: {self.rect_w:.0f} m x "
                    f"{self.rect_h:.0f} m")
                self.wp_count.set(
                    f"{len(self.waypoints)} WP | "
                    f"{self.line_length:.0f} m/line")
                self.export_btn.config(state="normal")
                area_mode = (
                    "hull" if self.hull_use_var.get()
                    else "polygon")
                if not self.poly_simple and \
                        not self.hull_use_var.get():
                    area_mode += "\u2192hull"
                dir_label = (
                    "ortho" if self.flight_dir_var.get() == "orthogonal"
                    else "par")
                self.status_var.set(
                    f"Route: {line_bearing:.0f} ({dir_label}) | "
                    f"{self.total_distance:.0f} m | "
                    f"{self.n_lines} lines ({area_mode})")
            else:
                self.wp_count.set("")
                self.rect_info_var.set("")
                self.export_btn.config(state="disabled")
                self.status_var.set("No route")
        except ValueError as e:
            messagebox.showerror("Error", str(e))
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def export_all(self):
        '''
        Exports the flight plan as KMZ + PDF
        to a user-selected folder.
        KMZ contains two files:
          template.kml  – display format (createTime, globalHeight, ellipsoidHeight)
          waylines.wpmz – execution format (distance, duration, executeHeight, waypointSpeed)
        '''
        if not self.waypoints:
            return
        base = self.filename_var.get().strip() or "mission"
        base = re.sub(r"[^A-Za-z0-9_\-]", "_", base)
        if not base:
            base = "mission"
        outdir = filedialog.askdirectory(
            title="Select folder for files")
        if not outdir:
            return
        kmz_path  = os.path.join(outdir, base + ".kmz")
        pdf_path  = os.path.join(outdir, base + ".pdf")
        try:
            speed_val = float(
                self.param_vars["speed"].get().strip() or "1")
            altitude_val = float(
                self.param_vars["altitude"].get().strip() or "30")
            v = self.hover_var.get()
            hover_s = v if v is not None else 0
            extension = int(self.extension_var.get() or 5)
            now_ms = int(time.time() * 1000)

            # Compute distance and duration
            total_dist = 0.0
            for i in range(len(self.waypoints) - 1):
                lat1, lon1, _ = self.waypoints[i]
                lat2, lon2, _ = self.waypoints[i + 1]
                dlat = math.radians(lat2 - lat1)
                dlon = math.radians(lon2 - lon1)
                a = (math.sin(dlat / 2) ** 2 +
                     math.cos(math.radians(lat1)) *
                     math.cos(math.radians(lat2)) *
                     math.sin(dlon / 2) ** 2)
                total_dist += 6371000 * 2 * math.atan2(
                    math.sqrt(a), math.sqrt(1 - a))
            hover_total = sum(h for _, h, _ in self.wp_actions)
            duration = (total_dist / speed_val + hover_total
                        if speed_val > 0 else hover_total)

            # ════════════════════════════════════════
            #  waylines.wpmz (execution format)
            # ════════════════════════════════════════
            
            wpmz = [
                '<?xml version="1.0" encoding="UTF-8"?>',
                '<kml xmlns="http://www.opengis.net/kml/2.2" '
                'xmlns:wpml="http://www.dji.com/wpmz/1.0.6">',
                '  <Document>',
                '    <wpml:missionConfig>',
                '      <wpml:flyToWaylineMode>safely</wpml:flyToWaylineMode>',
                '      <wpml:finishAction>goHome</wpml:finishAction>',
                '      <wpml:exitOnRCLost>executeLostAction</wpml:exitOnRCLost>',
                '      <wpml:executeRCLostAction>goBack</wpml:executeRCLostAction>',
                '      <wpml:takeOffSecurityHeight>20</wpml:takeOffSecurityHeight>',
                f'      <wpml:globalTransitionalSpeed>{int(speed_val)}</wpml:globalTransitionalSpeed>',
                '      <wpml:droneInfo>',
                '        <wpml:droneEnumValue>103</wpml:droneEnumValue>',
                '        <wpml:droneSubEnumValue>0</wpml:droneSubEnumValue>',
                '      </wpml:droneInfo>',
                '      <wpml:waylineAvoidLimitAreaMode>0</wpml:waylineAvoidLimitAreaMode>',
                '    </wpml:missionConfig>',
                '    <Folder>',
                '      <wpml:templateId>0</wpml:templateId>',
                '      <wpml:executeHeightMode>relativeToStartPoint</wpml:executeHeightMode>',
                '      <wpml:waylineId>0</wpml:waylineId>',
                f'      <wpml:distance>{total_dist:.15f}</wpml:distance>',
                f'      <wpml:duration>{duration:.15f}</wpml:duration>',
                f'      <wpml:autoFlightSpeed>{int(speed_val)}</wpml:autoFlightSpeed>',
                '      <wpml:realTimeFollowSurfaceByFov>0</wpml:realTimeFollowSurfaceByFov>',
            ]

            for i, (wp, action) in enumerate(
                    zip(self.waypoints, self.wp_actions)):
                lat, lon, alt = wp
                turn_ccw, hover, wp_type = action

                if wp_type in ("line_start", "line_end", "end"):
                    turn_mode = "toPointAndStopWithDiscontinuityCurvature"
                    damping = "0"
                else:
                    turn_mode = "toPointAndPassWithContinuityCurvature"
                    damping = "5.0"

                wpmz += [
                    '      <Placemark>',
                    '        <Point>',
                    '          <coordinates>',
                    f'            {lon},{lat}',
                    '          </coordinates>',
                    '        </Point>',
                    f'        <wpml:index>{i}</wpml:index>',
                    f'        <wpml:executeHeight>{int(alt)}</wpml:executeHeight>',
                    f'        <wpml:waypointSpeed>{int(speed_val)}</wpml:waypointSpeed>',
                    '        <wpml:waypointHeadingParam>',
                    '          <wpml:waypointHeadingMode>followWayline</wpml:waypointHeadingMode>',
                    '          <wpml:waypointHeadingAngle>0</wpml:waypointHeadingAngle>',
                    '          <wpml:waypointPoiPoint>0.000000,0.000000,0.000000</wpml:waypointPoiPoint>',
                    '          <wpml:waypointHeadingAngleEnable>0</wpml:waypointHeadingAngleEnable>',
                    '          <wpml:waypointHeadingPoiIndex>0</wpml:waypointHeadingPoiIndex>',
                    '        </wpml:waypointHeadingParam>',
                    '        <wpml:waypointTurnParam>',
                    f'          <wpml:waypointTurnMode>{turn_mode}</wpml:waypointTurnMode>',
                    f'          <wpml:waypointTurnDampingDist>{damping}</wpml:waypointTurnDampingDist>',
                    '        </wpml:waypointTurnParam>',
                    '        <wpml:useStraightLine>0</wpml:useStraightLine>',
                    '        <wpml:waypointGimbalHeadingParam>',
                    '          <wpml:waypointGimbalPitchAngle>0</wpml:waypointGimbalPitchAngle>',
                    '          <wpml:waypointGimbalYawAngle>0</wpml:waypointGimbalYawAngle>',
                    '        </wpml:waypointGimbalHeadingParam>',
                    '        <wpml:isRisky>0</wpml:isRisky>',
                    '        <wpml:waypointWorkType>0</wpml:waypointWorkType>',
                    '      </Placemark>',
                ]

            wpmz += [
                '    </Folder>',
                '  </Document>',
                '</kml>',
            ]
            wpmz_str = "\n".join(wpmz)

            # ════════════════════════════════════════
            #  template.kml (display format)
            # ════════════════════════════════════════
            
            template_kml = (
                '<?xml version="1.0" encoding="UTF-8"?>\n'
                '<kml xmlns="http://www.opengis.net/kml/2.2" '
                'xmlns:wpml="http://www.dji.com/wpmz/1.0.6">\n'
                '  <Document>\n'
                f'    <wpml:createTime>{now_ms}</wpml:createTime>\n'
                f'    <wpml:updateTime>{now_ms}</wpml:updateTime>\n'
                '    <wpml:missionConfig>\n'
                '      <wpml:flyToWaylineMode>safely</wpml:flyToWaylineMode>\n'
                '      <wpml:finishAction>goHome</wpml:finishAction>\n'
                '      <wpml:exitOnRCLost>executeLostAction</wpml:exitOnRCLost>\n'
                '      <wpml:executeRCLostAction>goBack</wpml:executeRCLostAction>\n'
                '      <wpml:takeOffSecurityHeight>20</wpml:takeOffSecurityHeight>\n'
                f'      <wpml:globalTransitionalSpeed>{int(speed_val)}</wpml:globalTransitionalSpeed>\n'
                '      <wpml:droneInfo>\n'
                '        <wpml:droneEnumValue>103</wpml:droneEnumValue>\n'
                '        <wpml:droneSubEnumValue>0</wpml:droneSubEnumValue>\n'
                '      </wpml:droneInfo>\n'
                '      <wpml:waylineAvoidLimitAreaMode>0</wpml:waylineAvoidLimitAreaMode>\n'
                '    </wpml:missionConfig>\n'
                '    <Folder>\n'
                '      <wpml:templateType>waypoint</wpml:templateType>\n'
                '      <wpml:templateId>0</wpml:templateId>\n'
                '      <wpml:waylineCoordinateSysParam>\n'
                '        <wpml:coordinateMode>WGS84</wpml:coordinateMode>\n'
                '        <wpml:heightMode>relativeToStartPoint</wpml:heightMode>\n'
                '        <wpml:positioningType>GPS</wpml:positioningType>\n'
                '      </wpml:waylineCoordinateSysParam>\n'
                f'      <wpml:autoFlightSpeed>{int(speed_val)}</wpml:autoFlightSpeed>\n'
                f'      <wpml:globalHeight>{int(altitude_val)}</wpml:globalHeight>\n'
                '      <wpml:caliFlightEnable>0</wpml:caliFlightEnable>\n'
                '      <wpml:gimbalPitchMode>manual</wpml:gimbalPitchMode>\n'
                '      <wpml:globalWaypointHeadingParam>\n'
                '        <wpml:waypointHeadingMode>followWayline</wpml:waypointHeadingMode>\n'
                '        <wpml:waypointHeadingAngle>0</wpml:waypointHeadingAngle>\n'
                '        <wpml:waypointPoiPoint>0.000000,0.000000,0.000000</wpml:waypointPoiPoint>\n'
                '        <wpml:waypointHeadingPoiIndex>0</wpml:waypointHeadingPoiIndex>\n'
                '      </wpml:globalWaypointHeadingParam>\n'
                '      <wpml:globalWaypointTurnMode>toPointAndStopWithDiscontinuityCurvature</wpml:globalWaypointTurnMode>\n'
                '      <wpml:globalUseStraightLine>0</wpml:globalUseStraightLine>\n'
            )

            for i, (wp, action) in enumerate(
                    zip(self.waypoints, self.wp_actions)):
                lat, lon, alt = wp
                turn_ccw, hover, wp_type = action

                if wp_type in ("line_start", "line_end", "end"):
                    turn_mode = "toPointAndStopWithDiscontinuityCurvature"
                else:
                    turn_mode = "toPointAndPassWithContinuityCurvature"

                template_kml += (
                    f'      <Placemark>\n'
                    f'        <Point>\n'
                    f'          <coordinates>\n'
                    f'            {lon},{lat}\n'
                    f'          </coordinates>\n'
                    f'        </Point>\n'
                    f'        <wpml:index>{i}</wpml:index>\n'
                    f'        <wpml:ellipsoidHeight>{int(alt)}</wpml:ellipsoidHeight>\n'
                    f'        <wpml:height>{int(alt)}</wpml:height>\n'
                    f'        <wpml:waypointTurnParam>\n'
                    f'          <wpml:waypointTurnMode>{turn_mode}</wpml:waypointTurnMode>\n'
                    f'          <wpml:waypointTurnDampingDist>0</wpml:waypointTurnDampingDist>\n'
                    f'        </wpml:waypointTurnParam>\n'
                    f'        <wpml:useGlobalSpeed>1</wpml:useGlobalSpeed>\n'
                    f'        <wpml:useGlobalHeadingParam>1</wpml:useGlobalHeadingParam>\n'
                    f'        <wpml:useStraightLine>0</wpml:useStraightLine>\n'
                    f'        <wpml:isRisky>0</wpml:isRisky>\n'
                    f'      </Placemark>\n'
                )

            template_kml += (
                '    </Folder>\n'
                '  </Document>\n'
                '</kml>'
            )

            # ════════════════════════════════════════
            #  KMZ
            # ════════════════════════════════════════
            
            with zipfile.ZipFile(kmz_path, 'w',
                                 zipfile.ZIP_DEFLATED) as zf:
                zf.writestr("template.kml", template_kml)
                zf.writestr("waylines.wpmz", wpmz_str)

            # ════════════════════════════════════════
            #  PDF
            # ════════════════════════════════════════

            date = (self.s1_vars["date"].get().strip()
                    or datetime.now().strftime("%d.%m.%Y"))
            pilot = self.s1_vars["pilot"].get().strip() or "-"
            spotter = self.s1_vars["spotter"].get().strip() or "-"
            altitude = (self.param_vars["altitude"].get().strip()
                        or "30")
            speed = self.param_vars["speed"].get().strip() or "1"
            line_spacing = (
                self.param_vars["line_spacing"].get().strip()
                or "5")
            wind_dir = self.wind_var.get().strip() or "270"
            wind_spd = self.wind_speed_var.get().strip() or "3"
            rain = self.wx_rain_var.get().strip() or "no"
            wind_meas = f"{wind_spd} m/s ({wind_dir} deg.)"
            wind_ok = "OK" if _ok_wind(wind_spd) else "NO"
            rain_ok = "OK" if rain.lower() in ("no", "n") else "NO"
            try:
                fly_time = self.total_distance / float(speed)
                hover_time = sum(h for _, h, _ in self.wp_actions)
                flight_time = f"{fly_time + hover_time:.0f}"
            except (ValueError, ZeroDivisionError):
                flight_time = "-"
            if self.weather_result and self.weather_result["points"]:
                pts = self.weather_result["points"]
                avg_temp = sum(p["temp_c"] for p in pts) / len(pts)
                avg_rh = sum(p["rh_pct"] for p in pts) / len(pts)
            else:
                avg_temp = self.wx_temp_var.get().strip() or "-"
                avg_rh = self.wx_rh_var.get().strip() or "-"
            fields = {
                "date": date,
                "kml_name": base + ".kmz",
                "pilot": pilot,
                "spotter": spotter,
                "altitude": altitude,
                "speed": speed,
                "line_spacing": line_spacing,
                "n_lines": self.n_lines,
                "total_dist": f"{self.total_distance:.0f}",
                "flight_time": flight_time,
                "rect_w": f"{self.rect_w:.0f}",
                "rect_h": f"{self.rect_h:.0f}",
                "hover_s": hover_s,
                "extension": extension,
                "wind_meas": wind_meas,
                "wind_ok": wind_ok,
                "rain": rain,
                "rain_ok": rain_ok,
                "temp_c": avg_temp,
                "rh_pct": avg_rh,
            }
            generate_pdf(pdf_path, fields)

            self.status_var.set(
                f"OK: Saved {base}.kmz + {base}.pdf")
        except Exception as e:
            import traceback
            messagebox.showerror(
                "Export error",
                f"{type(e).__name__}: {e}\n\n{traceback.format_exc()}")

    def _parse_single_coord(self, text):
        '''
        Parses a single coordinate value (latitude or longitude).

        Accepts:
          Decimal:  56.493029
          DMS (min): 56°29.6055'N
          DMS (sec): 56°29'34.19"N

        Returns:
            float in decimal degrees, or None if parsing fails.
        '''
        text = text.strip().replace(",", ".")
        if not text:
            return None
        if "°" in text:
            m = re.match(
                r"(\d+)\s*°\s*(\d+)\s*'\s*(\d+(?:\.\d+)?)\s*\"?\s*([NSEW])?",
                text)
            if m:
                val = (int(m.group(1))
                       + int(m.group(2)) / 60.0
                       + float(m.group(3)) / 3600.0)
                if m.group(4) in ("S", "W"):
                    val = -val
                return val
            m = re.match(
                r"(\d+)\s*°\s*(\d+(?:\.\d+)?)\s*'\s*([NSEW])?",
                text)
            if m:
                val = int(m.group(1)) + float(m.group(2)) / 60.0
                if m.group(3) in ("S", "W"):
                    val = -val
                return val
            return None
        try:
            return float(text)
        except ValueError:
            return None

    def _parse_coord_line(self, line):
        '''
        Parses a single line containing two coordinates.
        Accepts:
          Decimal:       56.493029 9.560055
          DMS (min):     56°29.6055'N 9°33.6107'E
          DMS (sec):     56°29'34.19"N 9°33'45.45"E
        Returns:
            (lat, lon) tuple in decimal degrees, or None if parsing fails.
        '''
        line = line.strip()
        if not line:
            return None
        if "°" in line:
            line = line.replace(",", ".")

            # DMS (sec): 56°29'34.19"N 9°33'45.45"E
            m = re.match(
                r"(\d+)\s*°\s*(\d+)\s*'\s*(\d+(?:\.\d+)?)\s*\"?\s*([NS])?\s*"
                r"(\d+)\s*°\s*(\d+)\s*'\s*(\d+(?:\.\d+)?)\s*\"?\s*([EW])?",
                line)
            if m:
                lat = (int(m.group(1))
                       + int(m.group(2)) / 60.0
                       + float(m.group(3)) / 3600.0)
                lon = (int(m.group(5))
                       + int(m.group(6)) / 60.0
                       + float(m.group(7)) / 3600.0)
                if m.group(4) == "S":
                    lat = -lat
                if m.group(8) == "W":
                    lon = -lon
                return (lat, lon)
            
            # DMS (min): 56°29.6055'N 9°33.6107'E
            m = re.match(
                r"(\d+)\s*°\s*(\d+(?:\.\d+)?)\s*'\s*([NS])?\s*"
                r"(\d+)\s*°\s*(\d+(?:\.\d+)?)\s*'\s*([EW])?",
                line)
            if m:
                lat = int(m.group(1)) + float(m.group(2)) / 60.0
                lon = int(m.group(4)) + float(m.group(5)) / 60.0
                if m.group(3) == "S":
                    lat = -lat
                if m.group(6) == "W":
                    lon = -lon
                return (lat, lon)
            return None
        
        # Decimal: 56.493029 9.560055
        parts = line.replace(",", " ").split()
        if len(parts) >= 2:
            try:
                return (float(parts[0]), float(parts[1]))
            except ValueError:
                return None
        return None

    def _to_dms(self, lat, lon):
        '''
        Converts decimal degrees to a DMS display string.
        '''
        ns = "N" if lat >= 0 else "S"
        ew = "E" if lon >= 0 else "W"
        lat_abs = abs(lat)
        lon_abs = abs(lon)
        lat_d = int(lat_abs)
        lat_m = (lat_abs - lat_d) * 60
        lon_d = int(lon_abs)
        lon_m = (lon_abs - lon_d) * 60
        return f"{lat_d}\u00b0{lat_m:.4f}'{ns} {lon_d}\u00b0{lon_m:.4f}'{ew}"

    def _get_polygon(self):
        '''
        Reads the polygon text field and returns a list of (lat, lon) tuples.
        '''
        text = self.poly_text.get("1.0", "end").strip()
        points = []
        for line in text.split("\n"):
            coord = self._parse_coord_line(line)
            if coord is not None:
                points.append(coord)
        if len(points) < 3:
            raise ValueError(
                "Polygon must have at least 3 vertices")
        return points

    def _get_params(self):
        '''
        Reads and validates the parameter fields from the GUI.

        Parses the three parameter inputs (line spacing, altitude,
        and speed) and ensures they are physically valid before
        passing them to the route generator.

        Args:
            (none – reads from self.param_vars StringVars)

        Returns:
            dict with keys:
                line_spacing – distance between flight lines (m), > 0
                altitude     – flight altitude (m), > 0
                speed        – flight speed (m/s), > 0

        Raises:
            ValueError:
                "Line spacing must be > 0 m"  – if line_spacing <= 0
                "Altitude must be > 0 m"     – if altitude <= 0
                "Speed must be > 0 m/s"      – if speed <= 0
                (also raised by float() if the field is empty or
                 contains non-numeric text)

        Note:
            All three values must be strictly positive. A value of
            zero or negative would cause an infinite loop (spacing),
            an invalid altitude, or a division by zero /
            negative flight time (speed) downstream.

            The caller (generate_route) catches ValueError and
            displays it in a messagebox, so no further handling
            is needed.
        '''
        ls = float(self.param_vars["line_spacing"].get())
        alt = float(self.param_vars["altitude"].get())
        spd = float(self.param_vars["speed"].get())
        if ls <= 0:
            raise ValueError("Line spacing must be > 0 m")
        if alt <= 0:
            raise ValueError("Altitude must be > 0 m")
        if spd <= 0:
            raise ValueError("Speed must be > 0 m/s")
        return {"line_spacing": ls, "altitude": alt, "speed": spd}
    
# ════════════════════════════════════════════
#  Helpers
# ════════════════════════════════════════════

def _ok_wind(speed):
    '''
    Checks if a wind speed value is within the acceptable range for drone flight.

    Acceptable range: 2 m/s (minimum) to 5 m/s (maximum).

    Args:
        speed – wind speed as string (e.g. "3" or "3.5")

    Returns:
        True if 2 <= speed <= 5, False otherwise (including parse errors).
    '''
    try:
        return 2 <= float(speed) <= 5
    except ValueError:
        return False

# ════════════════════════════════════════════
#  Entry Point
# ════════════════════════════════════════════

if __name__ == "__main__":
    root = tk.Tk()            # Create main window
    app = FlightPlanApp(root) # Build the GUI
    root.mainloop()           # Start event loop