#!/usr/bin/env python3
"""
heprep_to_vtp_clean.py
Convert HepRep (Geant4) -> two VTPs:
  - geometry.vtp : polygons from Detector Geometry
  - event_data_filtered.vtp : merged/filtered polylines from Event Data

Features:
 - filters out coordinates whose absolute component exceeds max_coord (default 1e7)
 - merges primitives whose endpoints are within `merge_tol` distance
 - drops segments with fewer than min_points
 - writes .vtp (VTK XML PolyData) without external deps
"""
import xml.etree.ElementTree as ET
import math, sys
import vtk
from collections import defaultdict

geometry_names = ["world_phys", "lab_phys", "LXe_phys"]

def localname(tag):
    return tag.split('}')[-1] if '}' in tag else tag

def parse_points_geom(prim, ns):
    pts = []
    for p in prim.findall(".//heprep:point", ns):
        try:
            pts.append(
                (float(p.get("x", "0")), float(p.get("y", "0")), float(p.get("z", "0")))
            )
        except ValueError:
            pass
    return pts

def parse_points(elem):
    pts = []
    for p in elem.findall(".//heprep:point", ns):
        x = p.get('x') or p.get('X'); y = p.get('y') or p.get('Y'); z = p.get('z') or p.get('Z')
        if x is None or y is None or z is None:
            txt = (p.text or "").strip()
            if txt:
                toks = txt.split()
                if len(toks) >= 3:
                    x,y,z = toks[0:3]
                else:
                    continue
            else:
                continue
        try:
            pts.append((float(x), float(y), float(z)))
        except:
            continue
    return pts

def distance(a,b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def parse_heprep(fname):
    tree = ET.parse(fname)
    root = tree.getroot()
    # namespace mapping for ElementTree find
    global ns
    ns = {'heprep': 'http://www.slac.stanford.edu/~perl/heprep/'}

    # Geometry primitives (collect polygons)
    geom_polys = []
    for prim in root.findall(".//heprep:type[@name='Detector Geometry']//heprep:primitive", ns):
        pts = parse_points(prim)
        if pts:
            geom_polys.append(pts)

    # Event primitives (polylines)
    event_primitives = []
    for prim in root.findall(".//heprep:type[@name='Event Data']//heprep:type[@name='TransientPolylines']//heprep:primitive", ns):
        pts = parse_points(prim)
        if pts:
            event_primitives.append(pts)

    return geom_polys, event_primitives

def filter_and_merge(primitives, max_coord=1e7, merge_tol=1e-3, min_points=2, drop_zero_length=True):
    """
    primitives: list of list-of-(x,y,z)
    Returns merged list-of-list-of-(x,y,z)
    """
    # first filter out any primitive that has coords exceeding max_coord
    good = []
    for pts in primitives:
        bad = False
        for x,y,z in pts:
            if abs(x) > max_coord or abs(y) > max_coord or abs(z) > max_coord:
                bad = True
                break
        if not bad and len(pts) >= min_points:
            good.append(pts)

    # merge primitives whose endpoints are within merge_tol
    # naive O(N^2) merge: repeatedly try to append primitives whose start is near end
    merged = []
    used = [False]*len(good)
    for i, pts in enumerate(good):
        if used[i]: continue
        cur = list(pts)
        used[i] = True
        changed = True
        while changed:
            changed = False
            for j, other in enumerate(good):
                if used[j]: continue
                # try append other to end
                if distance(cur[-1], other[0]) <= merge_tol:
                    cur.extend(other[1:])  # avoid duplicating shared point
                    used[j] = True
                    changed = True
                    break
                # try prepend other if other end near cur start
                if distance(other[-1], cur[0]) <= merge_tol:
                    cur = list(other[:-1]) + cur
                    used[j] = True
                    changed = True
                    break
            # end for
        # optionally drop segments with zero-length or extremely short extent
        if drop_zero_length:
            ext = math.sqrt((cur[-1][0]-cur[0][0])**2 + (cur[-1][1]-cur[0][1])**2 + (cur[-1][2]-cur[0][2])**2)
            if ext < 1e-6:
                continue
        merged.append(cur)
    return merged

def write_vtp_polydata_points_lines(fname, points, lines):
    # points: list of (x,y,z)
    # lines: list of lists of point indices (ints)
    N = len(points); M = len(lines)
    pts_flat = []
    for x,y,z in points:
        pts_flat.extend([repr(float(x)), repr(float(y)), repr(float(z))])
    connectivity = []
    offsets = []
    off = 0
    for cell in lines:
        for idx in cell:
            connectivity.append(str(int(idx)))
        off += len(cell)
        offsets.append(str(off))

    xml = '<?xml version="1.0"?>\n'
    xml += '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n'
    xml += '  <PolyData>\n'
    xml += f'    <Piece NumberOfPoints="{N}" NumberOfLines="{M}">\n'
    xml += '      <PointData />\n'
    xml += '      <Points>\n'
    xml += '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n'
    for i in range(0, len(pts_flat), 12):
        xml += '          ' + ' '.join(pts_flat[i:i+12]) + '\n'
    xml += '        </DataArray>\n'
    xml += '      </Points>\n'
    xml += '      <Lines>\n'
    xml += '        <DataArray type="Int32" Name="connectivity" format="ascii">\n'
    xml += '          ' + (' '.join(connectivity) if connectivity else '') + '\n'
    xml += '        </DataArray>\n'
    xml += '        <DataArray type="Int32" Name="offsets" format="ascii">\n'
    xml += '          ' + (' '.join(offsets) if offsets else '') + '\n'
    xml += '        </DataArray>\n'
    xml += '      </Lines>\n'
    xml += '    </Piece>\n'
    xml += '  </PolyData>\n'
    xml += '</VTKFile>\n'
    with open(fname, 'w') as f:
        f.write(xml)
    print(f"Wrote: {fname}  (points={N}, lines={M})")

def write_vtp_polys(fname, polys):
    # flatten points and create cell lists (polygons)
    points = []
    lines = []
    for poly in polys:
        idxs = []
        for pt in poly:
            idxs.append(len(points))
            points.append(pt)
        # close polygon? For ParaView a polygon can be given as poly; we store as polygon cells via Polys in separate writer if needed.
        lines.append(idxs)
    # we'll reuse the lines writer but name as polys in header (ParaView still reads)
    write_vtp_polydata_points_lines(fname, points, lines)

def write_vtp(filename, primitives, line_mode=False):
    """Write a list of primitives (each = list of (x,y,z)) to a .vtp file."""
    points = vtk.vtkPoints()
    polys = vtk.vtkCellArray()

    for prim in primitives:
        ids = []
        for x, y, z in prim:
            ids.append(points.InsertNextPoint(x, y, z))
        if line_mode:
            cell = vtk.vtkPolyLine()
            cell.GetPointIds().SetNumberOfIds(len(ids))
            for i, pid in enumerate(ids):
                cell.GetPointIds().SetId(i, pid)
            polys.InsertNextCell(cell)
        else:
            poly = vtk.vtkPolygon()
            poly.GetPointIds().SetNumberOfIds(len(ids))
            for i, pid in enumerate(ids):
                poly.GetPointIds().SetId(i, pid)
            polys.InsertNextCell(poly)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    if line_mode:
        polydata.SetLines(polys)
    else:
        polydata.SetPolys(polys)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(polydata)
    writer.Write()
    print(f"Wrote {filename} with {len(primitives)} primitives")

def main():
    if len(sys.argv) < 2:
        print("Usage: python heprep_to_vtp_clean.py input.heprep")
        sys.exit(1)
    infile = sys.argv[1]

    tree = ET.parse(infile)
    root = tree.getroot()

    geom_polys, event_prims = parse_heprep(infile)
#    print("Found geometry primitives:", len(geom_polys))
#    print("Found event primitives:", len(event_prims))

    for geom_name in geometry_names:
        path = f".//heprep:type[@name='Detector Geometry']//heprep:type[@name='{geom_name}']//heprep:primitive"
        geoms = root.findall(path, ns)
        primitives = [parse_points_geom(g, ns) for g in geoms if parse_points_geom(g, ns)]
        if primitives:
            outname = f"geometry_{geom_name}.vtp"
            write_vtp(outname, primitives, line_mode=False)
        else:
            print(f"No primitives found for {geom_name}")

    # parameters you can tweak:
    max_coord = 1e7    # filter coordinates exceeding this absolute value
    merge_tol = 1e-3   # distance tolerance to consider endpoints identical (units as in file)
    min_points = 2

    merged = filter_and_merge(event_prims, max_coord=max_coord, merge_tol=merge_tol, min_points=min_points)

    # Build point list and connectivity for merged polylines
    points = []
    lines = []
    for poly in merged:
        idxs = []
        for pt in poly:
            idxs.append(len(points))
            points.append(pt)
        lines.append(idxs)

    # write outputs
    if points and lines:
        write_vtp_polydata_points_lines("event_data_filtered.vtp", points, lines)
    else:
        print("No event primitives left after filtering/merging.")

if __name__ == "__main__":
    main()

