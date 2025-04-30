import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.geom

settings = ifcopenshell.geom.settings()
# settings.set(settings.USE_PYTHON_OPENCASCADE, False)  # Nie używamy pythonOCC
# Stałe powietrza
dENSITY_AIR = 1.2      # kg/m3
SPECIFIC_HEAT_AIR = 1005  # J/(kg·K)

# Pomocnicze odczyty PSet
def get_pset_value(ifc_elem, pset_name, prop_name):
    psets = ifcopenshell.util.element.get_psets(ifc_elem)
    pset = psets.get(pset_name, {})
    return float(pset.get(prop_name, 0.0)) if prop_name in pset else 0.0

# Wartość U z IFC
def get_u_value(ifc_element):
    # najpierw z PSet_ThermalTransmittance
    u = get_pset_value(ifc_element, 'PSet_ThermalTransmittance', 'U')
    if u > 0:
        return u
    # zdefiniowane typy
    predefined = getattr(ifc_element, 'PredefinedType', '')
    mapping = {
        'WALL': 0.3,
        'WINDOW': 1.1,
        'ROOF': 0.25,
        'SLAB': 0.4
    }
    return mapping.get(predefined, 0.35)

# Kubatura przestrzeni z atrybutów lub PSet
def get_space_volume(space):
    vol = getattr(space, 'GrossVolume', None) or getattr(space, 'NetVolume', None)
    if vol and vol > 0:
        return vol
    # z PSet_SpaceThermalProperties
    return get_pset_value(space, 'PSet_SpaceThermalProperties', 'GrossVolume')

# Obliczenia strat
def transmission_heat_loss(u_values, areas, delta_t):
    return sum(u * a for u, a in zip(u_values, areas)) * delta_t

def ventilation_heat_loss(ach, volume, delta_t):
    v_dot = ach * volume / 3600.0
    return v_dot * dENSITY_AIR * SPECIFIC_HEAT_AIR * delta_t

def compute_volume(verts, faces):
    volume = 0.0
    for i in range(0, len(faces), 3):
        idx0 = faces[i] * 3
        idx1 = faces[i+1] * 3
        idx2 = faces[i+2] * 3

        v0 = verts[idx0:idx0+3]
        v1 = verts[idx1:idx1+3]
        v2 = verts[idx2:idx2+3]

        # Obliczanie objętości czworościanu
        volume += (v0[0]*v1[1]*v2[2] + v1[0]*v2[1]*v0[2] + v2[0]*v0[1]*v1[2]
                   - v0[0]*v2[1]*v1[2] - v1[0]*v0[1]*v2[2] - v2[0]*v1[1]*v0[2]) / 6.0
    return abs(volume)

# Główna funkcja
def calculate_heat_loss_from_ifc(ifc_path, indoor_temp, outdoor_temp, ach=0.5):
    model = ifcopenshell.open(ifc_path)
    spaces = model.by_type('IfcSpace')
    results = []
    shapes = []

    for space in spaces:
        if space.Representation:
            shape = ifcopenshell.geom.create_shape(settings, space).geometry
            shapes.append((space, shape))
            # print volume calcualted with shapely
        if not space.Representation:
            continue
        shape = ifcopenshell.geom.create_shape(settings, space)
        verts = shape.geometry.verts  # Lista współrzędnych w formacie [x1, y1, z1, x2, y2, z2, ...]
        faces = shape.geometry.faces  # Lista indeksów wierzchołków tworzących trójkąty
        volume = compute_volume(verts, faces)
        print(f"Pomieszczenie '{space.Name}' ma objętość: {volume:.3f} jednostek³")

            
    print(shapes)
    #     vol = get_space_volume(space)
    #     delta_t = indoor_temp - outdoor_temp
    #
    #     # graniczne elementy
    #     rels = model.get_inverse(space, 'IfcRelSpaceBoundary')
    #     u_vals, areas = [], []
    #     for rel in rels:
    #         if not rel.is_a('IfcRelSpaceBoundary'):
    #             continue
    #         if not getattr(rel, 'AccessBoundary', False) or getattr(rel, 'PhysicalOrVirtualBoundary', '').upper() != 'PHYSICAL':
    #             continue
    #         elem = rel.RelatingBuildingElement
    #         # obszar granicy
    #         area = getattr(rel, 'CalculatedArea', 0.0)
    #         if area <= 0:
    #             # spróbuj z PSet_WallCommon lub PSet_WindowCommon
    #             area = get_pset_value(elem, 'PSet_WallCommon', 'GrossArea') or get_pset_value(elem, 'PSet_WindowCommon', 'GrossArea')
    #         u = get_u_value(elem)
    #         u_vals.append(u)
    #         areas.append(area)
    #
    #     q_trans = transmission_heat_loss(u_vals, areas, delta_t)
    #     q_vent = ventilation_heat_loss(ach, vol, delta_t)
    #     q_total = q_trans + q_vent
    #
    #     results.append({
    #         'space_name': space.Name or space.GlobalId,
    #         'volume_m3': vol,
    #         'surface_m2': sum(areas),
    #         'Q_trans_W': q_trans,
    #         'Q_vent_W': q_vent,
    #         'Q_total_W': q_total
    #     })
    #
    # return results

if __name__ == '__main__':
    IFC_FILE = 'model.ifc'
    T_IN, T_OUT, ACH = 20.0, 0.0, 0.5
    data = calculate_heat_loss_from_ifc(IFC_FILE, T_IN, T_OUT, ACH)
    # for r in data:
    #     print(f"Przestrzeń: {r['space_name']}")
    #     print(f"  Kubatura: {r['volume_m3']:.1f} m3")
    #     print(f"  Pow. przegród: {r['surface_m2']:.1f} m2")
    #     print(f"  Straty transmisyjne: {r['Q_trans_W']:.2f} W")
    #     print(f"  Straty wentylacyjne: {r['Q_vent_W']:.2f} W")
    #     print(f"  Moc całkowita: {r['Q_total_W']:.2f} W\n")

