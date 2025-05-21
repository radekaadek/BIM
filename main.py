import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.geom
import click
from datetime import datetime, timedelta
from meteostat import Point, Daily, Hourly
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# --- Constants ---
DENSITY_AIR = 1.2  # kg/m3
SPECIFIC_HEAT_AIR = 1005  # J/(kg·K)

# --- Helper functions for geometric calculations ---
def compute_volume_from_verts_faces(verts, faces):
    volume = 0.0
    for i in range(0, len(faces), 3):
        idx0, idx1, idx2 = faces[i], faces[i+1], faces[i+2]
        v0 = verts[idx0*3 : idx0*3+3]
        v1 = verts[idx1*3 : idx1*3+3]
        v2 = verts[idx2*3 : idx2*3+3]
        volume += (v0[0]*(v1[1]*v2[2] - v1[2]*v2[1]) +
                   v1[0]*(v2[1]*v0[2] - v2[2]*v0[1]) +
                   v2[0]*(v0[1]*v1[2] - v0[2]*v1[1]))
    return abs(volume / 6.0)

# --- BuildingModel Class ---
class BuildingModel:
    def __init__(self, ifc_path):
        self.ifc_path = ifc_path
        self.model = ifcopenshell.open(ifc_path)
        self.settings = ifcopenshell.geom.settings()
        self.spaces_data = []
        self._preprocess_spaces()

    def _get_pset_value(self, ifc_elem, pset_name, prop_name):
        if not ifc_elem:
            return 0.0
        psets = ifcopenshell.util.element.get_psets(ifc_elem)
        pset = psets.get(pset_name, {})
        value = pset.get(prop_name)
        if isinstance(value, (int, float)):
            return float(value)
        if hasattr(value, 'wrappedValue') and isinstance(value.wrappedValue, (int, float)):
            return float(value.wrappedValue)
        return 0.0

    def _get_u_value(self, ifc_element):
        if not ifc_element:
            return 0.35
        u = self._get_pset_value(ifc_element, 'PSet_ThermalTransmittance', 'ThermalTransmittance')
        if u > 0:
            return u
        u = self._get_pset_value(ifc_element, 'PSet_ThermalTransmittance', 'UValue')
        if u > 0:
            return u
        mapping = {
            'IfcWall': 0.3, 'IfcWallStandardCase': 0.3, 'IfcWindow': 1.1,
            'IfcRoof': 0.25, 'IfcSlab': 0.4, 'IfcDoor': 1.5, 'IfcCovering': 0.35
        }
        for pset, prop in {
            'PSet_WallCommon': 'ThermalTransmittance', 'PSet_WindowCommon': 'ThermalTransmittance',
            'PSet_DoorCommon': 'ThermalTransmittance', 'PSet_RoofCommon': 'ThermalTransmittance',
            'PSet_SlabCommon': 'ThermalTransmittance'
        }.items():
            u = self._get_pset_value(ifc_element, pset, prop)
            if u > 0:
                return u
        return mapping.get(ifc_element.is_a(), 0.35)

    def _get_space_volume(self, space_entity):
        vol = getattr(space_entity, 'GrossVolume', None)
        if vol and vol > 0:
            return float(vol)
        vol = getattr(space_entity, 'NetVolume', None)
        if vol and vol > 0:
            return float(vol)
        for pset in ['PSet_SpaceCommon', 'PSet_SpaceThermalDesign']:
            vol = self._get_pset_value(space_entity, pset, 'GrossVolume')
            if vol > 0:
                return vol
            vol = self._get_pset_value(space_entity, pset, 'NetVolume')
            if vol > 0:
                return vol
        if space_entity.Representation:
            try:
                shape = ifcopenshell.geom.create_shape(self.settings, space_entity)
                verts, faces = shape.geometry.verts, shape.geometry.faces
                if verts and faces:
                    return compute_volume_from_verts_faces(verts, faces)
            except Exception:
                pass
        return 0.0

    def _preprocess_spaces(self):
        ifc_spaces = self.model.by_type('IfcSpace')
        if not ifc_spaces:
            print('ERROR: No IfcSpace entities found.')
            return
        for space in ifc_spaces:
            vol = self._get_space_volume(space)
            if vol <= 0:
                continue
            boundaries = []
            total_area = 0.0
            for rel in self.model.get_inverse(space):
                if rel.is_a('IfcRelSpaceBoundary') and getattr(rel, 'PhysicalOrVirtualBoundary', 'PHYSICAL') == 'PHYSICAL':
                    elem = rel.RelatedBuildingElement
                    if not elem:
                        continue
                    area = float(getattr(rel, 'CalculatedArea', 0.0))
                    if area <= 0:
                        for pset, prop in {
                            'PSet_WallCommon': 'GrossArea', 'PSet_WindowCommon': 'Area',
                            'PSet_DoorCommon': 'Area', 'PSet_RoofCommon': 'GrossArea',
                            'PSet_SlabCommon': 'GrossArea'
                        }.items():
                            val = self._get_pset_value(elem, pset, prop)
                            if val > 0:
                                area = val
                                break
                    if area <= 0:
                        continue
                    u = self._get_u_value(elem)
                    boundaries.append({'u_value': u, 'area': area})
                    total_area += area
            self.spaces_data.append({
                'space_name': space.Name or space.GlobalId,
                'volume_m3': vol,
                'boundaries': boundaries,
                'total_surface_m2': total_area
            })

# --- Heat Loss Calculations ---
def _transmission_heat_loss(boundaries, delta_t):
    return sum(b['u_value'] * b['area'] for b in boundaries) * delta_t

def _ventilation_heat_loss(ach, volume, delta_t):
    v_dot = ach * volume / 3600.0
    m_dot = v_dot * DENSITY_AIR
    return m_dot * SPECIFIC_HEAT_AIR * delta_t

def calculate_total_heat_loss(model, indoor_temp, outdoor_temp, ach=0.5):
    delta_t = indoor_temp - outdoor_temp
    results = []
    for s in model.spaces_data:
        q_trans = _transmission_heat_loss(s['boundaries'], delta_t)
        q_vent = _ventilation_heat_loss(ach, s['volume_m3'], delta_t)
        q_tot = (q_trans + q_vent) if delta_t > 0 else 0.0
        results.append({'space_name': s['space_name'], 'Q_total_W': max(0, q_tot)})
    total = sum(r['Q_total_W'] for r in results)
    results.append({'space_name': '--- TOTAL BUILDING ---', 'Q_total_W': total})
    return results

# --- Plotting Functions ---
def plot_yearly_energy(model, lat, lon, alt, year, indoor_temp, ach, output_file):
    start = datetime(year, 1, 1)
    end = datetime(year, 12, 31)
    loc = Point(lat, lon, alt)
    df = Daily(loc, start=start, end=end).fetch().dropna(subset=['tavg'])
    dates, energy = [], []
    for date, row in df.iterrows():
        tm = row['tavg']
        res = calculate_total_heat_loss(model, indoor_temp, tm, ach)
        total = next(r['Q_total_W'] for r in res if r['space_name'] == '--- TOTAL BUILDING ---')
        dates.append(date)
        energy.append((total * 24) / 1000.0)
    plt.figure(figsize=(12, 6))
    plt.plot(dates, energy)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator())
    plt.title(f'Roczne zużycie energii {year} (wew.:{indoor_temp}°C, ACH:{ach})')
    plt.xlabel('Miesiąc')
    plt.ylabel('kWh/dzień')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_file)
    print(f'Zapisano wykres rocznego zużycia: {output_file}')

def plot_daily_hourly(model, lat, lon, alt, target_date, indoor_temp, ach, output_file):
    loc = Point(lat, lon, alt)
    start = datetime(target_date.year, target_date.month, target_date.day)
    end = start + timedelta(days=1)
    df = Hourly(loc, start=start, end=end).fetch().dropna(subset=['temp'])
    times, losses = [], []
    for ts, row in df.iterrows():
        tm = row['temp']
        res = calculate_total_heat_loss(model, indoor_temp, tm, ach)
        total = next(r['Q_total_W'] for r in res if r['space_name'] == '--- TOTAL BUILDING ---')
        times.append(ts)
        losses.append(total)
    plt.figure(figsize=(12, 6))
    plt.plot(times, losses)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=2))
    plt.title(f"Godzinowe straty ciepła {target_date.strftime('%Y-%m-%d')} (wew.:{indoor_temp}°C, ACH:{ach})")
    plt.xlabel('Godzina')
    plt.ylabel('Moc strat [W]')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_file)
    print(f'Zapisano godzinowy wykres strat ciepła: {output_file}')

# --- CLI with Click ---
@click.command()
@click.option('--ifc', 'ifc_path', type=click.Path(exists=True), required=True, help='Ścieżka do pliku IFC')
@click.option('--year', type=click.INT, help='Rok do analizy rocznej')
@click.option('--date', 'date_str', type=click.STRING, help='Data do analizy dziennej (YYYY-MM-DD)')
@click.option('--lat', type=click.FLOAT, default=52.2297, help='Szerokość geogr.')
@click.option('--lon', type=click.FLOAT, default=21.0122, help='Długość geogr.')
@click.option('--alt', type=click.INT, default=110, help='Wysokość npm [m]')
@click.option('--ach', type=click.FLOAT, default=0.5, help='Wymiany powietrza [ACH]')
@click.option('--tmp', 'indoor_temp', type=click.FLOAT, default=20.0, help='Temperatura wewnętrzna [°C]')
def main(ifc_path, year, date_str, lat, lon, alt, ach, indoor_temp):
    """Oblicz i generuj wykresy strat ciepła"""
    building = BuildingModel(ifc_path)
    if not building.spaces_data:
        click.echo('Błąd: Nie przetworzono żadnych przestrzeni. Kończę.', err=True)
        raise SystemExit(1)

    if year:
        output_year = f'roczne_zuzycie_{year}.png'
        plot_yearly_energy(building, lat, lon, alt, year, indoor_temp, ach, output_year)

    if date_str:
        try:
            target_date = datetime.strptime(date_str, '%Y-%m-%d')
        except ValueError:
            click.echo('Nieprawidłowy format daty. Użyj YYYY-MM-DD.', err=True)
            raise SystemExit(1)
        output_day = f'hourly_loss_{target_date.strftime("%Y-%m-%d")}.png'
        plot_daily_hourly(building, lat, lon, alt, target_date, indoor_temp, ach, output_day)

    if not year and not date_str:
        click.echo('Podaj przynajmniej --year lub --date aby wygenerować wykres.', err=True)
        raise SystemExit(1)

if __name__ == '__main__':
    main()

