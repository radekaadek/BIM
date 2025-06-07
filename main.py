import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.geom
import tkinter as tk
from tkinter import filedialog, messagebox
from datetime import datetime, timedelta
from meteostat import Point, Daily, Hourly
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
from PIL import Image, ImageTk

# --- Constants ---
DENSITY_AIR = 1.2
SPECIFIC_HEAT_AIR = 1005

HEAT_GAIN_PEOPLE_SENSIBLE = 75
HEAT_GAIN_LIGHTING_DENSITY = 8
HEAT_GAIN_EQUIPMENT_DENSITY = 10

DEFAULT_WINDOW_SHGC = 0.70


# --- Helper functions for geometric calculations ---
def compute_volume_from_verts_faces(verts, faces):
    volume = 0.0
    for i in range(0, len(faces), 3):
        idx0, idx1, idx2 = faces[i], faces[i + 1], faces[i + 2]
        v0 = verts[idx0 * 3: idx0 * 3 + 3]
        v1 = verts[idx1 * 3: idx1 * 3 + 3]
        v2 = verts[idx2 * 3: idx2 * 3 + 3]
        volume += (v0[0] * (v1[1] * v2[2] - v1[2] * v2[1]) +
                   v1[0] * (v2[1] * v0[2] - v2[2] * v0[1]) +
                   v2[0] * (v0[1] * v1[2] - v0[2] * v1[1]))
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

    def _get_shgc_value(self, ifc_element):
        if not ifc_element:
            return DEFAULT_WINDOW_SHGC

        shgc = self._get_pset_value(ifc_element, 'PSet_WindowCommon', 'SolarHeatGainCoefficient')
        if shgc > 0:
            return shgc
        shgc = self._get_pset_value(ifc_element, 'PSet_WindowCommon', 'GValue')
        if shgc > 0:
            return shgc

        shgc = self._get_pset_value(ifc_element, 'PSet_GlazingCommon', 'SolarHeatGainCoefficient')
        if shgc > 0:
            return shgc

        if ifc_element.is_a('IfcWindow'):
            return DEFAULT_WINDOW_SHGC

        return 0.0

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
            total_surface_area = 0.0
            total_window_area = 0.0

            for rel in self.model.get_inverse(space):
                if rel.is_a('IfcRelSpaceBoundary') and getattr(rel, 'PhysicalOrVirtualBoundary',
                                                               'PHYSICAL') == 'PHYSICAL':
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
                    shgc = self._get_shgc_value(elem)

                    boundaries.append({'u_value': u, 'area': area, 'element_type': elem.is_a(), 'shgc': shgc})
                    total_surface_area += area

                    if 'Window' in elem.is_a():
                        total_window_area += area

            gross_floor_area = self._get_pset_value(space, 'PSet_SpaceCommon', 'GrossFloorArea')

            self.spaces_data.append({
                'space_name': space.Name or space.GlobalId,
                'volume_m3': vol,
                'boundaries': boundaries,
                'total_surface_m2': total_surface_area,
                'total_window_area_m2': total_window_area,
                'gross_floor_area_m2': gross_floor_area if gross_floor_area > 0 else (vol ** (2 / 3))
            })


# --- Heat Load Calculations ---
def _transmission_heat_loss(boundaries, delta_t):
    return sum(b['u_value'] * b['area'] for b in boundaries) * delta_t


def _ventilation_heat_loss(ach, volume, delta_t):
    v_dot = ach * volume / 3600.0
    m_dot = v_dot * DENSITY_AIR
    return m_dot * SPECIFIC_HEAT_AIR * delta_t


def _internal_heat_gains(occupancy_count, floor_area_m2, lighting_density_w_m2, equipment_density_w_m2):
    q_people = occupancy_count * HEAT_GAIN_PEOPLE_SENSIBLE
    q_lighting = floor_area_m2 * lighting_density_w_m2
    q_equipment = floor_area_m2 * equipment_density_w_m2
    return q_people + q_lighting + q_equipment


def _solar_heat_gain(window_area, shgc, solar_irradiance_w_m2):
    return window_area * shgc * solar_irradiance_w_m2


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


def calculate_total_cooling_load(
        model,
        indoor_temp,
        outdoor_temp,
        ach=0.5,
        occupancy_count_per_space=None,
        lighting_density_w_m2=HEAT_GAIN_LIGHTING_DENSITY,
        equipment_density_w_m2=HEAT_GAIN_EQUIPMENT_DENSITY,
        average_solar_irradiance_w_m2=200  # Default/simplified for yearly
):
    delta_t = outdoor_temp - indoor_temp

    results = []
    for s in model.spaces_data:
        space_name = s['space_name']
        volume = s['volume_m3']
        boundaries = s['boundaries']
        gross_floor_area = s['gross_floor_area_m2']
        total_window_area = s['total_window_area_m2']

        q_trans = _transmission_heat_loss(boundaries, delta_t) if delta_t > 0 else 0.0
        q_vent = _ventilation_heat_loss(ach, volume, delta_t) if delta_t > 0 else 0.0

        occupancy = occupancy_count_per_space.get(space_name, 0) if occupancy_count_per_space else 0
        q_internal = _internal_heat_gains(occupancy, gross_floor_area,
                                          lighting_density_w_m2, equipment_density_w_m2)

        q_solar = _solar_heat_gain(total_window_area, DEFAULT_WINDOW_SHGC, average_solar_irradiance_w_m2)

        q_tot = q_trans + q_vent + q_internal + q_solar

        results.append({'space_name': space_name, 'Q_total_W': max(0, q_tot)})

    total_building_cooling_load = sum(r['Q_total_W'] for r in results)
    results.append({'space_name': '--- TOTAL BUILDING ---', 'Q_total_W': total_building_cooling_load})
    return results


# --- Plotting Functions ---
def plot_yearly_energy_heating(model, lat, lon, alt, year, indoor_temp, ach, output_file):
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
    plt.title(f'Roczne zapotrzebowanie mocy na ogrzewanie {year} (wew.:{indoor_temp}°C, ACH:{ach})')
    plt.xlabel('Miesiąc')
    plt.ylabel('kWh/dzień')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    return f'Zapisano wykres rocznego zapotrzebowania na ogrzewanie: {output_file}'


def plot_daily_hourly_heating(model, lat, lon, alt, target_date, indoor_temp, ach, output_file):
    loc = Point(lat, lon, alt)
    start = datetime(target_date.year, target_date.month, target_date.day)
    end = start + timedelta(days=1)
    df = Hourly(loc, start=start, end=end).fetch().dropna(subset=['temp'])

    if df.empty:
        messagebox.showerror("Błąd danych pogodowych",
                             f"Brak danych godzinowych dla daty {target_date.strftime('%Y-%m-%d')} i lokalizacji.")
        return None  # Indicate failure

    times, losses = [], []
    for ts, row in df.iterrows():
        tm = row['temp']
        res = calculate_total_heat_loss(model, indoor_temp, tm, ach)
        total = next(r['Q_total_W'] for r in res if r['space_name'] == '--- TOTAL BUILDING ---')
        times.append(ts)
        losses.append(total)

    if not times:
        messagebox.showerror("Błąd danych",
                             f"Brak wystarczających danych do wygenerowania wykresu dla {target_date.strftime('%Y-%m-%d')}.")
        return None

    plt.figure(figsize=(12, 6))
    plt.plot(times, losses)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=2))
    plt.title(
        f"Godzinowe zapotrzebowanie mocy na ogrzewanie {target_date.strftime('%Y-%m-%d')} (wew.:{indoor_temp}°C, ACH:{ach})")
    plt.xlabel('Godzina')
    plt.ylabel('Moc [W]')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    return f'Zapisano godzinowy wykres zapotrzebowania na ogrzewanie: {output_file}'


def plot_yearly_energy_cooling(model, lat, lon, alt, year, indoor_temp, ach,
                               occupancy_count_per_space, lighting_density, equipment_density, output_file):
    start = datetime(year, 1, 1)
    end = datetime(year, 12, 31)
    loc = Point(lat, lon, alt)
    df = Daily(loc, start=start, end=end).fetch().dropna(subset=['tavg'])

    dates, energy = [], []
    for date, row in df.iterrows():
        tm = row['tavg']

        solar_irradiance_for_cooling_yearly = 200 if tm > indoor_temp else 0

        res = calculate_total_cooling_load(
            model, indoor_temp, tm, ach,
            occupancy_count_per_space=occupancy_count_per_space,
            lighting_density_w_m2=lighting_density,
            equipment_density_w_m2=equipment_density,
            average_solar_irradiance_w_m2=solar_irradiance_for_cooling_yearly
        )
        total = next(r['Q_total_W'] for r in res if r['space_name'] == '--- TOTAL BUILDING ---')
        dates.append(date)
        energy.append((total * 24) / 1000.0)  # kWh/day

    plt.figure(figsize=(12, 6))
    plt.plot(dates, energy, color='red')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator())
    plt.title(f'Roczne zapotrzebowanie mocy na chłodzenie {year} (wew.:{indoor_temp}°C, ACH:{ach})')
    plt.xlabel('Miesiąc')
    plt.ylabel('kWh/dzień')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    return f'Zapisano wykres rocznego zapotrzebowania na chłodzenie: {output_file}'


def plot_daily_hourly_cooling(model, lat, lon, alt, target_date, indoor_temp, ach,
                              occupancy_count_per_space, lighting_density, equipment_density, output_file):
    loc = Point(lat, lon, alt)
    start = datetime(target_date.year, target_date.month, target_date.day)
    end = start + timedelta(days=1)

    df = Hourly(loc, start=start, end=end).fetch()

    if df.empty:
        messagebox.showerror("Błąd danych pogodowych",
                             f"Brak danych godzinowych dla daty {target_date.strftime('%Y-%m-%d')} i lokalizacji.")
        return None  # Indicate failure

    df = df.dropna(subset=['temp'])
    df['swsf'] = df['swsf'].fillna(0)

    times, loads = [], []
    for ts, row in df.iterrows():
        tm = row['temp']
        solar_irradiance_w_m2 = row['swsf']

        res = calculate_total_cooling_load(
            model, indoor_temp, tm, ach,
            occupancy_count_per_space=occupancy_count_per_space,
            lighting_density_w_m2=lighting_density,
            equipment_density_w_m2=equipment_density,
            average_solar_irradiance_w_m2=solar_irradiance_w_m2
        )
        total = next(r['Q_total_W'] for r in res if r['space_name'] == '--- TOTAL BUILDING ---')
        times.append(ts)
        loads.append(total)

    if not times:
        messagebox.showerror("Błąd danych",
                             f"Brak wystarczających danych do wygenerowania wykresu dla {target_date.strftime('%Y-%m-%d')}.")
        return None

    plt.figure(figsize=(12, 6))
    plt.plot(times, loads, color='red')
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.gca().xaxis.set_major_locator(mdates.HourLocator(interval=2))
    plt.title(
        f"Godzinowe zapotrzebowanie mocy na chłodzenie {target_date.strftime('%Y-%m-%d')} (wew.:{indoor_temp}°C, ACH:{ach})")
    plt.xlabel('Godzina')
    plt.ylabel('Moc [W]')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    return f'Zapisano godzinowy wykres zapotrzebowania na chłodzenie: {output_file}'


# --- GUI Application ---
class HeatLossApp:
    def __init__(self, master):
        self.master = master
        master.title("Analiza Zapotrzebowania Energetycznego Budynku")

        # Set custom favicon
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            favicon_filename_ico = 'favicon.ico'
            favicon_filename_png = 'micon.png'
            icon_path_ico = os.path.join(script_dir, favicon_filename_ico)
            icon_path_png = os.path.join(script_dir, favicon_filename_png)

            if os.path.exists(icon_path_ico):
                self.master.iconbitmap(icon_path_ico)
            elif os.path.exists(icon_path_png):
                original_image = Image.open(icon_path_png)
                sizes = [(16, 16), (24, 24), (32, 32), (48, 48)]
                photo_images = []
                for size in sizes:
                    resized_image = original_image.copy()
                    resized_image.thumbnail(size, Image.LANCZOS)
                    photo_image = ImageTk.PhotoImage(resized_image)
                    photo_images.append(photo_image)
                self.master.iconphoto(False, *photo_images)
                self.master.icon_references = photo_images
            else:
                print(
                    f"Warning: No favicon file found (looked for {favicon_filename_ico} and {favicon_filename_png}). Using default icon.")
        except Exception as e:
            print(f"Error setting favicon: {e}")

        self.ifc_path = tk.StringVar()
        self.analysis_mode = tk.StringVar(value="Ogrzewanie")
        self.analysis_period = tk.StringVar(value="Roczny")
        self.analysis_input_value = tk.StringVar()
        self.lat = tk.DoubleVar(value=52.2297)
        self.lon = tk.DoubleVar(value=21.0122)
        self.alt = tk.IntVar(value=110)
        self.ach = tk.DoubleVar(value=0.5)
        self.indoor_temp = tk.DoubleVar(value=20.0)
        self.cooling_setpoint_temp = tk.DoubleVar(value=24.0)
        self.occupancy_default_count = tk.IntVar(value=1)
        self.lighting_density_w_m2 = tk.DoubleVar(value=HEAT_GAIN_LIGHTING_DENSITY)
        self.equipment_density_w_m2 = tk.DoubleVar(value=HEAT_GAIN_EQUIPMENT_DENSITY)

        self._create_widgets()
        self._on_analysis_option_change()  # Call initial update to set visibility

        self.master.update_idletasks()  # Force update to calculate minimum size
        # Optional: Set a minimum window size, allow it to expand if content requires more
        # min_width = self.master.winfo_reqwidth()
        # min_height = self.master.winfo_reqheight()
        # self.master.geometry(f"{max(600, min_width)}x{max(400, min_height)}")
        self.master.geometry("600x400")  # Attempt to set, but pack() may override

    def _create_widgets(self):
        input_frame = tk.Frame(self.master, padx=5, pady=2)
        input_frame.pack(side=tk.TOP, fill=tk.X)

        self._create_labeled_entry(input_frame, "Ścieżka do pliku IFC:", self.ifc_path, width=40, is_file_path=True,
                                   help_text="Ścieżka do pliku modelu budynku w formacie IFC.")

        self._create_labeled_entry(input_frame, "Szerokość geograficzna (Lat):", self.lat, width=10,
                                   help_text="Szerokość geograficzna lokalizacji budynku w stopniach dziesiętnych. Używana do pobierania danych pogodowych.")
        self._create_labeled_entry(input_frame, "Długość geograficzna (Lon):", self.lon, width=10,
                                   help_text="Długość geograficzna lokalizacji budynku w stopniach dziesiętnych. Używana do pobierania danych pogodowych.")
        self._create_labeled_entry(input_frame, "Wysokość n.p.m. [m] (Alt):", self.alt, width=10,
                                   help_text="Wysokość nad poziomem morza lokalizacji budynku w metrach. Używana do pobierania danych pogodowych.")

        self._create_labeled_entry(input_frame, "Wymiany powietrza [ACH]:", self.ach, width=10,
                                   help_text="Współczynnik wymiany powietrza na godzinę. Reprezentuje infiltrację powietrza z zewnątrz i/lub wentylację. Jest to uproszczone założenie, które nie uwzględnia złożonych systemów wentylacji ani dynamicznych warunków infiltracji.")

        # Heating Temperature Widgets (dynamically controlled)
        self.heating_temp_frame = self._create_labeled_entry(input_frame, "Temperatura wewnętrzna (ogrzewanie) [°C]:",
                                                             self.indoor_temp, width=10,
                                                             help_text="Żądana temperatura wewnętrzna do obliczeń zapotrzebowania na ciepło. Stała wartość przyjęta dla wszystkich pomieszczeń.",
                                                             return_frame=True)

        # Cooling Temperature Widgets (dynamically controlled)
        self.cooling_temp_frame = self._create_labeled_entry(input_frame, "Temperatura wewnętrzna (chłodzenie) [°C]:",
                                                             self.cooling_setpoint_temp, width=10,
                                                             help_text="Żądana temperatura wewnętrzna do obliczeń zapotrzebowania na chłód. Stała wartość przyjęta dla wszystkich pomieszczeń.",
                                                             return_frame=True)

        # Cooling Specific Inputs (dynamically controlled - will be hidden for heating)
        self.occupancy_frame = self._create_labeled_entry(input_frame, "Domyślna liczba osób na przestrzeń:",
                                                          self.occupancy_default_count, width=10,
                                                          help_text=f"Uproszczona liczba osób w każdej przestrzeni. Każda osoba generuje {HEAT_GAIN_PEOPLE_SENSIBLE} W ciepła jawnego.",
                                                          return_frame=True)
        self.lighting_frame = self._create_labeled_entry(input_frame, "Gęstość oświetlenia [W/m²]:",
                                                         self.lighting_density_w_m2, width=10,
                                                         help_text="Średnia gęstość mocy oświetlenia w Watach na metr kwadratowy powierzchni podłogi. Jest to uproszczone założenie dla zysków ciepła wewnętrznego.",
                                                         return_frame=True)
        self.equipment_frame = self._create_labeled_entry(input_frame, "Gęstość urządzeń [W/m²]:",
                                                          self.equipment_density_w_m2, width=10,
                                                          help_text="Średnia gęstość mocy urządzeń (np. komputerów, elektroniki) w Watach na metr kwadratowy powierzchni podłogi. Jest to uproszczone założenie dla zysków ciepła wewnętrznego.",
                                                          return_frame=True)

        # Analysis Options
        analysis_frame = tk.Frame(self.master, padx=5, pady=2)
        analysis_frame.pack(side=tk.TOP, fill=tk.X)

        # --- Analysis Mode Row ---
        mode_option_frame = tk.Frame(analysis_frame)
        mode_option_frame.pack(fill=tk.X, pady=1, anchor="w")
        tk.Label(mode_option_frame, text="Tryb analizy:").pack(side=tk.LEFT, padx=5)
        mode_menu = tk.OptionMenu(mode_option_frame, self.analysis_mode, "Ogrzewanie", "Chłodzenie",
                                  command=self._on_analysis_option_change)
        mode_menu.pack(side=tk.LEFT, padx=5)
        self._add_help_button_pack(mode_option_frame,
                                   "Wybierz, czy przeprowadzasz analizę ogrzewania, czy chłodzenia. Ma to wpływ na rodzaj pobieranych danych pogodowych i uwzględniane obciążenia wewnętrzne/słoneczne.")

        # --- Analysis Period Row ---
        period_option_frame = tk.Frame(analysis_frame)
        period_option_frame.pack(fill=tk.X, pady=1, anchor="w")
        tk.Label(period_option_frame, text="Zakres czasowy:").pack(side=tk.LEFT, padx=5)
        period_menu = tk.OptionMenu(period_option_frame, self.analysis_period, "Roczny", "Dzienny",
                                    command=self._on_analysis_option_change)
        period_menu.pack(side=tk.LEFT, padx=5)
        self._add_help_button_pack(period_option_frame,
                                   "Wybierz zakres czasowy analizy. **Roczna** analiza wykorzystuje średnie dobowe temperatury i daje zużycie kWh/dzień w ciągu roku. **Dzienna** analiza wykorzystuje godzinowe temperatury i promieniowanie słoneczne dla wybranej daty, przedstawiając moc [W] w ciągu dnia.")

        # Dynamic Analysis Input
        dynamic_input_frame = tk.Frame(self.master, padx=5, pady=2)
        dynamic_input_frame.pack(side=tk.TOP, fill=tk.X)

        self.analysis_input_label = tk.Label(dynamic_input_frame, text="Wprowadź rok:")
        self.analysis_input_label.pack(side=tk.LEFT, padx=5)
        self.analysis_input_entry = tk.Entry(dynamic_input_frame, textvariable=self.analysis_input_value, width=15)
        self.analysis_input_entry.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
        self.analysis_input_help_button = self._add_help_button_pack(dynamic_input_frame,
                                                                     "Wprowadź rok (np. 2023) dla analizy rocznej lub datę (np. 2023-01-15) dla analizy dziennej, zgodnie z wybranym trybem i zakresem.",
                                                                     return_button=True)

        # Buttons
        button_frame = tk.Frame(self.master, pady=5)
        button_frame.pack(side=tk.TOP, fill=tk.X)
        tk.Button(button_frame, text="Generuj wykresy", command=self._run_analysis).pack(pady=2)

        # Status Message
        self.status_label = tk.Label(self.master, text="", fg="blue")
        self.status_label.pack(side=tk.BOTTOM, fill=tk.X, padx=5, pady=2)

    def _create_labeled_entry(self, parent_frame, label_text, textvariable, width=20, is_file_path=False, help_text="",
                              return_frame=False):
        frame = tk.Frame(parent_frame)
        frame.pack(fill=tk.X, pady=1, padx=5)

        label = tk.Label(frame, text=label_text)
        label.pack(side=tk.LEFT, anchor="w", padx=(0, 5))

        entry = tk.Entry(frame, textvariable=textvariable, width=width)
        entry.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=(0, 5))

        if is_file_path:
            tk.Button(frame, text="Wybierz plik", command=self._browse_ifc).pack(side=tk.LEFT, padx=(0, 0))

        self._add_help_button_pack(frame, help_text)

        if return_frame:
            return frame

    def _add_help_button_pack(self, parent_widget, text, return_button=False):
        button = tk.Button(parent_widget, text="?", width=2, command=lambda: messagebox.showinfo("Pomoc", text))
        button.pack(side=tk.LEFT, padx=(5, 0))
        if return_button:
            return button

    def _browse_ifc(self):
        file_path = filedialog.askopenfilename(
            filetypes=[("IFC files", "*.ifc"), ("All files", "*.*")]
        )
        if file_path:
            self.ifc_path.set(file_path)

    def _on_analysis_option_change(self, *args):
        # Update dynamic input label based on period
        selected_period = self.analysis_period.get()
        if selected_period == "Roczny":
            self.analysis_input_label.config(text="Wprowadź rok (np. 2023):")
            self.analysis_input_value.set("")
        elif selected_period == "Dzienny":
            self.analysis_input_label.config(text="Wprowadź datę (RRRR-MM-DD):")
            self.analysis_input_value.set("")

        # Update temperature and internal gains input visibility based on mode
        selected_mode = self.analysis_mode.get()
        if selected_mode == "Ogrzewanie":
            self.heating_temp_frame.pack(fill=tk.X, pady=1, padx=5)
            self.cooling_temp_frame.pack_forget()

            # Hide fields related to internal gains for heating
            self.occupancy_frame.pack_forget()
            self.lighting_frame.pack_forget()
            self.equipment_frame.pack_forget()

        elif selected_mode == "Chłodzenie":
            self.heating_temp_frame.pack_forget()
            self.cooling_temp_frame.pack(fill=tk.X, pady=1, padx=5)

            # Show fields related to internal gains for cooling
            # Re-pack them to ensure they appear in the correct order
            self.occupancy_frame.pack(fill=tk.X, pady=1, padx=5)
            self.lighting_frame.pack(fill=tk.X, pady=1, padx=5)
            self.equipment_frame.pack(fill=tk.X, pady=1, padx=5)

    def _run_analysis(self):
        ifc_path = self.ifc_path.get()
        if not ifc_path:
            messagebox.showerror("Błąd", "Proszę wybrać plik IFC.")
            return

        try:
            building = BuildingModel(ifc_path)
            if not building.spaces_data:
                messagebox.showerror("Błąd",
                                     "Nie przetworzono żadnych przestrzeni z pliku IFC. Upewnij się, że plik IFC zawiera elementy IfcSpace z poprawną geometrią lub właściwościami.")
                return

            lat = self.lat.get()
            lon = self.lon.get()
            alt = self.alt.get()
            ach = self.ach.get()

            analysis_mode = self.analysis_mode.get()
            analysis_period = self.analysis_period.get()
            analysis_input = self.analysis_input_value.get()

            if not analysis_input:
                messagebox.showerror("Błąd", f"Proszę wprowadzić wartość dla wybranej analizy.")
                return

            occupancy_per_space = {
                space['space_name']: self.occupancy_default_count.get()
                for space in building.spaces_data
            }
            lighting_density = self.lighting_density_w_m2.get()
            equipment_density = self.equipment_density_w_m2.get()

            self.status_label.config(text="Rozpoczynanie analizy...", fg="blue")
            self.master.update_idletasks()

            if analysis_mode == "Ogrzewanie":
                indoor_temp = self.indoor_temp.get()
                if analysis_period == "Roczny":
                    try:
                        year = int(analysis_input)
                        output_file = f'roczne_zapotrzebowanie_ogrzewanie_{year}.png'
                        msg = plot_yearly_energy_heating(building, lat, lon, alt, year, indoor_temp, ach, output_file)
                        if msg:
                            self.status_label.config(text=msg, fg="green")
                        else:
                            self.status_label.config(text="Błąd podczas generowania wykresu.", fg="red")

                    except ValueError:
                        messagebox.showerror("Błąd", "Nieprawidłowy format roku. Wprowadź liczbę całkowitą (np. 2023).")
                        self.status_label.config(text="Błąd wejścia.", fg="red")
                    except Exception as e:
                        messagebox.showerror("Błąd podczas analizy rocznej ogrzewania", str(e))
                        self.status_label.config(text="Błąd podczas analizy rocznej ogrzewania.", fg="red")
                elif analysis_period == "Dzienny":
                    try:
                        target_date = datetime.strptime(analysis_input, '%Y-%m-%d')
                        output_file = f'godzinowe_zapotrzebowanie_ogrzewanie_{target_date.strftime("%Y-%m-%d")}.png'
                        msg = plot_daily_hourly_heating(building, lat, lon, alt, target_date, indoor_temp, ach,
                                                        output_file)
                        if msg:
                            self.status_label.config(text=msg, fg="green")
                        else:
                            self.status_label.config(text="Błąd podczas generowania wykresu.", fg="red")
                    except ValueError:
                        messagebox.showerror("Błąd", "Nieprawidłowy format daty. Użyj RRRR-MM-DD (np. 2023-01-15).")
                        self.status_label.config(text="Błąd wejścia.", fg="red")
                    except Exception as e:
                        messagebox.showerror("Błąd podczas analizy dziennej ogrzewania", str(e))
                        self.status_label.config(text="Błąd podczas analizy dziennej ogrzewania.", fg="red")

            elif analysis_mode == "Chłodzenie":
                indoor_temp = self.cooling_setpoint_temp.get()
                if analysis_period == "Roczny":
                    try:
                        year = int(analysis_input)
                        output_file = f'roczne_zapotrzebowanie_chlodzenie_{year}.png'
                        msg = plot_yearly_energy_cooling(building, lat, lon, alt, year, indoor_temp, ach,
                                                         occupancy_per_space, lighting_density, equipment_density,
                                                         output_file)
                        if msg:
                            self.status_label.config(text=msg, fg="green")
                        else:
                            self.status_label.config(text="Błąd podczas generowania wykresu.", fg="red")
                    except ValueError:
                        messagebox.showerror("Błąd", "Nieprawidłowy format roku. Wprowadź liczbę całkowitą (np. 2023).")
                        self.status_label.config(text="Błąd wejścia.", fg="red")
                    except Exception as e:
                        messagebox.showerror("Błąd podczas analizy rocznej chłodzenia", str(e))
                        self.status_label.config(text="Błąd podczas analizy rocznej chłodzenia.", fg="red")
                elif analysis_period == "Dzienny":
                    try:
                        target_date = datetime.strptime(analysis_input, '%Y-%m-%d')
                        output_file = f'godzinowe_zapotrzebowanie_chlodzenie_{target_date.strftime("%Y-%m-%d")}.png'
                        msg = plot_daily_hourly_cooling(building, lat, lon, alt, target_date, indoor_temp, ach,
                                                        occupancy_per_space, lighting_density, equipment_density,
                                                        output_file)
                        if msg:
                            self.status_label.config(text=msg, fg="green")
                        else:
                            self.status_label.config(text="Błąd podczas generowania wykresu.", fg="red")
                    except ValueError:
                        messagebox.showerror("Błąd", "Nieprawidłowy format daty. Użyj RRRR-MM-DD (np. 2023-07-15).")
                        self.status_label.config(text="Błąd wejścia.", fg="red")
                    except Exception as e:
                        messagebox.showerror("Błąd podczas analizy dziennej chłodzenia", str(e))
                        self.status_label.config(text="Błąd podczas analizy dziennej chłodzenia.", fg="red")

        except Exception as e:
            messagebox.showerror("Ogólny błąd", f"Wystąpił nieoczekiwany błąd: {str(e)}")
            self.status_label.config(text="Wystąpił ogólny błąd.", fg="red")


if __name__ == '__main__':
    root = tk.Tk()
    app = HeatLossApp(root)
    root.mainloop()
