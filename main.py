import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.geom
import numpy as np
from datetime import datetime, timedelta
from meteostat import Point, Daily
import matplotlib.pyplot as plt # Import matplotlib
import matplotlib.dates as mdates # For formatting dates on the plot axis

# --- Constants ---
DENSITY_AIR = 1.2  # kg/m3
SPECIFIC_HEAT_AIR = 1005  # J/(kg·K)

# --- Helper functions for geometric calculations (can be outside class) ---
def compute_volume_from_verts_faces(verts, faces):
    """Computes the volume of a mesh given its vertices and faces."""
    volume = 0.0
    # Iterate over faces, assuming they are triangles
    for i in range(0, len(faces), 3):
        # Get vertex indices for the current triangle
        idx0, idx1, idx2 = faces[i], faces[i+1], faces[i+2]
        
        # Get coordinates of the vertices
        # Vertex coordinates are stored sequentially in verts: [x0,y0,z0, x1,y1,z1, ...]
        v0 = verts[idx0*3 : idx0*3+3]
        v1 = verts[idx1*3 : idx1*3+3]
        v2 = verts[idx2*3 : idx2*3+3]

        # Calculate the signed volume of the tetrahedron formed by the triangle and the origin (0,0,0)
        volume += (v0[0]*(v1[1]*v2[2] - v1[2]*v2[1]) +
                   v1[0]*(v2[1]*v0[2] - v2[2]*v0[1]) +
                   v2[0]*(v0[1]*v1[2] - v0[2]*v1[1]))
    return abs(volume / 6.0)

# --- BuildingModel Class ---
class BuildingModel:
    """
    Represents a building model loaded from an IFC file,
    with preprocessed data for heat loss calculations.
    """
    def __init__(self, ifc_path):
        """
        Initializes the BuildingModel by loading and preprocessing the IFC file.
        
        Args:
            ifc_path (str): Path to the IFC file.
        """
        self.ifc_path = ifc_path
        self.model = ifcopenshell.open(ifc_path)
        self.settings = ifcopenshell.geom.settings()
        # Example: self.settings.set(self.settings.USE_PYTHON_OPENCASCADE, True) # If needed
        
        self.spaces_data = [] # Stores preprocessed data for each IfcSpace
        self._preprocess_spaces()

    def _get_pset_value(self, ifc_elem, pset_name, prop_name):
        """Helper to get a numeric value from a PropertySet."""
        if not ifc_elem:
            return 0.0
        psets = ifcopenshell.util.element.get_psets(ifc_elem)
        pset = psets.get(pset_name, {})
        value = pset.get(prop_name)
        if isinstance(value, (int, float)):
            return float(value)
        # Handle IfcValue types if necessary, e.g. value.wrappedValue
        if hasattr(value, 'wrappedValue'):
             return float(value.wrappedValue) if isinstance(value.wrappedValue, (int, float)) else 0.0
        return 0.0


    def _get_u_value(self, ifc_element):
        """Gets U-value for a building element."""
        if not ifc_element:
            return 0.35 # Default U-value if element is None
            
        # 1. Try PSet_ThermalTransmittance
        u = self._get_pset_value(ifc_element, 'PSet_ThermalTransmittance', 'ThermalTransmittance') # 'U' or 'ThermalTransmittance'
        if u > 0:
            return u
        u = self._get_pset_value(ifc_element, 'PSet_ThermalTransmittance', 'UValue') # Common alternative name
        if u > 0:
            return u

        # 2. Fallback to predefined types (can be expanded)
        element_type_name = ifc_element.is_a() # e.g., 'IfcWall', 'IfcWindow'
        
        mapping = {
            'IfcWall': 0.3, 'IfcWallStandardCase': 0.3,
            'IfcWindow': 1.1,
            'IfcRoof': 0.25,
            'IfcSlab': 0.4,  # For slabs (floors, ceilings if external)
            'IfcDoor': 1.5,  # Example for doors
            'IfcCovering': 0.35, # Generic covering
        }
        
        # Check common PSet names for U-value if specific PSet_ThermalTransmittance is missing
        common_psets_u_props = {
            'PSet_WallCommon': 'ThermalTransmittance',
            'PSet_WindowCommon': 'ThermalTransmittance',
            'PSet_DoorCommon': 'ThermalTransmittance',
            'PSet_RoofCommon': 'ThermalTransmittance',
            'PSet_SlabCommon': 'ThermalTransmittance',
        }
        for pset_name, prop_name in common_psets_u_props.items():
            u = self._get_pset_value(ifc_element, pset_name, prop_name)
            if u > 0: return u

        return mapping.get(element_type_name, 0.35) # Default if no specific rule

    def _get_space_volume(self, space_entity):
        """Calculates or retrieves the volume of an IfcSpace."""
        # 1. Try direct attributes
        vol = getattr(space_entity, 'GrossVolume', None)
        if vol and vol > 0: return float(vol)
        vol = getattr(space_entity, 'NetVolume', None)
        if vol and vol > 0: return float(vol)

        # 2. Try PSet_SpaceCommon or PSet_SpaceThermalProperties
        vol = self._get_pset_value(space_entity, 'PSet_SpaceCommon', 'GrossVolume')
        if vol > 0: return vol
        vol = self._get_pset_value(space_entity, 'PSet_SpaceCommon', 'NetVolume')
        if vol > 0: return vol
        vol = self._get_pset_value(space_entity, 'PSet_SpaceThermalProperties', 'GrossVolume')
        if vol > 0: return vol

        # 3. Fallback to geometric calculation if representation exists
        if space_entity.Representation:
            try:
                shape = ifcopenshell.geom.create_shape(self.settings, space_entity)
                if shape and hasattr(shape.geometry, 'verts') and hasattr(shape.geometry, 'faces'):
                    verts = shape.geometry.verts
                    faces = shape.geometry.faces
                    if verts and faces and len(verts) > 0 and len(faces) > 0:
                        return compute_volume_from_verts_faces(verts, faces)
                    else:
                        print(f"Warning: Space {space_entity.Name or space_entity.GlobalId} has geometry but verts/faces are empty.")
                else:
                    print(f"Warning: Could not generate valid shape or geometry for volume calculation of space {space_entity.Name or space_entity.GlobalId}.")
            except Exception as e:
                print(f"Error computing volume for space {space_entity.Name or space_entity.GlobalId} from geometry: {e}")
        
        print(f"Warning: Volume for space {space_entity.Name or space_entity.GlobalId} could not be determined, defaulting to 0.")
        return 0.0

    def _preprocess_spaces(self):
        """
        Extracts and stores relevant data for each IfcSpace in the model.
        This includes volume and properties of boundary elements (U-value, area).
        """
        ifc_spaces = self.model.by_type('IfcSpace')
        if not ifc_spaces:
            print("No IfcSpace entities found in the model.")
            return

        for space_entity in ifc_spaces:
            space_name = space_entity.Name or space_entity.LongName or space_entity.GlobalId
            volume = self._get_space_volume(space_entity)

            if volume <= 0:
                print(f"Skipping space '{space_name}' due to zero or invalid volume ({volume:.2f} m³).")
                continue

            boundary_elements_data = []
            total_boundary_surface_m2 = 0.0
            
            rels = self.model.get_inverse(space_entity) 
            space_boundaries = [rel for rel in rels if rel.is_a('IfcRelSpaceBoundary')]

            for rel_sb in space_boundaries:
                physical_boundary = getattr(rel_sb, 'PhysicalOrVirtualBoundary', 'PHYSICAL') 
                if physical_boundary != 'PHYSICAL':
                    continue

                related_element = rel_sb.RelatedBuildingElement
                if not related_element: 
                    continue
                
                area = 0.0
                if hasattr(rel_sb, 'ConnectionGeometry') and rel_sb.ConnectionGeometry:
                    # Attempt to get area from ConnectionGeometry if available and represents a surface
                    # This is a more complex part, for now, we'll rely on CalculatedArea or PSet
                    # For simplicity, we'll prioritize CalculatedArea and PSet values.
                    # If you need more precise area from geometry, this part would need expansion.
                    pass

                if hasattr(rel_sb, 'CalculatedArea') and rel_sb.CalculatedArea is not None:
                    area = float(rel_sb.CalculatedArea)
                
                if area <= 0: # If CalculatedArea is not available or zero, try PSets
                    # Try to get area from PSet of the related element
                    # This is a simplified approach; specific PSets and properties might vary
                    psets = ifcopenshell.util.element.get_psets(related_element)
                    area_psets = { # More specific PSet names based on common practice
                        'PSet_WallCommon': 'GrossArea', 'PSet_WindowCommon': 'Area', 
                        'PSet_DoorCommon': 'Area', 'PSet_RoofCommon': 'GrossArea',
                        'PSet_SlabCommon': 'GrossArea'
                        # Add other relevant Psets and properties for area
                    }
                    for pset_name, prop_name in area_psets.items():
                        elem_area = self._get_pset_value(related_element, pset_name, prop_name)
                        if elem_area > 0:
                            area = elem_area
                            break # Found an area
                
                if area <= 0: # If still no area, skip this boundary element
                    # print(f"Warning: Could not determine area for boundary {related_element.Name if related_element else 'N/A'} of space {space_name}")
                    continue 

                u_value = self._get_u_value(related_element)
                
                boundary_elements_data.append({'u_value': u_value, 'area': area, 'element_type': related_element.is_a()})
                total_boundary_surface_m2 += area
            
            if not boundary_elements_data:
                print(f"Space '{space_name}' has volume but no valid physical boundaries found or processed. Heat transmission will be zero.")

            self.spaces_data.append({
                'space_name': space_name,
                'volume_m3': volume,
                'boundaries': boundary_elements_data,
                'total_boundary_surface_m2': total_boundary_surface_m2
            })
        
        if not self.spaces_data:
            print("No spaces were successfully processed for heat loss calculation.")


# --- Heat Loss Calculation Functions (using preprocessed data) ---
def _transmission_heat_loss(space_boundary_data, delta_t):
    """Calculates transmission heat loss for a single space."""
    q_trans = 0.0
    for boundary in space_boundary_data:
        q_trans += boundary['u_value'] * boundary['area']
    return q_trans * delta_t 

def _ventilation_heat_loss(ach, volume, delta_t):
    """Calculates ventilation heat loss."""
    if volume <= 0: return 0.0
    v_dot_air = (ach * volume) / 3600.0  
    m_dot_air = v_dot_air * DENSITY_AIR
    q_vent = m_dot_air * SPECIFIC_HEAT_AIR * delta_t
    return q_vent

def calculate_total_heat_loss(building_model: BuildingModel, indoor_temp: float, outdoor_temp: float, ach: float = 0.5):
    """
    Calculates heat loss for each space in the building model given temperatures.
    """
    results = []
    if not building_model.spaces_data:
        print("Cannot calculate heat loss: No space data available in the building model.")
        return results

    delta_t = indoor_temp - outdoor_temp

    if delta_t <= 0:
        if indoor_temp < outdoor_temp:
            # print(f"Indoor temperature ({indoor_temp}°C) is lower than outdoor temperature ({outdoor_temp}°C). Heat GAIN occurs. Heat loss set to 0.")
            pass # Silently handle, Q will be 0
        elif indoor_temp == outdoor_temp:
            # print(f"Indoor temperature ({indoor_temp}°C) is equal to outdoor temperature ({outdoor_temp}°C). No temperature difference, heat loss is zero.")
            pass # Silently handle, Q will be 0
        
        # Create results structure with 0 heat loss
        for space_data in building_model.spaces_data:
            results.append({
                'space_name': space_data['space_name'],
                'volume_m3': space_data['volume_m3'],
                'total_boundary_surface_m2': space_data['total_boundary_surface_m2'],
                'Q_trans_W': 0.0,
                'Q_vent_W': 0.0,
                'Q_total_W': 0.0,
                'delta_t_K': delta_t
            })
        if results: # Add total building summary
            results.append({
                'space_name': '--- TOTAL BUILDING ---',
                'volume_m3': sum(s['volume_m3'] for s in building_model.spaces_data if 'volume_m3' in s),
                'total_boundary_surface_m2': sum(s['total_boundary_surface_m2'] for s in building_model.spaces_data if 'total_boundary_surface_m2' in s),
                'Q_trans_W': 0.0,
                'Q_vent_W': 0.0,
                'Q_total_W': 0.0,
                'delta_t_K': delta_t
            })
        return results


    total_building_q_trans = 0
    total_building_q_vent = 0
    total_building_q_total = 0

    for space_data in building_model.spaces_data:
        q_trans = _transmission_heat_loss(space_data['boundaries'], delta_t)
        q_vent = _ventilation_heat_loss(ach, space_data['volume_m3'], delta_t)
        q_total = q_trans + q_vent

        results.append({
            'space_name': space_data['space_name'],
            'volume_m3': space_data['volume_m3'],
            'total_boundary_surface_m2': space_data['total_boundary_surface_m2'],
            'Q_trans_W': q_trans,
            'Q_vent_W': q_vent,
            'Q_total_W': q_total,
            'delta_t_K': delta_t
        })
        total_building_q_trans += q_trans
        total_building_q_vent += q_vent
        total_building_q_total += q_total
    
    if results: 
        results.append({
            'space_name': '--- TOTAL BUILDING ---',
            'volume_m3': sum(s['volume_m3'] for s in building_model.spaces_data if 'volume_m3' in s),
            'total_boundary_surface_m2': sum(s['total_boundary_surface_m2'] for s in building_model.spaces_data if 'total_boundary_surface_m2' in s),
            'Q_trans_W': total_building_q_trans,
            'Q_vent_W': total_building_q_vent,
            'Q_total_W': total_building_q_total,
            'delta_t_K': delta_t
        })
    return results

# --- Weather Data Function (single day, kept for potential other uses) ---
def get_outdoor_temperature(lat, lon, alt, target_date):
    """
    Fetches the average outdoor temperature for a given location and date.
    """
    location = Point(lat, lon, alt)
    try:
        data = Daily(location, start=target_date, end=target_date)
        df = data.fetch()
        if not df.empty and 'tavg' in df.columns and not np.isnan(df['tavg'].iloc[0]):
            tavg = df['tavg'].iloc[0]
            # print(f"Successfully fetched outdoor temperature for {target_date.strftime('%Y-%m-%d')}: {tavg:.1f} °C") # Less verbose
            return float(tavg)
        else:
            print(f"Warning: Could not fetch valid outdoor temperature for {target_date.strftime('%Y-%m-%d')}. Using default 0°C for this day.")
            return 0.0 
    except Exception as e:
        print(f"Error fetching weather data for {target_date.strftime('%Y-%m-%d')}: {e}. Using default 0°C for this day.")
        return 0.0

# --- Yearly Energy Consumption Plot Function ---
def plot_yearly_energy_consumption(building_model: BuildingModel, lat: float, lon: float, alt: int, year: int, indoor_temp_celsius: float, ach: float):
    """
    Fetches daily outdoor temperatures, calculates daily heat loss for the building,
    converts it to energy consumption (kWh), and plots it for a given year.
    Saves the plot to a file.

    Args:
        building_model (BuildingModel): The preprocessed building model.
        lat (float): Latitude.
        lon (float): Longitude.
        alt (int): Altitude in meters.
        year (int): The year for which to calculate and plot energy consumption.
        indoor_temp_celsius (float): The constant indoor temperature.
        ach (float): Air changes per hour for ventilation loss.
    """
    if not building_model or not building_model.spaces_data:
        print("Cannot plot yearly energy consumption: Building model is not loaded or has no space data.")
        return

    print(f"\nFetching daily temperatures and calculating energy consumption for {year} at Lat: {lat}, Lon: {lon}...")
    start_date = datetime(year, 1, 1)
    end_date = datetime(year, 12, 31)
    location = Point(lat, lon, alt)

    weather_data = Daily(location, start=start_date, end=end_date)
    df_weather = weather_data.fetch()

    if df_weather.empty or 'tavg' not in df_weather.columns:
        print(f"No temperature data found for {year}.")
        return

    df_cleaned = df_weather.dropna(subset=['tavg'])
    if df_cleaned.empty:
        print(f"All temperature data for {year} was invalid (NaN). Cannot plot energy consumption.")
        return

    dates = []
    daily_energy_kwh_list = []
    
    print(f"Calculating daily energy consumption for {year} (Indoor Temp: {indoor_temp_celsius}°C, ACH: {ach})...")
    for date_index, row in df_cleaned.iterrows():
        outdoor_temp_celsius_today = row['tavg']
        
        heat_loss_results_today = calculate_total_heat_loss(building_model, indoor_temp_celsius, outdoor_temp_celsius_today, ach)
        
        total_building_heat_loss_w = 0
        if heat_loss_results_today:
            for res in heat_loss_results_today:
                if res['space_name'] == '--- TOTAL BUILDING ---':
                    total_building_heat_loss_w = res['Q_total_W']
                    break
        
        # Convert average power (Watts) over a day to energy (kWh)
        # Energy (kWh) = Power (W) * hours / 1000
        # Since Q_total_W is already an average power, we multiply by 24 hours
        energy_kwh_today = (total_building_heat_loss_w * 24) / 1000.0
        
        # Ensure non-negative energy consumption (heating demand)
        if energy_kwh_today < 0:
            energy_kwh_today = 0

        dates.append(date_index) # Use the DataFrame index which is a DatetimeIndex
        daily_energy_kwh_list.append(energy_kwh_today)

    if not daily_energy_kwh_list:
        print("No energy consumption data could be calculated to plot.")
        return

    fig = plt.figure(figsize=(12, 6))
    plt.plot(dates, daily_energy_kwh_list, label=f'Szacowane dzienne zużycie energii na ogrzewanie (kWh) w {year}', color='orangered')
    
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b')) 
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator()) 
    plt.gcf().autofmt_xdate() 

    plt.title(f'Szacowane dzienne zużycie energii na ogrzewanie dla lokalizacji ({lat:.2f}, {lon:.2f}) w roku {year}\n(Temp. wewn.: {indoor_temp_celsius}°C, ACH: {ach})')
    plt.xlabel('Miesiąc')
    plt.ylabel('Dzienne zużycie energii na ogrzewanie (kWh)')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    
    plot_filename = f"yearly_energy_consumption_{year}.png"
    plt.savefig(plot_filename)
    print(f"Wykres zużycia energii pomyślnie zapisany jako {plot_filename}")


# --- Main Execution ---
if __name__ == '__main__':
    IFC_FILE_PATH = 'model.ifc' # Upewnij się, że ten plik istnieje lub podaj poprawną ścieżkę
    
    LATITUDE = 52.2297   # Warszawa
    LONGITUDE = 21.0122  # Warszawa
    ALTITUDE = 110       # Przybliżona wysokość dla Warszawy (w metrach)
    YEAR_FOR_PLOT = 2023 # Rok do analizy
    
    DEFAULT_ACH = 0.5             # Wymiany powietrza na godzinę (Air Changes per Hour)
    INDOOR_TEMPERATURE_CELSIUS = 20.0 # Stała temperatura wewnętrzna

    print(f"Ładowanie i przetwarzanie modelu budynku z: {IFC_FILE_PATH}")
    building = None # Initialize building to None
    try:
        building = BuildingModel(IFC_FILE_PATH)
        if not building.spaces_data:
            print("ERROR: Model IFC został załadowany, ale nie znaleziono lub nie przetworzono żadnych przestrzeni (IfcSpace) do analizy.")
            print("Sprawdź model IFC lub logi preprocesora.")
            building = None # Set to None if no spaces processed
    except FileNotFoundError:
        print(f"BŁĄD KRYTYCZNY: Plik IFC nie znaleziony w '{IFC_FILE_PATH}'. Sprawdź ścieżkę.")
        # building remains None
    except Exception as e:
        print(f"BŁĄD KRYTYCZNY podczas ładowania modelu IFC: {e}")
        # building remains None

    if building:
        print(f"\n--- Generowanie wykresu rocznego zużycia energii ---")
        plot_yearly_energy_consumption(
            building_model=building,
            lat=LATITUDE,
            lon=LONGITUDE,
            alt=ALTITUDE,
            year=YEAR_FOR_PLOT,
            indoor_temp_celsius=INDOOR_TEMPERATURE_CELSIUS,
            ach=DEFAULT_ACH
        )
    else:
        print("\nNie można wygenerować wykresu zużycia energii, ponieważ model budynku nie został pomyślnie załadowany lub przetworzony.")
    
    print("\n--- Skrypt zakończył działanie ---")
