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
                if hasattr(rel_sb, 'CalculatedArea') and rel_sb.CalculatedArea is not None:
                    area = float(rel_sb.CalculatedArea)
                
                if area <= 0:
                    area_psets = {
                        'PSet_WallCommon': 'GrossArea', 'PSet_WindowCommon': 'Area', 
                        'PSet_DoorCommon': 'Area', 'PSet_RoofCommon': 'GrossArea',
                        'PSet_SlabCommon': 'GrossArea'
                    }
                    for pset_name, prop_name in area_psets.items():
                        elem_area = self._get_pset_value(related_element, pset_name, prop_name)
                        if elem_area > 0:
                            area = elem_area
                            break 
                
                if area <= 0:
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
        # print(f"Indoor temperature ({indoor_temp}°C) is not higher than outdoor temperature ({outdoor_temp}°C). Heat loss is assumed to be zero.") # Original print
        # More informative print for this case
        if indoor_temp < outdoor_temp:
            print(f"Indoor temperature ({indoor_temp}°C) is lower than outdoor temperature ({outdoor_temp}°C). This indicates heat GAIN, not loss. Calculations will show negative loss (gain) or zero if capped.")
        elif indoor_temp == outdoor_temp:
            print(f"Indoor temperature ({indoor_temp}°C) is equal to outdoor temperature ({outdoor_temp}°C). No temperature difference, so heat loss is zero.")
        # Proceed with calculation; delta_t will be <=0, leading to zero or negative (gain) heat loss values.
        
        for space_data in building_model.spaces_data: # Still provide structure if delta_t is 0
            results.append({
                'space_name': space_data['space_name'],
                'volume_m3': space_data['volume_m3'],
                'total_boundary_surface_m2': space_data['total_boundary_surface_m2'],
                'Q_trans_W': 0.0, # Explicitly zero if no delta_t
                'Q_vent_W': 0.0,  # Explicitly zero if no delta_t
                'Q_total_W': 0.0, # Explicitly zero if no delta_t
                'delta_t_K': delta_t
            })
        # Add total building summary even if delta_t is zero
        if results:
             results.append({
                'space_name': '--- TOTAL BUILDING ---',
                'volume_m3': sum(s['volume_m3'] for s in building_model.spaces_data),
                'total_boundary_surface_m2': sum(s['total_boundary_surface_m2'] for s in building_model.spaces_data),
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
            'volume_m3': sum(s['volume_m3'] for s in building_model.spaces_data),
            'total_boundary_surface_m2': sum(s['total_boundary_surface_m2'] for s in building_model.spaces_data),
            'Q_trans_W': total_building_q_trans,
            'Q_vent_W': total_building_q_vent,
            'Q_total_W': total_building_q_total,
            'delta_t_K': delta_t
        })
    return results

# --- Weather Data Function ---
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
            print(f"Successfully fetched outdoor temperature for {target_date.strftime('%Y-%m-%d')}: {tavg:.1f} °C")
            return float(tavg)
        else:
            print(f"Warning: Could not fetch valid outdoor temperature for {target_date.strftime('%Y-%m-%d')}. Using default 0°C.")
            return 0.0 
    except Exception as e:
        print(f"Error fetching weather data for {target_date.strftime('%Y-%m-%d')}: {e}. Using default 0°C.")
        return 0.0

# --- Yearly Temperature Plot Function ---
def plot_yearly_temperatures(lat, lon, alt, year):
    """
    Fetches daily average temperatures for a given year and location,
    and plots them using matplotlib. Saves the plot to a file.
    
    Args:
        lat (float): Latitude.
        lon (float): Longitude.
        alt (int): Altitude in meters.
        year (int): The year for which to fetch and plot temperatures.
    """
    print(f"\nFetching daily temperatures for {year} at Lat: {lat}, Lon: {lon}...")
    start_date = datetime(year, 1, 1)
    end_date = datetime(year, 12, 31)
    location = Point(lat, lon, alt)
    
    # Get the current figure manager to check the backend, if needed for debugging
    # current_backend = plt.get_backend()
    # print(f"Matplotlib backend in use: {current_backend}")

    try:
        data = Daily(location, start=start_date, end=end_date)
        df = data.fetch()

        if df.empty or 'tavg' not in df.columns:
            print(f"No temperature data found for {year}.")
            return

        df_cleaned = df.dropna(subset=['tavg'])
        if df_cleaned.empty:
            print(f"All temperature data for {year} was invalid (NaN). Cannot plot.")
            return

        dates = df_cleaned.index
        temps = df_cleaned['tavg']

        fig = plt.figure(figsize=(12, 6)) # Get a reference to the figure
        plt.plot(dates, temps, label=f'Średnia temperatura dzienna (°C) w {year}')
        
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b')) 
        plt.gca().xaxis.set_major_locator(mdates.MonthLocator()) 
        plt.gcf().autofmt_xdate() 

        plt.title(f'Dzienne średnie temperatury dla lokalizacji ({lat:.2f}, {lon:.2f}) w roku {year}')
        plt.xlabel('Data')
        plt.ylabel('Średnia temperatura (°C)')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        
        # Save the plot to a file
        plot_filename = f"yearly_temperatures_{year}.png"
        plt.savefig(plot_filename)
        print(f"Plot successfully saved as {plot_filename}")
        
        # Attempt to show the plot. This might still show a warning if the backend is non-interactive,
        # but the plot has been saved.
        print(f"Attempting to display plot for {year}...")
        plt.show()
        
        # Close the figure to free up memory
        plt.close(fig) # Use the figure reference

    except Exception as e:
        print(f"Error generating temperature plot for {year}: {e}")
        # Ensure the figure is closed even if an error occurs mid-plotting
        if 'fig' in locals() and plt.fignum_exists(fig.number):
            plt.close(fig)


# --- Main Execution ---
if __name__ == '__main__':
    IFC_FILE_PATH = 'model.ifc' 
    
    LATITUDE = 52.2297  
    LONGITUDE = 21.0122 
    ALTITUDE = 110     
    YEAR_FOR_PLOT = 2024 
    DEFAULT_ACH = 0.5  

    print(f"Loading and preprocessing building model from: {IFC_FILE_PATH}")
    try:
        building = BuildingModel(IFC_FILE_PATH)
    except FileNotFoundError:
        print(f"ERROR: IFC file not found at '{IFC_FILE_PATH}'. Please check the path.")
        building = None
    except Exception as e:
        print(f"Critical error loading IFC model: {e}")
        building = None 

    print(f"\n--- Generating Yearly Temperature Plot ---")
    plot_yearly_temperatures(LATITUDE, LONGITUDE, ALTITUDE, YEAR_FOR_PLOT)
    
    print("\n--- Script Finished ---")

