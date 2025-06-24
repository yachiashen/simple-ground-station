import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt 
from cartopy.geodesic import Geodesic
from skyfield.api import wgs84, load
import numpy as np
from PIL import Image
import io
from datetime import timedelta
import datetime
from shapely.geometry import Polygon 
import imageio
import os
import folium 
import logging

logger = logging.getLogger(__name__)

# Skyfield timescale
ts = load.timescale()
# Mean Earth radius (km) - used for footprint calculation
R_EARTH_KM_MEAN = 6371.0
# Define available map tile options for Folium
ALL_TILE_OPTIONS = {
    "OpenStreetMap": {
        "tiles_param": "OpenStreetMap",
        "name": "OpenStreetMap",
        "attr": '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
    },
    "Stamen Terrain": {
        "tiles_param": "https://tiles.stadiamaps.com/tiles/stamen_terrain/{z}/{x}/{y}.png",
        "name": "Stamen Terrain (via Stadia Maps)",
        "attr": 'Map tiles by <a href="http://stamen.com">Stamen Design</a> (<a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a>). Data by <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors. Tiles hosted by <a href="https://stadiamaps.com/">Stadia Maps</a>.'
    },
    "Stamen Toner": {
        "tiles_param": "https://tiles.stadiamaps.com/tiles/stamen_toner/{z}/{x}/{y}.png",
        "name": "Stamen Toner (via Stadia Maps)",
        "attr": 'Map tiles by <a href="http://stamen.com">Stamen Design</a> (<a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a>). Data by <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors. Tiles hosted by <a href="https://stadiamaps.com/">Stadia Maps</a>.'
    },
    "CartoDB positron": {
        "tiles_param": "cartodbpositron",
        "name": "CartoDB Positron",
        "attr": '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors &copy; <a href="https://carto.com/attributions">CARTO</a>'
    },
    "Esri World Imagery": { 
        "tiles_param": 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        "name": "Esri World Imagery",
        "attr": 'Tiles &copy; Esri &mdash; Source: Esri, Maxar, Earthstar Geographics, and the GIS User Community'
    }
}

AVAILABLE_MAP_STYLES = list(ALL_TILE_OPTIONS.keys())
DEFAULT_MAP_STYLE = "OpenStreetMap"


# --- Folium Interactive Map Function ---
def generate_interactive_map_folium(satellite, 
                                    observer_location, 
                                    track_reference_time_utc: datetime.datetime, 
                                    display_time_utc: datetime.datetime, 
                                    track_duration_hours: float = 3, 
                                    points_per_hour: int = 60,
                                    min_elevation_for_footprint_deg: float = 5.0,
                                    map_tiles: str = 'OpenStreetMap',
                                    initial_zoom: int = 2,
                                    cached_track_data: dict | None = None,
                                    perform_fit_bounds: bool = True,
                                    center_location: tuple[float, float] | None = None,
                                    current_zoom_level: int | None = None
                                    ):
    """
    Generate an interactive Folium map showing:
      - The satellite ground track from `track_reference_time_utc` for `track_duration_hours`
      - The observer location
      - The satellite position and footprint at `display_time_utc`

    It can optionally reuse `cached_track_data` to avoid recomputing the track.
    Returns a tuple (map_html, newly_computed_track_data_or_none).

    Parameters
    ----------
    satellite : EarthSatellite
        Skyfield satellite object.
    observer_location : Topos
        Skyfield observer location object.
    track_reference_time_utc : datetime.datetime or Skyfield Time
        Start/reference time for the ground track.
    display_time_utc : datetime.datetime or Skyfield Time
        Time at which to render the satellite marker and footprint.
    track_duration_hours : float, default 3
        Duration (hours) for the ground track.
    points_per_hour : int, default 60
        Number of samples per hour along the track.
    min_elevation_for_footprint_deg : float, default 5.0
        Minimum elevation used to compute the footprint radius.
    map_tiles : str, default 'OpenStreetMap'
        Base tile key (must exist in ALL_TILE_OPTIONS).
    initial_zoom : int, default 2
        Initial zoom level if no fit_bounds is applied.
    cached_track_data : dict | None
        Optional cache containing "coords", "times_iso", "altitudes_km", etc.
    perform_fit_bounds : bool, default True
        If True, fit map bounds to include track, satellite and observer.
    center_location : tuple[float, float] | None
        Explicit map center (lat, lon) if provided.
    current_zoom_level : int | None
        Explicit zoom level if provided.

    Returns
    -------
    tuple[str, dict | None]
        HTML representation of the Folium map and optionally newly computed
        track data (or None if cache was used).
    """
    try:
        newly_computed_track_data = None

        # Try to reuse cached track data if provided and complete
        if cached_track_data:
            logger.info("Using cached track data for Folium map.")
            track_coordinates = cached_track_data.get("coords")
            track_times_iso = cached_track_data.get("times_iso")
            track_altitudes_km = cached_track_data.get("altitudes_km")
            if not all([track_coordinates, track_times_iso, track_altitudes_km]):
                logger.warning("Cached track data is incomplete. Recomputing.")
                cached_track_data = None

        # Compute track if no usable cache
        if not cached_track_data:
            logger.info("No valid cached track data. Computing new track for Folium map.")
            
            track_ref_dt = track_reference_time_utc.utc_datetime()
            t0_track_skyfield = ts.utc(track_ref_dt.year, track_ref_dt.month, track_ref_dt.day,
                                       track_ref_dt.hour, track_ref_dt.minute, track_ref_dt.second)
            
            num_track_points = int(track_duration_hours * points_per_hour)
            if num_track_points <= 0: num_track_points = 1

            current_track_end_skyfield_time_obj = track_reference_time_utc + datetime.timedelta(hours=track_duration_hours)
            current_track_end_dt = current_track_end_skyfield_time_obj.utc_datetime()
            track_end_time_skyfield = ts.utc(current_track_end_dt.year, current_track_end_dt.month, current_track_end_dt.day,
                                             current_track_end_dt.hour, current_track_end_dt.minute, current_track_end_dt.second)

            track_times_skyfield_objects = ts.linspace(t0_track_skyfield, track_end_time_skyfield, num_track_points)
            geocentric_track = satellite.at(track_times_skyfield_objects)
            subpoints_track = wgs84.subpoint(geocentric_track)
            track_coordinates = list(zip(subpoints_track.latitude.degrees, subpoints_track.longitude.degrees))
            track_altitudes_km = list(subpoints_track.elevation.km)
            track_times_iso = [t.utc_iso() for t in track_times_skyfield_objects]

            newly_computed_track_data = {
                "coords": track_coordinates,
                "times_iso": track_times_iso,
                "altitudes_km": track_altitudes_km,
                "track_reference_time_utc_iso": track_reference_time_utc.utc_iso(),
                "track_duration_hours": track_duration_hours,
                "points_per_hour": points_per_hour
            }
        
        actual_track_end_time_utc = track_reference_time_utc + datetime.timedelta(hours=track_duration_hours)
        
        # Satellite position at display_time_utc (always computed)
        logger.info(f"generate_interactive_map_folium: Received display_time_utc: {display_time_utc.utc_iso() if display_time_utc is not None else 'None'}")
        display_time_dt = display_time_utc.utc_datetime()
        t_display_skyfield = ts.utc(display_time_dt.year, display_time_dt.month, display_time_dt.day,
                                    display_time_dt.hour, display_time_dt.minute,
                                    display_time_dt.second + display_time_dt.microsecond/1e6)
        logger.info(f"generate_interactive_map_folium: Converted to t_display_skyfield: {t_display_skyfield.utc_iso()}")
        
        geocentric_at_display_time = satellite.at(t_display_skyfield)
        subpoint_at_display_time = wgs84.subpoint(geocentric_at_display_time)
        current_sat_lat = subpoint_at_display_time.latitude.degrees
        current_sat_lon = subpoint_at_display_time.longitude.degrees
        current_sat_alt_km = subpoint_at_display_time.elevation.km
        logger.info(f"generate_interactive_map_folium: Calculated satellite position at {t_display_skyfield.utc_iso()}: Lat={current_sat_lat:.4f}, Lon={current_sat_lon:.4f}, Alt={current_sat_alt_km:.2f} km")

        # Observer lat/lon
        obs_lat = observer_location.latitude.degrees
        obs_lon = observer_location.longitude.degrees

        # Map center and zoom
        map_loc = None
        zoom_lvl = initial_zoom
        if center_location:
            map_loc = center_location
        else:
            track_coord_defined = track_coordinates and len(track_coordinates) > 0 and track_coordinates[0] is not None
            map_loc = [
                current_sat_lat if current_sat_lat is not None else (track_coordinates[0][0] if track_coord_defined else (obs_lat if obs_lat is not None else 0)),
                current_sat_lon if current_sat_lon is not None else (track_coordinates[0][1] if track_coord_defined else (obs_lon if obs_lon is not None else 0))
            ]
        if current_zoom_level is not None:
            zoom_lvl = current_zoom_level
        
        m = folium.Map(
            location=map_loc, 
            zoom_start=zoom_lvl, 
            tiles=None,
            max_bounds=[[-85.05112878, -180.0], [85.05112878, 180.0]],
            world_copy_jump=False
        )

        # Base layers
        if map_tiles not in ALL_TILE_OPTIONS:
            logger.warning(f"Provided map_tiles key '{map_tiles}' not found in ALL_TILE_OPTIONS. Defaulting to OpenStreetMap.")
            map_tiles = DEFAULT_MAP_STYLE

        for layer_key, config in ALL_TILE_OPTIONS.items():
            is_default_layer = (layer_key == map_tiles)
            folium.TileLayer(
                tiles=config["tiles_param"], 
                name=config["name"],
                attr=config["attr"],
                overlay=False,
                control=True,
                show=is_default_layer,
                no_wrap=True
            ).add_to(m)

        # Ground track with anti-dateline splitting
        if track_coordinates:
            segments = []
            current_segment = []
            if len(track_coordinates) > 0:
                current_segment.append(track_coordinates[0])
                for i in range(1, len(track_coordinates)):
                    prev_point = track_coordinates[i-1]
                    curr_point = track_coordinates[i]
                    prev_lon = prev_point[1]
                    curr_lon = curr_point[1]
                    if abs(curr_lon - prev_lon) > 180:
                        segments.append(list(current_segment))
                        current_segment = [curr_point]
                    else:
                        current_segment.append(curr_point)
                if current_segment:
                    segments.append(list(current_segment))
            
            for i, segment_coords in enumerate(segments):
                if len(segment_coords) > 1:
                    folium.PolyLine(
                        locations=segment_coords, color='red', weight=2.5, opacity=1,
                        tooltip=f'{satellite.name} Ground Track (Segment {i+1})'
                    ).add_to(m)

            # Start/end markers
            if len(track_coordinates) > 0:
                folium.Marker(
                    location=track_coordinates[0],
                    popup=f'Track Start: {track_reference_time_utc.utc_datetime().strftime("%Y-%m-%d %H:%M:%S")} UTC',
                    icon=folium.Icon(color='green', icon='play')
                ).add_to(m)
            if len(track_coordinates) > 1:
                folium.Marker(
                    location=track_coordinates[-1],
                    popup=f'Track End: {actual_track_end_time_utc.utc_datetime().strftime("%Y-%m-%d %H:%M:%S")} UTC',
                    icon=folium.Icon(color='orange', icon='stop')
                ).add_to(m)
            
            # Sparse markers along the track for context
            if track_times_iso and track_altitudes_km and len(track_times_iso) == len(track_coordinates) and len(track_altitudes_km) == len(track_coordinates):
                num_display_track_points = len(track_coordinates)
                marker_frequency = max(1, num_display_track_points // 20 if num_display_track_points >= 20 else 1)
                for i, coord in enumerate(track_coordinates):
                    if i % marker_frequency == 0:
                        try:
                            point_time_utc_obj = datetime.datetime.fromisoformat(track_times_iso[i].replace('Z', '+00:00'))
                            altitude_km_pt = track_altitudes_km[i]
                            tooltip_text = (
                                f"<b>Satellite: {satellite.name}</b><br>"
                                f"Time: {point_time_utc_obj.strftime('%Y-%m-%d %H:%M:%S')} UTC<br>"
                                f"Lat: {coord[0]:.2f}°, Lon: {coord[1]:.2f}°<br>"
                                f"Alt: {altitude_km_pt:.2f} km"
                            )
                            folium.CircleMarker(
                                location=coord, radius=3, color='darkred', weight=1,
                                fill=True, fill_color='red', fill_opacity=0.7, tooltip=tooltip_text
                            ).add_to(m)
                        except Exception as e_marker:
                            logger.error(f"Error creating track point marker at index {i}: {e_marker}")
            else:
                logger.warning("Mismatch in track data lengths (times_iso, altitudes_km, coords) or data missing. Skipping intermittent track markers.")


        # Observer marker
        if obs_lat is not None and obs_lon is not None:
            folium.Marker(
                location=[obs_lat, obs_lon],
                popup=f'Observer<br>Lat: {obs_lat:.2f}°, Lon: {obs_lon:.2f}°',
                tooltip='Observer Location', icon=folium.Icon(color='blue', icon='user')
            ).add_to(m)

        # Satellite marker at selected time and its footprint
        if current_sat_lat is not None and current_sat_lon is not None:
            popup_html_current_sat = (
                f"<b>{satellite.name} (Selected Time)</b><br>"
                f"Time: {display_time_utc.utc_datetime().strftime('%Y-%m-%d %H:%M:%S')} UTC<br>"
                f"Lat: {current_sat_lat:.2f}°, Lon: {current_sat_lon:.2f}°<br>"
                f"Altitude: {current_sat_alt_km:.2f} km"
            )
            folium.Marker(
                location=[current_sat_lat, current_sat_lon],
                popup=folium.Popup(popup_html_current_sat, max_width=300),
                tooltip=f'{satellite.name} (at selected time)',
                icon=folium.Icon(color='purple', icon='satellite', prefix='fa')
            ).add_to(m)

            # Footprint (simple great-circle radius approximation)
            if current_sat_alt_km > 0:
                try:
                    el_min_rad = np.deg2rad(min_elevation_for_footprint_deg)
                    cos_gamma_arg = (R_EARTH_KM_MEAN / (R_EARTH_KM_MEAN + current_sat_alt_km)) * np.cos(el_min_rad)
                    cos_gamma_arg = np.clip(cos_gamma_arg, -1.0, 1.0)
                    gamma = np.arccos(cos_gamma_arg) - el_min_rad
                    if gamma > 0:
                        footprint_radius_meters = R_EARTH_KM_MEAN * 1000 * gamma
                        folium.Circle(
                            location=[current_sat_lat, current_sat_lon], radius=footprint_radius_meters,
                            color='magenta', fill=True, fill_color='magenta', fill_opacity=0.15,
                            tooltip=(f'Footprint (El ≥ {min_elevation_for_footprint_deg}°)<br>'
                                     f'Radius: {footprint_radius_meters/1000:.0f} km<br>'
                                     f'Time: {display_time_utc.utc_datetime().strftime("%H:%M:%S")} UTC')
                        ).add_to(m)
                except Exception as e_fp:
                    logger.error(f"Error calculating or drawing Folium footprint at display_time_utc: {e_fp}", exc_info=True)
        
        # Fit bounds if requested and enough points exist
        if perform_fit_bounds and track_coordinates and len(track_coordinates) > 0:
            all_lats = [p[0] for p in track_coordinates]
            all_lons = [p[1] for p in track_coordinates]
            if current_sat_lat is not None: all_lats.append(current_sat_lat)
            if current_sat_lon is not None: all_lons.append(current_sat_lon)
            if obs_lat is not None: all_lats.append(obs_lat)
            if obs_lon is not None: all_lons.append(obs_lon)
            
            all_lats = [lat for lat in all_lats if lat is not None]
            all_lons = [lon for lon in all_lons if lon is not None]

            if all_lats and all_lons: 
                min_lat, max_lat = min(all_lats), max(all_lats)
                min_lon, max_lon = min(all_lons), max(all_lons)
                padding_lat = (max_lat - min_lat) * 0.1 if (max_lat - min_lat) > 0 else 0.1
                padding_lon = (max_lon - min_lon) * 0.1 if (max_lon - min_lon) > 0 else 0.1
                bounds = [[min_lat - padding_lat, min_lon - padding_lon], [max_lat + padding_lat, max_lon + padding_lon]]
                bounds[0][0] = max(bounds[0][0], -85)
                bounds[0][1] = max(bounds[0][1], -180)
                bounds[1][0] = min(bounds[1][0], 85)
                bounds[1][1] = min(bounds[1][1], 180)
                
                if bounds[0][0] < bounds[1][0] and bounds[0][1] < bounds[1][1]:
                    logger.info(f"Performing fit_bounds with bounds: {bounds}")
                    m.fit_bounds(bounds)
                else:
                    logger.warning(f"Skipping fit_bounds due to invalid calculated bounds: {bounds}. Map will use explicit center/zoom or defaults.")
            else:
                logger.info("Not enough points to calculate bounds for fit_bounds. Map will use explicit center/zoom or defaults.")
        elif not perform_fit_bounds:
            logger.info("perform_fit_bounds is False. Skipping m.fit_bounds(). Map will use explicit center/zoom.")
        else: 
            logger.info("perform_fit_bounds is True, but no track_coordinates available or other issue. Skipping m.fit_bounds(). Map will use explicit center/zoom or defaults.")

        if len(ALL_TILE_OPTIONS) > 1:
            folium.LayerControl().add_to(m)

        map_html = m._repr_html_()
        return map_html, newly_computed_track_data
    except Exception as e:
        logger.error(f"Error generating Folium interactive map: {e}", exc_info=True)
        error_map_html = f"<html><body><p>Error generating map: {e}</p></body></html>"
        return error_map_html, None

# Skyfield timescale
ts = load.timescale()
# Mean Earth radius (km) - used for footprint calculation
R_EARTH_KM_MEAN = 6371.0

def _plot_single_frame_for_gif(satellite, observer_location, 
                               current_time_utc: datetime.datetime, 
                               base_track_start_time_utc: datetime.datetime,
                               base_track_duration_hours=3, 
                               points_per_hour=60,
                               min_elevation_for_footprint_deg=0.0,
                               map_style='Default',
                               fig_size=(10, 8)):
    """
    Render a single matplotlib/cartopy frame for the GIF animation.

    The frame includes:
      - Full ground track for a base window starting at `base_track_start_time_utc`
      - The observer location
      - The footprint circle and satellite marker at `current_time_utc`

    Parameters
    ----------
    satellite : EarthSatellite
        Skyfield satellite object.
    observer_location : Topos
        Skyfield observer location.
    current_time_utc : datetime.datetime
        Simulation time for the current frame (UTC).
    base_track_start_time_utc : datetime.datetime
        Start time for the background ground track.
    base_track_duration_hours : int, default 3
        Hours of the ground track to show in the background.
    points_per_hour : int, default 60
        Number of track samples per hour.
    min_elevation_for_footprint_deg : float, default 0.0
        Minimum elevation for footprint calculation.
    map_style : str, default 'Default'
        One of ('Default', 'Stamen Terrain', 'OpenStreetMap').
    fig_size : tuple, default (10, 8)
        Matplotlib figure size.

    Returns
    -------
    PIL.Image.Image
        The rendered frame as a PIL Image. On error, returns a blank
        white image sized according to fig_size.
    """
    try:
        # 1) Build time range for full ground track
        t0_track = ts.utc(base_track_start_time_utc.year, base_track_start_time_utc.month, base_track_start_time_utc.day,
                          base_track_start_time_utc.hour, base_track_start_time_utc.minute, base_track_start_time_utc.second)
        end_time_track_utc = base_track_start_time_utc + datetime.timedelta(hours=base_track_duration_hours)
        t1_track = ts.utc(end_time_track_utc.year, end_time_track_utc.month, end_time_track_utc.day,
                          end_time_track_utc.hour, end_time_track_utc.minute, end_time_track_utc.second)

        num_points_track = int(base_track_duration_hours * points_per_hour)
        times_track = ts.linspace(t0_track, t1_track, num_points_track)

        # 2) Compute ground track subpoints
        geocentric_positions_track = satellite.at(times_track)
        subpoints_track = wgs84.subpoint(geocentric_positions_track)
        subpoint_lat_track = subpoints_track.latitude.degrees
        subpoint_lon_track = subpoints_track.longitude.degrees

        # 3) Initialize map
        fig = Figure(figsize=fig_size) # Use passed figsize
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.set_global()

        if map_style == 'Stamen Terrain':
            stamen_terrain = cimgt.Stamen('terrain-background')
            ax.add_image(stamen_terrain, 8)
        elif map_style == 'OpenStreetMap':
            osm_tiles = cimgt.OSM()
            ax.add_image(osm_tiles, 8)
        else:  # Default style
            ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black', facecolor='#BEB88A')
            ax.add_feature(cfeature.OCEAN, zorder=0, facecolor='#D3E5E5')
        
        if map_style not in ['Stamen Terrain', 'OpenStreetMap'] or map_style == 'Default':
            ax.add_feature(cfeature.COASTLINE, zorder=1)
            ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=1)
        else:
            ax.add_feature(cfeature.COASTLINE, zorder=3, edgecolor='black', linewidth=0.5)
            ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=3, edgecolor='gray')
        ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, zorder=4)

        # 4) Plot ground track
        ax.plot(subpoint_lon_track, subpoint_lat_track, 'r-', transform=ccrs.Geodetic(), label=f'{satellite.name} Ground Track', zorder=5)

        # 5) Plot observer location
        obs_lat = observer_location.latitude.degrees
        obs_lon = observer_location.longitude.degrees
        ax.plot(obs_lon, obs_lat, 'bo', markersize=7, transform=ccrs.Geodetic(), label='Observer', zorder=6)

        # 6) Footprint at current_time_utc
        current_time_skyfield = ts.utc(current_time_utc.year, current_time_utc.month, current_time_utc.day,
                                       current_time_utc.hour, current_time_utc.minute, 
                                       current_time_utc.second + current_time_utc.microsecond/1e6)
        geocentric_current = satellite.at(current_time_skyfield)
        subpoint_current = wgs84.subpoint(geocentric_current)
        current_lat_deg = subpoint_current.latitude.degrees
        current_lon_deg = subpoint_current.longitude.degrees
        current_alt_km = subpoint_current.elevation.km

        if current_alt_km > 0:
            try:
                if min_elevation_for_footprint_deg > 0.0:
                    el_min_rad = np.deg2rad(min_elevation_for_footprint_deg)
                    cos_gamma_term = (R_EARTH_KM_MEAN / (R_EARTH_KM_MEAN + current_alt_km)) * np.cos(el_min_rad)
                    if cos_gamma_term > 1.0: cos_gamma_term = 1.0 
                    elif cos_gamma_term < -1.0: cos_gamma_term = -1.0
                    gamma = np.arccos(cos_gamma_term) - el_min_rad
                else:
                    gamma = np.arccos(R_EARTH_KM_MEAN / (R_EARTH_KM_MEAN + current_alt_km))
                
                if gamma > 0:
                    footprint_radius_meters = R_EARTH_KM_MEAN * 1000 * gamma
                    circle_points = Geodesic().circle(lon=current_lon_deg, lat=current_lat_deg, 
                                                    radius=footprint_radius_meters, n_samples=180)
                    footprint_polygon = Polygon(circle_points)
                    ax.add_geometries([footprint_polygon], ccrs.Geodetic(),
                                      facecolor='blue', alpha=0.25, edgecolor='blue', 
                                      linewidth=0.8, zorder=2)
                    ax.plot(current_lon_deg, current_lat_deg, 
                            marker='o', color='cyan', markersize=5, markeredgecolor='black',
                            transform=ccrs.Geodetic(), 
                            label=f'Footprint @ {current_time_utc.strftime("%H:%M:%S")}', 
                            zorder=6)
            except Exception:
                pass
        
        plt.title(f"Satellite: {satellite.name} - Footprint @ {current_time_utc.strftime('%Y-%m-%d %H:%M:%S')} UTC")

        img_buffer = io.BytesIO()
        canvas.print_png(img_buffer)
        img_buffer.seek(0)
        pil_image = Image.open(img_buffer)
        plt.close(fig)
        return pil_image

    except Exception as e:
        print(f"Error plotting single GIF frame for time {current_time_utc}: {e}")
        return Image.new('RGB', (int(fig_size[0]*100), int(fig_size[1]*100)), color = 'white')


def generate_footprint_animation_gif(
    satellite, 
    observer_location, 
    animation_start_time_utc: datetime.datetime, 
    animation_duration_minutes: int = 60, 
    simulation_time_step_seconds: int = 30,
    gif_playback_fps: int = 2,
    output_gif_path: str = "data/footprint_animation.gif",
    base_track_duration_hours: int = 3, 
    min_elevation_for_footprint_deg: float = 5.0,
    map_style: str = 'Default',
    fig_size_inches: tuple = (8, 6) 
):
    """
    Generate a GIF animation of the satellite footprint over time.

    Each frame shows the full background ground track (for a base window)
    plus the footprint and satellite marker at the frame's time.

    Parameters
    ----------
    satellite : EarthSatellite
        Skyfield satellite object.
    observer_location : Topos
        Skyfield observer location object.
    animation_start_time_utc : datetime.datetime
        Start time of the animation (UTC).
    animation_duration_minutes : int, default 60
        Total duration of the animation (minutes).
    simulation_time_step_seconds : int, default 30
        Simulation seconds per frame (controls frame count).
    gif_playback_fps : int, default 2
        Output GIF playback frames per second.
    output_gif_path : str, default "data/footprint_animation.gif"
        Output path for the generated GIF.
    base_track_duration_hours : int, default 3
        Hours of ground track to show in background per frame.
    min_elevation_for_footprint_deg : float, default 5.0
        Minimum elevation for footprint calculation.
    map_style : str, default 'Default'
        One of ('Default', 'Stamen Terrain', 'OpenStreetMap').
    fig_size_inches : tuple, default (8, 6)
        Figure size in inches for each frame.

    Returns
    -------
    str | None
        Path to the generated GIF on success, otherwise None.
    """
    frames = []
    
    if simulation_time_step_seconds <= 0:
        print("Error: simulation_time_step_seconds must be positive.")
        return None
    
    total_animation_seconds = animation_duration_minutes * 60
    total_frames = int(np.ceil(total_animation_seconds / simulation_time_step_seconds)) if total_animation_seconds > 0 else 0

    if total_frames <= 0:
        print(f"Error: With duration {animation_duration_minutes} mins and step {simulation_time_step_seconds}s, no frames would be generated.")
        return None
    MAX_FRAMES = 2000 
    if total_frames > MAX_FRAMES:
        print(f"Error: Number of frames ({total_frames}) exceeds the maximum limit of {MAX_FRAMES}. Please increase simulation_time_step_seconds or reduce animation_duration_minutes.")
        return None
    elif total_frames > 1000:
        print(f"Warning: Number of frames ({total_frames}) is large. Generation may take some time.")

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_gif_path), exist_ok=True)
    
    print(f"Generating GIF: {total_frames} frames, each frame represents {simulation_time_step_seconds}s of simulation.")
    print(f"Animation duration: {animation_duration_minutes} mins. GIF playback FPS: {gif_playback_fps}.")

    current_animation_time_utc = animation_start_time_utc
    for i in range(total_frames): 
        if hasattr(current_animation_time_utc, 'utc_datetime'):
            time_str = current_animation_time_utc.utc_datetime().strftime('%Y-%m-%d %H:%M:%S')
        else:
            time_str = current_animation_time_utc.strftime('%Y-%m-%d %H:%M:%S')
        print(f"  Generating frame {i+1}/{total_frames} for time {time_str} UTC...")
        
        frame_pil = _plot_single_frame_for_gif(
            satellite=satellite,
            observer_location=observer_location,
            current_time_utc=current_animation_time_utc, 
            base_track_start_time_utc=animation_start_time_utc, 
            base_track_duration_hours=base_track_duration_hours,
            points_per_hour=30, 
            min_elevation_for_footprint_deg=min_elevation_for_footprint_deg,
            map_style=map_style,
            fig_size=fig_size_inches
        )
        if frame_pil:
            frames.append(frame_pil)
        
        current_animation_time_utc += datetime.timedelta(seconds=simulation_time_step_seconds)

    if frames:
        try:
            print(f"Saving GIF to {output_gif_path} ({len(frames)} frames) with playback FPS: {gif_playback_fps}...")
            imageio.mimsave(output_gif_path, frames, fps=gif_playback_fps, loop=0) 
            print(f"GIF saved successfully to {output_gif_path}")
            return output_gif_path
        except Exception as e:
            print(f"Error saving GIF: {e}")
            import traceback
            traceback.print_exc()
            return None
    else:
        print("No frames generated for GIF.")
        return None

def plot_satellite_track(satellite, observer_location, 
                         start_time_utc: datetime.datetime, 
                         duration_hours=3, 
                         points_per_hour=60,
                         min_elevation_for_footprint_deg=0.0,
                         map_style='Default',
                         footprint_times_utc: list[datetime.datetime] | None = None):
    """
    Plot the satellite ground track for a given time window and mark the observer.

    Optionally draw one or more footprint circles at times in `footprint_times_utc`.
    If `footprint_times_utc` is None, a single footprint is drawn at `start_time_utc`.

    Parameters
    ----------
    satellite : EarthSatellite
        Skyfield satellite object.
    observer_location : Topos
        Skyfield observer location object.
    start_time_utc : datetime.datetime
        Start time (UTC) of the plotted track.
    duration_hours : float, default 3
        Duration of the track (hours).
    points_per_hour : int, default 60
        Number of points per hour to sample along the track.
    min_elevation_for_footprint_deg : float, default 0.0
        Minimum elevation for footprint calculation.
    map_style : str, default 'Default'
        One of ('Default', 'Stamen Terrain', 'OpenStreetMap').
    footprint_times_utc : list[datetime.datetime] | None, default None
        Times at which to draw footprints.

    Returns
    -------
    tuple[PIL.Image.Image | None, str]
        (image, footprint_status_message_html). Image is None on error.
    """
    footprint_status_messages = []
    try:
        # 1) Build time range for the track
        t0 = ts.utc(start_time_utc.year, start_time_utc.month, start_time_utc.day,
                    start_time_utc.hour, start_time_utc.minute, start_time_utc.second)
        end_time_utc = start_time_utc + datetime.timedelta(hours=duration_hours)
        t1 = ts.utc(end_time_utc.year, end_time_utc.month, end_time_utc.day,
                    end_time_utc.hour, end_time_utc.minute, end_time_utc.second)

        num_points = int(duration_hours * points_per_hour)
        times = ts.linspace(t0, t1, num_points)

        # 2) Compute ground track subpoints
        geocentric_positions = satellite.at(times)
        subpoints = wgs84.subpoint(geocentric_positions)
        subpoint_lat = subpoints.latitude.degrees
        subpoint_lon = subpoints.longitude.degrees

        # 3) Initialize map
        fig = Figure(figsize=(10, 8))
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        
        ax.set_global()

        if map_style == 'Stamen Terrain':
            stamen_terrain = cimgt.Stamen('terrain-background')
            ax.add_image(stamen_terrain, 8)
        elif map_style == 'OpenStreetMap':
            osm_tiles = cimgt.OSM()
            ax.add_image(osm_tiles, 8)
        else:
            ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black', facecolor='#BEB88A')
            ax.add_feature(cfeature.OCEAN, zorder=0, facecolor='#D3E5E5')
        
        if map_style not in ['Stamen Terrain', 'OpenStreetMap'] or map_style == 'Default':
            ax.add_feature(cfeature.COASTLINE, zorder=1)
            ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=1)
        else:
            ax.add_feature(cfeature.COASTLINE, zorder=3, edgecolor='black', linewidth=0.5)
            ax.add_feature(cfeature.BORDERS, linestyle=':', zorder=3, edgecolor='gray')

        ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, zorder=4)

        # 4) Plot ground track
        ax.plot(subpoint_lon, subpoint_lat, 'r-', transform=ccrs.Geodetic(), label=f'{satellite.name} Ground Track', zorder=5)

        # 5) Mark observer
        obs_lat = observer_location.latitude.degrees
        obs_lon = observer_location.longitude.degrees
        ax.plot(obs_lon, obs_lat, 'bo', markersize=7, transform=ccrs.Geodetic(), label='Observer', zorder=6)
        ax.text(obs_lon + 2, obs_lat + 2, 'Observer', transform=ccrs.Geodetic(), zorder=7)
        
        # 6) Start/end markers
        ax.plot(subpoint_lon[0], subpoint_lat[0], 'go', transform=ccrs.Geodetic(), label='Track Start', zorder=5)
        ax.plot(subpoint_lon[-1], subpoint_lat[-1], 'mo', transform=ccrs.Geodetic(), label='Track End', zorder=5)


        # --- Footprint computation and drawing ---
        times_to_plot_footprint = []
        if footprint_times_utc and len(footprint_times_utc) > 0:
            times_to_plot_footprint = footprint_times_utc
        else:
            times_to_plot_footprint = [start_time_utc]

        footprint_colors = ['blue', 'green', 'purple', 'orange']
        footprint_alphas = [0.15, 0.12, 0.10, 0.08]

        for i, fp_time_utc in enumerate(times_to_plot_footprint):
            current_time_skyfield = ts.utc(fp_time_utc.year, fp_time_utc.month, fp_time_utc.day,
                                           fp_time_utc.hour, fp_time_utc.minute, 
                                           fp_time_utc.second + fp_time_utc.microsecond/1e6)
            geocentric_current = satellite.at(current_time_skyfield)
            subpoint_current = wgs84.subpoint(geocentric_current)
            current_lat_deg = subpoint_current.latitude.degrees
            current_lon_deg = subpoint_current.longitude.degrees
            current_alt_km = subpoint_current.elevation.km

            loop_status_message = f"時刻 {fp_time_utc.utc_datetime().strftime('%H:%M:%S')} UTC: "

            if current_alt_km > 0:
                try:
                    if min_elevation_for_footprint_deg > 0.0:
                        el_min_rad = np.deg2rad(min_elevation_for_footprint_deg)
                        cos_gamma_term = (R_EARTH_KM_MEAN / (R_EARTH_KM_MEAN + current_alt_km)) * np.cos(el_min_rad)
                        if cos_gamma_term > 1.0: cos_gamma_term = 1.0 
                        elif cos_gamma_term < -1.0: cos_gamma_term = -1.0
                        gamma = np.arccos(cos_gamma_term) - el_min_rad
                    else:
                        gamma = np.arccos(R_EARTH_KM_MEAN / (R_EARTH_KM_MEAN + current_alt_km))
                    
                    if gamma > 0:
                        footprint_radius_meters = R_EARTH_KM_MEAN * 1000 * gamma
                        circle_points = Geodesic().circle(lon=current_lon_deg, lat=current_lat_deg, 
                                                        radius=footprint_radius_meters, n_samples=180)
                        footprint_polygon = Polygon(circle_points)
                        
                        color_idx = i % len(footprint_colors)
                        alpha_idx = i % len(footprint_alphas)
                        
                        ax.add_geometries([footprint_polygon], ccrs.Geodetic(),
                                          facecolor=footprint_colors[color_idx], 
                                          alpha=footprint_alphas[alpha_idx], 
                                          edgecolor=footprint_colors[color_idx], 
                                          linewidth=0.8, zorder=2)
                        ax.plot(current_lon_deg, current_lat_deg, 
                                marker='o', color=footprint_colors[color_idx], markersize=5, markeredgecolor='black',
                                transform=ccrs.Geodetic(), 
                                label=f'Footprint @ {fp_time_utc.utc_datetime().strftime("%H:%M")}, El>={min_elevation_for_footprint_deg}°', 
                                zorder=6)
                        loop_status_message += "Footprint 繪製成功。"
                    else:
                        loop_status_message += f"Footprint 未繪製：衛星可能太低或仰角設定 ({min_elevation_for_footprint_deg}°) 過高。"
                except Exception as e_fp:
                    loop_status_message += f"計算或繪製 footprint 時發生錯誤: {e_fp}"
                    print(f"Error calculating or plotting footprint for time {fp_time_utc}: {e_fp}")
            else:
                loop_status_message += "Footprint 未繪製：衛星高度為零或負值。"
            footprint_status_messages.append(loop_status_message)
        # --- End of footprint ---

        plt.title(f"Satellite: {satellite.name} - Ground Track ({duration_hours} hours from {start_time_utc.utc_datetime().strftime('%Y-%m-%d %H:%M:%S')} UTC)")
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels)

        img_buffer = io.BytesIO()
        canvas.print_png(img_buffer)
        img_buffer.seek(0)
        pil_image = Image.open(img_buffer)
        plt.close(fig)

        combined_footprint_status = "<br>".join(footprint_status_messages) + "<br>"
        return pil_image, combined_footprint_status

    except Exception as e:
        error_message = f"繪製衛星軌跡圖時發生嚴重錯誤: {e}<br>"
        print(error_message)
        import traceback
        traceback.print_exc()
        return None, error_message

# --- Test code (kept as-is; comments translated above) ---
if __name__ == '__main__':
    from src import tle_fetcher
    from src import satellite_tracker
    import os

    # 1) Prepare TLE data
    tle_file_path = os.path.join(tle_fetcher.DEFAULT_TLE_SAVE_DIR, tle_fetcher.DEFAULT_TLE_FILENAME)
    if not os.path.exists(tle_file_path):
        print(f"TLE file not found at {tle_file_path}, downloading...")
        tle_fetcher.download_tle_data()

    # 2) Load satellites
    sats = satellite_tracker.load_satellites_from_tle(tle_file_path)
    if not sats:
        print("Failed to load satellites. Exiting test.")
        exit()

    iss = satellite_tracker.get_satellite_by_name(sats, "ISS (ZARYA)")
    if not iss:
        print("ISS not found in TLE data. Exiting test.")
        exit()
    print(f"Found satellite: {iss.name}")

    # 3) Set observer (e.g., Taipei)
    observer = satellite_tracker.get_observer_location(25.0330, 121.5654)
    print(f"Observer at: Lat {observer.latitude.degrees:.2f}, Lon {observer.longitude.degrees:.2f}")

    # 4) Set time
    current_time_utc = datetime.datetime.now(datetime.timezone.utc)
    print(f"Current UTC time: {current_time_utc.strftime('%Y-%m-%d %H:%M:%S')}")

    # 5) Plot ground track
    print("Generating ground track plot (Default Style)...")
    footprint_test_times = [
        current_time_utc,
        current_time_utc + datetime.timedelta(minutes=20),
        current_time_utc + datetime.timedelta(minutes=40)
    ]

    track_image_default, footprint_status_default = plot_satellite_track(
        iss, observer, current_time_utc, 
        duration_hours=3, points_per_hour=60, map_style='Default',
        footprint_times_utc=footprint_test_times
    )
    if track_image_default:
        print("Default plot generated. Saving to test_track_plot_default.png")
        track_image_default.save("test_track_plot_default.png")
    print(footprint_status_default)

    print("\\nGenerating ground track plot (Stamen Terrain Style)...")
    track_image_stamen, footprint_status_stamen = plot_satellite_track(
        iss, observer, current_time_utc, 
        duration_hours=3, points_per_hour=60, map_style='Stamen Terrain',
        footprint_times_utc=footprint_test_times
    )
    if track_image_stamen:
        print("Stamen Terrain plot generated. Saving to test_track_plot_stamen.png")
        track_image_stamen.save("test_track_plot_stamen.png")
    print(footprint_status_stamen)

    print("\\nGenerating ground track plot (OpenStreetMap Style)...")
    track_image_osm, footprint_status_osm = plot_satellite_track(
        iss, observer, current_time_utc, 
        duration_hours=3, points_per_hour=60, map_style='OpenStreetMap',
        footprint_times_utc=footprint_test_times
    )
    if track_image_osm:
        print("OpenStreetMap plot generated. Saving to test_track_plot_osm.png")
        track_image_osm.save("test_track_plot_osm.png")
    print(footprint_status_osm)

    # --- Test GIF Generation ---
    print("\nGenerating Footprint Animation GIF...")
    gif_start_time = current_time_utc
    
    if not os.path.exists('data'):
        os.makedirs('data')

    gif_path = generate_footprint_animation_gif(
        satellite=iss,
        observer_location=observer,
        animation_start_time_utc=gif_start_time,
        animation_duration_minutes=1,
        simulation_time_step_seconds=10,
        gif_playback_fps=2,
        output_gif_path="data/test_footprint_animation_detailed.gif",
        base_track_duration_hours=3, 
        min_elevation_for_footprint_deg=5.0,
        map_style='Default',
        fig_size_inches=(7,5) 
    )

    if gif_path:
        print(f"Footprint animation GIF generated: {gif_path}")
    else:
        print("Failed to generate footprint animation GIF.")

    # --- Test Folium Interactive Map Generation ---
    print("\\nGenerating Folium Interactive Map...")
    folium_map_html, _ = generate_interactive_map_folium(
        satellite=iss,
        observer_location=observer,
        track_reference_time_utc=current_time_utc,
        display_time_utc=current_time_utc,
        track_duration_hours=3,
        points_per_hour=60,
        min_elevation_for_footprint_deg=5.0,
        map_tiles='OpenStreetMap' 
    )

    if folium_map_html:
        map_file_path = "data/test_folium_map.html"
        os.makedirs(os.path.dirname(map_file_path), exist_ok=True)
        with open(map_file_path, "w", encoding="utf-8") as f:
            f.write(folium_map_html)
        print(f"Folium interactive map generated and saved to: {map_file_path}")
    else:
        print("Failed to generate Folium interactive map.")

    print("\\nVisualizer tests completed.")
