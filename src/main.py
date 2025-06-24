import gradio as gr
import os
from datetime import datetime, timedelta, timezone
import pytz
import pandas as pd
from skyfield.api import load, wgs84
from timezonefinder import TimezoneFinder
import tempfile
import csv
import io
import json
import openpyxl
import logging
import traceback
import hashlib

from src.tle_fetcher import download_tle_file, DEFAULT_TLE_FILE, get_tle_sources, TLEFetchError
from src.satellite_tracker import TLE_FORMAT_ERROR_MESSAGE
from src import tle_fetcher, satellite_tracker, visualizer

# =========================
# Constants
# =========================
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
DEFAULT_TLE_URL = "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle"

# Cartopy-based static map styles (used by GIF generation)
STATIC_MAP_STYLE_CHOICES = ["cartopy.mpl.geoaxes.GeoAxes.stock_img", "cartopy.mpl.geoaxes.GeoAxes.nightshade", "cartopy.mpl.geoaxes.GeoAxes.blue_marble"]
DEFAULT_STATIC_MAP_STYLE = "cartopy.mpl.geoaxes.GeoAxes.stock_img"

# Folium tiles for interactive map
FOLIUM_MAP_TILES_CHOICES = ["OpenStreetMap", "Stamen Terrain", "Stamen Toner", "Stamen Watercolor", "CartoDB positron", "CartoDB dark_matter"]
DEFAULT_FOLIUM_MAP_TILES = "OpenStreetMap"

DEFAULT_MIN_ELEVATION = 10.0
DEFAULT_FOLIUM_INITIAL_ZOOM = 1

# Skyfield timescale
ts = load.timescale()

# Default empty HTML for the map output
EMPTY_MAP_CONTENT = "<html><body><p style=\\'text-align: center; margin-top: 20px;\\'>地圖將在此處顯示。</p></body></html>"

# Fallback Folium tile list if visualizer does not export AVAILABLE_MAP_STYLES
DEFAULT_FOLIUM_TILE_LIST = [
    "OpenStreetMap",
    "CartoDB positron",
    "CartoDB dark_matter",
    "Stamen Terrain",
    "Stamen Toner",
    "Stamen Watercolor"
]

# If visualizer provides styles, use them; otherwise fallback
if hasattr(visualizer, 'AVAILABLE_MAP_STYLES') and visualizer.AVAILABLE_MAP_STYLES:
    FOLIUM_MAP_TILES_CHOICES = visualizer.AVAILABLE_MAP_STYLES
    DEFAULT_FOLIUM_MAP_TILES = FOLIUM_MAP_TILES_CHOICES[0]
else:
    FOLIUM_MAP_TILES_CHOICES = DEFAULT_FOLIUM_TILE_LIST
    DEFAULT_FOLIUM_MAP_TILES = FOLIUM_MAP_TILES_CHOICES[0]
    print("Warning: visualizer.AVAILABLE_MAP_STYLES not found or empty. Using default Folium tile list for main.py.")

# Static map styles for GIF generation (Cartopy based) — simplified labels for end users
STATIC_MAP_STYLE_CHOICES = ["Default", "OpenStreetMap", "Stamen Terrain"]
DEFAULT_STATIC_MAP_STYLE = STATIC_MAP_STYLE_CHOICES[0]

# Default inputs
DEFAULT_SATELLITE_NAME = "ISS (ZARYA)"
DEFAULT_OBSERVER_LAT = 25.0330
DEFAULT_OBSERVER_LON = 121.5654
DEFAULT_MIN_ELEVATION_FOOTPRINT = 10.0
DEFAULT_GIF_DURATION_MINUTES = 5
DEFAULT_GIF_FPS = 10

# =========================
# Logger Setup
# =========================
LOG_FILE_PATH = os.path.join(DATA_DIR, "app.log")

# Ensure data directory exists
if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)

logger = logging.getLogger(__name__)

# =========================
# Global App State
# =========================
tle_data = []
app_state = {}    # shared state across Gradio callbacks

# =========================
# Helper Functions
# =========================
def format_current_position(pos_data, sat_name):
    """
    Build an HTML snippet for the current satellite position.

    Parameters
    ----------
    pos_data : dict | None
        Position dictionary returned by satellite_tracker.get_satellite_current_position.
    sat_name : str
        Satellite name for display.

    Returns
    -------
    str
        HTML-formatted string describing the current position or an error message.
    """
    if not pos_data:
        return f"無法計算衛星 {sat_name} 的目前位置。"
    # Use <br> for line breaks because Gradio Markdown respects basic HTML
    return (
        f"衛星 {sat_name} 目前位置:<br>"
        f"  星下點: 緯度 {pos_data['subpoint_lat_deg']:.2f}°, 經度 {pos_data['subpoint_lon_deg']:.2f}°<br>"
        f"  軌道高度: {pos_data['altitude_km']:.2f} km<br>"
        f"  相對於觀測點:<br>"
        f"    方位角: {pos_data['azimuth_deg']:.2f}°<br>"
        f"    仰角: {pos_data['elevation_deg']:.2f}°<br>"
        f"    距離: {pos_data['distance_km']:.2f} km"
    )

def format_pass_predictions(pass_data, sat_name, observer_lat, observer_lon):
    """
    Build an HTML table for pass predictions with local time conversion.

    Parameters
    ----------
    pass_data : list[dict]
        List of predicted pass events with UTC timestamps and angles.
    sat_name : str
        Satellite name for header text.
    observer_lat : float
        Observer latitude in degrees.
    observer_lon : float
        Observer longitude in degrees.

    Returns
    -------
    str
        HTML content with a header and a table of pass predictions, or a message if none exist.
    """
    if not pass_data:
        return f"在接下來的24小時內，未找到 {sat_name} 高於10°仰角的過境事件。"

    tf = TimezoneFinder()
    observer_tz_str = tf.timezone_at(lng=observer_lon, lat=observer_lat)
    observer_tz = pytz.utc
    if observer_tz_str:
        try:
            observer_tz = pytz.timezone(observer_tz_str)
        except pytz.exceptions.UnknownTimeZoneError:
            observer_tz_str = "UTC"

    header_message = f"找到 {len(pass_data)} 次 {sat_name} 的過境事件 (時間已轉換為 {observer_tz_str}):<br><br>"

    table_data = []
    column_headers = [
        f"<span style='white-space: nowrap;'>上升時間</span>", 
        f"<span style='white-space: nowrap;'>上升方位</span>", 
        f"<span style='white-space: nowrap;'>最高點時間</span>", 
        f"<span style='white-space: nowrap;'>最高點方位</span>", 
        f"<span style='white-space: nowrap;'>最高點仰角</span>", 
        f"<span style='white-space: nowrap;'>下降時間</span>", 
        f"<span style='white-space: nowrap;'>下降方位</span>", 
        f"<span style='white-space: nowrap;'>持續時間</span>"
    ]
    for p_info in pass_data:
        rise_time_local = p_info['rise_time_utc'].astimezone(observer_tz)
        cul_time_local = p_info['culmination_time_utc'].astimezone(observer_tz)
        set_time_local = p_info['set_time_utc'].astimezone(observer_tz)
        table_data.append([
            f"<span style='white-space: nowrap;'>{rise_time_local.strftime('%Y-%m-%d %H:%M:%S')}</span>",
            f"<span style='white-space: nowrap;'>{p_info['rise_azimuth_deg']:.1f}&nbsp;&deg;</span>",
            f"<span style='white-space: nowrap;'>{cul_time_local.strftime('%Y-%m-%d %H:%M:%S')}</span>",
            f"<span style='white-space: nowrap;'>{p_info['culmination_azimuth_deg']:.1f}&nbsp;&deg;</span>",
            f"<span style='white-space: nowrap;'>{p_info['culmination_elevation_deg']:.1f}&nbsp;&deg;</span>",
            f"<span style='white-space: nowrap;'>{set_time_local.strftime('%Y-%m-%d %H:%M:%S')}</span>",
            f"<span style='white-space: nowrap;'>{p_info['set_azimuth_deg']:.1f}&nbsp;&deg;</span>",
            f"<span style='white-space: nowrap;'>{p_info['duration_minutes']:.1f} 分鐘</span>"
        ])

    if not table_data:
        return header_message

    df = pd.DataFrame(table_data, columns=column_headers)
    # escape=False is required so the span/HTML is rendered properly by Gradio
    html_table = df.to_html(index=False, border=1, escape=False, classes='passes-table') 
    return header_message + html_table

def update_map_for_time_slider(slider_minutes_offset, current_app_state):
    """
    Update the interactive Folium map when the time slider changes.

    The function uses cached state (satellite object, observer, track reference time,
    map style, zoom, center) to render the map at a specific time offset.

    Parameters
    ----------
    slider_minutes_offset : int | float
        Minutes offset from the stored `track_reference_time_utc`.
    current_app_state : dict
        The application state persisted across callbacks.

    Returns
    -------
    tuple[str, str, dict]
        (folium_map_html_output, status_log, updated_app_state)
    """
    status_log = ""
    folium_map_html_output = None
    updated_app_state = current_app_state.copy() if current_app_state is not None else {}

    try:
        logger.info(f"update_map_for_time_slider: Called. Slider offset: {slider_minutes_offset} min.")
        status_log += f"DEBUG: Slider triggered with offset: {slider_minutes_offset} minutes.<br>"
        logger.debug(f"update_map_for_time_slider: Current app_state keys: {list(updated_app_state.keys())}")

        ts_from_state = updated_app_state.get("ts")
        if not ts_from_state:
            status_log += "錯誤：內部時間設定(ts)遺失，無法更新地圖。<br>"
            logger.error("update_map_for_time_slider: Timescale 'ts' not found in app_state.")
            folium_map_html_output = updated_app_state.get("last_folium_map_html", "<p>Error: Timescale missing in state.</p>")
            return folium_map_html_output, status_log, updated_app_state

        satellite_obj = updated_app_state.get("current_satellite_skyfield_obj")
        if not satellite_obj:
            status_log += "錯誤：衛星物件(current_satellite_skyfield_obj)遺失，無法更新地圖。<br>"
            logger.error("update_map_for_time_slider: 'current_satellite_skyfield_obj' not found in app_state.")
            folium_map_html_output = updated_app_state.get("last_folium_map_html", "<p>Error: Satellite object missing in state.</p>")
            return folium_map_html_output, status_log, updated_app_state
        
        logger.debug(f"update_map_for_time_slider: ts_from_state type: {type(ts_from_state)}, satellite_obj name: {getattr(satellite_obj, 'name', 'N/A')}")

        observer_obj = updated_app_state.get("current_observer_topos_obj")
        track_ref_time_utc_obj = updated_app_state.get("track_reference_time_utc")      # Skyfield Time expected
        min_elev = updated_app_state.get("min_elevation_for_footprint_deg", DEFAULT_MIN_ELEVATION_FOOTPRINT)
        map_tiles = updated_app_state.get("folium_map_tiles", DEFAULT_FOLIUM_MAP_TILES)
        track_duration_hrs = updated_app_state.get("track_duration_hours", 3.0)

        if observer_obj is None or track_ref_time_utc_obj is None:
            status_log += "錯誤：應用程式狀態不完整 (observer/track_ref_time)，無法更新地圖。<br>"
            logger.warning("update_map_for_time_slider: Incomplete state (observer_obj or track_ref_time_utc_obj) for map update.")
            folium_map_html_output = updated_app_state.get("last_folium_map_html", "<p>Error: Observer or track reference time missing.</p>")
            return folium_map_html_output, status_log, updated_app_state

        # Reference time is already a Skyfield Time object.
        track_reference_time_utc = track_ref_time_utc_obj

        # Compute the display time from reference + slider offset
        current_display_time_utc = track_reference_time_utc + timedelta(minutes=slider_minutes_offset)
        
        status_log += f"DEBUG: Map Update Triggered.<br>"
        status_log += f"DEBUG: Slider Offset: {slider_minutes_offset} minutes.<br>"
        status_log += f"DEBUG: Track Reference UTC: {track_reference_time_utc.utc_iso()}<br>"
        status_log += f"DEBUG: Calculated Display UTC for Map: {current_display_time_utc.utc_iso()}<br>"
        logger.info(f"update_map_for_time_slider: Calculated Display UTC for Map: {current_display_time_utc.utc_iso()}")

        new_map_center = None
        if hasattr(satellite_obj, 'at'): 
            logger.debug(f"update_map_for_time_slider: Calculating new subpoint for satellite {satellite_obj.name} at {current_display_time_utc.utc_iso()}")
            # Compute satellite subpoint to center the map
            geocentric_at_display = satellite_obj.at(current_display_time_utc)
            subpoint_at_display = wgs84.subpoint(geocentric_at_display)
            new_center_lat = subpoint_at_display.latitude.degrees
            new_center_lon = subpoint_at_display.longitude.degrees
            new_map_center = (new_center_lat, new_center_lon)
            status_log += f"DEBUG: New Map Center (Sat Subpoint): ({new_center_lat:.2f}, {new_center_lon:.2f}).<br>"
            logger.info(f"Slider: New map center at satellite subpoint: ({new_center_lat:.2f}, {new_center_lon:.2f})")
        else:
            new_map_center = updated_app_state.get("current_folium_center") 
            status_log += "DEBUG: Could not calculate new satellite subpoint for centering, using previous center.<br>"
            logger.warning("Slider: Could not get satellite_obj or it's not a valid Skyfield object to calculate new subpoint for centering.")

        persistent_zoom_level = updated_app_state.get("current_folium_zoom", DEFAULT_FOLIUM_INITIAL_ZOOM)
        status_log += f"DEBUG: Using Zoom Level: {persistent_zoom_level}.<br>"
        logger.info(f"Slider: Using zoom level: {persistent_zoom_level}")
        
        logger.info(f"update_map_for_time_slider: Calling generate_interactive_map_folium with display_time_utc: {current_display_time_utc.utc_iso()}, center: {new_map_center}, zoom: {persistent_zoom_level}")
        
        folium_map_html_output, newly_computed_track_slider = visualizer.generate_interactive_map_folium(
            satellite=satellite_obj,
            observer_location=observer_obj,
            track_reference_time_utc=track_reference_time_utc,  # global reference for the track
            display_time_utc=current_display_time_utc,          # current marker/footprint time
            track_duration_hours=track_duration_hrs,
            min_elevation_for_footprint_deg=min_elev,
            map_tiles=map_tiles,
            cached_track_data=updated_app_state.get('cached_track_data'),
            perform_fit_bounds=False, 
            center_location=new_map_center,
            current_zoom_level=persistent_zoom_level
        )
        
        if newly_computed_track_slider:
            # Track cache recomputed or refreshed
            logger.info("Slider update recomputed track data. Cache might not have been used as expected or was refreshed.")
            updated_app_state['cached_track_data'] = newly_computed_track_slider
        
        if folium_map_html_output and len(folium_map_html_output) > 100:
            html_hash = hashlib.md5(folium_map_html_output.encode()).hexdigest()
            status_log += f"DEBUG: Generated Map HTML. Hash: {html_hash}. Length: {len(folium_map_html_output)}.<br>"
            logger.info(f"update_map_for_time_slider: Generated Map HTML. Hash: {html_hash}, Length: {len(folium_map_html_output)}")
        else:
            status_log += "DEBUG: Generated Map HTML is empty or very short.<br>"
            logger.warning("update_map_for_time_slider: Generated Map HTML is empty or very short.")

        status_log += "地圖已根據時間滑桿更新（固定縮放/中心更新）。<br>"
        logger.info("update_map_for_time_slider: Folium map processing complete for slider.")
        updated_app_state["last_folium_map_html"] = folium_map_html_output
        if new_map_center: 
            updated_app_state["current_folium_center"] = new_map_center
        
    except Exception as e:
        error_msg = f"更新地圖時發生嚴重錯誤: {type(e).__name__}: {e}"
        status_log += f"嚴重錯誤：{error_msg}<br>"
        logger.error(f"update_map_for_time_slider: {error_msg}", exc_info=True) 
        folium_map_html_output = updated_app_state.get("last_folium_map_html", f"<p>Error updating map: {e}</p>")

    return folium_map_html_output, status_log, updated_app_state

def ground_station_interface(
    satellite_name_query,
    observer_lat,
    observer_lon,
    tle_source,             # category, e.g., Active Satellites, Space Stations
    tle_mode,               # how to get TLE: auto download / upload file / manual paste
    uploaded_tle_file_obj, 
    manual_tle_lines, 
    selected_map_style, 
    min_elevation_deg, 
    app_state_param
):
    """
    Main entry point for the Gradio UI: validates inputs, resolves TLE,
    finds the target satellite, computes current position and next passes,
    and renders the interactive Folium map.

    Returns
    -------
    tuple
        (current_pos_str_val, passes_str_val, status_log, folium_map_html_output,
         slider_update, app_state)
    """
    status_log = ""
    app_state = app_state_param.copy() if app_state_param is not None else {}
    temp_manual_tle_path = None

    # Initialize outputs for error paths
    current_pos_str_val = ""
    passes_str_val = ""
    folium_map_html_output = EMPTY_MAP_CONTENT
    slider_update = gr.update(visible=False, label="時間軸")

    try:
        logger.info(f"ground_station_interface: Request for '{satellite_name_query}', Lat: '{observer_lat}', Lon: '{observer_lon}', TLE Src: '{tle_source}', Map: '{selected_map_style}', MinEl: {min_elevation_deg}°")
        # Load timescale and store into app_state for other handlers (e.g., slider)
        ts = load.timescale()
        app_state["ts"] = ts

        # ---- Input validation
        if not satellite_name_query:
            status_log += "錯誤：請輸入衛星名稱。<br>"
            logger.warning("ground_station_interface: Satellite name missing.")
            return current_pos_str_val, passes_str_val, status_log, folium_map_html_output, slider_update, app_state

        if observer_lat is None or observer_lon is None:
            status_log += "錯誤：請提供觀測者緯度和經度。<br>"
            logger.warning(f"ground_station_interface: Observer lat or lon is None.")
            return current_pos_str_val, passes_str_val, status_log, folium_map_html_output, slider_update, app_state

        try:
            valid_observer_lat = float(observer_lat)
            valid_observer_lon = float(observer_lon)
            if not (-90 <= valid_observer_lat <= 90 and -180 <= valid_observer_lon <= 180):
                raise ValueError("緯度或經度超出範圍")
            # Build Skyfield Topos and keep both numeric coords and object in app_state
            observer_topos_obj = wgs84.latlon(valid_observer_lat, valid_observer_lon)
            app_state["current_observer_lat"] = valid_observer_lat
            app_state["current_observer_lon"] = valid_observer_lon
            app_state['current_observer_topos_obj'] = observer_topos_obj
        except ValueError as e:
            status_log += f"錯誤：無效的觀測者座標：{e}<br>"
            logger.error(f"ground_station_interface: Invalid observer coordinates: {e}")
            return current_pos_str_val, passes_str_val, status_log, folium_map_html_output, slider_update, app_state

        try:
            min_elevation_val = float(min_elevation_deg if min_elevation_deg is not None else 0.0)
            if not (0 <= min_elevation_val <= 90):
                raise ValueError("最小仰角必須介於 0 和 90 度之間。")
        except ValueError as e:
            status_log += f"錯誤：無效的最小仰角設定：{e}<br>"
            logger.error(f"ground_station_interface: Invalid min_elevation_deg: {e}")
            return current_pos_str_val, passes_str_val, status_log, folium_map_html_output, slider_update, app_state

        # ---- TLE handling
        tle_file_to_use = None
        status_message_tle = ""

        if tle_mode == "自動下載最新 TLE":
            try:
                # Ensure a download directory exists
                download_dir = os.path.join(DATA_DIR, "downloaded_tle")
                os.makedirs(download_dir, exist_ok=True)
                # Download the TLE for the selected category
                safe_source_name = "".join(c if c.isalnum() else "_" for c in tle_source).lower()
                download_path = os.path.join(download_dir, f"{safe_source_name}.tle")
                tle_file_to_use = tle_fetcher.download_tle_file(tle_source_name=tle_source, tle_file_path=download_path)
                status_message_tle = f"已從網路下載 {tle_source} TLE 檔案: {os.path.basename(tle_file_to_use)}<br>"
            except tle_fetcher.TLEFetchError as e:
                status_message_tle = f"下載 {tle_source} TLE 失敗: {e}<br>"
                logger.error(f"TLE download failed for {tle_source}: {e}")
                # Fallback: try local cached file; if none, fallback to Active Satellites
                safe_source_name = "".join(c if c.isalnum() else "_" for c in tle_source).lower()
                local_tle_path = os.path.join(DATA_DIR, "downloaded_tle", f"{safe_source_name}.tle")
                if os.path.exists(local_tle_path):
                    tle_file_to_use = local_tle_path
                    status_message_tle += f"使用本地快取的 TLE: {os.path.basename(tle_file_to_use)}<br>"
                else:
                    try:
                        download_path = os.path.join(DATA_DIR, "downloaded_tle", "active_satellites.tle")
                        tle_file_to_use = tle_fetcher.download_tle_file(tle_source_name="Active Satellites", tle_file_path=download_path)
                        status_message_tle += f"回退下載 Active Satellites TLE: {os.path.basename(tle_file_to_use)}<br>"
                    except tle_fetcher.TLEFetchError:
                        tle_file_to_use = None
                        status_message_tle += "無法下載任何 TLE 檔案。<br>"
        elif tle_mode == "上傳 TLE 檔案" and uploaded_tle_file_obj is not None:
            # The uploaded file path is provided by Gradio
            tle_file_to_use = uploaded_tle_file_obj.name
            status_message_tle = f"使用上傳的 TLE 檔案: {os.path.basename(tle_file_to_use)}<br>"
        elif tle_mode == "手動輸入 TLE" and manual_tle_lines and manual_tle_lines.strip():
            # Persist manual lines to a temporary file so downstream code can treat it uniformly
            try:
                with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tle", dir=DATA_DIR) as tmpfile:
                    tmpfile.write(manual_tle_lines)
                    temp_manual_tle_path = tmpfile.name
                tle_file_to_use = temp_manual_tle_path
                status_message_tle = "使用手動輸入的 TLE 資料。<br>"
            except Exception as e:
                status_message_tle = f"處理手動 TLE 輸入時發生錯誤: {e}<br>"
                logger.error(f"Error processing manual TLE: {e}")
                tle_file_to_use = None
        else:
            # If no explicit mode provided, try auto download for the selected category
            try:
                status_message_tle += f"嘗試自動下載 {tle_source} TLE 檔案...<br>"
                download_dir = os.path.join(DATA_DIR, "downloaded_tle")
                os.makedirs(download_dir, exist_ok=True)
                safe_source_name = "".join(c if c.isalnum() else "_" for c in tle_source).lower()
                download_path = os.path.join(download_dir, f"{safe_source_name}.tle")
                tle_file_to_use = tle_fetcher.download_tle_file(tle_source_name=tle_source, tle_file_path=download_path)
                status_message_tle += f"已自動下載 {tle_source} TLE 檔案: {os.path.basename(tle_file_to_use)}<br>"
            except tle_fetcher.TLEFetchError as e:
                status_message_tle += f"自動下載 {tle_source} TLE 失敗: {e}<br>"
                logger.error(f"Auto-download TLE failed for {tle_source}: {e}")
                tle_file_to_use = None
        
        status_log += status_message_tle
        if not tle_file_to_use or not os.path.exists(tle_file_to_use):
            status_log += f"錯誤：TLE 檔案 '{tle_file_to_use}' 不存在或無法存取。<br>"
            logger.error(f"TLE file not accessible: {tle_file_to_use}")
            if temp_manual_tle_path and os.path.exists(temp_manual_tle_path): os.remove(temp_manual_tle_path)
            return current_pos_str_val, passes_str_val, status_log, folium_map_html_output, slider_update, app_state
        
        app_state["current_tle_file"] = tle_file_to_use
        logger.info(f"Using TLE file: {tle_file_to_use}")

        try:
            # Load satellites from TLE
            sats = satellite_tracker.load_satellites_from_tle(tle_file_to_use)
            status_log += f"已從 {os.path.basename(tle_file_to_use)} 成功載入 {len(sats)} 顆衛星。<br>"
            logger.info(f"Loaded {len(sats)} satellites from {tle_file_to_use}")
        except ValueError as e:
            if TLE_FORMAT_ERROR_MESSAGE in str(e):
                status_log += f"錯誤：TLE 檔案格式無效。{TLE_FORMAT_ERROR_MESSAGE}<br>"
                logger.error(f"Invalid TLE file format: {tle_file_to_use}. Error: {e}")
            else:
                status_log += f"載入 TLE 檔案時發生錯誤: {e}<br>"
                logger.error(f"Error loading TLE file {tle_file_to_use}: {e}")
            if temp_manual_tle_path and os.path.exists(temp_manual_tle_path): os.remove(temp_manual_tle_path)
            return current_pos_str_val, passes_str_val, status_log, folium_map_html_output, slider_update, app_state

        # Find the target satellite by name (fuzzy or exact depending on your implementation)
        target_satellite = satellite_tracker.get_satellite_by_name(sats, satellite_name_query)
        if not target_satellite:
            status_log += f"錯誤：在 TLE 檔案 {os.path.basename(tle_file_to_use)} 中找不到名為 '{satellite_name_query}' 的衛星。<br>"
            logger.warning(f"ground_station_interface: Satellite '{satellite_name_query}' not found in {os.path.basename(tle_file_to_use)}.")
            available_sats = [s.name for s in sats[:20]]
            status_log += f"可用的衛星名稱範例: {', '.join(available_sats)}{'...' if len(sats) > 20 else ''}<br>"
            if temp_manual_tle_path and os.path.exists(temp_manual_tle_path): os.remove(temp_manual_tle_path)
            return current_pos_str_val, passes_str_val, status_log, folium_map_html_output, slider_update, app_state
        
        status_log += f"成功找到衛星: {target_satellite.name}。<br>"
        logger.info(f"ground_station_interface: Found satellite: {target_satellite.name}")

        # Reference time for tracking/slider
        current_time_utc = ts.now()
        app_state["track_reference_time_utc"] = current_time_utc

        # Cache invalidation if key parameters changed
        previous_satellite_name = app_state.get('satellite_name')
        previous_observer_lat = app_state.get('current_observer_lat')
        previous_observer_lon = app_state.get('current_observer_lon')
        previous_tle_file = app_state.get('current_tle_file')
        
        cache_invalidation_reasons = []
        if previous_satellite_name and previous_satellite_name != satellite_name_query:
            cache_invalidation_reasons.append(f"衛星變更: {previous_satellite_name} -> {satellite_name_query}")
        if previous_observer_lat is not None and abs(previous_observer_lat - valid_observer_lat) > 0.001:
            cache_invalidation_reasons.append(f"觀測者緯度變更: {previous_observer_lat:.3f}° -> {valid_observer_lat:.3f}°")
        if previous_observer_lon is not None and abs(previous_observer_lon - valid_observer_lon) > 0.001:
            cache_invalidation_reasons.append(f"觀測者經度變更: {previous_observer_lon:.3f}° -> {valid_observer_lon:.3f}°")
        if previous_tle_file and previous_tle_file != tle_file_to_use:
            cache_invalidation_reasons.append(f"TLE 檔案變更: {os.path.basename(previous_tle_file)} -> {os.path.basename(tle_file_to_use)}")
        
        if cache_invalidation_reasons:
            if 'cached_track_data' in app_state:
                del app_state['cached_track_data']
                logger.info(f"Cleared cached track data due to parameter changes: {'; '.join(cache_invalidation_reasons)}")
                status_log += f"已清除軌跡快取，原因: {'; '.join(cache_invalidation_reasons)}。<br>"

        # Persist objects and metadata in app_state for later handlers (slider/GIF)
        app_state['current_satellite_skyfield_obj'] = target_satellite
        app_state['current_observer_topos_obj'] = observer_topos_obj
        
        # Store TLE lines (for debugging/export)
        tle_lines_content = None
        if hasattr(target_satellite, 'tle') and isinstance(target_satellite.tle, (list, tuple)) and len(target_satellite.tle) >= 2:
            tle_lines_content = "\\n".join(target_satellite.tle)
        elif hasattr(target_satellite, 'tle_lines') and isinstance(target_satellite.tle_lines, str):
             tle_lines_content = target_satellite.tle_lines
        elif tle_file_to_use and os.path.exists(tle_file_to_use):
            try:
                with open(tle_file_to_use, 'r') as f_tle:
                    tle_lines_content = f_tle.read()
            except Exception as e_read:
                logger.warning(f"Could not read TLE lines from {tle_file_to_use}: {e_read}")
        app_state['tle_lines'] = tle_lines_content

        app_state['satellite_name'] = satellite_name_query

        # ---- Current position
        logger.info("Calculating current satellite position...")
        current_position = satellite_tracker.get_satellite_current_position(target_satellite, observer_topos_obj)
        if current_position:
            current_pos_str_val = (
                f"<b>目前衛星位置 ({target_satellite.name})</b><br>"
                f"方位角: {current_position['azimuth_deg']:.1f}°<br>"
                f"仰角: {current_position['elevation_deg']:.1f}°<br>"
                f"距離: {current_position['distance_km']:.1f} km<br>"
                f"星下點緯度: {current_position['subpoint_lat_deg']:.4f}°<br>"
                f"星下點經度: {current_position['subpoint_lon_deg']:.4f}°<br>"
                f"衛星高度: {current_position['altitude_km']:.1f} km<br>"
                f"<small>更新時間: {current_time_utc.utc_datetime().strftime('%H:%M:%S')} UTC</small>"
            )
            logger.info(f"Current satellite position calculated: Az={current_position['azimuth_deg']:.1f}°, El={current_position['elevation_deg']:.1f}°")
        else:
            current_pos_str_val = "<b>目前衛星位置</b><br>無法計算衛星位置"
            logger.warning("Failed to calculate current satellite position")

        # ---- Pass predictions (next 24 hours)
        logger.info("Calculating satellite pass predictions...")
        prediction_start_time = current_time_utc.utc_datetime()
        prediction_duration_hours = 24.0
        passes = satellite_tracker.predict_satellite_passes(
            target_satellite, 
            observer_topos_obj, 
            prediction_start_time, 
            prediction_duration_hours, 
            min_elevation_val
        )
        
        if passes:
            passes_str_val = f"<b>衛星過境預報 ({target_satellite.name})</b><br>"
            passes_str_val += f"<small>預報時間: 未來 {prediction_duration_hours:.0f} 小時 (最小仰角: {min_elevation_val}°)</small><br><br>"
            
            passes_str_val += """
            <table border="1" cellpadding="4" cellspacing="0" style="border-collapse: collapse; font-size: 12px; width: 100%;">
                <thead>
                    <tr>
                        <th style="white-space: nowrap; text-align: center;">#</th>
                        <th style="white-space: nowrap; text-align: center;">上升時間 (UTC)</th>
                        <th style="white-space: nowrap; text-align: center;">上升方位角</th>
                        <th style="white-space: nowrap; text-align: center;">最高時間 (UTC)</th>
                        <th style="white-space: nowrap; text-align: center;">最高仰角</th>
                        <th style="white-space: nowrap; text-align: center;">下降時間 (UTC)</th>
                        <th style="white-space: nowrap; text-align: center;">下降方位角</th>
                        <th style="white-space: nowrap; text-align: center;">持續時間</th>
                    </tr>
                </thead>
                <tbody>
            """
            
            for i, pass_info in enumerate(passes[:10]):
                passes_str_val += f"""
                    <tr>
                        <td style="white-space: nowrap; text-align: center;">{i+1}</td>
                        <td style="white-space: nowrap; text-align: center;">{pass_info['rise_time_utc'].strftime('%m/%d %H:%M')}</td>
                        <td style="white-space: nowrap; text-align: center;">{pass_info['rise_azimuth_deg']:.0f}°</td>
                        <td style="white-space: nowrap; text-align: center;">{pass_info['culmination_time_utc'].strftime('%m/%d %H:%M')}</td>
                        <td style="white-space: nowrap; text-align: center;">{pass_info['culmination_elevation_deg']:.0f}°</td>
                        <td style="white-space: nowrap; text-align: center;">{pass_info['set_time_utc'].strftime('%m/%d %H:%M')}</td>
                        <td style="white-space: nowrap; text-align: center;">{pass_info['set_azimuth_deg']:.0f}°</td>
                        <td style="white-space: nowrap; text-align: center;">{pass_info['duration_minutes']:.1f} 分鐘</td>
                    </tr>
                """
            
            passes_str_val += """
                </tbody>
            </table>
            """
            
            if len(passes) > 10:
                passes_str_val += f"<br><small>... 還有 {len(passes) - 10} 次過境 (下載完整資料以查看全部)</small>"
            
            logger.info(f"Found {len(passes)} satellite passes in next {prediction_duration_hours} hours")
        else:
            passes_str_val = f"<b>衛星過境預報 ({target_satellite.name})</b><br>"
            passes_str_val += f"<small>未來 {prediction_duration_hours:.0f} 小時內無可見過境</small><br>"
            passes_str_val += f"<small>(最小仰角: {min_elevation_val}°)</small>"
            logger.info(f"No satellite passes found in next {prediction_duration_hours} hours with min elevation {min_elevation_val}°")

        # ---- Interactive Folium map (initial render)
        folium_map_html_output, newly_computed_track_data_main = visualizer.generate_interactive_map_folium(
            satellite=target_satellite,
            observer_location=observer_topos_obj,
            track_reference_time_utc=current_time_utc,
            display_time_utc=current_time_utc,
            track_duration_hours=app_state.get("track_duration_hours", 3.0),
            min_elevation_for_footprint_deg=min_elevation_val,
            map_tiles=selected_map_style,
            cached_track_data=app_state.get('cached_track_data'),
            perform_fit_bounds=True,
            center_location=app_state.get("current_folium_center"),
            current_zoom_level=app_state.get("current_folium_zoom")
        )
        if newly_computed_track_data_main:
            app_state['cached_track_data'] = newly_computed_track_data_main
        app_state["last_folium_map_html"] = folium_map_html_output
        app_state["folium_map_tiles"] = selected_map_style
        app_state["min_elevation_for_footprint_deg"] = min_elevation_val

        # Slider range: 0 to (track_duration_hours * 60) minutes
        slider_max_minutes = int(app_state.get("track_duration_hours", 3.0) * 60)
        slider_update = gr.Slider(
            minimum=0,
            maximum=slider_max_minutes,
            value=0,
            step=5,
            label=f"時間軸 (相對於 {current_time_utc.utc_datetime().strftime('%H:%M')} UTC, 0 至 {slider_max_minutes} 分鐘)",
            visible=True
        )
        status_log += "互動式地圖已產生。<br>"
        
        if passes:
            app_state['passes_data'] = passes
        else:
            status_log += "無過境預報資料可供下載。<br>"
        
    except Exception as e:
        error_msg_detail = f"處理請求時發生未預期的錯誤: {type(e).__name__}: {e}"
        status_log += f"嚴重錯誤：{error_msg_detail}。請查看日誌。<br>"
        logger.error(f"ground_station_interface: UNEXPECTED ERROR: {error_msg_detail}", exc_info=True)
        
        if temp_manual_tle_path and os.path.exists(temp_manual_tle_path):
            try:
                os.remove(temp_manual_tle_path)
                logger.info(f"Cleaned up temp manual TLE due to error: {temp_manual_tle_path}")
            except Exception as e_clean_err:
                logger.warning(f"Could not clean up temp TLE {temp_manual_tle_path} during error handling: {e_clean_err}")
        
        current_pos_str_val = current_pos_str_val if 'current_pos_str_val' in locals() else ""
        passes_str_val = passes_str_val if 'passes_str_val' in locals() else ""
        folium_map_html_output = folium_map_html_output if 'folium_map_html_output' in locals() else EMPTY_MAP_CONTENT
        slider_update = slider_update if 'slider_update' in locals() else gr.update(visible=False)
        app_state_to_return = app_state if 'app_state' in locals() and isinstance(app_state, dict) else {}
        return (current_pos_str_val, passes_str_val, status_log, folium_map_html_output, 
                slider_update, app_state_to_return)

    return (current_pos_str_val, passes_str_val, status_log, folium_map_html_output, 
            slider_update, app_state)

def handle_gif_generation(current_app_state, gif_duration_mins, gif_fps, gif_sim_step_secs, gif_map_style, min_elevation_gif):
    """
    Generate a footprint animation GIF using Cartopy-based static rendering.

    Parameters
    ----------
    current_app_state : dict
        Global state containing satellite, observer, and reference time.
    gif_duration_mins : int | float
        Total animation duration in minutes.
    gif_fps : int
        Output GIF playback frames per second.
    gif_sim_step_secs : int
        Simulation time step in seconds per frame.
    gif_map_style : str
        Static map background style name.
    min_elevation_gif : float
        Minimum elevation (deg) for footprint rendering.

    Returns
    -------
    tuple[str | None, str]
        (gif_output_path, gif_status_log). Path may be None on failure.
    """
    gif_status_log = ""
    gif_output_path = None

    try:
        logger.info(f"handle_gif_generation: Request to generate GIF. Duration: {gif_duration_mins}m, FPS: {gif_fps}, SimStep: {gif_sim_step_secs}s, Style: {gif_map_style}, MinElev: {min_elevation_gif}°")

        satellite_obj = current_app_state.get("current_satellite_skyfield_obj")
        observer_obj = current_app_state.get("current_observer_topos_obj")
        track_ref_time_utc_from_state = current_app_state.get("track_reference_time_utc")

        if not (satellite_obj and observer_obj and track_ref_time_utc_from_state is not None):
            gif_status_log = "錯誤：無法產生 GIF。缺少必要的衛星、觀測者或時間資訊。請先執行主分析。<br>"
            logger.warning("handle_gif_generation: Missing necessary data in app_state for GIF generation.")
            return None, gif_status_log

        # Normalize reference time to a Python datetime (UTC)
        animation_start_time_utc = track_ref_time_utc_from_state
        if isinstance(track_ref_time_utc_from_state, str):
            try:
                animation_start_time_utc = datetime.fromisoformat(track_ref_time_utc_from_state.replace("Z", "+00:00"))
                if animation_start_time_utc.tzinfo is None:
                    animation_start_time_utc = animation_start_time_utc.replace(tzinfo=timezone.utc)
            except ValueError:
                logger.error(f"Could not parse track_reference_time_utc from string for GIF: {track_ref_time_utc_from_state}")
                gif_status_log = "錯誤: 無法解析軌跡參考時間以產生 GIF。<br>"
                return None, gif_status_log
        elif hasattr(track_ref_time_utc_from_state, 'utc_datetime'):
            try:
                animation_start_time_utc = track_ref_time_utc_from_state.utc_datetime()
            except Exception as e:
                logger.error(f"Could not convert Skyfield Time to datetime for GIF: {e}")
                gif_status_log = "錯誤: 無法轉換 Skyfield 時間物件以產生 GIF。<br>"
                return None, gif_status_log
        elif not isinstance(animation_start_time_utc, datetime):
            logger.error(f"track_reference_time_utc for GIF is not a datetime object, string, or Skyfield Time: {type(animation_start_time_utc)}")
            gif_status_log = "錯誤: 軌跡參考時間格式不正確以產生 GIF。<br>"
            return None, gif_status_log

        gif_status_log += f"正在準備產生 GIF 動畫...<br>動畫開始時間: {animation_start_time_utc.strftime('%Y-%m-%d %H:%M:%S')} UTC<br>"
        
        # Unique filename for output (into data/ directory)
        timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_gif_filename = f"footprint_anim_{timestamp_str}.gif"
        gif_save_path = os.path.join("data", output_gif_filename)

        gif_output_path = visualizer.generate_footprint_animation_gif(
            satellite=satellite_obj,
            observer_location=observer_obj,
            animation_start_time_utc=animation_start_time_utc,
            animation_duration_minutes=int(gif_duration_mins),
            simulation_time_step_seconds=int(gif_sim_step_secs),
            gif_playback_fps=int(gif_fps),
            output_gif_path=gif_save_path,
            min_elevation_for_footprint_deg=float(min_elevation_gif),
            map_style=gif_map_style,
            base_track_duration_hours=float(current_app_state.get("track_duration_hours", 3.0))
        )

        if gif_output_path:
            gif_status_log += f"GIF 動畫已成功產生並儲存至: {os.path.basename(gif_output_path)}<br>"
            logger.info(f"handle_gif_generation: GIF generated successfully at {gif_output_path}")
        else:
            gif_status_log += "GIF 動畫產生失敗。請查看應用程式日誌 (data/app.log) 以獲取詳細資訊。<br>"
            logger.error("handle_gif_generation: GIF generation failed (visualizer returned None).")

    except Exception as e:
        error_msg = f"產生 GIF 動畫時發生未預期錯誤: {e}"
        gif_status_log += f"嚴重錯誤：{error_msg}<br>"
        logger.error(f"handle_gif_generation: Unexpected error: {error_msg}", exc_info=True)
        gif_output_path = None

    return gif_output_path, gif_status_log

# =========================
# Gradio UI
# =========================
if __name__ == "__main__":
    with gr.Blocks(theme=gr.themes.Soft()) as demo:
        app_state = gr.State({})
        gr.Markdown("""# 簡易衛星地面站""", elem_id="main-title")

        with gr.Row():
            with gr.Column(scale=1):
                gr.Markdown("### 輸入參數")
                satellite_name_input = gr.Textbox(label="衛星名稱 (例如：ISS (ZARYA))", value=DEFAULT_SATELLITE_NAME)
                observer_lat_input = gr.Number(label="觀測者緯度 (°)", value=DEFAULT_OBSERVER_LAT)
                observer_lon_input = gr.Number(label="觀測者經度 (°)", value=DEFAULT_OBSERVER_LON)
                
                # TLE category dropdown (source list only)
                base_tle_sources = get_tle_sources()
                tle_source_dropdown = gr.Dropdown(label="TLE 衛星分類", choices=base_tle_sources, value=base_tle_sources[0] if base_tle_sources else "Active Satellites")
                
                with gr.Accordion("進階 TLE 選項", open=False):
                    tle_mode = gr.Radio(
                        label="TLE 取得方式", 
                        choices=["自動下載最新 TLE", "上傳 TLE 檔案", "手動輸入 TLE"], 
                        value="自動下載最新 TLE"
                    )
                    # Only visible when "Upload" mode is selected
                    tle_upload_button = gr.File(label="上傳 TLE 檔案 (.tle)", file_types=[".tle"], visible=False)
                    # Only visible when "Manual" mode is selected
                    manual_tle_input_textbox = gr.Textbox(label="手動輸入 TLE (兩行格式)", lines=3, placeholder="LINE 1\\\\nLINE 2", visible=False)

                min_elevation_input = gr.Slider(label="Footprint 最小仰角 (°)", minimum=0, maximum=90, value=DEFAULT_MIN_ELEVATION_FOOTPRINT, step=1)
                map_style_dropdown_folium = gr.Dropdown(label="地圖底圖樣式 (互動圖)", choices=FOLIUM_MAP_TILES_CHOICES, value=DEFAULT_FOLIUM_MAP_TILES)
                
                run_button = gr.Button("執行分析", variant="primary")

                gr.Markdown("### 狀態與日誌")
                status_log_output = gr.HTML(label="狀態日誌", value="請點擊執行分析。")
                
            with gr.Column(scale=2):
                gr.Markdown("### 目前衛星位置")
                current_position_output = gr.Markdown()

                gr.Markdown("### 衛星過境預報")
                pass_predictions_output = gr.HTML()
                
                gr.Markdown("### 互動式軌跡地圖")
                folium_map_output = gr.HTML(label="軌跡地圖")
                
                # Hidden by default; becomes visible after first successful run
                time_slider = gr.Slider(label="時間軸 (相對軌跡開始時間, 分鐘)", minimum=0, maximum=180, step=1, visible=False, interactive=True)

        # --- UI event handlers
        
        def update_tle_options_visibility(tle_mode):
            """
            Toggle visibility of the TLE file upload and manual input widgets.

            Parameters
            ----------
            tle_mode : str
                Selected mode for TLE acquisition.

            Returns
            -------
            tuple
                (upload_visible_update, manual_visible_update)
            """
            upload_visible = tle_mode == "上傳 TLE 檔案"
            manual_visible = tle_mode == "手動輸入 TLE"
            return gr.update(visible=upload_visible), gr.update(visible=manual_visible)
        
        tle_mode.change(
            fn=update_tle_options_visibility,
            inputs=[tle_mode],
            outputs=[tle_upload_button, manual_tle_input_textbox]
        )
        
        run_button.click(
            fn=ground_station_interface,
            inputs=[satellite_name_input, observer_lat_input, observer_lon_input, tle_source_dropdown, tle_mode, tle_upload_button, manual_tle_input_textbox, map_style_dropdown_folium, min_elevation_input, app_state],
            outputs=[
                current_position_output, 
                pass_predictions_output, 
                status_log_output, 
                folium_map_output,
                time_slider,
                app_state
            ]
        )
        
        time_slider.change(
            fn=update_map_for_time_slider,
            inputs=[time_slider, app_state],
            outputs=[folium_map_output, status_log_output, app_state]
        )
        
        # GIF generator panel
        with gr.Accordion("GIF 動畫產生器", open=False):
            with gr.Row():
                with gr.Column(scale=1):
                    gif_duration_input = gr.Number(label="動畫時長 (分鐘)", value=DEFAULT_GIF_DURATION_MINUTES, minimum=1)
                    gif_fps_input = gr.Number(label="GIF 播放 FPS", value=DEFAULT_GIF_FPS, minimum=1, step=1)
                    gif_sim_step_input = gr.Number(label="模擬時間間隔 (秒/幀)", value=30, minimum=1, step=1)
                    gif_map_style_dropdown = gr.Dropdown(label="地圖樣式 (GIF)", choices=STATIC_MAP_STYLE_CHOICES, value=DEFAULT_STATIC_MAP_STYLE)
                with gr.Column(scale=1):
                    generate_gif_button = gr.Button("產生軌跡動畫 (GIF)")
                    cancel_gif_button = gr.Button("取消產生")
            with gr.Row():        
                with gr.Column(scale=2):
                    gif_output_image = gr.Image(label="GIF 動畫", type="filepath", interactive=False, height=400)
                    gif_status_log_output = gr.Markdown()
            
            generate_gif_button.click(
                fn=handle_gif_generation,
                inputs=[
                    app_state,
                    gif_duration_input,
                    gif_fps_input,
                    gif_sim_step_input,
                    gif_map_style_dropdown,
                    min_elevation_input
                ],
                outputs=[gif_output_image, gif_status_log_output]
            )

    demo.queue()
    demo.launch(share=False)
