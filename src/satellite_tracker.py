from skyfield.api import load, wgs84, Topos, N, E, S, W
from skyfield.sgp4lib import EarthSatellite
from skyfield.framelib import itrs
import os
import datetime

# Constant for TLE format error message
TLE_FORMAT_ERROR_MESSAGE = "TLE file format error or expired"
# Default TLE file path
DEFAULT_TLE_FILE = os.path.join(os.path.dirname(__file__), '..', 'data', 'active.tle')

# Skyfield timescale (global)
ts = load.timescale()

def load_satellites_from_tle(tle_file=DEFAULT_TLE_FILE):
    """
    Load satellites from a TLE file.

    Parameters
    ----------
    tle_file : str
        Path to the TLE file.

    Returns
    -------
    list[EarthSatellite]
        A list of EarthSatellite objects. Returns an empty list if file does not
        exist or if loading fails.
    """
    if not os.path.exists(tle_file):
        print(f"錯誤: TLE 檔案 {tle_file} 不存在。")
        return []
    try:
        satellites = load.tle_file(tle_file)
        print(f"已從 {tle_file} 成功載入 {len(satellites)} 顆衛星。")
        return satellites
    except Exception as e:
        print(f"載入 TLE 檔案時發生錯誤: {e}")
        return []

def get_satellite_by_name(satellites, name_query):
    """
    Find a satellite by (partial) name within a list.

    Parameters
    ----------
    satellites : list[EarthSatellite]
        List of EarthSatellite objects to search.
    name_query : str
        Full or partial satellite name (case-insensitive).

    Returns
    -------
    EarthSatellite | None
        Matching satellite object if found, otherwise None.
    """
    name_query = name_query.upper()
    for sat in satellites:
        if name_query in sat.name.upper():
            return sat
    print(f"未找到名稱包含 '{name_query}' 的衛星。")
    return None

def get_observer_location(latitude_deg, longitude_deg, elevation_m=0):
    """
    Build a Skyfield Topos object for the observer location.

    Parameters
    ----------
    latitude_deg : float
        Observer latitude in degrees.
    longitude_deg : float
        Observer longitude in degrees.
    elevation_m : float, optional
        Observer elevation in meters (default is 0).

    Returns
    -------
    Topos
        Skyfield Topos object representing the observer location.
    """
    return wgs84.latlon(latitude_deg, longitude_deg, elevation_m)

# --- Core computation functions ---

def get_satellite_current_position(satellite: EarthSatellite, observer_location: Topos):
    """
    Compute the satellite's current topocentric position (relative to observer)
    and subpoint on Earth.

    Parameters
    ----------
    satellite : EarthSatellite
        Satellite to compute.
    observer_location : Topos
        Skyfield Topos for the observer.

    Returns
    -------
    dict | None
        Dictionary containing:
        {
            "azimuth_deg": float,       # azimuth (deg)
            "elevation_deg": float,     # elevation (deg)
            "distance_km": float,       # range (km)
            "subpoint_lat_deg": float,  # subpoint latitude (deg)
            "subpoint_lon_deg": float,  # subpoint longitude (deg)
            "altitude_km": float        # satellite altitude above ground (km)
        }
        Returns None on failure.
    """
    try:
        # Current time
        now = ts.now()
        
        # Topocentric vector from observer to satellite
        difference = satellite - observer_location
        topocentric = difference.at(now)

        # Extract elevation, azimuth, and range
        alt, az, distance = topocentric.altaz()

        # Satellite subpoint (projection on Earth surface)
        geocentric = satellite.at(now)
        subpoint = wgs84.subpoint(geocentric)

        return {
            "azimuth_deg": az.degrees,
            "elevation_deg": alt.degrees,
            "distance_km": distance.km,
            "subpoint_lat_deg": subpoint.latitude.degrees,
            "subpoint_lon_deg": subpoint.longitude.degrees,
            "altitude_km": subpoint.elevation.km
        }
    except Exception as e:
        print(f"計算衛星目前位置時發生錯誤: {e}")
        return None

def predict_satellite_passes(satellite: EarthSatellite, observer_location: Topos,
                             start_time_utc: datetime.datetime,
                             duration_hours: float,
                             min_elevation_deg=10.0):
    """
    Predict satellite passes over the observer within a time window.

    This stitches together events from Skyfield's `find_events`:
    0 = rise, 1 = culminate (highest elevation), 2 = set,
    and returns a list of complete passes (rise -> culminate -> set)
    whose elevation exceeds the given threshold.

    Parameters
    ----------
    satellite : EarthSatellite
        Satellite to predict.
    observer_location : Topos
        Skyfield Topos for the observer.
    start_time_utc : datetime.datetime
        Start time (UTC).
    duration_hours : float
        Duration from start time in hours.
    min_elevation_deg : float, optional
        Minimum elevation (deg) for visibility (default 10.0).

    Returns
    -------
    list[dict]
        Each pass dictionary contains:
        {
            "rise_time_utc": datetime.datetime,
            "rise_azimuth_deg": float,
            "culmination_time_utc": datetime.datetime,
            "culmination_azimuth_deg": float,
            "culmination_elevation_deg": float,
            "set_time_utc": datetime.datetime,
            "set_azimuth_deg": float,
            "duration_minutes": float
        }
        Returns an empty list if no passes or on failure.
    """
    try:
        # Compute end time
        end_time_utc = start_time_utc + datetime.timedelta(hours=duration_hours)

        # Convert datetimes to Skyfield Time objects (UTC)
        t0 = ts.utc(start_time_utc.year, start_time_utc.month, start_time_utc.day,
                    start_time_utc.hour, start_time_utc.minute, start_time_utc.second)
        t1 = ts.utc(end_time_utc.year, end_time_utc.month, end_time_utc.day,
                    end_time_utc.hour, end_time_utc.minute, end_time_utc.second)

        # Find rise/culmination/set events above elevation threshold
        times, events = satellite.find_events(observer_location, t0, t1, altitude_degrees=min_elevation_deg)

        passes = []
        # Event codes: 0 = rise, 1 = culminate, 2 = set
        # Stitch sequential events into complete passes
        i = 0
        while i < len(events):
            if events[i] == 0:        # rise
                rise_t = times[i]                
                cul_t = None
                set_t = None
                
                # Search for subsequent culmination and set
                j = i + 1
                while j < len(events):
                    if events[j] == 1 and cul_t is None:
                        cul_t = times[j]
                    elif events[j] == 2:
                        set_t = times[j]
                        break
                    j += 1
                
                if cul_t is not None and set_t is not None:
                    # Azimuth at rise
                    _, rise_az, _ = (satellite - observer_location).at(rise_t).altaz()
                    # Azimuth/elevation at culmination
                    cul_alt, cul_az, _ = (satellite - observer_location).at(cul_t).altaz()
                    # Azimuth at set
                    _, set_az, _ = (satellite - observer_location).at(set_t).altaz()

                    duration = (set_t.utc_datetime() - rise_t.utc_datetime()).total_seconds() / 60.0

                    passes.append({
                        "rise_time_utc": rise_t.utc_datetime().replace(tzinfo=datetime.timezone.utc),
                        "rise_azimuth_deg": rise_az.degrees,
                        "culmination_time_utc": cul_t.utc_datetime().replace(tzinfo=datetime.timezone.utc),
                        "culmination_azimuth_deg": cul_az.degrees,
                        "culmination_elevation_deg": cul_alt.degrees,
                        "set_time_utc": set_t.utc_datetime().replace(tzinfo=datetime.timezone.utc),
                        "set_azimuth_deg": set_az.degrees,
                        "duration_minutes": duration
                    })
                    # Move i to the set event index to continue
                    i = j
                else:
                    # Incomplete sequence; skip this rise
                    i += 1
            else:
                # Not a rise; advance
                i += 1

        return passes

    except Exception as e:
        print(f"預測衛星過境時發生錯誤: {e}")
        import traceback
        traceback.print_exc()
        return []

if __name__ == "__main__":
    # Test loading TLE and selecting a satellite
    sats = load_satellites_from_tle() 
    if sats:
        # Try finding the ISS (multiple name variants just in case)
        iss_name_parts = ["ISS (ZARYA)", "ISS", "ZARYA"]
        iss = None
        for name_part in iss_name_parts:
            iss = get_satellite_by_name(sats, name_part)
            if iss:
                print(f"成功找到衛星: {iss.name} ({iss.model.satnum})")
                break
        if not iss:
            print("在預設 TLE 檔案中未找到 ISS。請確認 active.tle 包含 ISS 資料，或嘗試其他衛星名稱。")

        # Example observer: Taipei
        taipei_lat = 25.0330
        taipei_lon = 121.5654
        observer = get_observer_location(taipei_lat, taipei_lon)
        print(f"觀測點設定: 緯度 {taipei_lat}, 經度 {taipei_lon}")

        if iss and observer:
            print("\n--- 測試目前位置計算 ---")
            current_pos = get_satellite_current_position(iss, observer)
            if current_pos:
                print(f"衛星 {iss.name} 目前位置:")
                print(f"  星下點: 緯度 {current_pos['subpoint_lat_deg']:.2f}°, 經度 {current_pos['subpoint_lon_deg']:.2f}°")
                print(f"  高度: {current_pos['altitude_km']:.2f} km")
                print(f"  相對於觀測點:")
                print(f"    方位角: {current_pos['azimuth_deg']:.2f}°")
                print(f"    仰角: {current_pos['elevation_deg']:.2f}°")
                print(f"    距離: {current_pos['distance_km']:.2f} km")
            else:
                print(f"無法計算衛星 {iss.name} 的目前位置。")

            print("\n--- 測試未來24小時過境預測 (最低仰角10°) ---")
            now_utc = datetime.datetime.now(datetime.timezone.utc)
            start_predict_time = now_utc
            duration_hours = 24.0
            
            print(f"預測時間範圍 (UTC): {start_predict_time.strftime('%Y-%m-%d %H:%M:%S %Z')} 持續 {duration_hours} 小時")

            passes = predict_satellite_passes(iss, observer, start_predict_time, duration_hours, min_elevation_deg=10.0)

            if passes:
                print(f"找到 {len(passes)} 次 {iss.name} 的過境事件:")
                for p_idx, p_info in enumerate(passes):
                    print(f"\n  過境事件 #{p_idx + 1}")
                    print(f"    上升時間: {p_info['rise_time_utc'].strftime('%Y-%m-%d %H:%M:%S %Z')} (方位角: {p_info['rise_azimuth_deg']:.1f}°)")
                    print(f"    最高點時間: {p_info['culmination_time_utc'].strftime('%Y-%m-%d %H:%M:%S %Z')} (方位角: {p_info['culmination_azimuth_deg']:.1f}°, 仰角: {p_info['culmination_elevation_deg']:.1f}°)")
                    print(f"    下降時間: {p_info['set_time_utc'].strftime('%Y-%m-%d %H:%M:%S %Z')} (方位角: {p_info['set_azimuth_deg']:.1f}°)")
                    print(f"    持續時間: {p_info['duration_minutes']:.1f} 分鐘")
            else:
                print(f"在接下來的24小時內，未找到 {iss.name} 高於10°仰角的過境事件。")
        
        elif not iss:
            print("由於未找到衛星，無法進行下一步計算。")
        elif not observer:
            print("觀測點設定失敗，無法進行下一步計算。")

    else:
        print("未能載入衛星資料，請先執行 tle_fetcher.py 以下載 TLE 檔案。")
