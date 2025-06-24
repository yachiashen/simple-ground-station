from skyfield.api import load, wgs84, Topos, N, E, S, W
from skyfield.sgp4lib import EarthSatellite
from skyfield.framelib import itrs
import os
import datetime

# TLE 格式錯誤訊息常數
TLE_FORMAT_ERROR_MESSAGE = "TLE file format error or expired"

# 預設 TLE 檔案路徑 (假設 tle_fetcher.py 會將其下載到此處)
DEFAULT_TLE_FILE = os.path.join(os.path.dirname(__file__), '..', 'data', 'active.tle')

# Skyfield 時間尺度
ts = load.timescale()

def load_satellites_from_tle(tle_file=DEFAULT_TLE_FILE):
    """
    從 TLE 檔案載入衛星資料。

    Args:
        tle_file (str): TLE 檔案的路徑。

    Returns:
        list: EarthSatellite 物件的列表，如果檔案不存在或為空則返回空列表。
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
    從衛星列表中按名稱搜尋衛星。

    Args:
        satellites (list): EarthSatellite 物件的列表。
        name_query (str): 要搜尋的衛星名稱 (部分或完整)。

    Returns:
        EarthSatellite: 找到的衛星物件，若未找到則返回 None。
    """
    name_query = name_query.upper()
    for sat in satellites:
        if name_query in sat.name.upper():
            return sat
    print(f"未找到名稱包含 '{name_query}' 的衛星。")
    return None

def get_observer_location(latitude_deg, longitude_deg, elevation_m=0):
    """
    建立觀測者位置的 Topos 物件。

    Args:
        latitude_deg (float): 觀測者緯度 (度)。
        longitude_deg (float): 觀測者經度 (度)。
        elevation_m (float): 觀測者海拔高度 (米)。

    Returns:
        Topos: 代表觀測者位置的 Skyfield Topos 物件。
    """
    return wgs84.latlon(latitude_deg, longitude_deg, elevation_m)

# --- 核心計算函數 ---

def get_satellite_current_position(satellite: EarthSatellite, observer_location: Topos):
    """
    計算衛星相對於觀測者的目前位置和星下點。

    Args:
        satellite (EarthSatellite): 要計算的衛星物件。
        observer_location (Topos): 觀測者位置的 Skyfield Topos 物件。

    Returns:
        dict: 包含衛星位置資訊的字典，格式如下：
              {
                  "azimuth_deg": float,      # 方位角 (度)
                  "elevation_deg": float,    # 仰角 (度)
                  "distance_km": float,      # 距離 (公里)
                  "subpoint_lat_deg": float, # 星下點緯度 (度)
                  "subpoint_lon_deg": float, # 星下點經度 (度)
                  "altitude_km": float       # 衛星高度 (公里)
              }
              如果計算失敗則返回 None。
    """
    try:
        now = ts.now()  # 取得目前時間
        difference = satellite - observer_location
        topocentric = difference.at(now)  # 計算衛星相對於觀測站的位置

        alt, az, distance = topocentric.altaz() # 取得仰角、方位角、距離

        # 計算衛星的星下點 (subpoint)
        geocentric = satellite.at(now)
        subpoint = wgs84.subpoint(geocentric)

        return {
            "azimuth_deg": az.degrees,
            "elevation_deg": alt.degrees,
            "distance_km": distance.km,
            "subpoint_lat_deg": subpoint.latitude.degrees,
            "subpoint_lon_deg": subpoint.longitude.degrees,
            "altitude_km": subpoint.elevation.km  # 這裡的 elevation 是指衛星離地表的高度
        }
    except Exception as e:
        print(f"計算衛星目前位置時發生錯誤: {e}")
        return None

def predict_satellite_passes(satellite: EarthSatellite, observer_location: Topos,
                             start_time_utc: datetime.datetime,
                             duration_hours: float,
                             min_elevation_deg=10.0):
    """
    預測衛星在指定時間段內經過觀測點的情況。

    Args:
        satellite (EarthSatellite): 要預測的衛星物件。
        observer_location (Topos): 觀測者位置的 Skyfield Topos 物件。
        start_time_utc (datetime.datetime): 預測開始時間 (UTC)。
        duration_hours (float): 從開始時間算起，要預測的時長（小時）。
        min_elevation_deg (float): 衛星可見的最小仰角 (度)。

    Returns:
        list: 包含過境事件資訊的字典列表。每個字典代表一次過境，包含：
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
              如果計算失敗或沒有過境事件則返回空列表。
    """
    try:
        # 計算結束時間
        end_time_utc = start_time_utc + datetime.timedelta(hours=duration_hours)

        # 將 datetime 物件轉換為 Skyfield Time 物件
        t0 = ts.utc(start_time_utc.year, start_time_utc.month, start_time_utc.day,
                    start_time_utc.hour, start_time_utc.minute, start_time_utc.second)
        t1 = ts.utc(end_time_utc.year, end_time_utc.month, end_time_utc.day,
                    end_time_utc.hour, end_time_utc.minute, end_time_utc.second)

        times, events = satellite.find_events(observer_location, t0, t1, altitude_degrees=min_elevation_deg)

        passes = []
        # 事件類型: 0 = 上升 (rise), 1 = 到達最高點 (culminate), 2 = 下降 (set)
        # 我們需要將這些事件組合成完整的過境
        
        # 使用一個迴圈來迭代事件，並建立完整的過境資訊
        i = 0
        while i < len(events):
            # 尋找上升事件
            if events[i] == 0:
                rise_t = times[i]
                
                # 尋找接下來的中天和下降事件
                cul_t = None
                set_t = None
                
                # 從上升事件之後開始尋找
                j = i + 1
                while j < len(events):
                    if events[j] == 1 and cul_t is None: # 找到中天
                        cul_t = times[j]
                    elif events[j] == 2: # 找到下降
                        set_t = times[j]
                        break # 找到完整的過境，跳出內部迴圈
                    j += 1
                
                if cul_t is not None and set_t is not None:
                    # 計算上升時的方位角
                    _, rise_az, _ = (satellite - observer_location).at(rise_t).altaz()
                    # 計算中天時的方位角和仰角
                    cul_alt, cul_az, _ = (satellite - observer_location).at(cul_t).altaz()
                    # 計算下降時的方位角
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
                    i = j # 移動外部迴圈的索引到下降事件之後
                else:
                    # 如果沒有找到完整的中天和下降，則跳過此上升事件
                    i += 1
            else:
                i += 1 # 如果不是上升事件，繼續尋找下一個

        return passes

    except Exception as e:
        print(f"預測衛星過境時發生錯誤: {e}")
        import traceback
        traceback.print_exc()
        return []

if __name__ == "__main__":
    # 測試載入 TLE 和選取衛星
    sats = load_satellites_from_tle() 
    if sats:
        # 嘗試尋找國際太空站 (ISS)
        iss_name_parts = ["ISS (ZARYA)", "ISS", "ZARYA"]
        iss = None
        for name_part in iss_name_parts:
            iss = get_satellite_by_name(sats, name_part)
            if iss:
                print(f"成功找到衛星: {iss.name} ({iss.model.satnum})")
                break
        if not iss:
            print("在預設 TLE 檔案中未找到 ISS。請確認 active.tle 包含 ISS 資料，或嘗試其他衛星名稱。")

        # 測試觀測點 (例如：台北)
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
