import requests
import os

# 自訂錯誤類別
class TLEFetchError(Exception):
    """Custom exception for TLE fetching errors."""
    pass

# 預設的 TLE 檔案儲存路徑
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
# 預設 TLE 檔案路徑，如果使用者未上傳檔案且未選擇其他來源，則使用此檔案
DEFAULT_TLE_FILE = os.path.join(DATA_DIR, "active_satellites.tle")

# 定義多個 TLE 來源
TLE_SOURCES = {
    "Active Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle",
    "Weather Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=weather&FORMAT=tle",
    "NOAA Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=noaa&FORMAT=tle",
    "GOES Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=goes&FORMAT=tle",
    "Earth Resources Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=resource&FORMAT=tle",
    "Search & Rescue (SARSAT) Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=sarsat&FORMAT=tle",
    "Disaster Monitoring Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=dmc&FORMAT=tle",
    "Tracking and Data Relay System (TDRSS) Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=tdrss&FORMAT=tle",
    "GPS Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=gps-ops&FORMAT=tle",
    "GLONASS Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=glo-ops&FORMAT=tle",
    "Galileo Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=galileo&FORMAT=tle",
    "Beidou Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=beidou&FORMAT=tle",
    "Satellite-Based Augmentation System (SBAS) Satellites": "https://celestrak.org/NORAD/elements/gp.php?GROUP=sbas&FORMAT=tle",
    "ISS": "https://celestrak.org/NORAD/elements/gp.php?NAME=ISS&FORMAT=tle",
    "Starlink": "https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle",
}

def ensure_data_dir_exists():
    """確保 data 資料夾存在"""
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)

def get_tle_sources():
    """返回可用的 TLE 來源名稱列表"""
    return list(TLE_SOURCES.keys())

def download_tle_file(tle_source_name="Active Satellites", tle_file_path=None):
    """
    從指定的 Celestrak 來源下載 TLE 檔案。

    Args:
        tle_source_name (str): TLE_SOURCES 字典中的來源名稱。
        tle_file_path (str, optional): 下載 TLE 檔案的本機路徑。
                                     如果為 None，則使用預設路徑 data/{source_name}.tle。

    Returns:
        str: 下載的 TLE 檔案路徑，如果下載失敗則返回 None。
    """
    if tle_source_name not in TLE_SOURCES:
        print(f"錯誤：未知的 TLE 來源 '{tle_source_name}'")
        return None

    url = TLE_SOURCES[tle_source_name]
    
    # 標準化來源名稱以用於檔案名，例如 "Active Satellites" -> "active_satellites.tle"
    # 並移除特殊字元，轉換為小寫
    safe_source_name = "".join(c if c.isalnum() else "_" for c in tle_source_name).lower()
    
    if tle_file_path is None:
        ensure_data_dir_exists()
        tle_file_path = os.path.join(DATA_DIR, f"{safe_source_name}.tle")

    try:
        print(f"正在從 {url} 下載 TLE 檔案到 {tle_file_path}...")
        response = requests.get(url, timeout=10) # 設定超時以避免無限等待
        response.raise_for_status()  # 如果請求失敗 (狀態碼 4xx 或 5xx)，則拋出 HTTPError

        with open(tle_file_path, "wb") as f:
            f.write(response.content)
        print(f"TLE 檔案已成功下載到: {tle_file_path}")
        return tle_file_path
    except requests.exceptions.RequestException as e:
        error_msg = f"下載 TLE 檔案失敗: {e}"
        print(f"錯誤：{error_msg}")
        raise TLEFetchError(error_msg) from e
    except IOError as e:
        error_msg = f"儲存 TLE 檔案失敗: {e}"
        print(f"錯誤：{error_msg}")
        raise TLEFetchError(error_msg) from e

if __name__ == "__main__":
    try:
        downloaded_file_path = download_tle_file()
        if downloaded_file_path:
            print(f"測試下載完成，檔案位於: {downloaded_file_path}")
        else:
            print("測試下載失敗 (download_tle_file 返回 None，但不應發生)。")
    except TLEFetchError as e:
        print(f"測試下載過程中發生 TLEFetchError: {e}")
    except Exception as e:
        print(f"測試下載過程中發生未預期錯誤: {e}")
