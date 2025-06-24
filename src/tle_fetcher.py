import requests
import os

# Custom error type for TLE fetching failures
class TLEFetchError(Exception):
    """Custom exception for TLE fetching errors."""
    pass
# Default data directory for storing downloaded TLE files
DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
# Default TLE file path used when the user does not upload/select another source
DEFAULT_TLE_FILE = os.path.join(DATA_DIR, "active_satellites.tle")
# Multiple available TLE sources on Celestrak (grouped endpoints)
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
    """
    Ensure the local data directory exists.

    Creates DATA_DIR if it does not already exist.
    """
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)

def get_tle_sources():
    """
    Return the list of available TLE source names.

    Returns
    -------
    list[str]
        A list of keys corresponding to available sources in TLE_SOURCES.
    """
    return list(TLE_SOURCES.keys())

def download_tle_file(tle_source_name="Active Satellites", tle_file_path=None):
    """
    Download a TLE file from a specified Celestrak source.

    Parameters
    ----------
    tle_source_name : str, optional
        The display name of the TLE source as defined in TLE_SOURCES.
        Defaults to "Active Satellites".
    tle_file_path : str | None, optional
        Local filesystem path to save the downloaded TLE file. If None,
        a default path under DATA_DIR will be used: data/{safe_source_name}.tle

    Returns
    -------
    str
        The local path to the downloaded TLE file on success.

    Raises
    ------
    TLEFetchError
        If the HTTP request fails or the file cannot be saved.
    """
    if tle_source_name not in TLE_SOURCES:
        print(f"錯誤：未知的 TLE 來源 '{tle_source_name}'")
        return None

    url = TLE_SOURCES[tle_source_name]
    
    # Normalize the source name for a safe filename (lowercase, non-alnum -> underscore)
    safe_source_name = "".join(c if c.isalnum() else "_" for c in tle_source_name).lower()
    
    if tle_file_path is None:
        ensure_data_dir_exists()
        tle_file_path = os.path.join(DATA_DIR, f"{safe_source_name}.tle")

    try:
        print(f"正在從 {url} 下載 TLE 檔案到 {tle_file_path}...")
        # Set a timeout to avoid indefinite waiting on network calls
        response = requests.get(url, timeout=10)
        # Raise HTTPError for 4xx/5xx responses
        response.raise_for_status()

        # Save content to the chosen path in binary mode
        with open(tle_file_path, "wb") as f:
            f.write(response.content)
        
        print(f"TLE 檔案已成功下載到: {tle_file_path}")
        return tle_file_path
    
    except requests.exceptions.RequestException as e:
        # Network / HTTP errors
        error_msg = f"下載 TLE 檔案失敗: {e}"
        print(f"錯誤：{error_msg}")
        raise TLEFetchError(error_msg) from e
    
    except IOError as e:
        # Local filesystem write errors
        error_msg = f"儲存 TLE 檔案失敗: {e}"
        print(f"錯誤：{error_msg}")
        raise TLEFetchError(error_msg) from e

if __name__ == "__main__":
    # Simple smoke test for downloading the default source
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
