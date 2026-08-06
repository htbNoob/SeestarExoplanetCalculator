"""
Observer site and tonight's dark-window helpers, shared by the "tonight's
best transits" and "tonight's best eclipsing binaries" pages.
"""
from dataclasses import dataclass
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

import numpy as np
from astropy.coordinates import AltAz, EarthLocation, get_sun
from astropy.time import Time
import astropy.units as u


@dataclass
class Site:
    name: str
    lat: float   # degrees, +N
    lon: float   # degrees, +E
    elev_m: float = 0.0
    tz: str = "UTC"   # IANA timezone name (e.g. "Europe/Berlin"), for display only


@dataclass
class DarkWindow:
    site: Site
    obs_loc: EarthLocation
    dark_start: Time
    dark_end: Time
    ds_jd: float
    de_jd: float
    placeholder: bool = False   # True if no true dark window was found in the scan


def compute_dark_window(site: Site, horizon_deg: float = -12.0, scan_hours: float = 36.0) -> DarkWindow:
    """Find tonight's dark window (sun below `horizon_deg`) for a site.

    Mirrors the notebook's nautical-twilight (-12 deg) scan used ahead of both
    the transit and eclipsing-binary searches.
    """
    obs_loc = EarthLocation(lat=site.lat * u.deg, lon=site.lon * u.deg, height=site.elev_m * u.m)
    now = Time.now()

    t_arr = now + np.arange(0, scan_hours, 0.1) * u.hour
    af_scan = AltAz(obstime=t_arr, location=obs_loc)
    sun_alt = get_sun(t_arr).transform_to(af_scan).alt.deg
    dark = sun_alt < horizon_deg

    if dark.sum() > 2:
        idx0 = int(np.argmax(dark))
        tail = dark[idx0:]
        idx1 = idx0 + (len(tail) - 1 if tail.all() else int(np.argmax(~tail)) - 1)
        dark_start = t_arr[idx0]
        dark_end = t_arr[idx1]
        placeholder = False
    else:
        dark_start, dark_end = now, now + 8 * u.hour
        placeholder = True

    return DarkWindow(
        site=site,
        obs_loc=obs_loc,
        dark_start=dark_start,
        dark_end=dark_end,
        ds_jd=dark_start.tdb.jd,
        de_jd=dark_end.tdb.jd,
        placeholder=placeholder,
    )


def min_altitude_deg(coord, t0_jd, t1_jd, obs_loc, n=10) -> float:
    """Minimum altitude (deg) of `coord` between two JD times, as seen from `obs_loc`."""
    tt = Time(np.linspace(t0_jd, t1_jd, n), format="jd", scale="tdb")
    return float(coord.transform_to(AltAz(obstime=tt, location=obs_loc)).alt.deg.min())


def is_valid_tz(tz_name: str) -> bool:
    """Whether `tz_name` is a recognized IANA timezone name."""
    try:
        ZoneInfo(tz_name)
        return True
    except (ZoneInfoNotFoundError, ValueError):
        return False


def local_iso(t: Time, tz_name: str) -> str:
    """Format an astropy Time as 'YYYY-MM-DD HH:MM' in the given IANA timezone.

    Falls back to UTC (with a trailing " UTC") if `tz_name` isn't recognized.
    """
    try:
        return t.to_datetime(timezone=ZoneInfo(tz_name)).strftime("%Y-%m-%d %H:%M")
    except (ZoneInfoNotFoundError, ValueError):
        return t.utc.iso[:16] + " UTC"
