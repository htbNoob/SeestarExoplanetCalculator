"""Cloud-cover / precipitation / temperature forecast for tonight's dark
window, via the free Open-Meteo API (no key required)."""
from dataclasses import dataclass
from typing import List, Optional

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy.time import Time

from seestar_site import DarkWindow

FORECAST_URL = "https://api.open-meteo.com/v1/forecast"


@dataclass
class HourlyForecast:
    time_utc: List[str]        # ISO8601 UTC hour strings
    cloudcover: List[float]    # %
    precip_prob: List[float]   # %
    temperature_c: List[float]
    dewpoint_c: List[float]
    windspeed_kmh: List[float]


def fetch_hourly_forecast(lat: float, lon: float, forecast_days: int = 2,
                           timeout: int = 15) -> HourlyForecast:
    """Hourly cloud cover / precipitation probability / temperature / dew
    point / wind speed forecast for a location, covering `forecast_days`
    ahead (default 2, enough to cover the 36 h dark-window scan)."""
    resp = requests.get(FORECAST_URL, params={
        "latitude": lat, "longitude": lon,
        "hourly": "cloudcover,precipitation_probability,temperature_2m,dewpoint_2m,windspeed_10m",
        "timezone": "UTC",
        "forecast_days": forecast_days,
    }, timeout=timeout)
    resp.raise_for_status()
    h = resp.json()["hourly"]
    return HourlyForecast(
        time_utc=h["time"],
        cloudcover=h["cloudcover"],
        precip_prob=h["precipitation_probability"],
        temperature_c=h["temperature_2m"],
        dewpoint_c=h["dewpoint_2m"],
        windspeed_kmh=h["windspeed_10m"],
    )


@dataclass
class DarkWindowWeather:
    mean_cloudcover: float
    max_cloudcover: float
    mean_precip_prob: float
    min_temperature_c: float
    max_windspeed_kmh: float
    rating: str   # "Good" / "Fair" / "Poor"


def summarize_dark_window_weather(forecast: HourlyForecast,
                                   dark_window: DarkWindow) -> Optional[DarkWindowWeather]:
    """Aggregate the hourly forecast over just tonight's dark-window hours."""
    jds = Time(forecast.time_utc, format="isot", scale="utc").jd
    mask = (jds >= dark_window.ds_jd) & (jds <= dark_window.de_jd)
    if not mask.any():
        return None

    cc = np.asarray(forecast.cloudcover, dtype=float)[mask]
    pp = np.asarray(forecast.precip_prob, dtype=float)[mask]
    temp = np.asarray(forecast.temperature_c, dtype=float)[mask]
    wind = np.asarray(forecast.windspeed_kmh, dtype=float)[mask]

    mean_cc = float(np.nanmean(cc))
    mean_pp = float(np.nanmean(pp))

    if mean_cc < 20 and mean_pp < 20:
        rating = "Good"
    elif mean_cc < 60 and mean_pp < 40:
        rating = "Fair"
    else:
        rating = "Poor"

    return DarkWindowWeather(
        mean_cloudcover=mean_cc,
        max_cloudcover=float(np.nanmax(cc)),
        mean_precip_prob=mean_pp,
        min_temperature_c=float(np.nanmin(temp)),
        max_windspeed_kmh=float(np.nanmax(wind)),
        rating=rating,
    )


def build_weather_figure(forecast: HourlyForecast, dark_window: DarkWindow) -> plt.Figure:
    """Two-panel forecast plot (cloud cover + precip. probability / temperature
    + dew point) with tonight's dark window shaded, in UTC."""
    dt = Time(forecast.time_utc, format="isot", scale="utc").datetime
    dark_start_dt = dark_window.dark_start.utc.datetime
    dark_end_dt = dark_window.dark_end.utc.datetime

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    fig.suptitle(f"Forecast - {dark_window.site.name} (UTC)", fontsize=12, fontweight="bold")

    ax1.axvspan(dark_start_dt, dark_end_dt, color="#dbe7f7", zorder=0,
                label="Tonight's dark window")
    ax1.plot(dt, forecast.cloudcover, color="steelblue", lw=2, label="Cloud cover (%)")
    ax1.plot(dt, forecast.precip_prob, color="tomato", lw=1.5, ls="--",
              label="Precipitation prob. (%)")
    ax1.set_ylabel("%")
    ax1.set_ylim(0, 100)
    ax1.grid(alpha=0.2)
    ax1.legend(fontsize=8, loc="upper right")

    ax2.axvspan(dark_start_dt, dark_end_dt, color="#dbe7f7", zorder=0)
    ax2.plot(dt, forecast.temperature_c, color="darkorange", lw=2, label="Temperature (C)")
    ax2.plot(dt, forecast.dewpoint_c, color="teal", lw=1.5, ls="--", label="Dew point (C)")
    ax2.set_ylabel("deg C")
    ax2.grid(alpha=0.2)
    ax2.legend(fontsize=8, loc="upper right")

    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %H:%M"))
    fig.autofmt_xdate()
    plt.tight_layout()
    return fig
