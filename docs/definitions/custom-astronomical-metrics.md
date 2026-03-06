# Custom Astronomical Metric Definitions

There are some customized metrics contained in the datasets that may need special reference of explanation in analysis descriptions:

**Solaration:** `solar_secs` --- Seconds elapsed since the preceding winter solstice. Tracks Earth's orbital position and maps directly to the Sun's declination angle, replacing the arbitrary civil calendar with an astronomically anchored annual cycle. Full metric description details can be found at https://github.com/jakeYeager/nornir-urd-two/blob/main/metric_info_solar.md

**Lunation:** `lunar_secs` --- Seconds elapsed since the preceding new moon. Simultaneously encodes the Moon's orbital position and the Sun-Moon phase alignment, capturing the spring/neap tidal cycle that drives both ocean tides and the lesser-known earth tides (periodic crustal deformation). Full metric description details can be found at https://github.com/jakeYeager/nornir-urd-two/blob/main/metric_info_lunar.md

**Midnight:** `midnight_secs` --- Seconds elapsed since local solar midnight at the event's longitude. Uses a pure longitude-based time offset rather than civil time zones, providing an exact measure of the Sun's rotational position relative to the earthquake location. Full metric description details can be found at https://github.com/jakeYeager/nornir-urd-two/blob/main/metric_info_midnight.md
