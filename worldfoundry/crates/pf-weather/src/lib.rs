//! pf-weather: WEATHER-V0 and COLDPOOL-1 (lane WEATHER-COLDPOOL, 2026-10-21).
//!
//! Operator ruling 2026-10-21 (model specification design review's perception mechanics,
//! `unpublished-design-record review-mechanics/REPORT.md` 1.3 and 1.4; the 21:54
//! STREAMING_DESIGN addendum names the visibility model as OPEN).
//!
//! WHAT THIS IS: an opt-in, pure function from canonical inputs to a weather
//! record. Inputs: the published T0 climate fields (`pf.climate/v0`:
//! temperature mean and seasonal half range, precipitation and its
//! orographic part, aridity index, prevailing wind u/v) and T0 elevation; the
//! world seed; the world hour; a planet point (lat, lon) and optionally an
//! elevation. Output: cloud cover, precipitation type and rate, wind, air
//! temperature, dewpoint and relative humidity, the extinction budget and
//! meteorological visibility V = 3.912 / beta (Koschmieder), plus the night
//! radiative-cooling factor COLDPOOL-1 integrates.
//!
//! CANON (spec decisions 4, 9, 11, 23; ADR 005): weather is ontic and
//! deterministic from (world seed, clock, place); this crate computes it as a
//! pure function, so the kernel's perception query and any presentation
//! (UE height fog, SkyAtmosphere, the terrain streamer's load radius) read the
//! same numbers. Renderer state is never an input. Nothing here writes canon;
//! a kernel adapter that mints weather records would cite law id
//! `pf-weather-v0` and the constants digest. World time is not compute time:
//! the hour index is supplied (tick -> hour mapping OPEN, see `sun`).
//!
//! RANDOMNESS: `vz.keyed-rng/v1` only, law-set scope `pf-weather-v0`, keyed on
//! planet identity (a planet-fixed 3D cube lattice for disturbances, the T0
//! cell id for hourly jitter) and the time step; never on a chart. Declared
//! purposes: `weather.strat.{count,cx,cy,cz,start,dur,radius,thin,intensity}`,
//! `weather.conv.{same}`, `weather.hour.{rate,cloud,wind-u,wind-v}`.
//!
//! MODEL (90 %, not 99.9 %): two marked Poisson disc processes (synoptic
//! systems 200-700 km, convective showers 8-40 km) drifting with a latitude
//! steering wind; each disc is kept at a place with a probability set from
//! that place's climatology so that the wet-hour fraction is exact
//! (Poisson void probability) and the long-run precipitation equals the T0
//! annual precipitation in expectation (see `disturb` and `state`). Daily and
//! seasonal temperature cycles from the T0 mean and half range; humidity from
//! the aridity index with vapour pressure held through the day; wind from the
//! T0 prevailing wind, quieter under high pressure and at night, plus the
//! cyclonic wind of passing systems; extinction = Rayleigh + humid haze
//! (hygroscopic growth) + rain/snow (NCEP CALVIS power laws) + dust (arid,
//! windy, dry) + open-ground radiation fog.
//!
//! DETERMINISM: no clock, no hash-map iteration in outputs, no threads.
//! Transcendental functions are the platform's (exp, ln, powf, sin, cos), so
//! byte identity is a guarantee on a pinned build (as pf-surface's climate).
//!
//! NOT MODELLED (V0): precipitation seasonality (the T0 tier has none),
//! interannual variability, sea/advection fog, hill fog (cloud base below a
//! summit), snow cover, rain-wetted dust suppression after the rain stops,
//! frontal temperature structure beyond a warm/cold swing.
//!
//! V1 (lane WEATHER-V1, 2026-10-21, all opt-in; V0 outputs unchanged when the
//! keys are off): `fronts` (frontal band families and convective clusters,
//! `fronts=1`, optional orographic modulation `oro=1`) and `coldflow`
//! (COLDPOOL cold-air drainage over the chart-wide river network, per-airshed
//! lids, `coldflow=1`).

pub mod climate;
pub mod coldflow;
pub mod coldpool;
pub mod constants;
pub mod disturb;
pub mod fronts;
pub mod state;
mod state_v1;
pub mod sun;

pub use climate::{ClimateGrid, ClimateSample};
pub use constants::WeatherConstantsV0;
pub use state::{PrecipType, WeatherModel, WeatherState, LAW_ID};
