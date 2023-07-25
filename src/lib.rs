// https://github.com/radcheb/Adhan/tree/master/C/adhan/src
use chrono::{DateTime, NaiveDate, Utc, Duration, Datelike, Timelike};
use anyhow::{ Result, anyhow };
use std::error::Error;

const PI: f64 = 3.14159265358979323846;

enum Prayers {
    None,
    Fajr,
    Sunrise,
    Dhur,
    Asr,
    Maghrib,
    Isha,
}

enum Madhab {
    Shafi,
    Hanafi,
}

enum CalculationMethod {
    MuslimWorldLeague,
    Egyptian,
    Karachi,
    UmmAlQura,
    Gulf,
    Dubai,
    MoonSightingCommittee,
    NorthAmerica,
    Kuwait,
    Qatar,
    Singapoor,
    Other,
}

enum HighLatitudeRule {
    MiddleOfTheNight,
    SeventhOfTheNight,
    TwilightAngle,
}

enum ShadowLength {
    Single,
    Double,
}

struct PrayerAdjustments {
    fajr: i32,
    sunrise: i32,
    dhuhr: i32,
    asr: i32,
    maghrib: i32,
    isha: i32,
}

impl PrayerAdjustments {
    fn init() -> Self {
        PrayerAdjustments {
            fajr: 0,
            sunrise: 0,
            dhuhr: 0,
            asr: 0,
            maghrib: 0,
            isha: 0,
        }
    }
}

fn get_shadow_length(madhab: Madhab) -> ShadowLength {
    return match madhab {
        Madhab::Shafi => ShadowLength::Single,
        Madhab::Hanafi => ShadowLength::Double,
    };
}

struct NightPortions {
    fajr: f64,
    isha: f64,
}

impl NightPortions {
    fn get(parameters: CalculationParameters) -> Self {
        match parameters.high_latitude_rule {
            HighLatitudeRule::MiddleOfTheNight => {
                return NightPortions {
                    fajr: (1.0 / 2.0),
                    isha: (1.0 / 2.0),
                }
            }
            HighLatitudeRule::SeventhOfTheNight => {
                return NightPortions {
                    fajr: (1.0 / 7.0),
                    isha: (1.0 / 7.0),
                }
            }
            HighLatitudeRule::TwilightAngle => {
                return NightPortions {
                    fajr: parameters.fajr_angle / 60.0,
                    isha: parameters.isha_angle / 60.0,
                }
            }
        }
    }
}

struct CalculationParameters {
    calculation_method: CalculationMethod,
    fajr_angle: f64,
    isha_angle: f64,
    isha_interval: i32,
    madhab: Madhab,
    high_latitude_rule: HighLatitudeRule,
    prayer_adjustments: PrayerAdjustments,
}

impl CalculationParameters {
    fn new(&self) -> Self {
        Self::with_parameters(
            CalculationMethod::Other,
            0.0,
            0.0,
            0,
            Madhab::Shafi,
            HighLatitudeRule::TwilightAngle,
            PrayerAdjustments::init(),
        )
    }

    fn with_parameters(
        calculation_method: CalculationMethod,
        fajr_angle: f64,
        isha_angle: f64,
        isha_interval: i32,
        madhab: Madhab,
        high_latitude_rule: HighLatitudeRule,
        prayer_adjustments: PrayerAdjustments,
    ) -> Self {
        CalculationParameters {
            calculation_method: calculation_method,
            fajr_angle: fajr_angle,
            isha_angle: isha_angle,
            isha_interval: isha_interval,
            madhab: madhab,
            high_latitude_rule: high_latitude_rule,
            prayer_adjustments: prayer_adjustments,
        }
    }

    fn by_method(method: CalculationMethod) -> Self {
        match method {
            CalculationMethod::MuslimWorldLeague => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 18.0,
                    isha_angle: 17.0,
                    isha_interval: 0,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
            CalculationMethod::Egyptian => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 20.0,
                    isha_angle: 18.0,
                    isha_interval: 0,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
            CalculationMethod::Karachi => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 18.0,
                    isha_angle: 18.0,
                    isha_interval: 0,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
            CalculationMethod::UmmAlQura => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 18.5,
                    isha_angle: 0.0,
                    isha_interval: 90,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
            CalculationMethod::Gulf => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 19.5,
                    isha_angle: 0.0,
                    isha_interval: 90,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
            CalculationMethod::Dubai => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 0.0,
                    isha_angle: 0.0,
                    isha_interval: 0,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
            CalculationMethod::MoonSightingCommittee => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 18.0,
                    isha_angle: 18.0,
                    isha_interval: 0,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
            CalculationMethod::NorthAmerica => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 15.0,
                    isha_angle: 15.0,
                    isha_interval: 0,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
            CalculationMethod::Kuwait => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 18.0,
                    isha_angle: 17.5,
                    isha_interval: 0,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
            CalculationMethod::Qatar => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 18.0,
                    isha_angle: 0.0,
                    isha_interval: 90,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
            CalculationMethod::Singapoor => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 0.0,
                    isha_angle: 0.0,
                    isha_interval: 0,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
            CalculationMethod::Other => {
                return CalculationParameters {
                    calculation_method: method,
                    fajr_angle: 0.0,
                    isha_angle: 0.0,
                    isha_interval: 0,
                    madhab: Madhab::Shafi,
                    high_latitude_rule: HighLatitudeRule::TwilightAngle,
                    prayer_adjustments: PrayerAdjustments::init(),
                }
            }
        }
    }
}

#[derive(Clone)]
struct Coordinates {
    latitude: f64,
    longitude: f64,
}

impl Coordinates {
    fn new(latitude: f64, longitude: f64) -> Self {
        Coordinates {
            latitude: latitude,
            longitude: longitude,
        }
    }
}

struct SolarCoordinates {
    declination: f64,
    right_ascension: f64,
    apparent_side_realtime: f64,
}


type JulianDay = f64;

impl SolarCoordinates {
    fn new(julian_day: JulianDay) -> Self {
        let T = self::julian_century(julian_day);
        let L0 = mean_solar_longitude(T);
        let Lp = mean_lunar_longitude(T);
        let omega = ascending_lunar_node_longitude(T);
        let lambda = to_radius(apparent_solar_longitude(T, L0));
        let theta0 = mean_side_real_time(T);
        let delta_psi = nutation_in_longitude(L0, Lp, omega);
        let delta_epsilon = nutation_in_obliquity(L0, Lp, omega);
        let epsilon0 = mean_obliquity_of_the_ecliptic(T);
        let epsilonapp = to_radius(apparent_obliquity_of_the_ecliptic(T, epsilon0));
        /* Equation from Astronomical Algorithms page 165 */
        let declination = to_degrees(((epsilonapp).sin() * (lambda).sin()).asin());
        /* Equation from Astronomical Algorithms page 165 */
        let right_ascension = unwind_angle(to_degrees(
            ((epsilonapp).cos() * (lambda).sin()).atan2((lambda).cos()),
        ));

        /* Equation from Astronomical Algorithms page 88 */
        let apparent_side_real_time0 = theta0
            + (((delta_psi * 3600.0) * (to_radius(epsilon0 + delta_epsilon)).cos()) / 3600.0);
        SolarCoordinates {
            declination: declination,
            right_ascension: right_ascension,
            apparent_side_realtime: apparent_side_real_time0,
        }
    }
}

struct SolarTime {
    transit: f64,
    sunrise: f64,
    sunset: f64,
    observer: Coordinates,
    solar: SolarCoordinates,
    prev_solar: SolarCoordinates,
    next_solar: SolarCoordinates,
    approximate_transit: f64,
}

impl SolarTime {
    fn new(time: TimeComponent, date: DateComponent, coordinates: Coordinates) -> Result<Self, Box<dyn Error>> {
        let maybe_datetime = NaiveDate::from_ymd_opt(date.year, date.month as u32, date.day as u32)
                                        .ok_or(anyhow!("oops"))?
                                        .and_hms_opt(time.hours as u32, time.minutes as u32, time.seconds as u32).ok_or(anyhow!("on no"))?;

        let today = DateTime::<Utc>::from_utc(maybe_datetime, Utc);


        let tomorrow = today + Duration::days(1);
        let yesterday = today + Duration::days(-1);

        let solar = SolarCoordinates::new(_julian_day(
            today.year(),
            today.month().try_into()?,
            today.day().try_into()?,
            today.hour().into(),
        ));

        let prev_solar = SolarCoordinates::new(_julian_day(
            yesterday.year(),
            yesterday.month().try_into()?,
            yesterday.day().try_into()?,
            yesterday.hour().into(),
        ));

        let next_solar = SolarCoordinates::new(_julian_day(
            tomorrow.year(),
            tomorrow.month().try_into()?,
            tomorrow.day().try_into()?,
            tomorrow.hour().into(),
        ));

        let approximate_transit = get_approximate_transit(coordinates.longitude, solar.apparent_side_realtime,
            solar.right_ascension);

        let solarAltitude = -50.0 / 60.0;

        let transit = corrected_transit(approximate_transit, coordinates.longitude,
        solar.apparent_side_realtime, solar.right_ascension, prev_solar.right_ascension,
        next_solar.right_ascension);

        let sunrise = corrected_hour_angle(approximate_transit, solarAltitude,
        coordinates.clone(), false, solar.apparent_side_realtime, solar.right_ascension,
        prev_solar.right_ascension, next_solar.right_ascension, solar.declination,
        prev_solar.declination, next_solar.declination);

        let sunset = corrected_hour_angle(approximate_transit, solarAltitude,
        coordinates.clone(), true, solar.apparent_side_realtime, solar.right_ascension,
        prev_solar.right_ascension, next_solar.right_ascension, solar.declination,
        prev_solar.declination, next_solar.declination);

        Ok(SolarTime {
            transit: transit,
            sunrise: sunrise,
            sunset: sunset,
            observer: coordinates.clone(),
            solar: solar,
            prev_solar: prev_solar,
            next_solar: next_solar,
            approximate_transit: approximate_transit,
        })
    }
}

fn normalize_with_bound(value: f64, max: f64) -> f64 {
    value - (max * (value / max).floor())
}

fn unwind_angle(value: f64) -> f64 {
    self::normalize_with_bound(value, 360.0)
}

fn closest_angle(angle: f64) -> f64 {
    if angle >= -180.0 && angle <= 180.0 {
        return angle;
    }
    return angle - (360.0 * (angle / 360.0).round());
}

struct TimeComponent {
    hours: i32,
    minutes: i32,
    seconds: i32,
}

struct DateComponent {
    day: i32,
    month: i32,
    year: i32,
}

impl TimeComponent {
    fn new(hours: i32, minutes: i32, seconds: i32) -> Self {
        TimeComponent {
            hours: hours,
            minutes: minutes,
            seconds: seconds
        }
    }

    fn is_valid_time(&self) -> bool {
        !(
            self.hours == -1 &&
            self.minutes == -1 &&
            self.seconds == -1
        )
    }

    fn from_f64(value: f64) -> Self {
        let hours = value.floor() as i32;
        let minutes = ((value - hours as f64) * 60.0).floor() as i32;
        let seconds = ((value - (hours as f64 + minutes as f64 / 60.0)) * 60.0 * 60.0).floor() as i32;
        TimeComponent {
            hours: hours,
            minutes: minutes,
            seconds: seconds
        }
    }
}


struct PrayerTimes {
    fajr: TimeComponent,
    sunrise: TimeComponent,
    dhuhr: TimeComponent,
    asr: TimeComponent,
    maghrib: TimeComponent,
    isha: TimeComponent,
}

/* -------------------------------------------------------------------------------------------------
    Astronomical calculations
--------------------------------------------------------------------------------------------------*/

fn to_radius(deg: f64) -> f64 {
    deg * (PI / 180.0)
}

fn to_degrees(radians: f64) -> f64 {
    radians * (180.0 / PI)
}

fn mean_solar_longitude(t: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 163 */
    let term1 = 280.4664567;
    let term2 = 36000.76983 * t;
    let term3 = 0.0003032 * f64::powf(t, 2.0);
    let l0 = term1 + term2 + term3;
    self::unwind_angle(l0)
}

fn mean_lunar_longitude(t: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 144 */
    let term1 = 218.3165;
    let term2 = 481267.8813 * t;
    let lp = term1 + term2;
    self::unwind_angle(lp)
}

fn apparent_solar_longitude(t: f64, L0: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 164 */
    let longitude = L0 + self::solar_equation_of_the_center(t, self::mean_solar_anomaly(t));
    let omega = 125.04 - (1934.136 * t);
    let lambda = longitude - 0.00569 - (0.00478 * (self::to_radius(omega)).sin());
    self::unwind_angle(lambda)
}

fn ascending_lunar_node_longitude(t: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 144 */
    let term1 = 125.04452;
    let term2 = 1934.136261 * t;
    let term3 = 0.0020708 * f64::powf(t, 2.0);
    let term4 = f64::powf(t, 3.0) / 450000.0;
    let omega = term1 - term2 + term3 + term4;
    self::unwind_angle(omega)
}

fn mean_solar_anomaly(t: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 163 */
    let term1 = 357.52911;
    let term2 = 35999.05029 * t;
    let term3 = 0.0001537 * f64::powf(t, 2.0);
    let m = term1 + term2 - term3;
    self::unwind_angle(m)
}

fn solar_equation_of_the_center(t: f64, m: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 164 */
    let mrad = self::to_radius(m);
    let term1 = (1.914602 - (0.004817 * t) - (0.000014 * f64::powf(t, 2.0))) * (mrad).sin();
    let term2 = (0.019993 - (0.000101 * t)) * (2.0 * mrad).sin();
    let term3 = 0.000289 * (3.0 * mrad).sin();
    term1 + term2 + term3
}

fn mean_obliquity_of_the_ecliptic(t: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 147 */
    let term1 = 23.439291;
    let term2 = 0.013004167 * t;
    let term3 = 0.0000001639 * f64::powf(t, 2.0);
    let term4 = 0.0000005036 * f64::powf(t, 3.0);
    term1 - term2 - term3 + term4
}

fn apparent_obliquity_of_the_ecliptic(t: f64, epsilon0: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 165 */
    let o = 125.04 - (1934.136 * t);
    epsilon0 + (0.00256 * (self::to_radius(o)).cos())
}

fn mean_side_real_time(t: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 165 */
    let jd = (t * 36525.0) + 2451545.0;
    let term1 = 280.46061837;
    let term2 = 360.98564736629 * (jd - 2451545.0);
    let term3 = 0.000387933 * f64::powf(t, 2.0);
    let term4 = f64::powf(t, 3.0) / 38710000.0;
    let theta = term1 + term2 + term3 - term4;
    self::unwind_angle(theta)
}

fn nutation_in_longitude(lo: f64, lp: f64, omega: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 144 */
    let term1 = (-17.2 / 3600.0) * (self::to_radius(omega)).sin();
    let term2 = (1.32 / 3600.0) * (2.0 * self::to_radius(lo)).sin();
    let term3 = (0.23 / 3600.0) * (2.0 * self::to_radius(lp)).sin();
    let term4 = (0.21 / 3600.0) * (2.0 * self::to_radius(omega)).sin();
    term1 - term2 - term3 + term4
}

fn nutation_in_obliquity(lo: f64, lp: f64, omega: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 144 */
    let term1 = (9.2 / 3600.0) * (self::to_radius(omega)).cos();
    let term2 = (0.57 / 3600.0) * (2.0 * self::to_radius(lo)).cos();
    let term3 = (0.10 / 3600.0) * (2.0 * self::to_radius(lp)).cos();
    let term4 = (0.09 / 3600.0) * (2.0 * self::to_radius(omega)).cos();
    term1 + term2 + term3 - term4
}

fn altitude_of_celestial_body(phi: f64, delta: f64, H: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 93 */
    let term1 = (self::to_radius(phi).sin()) * (self::to_radius(delta).sin());
    let term2 =
        (self::to_radius(phi).cos()) * (self::to_radius(delta).cos()) * (self::to_radius(H).cos());
    self::to_degrees((term1 + term2).asin())
}

fn get_approximate_transit(l: f64, theta0: f64, alpha2: f64) -> f64 {
    /* Equation from page Astronomical Algorithms 102 */
    let lw = l * -1.0;
    self::normalize_with_bound((alpha2 + lw - theta0) / 360.0, 1.0)
}

fn corrected_transit(m0: f64, l: f64, theta0: f64, alpha2: f64, alpha1: f64, alpha3: f64) -> f64 {
    /* Equation from page Astronomical Algorithms 102 */
    let lw = l * -1.0;
    let theta = self::unwind_angle(theta0 + (360.985647 * m0));
    let alpha = self::unwind_angle(self::interpolate_angles(
        /* value */ alpha2, /* previousValue */ alpha1, /* nextValue */ alpha3,
        /* factor */ m0,
    ));
    let h = self::closest_angle(theta - lw - alpha);
    let deltam = h / -360.0;
    (m0 + deltam) * 24.0
}

fn corrected_hour_angle(
    m0: f64,
    h0: f64,
    coordinates: Coordinates,
    after_transit: bool,
    theta0: f64,
    alpha2: f64,
    alpha1: f64,
    alpha3: f64,
    delta2: f64,
    delta1: f64,
    delta3: f64,
) -> f64 {
    /* Equation from page Astronomical Algorithms 102 */
    let Lw = coordinates.longitude * -1.0;
    let term1 = (self::to_radius(h0).sin())
        - ((self::to_radius(coordinates.latitude).sin()) * (self::to_radius(delta2).sin()));
    let term2 = (self::to_radius(coordinates.latitude).cos()) * (self::to_radius(delta2).cos());
    let H0 = self::to_degrees((term1 / term2).acos());
    let m = if after_transit {
        m0 + (H0 / 360.0)
    } else {
        m0 - (H0 / 360.0)
    };
    let theta = self::unwind_angle(theta0 + (360.985647 * m));
    let alpha = self::unwind_angle(self::interpolate_angles(
        /* value */ alpha2, /* previousValue */ alpha1, /* nextValue */ alpha3,
        /* factor */ m,
    ));
    let delta = self::interpolate(
        /* value */ delta2, /* previousValue */ delta1, /* nextValue */ delta3,
        /* factor */ m,
    );
    let H = theta - Lw - alpha;
    let h = self::altitude_of_celestial_body(
        /* observerLatitude */ coordinates.latitude,
        /* declination */ delta,
        /* localHourAngle */ H,
    );
    let term3 = h - h0;
    let term4 = 360.0
        * (self::to_radius(delta)).cos()
        * (self::to_radius(coordinates.latitude)).cos()
        * (self::to_radius(H)).sin();
    let deltam = term3 / term4;
    (m + deltam) * 24.0
}

fn interpolate(y2: f64, y1: f64, y3: f64, n: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 24 */
    let a = y2 - y1;
    let b = y3 - y2;
    let c = b - a;
    y2 + ((n / 2.0) * (a + b + (n * c)))
}

fn interpolate_angles(y2: f64, y1: f64, y3: f64, n: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 24 */
    let a = self::unwind_angle(y2 - y1);
    let b = self::unwind_angle(y3 - y2);
    let c = b - a;
    y2 + ((n / 2.0) * (a + b + (n * c)))
}

/* -------------------------------------------------------------------------------------------------
    Date utils
--------------------------------------------------------------------------------------------------*/

fn _julian_day(year: i32, month: i32, day: i32, hours: f64) -> JulianDay {
    /* Equation from Astronomical Algorithms page 60 */
    // NOTE: Integer conversion is done intentionally for the purpose of decimal truncation
    let Y = if month > 2 { year } else { year - 1 };
    let M = if month > 2 { month } else { month + 12 };
    let D = day as f64 + (hours / 24.0);

    let A = Y / 100;
    let B = 2 - A + (A / 4);

    let i0 = (365.25 * (Y + 4716) as f64) as i32;
    let i1 = (30.6001 * (M + 1) as f64) as i32;
    return i0 as f64 + i1 as f64 + D as f64 + B as f64 - 1524.5;
}

fn julian_day(year: i32, month: i32, day: i32) -> JulianDay {
    self::_julian_day(year, month, day, 0.0)
}

// fn julianDay2(const struct tm* date) -> f64 {
//     return _julianDay(date->tm_year + 1900,
//                      date->tm_mon + 1, date->tm_mday,
//                      date->tm_hour + date->tm_min / 60.0);
//}

fn julian_century(jd: f64) -> f64 {
    /* Equation from Astronomical Algorithms page 163 */
    (jd - 2451545.0) / 36525.0
}

/* -------------------------------------------------------------------------------------------------
    lib exports
--------------------------------------------------------------------------------------------------*/
pub fn ComputePrayers() -> String {
    "TODO".to_string()
}