# Prayer

_High precision prayer time library_

Adhan is a well tested and well documented library for calculating Islamic prayer times. All astronomical calculations
are high precision equations directly from the book
[“Astronomical Algorithms” by Jean Meeus](http://www.willbell.com/math/mc1.htm). This book is recommended
by the Astronomical Applications Department of the U.S. Naval Observatory and the Earth System Research Laboratory
of the National Oceanic and Atmospheric Administration.

## Languages

Adhan is based on the following port:

- C (C99) [(Usage and examples)](https://github.com/radcheb/Adhan/tree/master/C/adhan)

## Usage

Add the library dependency into your `Cargo.toml` file for example:

```toml
[dependencies]
prayer = { git "https://github.com/AnharHussainMiah/prayer.git" }
```

Next import the library:

```rust
use prayer::{
    CalculationMethod, CalculationParameters, Coordinates, DateComponent, PrayerTimes,
    TimeComponent,
};
```

A basic example of computing the prayer times using a predefined `CalculationMethod` :

```rust
if let Ok(prayers) = prayer::PrayerTimes::new(
        Coordinates::new(51.509865, -0.118092), // latitude, longitude
        TimeComponent::new(0, 0, 0),            // hours, minutes, seconds
        DateComponent::new(1, 8, 2023),         // day, month, year
        CalculationParameters::by_method(CalculationMethod::MuslimWorldLeague),
    ) {
        println!(
            "Fajr -> {}:{}:{}",
            prayers.fajr.hours, prayers.fajr.minutes, prayers.fajr.seconds
        );
        println!(
            "Sunrise -> {}:{}:{}",
            prayers.sunrise.hours, prayers.sunrise.minutes, prayers.sunrise.seconds
        );
        println!(
            "Dhuhr -> {}:{}:{}",
            prayers.dhuhr.hours, prayers.dhuhr.minutes, prayers.dhuhr.seconds
        );
        println!(
            "Asr -> {}:{}:{}",
            prayers.asr.hours, prayers.asr.minutes, prayers.asr.seconds
        );
        println!(
            "Maghrib -> {}:{}:{}",
            prayers.maghrib.hours, prayers.maghrib.minutes, prayers.maghrib.seconds
        );
        println!(
            "isha -> {}:{}:{}",
            prayers.isha.hours, prayers.isha.minutes, prayers.isha.seconds
        );
    }
```

### Using custom calculation parameters

[TODO]

## Contributing

Prayer is made publicly available to provide a well tested and well documented library for Islamic prayer times to all
developers. We accept feature contributions provided that they are properly documented and include the appropriate
unit tests. We are also looking for contributions in the form of unit tests of of prayer times for different
locations, we do ask that the source of the comparison values be properly documented.

## License

Prayer is available under the MIT license. See the LICENSE file for more info.
