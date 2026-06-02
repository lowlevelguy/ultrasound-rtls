# Ultrasound-Based Indoors Real-Time Location System
This repository contains all the necessary code for the implementation of a prototype 2D Real-Time Location System (RTLS) with time-position history storage on a cloud database and a web interface for playback.

## Prerequisites
This prototype consists of four different components:
- a position tracker: running on a Teensy 4.0
- a history logger: running on an ESP32
- an MQTT and InfluxDB server
- an HTML/JS interface (`dashboard.html`)

Moreover, the position tracker expects at least 3 ultrasound ranging modules, communicating over UART, to function.

## Running
To build, upload and monitor the tracker component:
```
pio run -t upload -t monitor -e tracker
```

To build, upload and monitor the logger component:
```
pio run -t upload -t monitor -e logger
```

To setup the MQTT-InfluxDB data pipeline, make sure you have Mosquitto, InfluxDB2 and Telegraf installed on your system. Then, copy `server/mosquitto.conf` into `/etc/mosquitto/conf.d/rtls.conf` (or whatever you wish to name it), and `server/telegraf.conf` into `/etc/telegraf/telegraf.d/rtls.conf` (or whatever you wish to name it). Once that's done, simply run the bridge script:
```
python3 server/bridge.py
```
That will serve as the translation unit between the ESP32's struct-packed dumps and telegraf.

To run the dashboard, simply open a python server:
```
python3 -m http.server 8080
```
and open `localhost:8080`.

## Description
### Position Tracker
This component, running on a [Teensy 4.0](https://www.pjrc.com/store/teensy40.html), is tasked with collecting distance information from 3 or 4 ultrasound ranging modules (depending on your configuration). The main reason for this choice of board is its 7 hardware UART interfaces, sufficient to drive our ranging modules.

#### Ultrasound Ranging Module
We make use of 3 [AJ-SR04M modules](https://www.fabian.com.mt/viewer/42585/pdf.pdf), along with 1 module named SR04M-2 of unknown origin; though online reviews indicate it's functionally identical to the [JSN SR04T 2.0](https://github.com/HamidSaffari/JSN-SR04T/blob/master/JSN-SR04T-2.0.pdf).

We the AJ-SR04M modules in mode 4: low power UART (9600 baud rate). The module sits idle, until it receives a trigger byte `0x01`, upon which it initiates its ranging function, and reports back its measurement in a 4-bytes binary format: `0xFF, high_byte, low_byte, checksum`. The checksum is computed as the sum of `high_byte` and `low_byte` with no carry. The AJ-SR04M is reported in the datasheet to be able to achieve ranging upto 8m, however testing in real life revealed the true distance of detection highly depends on the size of the object. In our case, it ended up only detecting our object upto 4m away from its probe.

The SR04M-2 module is similar to the AJ-SR04M module. Though it is called mode 3, it also presents a low power UART mode, also with a 9600 baud rate. The only differences are that its trigger byte is `0x55`, and the checksum is computed as the sum of high and low bytes, as well as the start byte `0xFF`; yet again with no carry. This is equivalent to `high_byte + low_byte - 1`. Moreover, its (or rather, the JSN SR04T 2.0's) documented maximum range is 6m, as opposed to 8m for the AJ-SR04M. Again, in practice, this highly depends on the dimensions of the target. For ours, we were only able to reliably do range detection upto 2.5m away from the probe.

One thing to note is that the AJ-SR04M modules do not seem to exclusively accept `0x01` trigger bytes; for example, a `0x55` trigger byte will function all the same. It seems that they only check for the LSB being set, rather than the value of the entire byte being 1.

#### Module Calibration
Upon testing, it seems as if each sensor has an innate error term that is linear in distance. We took 5 to 8 measurements at known positions, 0.5m apart each, for each module. After coming to a convention on the module labeling, we approximated the following linear bias models:

```
measured[0] = 0.9896 * dist - 21.5
measured[1] = 0.9876 * dist - 12.6
measured[2] = 0.9713 * dist - 19.9
measured[3] = 1.0700 * dist - 43.9     // SR04M-2
```

As a result, we made sure to correct the bias for each reading in our implementation in accordance with this model.

#### Positioning Algorithms
We implement five positioning algorithms: trilateration, Ordinary Least Squares (OLS), Feasible Generalized Least Squares (FGLS), and Maximum Likelihood Estimation (MLE) with 4 anchors, and MLE with 3 anchors. Below are the execution time benchmarks for each, taken over 300 iterations on the Teensy 4.0, excluding a first discarded iteration.

|             Algorithm             | Min Cycle Count | Max Cycle Count |
| :-------------------------------: | :-------------: | :-------------: |
|           Trilateration           |   34 (56.7 ns)  |   35 (58.3 ns)  |
|                OLS                |   270 (450 ns)  |   270 (450 ns)  |
|               FGLS                |  1074 (1.79 us) |  1102 (1.84 us) |
| MLE-4 (4 Gauss-Newton iterations) |  3974 (6.6 us)  |  4061 (6.8 us)  |
| MLE-3 (8 Guass-Newton iterations) |   6616 (11 us)  |  6795 (11.3 us) |

Likewise, below are the error benchmarks, taken over 500 iterations, for simulated additive white noise-d input distance measurements with a standard deviation of 5 mm, in terms of distance between the true position and the estimated one:

|             Algorithm             | Mean Error (mm) | Standard Deviation (mm) |
| :-------------------------------: | :-------------: | :---------------------: |
|           Trilateration           |       3.47      |          0.02           |
|                OLS                |       2.68      |          2.32           |
|               FGLS                |       2.40      |          1.88           |
| MLE-4 (4 Gauss-Newton iterations) |       4.17      |          2.28           |
| MLE-3 (8 Guass-Newton iterations) |       2.86      |          2.57           |

### ESP32
We use an off-brand ESP32 board built around the ESP32 Wrover SoC. As such, it comes natively with Bluetooth and WiFi connectivity, the latter of which we need to make use of cloud databases.

