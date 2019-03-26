# COVIS Data Description

## Description of COVIS

COVIS (Cabled Observatory Vent Imaging Sonar) employs a modified Reson 7125
SeaBat multibeam sonar to image and quantify hydrothermal discharge.  The COVIS
sonar transducers are mounted atop a 4.2-m tower and have three independent
degrees of freedom (pan, tilt, roll). The electronics and data-handling hardware
are contained in two pressure vessels near the base of the tower. COVIS is
specifically designed to (a) acquire water-column multi-beam backscatter data
for imaging hydrothermal plumes and measuring their vertical velocity and heat
flux and (b) acquire acoustic backscatter from the seafloor covered by diffuse
hydrothermal discharge.

COVIS makes three types of sonar measurements: (1) **Imaging,** which provides
3D images of volume scattering strength with relatively little processing, (2)
**Doppler,** which uses the full waveform data to extract fluid line-of-sight
velocities, and (3) **Diffuse-flow imaging** in which the change in seafloor
echos from separate transmission is used to map and quantify diffuse-flow
activity.   

The Reson Seabat has a single EM7216 dual-frequency receive array,
complemented by two projectors:  the TC2162 wide-beam projector operates at both
200 and 400 kHz, while the TC2160 fan-beam projector operates only at 400kHz.
The receiver array, with digital beamforming, provides 128 beams with a 3-dB
azimuthal width of 1.0º at 200 kHz and 256 beams of 0.5º width at 400 kHz. In
both cases, a 128º sector is covered. If necessary, the sensor head on COVIS can
be panned to provide a total coverage of up to 256º of azimuth by taking
multiple 128º sector scans in sequence.  The system’s total range of motion in
yaw is limited by the system cabling.

For the diffuse-flow measurements, COVIS
operates at 200 kHz with the TC2162 projector covering a wide azimuthal sector
(130º 3-dB beamwidth) with 20º 3-dB beamwidth in elevation. The sonar head is
held in a fixed, negative (tilted downward) elevation.    Both Doppler and
Imaging measurements (which target the plumes rising above focused discharges)
use the TC2160 projector at 400 kHz with a transmit beam that is wide in the
azimuthal coordinate (130º 3-dB beamwidth) and narrow in the elevation
coordinate (1º 3-dB beamwidth).  In these cases, COVIS’ sensor array is swept
through a range of elevations in 1º increments. Overall, the Imaging
measurements use a similar setup to the Doppler, but with a narrower pulse and
fewer pings per elevation to get a higher resolution image.

In the previous deployment at Grotto vent on the Ocean Network Canada Neptune
cabled observatory, a single data collection sequence included one set of
each of the three modes (imaging, Doppler, diffuse) with the overall sequence
requiring approximately 1 hour of active sonar pinging. COVIS ran this complete
sequence once every three hours (8 times / day), with the minimum interval
between sequences determined by the time required for packaging and compressing
the sonar data for upload.  A similar schedule of 8 acquisitions per day is
planned for the OOI deployment.

Under normal operation, COVIS follows a pre-programmed schedule performing sonar
collections as defined in a configuration file, then packaging the resulting
sonar data files for upload.   After each sonar data collection, the control
computer retrieves the sonar files
from the sonar, packages them (see archive file description below), and stores
them to an onboard disk for upload.  Based on the data file sizes and
schedule while deployed on Neptune, the control computer has storage space for
at most one month of data collection; this permits the off-loading of data at a
different schedule than the sonar operations schedule.

## COVIS Data Format

Terminology:

 * **Mode:**  One mode of sonar operation: plume imaging, plume doppler or diffuse flow.
 * **Ping:** One sonar transmission.
 * **Scan:** Set of pings are one elevation angle.
 * **Sweep:**  Complete sweep across a range of elevations, comprised of N scans (except for diffuse flow, where elevation is held fixed).
 * **Run:**  Complete data acquisition in one mode, comprised of 1 to N sweeps.
 * **Sequence:**  Complete set of runs spanning all modes.


For each sonar run (one period of collection of one sonar mode), the
data format as seen by the end user is a compressed archive
(in the Unix-standard .tar.gz format) which contains the raw Reson IQ data
files and metadata about rotator position during the sonar collection and
other ping metadata.    This archive is generated onboard COVIS by the SIC
at the conclusion of the sweep and stored on disk in the SIC, where it can
then be retrieved over the OOI network using standard network file transfer
tools (e.g., rsync, scp, ftp, etc).    As each sonar archive file corresponds
to a single sonar mode, for the “standard” 8 sequences daily COVIS sampling
routine used on Neptune, a total of 24 archive files are produced per day.  
As shown in Table 2, the archive file size varies greatly between sonar
operating modes.

Each archive file contains the following files

  1. *sweep.json*, a json-formatted file that for each sweep describes start
   time, mode, and possibly additional annotations;
  1. *transducer.xml* gives a list of configuration values from the Reson
    control software.
  1.  *index.csv*, which contains a list of all the pings and the
  corresponding sonar orientations;
  1.  a number of record files (*rec_7038_xxxxxx.bin* and *rec_7000_xxxxxx.json*)
  equal to the number of pings (where “xxxxxx” is the zero-padded ping number).  
  The *rec_7038_xxxxxx.bin* files are the binary sonar data (see Reson
  documentation), and the *rec_7000_xxxxxx.json* files contain metadata
  describing the ping and the sonar.  The number of pings varies depending
  on the mode, sonar operations (the setup is by time not ping due to
  specifics of Reson control commands), and other operations parameters.

Each file is described in more detail below:

### sweep.json

A sample `sweep.json` is shown below:

    {
        "mode": "imageleft",
        "_id": "COVIS-20181109T030002-imageleft",
        "starttime": [
            1541732527,
            313958
        ],
        "endtime": [
            1541733603,
            73036
        ],
        "motion": {
            "start": [
                25,
                0,
                70
            ],
            "inc": [
                0,
                0,
                1
            ],
            "steps": 105,
            "roundtrip": false,
            "mask": 2
        },
        "settings": {
            "rxgain": 53,
            "txpower": 220,
            "pulse_width": 0.0005,
            "range": 75
        }
    }



### index.csv

The `index.csv` file lists the ping number, Unix system time (in seconds and microseconds) for each ping, and the sonar attitude during the ping.  The attitude is given in two triplets:

  * pitch, roll, yaw:
  * kPAngle,kRAngle,kHeading:

The data is stored in CSV format with a single header line:

    ping,seconds,usecs,pitch,roll,yaw,kPAngle,kRAngle,kHeading
    1,1541732528,233184,-1.0,25.1,70.6,7.4,-10.9,216.7
    2,1541732528,688452,-1.0,25.1,70.6,7.4,-10.9,216.7
    3,1541732529,188446,-1.0,25.1,70.6,7.4,-10.9,216.7
    4,1541732529,688427,-1.0,25.1,70.6,7.4,-10.9,216.7
    6,1541732534,887476,-1.0,25.1,71.0,7.4,-11.0,216.6
    7,1541732535,387463,-1.0,25.1,71.0,7.4,-11.0,216.6
    8,1541732535,887448,-1.0,25.1,71.0,7.4,-11.0,216.6

### transducer.xml

Transducer.xml is an XML configuration file generated by the Reson sonar.  It is
generally not required for post-processing and is not described here.

### sonar data files

The sonar data itself is stored as pairs of files:

  * `rec_7000_XXXXXX.json` files contain a ping description, with XXXXXX being the ping number
  * `rec_7038_XXXXXX.bin` files contain binary Reson sonar data.

A sample `.json` file contains:

    {
        "hdr": {
            "sonar_id": 0,
            "ping_num": 1,
            "multi_ping": 0,
            "xmit_freq": 396000,
            "sample_rate": 34482.758,
            "rcvr_bandwidth": 0,
            "pulse_width": 0.0005,
            "pulse_type": 0,
            "envelope_type": 0,
            "envelope_param": 0,
            "pulse_extra": 0,
            "max_ping_rate": 2,
            "ping_period": 1541732500,
            "range_sel": 75,
            "power_sel": 220,
            "gain_sel": 53,
            "control": {
                "auto_bd_filter_method": 0,
                "auto_gain_method": 0,
                "auto_range_method": 0,
                "bd_depth": false,
                "bd_range": false
            },
            "prj": {
                "id": 1,
                "vert_angle": 0,
                "horiz_angle": 0,
                "vert_width": 0.0174533,
                "horiz_width": 2.0943952,
                "focal_point": 10000000,
                "window_type": 0,
                "window_param": 0
            },
            "transmit": {
                "pitch_stabilization": 0,
                "yaw_stabilization": 0
            },
            "hydrophone_id": 0,
            "recv": {
                "window_type": 0,
                "window_param": 0,
                "flags": 1,
                "beam_width": 0.008
            },
            "bd": {
                "min_range": 1,
                "max_range": 5,
                "min_depth": 1,
                "max_depth": 5
            },
            "absorption": 0,
            "sound_speed": 1468,
            "spreading_loss": 0
      }
    }
