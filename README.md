[![Left: station-to-station path of a single train through Switzerlanad, as contained in official timetable data. Right: path of the same train map-matched by pfaedle.](geo/schweiz_ex_res.png?raw=true)](geo/schweiz_ex.png?raw=true)
[![Left: station-to-station path of a single bus through Stuttgart, as contained in official timetable data. Right: path of the same bus map-matched by pfaedle.](geo/stuttgart_ex_res.png?raw=true)](geo/stuttgart_ex.png?raw=true)

[![Build
Status](https://travis-ci.org/ad-freiburg/pfaedle.svg?branch=master)](https://travis-ci.org/ad-freiburg/pfaedle)

# pfaedle

Precise map-matching for public transit schedules (GTFS data).

## Requirements

 * `cmake`
 * `gcc >= 4.8` (or `clang >= 5.0`)

## Building and Installation

Fetch this repository and init submodules:

```
git clone --recurse-submodules https://github.com/ad-freiburg/pfaedle
```

```
mkdir build && cd build
cmake ..
make -j
```

To install, type
```
make install
```

# General Usage

## Generating shapes for a GTFS feed

```
pfaedle -c <CFG FILE> -x <OSM FILE> <GTFS INPUT FOLDER>
```

A shape'd version of the input GTFS feed will be written to `./gtfs-out`.

By default, shapes are only calculated for trips that don't have a shape in the
input feed. To drop all existing shapes, use the `-D` flag.

For example, you may generate (and replace existing, see -D parameter) shapes for the GTFS dataset for Freiburg like this:

```
$ mkdir freiburg_gtfs && cd freiburg_gtfs
$ wget https://fritz.freiburg.de/csv_Downloads/VAGFR.zip
$ unzip VAGFR.zip
$ wget http://download.geofabrik.de/europe/germany/baden-wuerttemberg/freiburg-regbez-latest.osm.bz2
$ bunzip2 freiburg-regbez-latest.osm.bz2
$ mkdir gtfs-out
$ pfaedle -D -c pfaedle.cfg -x freiburg-regbez-latest.osm .
```

A default configuration file `pfaedle.cfg` can be found in this repo.


## Generating shapes for a specific MOT

To generate shapes only for a specific mot, use the `-m` option. Possible
values are either `tram`, `bus`, `rail`, `subway`, `ferry`, `funicular`,
`gondola`, `all` (default) or GTFS vehicle type codes (0, 1, 2, 3, 4, 5, 6, 7).

Multiple values can be specified (comma separated).

## OSM filtering

`pfaedle` comes with the ability to filter OpenStreetMap data. If you specify
the `-X` flag, `pfaedle` will filter the input OSM file and output a new OSM
file which contains *exactly* the data needed to calculate the shapes for the
input GTFS feed and the input configuration.

This can be used to avoid parsing (for example) the entire world.osm on each
run.

## Debugging

The following flags may be useful for debugging:

 * `-T <GTFS TRIP ID>` only calculate shape for a single trip (specified via its GTFS trip id) and output it as GeoJSON to
   `<dbg-path>/path.json`
 * `--write-graph` write the graph used for routing as GeoJSON to
   `<dbg-path>/graph.json`
 * `--write-cgraph` if `-T` is set, write the combination graph used for
   routing as GeoJSON to `<dbg-path>/combgraph.json`
 * `--write-trgraph` write the complete network graph to `<dbg-path>/trgraph.json`

# Configuration

The main config file distributed with this repository is `pfaedle.cfg`. The
config file has some comments which hopefully explain the meaning behind the
parameters.

# Evaluation

You may run an entire evaluation of our testing datasets Vitoria-Gasteiz, Paris, Switzerland and
Stuttgart with

```
mkdir build && cd build
cmake ..
make -j
make eval
```

*Notes:*
 * this will download, and filter, the entire OSM files for Spain, the
Stuttgart region. Make sure you have enough space left on your hard drive.
 * in evaluation mode, pfaedle needs significantly more time, because the
   calculation of the similarity measurements between shapes are expensive
 * if you are only interested in the end results of a single dataset, run
   `make <dataset>.lighteval` in `/eval`. For example, `make paris.lighteval`
   generates a shaped version of the paris dataset, without doing extensive
   comparisons to the ground truth.
 * similarily, if you want to run the extensive evaluation for a single dataset,
   run `make <dataset>.eval` in `/eval`.


## Evaluation requirements

 * zlib

On Debianesque systems, type

```
sudo apt-get install zlib1g-dev
```

to install the dependencies.
