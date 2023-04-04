[![Left: station-to-station path of a single train through Switzerland obtained from schedule timetable data. Right: path of the same train map-matched by pfaedle.](geo/schweiz_ex_res.png?raw=true)](geo/schweiz_ex.png?raw=true)
*Left: station-to-station path of a single train through Switzerland obtained from official schedule data. Right: path of the same train map-matched by pfaedle.*

[![Left: station-to-station path of a single bus through Stuttgart obtained from official schedule data. Right: path of the same bus map-matched by pfaedle.](geo/stuttgart_ex_res.png?raw=true)](geo/stuttgart_ex.png?raw=true)
*Left: station-to-station path of a single bus through Stuttgart obtained from official schedule data. Right: path of the same bus map-matched by pfaedle.*

[![Build](https://github.com/ad-freiburg/pfaedle/actions/workflows/build.yml/badge.svg)](https://github.com/ad-freiburg/pfaedle/actions/workflows/build.yml)

# pfaedle

Precise OpenStreetMap (OSM) map-matching for public transit schedules ([GTFS](https://developers.google.com/transit/gtfs/reference/) data).
First described in [this 2018 SIGSPATIAL paper](http://ad-publications.informatik.uni-freiburg.de/SIGSPATIAL_Sparse%20map%20matching%202018.pdf).

For a quick visual inspection of the shape quality, see for example the schedule data for Germany or Switzerland in our tool [TRAVIC](https://travic.app/?z=7&x=1261608.6&y=6430601.6).

## Requirements

 * `cmake`
 * `gcc >= 5.0` (or `clang >= 3.9`)
 * `libzip` (optional, for ZIP support)

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
pfaedle -x <OSM FILE> <GTFS INPUT FEED>
```

A shape'd version of the input GTFS feed will be written to `./gtfs-out`.

By default, shapes are only calculated for trips that don't have a shape in the
input feed. To drop all existing shapes, use the `-D` flag.

For example, you may generate (and replace existing, see -D parameter) shapes for the GTFS dataset for Freiburg like this:

```
$ wget https://fritz.freiburg.de/csv_Downloads/VAGFR.zip
$ wget http://download.geofabrik.de/europe/germany/baden-wuerttemberg/freiburg-regbez-latest.osm.bz2 && bunzip2 freiburg-regbez-latest.osm.bz2
$ pfaedle -D -x freiburg-regbez-latest.osm VAGFR.zip
```

## Generating shapes for a specific MOT

To generate shapes for a specific mot only, use the `-m` option. Possible
values are either `tram`, `bus`, `coach`, `rail`, `subway`, `ferry`, `funicular`,
`gondola`, `all` (default) or GTFS route type codes (0, 1, 2, 3, 4, 5, 6, 7, or [extended route types](https://developers.google.com/transit/gtfs/reference/extended-route-types)). Integer codes will only match the specific route type, while string codes will match classes of route types. For example, `-m 101` will only match routes with `route_type` `101` (high speed rail), while `-m rail` will match any rail service encoded via a standard `route_type` `2` or an extended `route_type` describing a rail service (e.g. `100`, `101`, `102`, ...).

Multiple values can be specified (comma separated).

## OSM filtering

`pfaedle` comes with the ability to filter OpenStreetMap data. If you specify
the `-X` flag, `pfaedle` will filter the input OSM file and output a new OSM
file which contains exactly the data needed to calculate the shapes for the
input GTFS feed and the input configuration.

This can be used to avoid parsing (for example) the entire `planet.osm` on each
run.

## via Docker

You can use the [`adfreiburg/pfaedle` Docker image](https://hub.docker.com/r/adfreiburg/pfaedle) by mounting the OSM & GTFS data into the container:

```shell
docker run -i --rm \
	# mount OSM data
	--volume /path/to/osm/data:/osm \
	# mount GTFS data
	--volume /path/to/gtfs/data:/gtfs \
	# tell pfaedle where to find the data
	adfreiburg/pfaedle -x /osm/osm-data.xml -i /gtfs
```

## Debugging

The following flags may be useful for debugging:

 * `-T <GTFS TRIP ID>` only calculate shape for a single trip (specified via its GTFS trip id) and output it as GeoJSON to
   `<dbg-path>/path.json`
 * `--write-graph` write the graph used for routing as GeoJSON to
 * `--write-trgraph` write the complete network graph to `<dbg-path>/trgraph.json`

# Configuration

A default configuration file `pfaedle.cfg` can be found in this repo and will be installed with `make install`. Custom configuration files can be specified with the `-c` flag. If no `-c` flag is set, `pfaedle` will parse and merge the following cfg files in the given order (if present): `<install prefix>/etc/pfaedle/pfaedle.cfg`, `$HOME/.config/pfaedle/pfaedle.cfg`, `<CWD>/pfaedle.cfg`. Values given in later files will overwrite earlier defined values.
