EVAL_DF_BINS=10,20,30,40,50,60,70,80,90,100

all: eval lighteval

lighteval: vitoria.lighteval stuttgart.lighteval paris.lighteval switzerland.lighteval

eval: vitoria.eval stuttgart.eval paris.eval switzerland.eval

clean:
	@rm -f *.eval
	@rm -rf gtfs
	@rm -rf osm
	@rm -rf evalout

osmconvert:
	@echo `date +"[%F %T.%3N]"` "EVAL : Fetching osmconvert..."
	@curl http://m.m.i24.cc/osmconvert.c | cc -x c - -lz -O3 -o osmconvert

%.lighteval: osm/%.osm gtfs/%/stops.txt gtfs/%/stop_times.txt gtfs/%/trips.txt gtfs/%/routes.txt eval.cfg
	@echo `date +"[%F %T.%3N]"` "EVAL : Running light (without stats) evaluation for '"$*"'..."
	@mkdir -p gtfs/$*/shaped
	@rm -f gtfs/$*/shaped/*
	@../build/pfaedle -x $< -i gtfs/$* -c eval.cfg -o gtfs/$*/shaped -D -m all 2>&1 | tee $@

%.eval: osm/%.osm gtfs/%/stops.txt gtfs/%/stop_times.txt gtfs/%/trips.txt gtfs/%/routes.txt eval.cfg eval-wo-osm-line-rels.cfg
	@echo `date +"[%F %T.%3N]"` "EVAL : Running evaluation for '"$*"'..."
	@mkdir -p gtfs/$*/shaped
	@rm -f gtfs/$*/shaped/*
	@mkdir -p evalout/
	@mkdir -p evalout/$*/
	@mkdir -p evalout/$*/hmm+osm
	@../build/pfaedle -x $< -i gtfs/$* -c eval.cfg --eval-path evalout/$*/hmm+osm -o gtfs/$*/shaped -D -m all --eval --eval-df-bins $(EVAL_DF_BINS) 2>&1 | tee $@
	@find evalout/$*/hmm+osm/ -name "*.json" -print0 | xargs -0 rm

	@mkdir -p evalout/$*/greedy
	@../build/pfaedle -x $< -i gtfs/$* -c eval.cfg --method greedy --eval-path evalout/$*/greedy -o gtfs/$*/shaped -D -m all --eval --eval-df-bins $(EVAL_DF_BINS) 2>&1 | tee $@
	@find evalout/$*/greedy/ -name "*.json" -print0 | xargs -0 rm

	@mkdir -p evalout/$*/greedy2
	@../build/pfaedle -x $< -i gtfs/$* -c eval.cfg --method greedy2 --eval-path evalout/$*/greedy2 -o gtfs/$*/shaped -D -m all --eval --eval-df-bins $(EVAL_DF_BINS) 2>&1 | tee $@
	@find evalout/$*/greedy2/ -name "*.json" -print0 | xargs -0 rm

	@mkdir -p evalout/$*/hmm
	@../build/pfaedle -x $< -i gtfs/$* -c eval-wo-osm-line-rels.cfg --eval-path evalout/$*/hmm -o gtfs/$*/shaped -D -m all --eval --eval-df-bins $(EVAL_DF_BINS) 2>&1 | tee $@
	@find evalout/$*/hmm/ -name "*.json" -print0 | xargs -0 rm

osm/spain-latest.osm.pbf:
	@mkdir -p osm
	@echo `date +"[%F %T.%3N]"` "EVAL : Downloading OSM data for Spain..."
	@curl --progress-bar http://download.geofabrik.de/europe/spain-latest.osm.pbf > $@

osm/spain-latest.osm: osm/spain-latest.osm.pbf osmconvert
	@# pre-filter to vitoria gasteiz
	@echo `date +"[%F %T.%3N]"` "EVAL : Pre-filtering OSM data to Vitoria-Gasteiz..."
	@osmconvert -b=-2.8661,42.7480,-2.4788,43.0237 $< > $@

osm/baden-wuerttemberg-latest.osm.pbf:
	@mkdir -p osm
	@echo `date +"[%F %T.%3N]"` "EVAL : Downloading OSM data for Baden-Württemberg..."
	@curl --progress-bar http://download.geofabrik.de/europe/germany/baden-wuerttemberg-latest.osm.pbf > $@

osm/baden-wuerttemberg-latest.osm: osm/baden-wuerttemberg-latest.osm.pbf osmconvert
	@echo `date +"[%F %T.%3N]"` "EVAL : Extracting OSM data..."
	@osmconvert $< > $@

osm/france-latest.osm.pbf:
	@mkdir -p osm
	@echo `date +"[%F %T.%3N]"` "EVAL : Downloading OSM data for France..."
	@curl --progress-bar http://download.geofabrik.de/europe/france-latest.osm.pbf > $@

osm/paris-latest.osm: osm/france-latest.osm.pbf osmconvert
	@# pre-filter to greater ile de france
	@echo `date +"[%F %T.%3N]"` "EVAL : Pre-filtering OSM data to Île-de-France..."
	@osmconvert -b=0.374,47.651,4.241,50.261 $< > $@

osm/europe-latest.osm.pbf:
	@mkdir -p osm
	@echo `date +"[%F %T.%3N]"` "EVAL : Downloading OSM data for Europe..."
	@curl --progress-bar http://download.geofabrik.de/europe-latest.osm.pbf > $@

osm/switzerland-latest.osm: osm/europe-latest.osm.pbf osmconvert
	@# pre-filter to greater switzerland
	@echo `date +"[%F %T.%3N]"` "EVAL : Pre-filtering OSM data to Switzerland..."
	@osmconvert -b=3.757,44.245,15.579,52.670 $< > $@

gtfs/vitoria/%.txt:
	@echo `date +"[%F %T.%3N]"` "EVAL : Downloading GTFS data for Vitoria-Gasteiz..."
	@mkdir -p gtfs
	@mkdir -p gtfs/vitoria
	@curl --progress-bar https://transitfeeds.com/p/tuvisa-euskotran/239/latest/download > gtfs/vitoria/gtfs.zip
	@cd gtfs/vitoria && unzip -qq -o gtfs.zip
	@rm gtfs/vitoria/gtfs.zip

gtfs/stuttgart/%.txt:
	@echo `date +"[%F %T.%3N]"` "EVAL : Downloading GTFS data for Stuttgart..."
	@mkdir -p gtfs
	@mkdir -p gtfs/stuttgart
	@echo "******************************************************************"
	@echo "* A password is required to access the VVS dataset. Send a mail  *"
	@echo "* to brosi@cs.informatik.uni-freiburg.de to receive the password. "
	@echo "******************************************************************"
	@curl --progress-bar http://www.vvs.de/download/opendata/VVS_GTFS.zip -su vvsopendata01 > gtfs/stuttgart/gtfs.zip
	@cd gtfs/stuttgart && unzip -qq -o gtfs.zip
	@rm gtfs/stuttgart/gtfs.zip

gtfs/paris/%.txt:
	@echo `date +"[%F %T.%3N]"` "EVAL : Downloading GTFS data for Paris..."
	@mkdir -p gtfs
	@mkdir -p gtfs/paris
	@curl --progress-bar https://transitfeeds.com/p/stif/822/latest/download > gtfs/paris/gtfs.zip
	@cd gtfs/paris && unzip -qq -o gtfs.zip
	@rm gtfs/paris/gtfs.zip

gtfs/switzerland/%.txt:
	@echo `date +"[%F %T.%3N]"` "EVAL : Downloading GTFS data for Switzerland..."
	@mkdir -p gtfs
	@mkdir -p gtfs/switzerland
	@curl --progress-bar http://gtfs.geops.ch/dl/gtfs_complete.zip > gtfs/switzerland/gtfs.zip
	@cd gtfs/switzerland && unzip -qq -o gtfs.zip
	@rm gtfs/switzerland/gtfs.zip


osm/vitoria.osm: osm/spain-latest.osm gtfs/vitoria/stops.txt gtfs/vitoria/trips.txt gtfs/vitoria/routes.txt gtfs/vitoria/stop_times.txt eval.cfg
	@../build/pfaedle -x $< -i gtfs/vitoria/ -c eval.cfg -m all -X $@

osm/stuttgart.osm: osm/baden-wuerttemberg-latest.osm gtfs/stuttgart/stops.txt gtfs/stuttgart/trips.txt gtfs/stuttgart/routes.txt gtfs/stuttgart/stop_times.txt eval.cfg
	@../build/pfaedle -x $< -i gtfs/stuttgart/ -c eval.cfg -m all -X $@

osm/paris.osm: osm/paris-latest.osm gtfs/paris/stops.txt gtfs/paris/trips.txt gtfs/paris/routes.txt gtfs/paris/stop_times.txt eval.cfg
	@../build/pfaedle -x $< -i gtfs/paris/ -c eval.cfg -m all -X $@

osm/switzerland.osm: osm/switzerland-latest.osm gtfs/switzerland/stops.txt eval.cfg
	@../build/pfaedle -x $< -i gtfs/switzerland/ -c eval.cfg -m all -X $@
