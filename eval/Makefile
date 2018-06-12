EVAL_DF_BINS=10,20,30,40,50,60,70,80,90,100

all: eval lighteval

lighteval: vitoria.lighteval stuttgart.lighteval paris.lighteval switzerland.lighteval

eval: vitoria.eval stuttgart.eval paris.eval switzerland.eval

clean:
	rm -f *.eval
	rm -rf gtfs
	rf -rf osm
	rm -rf evalout

osmconvert:
	wget -O - http://m.m.i24.cc/osmconvert.c | cc -x c - -lz -O3 -o osmconvert

%.lighteval: osm/%.osm gtfs/%/stops.txt gtfs/%/stop_times.txt gtfs/%/trips.txt gtfs/%/routes.txt eval.cfg eval-wo-osm.cfg
	mkdir -p gtfs/$*/shaped
	rm -f gtfs/$*/shaped/*
	../build/pfaedle -x $< -i gtfs/$* -c eval.cfg -o gtfs/$*/shaped -D -m all 2>&1 | tee $@

%.eval: osm/%.osm gtfs/%/stops.txt gtfs/%/stop_times.txt gtfs/%/trips.txt gtfs/%/routes.txt eval.cfg eval-wo-osm.cfg
	mkdir -p gtfs/$*/shaped
	rm -f gtfs/$*/shaped/*
	mkdir -p evalout/
	mkdir -p evalout/$*/
	mkdir -p evalout/$*/hmm+osm
	../build/pfaedle -x $< -i gtfs/$* -c eval.cfg --eval-path evalout/$*/hmm+osm -o gtfs/$*/shaped -D -m all --eval --eval-df-bins $(EVAL_DF_BINS) 2>&1 | tee $@
	find evalout/$*/hmm+osm/ -name "*.json" -print0 | xargs -0 rm

	mkdir -p evalout/$*/hmm
	../build/pfaedle -x $< -i gtfs/$* -c eval-wo-osm.cfg --eval-path evalout/$*/hmm -o gtfs/$*/shaped -D -m all --eval --eval-df-bins $(EVAL_DF_BINS) 2>&1 | tee $@
	find evalout/$*/hmm/ -name "*.json" -print0 | xargs -0 rm

	mkdir -p evalout/$*/greedy
	../build/pfaedle -x $< -i gtfs/$* -c eval.cfg --method greedy --eval-path evalout/$*/greedy -o gtfs/$*/shaped -D -m all --eval --eval-df-bins $(EVAL_DF_BINS) 2>&1 | tee $@
	find evalout/$*/greedy/ -name "*.json" -print0 | xargs -0 rm

	mkdir -p evalout/$*/greedy2
	../build/pfaedle -x $< -i gtfs/$* -c eval.cfg --method greedy2 --eval-path evalout/$*/greedy2 -o gtfs/$*/shaped -D -m all --eval --eval-df-bins $(EVAL_DF_BINS) 2>&1 | tee $@
	find evalout/$*/greedy2/ -name "*.json" -print0 | xargs -0 rm

osm/spain-latest.osm.pbf:
	mkdir -p osm
	wget http://download.geofabrik.de/europe/spain-latest.osm.pbf -O $@

osm/spain-latest.osm: osm/spain-latest.osm.pbf osmconvert 
	@# pre-filter to vitoria gasteiz
	osmconvert -b=-2.8661,42.7480,-2.4788,43.0237 $< > $@

osm/baden-wuerttemberg-latest.osm.pbf:
	mkdir -p osm
	wget http://download.geofabrik.de/europe/germany/baden-wuerttemberg-latest.osm.pbf -O $@

osm/baden-wuerttemberg-latest.osm: osm/baden-wuerttemberg-latest.osm.pbf osmconvert 
	osmconvert $< > $@

osm/france-latest.osm.pbf:
	mkdir -p osm
	wget http://download.geofabrik.de/europe/france-latest.osm.pbf -O $@

osm/paris-latest.osm: osm/france-latest.osm.pbf osmconvert 
	@# pre-filter to greater ile de france
	osmconvert -b=0.374,47.651,4.241,50.261 $< > $@

osm/europe-latest.osm.pbf:
	mkdir -p osm
	wget http://download.geofabrik.de/europe-latest.osm.pbf -O $@

osm/switzerland-latest.osm: osm/europe-latest.osm.pbf osmconvert 
	@# pre-filter to greater switzerland
	osmconvert -b=3.757,44.245,15.579,52.670 $< > $@

gtfs/vitoria/%.txt:
	mkdir -p gtfs
	mkdir -p gtfs/vitoria
	wget https://transitfeeds.com/p/tuvisa-euskotran/239/latest/download -O gtfs/vitoria/gtfs.zip
	cd gtfs/vitoria && unzip -o gtfs.zip
	rm gtfs/vitoria/gtfs.zip

gtfs/stuttgart/%.txt:
	mkdir -p gtfs
	mkdir -p gtfs/stuttgart
	echo "******************************************************************"
	echo "* A password is required to access the VVS dataset. Send a mail  *"
	echo "* to brosi@cs.informatik.uni-freiburg.de to receive the password. "
	echo "******************************************************************"
	wget http://www.vvs.de/download/opendata/VVS_GTFS.zip --ask-password --user vvsopendata01 -O gtfs/stuttgart/gtfs.zip
	cd gtfs/stuttgart && unzip gtfs.zip
	rm gtfs/stuttgart/gtfs.zip

gtfs/paris/%.txt:
	mkdir -p gtfs
	mkdir -p gtfs/paris
	wget https://transitfeeds.com/p/stif/822/latest/download -O gtfs/paris/gtfs.zip
	cd gtfs/paris && unzip gtfs.zip
	rm gtfs/paris/gtfs.zip

gtfs/switzerland/%.txt:
	mkdir -p gtfs
	mkdir -p gtfs/switzerland
	wget http://gtfs.geops.ch/dl/gtfs_complete.zip -O gtfs/switzerland/gtfs.zip
	cd gtfs/switzerland && unzip gtfs.zip
	rm gtfs/switzerland/gtfs.zip

osm/vitoria.osm: osm/spain-latest.osm gtfs/vitoria/stops.txt eval.cfg
	../build/pfaedle -x $< -i gtfs/vitoria/ -c eval.cfg -m all -X $@

osm/stuttgart.osm: osm/baden-wuerttemberg-latest.osm gtfs/stuttgart/stops.txt eval.cfg
	../build/pfaedle -x $< -i gtfs/stuttgart/ -c eval.cfg -m all -X $@

osm/paris.osm: osm/paris-latest.osm gtfs/paris/stops.txt eval.cfg
	../build/pfaedle -x $< -i gtfs/paris/ -c eval.cfg -m all -X $@

osm/switzerland.osm: osm/switzerland-latest.osm gtfs/switzerland/stops.txt eval.cfg
	../build/pfaedle -x $< -i gtfs/switzerland/ -c eval.cfg -m all -X $@
