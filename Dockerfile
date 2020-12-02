FROM debian:buster-slim AS builder

WORKDIR /app

RUN apt-get update && \
	apt-get install -y g++ cmake git

ADD . /app
RUN mkdir build && \
	cd build && \
	cmake .. && \
	make -j && \
	pwd && \
	make install

FROM debian:buster-slim

RUN apt-get update && \
	apt-get install -y libgomp1 && \
	rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/etc/pfaedle /usr/local/etc/pfaedle
COPY --from=builder /usr/local/bin/pfaedle /usr/local/bin/pfaedle

ENTRYPOINT ["/usr/local/bin/pfaedle"]
