FROM debian:bookworm-slim AS builder

WORKDIR /app

RUN apt-get update && \
	apt-get install -y g++ cmake git libzip-dev zlib1g-dev libbz2-dev

ADD . /app
RUN mkdir build && \
	cd build && \
	cmake .. && \
	make -j && \
	pwd && \
	make install

FROM debian:bookworm-slim

RUN apt-get update && \
	apt-get install -y libzip4 zlib1g libbz2-1.0 && \
	rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/etc/pfaedle /usr/local/etc/pfaedle
COPY --from=builder /usr/local/bin/pfaedle /usr/local/bin/pfaedle

ENTRYPOINT ["/usr/local/bin/pfaedle"]
