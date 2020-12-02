#!/bin/sh
set -ex
set -o pipefail

echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin
docker push ad-freiburg/pfaedle