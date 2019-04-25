#!/bin/bash
echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin
docker build -t mgnify/fetch_tool .
docker push mgnify/fetch_tool:latest