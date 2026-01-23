#!/bin/bash

PORT=3838

echo "🔍 Checking if Docker is using port $PORT..."

CONTAINER_ID=$(docker ps --filter "publish=$PORT" -q)

if [ -n "$CONTAINER_ID" ]; then
  echo "⚠️  Port $PORT is used by container: $CONTAINER_ID"
  echo "🛑 Stopping container..."
  docker stop "$CONTAINER_ID"
  echo "🧹 Removing container..."
  docker rm "$CONTAINER_ID"
else
  echo "✅ No Docker container is using port $PORT"
fi

echo "🚀 You can now safely start your container"

cp -r /home/os/thermochemicalDenaturationApp/* /home/os/eSPC_biophysics_platform/thermochemicalDenaturationApp/
docker build -t chemelt -f ./dockerFiles/Dockerfile_thermoChemical .
docker run -p 3838:3838 chemelt
