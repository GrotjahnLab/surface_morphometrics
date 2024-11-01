#!/bin/bash
set -e


# Start container
echo "Starting container..."
docker-compose -f docker-compose.yml up -d

# Wait for container to be ready
echo "Waiting for container..."
sleep 5

# Get the container ID from docker-compose
CONTAINER_NAME=$(docker-compose -f docker-compose.yml ps -q surface_morphometrics)

# Enter container if it exists
if [ -n "$CONTAINER_NAME" ]; then
    echo "Entering container..."
    docker exec -it $CONTAINER_NAME bash
else
    echo "Error: Container not found."
    exit 1
fi