#!/bin/bash

# Exit on error
set -e

# Container name/ID to stop 
CONTAINER_NAME="surface_morphometrics_container"

# Timeout in seconds
TIMEOUT=5

# Handle Ctrl+C
trap 'echo "Received interrupt signal..."; exit 1' INT

echo "Stopping container $CONTAINER_NAME..."

# Check if container is running and stop it
if docker ps -q --filter "name=$CONTAINER_NAME" | grep -q .; then
    docker stop -t "$TIMEOUT" "$CONTAINER_NAME" || {
        echo "Failed to stop container $CONTAINER_NAME"
        exit 1
    }
    echo "Container stopped successfully"
else
    echo "Container $CONTAINER_NAME is not running"
    exit 0
fi