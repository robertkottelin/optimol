# Docker Image Update Procedure for Molecular Optimization Backend

## Environment Specifications

- **Deployment Architecture**: Docker Compose
- **Backend Image**: robertkottelin/optimize-molecule:backend-latest
- **Location**: `/opt/optimize-molecule/`
- **Containers**: optimize-molecule (backend), nginx-proxy (reverse proxy)

## 1. Image Update Workflow

### 1.1. Manual Update Procedure

```bash
# SSH into droplet
ssh root@64.227.122.193

# Navigate to application directory
cd /opt/optimize-molecule

# prune af if needed 
docker image prune -af

# Pull latest images from Docker Hub
docker-compose pull

# View available images (verify new image exists)
docker images

# Stop and remove existing containers
docker-compose down

# Start containers with updated images
docker-compose up -d

# Verify container status
docker-compose ps

# Check container logs for errors
docker-compose logs --tail=50
```

### 1.2. Verification Procedure

```bash
# Verify API health endpoint
curl -I https://optimizemolecule.com/health

# Check container operational metrics
docker stats --no-stream

# Verify container is using updated image
docker inspect optimize-molecule | grep Image
```

## 2. Targeted Single-Service Updates

For updating only the backend while maintaining nginx configuration:

```bash
# Pull only the backend image
docker-compose pull backend

# Restart only backend service
docker-compose up -d --no-deps backend

# Verify backend status
docker ps --filter name=optimize-molecule
```

## 3. Rollback Procedure

If issues are detected with the updated image:

```bash
# List available images
docker images robertkottelin/optimize-molecule

# Modify docker-compose.yml to use specific version
sed -i 's|robertkottelin/optimize-molecule:backend-latest|robertkottelin/optimize-molecule:backend-previous|g' docker-compose.yml

# Restart with previous version
docker-compose up -d

# Verify rollback was successful
docker inspect optimize-molecule | grep Image
```

## 4. Automated Update Implementation

### 4.1. Watchtower Configuration

Install and configure Watchtower for automated updates:

```bash
# Create watchtower configuration directory
mkdir -p /opt/watchtower

# Create docker-compose.yml for watchtower
cat > /opt/watchtower/docker-compose.yml << 'EOF'
version: '3'

services:
  watchtower:
    image: containrrr/watchtower
    container_name: watchtower
    volumes:
      - /var/run/docker.sock:/var/run/docker.sock
    command: --interval 86400 --cleanup optimize-molecule
    restart: unless-stopped
EOF

# Start watchtower
cd /opt/watchtower
docker-compose up -d

# Verify watchtower is running
docker ps | grep watchtower
```

### 4.2. Webhook-based Updates

Configure webhook receiver for CI/CD pipeline integration:

```bash
# Install webhook
apt update && apt install -y webhook

# Create webhook configuration
mkdir -p /etc/webhook
cat > /etc/webhook/hooks.json << 'EOF'
[
  {
    "id": "deploy-backend",
    "execute-command": "/opt/scripts/update-backend.sh",
    "command-working-directory": "/opt/optimize-molecule",
    "trigger-rule": {
      "match": {
        "type": "payload-hash-sha1",
        "secret": "YOUR_SECRET_TOKEN",
        "parameter": {
          "source": "header",
          "name": "X-Hub-Signature"
        }
      }
    }
  }
]
EOF

# Create update script
mkdir -p /opt/scripts
cat > /opt/scripts/update-backend.sh << 'EOF'
#!/bin/bash
set -e

cd /opt/optimize-molecule
docker-compose pull backend
docker-compose up -d --no-deps backend

# Verify health endpoint
sleep 5
curl -s https://optimizemolecule.com/health | grep -q '"status":"healthy"'
exit_code=$?

# Log update result
if [ $exit_code -eq 0 ]; then
  echo "$(date): Backend update successful" >> /var/log/backend-updates.log
else
  echo "$(date): Backend update failed, health check returned non-zero" >> /var/log/backend-updates.log
fi

exit $exit_code
EOF

# Set executable permissions
chmod +x /opt/scripts/update-backend.sh

# Configure webhook service
cat > /etc/systemd/system/webhook.service << 'EOF'
[Unit]
Description=Webhook Service
After=network.target

[Service]
ExecStart=/usr/bin/webhook -hooks /etc/webhook/hooks.json -verbose
Restart=always
User=root
Group=root

[Install]
WantedBy=multi-user.target
EOF

# Enable and start webhook service
systemctl daemon-reload
systemctl enable webhook
systemctl start webhook

# Open firewall for webhook (default port 9000)
ufw allow 9000/tcp
```

## 5. Blue-Green Deployment Configuration

For zero-downtime updates:

```bash
# Create blue-green deployment script
cat > /opt/scripts/blue-green-deploy.sh << 'EOF'
#!/bin/bash
set -e

# Pull latest image
docker pull robertkottelin/optimize-molecule:backend-latest

# Determine current active container (blue or green)
if docker ps --filter name=optimize-molecule-blue | grep -q optimize-molecule-blue; then
  ACTIVE="blue"
  INACTIVE="green"
else
  ACTIVE="green"
  INACTIVE="blue"
fi

echo "Current active container: optimize-molecule-$ACTIVE"
echo "Deploying to: optimize-molecule-$INACTIVE"

# Start inactive container with new image
docker run -d \
  --name optimize-molecule-$INACTIVE \
  -p 500$([[ "$INACTIVE" == "blue" ]] && echo "0" || echo "1")":5000 \
  --restart unless-stopped \
  robertkottelin/optimize-molecule:backend-latest

# Wait for container to be healthy
echo "Waiting for container to initialize..."
sleep 10

# Test health of new container
HEALTH_PORT=$([[ "$INACTIVE" == "blue" ]] && echo "5000" || echo "5001")
if curl -s http://localhost:$HEALTH_PORT/health | grep -q '"status":"healthy"'; then
  echo "New container is healthy, updating Nginx configuration..."
  
  # Update Nginx config to point to new container
  sed -i "s/backend:5000/optimize-molecule-$INACTIVE:5000/g" /opt/optimize-molecule/nginx/conf/default.conf
  
  # Reload Nginx
  docker exec nginx-proxy nginx -s reload
  
  # Stop and remove old container
  docker stop optimize-molecule-$ACTIVE
  docker rm optimize-molecule-$ACTIVE
  
  echo "Deployment complete. Active container is now: optimize-molecule-$INACTIVE"
else
  echo "Health check failed for new container. Aborting deployment."
  docker stop optimize-molecule-$INACTIVE
  docker rm optimize-molecule-$INACTIVE
  exit 1
fi
EOF

# Set executable permissions
chmod +x /opt/scripts/blue-green-deploy.sh
```

## 6. Security Considerations

1. **Image Scanning**:
   ```bash
   # Install trivy scanner
   apt install -y wget apt-transport-https gnupg lsb-release
   wget -qO - https://aquasecurity.github.io/trivy-repo/deb/public.key | apt-key add -
   echo deb https://aquasecurity.github.io/trivy-repo/deb $(lsb_release -sc) main | tee -a /etc/apt/sources.list.d/trivy.list
   apt update && apt install -y trivy
   
   # Scan image for vulnerabilities
   trivy image robertkottelin/optimize-molecule:backend-latest
   ```

2. **Image Integrity Verification**:
   ```bash
   # Pull image with digest for immutability
   docker pull robertkottelin/optimize-molecule@sha256:f3beb287b79bfadce0ea07d903d1c356e7eca70680f5437a6c6beea1d32d6140
   ```

## 7. Monitoring and Alerting

Set up update notifications:

```bash
# Create notification script for successful updates
cat > /opt/scripts/update-notification.sh << 'EOF'
#!/bin/bash

# Send notification via curl to webhook/Slack/Discord
curl -X POST \
  -H "Content-Type: application/json" \
  -d "{\"text\":\"Backend image updated successfully on $(hostname) at $(date)\"}" \
  https://hooks.slack.com/services/YOUR_WEBHOOK_URL
EOF

# Append notification to update script
echo "/opt/scripts/update-notification.sh" >> /opt/scripts/update-backend.sh
chmod +x /opt/scripts/update-notification.sh
```