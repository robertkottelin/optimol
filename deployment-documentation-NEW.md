### Deployment Script
```bash
#!/bin/bash

# NEW:

# SSH in to droplet
ssh root@64.227.122.193

# Stop and remove existing containers
docker stop optimol-frontend optimol-backend optimol-redis || true
docker rm optimol-frontend optimol-backend optimol-redis || true

docker image prune -af

# Pull latest images
docker pull robertkottelin/optimize-molecule:frontend-latest
docker pull robertkottelin/optimize-molecule:backend-latest
docker pull redis:7-alpine

# Create Docker network if not exists
docker network create optimol-network || true

# Run Redis container
docker run -d \
  --name optimol-redis \
  --restart unless-stopped \
  --network optimol-network \
  -v optimol_redis_data:/data \
  redis:7-alpine redis-server --appendonly yes

# Run backend container
docker run -d \
  --name optimol-backend \
  --restart unless-stopped \
  --network optimol-network \
  -p 5000:5000 \
  -v /opt/optimol/.env:/app/.env \
  -v optimol_data:/app/instance \
  -e REDIS_URL=redis://optimol-redis:6379/0 \
  robertkottelin/optimize-molecule:backend-latest

# Run frontend container
docker run -d \
  --name optimol-frontend \
  --restart unless-stopped \
  -p 3000:3000 \
  --network optimol-network \
  robertkottelin/optimize-molecule:frontend-latest

# Check logs
docker logs optimol-backend
docker logs optimol-frontend
docker logs optimol-redis

# NEW AUTOMATION
cat > deploy.sh << 'EOF'
#!/bin/bash

# Stop and remove existing containers
docker stop optimol-frontend optimol-backend optimol-redis || true
docker rm optimol-frontend optimol-backend optimol-redis || true

docker image prune -af

# Pull latest images
docker pull robertkottelin/optimize-molecule:frontend-latest
docker pull robertkottelin/optimize-molecule:backend-latest
docker pull redis:7-alpine

# Create Docker network if not exists
docker network create optimol-network || true

# Run Redis container
docker run -d \
  --name optimol-redis \
  --restart unless-stopped \
  --network optimol-network \
  -v optimol_redis_data:/data \
  redis:7-alpine redis-server --appendonly yes

# Run backend container
docker run -d \
  --name optimol-backend \
  --restart unless-stopped \
  --network optimol-network \
  -p 5000:5000 \
  -v /opt/optimol/.env:/app/.env \
  -v optimol_data:/app/instance \
  -e REDIS_URL=redis://optimol-redis:6379/0 \
  robertkottelin/optimize-molecule:backend-latest

# Run frontend container
docker run -d \
  --name optimol-frontend \
  --restart unless-stopped \
  -p 3000:3000 \
  --network optimol-network \
  robertkottelin/optimize-molecule:frontend-latest

# Check logs
docker logs optimol-backend
docker logs optimol-frontend
docker logs optimol-redis
EOF

# make executable:
chmod +x deploy.sh

# run automation script
./deploy.sh

cat > /opt/optimol/.env << 'EOF'
STRIPE_SECRET=""
JWT_SECRET_KEY=""
FLASK_ENV=""
PRICE_ID=""
EOF

# Update Nginx configuration
cat > /etc/nginx/sites-available/optimizemolecule.com << 'EOF'
server {
    listen 80;
    listen 443 ssl;
    server_name optimizemolecule.com www.optimizemolecule.com;

    # SSL configuration
    ssl_certificate /etc/letsencrypt/live/optimizemolecule.com/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/optimizemolecule.com/privkey.pem;
    ssl_session_cache shared:SSL:10m;
    ssl_session_timeout 10m;
    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_prefer_server_ciphers on;
    ssl_ciphers "EECDH+AESGCM:EDH+AESGCM:AES256+EECDH:AES256+EDH";

    # Redirect non-https traffic to https
    if ($scheme != "https") {
        return 301 https://$host$request_uri;
    }

    # Direct API endpoints (must be defined before the catch-all location)
    location = /subscribe {
        proxy_pass http://localhost:5000/subscribe;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    location = /login {
        proxy_pass http://localhost:5000/login;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    location = /register {
        proxy_pass http://localhost:5000/register;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    location = /me {
        proxy_pass http://localhost:5000/me;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    location = /cancel-subscription {
        proxy_pass http://localhost:5000/cancel-subscription;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    location = /health {
        proxy_pass http://localhost:5000/health;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    location = /optimize-molecule {
        proxy_pass http://localhost:5000/optimize-molecule;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # API endpoint with /api prefix
    location /api/ {
        rewrite ^/api(/.*)$ $1 break;
        proxy_pass http://localhost:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # Frontend application (catch-all, must be last)
    location / {
        proxy_pass http://localhost:3000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_cache_bypass $http_upgrade;
        
        # Static asset optimization
        add_header Cache-Control "public, max-age=3600";
        add_header X-Content-Type-Options "nosniff";
    }
}
EOF

# Apply Nginx configuration
nginx -t && systemctl restart nginx

# Test endpoints
curl -s https://optimizemolecule.com/api/health
```
