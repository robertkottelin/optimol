# Molecular Optimization Backend Deployment Documentation

## Environment Specifications

- **Host**: DigitalOcean Droplet
- **IP**: 134.209.226.93
- **Domain**: optimizemolecule.com
- **Backend Image**: robertkottelin/optimize-molecule:backend-latest
- **Frontend Repository**: https://github.com/robertkottelin/optimol
- **Frontend Host**: GitHub Pages (https://robertkottelin.github.io/optimol/)

## 1. Initial Server Configuration

SSH into the Digital Ocean droplet:

```bash
ssh root@134.209.226.93
```

Install Docker and dependencies:

```bash
apt update
apt install -y apt-transport-https ca-certificates curl software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo "deb [arch=amd64 signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null
apt update
apt install -y docker-ce docker-ce-cli containerd.io
```

Verify Docker installation:

```bash
docker --version
```

## 2. Docker Image Management

Pull the Docker image directly from Docker Hub:

```bash
docker pull robertkottelin/optimize-molecule:backend-latest
```

Test running the container directly:

```bash
docker run -d \
  --name optimize-molecule \
  -p 5000:5000 \
  --restart unless-stopped \
  robertkottelin/optimize-molecule:backend-latest
```

## 3. Initial Nginx Configuration (System-level)

Install and configure Nginx system service:

```bash
apt install -y nginx
```

Generate self-signed certificate for initial testing:

```bash
mkdir -p /etc/nginx/ssl
openssl req -x509 -nodes -days 365 -newkey rsa:2048 \
  -keyout /etc/nginx/ssl/nginx.key \
  -out /etc/nginx/ssl/nginx.crt \
  -subj "/C=US/ST=State/L=City/O=Organization/CN=134.209.226.93"
```

Create Nginx configuration:

```bash
cat > /etc/nginx/sites-available/optimize-molecule << 'EOF'
server {
    listen 80;
    server_name 134.209.226.93;
    
    location / {
        return 301 https://$host$request_uri;
    }
}

server {
    listen 443 ssl;
    server_name 134.209.226.93;
    
    ssl_certificate /etc/nginx/ssl/nginx.crt;
    ssl_certificate_key /etc/nginx/ssl/nginx.key;
    
    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_prefer_server_ciphers on;
    ssl_ciphers ECDHE-RSA-AES256-GCM-SHA512:DHE-RSA-AES256-GCM-SHA512:ECDHE-RSA-AES256-GCM-SHA384:DHE-RSA-AES256-GCM-SHA384;
    
    add_header 'Access-Control-Allow-Origin' 'https://robertkottelin.github.io' always;
    add_header 'Access-Control-Allow-Methods' 'GET, POST, OPTIONS' always;
    add_header 'Access-Control-Allow-Headers' 'Content-Type, Authorization' always;
    
    location / {
        if ($request_method = 'OPTIONS') {
            add_header 'Access-Control-Allow-Origin' 'https://robertkottelin.github.io' always;
            add_header 'Access-Control-Allow-Methods' 'GET, POST, OPTIONS' always;
            add_header 'Access-Control-Allow-Headers' 'Content-Type, Authorization' always;
            add_header 'Access-Control-Max-Age' 1728000;
            add_header 'Content-Type' 'text/plain charset=UTF-8';
            add_header 'Content-Length' 0;
            return 204;
        }
        
        proxy_pass http://localhost:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 90;
    }
}
EOF

ln -s /etc/nginx/sites-available/optimize-molecule /etc/nginx/sites-enabled/
nginx -t
systemctl restart nginx
```

## 4. DNS and Domain Configuration

After acquiring optimizemolecule.com domain through Cloudflare:

1. Configured A record in Cloudflare:
   - Type: A
   - Name: @
   - Content: 134.209.226.93
   - Proxy status: Proxied

2. Configured CNAME for www subdomain:
   - Type: CNAME
   - Name: www
   - Content: optimizemolecule.com
   - Proxy status: Proxied

3. Configured SSL/TLS in Cloudflare:
   - Set Encryption mode to "Full"
   - Enabled Always Use HTTPS
   - Created Origin Certificate for the server

## 5. Docker Compose Implementation

Created Docker Compose structure:

```bash
mkdir -p /opt/optimize-molecule
cd /opt/optimize-molecule

# Create docker-compose.yml
cat > docker-compose.yml << 'EOF'
version: '3'

services:
  backend:
    image: robertkottelin/optimize-molecule:backend-latest
    container_name: optimize-molecule
    restart: unless-stopped
    ports:
      - "5000:5000"
    volumes:
      - ./data:/app/data
      - ./config:/app/config

  nginx:
    image: nginx:1.25
    container_name: nginx-proxy
    restart: unless-stopped
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - ./nginx/conf:/etc/nginx/conf.d
      - ./nginx/ssl:/etc/nginx/ssl
      - ./nginx/www:/var/www/html
    depends_on:
      - backend
EOF

# Create required directories
mkdir -p data config nginx/conf nginx/ssl nginx/www
```

## 6. Nginx Configuration for Production Domain

Created Nginx configuration for Docker-based implementation:

```bash
cat > /opt/optimize-molecule/nginx/conf/default.conf << 'EOF'
server {
    listen 80;
    server_name optimizemolecule.com www.optimizemolecule.com;
    
    location / {
        return 301 https://$host$request_uri;
    }
}

server {
    listen 443 ssl;
    server_name optimizemolecule.com www.optimizemolecule.com;
    
    ssl_certificate /etc/nginx/ssl/origin.pem;
    ssl_certificate_key /etc/nginx/ssl/origin.key;
    
    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_prefer_server_ciphers on;
    ssl_ciphers ECDHE-ECDSA-AES128-GCM-SHA256:ECDHE-RSA-AES128-GCM-SHA256;
    
    real_ip_header CF-Connecting-IP;
    
    location / {
        proxy_pass http://backend:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 90;
    }
}
EOF
```

## 7. SSL Certificate Implementation

Added Cloudflare Origin Certificate to Docker container:

```bash
# Transfer certificate files to server
cd /opt/optimize-molecule/nginx/ssl

# Created certificate files with appropriate content (user adds content manually)
touch origin.pem
touch origin.key
chmod 644 origin.pem
chmod 600 origin.key
```

## 8. Port Conflict Resolution

Stopped and disabled system Nginx to free ports for Docker services:

```bash
systemctl stop nginx
systemctl disable nginx
```

Verified no port conflicts:

```bash
netstat -tlpn | grep -E ':80|:443'
```

## 9. Container Launch and Validation

Started the container stack:

```bash
cd /opt/optimize-molecule
docker-compose up -d
```

Verified container status:

```bash
docker-compose ps
```

Expected output:
```
      Name                     Command               State                             Ports
----------------------------------------------------------------------------------------------------------------------
nginx-proxy         /docker-entrypoint.sh ngin ...   Up      0.0.0.0:443->443/tcp,:::443->443/tcp,
                                                             0.0.0.0:80->80/tcp,:::80->80/tcp
optimize-molecule   gunicorn --bind 0.0.0.0:50 ...   Up      0.0.0.0:5000->5000/tcp,:::5000->5000/tcp
```

Tested the configuration:

```bash
# Test HTTP redirect
curl -I http://localhost:80

# Test HTTPS functionality
curl -I https://localhost:443 -k
```

## 10. Backend CORS Configuration

Updated CORS configuration in the API application:

```bash
docker exec -it optimize-molecule bash -c "cat > /app/config.json << 'EOF'
{
    \"SQLALCHEMY_DATABASE_URI\": \"sqlite:///users.db\",
    \"SQLALCHEMY_TRACK_MODIFICATIONS\": false,
    \"CORS\": {
        \"origins\": [\"https://robertkottelin.github.io\", \"https://optimizemolecule.com\"],
        \"methods\": [\"GET\", \"POST\", \"OPTIONS\"],
        \"allow_headers\": [\"Content-Type\", \"Authorization\"]
    }
}
EOF"
```

## 11. Frontend Configuration

Updated the frontend to use the domain-based API endpoint:

```javascript
// In App.js
const apiBaseUrl = "https://optimizemolecule.com";
```

## 12. Validation & Troubleshooting

Testing API endpoints:

```bash
# Health check
curl -I https://optimizemolecule.com/health

# CORS pre-flight verification
curl -H "Origin: https://robertkottelin.github.io" \
     -H "Access-Control-Request-Method: POST" \
     -H "Access-Control-Request-Headers: Content-Type" \
     -X OPTIONS \
     -I https://optimizemolecule.com/health
```

## Security Considerations

1. **Firewall Configuration**: Implemented minimal port exposure (80, 443, 22)
2. **SSL/TLS**: Enforced TLS 1.2+ and strong cipher suites
3. **Cloudflare Protection**: Added edge security via Cloudflare proxying
4. **Docker Isolation**: Used container segmentation for application components

## Maintenance Procedures

1. **Certificate Renewal**: Cloudflare Origin Certificates have 15-year validity
2. **Container Updates**:
   ```bash
   cd /opt/optimize-molecule
   docker-compose pull
   docker-compose down
   docker-compose up -d
   ```
3. **Configuration Updates**:
   ```bash
   docker exec -it optimize-molecule bash
   # Edit configuration files as needed
   exit
   docker restart optimize-molecule
   ```