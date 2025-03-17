ssh root@64.227.122.193

# Pull frontend image
docker pull robertkottelin/optimize-molecule:frontend-latest

# Stop existing frontend container if present
docker stop optimol-frontend || true
docker rm optimol-frontend || true

# Run frontend container
docker run -d \
  --name optimol-frontend \
  --restart unless-stopped \
  -p 3000:3000 \
  robertkottelin/optimize-molecule:frontend-latest

  # Restart Nginx
systemctl restart nginx

# Check container logs
docker logs optimol-frontend

# Test the API endpoint
curl -k https://optimizemolecule.com/api/health

# List images
docker images
# Prune af if needed 
docker image prune -af

# Nginx:
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

    # Frontend (React application)
    location / {
        proxy_pass http://localhost:3000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    

    # Backend API (Flask application)
    location /api/ {
        # Strip /api prefix before forwarding to backend
        rewrite ^/api(/.*)$ $1 break;
        
        # Proxy to backend server
        proxy_pass http://localhost:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}
EOF