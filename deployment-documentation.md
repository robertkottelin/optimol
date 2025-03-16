# Digital Ocean Droplet Setup for Docker Backend Deployment

## Initial Server Configuration

Log into your Digital Ocean droplet and perform initial setup:

```bash
# Update package lists and install essential packages
apt update && apt upgrade -y
apt install -y ca-certificates curl gnupg lsb-release vim

# Configure proper hostname
hostnamectl set-hostname optimizemolecule
echo "64.227.122.193 optimizemolecule.com" >> /etc/hosts
```

## Docker Installation

Install Docker Engine:

```bash
# Add Docker's official GPG key
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg

# Add Docker repository
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null

# Install Docker Engine
apt update
apt install -y docker-ce docker-ce-cli containerd.io

# Verify Docker installation
docker --version
```

## Nginx Installation and Configuration

Install and configure Nginx as a reverse proxy:

```bash
# Install Nginx
apt install -y nginx

# Create Nginx configuration for the domain
cat > /etc/nginx/sites-available/optimizemolecule.com << 'EOF'
server {
    listen 80;
    server_name optimizemolecule.com www.optimizemolecule.com;

    location / {
        proxy_pass http://localhost:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        
        # CORS headers
        add_header 'Access-Control-Allow-Origin' 'https://robertkottelin.github.io' always;
        add_header 'Access-Control-Allow-Methods' 'GET, POST, OPTIONS' always;
        add_header 'Access-Control-Allow-Headers' 'Origin, X-Requested-With, Content-Type, Accept, Authorization, Cookie, Set-Cookie, x-requested-with' always;
        add_header 'Access-Control-Allow-Credentials' 'true' always;
        add_header 'Access-Control-Expose-Headers' 'Set-Cookie' always;

        # Handle preflight requests
        if ($request_method = 'OPTIONS') {
            add_header 'Access-Control-Allow-Origin' 'https://robertkottelin.github.io' always;
            add_header 'Access-Control-Allow-Methods' 'GET, POST, OPTIONS' always;
            add_header 'Access-Control-Allow-Headers' 'Origin, X-Requested-With, Content-Type, Accept, Authorization, Cookie, Set-Cookie, x-requested-with' always;
            add_header 'Access-Control-Allow-Credentials' 'true' always;
            add_header 'Access-Control-Expose-Headers' 'Set-Cookie' always;
            add_header 'Access-Control-Max-Age' '600' always;
            add_header 'Content-Type' 'text/plain charset=UTF-8';
            add_header 'Content-Length' '0';
            return 204;
        }
    }
}
EOF

# Enable the site configuration
ln -s /etc/nginx/sites-available/optimizemolecule.com /etc/nginx/sites-enabled/
rm -f /etc/nginx/sites-enabled/default

# Test Nginx configuration
nginx -t

# Restart Nginx
systemctl restart nginx
```

## SSL Certificate Configuration with Certbot

Install Certbot and obtain SSL certificates:

```bash
# Install Certbot
apt install -y certbot python3-certbot-nginx

# Obtain SSL certificates
certbot --nginx -d optimizemolecule.com -d www.optimizemolecule.com --non-interactive --agree-tos --email your-email@example.com

# Verify auto-renewal
certbot renew --dry-run
```

## Environment Variables Setup

Create a `.env` file for the application:

```bash
mkdir -p /opt/optimol
cat > /opt/optimol/.env << 'EOF'
STRIPE_SECRET=your_stripe_secret_key
JWT_SECRET_KEY=a_random_secure_key_for_jwt_tokens
FLASK_ENV=production
EOF

# Set proper permissions
chmod 600 /opt/optimol/.env
```

## Docker Image Deployment

Pull and run the Docker image:

```bash
# Pull the latest image
docker pull robertkottelin/optimize-molecule:backend-latest

# Run the container with proper configuration
docker run -d \
  --name optimol-backend \
  --restart unless-stopped \
  -p 5000:5000 \
  -v /opt/optimol/.env:/app/.env \
  -v optimol_data:/app/instance \
  robertkottelin/optimize-molecule:backend-latest
```

## Firewall Configuration

Configure UFW (Uncomplicated Firewall):

```bash
# Install UFW if not already installed
apt install -y ufw

# Set default policies
ufw default deny incoming
ufw default allow outgoing

# Allow SSH, HTTP, and HTTPS
ufw allow ssh
ufw allow http
ufw allow https

# Enable the firewall
ufw --force enable

# Verify the firewall status
ufw status
```

## Verify the Setup

Test that your API is accessible:

```bash
# Check container logs
docker logs optimol-backend

# Test the API endpoint
curl -k https://optimizemolecule.com/health
```

## Configure Cloudflare

Ensure Cloudflare is correctly set up:

1. Confirm DNS records are correct:
   - A record for optimizemolecule.com pointing to 64.227.122.193
   - CNAME record for www pointing to optimizemolecule.com

2. Configure SSL/TLS settings:
   - Set SSL/TLS mode to "Full (strict)" in Cloudflare dashboard
   - Enable "Always Use HTTPS" option

3. Verify CORS settings:
   - Confirm Page Rules or Transform Rules allowing cross-origin requests from GitHub Pages

## Deploying Updates

Create a deployment script for future updates:

```bash
cat > /opt/deploy.sh << 'EOF'
#!/bin/bash
set -e

# Pull the latest image
docker pull robertkottelin/optimize-molecule:backend-latest

# Stop and remove the existing container
docker stop optimol-backend || true
docker rm optimol-backend || true

# Run the new container
docker run -d \
  --name optimol-backend \
  --restart unless-stopped \
  -p 5000:5000 \
  -v /opt/optimol/.env:/app/.env \
  -v optimol_data:/app/instance \
  robertkottelin/optimize-molecule:backend-latest

# Clean up unused images
docker image prune -af
EOF

chmod +x /opt/deploy.sh
```

## Monitoring and Logging

Set up basic monitoring:

```bash
# Install monitoring tools
apt install -y htop glances

# Set up log rotation for Docker
cat > /etc/logrotate.d/docker << 'EOF'
/var/lib/docker/containers/*/*.log {
    rotate 7
    daily
    compress
    missingok
    delaycompress
    copytruncate
}
EOF
```

The setup is now complete. Your Digital Ocean droplet is configured to run your containerized backend application with proper SSL certificates, reverse proxy, and CORS settings for your GitHub Pages frontend.