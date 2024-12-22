To host your backend and frontend applications on your Ubuntu machine using Docker, follow these steps:

1. Install Docker and Docker Compose
Make sure Docker and Docker Compose are installed on your Ubuntu machine.

bash
Copy code
# Update package index and install prerequisites
sudo apt update
sudo apt install -y apt-transport-https ca-certificates curl software-properties-common

# Add Docker's official GPG key
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg

# Add Docker's stable repository
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

# Install Docker
sudo apt update
sudo apt install -y docker-ce docker-ce-cli containerd.io

# Verify Docker installation
docker --version

# Install Docker Compose
sudo curl -L "https://github.com/docker/compose/releases/latest/download/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose

# Verify Docker Compose installation
docker-compose --version
2. Pull Backend and Frontend Docker Images
Use the images you pushed to Docker Hub.

bash
Copy code
# Pull backend Docker image
sudo docker pull robertkottelin/optimize-molecule:backend-latest

# Pull frontend Docker image
sudo docker pull robertkottelin/optimize-molecule:frontend-latest
3. Create a docker-compose.yml File
Create a docker-compose.yml file to define how to run both containers.

bash
Copy code
# Create and edit docker-compose.yml
nano docker-compose.yml
Add the following configuration:

yaml
Copy code
version: '3.8'
services:
  backend:
    image: robertkottelin/optimize-molecule:backend-latest
    container_name: optimize-backend
    ports:
      - "5000:5000" # Maps container port 5000 to localhost port 5000
    environment:
      - STRIPE_SECRET=your_stripe_secret_key_here
    restart: always

  frontend:
    image: robertkottelin/optimize-molecule:frontend-latest
    container_name: optimize-frontend
    ports:
      - "3000:3000" # Maps container port 3000 to localhost port 3000
    restart: always
4. Start the Services
Run both containers using Docker Compose.

bash
Copy code
# Start the services
sudo docker-compose up -d

# Verify containers are running
sudo docker ps
You should see two running containers: optimize-backend and optimize-frontend.

5. Access the Applications
Backend: Open http://<your-server-ip>:5000 to verify the backend API is accessible.
Frontend: Open http://<your-server-ip>:3000 to access the frontend application.
If you're running this on your local machine, use http://localhost:5000 and http://localhost:3000.

6. Ensure Containers Restart on Reboot
To make sure the containers start automatically after a reboot, you can set the restart policy in docker-compose.yml (already added as restart: always) or enable the Docker service:

bash
Copy code
sudo systemctl enable docker
7. Monitor Logs
To view logs for debugging or monitoring:

bash
Copy code
# Backend logs
sudo docker logs -f optimize-backend

# Frontend logs
sudo docker logs -f optimize-frontend
8. Optional: Accessing the Backend from Frontend
Make sure the frontend is configured to communicate with the backend. Since both are hosted on the same server, update the frontend API base URL to http://<your-server-ip>:5000.

If this isn't already set, you may need to modify the frontend code or set the apiBaseUrl dynamically through environment variables or a .env file.

To set up domains so users can access your frontend (e.g., optimizemolecule.com) and the frontend can call the backend, you can follow these steps:

1. Obtain a Domain Name
Purchase a domain from a domain registrar like Namecheap, Google Domains, or GoDaddy.
For this guide, we’ll assume your domain is optimizemolecule.com.
2. Point Your Domain to Your Server
You need to associate your domain name with your Ubuntu server's public IP address.

Find Your Server's Public IP Address:

bash
Copy code
curl ifconfig.me
This will return your public IP address.

Set Up an A Record in Your Domain Registrar:

Log in to your domain registrar.
Go to the DNS settings for your domain.
Add an A Record pointing @ (or your root domain) to your server's public IP address.
Example:

makefile
Copy code
Type: A
Host: @
Value: <YOUR_SERVER_IP>
TTL: 1 Hour
To make the backend accessible via api.optimizemolecule.com, add another A Record:

makefile
Copy code
Type: A
Host: api
Value: <YOUR_SERVER_IP>
TTL: 1 Hour
3. Install Nginx as a Reverse Proxy
Nginx will act as a reverse proxy to route requests to your frontend and backend containers.

Install Nginx:

bash
Copy code
sudo apt update
sudo apt install -y nginx
Configure Nginx: Create a configuration file for your domain.

bash
Copy code
sudo nano /etc/nginx/sites-available/optimizemolecule
Add the following content to route traffic to the frontend and backend:

nginx
Copy code
server {
    listen 80;
    server_name optimizemolecule.com www.optimizemolecule.com;

    location / {
        proxy_pass http://localhost:3000; # Frontend container
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_cache_bypass $http_upgrade;
    }

    location /api/ {
        proxy_pass http://localhost:5000/; # Backend container
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_cache_bypass $http_upgrade;
    }
}
Enable the Configuration: Link the configuration file to the sites-enabled directory.

bash
Copy code
sudo ln -s /etc/nginx/sites-available/optimizemolecule /etc/nginx/sites-enabled/
Restart Nginx:

bash
Copy code
sudo nginx -t # Test the configuration
sudo systemctl restart nginx
4. Set Up SSL with Let's Encrypt
Secure your domain with HTTPS using Let's Encrypt.

Install Certbot:

bash
Copy code
sudo apt install -y certbot python3-certbot-nginx
Obtain and Install SSL Certificates:

bash
Copy code
sudo certbot --nginx -d optimizemolecule.com -d www.optimizemolecule.com
For the backend, if you're using a subdomain (e.g., api.optimizemolecule.com), add it as well:

bash
Copy code
sudo certbot --nginx -d api.optimizemolecule.com
Automatically Renew Certificates: Certbot automatically adds a cron job for renewal. You can test it with:

bash
Copy code
sudo certbot renew --dry-run
5. Update Frontend API Base URL
Make sure your frontend calls the backend using the correct domain. Update the apiBaseUrl in your frontend code to:

javascript
Copy code
const apiBaseUrl = "https://api.optimizemolecule.com";
Rebuild and redeploy the frontend container if needed:

bash
Copy code
sudo docker-compose down
sudo docker-compose up -d
6. Verify the Setup
Frontend: Visit https://optimizemolecule.com to confirm the frontend is accessible.
Backend: Test the backend by calling an endpoint (e.g., https://api.optimizemolecule.com/optimize) from the browser or a tool like Postman.
7. Optional: Use Cloudflare for Additional Features
If desired, use Cloudflare for:

Faster DNS resolution.
Free SSL certificates (if you don't want to use Let's Encrypt).
Protection against DDoS attacks.
Now, your application should be accessible at https://optimizemolecule.com for the frontend and https://api.optimizemolecule.com for the backend! 🎉
