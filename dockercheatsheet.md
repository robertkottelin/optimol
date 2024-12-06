BACKEND:
build:
sudo docker build -t optimize-molecule .

run:
sudo docker run -p 5000:5000 optimize-molecule

tag:
sudo docker tag optimize-molecule robertkottelin/optimize-molecule:backend-latest

login:
sudo docker login

push:
sudo docker push robertkottelin/optimize-molecule:backend-latest


FRONTEND:
build:
sudo docker build -t optimize-molecule .

run:
sudo docker run -p 8080:8080 optimize-molecule

tag:
sudo docker tag optimize-molecule robertkottelin/optimize-molecule:frontend-latest

login:
sudo docker login

push:
sudo docker push robertkottelin/optimize-molecule:frontend-latest


Stop all containers:
sudo docker stop $(sudo docker ps -q)
sudo docker rm $(sudo docker ps -aq)
