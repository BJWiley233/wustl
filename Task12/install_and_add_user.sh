#! /bin/bash

# https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository
# this already builds ubuntu:16.04 as separate image
sudo apt-get update

sudo apt-get install \
	apt-transport-https \
	ca-certificates \
	gnupg-agent \
	software-properties-common
    
   
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo apt-key fingerprint 0EBFCD88

sudo add-apt-repository \
	"deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
   
sudo apt-get update
apt-cache madison docker-ce

VERSION_STRING=5:20.10.1~3-0~ubuntu-focal
sudo apt-get install docker-ce=$VERSION_STRING docker-ce-cli=$VERSION_STRING containerd.io


# https://askubuntu.com/questions/477551/how-can-i-use-docker-without-sudo
sudo gpasswd -a $USER docker
grep /etc/group -e "docker"
newgrp docker
docker run hello-world
