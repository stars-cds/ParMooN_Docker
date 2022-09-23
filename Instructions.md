# ParMooN Execution Instructions
---

## 1. Prerequisites
---

Check for all the installed libraries 

- `gcc --version`
- `git --version`
- `cmake --version`

Check for installation of docker and paraview

- `whereis paraview`
- `whereis docker`

## 2. Pull the code from the github repo
---

### Create a Directory 
Before clonning the code, create a directory called `PARMOON_CODES_NEW` in your home directory 
`/home/user/PARMOON_CODES_NEW`

Either do it from command line, such as 

`mkdir /home/user/PARMOON_CODES_NEW` 

or do it from the terminal

### Navigate to Directory
Then, get into the folder using 

`cd /home/user/PARMOON_CODES_NEW`


### Clone the Git repo 

`git clone https://github.com/thivinanandh/ParMooN_Docker.git`


Once its done, you will see a folder called ParMooN_Docker inside your current terminal 



## 2. Configuring Docker container 
---

### Check for access to docker 

Please try on the following commands to check the usage of docker on your system 
* `docker ps`
* `docker images`


If this process returns an error due to any sudo access issue, then run the following command on your terminal

`su - user`

and try the docker commands again, If you are still not able to run the docker commands , then please restart the system 



### Downloading Docker Images

Check for the list of installed docker images in your system by typing

`docker images `

and check if there is an image called "parmoon/ParMooN_v2.0". If the image does not exists, then pull the image from dockerhub using 

`docker pull parmoon/parmoon_v2.0`


Once its pulled, we check for the image on `docker ps` to validate whether it has been pulled or not. 

## 3. Running the Docker container
---


Now, we need to run the docker container and mount the folder in which our codes are present into our main system 


`docker run -it -v /home/user/PARMOON_CODES_NEW:/home/parmoon parmoon/parmoon_v2.0`

The container should immediately give a terminal with root access



You can check the docker process running on your system with 
`docker ps`


## 4. Connecting VS code to docker instance
---

Open vscode by either typing `code` on terminal or using the GUI

Then, Navigate to the extensions section ( the one with 4 boxes ) and then search for the following extensions and install the extensions


* docker
* Remote- Container

Note: Both the extensions are from microsoft

Now, press `cntrl + shift + p` on your keyboard to open command panel and type "Attach Container"

You will find an option called `remote containers: Attach to running containers`. Click on that and select the running parmoon container which will be visible.


## 5. Navigating to the Folder
---
Once the VS code window is open, press on open folder and select `/home/parmoon/`

Now, you should be able to see the codes being available in the left panel .




## 6. Running the PARMOON CODE 
---

Create a directory inside the `ParMooN_Docker` directory clled BUILD 

```
cd /home/parmoon/ParMooN_Docker/
mkdir BUILD
cd BUILD
```


Once you are inside , then run 

```
cmake ../
```

Which will generate the make file


then run 

```
make -j
```


This will have generated an output file in the given 







