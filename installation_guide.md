
Last time updated: August 2024

## General description

The eSPC webserver is currently hosted through ShinyProxy, with each tool contained in a docker image. To modify an existing app, you need to have proficiency in Python and R. 

0) Clone the Repo and the submodules

```
git clone https://github.com/SPC-Facility-EMBL-Hamburg/eSPC_biophysics_platform
cd eSPC_biophysics_platform
git submodule update --init --recursive
```

## If you want to try the tools but don't plan to edit the code frequently

1) Install docker (https://docs.docker.com/engine/install/)

```
sudo bash updateDocker.bash
``` 
 
## To install individual tools (for developers)

Follow the README.md document:

| Tool                  | Path                                          |
| :---                  |:----------------------------------------------|
| MoltenProt            | differentialScanningFluorimetryApps/README.md |
| FoldAffinity          | differentialScanningFluorimetryApps/README.md |
| ThermoAffinity        | microscaleThermophoresisApp/README.md         |
| PhotoMol              | massPhotometryApp/README.md                   | 
| Raynals               | dynamicLightScatteringApp/README.md           | 
| ChiraKit              | circularDichroismApp/README.md                | 

## Installation of the whole eSPC platform 

To install the eSPC platform in the development machine, follow these instructions:

1) Clone or download this repo
2) Install docker       (tested with Docker version 27.1.1, https://docs.docker.com/engine/install/)
3) Install java         (tested with openjdk 22.0.2)
4) Install shinyproxy   (tested with shinyproxy3.1.1, https://www.shinyproxy.io/downloads/)
5) Create the file '/etc/systemd/system/docker.service.d/override.conf' and fill it with the following text:

    [Service]  
    ExecStart=  
    ExecStart=/usr/bin/dockerd -H unix:// -D -H tcp://127.0.0.1:2375  

This file will be used for setting up the Docker systemd service

6) Modify the file paths shown in application.yml according to your machine
7) Start the docker service
```
systemctl daemon-reload
systemctl start docker
```
8) Build the docker images
```
sudo bash updateDocker.bash
```

9) To test the correct installation run shinyproxy
```
java -jar shinyproxy-x.y.z.jar
```    
10) Open localhost:8080 inside any browser

Additionally, you can enter and inspect each of the docker images by running

```
docker run -t -i docker_container_name /bin/bash
```  

## General code structure

The following tree represents how the files are organised and their purpose:

    README.md                                   # Information about the eSPC platform
    fresh_install_guide.md                      # Manual to install from scratch the eSPC platform in a new host
    updateDocker.bash                           # Script to update the eSPC docker images
    application.yml                             # Example configuration file for ShinyProxy
    dockerFiles/                                # Contains the docker files to build the docker images
    dynamicLightScatteringApp/                  # Files for the Raynals app
    circularDichroismApp/                       # Files for the ChiraKit app 
    differentialScanningFluorimetryApps/        # Files for the FoldAffinity and Moltenprot apps (they share the same docker image)
    massPhotometryApp/                          # Files for the PhotoMol app
    microscaleThermophoresisApp/                # Files for the ThermoAffinity app
    templates                                   # HTML and CSS files, user documentation files, and images 
       |-- custom                               
       |   |-- app.html                         # Template file used by ShinyProxy to render the user interface of a Shiny web application
       |   |-- assets                           
       |   |   |-- apps_user_documentation      # User documentation PDFs
       |   |   |-- css                          
       |   |   |   |-- custom.css               # Used to customize the visual appearance of a Shiny web application
       |   |   |-- img                          # Images of the apps (now presented at spc.embl-hamburg.de)
       |   |-- index.html                       # Front page of the eSPC platform
