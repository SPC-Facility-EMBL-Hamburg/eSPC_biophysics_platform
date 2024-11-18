
Last time updated: August 2024

## General description

The eSPC webserver is currently hosted through ShinyProxy, with each tool contained in a docker image. To modify an existing app, you need to have proficiency in Python and R. 

0) Clone the Repo and the submodules

```
git clone https://github.com/SPC-Facility-EMBL-Hamburg/eSPC_biophysics_platform
cd eSPC_biophysics_platform
git submodule update --init --recursive
```

## If you want to try the tools but don't plan to edit the code

1) Install docker (https://docs.docker.com/engine/install/)


```
sudo bash updateDocker.bash
``` 
 
## To install individual tools (for developers)

Follow the README.md document:

| Tool                  | Path                                    |
| :---                  |     :---                                |
| MoltenProt            | foldA_moltenP_apps/README.md            |
| FoldAffinity          | foldA_moltenP_apps/README.md            |
| ThermoAffinity        | thermoA_app/README.md                   |
| PhotoMol              | refeynApp/README.md                     | 
| Raynals               | dynamicLightScatteringApp/README.md     | 
| ChiraKit              | circularDichroismApp/README.md          | 

## Installation of the whole eSPC platform 

To install the eSPC platform in the development machine, follow this instructions:

1) Clone or download this repo
2) Install docker       (tested with Docker version 23.0.1, https://docs.docker.com/engine/install/)
3) Install java         (tested with openjdk 11.0.18)
4) Install shinyproxy   (tested with shinyproxy3.0, https://www.shinyproxy.io/downloads/)
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
    developer_manual.md                         # Manual to understand how to modify or develop new tools
    fresh_install_guide.md                      # Manual to install from scratch the eSPC platform in a new host
    updateDocker.bash                           # Script to update the eSPC server in the development/ production machine
    application.yml                             # Example configuration file for ShinyProxy
    dockerFiles/                                # Contains the docker files to build the docker images
    dynamicLightScatteringApp/                  # Contains the files for the Raynals app
    circularDichroismApp/                       # Contains the files for the circular dichroism app 
    foldA_moltenP_apps/                         # Contains the files for the FoldAffinity and Moltenprot apps (they share the same docker image)
    refeynApp/                                  # Contains the files for the PhotoMol app
    thermoA_app/                                # Contains the files for the ThermoAffinity app
    templates                                   # Contains the html and css files, user documentation files, and images 
       |-- custom                               
       |   |-- app.html                         # Template file used by ShinyProxy to render the user interface of a Shiny web application
       |   |-- assets                           
       |   |   |-- apps_user_documentation      # User documentation PDFs
       |   |   |-- css                          
       |   |   |   |-- custom.css               # Used to customize the visual appearance of a Shiny web application
       |   |   |-- img                          # Images of the apps (now presented at spc.embl-hamburg.de)
       |   |-- index.html                       # HTML code of the webpage the users are going to see when they enter the eSPC platform
    testing                                     # Contains testing scripts: load of input files and installation of neccesary libraries in docker
       |-- data                                 # Input data to be tested
       |-- foldAffinity                         # Tests for FoldAffinity  
       |-- moltenProt                           # Tests for MoltenProt 
       |-- photoMol                             # Tests for PhotoMol 
       |-- raynals                              # Tests for Raynals 
       |-- thermoAffinity                       # Tests for ThermoAffinity 

