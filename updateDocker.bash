#!/bin/bash

if [[ $EUID -ne 0 ]]; then
   echo ''	
   echo "This script must be run with sudo." 
   echo ''
   exit 1
fi

echo "Which Docker image do you want to build?"
echo "1. FoldAffinity and MoltenProt"
echo "2. ThermoAffinity"
echo "3. PhotoMol"
echo "4. Raynals"
echo "5. ChiraKit"

read -p "Enter your choice [1-5]: " choice

askIfTestsWereDone () {

    echo ''
    read -p "Did you run the test $1 from the testing folder (y/n)?" CONT

    if [ "$CONT" != "y" ]; then echo 'update aborted' && exit ; fi

}

askIfInfoWasUpdated () {

    read -p "Did you add a new function to load a different kind of input file (y/n)?" CONT2

    if [ "$CONT2" = "y" ]; then

        echo  -e '     \nOkay... Did you check that'
        echo  -e '     1) There is an example file available in the corresponding www folder?'
        echo  -e '     2) The user documentation was updated?'
        echo  -e '     3) The user guide was updated?\n'

        read -p "Would you like to continue (y/n)?" CONT3
    
        if [ "$CONT3" != "y" ]; then exit ; fi

    fi    
    
}

printMessageHowToTestDocker () {

    echo -e "\nTo test $1 in the developer machine run $ sudo docker run -p 3838:3838 $2"
    echo -e  "and then open http://localhost:3838/ in your browser"

}
 
case $choice in

    1)
        
        askIfTestsWereDone 'testFoldAffinity.py and testMoltenProt.py'
        askIfInfoWasUpdated

        echo -e "\nBuilding Image for FoldAffinity and MoltenProt..."
        docker build -t spc_apps_docker_container -f ./dockerFiles/Dockerfile_foldAffinity_moltenProt_multiStageBuild .
       
        printMessageHowToTestDocker 'MoltenProt' 'spc_apps_docker_container'

        ;;
    2)

        askIfTestsWereDone 'testThermoAffinity.py'
        askIfInfoWasUpdated

        echo -e "\nBuilding Image for ThermoAffinity..."
        docker build -t thermo_affinity -f ./dockerFiles/Dockerfile_thermoAffinity_multiStageBuild .
       
        printMessageHowToTestDocker 'ThermoAffinity' 'thermo_affinity'

        ;;
    3)

        askIfTestsWereDone 'testPhotoMol.py'
        askIfInfoWasUpdated

        echo -e "\nBuilding Image for PhotoMol..."
        docker build -t photo_mol -f ./dockerFiles/Dockerfile_refeynApp_multiStageBuild .
       
        printMessageHowToTestDocker 'PhotoMol' 'photo_mol'

        ;;
    4)

        askIfTestsWereDone 'testRaynals.py'
        askIfInfoWasUpdated

        echo -e "\nBuilding Image for Raynals..."
        docker build -t raynals -f ./dockerFiles/Dockerfile_dlsApp_multiStageBuild .
       
        printMessageHowToTestDocker 'Raynals' 'raynals'

        ;;

    5)

        #askIfTestsWereDone 'testRaynals.py'
        askIfInfoWasUpdated

        echo -e "\nBuilding Image for Circular dichroism..."
        docker build -t chirakit -f ./dockerFiles/Dockerfile_cdApp_multiStageBuild .
       
        printMessageHowToTestDocker 'ChiraKit' 'chirakit'

        ;;
    *)
        echo "Invalid choice. Please select a valid option [1-5]."
        ;;

esac

