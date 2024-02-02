import unittest
import sys
import inspect
"""

To run the tests $ python testRaynals.py

We are testing that the Raynals class can load all the example files
and that the required R and python packages are going to be installed in the docker image

"""
# Change variables if needed.

raynals_relative_path = "../../dynamicLightScatteringApp/appFiles/Raynals"
docker_file_relative_path  = "../../dockerFiles/Dockerfile_dlsApp_multiStageBuild"
r_packages_relative_path   = "../../dynamicLightScatteringApp/appFiles/install_r_packages.R"

# End of variables to be changed

sys.path.append(raynals_relative_path)

from testHelpers  import *
from dlsAnalyzer  import dls_experiment

class testLoadFiles(unittest.TestCase):

    def test_load_raynals(self):

        dls = dls_experiment()
        dls.loadWyatFile('../data/testDLS.csv')

        self.assertTrue(hasattr(dls,'sampleInfo'),        'check  that we loaded the sampleInfo data')
        self.assertTrue(hasattr(dls,'time'),       'check  that we loaded the time data')
        self.assertEqual(len(dls.autocorrelationOriginal.shape), 2,'the signal data is not a 2D matrix')

        return None

class testRequiredPackages(unittest.TestCase):

    def test_R_packages(self):

        """
        Test that all the R required packaged are going to be installed

        """

        _ , files = run_fast_scandir(raynals_relative_path,".r")
        files = [f for f in files if "git" not in f and ".Rhistory" not in f]
        pkgs =  []

        for f in files:
            with open(f,"r") as rf:
                ls = rf.read().splitlines()
                
                for l in ls:

                    if "library" in l and len(l.split(")")) == 2:
                        pkgs.append(l.split(")")[0].split("(")[1])

                    if "packages <-" in l:
                        l2 = l.split("c(")[1]
                        l3 = l2.split("\"")

                        l3 = [x for x in l3 if x not in ["",",",")"]]

                        pkgs += l3
                
                    if "::" in l:
                        l2  = l.strip()
                        
                        # l3 returns the string before ::
                        l3  = l2.split("::")[0]
                        # Replace "(" with spaces
                        l4  = l3.replace("("," ")
                        # Get last element
                        pkg = l4.split(" ")[-1]

        pkgs = list(set(pkgs))

        pkgs_to_be_installed = []
        with open(r_packages_relative_path,"r") as rf:
            ls = rf.read().splitlines()
            for l in ls:
                if "install.packages" in l:
                   l2 = l.split("packages(\"")[1]
                   l3 = l2.split("\"")[0]
                   pkgs_to_be_installed.append(l3)

        pkgs_to_be_installed.append("htmltools") # Installed automatically

        self.assertTrue(all(elem in pkgs_to_be_installed  for elem in pkgs))
      
        return None

    def test_python_packages(self):

        """
        Test that all the python required packages are going to be installed

        """

        _ , files = run_fast_scandir(raynals_relative_path,".py")
        pkgs =  []

        for f in files:

            if 'test_reticulate.R' in f:
                continue

            with open(f,"r") as rf:
                ls = rf.read().splitlines()
                
                for l in ls:

                    if "import " in l:
                        if "from" in l:
                            pkg = l.split(" ")[1]
                        else:
                            pkg = l.split("import ")[1]

                        if " as" in pkg:
                            pkg = pkg.split(" as")[0]
                        if "." in pkg:
                            pkg = pkg.split(".")[0]

                        pkgs.append(pkg.replace(" ",""))
                 
        pkgs = list(set(pkgs))
        
        pkgs_to_be_installed = []
        with open(docker_file_relative_path,"r") as rf:
            ls = rf.read().splitlines()
            for l in ls:

                if "r-reticulate -c" in l:
                    
                    pkg = l.split("r-reticulate -c")[1].split()[1]
                    if "=" in pkg: 
                        pkg = pkg.split("=")[0]
                    else:
                        pkg = pkg.split("&&")[0]

                    pkgs_to_be_installed.append(pkg.replace(" ","")) 

                if 'pip install' in l:
                    pkg = l.split("install")[1].split()[0]
                    pkgs_to_be_installed.append(pkg.replace(" ","")) 

        pkgs_to_be_installed = list(set(pkgs_to_be_installed))

        # python scripts 
        pkgs.remove("simulation_helpers")
        pkgs.remove("helpers")
        pkgs.remove("loadDLSdataHelpers")

        # installed by default
        pkgs.remove("math")

        self.assertTrue(all(elem in pkgs_to_be_installed  for elem in pkgs))
        
        return None

if __name__ == '__main__':
    unittest.main()