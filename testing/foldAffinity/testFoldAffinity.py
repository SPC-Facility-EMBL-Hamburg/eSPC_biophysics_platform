import unittest
import sys
import inspect
"""

To run the tests $ python testFoldAffinity.py

We are testing that the FoldAffinity class can load all the example files
and that the required R and python packages are going to be installed in the docker image

"""

# Change variables if needed.

foldAffinity_relative_path = "../../foldA_moltenP_apps/foldAffinity"
docker_file_relative_path  = "../../dockerFiles/Dockerfile_foldAffinity_moltenProt_multiStageBuild"
r_packages_relative_path   = "../../foldA_moltenP_apps/install_r_packages.R"

# End of variables to be changed

sys.path.append(foldAffinity_relative_path)

from testHelpers  import *
from foldAffinity import DSF_binding

class testLoadFiles(unittest.TestCase):

    def test_load_nanoDSF_nanothemper_processed(self):

        dsf = DSF_binding()
        dsf.load_nanoDSF_xlsx("../data/nanoThemperProcessed.xlsx")
        dsf.set_signal("330nm")
        self.testLoad(dsf)

        return None

    def test_load_panta_raw(self):

        dsf = DSF_binding()
        dsf.load_panta_xlsx("../data/panta.xlsx")
        dsf.set_signal("330nm")
        self.testLoad(dsf)

        return None

    def test_load_QuantStudio_txt(self):

        dsf = DSF_binding()
        dsf.load_QuantStudio_txt("../data/quantStudio.txt")
        dsf.set_signal("Fluorescence")
        self.testLoad(dsf)

        return None

    def test_load_Thermofluor_xlsx(self):

        dsf = DSF_binding()
        dsf.load_Thermofluor_xlsx("../data/qPCRdemoFile.xls")
        dsf.set_signal("DSF_RFU")
        self.testLoad(dsf)

        return None

    def test_load_tycho_xlsx(self):

        dsf = DSF_binding()
        dsf.load_tycho_xlsx("../data/nanoThemperTycho.xlsx")
        dsf.set_signal("330nm")
        self.testLoad(dsf)

        return None

    def test_load_Agilents_MX3005P_qPCR_txt(self):

        dsf = DSF_binding()
        dsf.load_Agilents_MX3005P_qPCR_txt("../data/MX3005P.txt")
        dsf.set_signal("Fluorescence")
        self.testLoad(dsf)

        return None

    def testLoad(self,dsf=None):

        #Don't test anything
        if dsf is None:
            return None

        self.assertTrue(max(dsf.temps)<200, 'check that the data is in centigrades' )    
        self.assertEqual(len(dsf.fluo.shape), 2,'the fluorescence data is not a 2D matrix')

        return None

class testRequiredPackages(unittest.TestCase):

    def test_R_packages(self):

        """
        Test that all the R required packaged are going to be installed

        """

        _ , files = run_fast_scandir(foldAffinity_relative_path,".r")
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

        _ , files = run_fast_scandir(foldAffinity_relative_path,".py")
        pkgs =  []

        for f in files:
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

        pkgs_to_be_installed = list(set(pkgs_to_be_installed))

        # Helpers, fitting_helpers_unfolded_fraction, and fitting_helpers_thermal_unfolding are python scripts 
        # copy and codecs are installed automatically

        pkgs.remove("fitting_helpers_unfolded_fraction")
        pkgs.remove("fitting_helpers_thermal_unfolding")
        pkgs.remove("helpers")
        pkgs.remove("codecs")
        pkgs.remove("copy")
        pkgs.remove("json")

        self.assertTrue(all(elem in pkgs_to_be_installed  for elem in pkgs))
        
        return None

if __name__ == '__main__':
    unittest.main()