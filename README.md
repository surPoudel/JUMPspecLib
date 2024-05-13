# JUMPlib #

 * [Introduction](#introduction)
 * [Basic Installation](#basic-installation)
 * [Input Data](#input-data)
 * [Sample Data](#sample-data)
 * [Basic Installation](#basic-installation)
 * [JUMPlib Commands](#jumplib-commands)
 * [Test Data Exercise](#test-data-exercise)
 * [Input and Output Data Organization ](#input-and-output-data-organization)

---

## Introduction ##

JUMPlib is a specialized tool designed for searching TMT-based proteomics data. The JUMPlib program enables TMT library generation, database search, identification filtering, and protein quantification. To evaluate the performance of the JUMPlib program, we conducted an assessment using a large-scale TMT data set. In addition, the JUMPlib program can readily be adapted for label-free library generation and database search. Moreover, we curated comprehensive 11-plex and 18-plex TMT libraries from human brain samples, providing valuable resources to the research community.


[Top of page](#JUMPlib)

----
## JUMPlib Publication:
  * The manuscript is submitted and this part will be updated later.
  * If you use JUMPlib as part of a publication, please include this reference.

[Top of page](#JUMPlib)

---

## Basic Installation ##
The installation is tested in the linux system and HPC servers but this should work properly in windows and mac too. We highly recommend installing JUMPspecLib in a virtual environment, for example using the anaconda or miniconda package manager
1. Create a virtual enviroment and install required packages.

Here are some commands that would create an `anaconda` python environment for
running JUMPspecLib:

```
conda create -n jumplib python=3.8
conda activate jumplib
conda install numpy pandas=1.5.3 matplotlib scipy seaborn statsmodels pyteomics rpy2
# we recommend pandas=1.5.3 particular version of pandas because the libraries were created using this version but if you are going to create your own library using JUMPspecLib you can install latest version

```

2. Place the JUMPlib distribution source in the desired location (call
this `<path to JUMPlib>`)


----

  * Obtaining JUMPspecLib source 
You can obtain the latest version of JUMPspecLib from git; simple clone the
git repository:

```
    git clone https://github.com/surPoudel/JUMPspecLib.git
```

in the directory _where you would like JUMPlib to be installed_ (call this directory `<path to JUMPlib>`).  Note
that JUMPlib does not support out-of-place installs; the JUMPlib git
repository _is_ the entire installation.  

[Top of page](#JUMPlib)

----

## JUMPlib Commands ##

Once the conda environment (JUMPlib) is activated

1. make a working directory
2. keep all the mzXML files in the same directory
3. copy the parameter file from [parameterFiles](./parameterFiles) to the same directory
4. make necessary changes for the parameters 
5. Run the command below

a. Library generation
```
Preprocessing 
jump_lib -pp jump_lib_preprocess.params *.mzXML 

Library generation
jump_lib -d jump_lib_gen.params

Library merging
jump_lib -d_merge jump_lib_specLibMerge.params
```
Note: We also provide the comprehensive TMT libraries so you may skip Library generation becasue it takes time.

b. Library searching (with presearch and without presearch)
```
jump_lib -pp jump_preprocess.params *.mzXML
jump_lib -s jumplib_search.params
jump_lib -pp_s jumplib_search.params
```

c. Filter the search results
```
jump_lib -f jumplib_filter.params
```

d. Quantification of filtered dataset
```
jump_lib -q jump_lib_q.params
```

[Top of page](#JUMPlib)

----

## Test Data Exercise ##

1. Make a directory called **Test**
2. Copy mzXML files from [sample_data](./sample_data) to **Test**
3. Copy --------- from [parameterFiles](./parameterFiles) to **Test**
4. Open the parameter file copied from step # 3 and change the following information
    * write the changes here
    #### Download [database](https://drive.google.com/drive/folders/1RMV537_0cBQER5BoP3GWRMR9HWC9MwDl?usp=sharing)
    * pitfile = replace the path with the folder that has the sample pitfile name for eg. /home/spoudel1/JUMP-ptm/resources/database/humanComprehensive_v1_ft_mc2_c57_TMT_K229.pit 
    #### Download [pitfile](https://drive.google.com/drive/folders/1RMV537_0cBQER5BoP3GWRMR9HWC9MwDl?usp=sharing)

    #### Download test database and pitfile here [database_pitfile](https://drive.google.com/drive/folders/1RMV537_0cBQER5BoP3GWRMR9HWC9MwDl?usp=sharing)


[Top of page](#JUMPlib)

----

Maintainers
----

* To submit bug reports and feature suggestions, please contact

  **Suresh Poudel (suresh.poudel@stjude.org)**

[Top of page](#JUMPlib)

----

