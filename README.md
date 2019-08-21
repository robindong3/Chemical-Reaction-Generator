# An Intellegent Chemical Reaction Generator 

##### Ver 1.0

This is a program to generate reactions and species in low-temperature plasmas. Those reactions and species are sorted with their enthalpies.

### Prerequisites

Please clone this project into a Windows computer with:

```
git clone https://github.com/robindong3/Chemical-Reaction-Generator.git
```

There are two ways of running this project:

#### 1. Run the main.exe file

Theoretically, directly running the main.exe file does not require pre-installed python and related packages. However, some enviroments require an installation of Visual C++, otherwise they will encounter an error with missing api-ms-win-crt-runtime|1-1-0.dll. The Visual C++ can be downloaded from the following link.

```
https://support.microsoft.com/en-gb/help/2977003/the-latest-supported-visual-c-downloads
```
#### 1. Directly run the main.py file

This program is running on a python 3.7 enviroment. If you want to directly run the main.py, please ensure you have the tqdm package and the numpy package in your python 3.7 enviroment. To download those packages you can use pip install as:

```
pip install numpy
pip install tqdm
```

And then, use

```
python .\main.py
```

to run the code.

## Manually Adding Enthalpy of Species

For enthalpies which are not included in the code right now, and you want to insert them manually, please type in the enthalpies into the EnthalpyInsertManually.txt file in the root directory. Please keep your mannully input species and enthalpies in the same style as the examples in the file.

## Tests

Here I will list three tests of this code.

### 1. The N2,H2 Chemistry

This test can give you a general overview of this program.

1. Run the program by the two methods above.
2. Insert N2,H2 when the program says 'Please insert input gases, separated by comma': 
3. Insert Y to get a full reaction list.
4. Insert 2 to set the maximum Nitrogen atoms in each molecular as 2.
5. Directly press Enter without inerting anything for setting the overall scale of molecules as the defult.
6. Wait for the calculation completing.
7. Directly press Enter without inerting anything to exit the program.

In conclusion, this test has the following steps.  
```
N2,H2
Y
2
[Enter]
[Enter]
```

The results are saved the results folder. The two most important result files are:

EquationWithEnth.txt saves the reactions sorted with their enthalpies.
SpiciesList1.txt saves the species sorted with their enthalpies.

Other files in results folder are also useful:

EquationDict.json saves all reactions in python dictionaries, this can be read by other programs in the future research.
EquationStr.txt shows reaction without enthalpies
SpiciesList2.txt saves the unphysical species which are detected by the program, these species are not included in the SpiciesList1.txt.
EquationStr2.txt saves the reactions relates to those unphysical species in SpiciesList2.txt

### 2. The CF4,SiO2 Chemistry
This test can give you the performance of this program on more complex gasses.

1. Run the program by the two methods above.
2. Insert CF4,SiO2 when the program says 'Please insert input gases, separated by comma': 
3. Insert Y to get a full reaction list.
4. Insert 1 to set the maximum Carbon atoms in each molecular as 1.
5. Insert 1 to set the maximum Silicon atoms in each molecular as 1.
6. Directly press Enter without inerting anything for setting the maximum O as the defult.
7. Insert 4 to set the overall scale of molecules as 4.
8. Wait for the calculation completing.
9. Directly press Enter without inerting anything to exit the program.

In conclusion, this test has the following steps.  
```
CF4,SiO2
Y
1
1
[Enter]
4
[Enter]
```

The defult O number are calculated from the sum of other atoms, so it is better to make it as defult while using the program. The overall scale of molecules includes the number of O. Results can be seen in the results folder (Same as the first test)


### 3. The SF6,O2,CF6 chemistry

This test can give you the performance of this program on a larger number of input gasses. This test will take longer time due to the large number of reactions.

1. Run the program by the two methods above.
2. Insert SF6,O2,CF6 when the program says 'Please insert input gases, separated by comma': 
3. Insert Y to get a full reaction list.
6. Directly press Enter without inerting anything for setting the maximum C as the defult.
6. Directly press Enter without inerting anything for setting the maximum S as the defult.
6. Directly press Enter without inerting anything for setting the maximum O as the defult.
5. Directly press Enter without inerting anything for setting the overall scale of molecules as the defult.
8. Wait for the calculation completing.
9. Directly press Enter without inerting anything to exit the program.

In conclusion, this test has the following steps.  
```
SF6,O2,CF6
Y
[Enter]
[Enter]
[Enter]
[Enter]
[Enter]
```

This chemistry takes 78.6 seconds to complete on my computer. Results can be seen in the results folder (Same as the first test)


## Authors

* **Hongyang Dong** - hongyang.dong.18@ucl.ac.uk

