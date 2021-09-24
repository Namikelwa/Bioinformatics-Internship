How to Install R and RStudio
16 comments
How to install R on different computers

Instructions for Mac Users:
Open your internet browser and go to https://www.r-project.org/
In the Getting Started section, follow the download R link available
You will be redirected to a webpage asking you to choose your preferred CRAN (Comprehensive R Archive Network) mirror. You will have to choose the closest mirror to your current location
Follow the Download R for mac OS link available
Select the link containing the R version that is most suited to your Mac OS Version in order to download it (.pkg file)
Open the file and follow the recommended instructions for installation
To check whether the installation was successful: Open your Terminal, simply write R, and the available R version will be displayed to you automatically
Screen you get by invoking R in your terminal window

8. You might also need to install XQuartz if it is not part of your Mac OS Version. Please read carefully the available instructions (in this same webpage or consult the https://www.xquartz.org/ website) on how to install it. Select the link containing the XQuartz version that is most suited to your Mac OS Version in order to download it (.dmg file), then open the file and follow the recommended instructions for installation.

Instructions for Linux Users
Open your internet browser to go to https://www.r-project.org/
In the Getting Started section, follow the download R link
You will be redirected to a webpage asking you to choose your preferred CRAN(Comprehensive R Archive Network) mirror. You will have to select the closest mirror to your current location
Click on the Download R for Linux link
Select your Linux OS system
Comprehensive R Archive Network mirror site screen

6. Follow carefully the installation instructions

7. Example for Ubuntu

To install the complete R system, use:

sudo apt-get update
sudo apt-get install r-base 
Users who need to compile R packages from source [e.g. package maintainers, or anyone installing packages with install.packages()] should also install the r-base-dev package:

sudo apt-get install r-base-dev
8. To check whether the installation was successful:

Open your Terminal
write R
the available R version will be displayed to you automatically
Instructions for Windows Users
Open your internet browser and go to https://www.r-project.org/
In the Getting Started section, follow the download R link
You will be redirected to a webpage asking you to choose your preferred CRAN (Comprehensive R Archive Network) mirror. You will have to select the closest mirror to your current location
Click on the Download R for Windows link
Click on the install R for the first time link
This will open the following webpage
Comprehensive R Archive Network mirror site screen 7. You can directly click on Download R 4.1.1 for Windows

8. If you need more information on how to open and use R under windows, please click on the Installation and other instructions link (immediately under the first link) and read carefully the recommendations given to you.

How to install RStudio - Instructions for Mac, Linux and Windows Users
Open your internet browser and go to https://www.rstudio.com
Select the Products tab and go to RStudio R studio main screen
Alternatively, you can go directly to https://www.rstudio.com/products/rstudio/
Click on the RStudio Desktop link which will redirect you down the same page to the Download RStudio Desktop blue button which is part of the Open Source Edition 5. Download the most suited RStudio Desktop to your OS (Mac, Linux or Windows) installer RStudio screen 6. Open the file and follow the recommended instructions for installation
Discussion
Please use the comments section below to discuss any problems you encounter – facilitators on this course will be checking your comments and some of your fellow learners might be able to help as well. If you still need help, there is also Week 3 Help area where you can also post your question.


Let’s Begin with R!
29 comments
What is R?

R is a “language and environment for statistical computing and graphics”. R is an integrated suite of facilities for data manipulation, calculation and graphical display. Among other things it has:

an effective data handling and storage facility
a suite of operators for calculations on arrays
a large, coherent, integrated collection of intermediate tools for data analysis
graphical facilities for data analysis and display either directly at the computer or on hard-copy
R is very much a vehicle for newly developing methods of interactive data analysis. It has developed rapidly, and has been extended by a large collection of packages. (Please note that these are official definitions taken from https://www.r-project.org/ and http://mercury.webster.edu/aleshunas/R_learning_infrastructure/Introduction_to_R_and_RStudio.html).

How to work under R
Step 1. In your bash window, create a new working subdirectory (we recommend you to use separate working directories for analyses conducted with R, here under Desktop), and move to it

$ cd Desktop
$ mkdir exerciseR 
$ cd exerciseR
Step 2. Start the R program by simply typing R

$ R
This command will open the R environment for you, and a new prompt will appear: “>”

Step 3. To work under R, write any command following the “>” prompt as in Unix

>  
Example of basic operations in R:

> x = 3
> y = 2
> x + y
[1] 5
Step 4. R being created for statistical purposes, it has very helpful built-in functions

> sqrt(3 * 4 + 2 * 5 + 3)
[1] 5
> log(5) + log(10)
[1] 3.912023
> log(50) 
[1] 3.912023
> sum(3 * 4 + 2 * 5 + 3)
[1] 25
Where sqrt=square root, …

Step 5. When working under your current R session, the entities (variables, functions, etc…) that R (you) creates and uses are called “objects”. To display the names of the objects, there are 2 options:

> objects() 
> ls()
Step 6. As in Unix, to remove the 2 “objects” we previously created (Step 3) named x and y, you can use rm

> rm(x, y) 
Step 7. To quit R

> q() 
You will be asked whether you want to save your workspace image

> q()
Save workspace image? [y/n/c]:
Where y=yes (save data and quit), n=no (quit without saving data) and c=cancel (abort the request and continue working under the current R session). If you choose “y”, your objects will be saved in a “.RData” file, and the command lines in a “.Rhistory” file

Important things to know when you work under R
Note 1. As in Unix (man), obtaining help for functions is possible in R (help or ?). Example to obtain information on a command called sum

>  help(sum)
> ?sum
To close this subsection type “q”

Note 2. A “+” symbol might appear when you try to execute a command

+  
This means R is expecting you to complete your command

Note 3. As Unix, R is case sensitive

> help(sum)
> HELP(sum) 
The first command will work, the second will output Error in HELP(sum): impossible to find the function “HELP”

Note 4. As for Unix, you can use the upwards arrow on your keyboard (↑) to go back to the previous command you used.

Note 5. A very detailed introduction on how to use R can be found in https://cran.r-project.org/doc/manuals/r-release/R-intro.html

Discussion
Now try it yourself and discuss in the comment area below:

Question 1. Did you manage to start R?

Question 2. Did you find the R documentation helpful in this document?

Question 3. Did you find the R documentation helpful in the links?

Question 4. Did you manage to get out of the help command as indicated?

