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



How to Create and Manipulate Variables and Vectors in R
25 comments
Variables and Vectors in R

As we just saw previously in Step 3.4, when working under your current R session, the entities that R creates and uses are called “objects”.

There is a wide range of objects that can be created and that can contain variables, arrays of numbers, character strings, functions….

They can be classified into different types: Vectors, Lists, Data frames, Matrices, Factors, or Functions. In this course we will show you how to create and manipulate Vectors, Lists, and Data frames. For the sake of time and course content we will not cover the rest of the object types, but we encourage you to find more information on these in the following link: https://cran.r-project.org/doc/manuals/r-release/R-intro.html Let’s see together how Variables and Vectors can be created and manipulated.

What is a Variable and what is a Vector
Definition of a variable

Variables are objects in R that you can use to store values. It can consist of a single value, basic or complex arithmetic operations, or even be more complex such as a column in a data matrix or a data frame. We will see these complex forms in the following steps of this course.

Definition of a vector

A vector is substantially a list of variables, and the simplest data structure in R. A vector consists of a collection of numbers, arithmetic expressions, logical values or character strings for example. However, each vector must have all components of the same mode, that are called numeric, logical, character, complex, raw.

How to create and manipulate Variables
Step 1. We recommend you to work in the same working sub-directory that you created previously (in Step 3.4 of the course) for analyses conducted with R. If the sub-directory is not created yet or mistakenly removed, please do it again, and launch R

$ mkdir exerciseR 
$ cd exerciseR
$ R
Step 2. If you forgot before launching R, there is another option you can use to make sure to set the correct working directory using the “setwd()” command, and then check your position using the “getwd()” command (this will also be helpful in RStudio). You can use “getwd()” in R as you used “pwd” in Unix

Launch R

$ R
Go to your working directory

> setwd("/Users/imac/Desktop/exerciseR")
> getwd()
[1] "/Users/imac/Desktop/exerciseR"
Step 3. Let’s create a simple variable called x. We need to assign elements to this variable. The assignment to a variable can be done in 2 different but equivalent ways, using either the “<-“ or “=” operators. You can retrieve the value of x simply by typing x

> x <- 3 * 4 + 2 * 5 + 3
> x = 3 * 4 + 2 * 5 + 3
> x
[1] 25
Step 4. Let’s create another variable called y that can either contain a new value or for example contain a basic or more complex operation on the first variable x

> y <- x^4 - 4*x + 5
> y
[1] 390530
Note 1. Naming a Variable is not trivial and must be done appropriately:

Variable names can contain letters, numbers, underscores and periods
Variable names cannot start with a number or an underscore
Variable names cannot contain spaces at all
> x.length <- 3*2
> x.length
[1] 6
> _x.length <- 3*2
Error : unexpected input in "_"
> 3x.length <- 3*2
Error : unexpected symbol in "3x.length"
Note 2. Long Variable names are allowed but must be formatted using:

Periods to separate words: x.y.z
Underscores to separate words: x_y_z
Camel Case to separate words: XxYyZz
> x.length <- 3*2
> x.length
[1] 6
> x_length <- 3*2
> x_length
[1] 6
> xLength <- 3*2
> xLength
[1] 6
How to create and manipulate Vectors
Step 1. A vector can be created using an in-built function in R called c(). Elements must be comma-separated.

> c(10, 20, 30)
[1] 10 20 30
Step 2. A vector can be of different modes: numeric (and arithmetic), logical, or can consist of characters

> c(1.1, 2.2, 3.5)     			# numeric
[1] 1.1 2.2 3.5
>
> c(FALSE, TRUE, FALSE)		# logical
[1] FALSE TRUE FALSE
>
> c("Darth Vader", "Luke Skywalker", "Han Solo")		# character
[1] "Darth Vader"   "Luke Skywalker"   "Han Solo"
Note. Please note that when the value is a character data type, quotations must be used around each value, such as in “Han Solo”

Step 3. A vector can be assigned to a variable name, using 3 methods: either using the “<-“ or “=” operators or the assign function. You will very rarely see the last method which is to revert the order of assignment

> assign("x", c(10, 20, 30))
> x
[1] 10 20 30
>
> x <- c(10, 20, 30)
> x
[1] 10 20 30
>
> x = c(10, 20, 30)
> x
[1] 10 20 30
>
> c(10, 20, 30)  ->  x
> x
[1] 10 20 30
Step 4. In R, an object must be defined by properties of its fundamental components, such as the mode, that can be retrieved by the function “mode()” and the length by the function “length()”. An empty vector can be created and may still have a mode

> v <- numeric()
> w <- character()
> mode(x)
[1] "numeric"
> mode(v)
[1] "numeric"
> mode(w)
[1] "character"
> 
> length(x)
[1] 3
Step 5. Basic operations with numeric vectors

> x <- c(10, 20, 30)
> x
[1] 10 20 30
> 1/x
[1] 0.10000000 0.05000000 0.03333333
Step 6. A vector can be used in arithmetic expressions and/or as a combination of existing vectors

> x <- c(10, 20, 30)
> y <- x*3+4
> y
[1] 34 64 94
> z <- c(x, 0, 0, 0, x)
> z
[1] 10 20 30  0  0  0 10 20 30
> w <- 2*x + y + z
> w
[1]  64 124 184  54 104 154  64 124 184
Note how the addition is sequential in this last case (value 1 of x is multiplied by 2, then added to value 1 of y then added to value 1 of z, etc…)

Step 7. A vector can use built-in functions in R, such as mean() to calculate the mean of a certain object (here x), var() to calculate its variance, and sort() to sort the content here of object z.

> mean(x)
[1] 20
> var(x)
[1] 100
> sort(z)
[1]  0  0  0 10 10 20 20 30 30
Step 8. R uses built-in functions and operators to generate regular sequences. Here are examples of how to use rep() to repeat items (arguments needed are the value to repeat and the number of repeats) and seq() (arguments needed are the start, the end, and the interval) to create a sequence of items.

> a <- c(1:10)
> a
 [1]  1  2  3  4  5  6  7  8  9 10
> b <- rep(a, times=2)
> b
 [1]  1  2  3  4  5  6  7  8  9 10  1  2  3  4  5  6  7  8  9 10
> b <- rep(a, each=2)
> b
 [1]  1  1  2  2  3  3  4  4  5  5  6  6  7  7  8  8  9  9 10 10
> c <- seq(-2, 2, by=.5)
> c
[1] -2.0 -1.5 -1.0 -0.5  0.0  0.5  1.0  1.5  2.0
Step 9. The content of a vector can be compared to another using basic operators

> x==x
[1] TRUE TRUE TRUE
> x==y
[1] FALSE FALSE FALSE
> x!=y
[1] TRUE TRUE TRUE
Step 10. The content of a Vector can be easily queried and modified. For this it is possible to use Index Vectors to subset some elements of an existing vector, using square brackets

> x
[1] 10 20 30
> x[3]
[1] 30
> x[3] <- 50
> x
[1] 10 20 50
> length(x)
[1] 3
Step 11. For Index Vectors of character strings, a “names” attribute can help identify components and query the data.

> dairy <- c(10, 20, 1, 40)
> names(dairy) <- c("milk", "butter", "cream", "yogurt")
> breakfast <- dairy[c("milk","yogurt")]
> breakfast
  milk yogurt 
    10     40
This might be useful to remember when you will manipulate data frames.

Discussion
Now try it yourself and discuss in the comment area below:

Question 1. Did you manage to create and manipulate Variables?

Question 2. Did you manage to create and manipulate Vectors?

Exercise
Let’s try it !

Question 1. Could you create 3 vectors:

a vector x containing the numbers 3, 10 and 30

a vector m containing the content of x repeated twice

a vector n containing two copies of x separated by a 0

Question 2. Is the content of m equal to the content of n?

Question 3. Note that you should also obtain a warning message because the 2 vectors are not of the same length. How can you check the length of both vectors?
