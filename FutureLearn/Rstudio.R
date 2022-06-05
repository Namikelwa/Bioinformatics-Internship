## How to Install R and RStudio

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


## Let’s Begin with R!

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



## How to Create and Manipulate Variables and Vectors in R

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



SOLUTION 1.
> x <- c(3, 10, 30)
> m <- rep(x, times=2)
> m
[1] 3 10 30 3 10 30
> n <- c(x,0,x)
> n
[1] 3 10 30 0 3 10 30
SOLUTION 2.
> m==n
[1] TRUE TRUE TRUE FALSE FALSE FALSE FALSE
SOLUTION 3.
> length(m)
[1] 6
> length(n)
[1] 7




## How to Create and Manipulate Lists and Data frames in R

Lists and Data frames in R

A List is an “object” in R that consist of a collection of other objects known as components.

Lists can have components of the same type or mode, or components of different types or modes. They can hence combine different components (numeric, logical…) in a single object.

A Data frame is simply a List of a specified class called “data.frame”, but the components of the list must be vectors (numeric, character, logical), factors, matrices (numeric), lists, or even other data frames. Other restrictions include the fact that the contents of a data frame must have the same length (vectors), or be of the same row size (matrices). A data frame can be considered as a simple matrix containing rows and columns, having potentially different modes and attributes.

We encourage you to find more information on Lists and Data frames in the following link: https://cran.r-project.org/doc/manuals/r-release/R-intro.html Let’s see together in this step how Lists and Data frames can be created and manipulated.

How to create and manipulate Lists
Step 1. We recommend you to work in the same working sub-directory that you created previously, using one of the following options

Before launching R

$ cd exerciseR
$ pwd
/Users/imac/Desktop/exerciseR
$ R
After launching R

> setwd("/Users/imac/Desktop/exerciseR")
> getwd()
[1] "/Users/imac/Desktop/exerciseR"
Step 2. Let’s create a simple List called L. We need to assign elements to this List.

> L<-list(dairy="milk",type="almond",form="liquid",contain.liter=c(0.5,1,2))
This will create a list called L with 4 elements.


> L<-list(dairy="milk",type="almond",form="liquid",contain.liter=c(0.5,1,2))
> 
> L
$dairy
[1] "milk"
 
$type
[1] "almond"
 
$form
[1] "liquid"
 
$contain.liter
[1] 0.5 1.0 2.0
Step 3. The function “length()” can allow you to easily retrieve the number of top level components of this List.

> length(L)
[1] 4
Step 4. Please note that components of a list are always numbered by default. If L is the list with 4 components we just created, then L[[1]] will be its first component, etc. We can also refer to the first entry of Component 4 independently, as L[[4]] is a vector itself and L[[4]][[1]] will refer to its first entry.

> L[[1]]
[1] "milk"
> L[[2]]
[1] "almond"
> L[[3]]
[1] "liquid"
> L[[4]]
[1] 0.5 1.0 2.0
> 
> L[[4]][[1]]
[1] 0.5
Step 5. More conveniently, we could also refer to the Component of a List by its name, instead of using its position number between double brackets. When using this option, you will need to provide the Component name by using the “$” symbol. As an example, you can use L$dairy or L[[1]] to refer to the first component of a List equally

> L$dairy
[1] "milk"
> L[[1]]
[1] "milk"
Step 6. Different Lists can be combined using the concatenation function that we used before, c(). This will result in an object of mode List also, because we gave this function arguments of the mode List.

> List.A <- list(dairy="milk", type="almond")
> List.B <- list(dairy="yogurt", type="frozen")
> List.AB <- c(List.A, List.B)
> List.AB
$dairy
[1] "milk"
 
$type
[1] "almond"
 
$dairy
[1] "yogurt"
 
$type
[1] "frozen"
Note. Please note that there are many more options to use the Lists. We only covered those that might be relevant to understand for the rest of the course. We encourage you to find more information on lists in the link given to you in the introductory part of this Step.

How to create and manipulate Data frames
Step 1. Now that we know how Variables, Vectors, and Lists can be created and manipulated, let’s use this knowledge to create sequentially a data frame called df.

> Name <- c("Lilly", "James", "Harry")
> Age <- c(30,  31, 11)
> Height <- c(168, 179, 139)
> Weight <- c(57, 69, 32)
> df <- data.frame (row.names = Name, Age, Height, Weight)
Step 2. Once created, the data frame can be called directly by simply typing its name.

> Name <- c("Lilly", "James", "Harry")
> Age <- c(30,  31, 11)
> Height <- c(168, 179, 139)
> Weight <- c(57, 69, 32)
> df <- data.frame (row.names = Name, Age, Height, Weight)
> df
      Age Height Weight
Lilly  30    168     57
James  31    179     69
Harry  11    139     32
Step 3. Additional information can be added to an existing data frame. We can create a new data frame containing the same names, to be able to make the correspondence (used as a primary key) and then combine both data frames using cbind(), a function used to combine objects (vectors, data frames,…) by columns.

> Name <- c("Lilly", "James", "Harry")
> Sex <- c("F", "M", "M")
> df.add <- data.frame(row.names = Name, Sex)
> df.add
      Sex
Lilly   F
James   M
Harry   M
> df.all <- cbind(df, df.add)
> df.all
      Age Height Weight Sex
Lilly  30    168     57   F
James  31    179     69   M
Harry  11    139     32   M
Step 4. The information added can be in the form of Factors that can be used to represent categorical data, and can help you using plotting functions later on. Let’s create again a new data frame (df.add.fact) with the information in the Sex vector added as a Factor, and combine both in a new data frame (df.all.fact)

> Name <- c("Lilly", "James", "Harry")
> Sex <- as.factor(c("F", "M", "M"))
> df.add.fact <- data.frame(row.names = Name, Sex)
> df.all.fact <- cbind(df, df.add.fact)
> df.all.fact
      Age Height Weight Sex
Lilly  30    168     57   F
James  31    179     69   M
Harry  11    139     32   M
Note 1. Note that we coerced the content of “Sex” to be a Factor by as.factor(). Factors are categorical data that can only take certain values such as “M” and “F”, which is the case of the field “Sex”. These distinct values are predefined and will be called Levels. This can be checked using the functions class() and levels()

> class(Sex)
[1] "factor"
> levels(Sex)
[1] "F" "M"
Note 2. At this stage, you will notice no difference. But using levels() you will be able to see how factors can now be recognised as such.

> levels(df.all$Sex) <- c("M", "F")
> df.all
      Age Height Weight Sex
Lilly  30    168     57   M
James  31    179     69   F
Harry  11    139     32   F
> 
> levels(df.all.fact$Sex) <- c("M", "F")
> df.all.fact
      Age Height Weight Sex
Lilly  30    168     57   M
James  31    179     69   F
Harry  11    139     32   F
To query the type of levels we can use “levels()” and to query the number of levels, you can use “nlevels()”

> levels(Sex)
[1] "F" "M"
> nlevels(Sex)
[1] 2
Step 5. Let’s see a set of useful functions to explore and manipulate a data frame

How many rows and columns are in the data frame df.all.fact ? you can use dim() to set or get the dimension of the data frame (rows, columns), or more specifically nrow() for the number of columns and ncol() for the number of columns.
> dim(df.all.fact)
[1] 3 4
> nrow(df.all.fact)
[1] 3
> ncol(df.all.fact)
[1] 4
2. What is the class of data in each column? Use sapply() which will output the result of a certain function to an object (here will output the classes in the data frame) or str() to display the structure of an object in R or all basic structures of a data frame (one line for each)

> sapply(df.all.fact, class)
      Age    Height    Weight       Sex 
"numeric" "numeric" "numeric"  "factor" 
> 
> str(df.all.fact)
'data.frame':	3 obs. of  4 variables:
 $ Age   : num  30 31 11
 $ Height: num  168 179 139
 $ Weight: num  57 69 32
 $ Sex   : Factor w/ 2 levels "M","F": 1 2 2
3. It is possible to subset or filter a data frame, as simply as we did it for Lists. For instance, let’s see here how to select one column or one row.

To select a column: [1] is column 1, [,1] is column 1 displayed as a vector

To select a row: [1,] is row 1

> df.all.fact[1]
      Age
Lilly  30
James  31
Harry  11
> df.all.fact[,1]
[1] 30 31 11
> 
> df.all.fact[1,]
      Age Height Weight Sex
Lilly  30    168     57   M
4. This is how to select a group of elements.

To select the element in column 1 and row 1: [1,1]

To select elements 1 to 3 in column 3

> df.all.fact[1,1]
[1] 30
> 
> df.all.fact[1:3,3]
[1] 57 69 32
5. To re-order data in a data frame, there are different options. We can use “order()”. Let’s try to re-order here based on the Height column

> df.all.fact
      Age Height Weight Sex
Lilly  30    168     57   M
James  31    179     69   F
Harry  11    139     32   F
> df.all.fact[order(df.all.fact$Height),]
      Age Height Weight Sex
Harry  11    139     32   F
Lilly  30    168     57   M
James  31    179     69   F
6. To filter the data, functions such as “unique()” and “sort()” can be used. This should remind you of the sort and uniq functions that can be used in Unix.

To obtain unique values of the column Age: unique(df.all.fact$Age)

To obtain sorted unique values of the column age: sort(unique(df.all.fact$Age))

> unique(df.all.fact$Age)
[1] 30 31 11
> 
> sort(unique(df.all.fact$Age))
[1] 11 30 31
Note. There are different types of data that can be considered, and treated differently according to their nature.

table explaining different types of data, qualitative and quantitative


Working on Data Frames - Exercises and Discussion
21 comments
Working on Data frames: here are some questions for you to practice

Try to answer the following questions, then compare your answers with your fellow learners in the comment area.

Question 1. Create a Data frame
Create a Data frame called df_fruits from the following Vectors:

A Vector called Fruits composed of Apple, Banana, Orange, Mango
A Vector called Price composed of 4, 3, 2, 8
A Vector called Nature composed of Local, Exotic, Local, Exotic as factors
Question 2. Check your Data frame
How would you verify that the newly Data frame called df_fruits

Is properly created
What are the levels of the factors in the column called Nature
How many different of these levels were created
Question 3. Manipulate your Data frame
How would you verify the number of rows and columns of the new Data frame called df_fruits?
What class of data is in each column?
How would you select all elements of the column called Price?


Solution 1.
> Fruits <- c("Apple", "Banana", "Orange", "Mango")
> Price <- c(4, 3, 2, 8)
> Nature<- as.factor(c("Local", "Exotic", "Local", "Exotic"))
> df_fruits <- data.frame(row.names = Fruits, Price, Nature)
Solution 2.
> df_fruits
 Price Nature
Apple 4 Local
Banana 3 Exotic
Orange 2 Local
Mango 8 Exotic
>
> levels(Nature)
[1] "Exotic" "Local"
> nlevels(Nature)
[1] 2
Solution 3.
1. Number of rows and columns
> df_fruits
 Price Nature
Apple 4 Local
Banana 3 Exotic
Orange 2 Local
Mango 8 Exotic
> dim(df_fruits)
[1] 4 2
> nrow(df_fruits)
[1] 4
> ncol(df_fruits)
[1] 2
2. Classes of data
> sapply(df_fruits, class)
 Price Nature
"numeric" "factor"
> str(df_fruits)
'data.frame': 4 obs. of 2 variables:
$ Price : num 4 3 2 8
$ Nature: Factor w/ 2 levels "Exotic","Local": 2 1 2 1
3. There are 2 Possible solutions depending on the desired output
> df_fruits[1]
 Price
Apple 4
Banana 3
Orange 2
Mango 8
> df_fruits[1:4,1]
[1] 4 3 2 8



## Reading Data from Files in R

Reading data from files

Introduction
Now that we have learned how to create, query and manipulate simple data frames, let’s see how the basic functions we covered can be useful with more complex datasets.

Indeed, you will often be willing to exploit the usefulness of built-in functions in R to manipulate large data sets. Although these data sets are organized as data frames, with information organized in rows and columns, it is obviously not easy to create it from scratch. Instead, there are specific built-in functions in R that allows us to import existing data contained in a file.

Let’s see together in this step how to read data from an existing file. We will also see how to simply query with basic functions in R this exiting data set.

We will use the iris dataset, another example of tab-delimited file such as the diamonds dataset that you have used before in this course. It contains information about 3 plant species (setosa, virginica, versicolor) and related measures about 4 features (Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) of these plants.

Note that this information is organized in a particular form:

Each column is named to indicate the features and species to consider
Each row contains a row label and information corresponding to features and species
If you want in the future to use another file you generated, make sure that no missing information is left empty. Instead, indicate it as NA (Not Applicable).

We encourage you to find more information on reading, querying and manipulating data frames from files in the following links: https://cran.r-project.org/doc/manuals/r-release/R-intro.html and https://rpubs.com/moeransm/intro-iris

How to import and read data from files in R
Step 1. We recommend you to work in the same working sub-directory that you created previously, using one of the following options

Before launching R

$ cd exerciseR
$ pwd
/Users/imac/Desktop/exerciseR
$ R
After launching R

> setwd("/Users/imac/Desktop/exerciseR")
> getwd()
[1] "/Users/imac/Desktop/exerciseR"
Step 2. There are 2 ways to use existing files and access the data sets they contain in R:

Option 1. If the file exists in your computer, use the “read.table()” function to read the file / data frame. To view part of its content, you can use many functions as the “head()” or “tail()” functions as in Unix

> Iris <- read.table("iris.txt")
> head(Iris)
  Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa
Option 2. R has also pre-built datasets that can be used for exercise purposes. We will use the “iris“ data set. This data set is available under the datasets package. To load a package under R, use the “library()” function. Once you loaded the datasets package, call the data set you want using the “data()” function

> library(datasets)
> data(iris)
> head(iris)
  Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa
To review all the pre-build data sets available in R (type q to quit)

> data()
Note 1. Note that we called the variable in which we imported the file “Iris”, whereas the data set called from R is named “iris”. You are free to choose the name of the variables you use in R, but if you call data sets from the existing datasets package you must use its proper nomenclature

Note 2. Note that once your data set of interest is loaded, all the commands and functions that we will use will be applicable to both “Iris” and “iris” equally. As an example, we will use “iris” in this Step.

How to make simple queries and data manipulation in R
Step 1. To view a summary statistics of the whole data set, use the “summary()” function. You can also view summary statistics of one of the variables using the “$” option we saw previously

> summary(iris)
  Sepal.Length    Sepal.Width     Petal.Length    Petal.Width   
 Min.   :4.300   Min.   :2.000   Min.   :1.000   Min.   :0.100  
 1st Qu.:5.100   1st Qu.:2.800   1st Qu.:1.600   1st Qu.:0.300  
 Median :5.800   Median :3.000   Median :4.350   Median :1.300  
 Mean   :5.843   Mean   :3.057   Mean   :3.758   Mean   :1.199  
 3rd Qu.:6.400   3rd Qu.:3.300   3rd Qu.:5.100   3rd Qu.:1.800  
 Max.   :7.900   Max.   :4.400   Max.   :6.900   Max.   :2.500  
       Species  
 setosa    :50  
 versicolor:50  
 virginica :50
> 
> summary(iris$Species)
    setosa versicolor  virginica 
        50         50         50 
> 
> summary(iris$Petal.Length)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.000   1.600   4.350   3.758   5.100   6.900 
Step 2. Let’s now query the names of columns using the “names()” function, and the data set content in terms of number of columns and rows, structure, etc…

> names(iris)
[1] "Sepal.Length" "Sepal.Width"  "Petal.Length" "Petal.Width"  "Species"     
>
> names(iris)
[1] "Sepal.Length" "Sepal.Width"  "Petal.Length" "Petal.Width"  "Species"     
> dim(iris)
[1] 150   5
> ncol(iris)
[1] 5
> nrow(iris)
[1] 150
>
> sapply(iris, class)
Sepal.Length  Sepal.Width Petal.Length  Petal.Width      Species 
   "numeric"    "numeric"    "numeric"    "numeric"     "factor" 
> str(iris)
'data.frame':	150 obs. of  5 variables:
 $ Sepal.Length: num  5.1 4.9 4.7 4.6 5 5.4 4.6 5 4.4 4.9 ...
 $ Sepal.Width : num  3.5 3 3.2 3.1 3.6 3.9 3.4 3.4 2.9 3.1 ...
 $ Petal.Length: num  1.4 1.4 1.3 1.5 1.4 1.7 1.4 1.5 1.4 1.5 ...
 $ Petal.Width : num  0.2 0.2 0.2 0.2 0.2 0.4 0.3 0.2 0.2 0.1 ...
 $ Species     : Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 ...
> 
Step 3. To query or manipulate this data set, it is possible to use basic operators in R

> setosa1 <- iris[iris$Species == "setosa",]
> head(setosa1)
  Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa
> nrow(setosa1)
[1] 50
Alternative option

> setosa2 <- iris[iris$Species %in% "setosa",]
> nrow(setosa2)
[1] 50
To select data related to the setosa species in which Sepal.Length > 5

> setosa3<- setosa2<-iris[iris$Species %in% "setosa" & iris$Sepal.Length>5,]
> nrow(setosa3)
[1] 22
Step 4. To avoid using operators, conditional subsetting is also possible with base functions in R that can ease the process and using the same principles. An example is the “subset()” function

To select only data related to the setosa species

> setosa.sub1 <- subset(iris, Species == "setosa")
> head(setosa.sub1)
  Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa
> nrow(setosa.sub1)
[1] 50
To select again data related to the setosa species in which Sepal.Length > 5

> setosa.sub2 <- subset(iris, Species == "setosa" & Sepal.Length > 5)
> nrow(setosa.sub2)
[1] 22
Exercise
Question 1. How would you check if the variables setosa1 and setosa.sub1 are equivalent?

Question 2. What is the structure of setosa.sub2 ?

SOLUTION 1.
> setosa.test <- setosa1 == setosa.sub1
> head(setosa.test)
 Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1 TRUE TRUE TRUE TRUE TRUE
2 TRUE TRUE TRUE TRUE TRUE
3 TRUE TRUE TRUE TRUE TRUE
4 TRUE TRUE TRUE TRUE TRUE
5 TRUE TRUE TRUE TRUE TRUE
6 TRUE TRUE TRUE TRUE TRUE
SOLUTION 2.
> str(setosa.sub2)
'data.frame': 22 obs. of 5 variables:
$ Sepal.Length: num 5.1 5.4 5.4 5.8 5.7 5.4 5.1 5.7 5.1 5.4 ...
$ Sepal.Width : num 3.5 3.9 3.7 4 4.4 3.9 3.5 3.8 3.8 3.4 ...
$ Petal.Length: num 1.4 1.7 1.5 1.2 1.5 1.3 1.4 1.7 1.5 1.7 ...
$ Petal.Width : num 0.2 0.4 0.2 0.2 0.4 0.4 0.3 0.3 0.3 0.2 ...
$ Species : Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 .


## Querying and Manipulating Data from Files Using Dedicated Packages in R

Querying and manipulating data from existing files
Introduction

Querying and manipulating data from files might require you to use advanced options. In R, packages have been developed for specific purposes such as data manipulation, data analysis, or plotting. Packages need first to be installed and loaded into R. Each package comes with a set of functions.

Let’s see together in this step how to query and manipulate data from an existing file using a dedicated package called dplyr. We will see how some of its functions can be helpful for more complex queries and manipulation of your data than basic R functions.

We will keep using the iris dataset, that you have used already.

We encourage you to find more information on reading, querying and manipulating data frames from files in the following links: https://cran.r-project.org/doc/manuals/r-release/R-intro.html and https://rpubs.com/moeransm/intro-iris

Import and read data from files in R
Step 1. We recommend you to work in the same working sub-directory that you created previously, using one of the following options

Before launching R

$ cd exerciseR
$ pwd
/Users/imac/Desktop/exerciseR
$ R
After launching R

> setwd("/Users/imac/Desktop/exerciseR")
> getwd()
[1] "/Users/imac/Desktop/exerciseR"
Step 2. You should be able now to call again the dataset we want you to work on, the iris data set.

To read the file from your computer

> Iris <- read.table("iris.txt")
> head(Iris)
  Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa
To read the file from the available data sets in R

> library(datasets)
> data(iris)
> head(iris)
  Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa
Using the dplyr package
Step 1. To use the dplyr package on the iris data set, we will need to call the package.

Because dplyr is not genuinely part of R, we will need to install it first using the “install.packages()” function

> install.packages("dplyr")
--- Please select a CRAN mirror for use in this session ---
Secure CRAN mirrors 
 
 1: 0-Cloud [https]
 2: Australia (Canberra) [https]
 3: Australia (Melbourne 1) [https]
……
76: USA (TX 1) [https]
77: Uruguay [https]
78: (other mirrors)
Once you select your CRAN mirror of interest (the closest to your location) as we did before, the installation will proceed. Once terminated, you will be able to load the package using the “library()” function

> library(dplyr)
Note 1. A package comes with specific functions that would not otherwise be recognized in R. Examples of basic verbs for data manipulation available with the dplyr package are “filter()”, “select()”, “mutate()”, “arrange()”, “rename()”, “relocate()”, “slice()”, “summarise()”. We will see how to use the first 3 verbs, but if you want information on other dplyr functions, or more advanced options, please refer to https://dplyr.tidyverse.org/articles/dplyr.html

Note 2. For all dplyr functions, as it will be the case for other packages in R, the first argument needs to be the data frame, also called tibble.

Step 2. Let’s see how you can now use the “filter()” function to filter specific data from this file, as we used the “subset()” function in base R. Let’s again filter only the data related to the Species “setosa”

Filtering the data after installing and loading the package

> setosa.filt <- filter(iris, Species == "setosa")
> head(setosa.filt)
  Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa
To check that the new variable setosa.filter generated contains only data related to setosa species, as it was also the case for the variable setosa generated with the “subset()” function in base R (previous Step9)

> nrow(setosa)
[1] 50
> nrow(setosa.filt)
[1] 50
To filter on multiple conditions: here based on setosa species having a Petal.Length smaller than 2, then > 2

> setosa.filt.pl2 <- filter(iris, Species == "setosa", Petal.Length < 2)
> head(setosa.filt.pl2)
  Sepal.Length Sepal.Width Petal.Length Petal.Width Species
1          5.1         3.5          1.4         0.2  setosa
2          4.9         3.0          1.4         0.2  setosa
3          4.7         3.2          1.3         0.2  setosa
4          4.6         3.1          1.5         0.2  setosa
5          5.0         3.6          1.4         0.2  setosa
6          5.4         3.9          1.7         0.4  setosa
> nrow(setosa.filt.pl2)
[1] 50
> setosa.filt.pl2 <- filter(iris, Species == "setosa", Petal.Length > 2)
> nrow(setosa.filt.pl2)
[1] 0
Note. Note that if the 2 variables are of the same name (data filtered based on > 2 or <2), the latter will replace the previous

Step 3. To select specific columns, the “select()” function can be very helpful To select specified columns that can be distant

> Iris.select <- select(iris, Sepal.Length, Petal.Length) 
> head(Iris.select)
  Sepal.Length Petal.Length
1          5.1          1.4
2          4.9          1.4
3          4.7          1.3
4          4.6          1.5
5          5.0          1.4
6          5.4          1.7
To select a group of consecutive columns

> Iris.select.2 <- select(iris, Sepal.Length:Petal.Width) 
> head(Iris.select.2)
  Sepal.Length Sepal.Width Petal.Length Petal.Width
1          5.1         3.5          1.4         0.2
2          4.9         3.0          1.4         0.2
3          4.7         3.2          1.3         0.2
4          4.6         3.1          1.5         0.2
5          5.0         3.6          1.4         0.2
6          5.4         3.9          1.7         0.4
Step 4. Now imagine you would like to add new columns to an existing data frame. The “mutate()” function can be used in the following example to add a new column called Test, with the information of whether Sepal.Length is greater than twice the size of Petal.Length (TRUE) or not (FALSE)

> test.col <- mutate(iris, Test = Sepal.Length > 2 * Petal.Length)
> head(test.col)
  Sepal.Length Sepal.Width Petal.Length Petal.Width Species Test
1          5.1         3.5          1.4         0.2  setosa TRUE
2          4.9         3.0          1.4         0.2  setosa TRUE
3          4.7         3.2          1.3         0.2  setosa TRUE
4          4.6         3.1          1.5         0.2  setosa TRUE
5          5.0         3.6          1.4         0.2  setosa TRUE
6          5.4         3.9          1.7         0.4  setosa TRUE
> tail(test.col)
    Sepal.Length Sepal.Width Petal.Length Petal.Width   Species  Test
145          6.7         3.3          5.7         2.5 virginica FALSE
146          6.7         3.0          5.2         2.3 virginica FALSE
147          6.3         2.5          5.0         1.9 virginica FALSE
148          6.5         3.0          5.2         2.0 virginica FALSE
149          6.2         3.4          5.4         2.3 virginica FALSE
150          5.9         3.0          5.1         1.8 virginica FALSE
Exercise
Question 1. Using this last example, how would you count the number of TRUE items in the newly created Test column?


SOLUTION 1.
> test.col.num <- filter(test.col, test == "TRUE")
> nrow(test.col.num)
[1] 50
You can double check the total number of rows in test.col and the complementary
number of FALSE items
> test.col.num2 <- filter(test.col, test == "FALSE")
> nrow(test.col.num2)
[1] 100
> nrow(test.col)
[1] 150


## Making Simple Plots in R

Introduction

The majority of us have long been in the habit of using Excel to make plots. This is still fine for some small data sets, but for large data sets, it is really complicated to load data in Excel and manipulate it. Moreover, Excel is prone to errors when manipulating the data when it comes to filtering, arranging, or querying complex data. This is without even mentioning the complexity of generating and arranging plots in a publishable way.

In a world where information is mainly visual, whether science, journalism or even social media publications, R can rapidly and efficiently meet our needs. We will see here how very basic functions in R can make your life easier !

We encourage you to find more information on plotting data frames from files in the following links: https://cran.r-project.org/doc/manuals/r-release/R-intro.html and/or https://rpubs.com/moeransm/intro-iris and/or https://hbctraining.github.io/Intro-to-R/lessons/basic_plots_in_r.html

Import and read data from files in R
Step 1. We recommend that you work in the same working sub-directory that you created previously, using one of the following options

Before launching R

$ cd exerciseR
$ pwd
/Users/imac/Desktop/exerciseR
$ R
After launching R

> setwd("/Users/imac/Desktop/exerciseR")
> getwd()
[1] "/Users/imac/Desktop/exerciseR"
Step 2. Import or load the iris dataset we want you to work on

From your computer

> Iris <- read.table("iris.txt")
From the available data sets in R

> library(datasets)
> data(iris)
Basics of plotting graphics in R
Introduction

The principle of plotting using R commands is to provide generally 2 main information types: (1) the data we want to use, and (2) preferences or options for display. This information should be provided as individual elements, that will be interpreted in R as arguments. They are interpreted as layers of information.

Simple basic graphics in R

Basic graph types found in typical spreadsheet software also exist in R, such as histograms, barplots, scatterplots or boxplots. These can be generated using commands or functions such as “hist()”, “barplot()”, “plot()”, “boxplot()” respectively. Many others exist, all as part of the ‘base’ graphics package of R, but we will only cover examples of graphs generated with “hist()” and “plot()”. For a full list of variants, use

> library(help = "graphics")
Making simple basic graphics in R
The structure usage is function(data, options) Whatever is specified after the command is called arguments. Each command or function has its own set of arguments, but they will all follow the same structure. Please note that:

Some plotting functions in R can be used with either a whole data set, or specific data from a data frame (such as “plot()”), but others need data to be specified (such as “hist()”)
you can also provide your data with no options. This will generate an automatic graph of the data and will use the default options of the command.
Histogram with function “hist()”

Step 1. Why choose histograms to represent your data? Generally because you want to show the distribution of numerical data. To see examples of graphs you can generate with “hist()”, use the function “example()”

> example("hist")
Step 2. Usage:

hist(x, …)

example of possible arguments (https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/hist)

hist(x, breaks =, freq =, probability =, include.lowest =, right , density =, angle , col =, border =, main =, xlim =, ylim =, xlab =, ylab =, axes =, plot =, labels =, nclass =, warn.unused =, …)

Note . We need to specify which specific data from the iris data set we want to represent. The x corresponds to the data to represent. The following arguments are generally parameters that impact the graphical output.

Step 3. Let’s generate a histogram of Sepal.Length using default (only the data is specified) or advanced arguments.

Default arguments will output a histogram in a simple format (default naming of axis, colors, font…), but note how the axes have been optimized for the data.

> hist(iris$Sepal.Length)
histogram of Sepal.Length using default arguments

Using advanced arguments can allow you to customize different features in your output. The following options will rename the x-axis (xlab), give a title to the graph (main), color the borders (border), color the bars (col), and modify the y-axis limits (ylim)

> hist(iris$Sepal.Length, xlab="Sepal Length", 
main="Histogram of Sepal Length", border="white", 
col="red3", ylim=c(1, 40))
same histogram using advanced arguments

Note. Colors can also be specified using the HEX (hexadecimal) color code. You can find more information on HEX color codes in https://www.color-hex.com/

Histogram of Sepal.Length with the same arguments as before, except that we will remove the borders and color the bars with the same red colour but using now the HEX color code

> hist(iris$Sepal.Length, xlab="Sepal Length", 
main="Histogram of Sepal Length", border=FALSE, 
col="#CD0000", ylim=c(1, 40))
same histogram but with removed the borders and colour of the bars red

Plot with function “plot()”

Step 1. As you can imagine, plot is a generic term to design a wide range of graphics. The function “plot()” allows us to create many different plots

> methods(plot)
Step 2. Usage:

plot(x, y,…)

example of possible arguments (https://www.rdocumentation.org/packages/ROCR/versions/1.0-11/topics/plot-methods)

plot(x, y, type=, main=, xlab=, ylab=, pch=, col=,…)

Note. Here type specifies the type of plot that can be generated and are of many types such as “p” (points), “l” (lines), “b” (both), “o”, (both overplotted), etc

Step 3. With “plot()”, you can either use by default the whole data set or specify which specific data from the iris data set we want to represent.

Scatterplot of the whole dataset

> plot(iris) 
Scatterplot of the whole iris dataset

Scatterplot using specified data (Sepal.Length vs. Petal.Length). Remember that these are continuous numeric data. You can test how to produce the same output with:

Option 1

> plot(iris$Sepal.Length, iris$Petal.Length)
Option 2

> plot(Petal.Length ~ Sepal.Length, data=iris)
Option 3

> plot(Petal.Length ~ Sepal.Length, iris)
Option 4

> with(iris, plot(Sepal.Length, Petal.Length))
Scatterplot using specified data (Sepal.Length vs. Petal.Length)

Scatterplot using specified data and options. We will shape the points with pch, change their size using cex, and their color using col

> plot(iris$Sepal.Length, iris$Petal.Length, 
main="Sepal vs Petal Lengths", xlab="Sepal.Length", 
ylab="Petal.Length", pch="*", cex=2.0, col="red3") 
Scatterplot using specified data and options.

Scatterplot using specified data and options to change the background color and margin sizes with “par()”, a function used to specify general graphical parameters such as bg (background color), or mai (margins in inches for bottom, left, top and right)

> par(bg="lightgrey", mai=c(2,1,2,1.5))
> plot(iris$Sepal.Length, iris$Petal.Length, 
main="Sepal vs Petal Lengths", xlab="Sepal.Length", 
ylab="Petal.Length", pch="*", cex=3.0, col="red3")
Scatterplot using specified data and options to change the background color and margin sizes 

Note. You can quit these global graphical options by using “dev.off()”, or closing the graphical display

Saving a plot
By default, any plot you generate will be displayed in your graphic device window. To save a plot, you will have different options.

Option 1. First choose the output format (such as jpeg, png, pdf…), name your plot, generate it, then escape by closing the file. You will find the saved file in your working directory.

To make and save a file using default options

> pdf('test_hist.pdf')
> hist(iris$Sepal.Length)
> dev.off()
To make and save a file using advanced options, such as width and height (in inches)

> pdf('test_hist.pdf', 7, 10)
> hist(iris$Sepal.Length)
> dev.off()
Option 2. If you already generated your plot, and forgot to create the output file first, you can still use the “dev.copy()” command, with both the default or advanced options

> hist2(iris$Sepal.Length)
> dev.copy(pdf,'test_hist2.pdf')
> dev.off()
Note. These options will work with any OS (Linux, Mac, Windows). Some OSes offer the possibility to save the graphic window that opens with the “Save” or “Save as” option.

Exercise and Discuss - Working on Simple Plots
21 comments
Practise simple plotting in R

Question 1
Based on the same plotting principles you saw with the “plot()” function, can you draw a boxplot of the whole iris dataset?

Question 2
Based on the same plotting principles you saw previously, can you draw a boxplot of the Sepal.Width (x-axis) and Petal.Width (y-axis) from the iris dataset?

Question 3
Can you add to the previous boxplot the title “Boxplot Petal.Width vs. Sepal.Width”, and colour the boxplot using the “LightBlue” colour with its HEX colour code, and colour borders in “DarkSlateGray” with its HEX colour code?

Hint: you can browse websites dedicated to color code conversions to find the correspondence between a colour name and its HEX colour code, there are many available.
Please try to answer the questions yourself first and then compare the results with other learners. When you have tried the exercise, you can find solutions in the download area.

Solution 1.
To draw the boxplot, use the following
> boxplot(iris)
Solution 2.
To draw the boxplot, use the following
> boxplot(Petal.Width ~ Sepal.Width, data=iris)
Solution 3.
To draw the boxplot, use the following
> boxplot(Petal.Width ~ Sepal.Width, iris, col="#add8e6", border="#2f4f4f", main="Boxplot Peal.Width
vs. Sepal.Width")
Expected output:



## What is RStudio
[00:00:00.47] [MUSIC PLAYING]
[00:00:07.52] At this stage of the course, you should be able now to easily create variables, or
read and manipulate data frames, as well as creating simple graphs with ease. This is because all
what we covered in basic R will be helpful for this upcoming part of the course where we will be
working with RStudio. So for the next two steps, let's start first seeing together what RStudio is
and then how to use RStudio.
[00:00:33.29] You should have now properly installed RStudio from the first steps of the course.
So open it by simply clicking on it. And let's have first a general presentation of RStudio. And
then we will see together how the RStudio interface works.
[00:00:49.88] So RStudio is a free and open source integrated development environment, or IDE,
for abbrv. RStudio runs on all the major operating systems, such as Windows, Linux, and Mac.
What you are seeing now in the screen is the RStudio interface that should contain generally four
different areas, or quadrants, each devoted to give you access to a certain type of information.
Now, what are all these areas?
[00:01:21.06] So let's start with the bottom left quadrant. So this one is, by default, the console,
as you can see it written here. And it is what reproduces the exact same terminal access to R that
we used before in R, except of course, that it will open directly in our session. You can see here
that you have the greater than, or superior sign as a prompt. But this area also allows you to
access other resources, such as your main terminal here through the terminal tab. But we will be
using the console for now.
[00:01:58.52] Now, the bottom right quadrant, which is this one, has many functions. It is what
allows you to access your working directory. So if I click here on Files, you will be able to
access either your working directory or the directories through that Files tab. You can also create
a new folder. You can delete or rename folders. And many more actions you can do here.
[00:02:25.10] It is also where you will be able to view the plots you generate through the Plots
tab. So of course, we don't have any plots now. But once you generate them, you will be able to
see them through this tab.
[00:02:39.38] You can also from instal and load packages through the Packages tab. And you
can also have access to other information, such as the Help tab here. So for example, it can allow
you to interrogate a certain function for its usage.
[00:02:56.31] So let's say that you want to click, for example, on a package. So these are all the
packages that you have installed. So some of them are installed by default. Some others will be
installed because you want them to be there. So we will see examples of that. 
[00:03:12.06] But here, let's imagine that we want to see what this base package is. So if you
click on it, you will have the help of the documentation and the help for the package that is called
a call base.
[00:03:30.10] Click and go to the top right quadrant now. It is by default the environment. And it
is where you will find information on the objects you're working with, such as the variable that
you are generating while you work. As long as you work here in your console, you will be seeing
your variables generated in here. So you can easily remember them.
[00:03:54.01] And I loaded for you here as an example the iris dataset. This quadrant should, of
course, be empty when you open RStudio by default. But just for the sake of giving you an
example, I just loaded the iris dataset before. So from this quadrant, you will also be able to
access the history of all your commands from this area when clicking on this tab.
[00:04:19.21] The final quadrant, which is the top left quadrant, is what we call the source area.
And it is where you will generally view the source content of a script, for example. It could be an
existing file that you directly open from files existing here. Or it could also be a new one that you
are creating, and where you can add comments that are successful as long as you are generating
commands, or typing commands in the console. The fact that commands are placed into files
would allow you, of course, to manipulate them very easily later on, and to retrieve them very
easily.
[00:05:01.49] So you can access all these areas of quadrants at the same time, as you can see it
here. Or you can simply click to reduce or expand some of them to reorganise your interface
according to your preferences, of course. So let's say for example, that I want to view mainly the
console and not necessarily the script. So I just click on the small icon, and I will have the
console taking this area.
[00:05:30.44] I can also reduce it completely, or I can restore the first version of it. I can also
click on any of these quadrants to do exactly the same. I am not going into the details of
modifying other basic settings, such as the background colour or the default working directory.
But there are many other settings you can adjust. So know that it is feasible from the Preferences
of RStudio. So if you open your Preferences, you will be able to access a lot of different settings,
such as the Default Working Directory, for example, or the appearance of your screen.
[00:06:10.69] In short, if we can summarise what we saw in this first video, the RStudio
interface is very intuitive. It is organised in a way that allows you, as a user, to clearly view in
one single interface many information, such as the code you are using here, the commands you're
writing or executing here, the variables you are creating, or the graphics you're generating. So it
basically offers the possibility to interact with an R based environment through a user interface
solution such as simply clicking on File, for example, to open it, without having to write or code
in R to do so.
[00:06:51.88] In the next video, or step, we will start learning together how to use RStudio. 

## Bioinformatics - Lets begin with RStudio
[00:00:07.12] In this new step let's see how to start working now with RStudio and how the
RStudio interface is easing a lt of the process of working in R-based environment again. As an
example, when we see together in the first part of the video, how to open a new project, and then
in the second part, how to instal and use packages. So opening a new project before you start
work will allow you to do many things.
[00:00:34.26] First, it will allow you to keep track of a certain analysis you did, so it will allow
you, for example, to retrieve very easily, substeps that you've been through and that worked for
you. And it will also allow you to keep a certain organisation in your files so you can easily work
on different projects and each project will have its own script or set of scripts.
[00:01:00.09] So how to do that? There are many ways to open a new project. For example, in
the top right of your screen, you have a blue icon here that is indicating which project you're
working on. We have no projects right now, so it is indicated as none, so let's create a new one.
If you pull down the project's options, you will be given the option to create a new one from
here, so let's click on a new project. .
[00:01:29.15] You will be asked where you want to create it. Let's select, for example, a new
directory, and we want to create it, let's select New Project. And let's name it, for example, I will
name it Project Test, and where we want to put it is the exercise R folder. So let's continue using
that one. And it is the same that we were using for the R sessions.
[00:01:59.94] Then click on Create Project. So you will now see a new path appearing here,
which says that we are working under this project that we just created. So by default, it will open
the console again in its bigger format. You can just simply click on it to reduce it and you can
also resize it. That's absolutely no problem to do that.
[00:02:29.55] You can then create, also, new folders and files under this newly created project to
contain scripts or data, for example, that will belong to this specific project. This is kind of a
good practise that we generally use, because you will have dedicated scripts and dedicated data
for each project, and that you can easily retrieve them.
[00:02:52.68] So you can also load an existing script, as it is the case here with the previous
script I created with some very basic commands that we saw in our previous R steps. The script
is called "Test R commands," so I put it in the same exercise R. It is called "Test R commands."
So imagine that this is an existing script where you put interesting commands that you want to
retrieve. If I just simply click on it, it will open the different commands I put in the script in here
directly, in your source quadrant.
[00:03:30.74] A good practise is to leave notes for yourself in your scripts. So basically,
everything that is written following the hashtag here will be a comment that you leave for
yourself or, of course, for someone else that will be using your script. So this allows you to very
easily use a script that you generated some time ago, because you will retrieve indications or
explanations for each step easily. 
[00:03:55.19] For example, you can see here that I indicated how to set my working directory for
where I want to be here, and how to check that I am in the right place after that using that getwd
command. I guess you remember these commands from the previous R steps.
[00:04:12.39] So these comments allow me to remember what the command is doing. And you
can also have comments written at the end of the line. It's perfectly fine. This will not be read as
part of the script as long as it is written after a hashtag.
[00:04:30.38] So if you want to run a command from a script, you can then place yourself either
at the end of the command or you can just select that command and then click on Run. You just
hit that Run button at the top-right of the source quadrant. And you can see that you will be
redirected to where you want to be. OK?
[00:04:56.81] And you can see that the path changed in here. So if I want to go again in Project
Test, I can just basically type the same command, and slash, and I will select for Project Test. So
if I do it, I can do it from here. I will select the Project Test, and then I can again run that
command from here, and it will redirect me to that project list.
[00:05:30.40] So before we move forward, one very general comment I would like to make
about using RStudio is to tell you to please remember that when you end an RStudio session
using the Q-command, as we saw it, you will generally be asked to save your workspace, as we
saw it together. This can also be changed by accessing the settings and ask RStudio to register
automatically your session in a certain directory.
[00:06:02.31] Now, because we created a new project that we want to work from, when you
open again an RStudio session, you can either open the app directly and then use the setwd to
access that specific project, or you can simply go to your desktop and to the folder that you're
working on, and then when you select this Project Test and click on it, this will open a new
RStudio session with that new project that is called Project Test opened directly. So you can
notice that here we have an environment that is empty, whereas in the first one, we had the Iris
data set that was pre-loaded.
[00:06:49.62] So again, we can resize our console and continue working on that. So if I type now
getwd here as a command-- this is by default kind of a new script, so if I click on Run, it tells me
that I am definitely in that Project Test. So I can also open by going into my exercise R. I can
open again my commands, and everything is again accessible from that same options.
[00:07:24.04] All right. Now, the second example we want to see together, and that, again, will
show you how easy it is to use RStudio, is for installing packages. When you are in RStudio, you
have different possibilities to instal packages. The first one is to write the command as we saw it
in our steps again-- so basically writing instal packages in the console, and then simply typing
the name of the package in quotes, as you can see it here.
[00:07:52.27] For dplyr or ggplot2 packages that I put you here as examples, because these are
two packages that we will be needing. So you remember the dplyr from the previous step that we
used. And for the next steps, we will be using ggplot2 a lot. 
[00:08:06.97] I could hit the Run button from here directly, but I want to show you what happens
if you type a command in the console directly. In fact, it opens a pop-up window. So if I start
typing instal packages, it opens that first pop-up window to give me all the possible commands.
And once I click on the one I want, it will open another pop-up window to show me the usage of
the command. And this can be very helpful sometimes if you forgot some usage of this
command. And this is also possible with other commands you type in RStudio, of course. It is
not something that is particular to instal packages.
[00:08:46.84] So let me just imagine that I'm trying to instal dplyr. You need to remember that,
by default, you have some packages that are pre-installed. And specifically, if you installed a
newer version of RStudio, there are many that should be pre-installed, including here the dplyr
package. They are organised in an alphabetical order. So in fact, you can easily retrieve any
package that you're looking for.
[00:09:14.45] So this shows me that I don't need to instal the package again. This is already done
from RStudio. However, if I look for ggplot2, I'm not finding it here at all. And this is totally fine
if you don't find it by default from the beginning, simply because it is not pre-installed as a preinstalled package in RStudio. So let's launch the installation for ggplot2.
[00:09:52.50] So you will see a stop sign in here saying that you will not be able to use your
console while it is installing. So you have to wait a few moments for it to be installed. So here,
the connection was quite rapid. And I did try it before, so it is quite rapid. And then I removed it
to show you this example. But it might take a couple of minutes, in fact, to instal it, with a lot of
different codes that is written, as long as you go for installing the package. Don't worry. Don't
touch anything. Leave it until it will give you the prompt again, that creator sign that is particular
to R sessions.
[00:10:36.48] There is another option to instal a package from RStudio. You can just go to the
package, hit the Instal button, and then it will allow you to search on CRAN for your package of
interest. If you remember, CRAN is a repository from which the packages can be accessed, and
that we used previously in this course.
[00:10:55.74] Then you will simply have to type or select the package you want and wait for the
installation to proceed. So if I start typing ggplot2, it will show me different options. I will only
have to select ggplot2 and then click on the Instal button. It will do exactly the same as it was
done here. So I will click on Cancel, but you can Instal it in here.
[00:11:17.85] All right. Now that we've installed the packages we want, before being able to use
them, we will need to load them into our session, as we saw it in our steps. So you have, again,
two ways of doing that in RStudio. One would be to do it directly from the console using the
library function, as we saw it previously in this course. The second option you have is to open
the available packages and simply click on the package name to select it.
[00:11:46.21] Now, because I showed you how to use the console to Instal the package, let me
show you how to easily load the package from the packages tab directly. So you basically just
have to select your package of interest. You can see now ggplot2. That is installed. And then 
simply click on it, and you can see that automatically it recognises that click as being your
willingness to load that package. And it will write library for you to call that package. And it is
as simple as that.
[00:12:21.06] So one important thing to remember, however, is that you will have to Instal the
package you want to use only once. You need to do it only once. When you will open new R
sessions, ggplot2 will be pre-installed now. However, you will need to load the packages you
want each time you open a new RStudio session.
[00:12:44.14] So if we take again the example of dplyr, dplyr was pre-installed, but it was not
loaded. So if I want to use dplyr, I just simply click on it, and it will load that package. So make
sure to keep that in mind. Otherwise, you will not be able to use the functions that comes with a
package.
[00:13:04.20] However, you might have noticed that some packages that are previously preinstalled are also pre-loaded-- such as, for example, the data sets package-- you remember that
was the one we used to select the Iris data set-- and also the graphics package that are pre-loaded,
and that we use in this R session
[00:13:25.89] So we can take some time now to get familiarised with this display in RStudio
before moving forward in the upcoming steps. But you should have noticed that it is very easy
and simple to use. And you have many things that you can do with either script-based options, as
in R, or interactive options, by just clicking and selecting things as you would do it classically in
your machine. Now you're ready to go with RStudio. 


## Data Visualisation Packages and Principles - a Focus on ggplot2

As we saw, much of what is done with R can also be done in RStudio.

Data Visualization in RStudio
Introduction

As we saw, many things you do with R can also be done in RStudio. We would advise you to take some time to test in RStudio a few commands we used previously in R, and put in a Script what you want to use later. We would also recommend that you take some time to get familiarized with the RStudio display and test the options you can use through the user interface. This will be useful for the upcoming steps, where we will need to install packages, load them, use them, and generate plots using ggplot2.

Data visualization packages in R/RStudio
For the rest of the Steps, we will see how to use a dedicated package in RStudio for data visualization, called ggplot2.

Many packages exist that are devoted to data visualization.

The base graphics comes for example with the graphics package, and allows you to generate simple plots and then to improve aspects of the plot possibly through series of functions (we saw a quick example with “par()” and “plot()”).

The lattice plotting system is implemented with other packages such as the lattice package, which supports generating trellis graphs. It is generally used with a single function call that would specify all graphical parameters, which allow R to automatically compute the necessary graphical display.

With the ggplot2 package, we are almost combining both concepts. This package is based on using the grammar of graphics concept. It became so popular, that many other data visualization packages can complement it or use it, or are based on the same concepts. Examples of other packages are ggforce or ggvis.

It is also possible now to create more advanced interactive graphs in R/RStudio using packages such as Plotly or Shiny.

The grammar of graphics: basics of ggplot2
The grammar of graphics is the concept of using a particular grammar or language to specify and create certain statistical and graphical displays for data visualization. Applied to ggplot2, the grammar of graphics is implemented in a layered approach, using layers of information (data complemented with statistical or graphical information) to build up step by step to a final display of a graph. Layers complement each other in ggplot2 with different information types such as the aesthetics, the geometries, the faceting, the scales, the themes, the coordinates, the labels, and many others.

Here are explanations of some of these layers :

 	 
DATA	The data you want to plot
AESTHETICS	The graphical properties of data on the plot (x, y, …)
GEOMETRIES	The graphical elements that determines the visual display of a plot (point, line, bar, area, …). Each geometric object is related to specific aesthetics. For example, a geometric object “point” is related to aesthetics shape, size, color and position
FACETING	reorganizes the variables of a data into subsets with a certain graphical arrangement of elements
SCALES	defines how data is mapped as related to aesthetics (item colors according to their class, …)
THEMES	customizes non-data related display by defining options not directly related to the data itself
COORDINATES	customizes the coordinates
LABELS	sets plot title, legend or axis names
layers of grammar presented as parallel layers of different elements

## Data Visualisation with ggplot2 - Setting Data, Aesthetics and Geometries

Making data visualisations in RStudio.

Data Visualisation with ggplot2
Introduction

Let’s load our data set of interest, install and load all the packages we need, and start making data visualizations in RStudio.

Good practice

A good practice is to load the packages you need before starting your analysis. It is also recommended to write the packages you need in the script you prepare for a project. This is a list of convenient packages to use with ggplot2

> library(ggplot2)
> library(RColorBrewer)
> install.packages("viridis")
> library(viridis)
Note. Another widely used package in data science is called tidyverse, and is a collection of packages including ggplot2, dplyr and many other helpful resources. It can be worth trying to use it on your own after this course.

Setting your working directory in RStudio
Step 1. We recommend you to work in the Project folder Project_Test that we created previously, either by clicking directly on the Project_Test or using the following command

> setwd("/Users/imac/Desktop/exerciseR/Project_Test")
> getwd()
[1] "/Users/imac/Desktop/exerciseR/Project_Test"
Step 2. As a reminder, you can create a specific script file to write your commands and related comments.

Setting your data
Step 1. Import or load the iris dataset we want you to work on in RStudio. All options are identically accessed in R, but the two final options are particular to RStudio.

From your computer, if you placed the iris dataset file in your working directory
> Iris <- read.table("iris.txt")
From your computer, if the iris dataset file is in the parent folder exerciseR
> Iris <- read.table("/Users/imac/Desktop/exerciseR/iris.txt")
From the available data sets in R
> data(iris)
From the “Import Dataset” tab in the Environment, by selecting the correct file with its type and parent folder.

From the File menu, by choosing the “Import Dataset” option.

Note 1. Be careful to choose the “iris” dataset as “Iris” would here correspond to the same data set but with changes that could impede the rest of the commands.

Note 2. As other functions in R, the “read.table()” function has different options that you can view in the following link, which also shows you other functions used to import data from other file formats (for example with the “read.csv()” function to read “.csv” files). https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/read.table.

Step 2. You can also display and work on specific data chosen from the iris data set

> iris_length <- iris %>% select(Sepal.Length, Petal.Length)
> head(iris_length)
  Sepal.Length Petal.Length
1          5.1          1.4
2          4.9          1.4
3          4.7          1.3
4          4.6          1.5
5          5.0          1.4
6          5.4          1.7 
Setting Aesthetics and Geometries
Step 1. Let’s use basic layers to plot Petal.Length vs. Sepal.Length. With ggplot2, “aes()” specifies aesthetics for x and y-axis, and “geom_point()” generates a scatterplot

> ggplot(data = iris,aes(x = Sepal.Length, y = Petal.Length)) + 
geom_point()
Note 1. Here is an example of how you should see the output in your “Plots” area in RStudio. Note that you have the possibility to save your plot using the “Export” button, with options related to file formats. Other R options we saw for saving plots remain possible.

Note 2. We will not show the whole area again, but remember that the plots you generate will appear here.

Step 2. Using the same previous plot options, let’s color the points according to the Species

> ggplot(data = iris, aes(x = Sepal.Length, y = Petal.Length, color = Species)) + geom_point()
scattergraph with different species in different colours

Step 3. There are other possible shorter ways for generating this same output

> ggplot(iris, aes(Sepal.Length, Petal.Length, color = Species)) + 
geom_point()
Note. However, for the sake of clarity, we will mainly keep the full details such when using data, x and y to ease the understanding

scattergraph mono colour species

Step 4. It is possible to create a variable with your base aesthetics and then simply call it to apply other layers. The following will create the same output as the previous graph

> key <- ggplot(data = iris, aes(x = Sepal.Length, y = Petal.Length, , 
color = Species))
> key + geom_point()
Step 5. Different geometries can also be used to complement each other. Here “geom_smooth()” adds a trend line and area to the points

> key + geom_point() + geom_smooth()
scattergraph with trend line

Step 6. You should have noticed how geometries are here added with default options. Each has a set of options, such as removing the trend area in the following with se=FALSE

> key + geom_point() + geom_smooth(se=FALSE)
scattergraph with smooth trends

Step 7. You can easily change the points size, shape and colour from “geom_point()” options, but see how it affects the display: if you force one colour, you will not have any more colors by Species, even if they are required in the key variable

> key + geom_point(size=4, shape=15, color="red3")
scatter graph with red dots

Step 8. Or the size, shape and color as dependent now on Sepal.Length values from aes

> ggplot(data = iris, aes(x = Sepal.Length, y = Petal.Length, 
color = Sepal.Length, size = Sepal.Length)) + geom_point()
scatter graph with blue dots

Note. We used here the default ggplot2 colors, but we will see later on how to use other color palettes

Other Functions and Plots
Step 1. Remember that we are only covering here the “ggplot()” usage, but other possibilities exist to generate the same output as in Step 6 of this Article, such as “qplot()” which is used to generate quick plots with ggplot2

> qplot(Sepal.Length, Petal.Length, data = iris, color = 
factor(Species)) + 
geom_point() + 
geom_smooth(se=FALSE)
Step 2. Generating different plots will require different geometries

Boxplot with default options
> ggplot(data = iris, aes(x = Sepal.Length, y = Petal.Length, 
color = Species)) + geom_boxplot()
box plot with default options

Bar plot with default options
> ggplot(data=Iris,aes(x=Sepal.Length)) + geom_bar()
bar plot with default options

Or more complex ones even with default options such as Density plot
> ggplot(data=Iris,aes(x=Sepal.Length, y = Petal.Length)) + 
geom_density_2d_filled()
density plot with default options

Step 3. An important thing to remember is that each plotting functions comes with its own set of option, that might not work for other functions. Let’s see how to generate and modify histograms

Default options
> ggplot(data=Iris,aes(x=Sepal.Length)) + geom_histogram()
black and white histogram

Filling histogram colurs by Species. Note how calling the colour option is different here
> ggplot(data=Iris, aes(x=Sepal.Length,fill=Species)) + 
geom_histogram()

colourful histogram

Use binwidth option with histograms
> ggplot(data=Iris,aes(x=Sepal.Length,fill=Species)) + 
geom_histogram(binwidth = 0.05)
binwidth histogram in colours

Note. A wide range of different plots can be generated with ggplot2 such as Bar plots, Boxplots, Violin Plots, Density Plots, Area Charts, Correlograms…and many many more !

Exercise and Discuss - Working on Simple Plots
21 comments
Practise simple plotting in R

Question 1
Based on the same plotting principles you saw with the “plot()” function, can you draw a boxplot of the whole iris dataset?

Question 2
Based on the same plotting principles you saw previously, can you draw a boxplot of the Sepal.Width (x-axis) and Petal.Width (y-axis) from the iris dataset?

Question 3
Can you add to the previous boxplot the title “Boxplot Petal.Width vs. Sepal.Width”, and colour the boxplot using the “LightBlue” colour with its HEX colour code, and colour borders in “DarkSlateGray” with its HEX colour code?

Solution 1.
To draw the boxplot, use the following
> boxplot(iris)
Solution 2.
To draw the boxplot, use the following
> boxplot(Petal.Width ~ Sepal.Width, data=iris)
Solution 3.
To draw the boxplot, use the following
> boxplot(Petal.Width ~ Sepal.Width, iris, col="#add8e6", border="#2f4f4f", main="Boxplot Peal.Width
vs. Sepal.Width")


## Data Visualisation with ggplot2: Setting Facets and Scales
4 comments
Now that we learned how to choose our data, and apply layers of aesthetics and geometries, let’s explore other possible layers of information such as facets and scales.

Good practice

A good practice is to load the packages you need before starting your analysis. It is also recommended to write the packages you need in the script you prepare for a project. This is a list of convenient packages to use with ggplot2:

> library(ggplot2)
> library(RColorBrewer)
> library(viridis)
> install.packages("ggsci")
> library(ggsci)
Note that we need an extra package compared to the previous step.

Setting your working directory in RStudio
Step 1. We recommend that you work in the Project folder Project_Test that we created previously, either by clicking directly on the Project_Test or using the following command

> setwd("/Users/imac/Desktop/exerciseR/Project_Test")
> getwd()
[1] "/Users/imac/Desktop/exerciseR/Project_Test"
Step 2. As a reminder, you can create a specific script file to write your commands and related comments.

Setting Faceting
Step 1. Here is an example of how to facet the output by Species, using “facet_wrap()”, that requires the facets argument to be specified, i.e. here the Species

> key <- ggplot(data = iris,aes( x= Sepal.Length,y = 
Petal.Length, color = Species))
> key + geom_point() + geom_smooth(se=FALSE) + 
facet_wrap(~Species)
faceted output of species

Step 2. Other arguments with “facet_wrap()” give the possibility of fitting the y-axis scales to the values in order to optimize the output

> key + geom_point() + geom_smooth(se=FALSE) + 
facet_wrap(~Species, scale='free_y')
faceted output of species optimised

Setting Scales
Scales for positions and axis

Arguments for continuous x and y aesthetics are by default “scale_x_continuous()” and “scale_y_continuous()”. Variants include reversing order or transforming to a log scale. Please see usage in https://ggplot2.tidyverse.org/reference/scale_continuous.html

It is also possible to plot discrete variables using “scale_x_discrete()” or “scale_y_discrete()”

Step 1. To set x-axis limits using “scale_x_continuous()” with limits option

> key + geom_point() + geom_smooth(se=FALSE) + 
facet_wrap(~Species, scale='free_y') + 
scale_x_continuous(limits = c(1, 10))
faceted output with scaled x and y axis

Step 2. To reverse the x-axis using “scale_x_reverse()”

> key + geom_point() + geom_smooth(se=FALSE) + 
facet_wrap(~Species, scale='free_y') + 
scale_x_reverse()
faceted output with reversed x axis

Scales for colours, sizes and shapes
Continuous colour scales can be specified using many options. Some are already pre-installed in RStudio, such as the RColorBrewer, but other specific colour palettes can be easily installed, loaded and used. Usage for “scale_colour_gradient()” or “scale_fill_gradient()” as examples can be found here https://ggplot2.tidyverse.org/reference/scale_gradient.html

Step 1. Scales can be manually set by choosing specific colours, sizes and shapes

> key + geom_point() + geom_smooth(se=FALSE) + 
facet_wrap(~Species, scale='free_y') + 
scale_shape_manual(values=c(3, 16, 17)) +
scale_size_manual(values=c(2,3,4)) + 
scale_color_manual(values=c('#669999','#a3c2c2', '#b30059'))
scales set by choosing colours graph

Step 2. Scales can be set using existing colour palettes from the RColorBrewer package

> key + geom_point() + geom_smooth(se=FALSE) + 
facet_wrap(~Species, scale='free_y') + 
scale_color_brewer(palette="RdYlBu")
Note. We use “scale_color_brewer()” to customize colours for lines or points, whereas we would use “scale_fill_brewer()” for filling colours of area, histogram bars, boxplots, etc

graph, scales set using RColorBrewer

Note. RColorBrewer palettes can be consulted with

> display.brewer.all()
RColorBrewer palette

Step 3. Scales can use different options with other color palettes from the viridis package

> new_key <- key + geom_point(aes(color = Species)) + 
geom_smooth(aes(color = Species, fill = Species), method = "lm") + 
facet_wrap(~Species, scale='free_y') 
 
> new_key + scale_shape_manual(values=c(3, 16, 17)) +
scale_size_manual(values=c(2,3,4)) + 
scale_color_viridis(discrete = TRUE, option = "D") + 
scale_fill_viridis(discrete = TRUE)
Note. We changed options for “geom_point()” and “geom_smooth()” to show you some other possible display variations

graph, scales in different options with other color palettes from the viridis package

Step 4. Scales can use packages designed to offer color palettes taken from sources such as highly accessed journals. Examples are “scale_color_npg()” (from Nature Publishing Group), “scale_color_lancet()” (from The Lancet journal) or even “scale_color_tron()” (from the film “Tron: Legacy”). Remember that scale_color functions have their scale_fill counterparts

> new_key + scale_color_tron() + scale_fill_tron()
Graph, scales use tron() packages designed to offer color palettes taken from sources

The particular case of missing values
Step 1. Missing values (NA) can exist in any data set, and need to be taken into account when plotting data. Let’s use this very simple data frame containing NA values

Plotting with default colors in ggplot2. By default, a grey colour will be used for NA
> df_test <- data.frame(x = 1:10, y = 1:10, 
z = c(1, 2, 3, NA, 5, 6, 7, NA, 8, NA))
> plot_test <- ggplot(df_test, aes(x, y)) + 
geom_tile(aes(fill = z), size = 10)
> plot_test
Plot with default colours

We can ask to have no colour of NA values
> plot_test + scale_fill_gradient(na.value = NA)
graph, no colours on missing values

Or color NA values in a chosen colour, such as “red3” here
> plot_test + scale_fill_gradient(na.value = "red3")
graph, red colour on no missing value

Or use another colour palette, instead of default ggplot2 colours, and you will still be able to specify the NA values. Because we need many colors in this palette we use “scale_fill_gradientn()” instead of “scale_fill_gradient()”
> plot_test + scale_fill_gradientn(colours = viridis(7), na.value = "white")
graph, viridis colour palette for NA



## Data Visualisation with ggplot2: Setting Themes, Coordinates and Labels
2 comments
Data Visualisation in RStudio - other possible layers of information such as facets and scales

Introduction

Now that we learned how to choose our data, and apply layers of aesthetics and geometries, let’s explore other possible layers of information such as facets and scales.

Good practice

> library(ggplot2)
> library(RColorBrewer)
> library(viridis)
> library(ggsci)
Setting your working directory in RStudio
Step 1. We recommend you to work in the Project folder Project_Test that we created previously, either by clicking directly on the Project_Test or using the following command

> setwd("/Users/imac/Desktop/exerciseR/Project_Test")
> getwd()
[1] "/Users/imac/Desktop/exerciseR/Project_Test"
Step 2. As a reminder, you can create a specific script file to write your commands and related comments.

Setting Themes
Themes encompass options for the general graphical display of a graph, by modifying non-data output such as the legends, the panel background, etc… You can refer to “theme()” usage including examples in https://www.rdocumentation.org/packages/ggplot2/versions/3.3.3/topics/theme

Step 1. We will use the same new_key variable that we created in the previous Step

> new_key <- key + geom_point(aes(color = Species)) + 
geom_smooth(aes(color = Species, fill = Species), 
method = "lm") + facet_wrap(~Species, scale='free_y') 
> new_key + scale_color_tron() + scale_fill_tron()
graph, using new_key variable from the step before

Step 2. By default the themes are set to “theme_grey()” (or gray). The same previous output will be displayed if you type:

> new_key + scale_color_tron() + scale_fill_tron() + theme_grey()
Step 3. Don’t forget that, as for other functions, you can learn about other variants and options of themes while typing the command in RStudio. If not, simply start typing the function and hit the TAB button on your keyboard. Here is different theme:

> new_key + scale_color_tron() + scale_fill_tron() + theme_minimal() 
graph with iron colours and minimal theme

Step 4. You can now customize the background color, and also set the legends positions. Here is an example on how to place the legend at the bottom.

> new_key + scale_color_tron() + scale_fill_tron() + 
theme_minimal() + theme(legend.position = "bottom", 
panel.background = element_rect(fill = "#e0ebeb"), 
legend.key = element_rect(fill = "#669999"))
graph, with background colour and legend set

Step 5. Note the difference hereafter if when we leave the default theme (just by removing “theme_minimal()” to come back to basic)

> new_key + scale_color_tron() + scale_fill_tron() +  
theme(legend.position = "bottom", panel.background = 
element_rect(fill = "#e0ebeb"), 
legend.key = element_rect(fill = "#669999"))
graph, back to basic theme

Step 6. It is also possible to customize the position or justification of your legend. To do this, you can use “theme_position()” to set the position in the whole panel, and “theme_justification()” to set the position in the legend area. They are defined with a vector of length 2, indicating x and y positions in terms of space coordinates, where c(1, 0) is the bottom-right position

> new_key + scale_color_tron() + scale_fill_tron() + 
theme(legend.position=c(1,0), legend.justification = c(1, 0))
 graph with customized position of legend

Setting Coordinates
The coordinate system controls the position of objects into the main panel, as well as axis and grid lines display. The most classically used are cartesian coordinates, but many other exist

Step 1. Let’s start with “coord_cartesian()” which is set by default. The following commands will have the same output

> new_key + scale_color_tron() + scale_fill_tron()
> new_key + scale_color_tron() + scale_fill_tron() 
+ coord_cartesian()
Step 2. We are now used to see specific options with each function. Let’s fix a ratio from x to y axis, and see how this affects the display in the main panel

> new_key2 <- ggplot(data = iris) + geom_point(aes(x = 
Petal.Length, y = Petal.Width, color = Species)) + 
facet_wrap(~Species) + scale_color_tron() + scale_fill_tron()
 
> new_key2 + coord_fixed(ratio = 2)
Note. The panel borders colored in grey show you the new display

graph illustrating cartesian coordinates

Step 3. A logarithmic transformation can be applied to x and y-axis

> new_key2 + coord_trans(x = "log2", y = "log2")
graph illustrating logarithmic coordinates

Step 4. It also possible to swap x-axis values to y-axis and vice-versa

> new_key2 + coord_flip()
graph illustrating x and y axis swapped

Step 5. Compare these 2 outputs to see how it changes the display for other plot types

Barplot
> new_key3 <- ggplot(iris, aes(x=Petal.Length, Petal.Width)) 
+ geom_bar(stat="identity", fill="white", color="red3")
> new_key3
Right Barplot
> new_key3 + coord_flip()
graph illustrating output compare

Setting Labels
The Labs or Labels system controls the legend, axis and plot labels. You can find examples and usage here https://ggplot2.tidyverse.org/reference/labs.html

Step 1. Let’s start with default “labs()” to set the legend title. The following commands will have the same output

> new_key + scale_color_tron() + scale_fill_tron()
> new_key + scale_color_tron() + scale_fill_tron() 
+ labs(color = "Species")
Step 2. Notice how we modify the legend title by:

> new_key + scale_color_tron() + scale_fill_tron() + 
labs(color = "Iris_Species")
graph illustrating modified legends

Step 3. To modify also the x-axis and y-axis names

> new_key + scale_color_tron() + scale_fill_tron() + 
labs(color = "Iris_Species", x = "Sepal Length values", 
y = "Petal Length values")
graph illustrating modified names

Step 4. Note that it might often be possible to do the same thing using different functions or options. Let’s see 2 options to add a title

> new_key + scale_color_tron() + scale_fill_tron() + 
labs(color = "Iris_Species", x = "Sepal Length values", 
y = "Petal Length values", title = "Petal vs. Sepal Legnth")

 
> new_key + scale_color_tron() + scale_fill_tron() + 
labs(color = "Iris_Species", x = "Sepal Length values", 
y = "Petal Length values") + ggtitle("Petal vs. Sepal Legnth")
graph illustrating added titles


## Working with ggplot2
10 comments
Here are some questions on Working with ggplot2 for you to check your learning so far. Please discuss your answers with other learners in the discussion section below:

Question 1.
What is the difference between the 3 following commands?

> ggplot(data=iris, aes(x=Sepal.Length, y=Petal.Length, 
color = Species)) + geom_point()
> ggplot(iris, aes(x=Sepal.Length, y=Petal.Length, 
color = Species)) + geom_point()
> ggplot(iris, aes(Sepal.Length, Petal.Length, 
color = Species)) + geom_point()
Each will produce a different output
They will produce the same output
Only the first one is correct
Only the second one is correct
Only the third one is correct
I don’t know
Question 2.
In your opinion, what is the correct command that should have generated the following plot (knowing that this type of plot is called a boxplot)

Option 1
> ggplot(data = iris, aes(x = Sepal.Length, y = Petal.Length))+ 
geom_boxplot()
Option 2
> ggplot(data = iris, aes(x = Sepal.Length, y = Petal.Length))+ 
geom_boxplot() + facet_wrap(~Species)
Option 3
> ggplot(data = iris, aes(x = Sepal.Length, y = Petal.Length, 
color = Species)) + geom_boxplot()
boxplot generated by the commands in question

Question 3.
Consider the following command that generated the associated plot below:

> ggplot(data = iris, aes(x = Sepal.Length, y = Petal.Length, 
color = Species)) + geom_boxplot() + geom_smooth(se=FALSE) + 
facet_wrap(~Species, scale='free_y') + 
scale_color_brewer(palette="RdYlBu")
plot generated  in question 3

Which layer(s) should we add in order to produce the following plot?

+ facet_wrap(~Species)
+ facet_wrap(~Species, scale='free_y')
+ facet_wrap(~Species, scale='free_y') + 
labs(x = "Sepal Length", title = "Petal vs Sepal Length")
+ facet_wrap(~Species) + labs(x = "Sepal Length", 
title = "Petal vs Sepal Length")
+ facet_wrap(~Species, scale='free_y') + labs(y = "Petal Length", 
title = "Petal vs Sepal Length")
plot generated in question 3 with layers

Please try to answer the questions yourself first and then compare the results with other learners. Finally, you can find solutions in the download area.


## Data Visualisation in R/RStudio - Conclusion and Resources
2 comments
Data visualisation in R and RStudio makes it possible to easily use basic plotting functions, or apply more advanced functions through packages.

As you must have noticed throughout this week, the undeniable added value of R/RStudio compared to the more classical resources such as Excel, is the ability to produce publication-ready graphics. For this, you can either use default functions and options, which already produce a highly controlled output quality, or you can also pre-define advanced options and use them in variables.

Moreover, you can easily ensure reproducibility of your output, simply by writing your commands in a script. This would allow you to apply changes in your data or display with the same final output.

Importantly, there are plenty of websites and online resources to easily take you through using packages, either using their official documentation, or personal contributions and forums to tackle certain aspects of the analysis.

Although using packages such as ggplot2 might seem cumbersome to understand in the beginning, I’m pretty confident that you will quickly understand that the game is worth the candle!

Resources for R/RStudio
There are obviously tons of online resources to guide you through introductory or advanced R or RStudio usage.

One I can suggest is http://www.sthda.com/english/wiki/ggplot2-essentials. This very famous website (STHDA: Statistical Tools for High-throughput Data Analysis) will give you a step by step explanation of how to use R and RStudio, taking you from very basic introductions of R principles to more advanced packages for data manipulation, statistical analysis, or data visualization with many detailed examples.

Resources for ggplot2
As I believe you understood at this point, there are still many functions and options to explore in ggplot2.

For a list of functions available in ggplot2, there are again many resources on the web. It is important to browse these resources by yourself, because each one of you will have specific display requirements and we obviously cannot cover them all.

However, now that I showed you some basic usage, and in order to help you start using ggplot2 with good references, my personal preferences are for the following resources:

https://ggplot2-book.org/ This is the on-line version of the book “ggplot2: elegant graphics for data analysis” published (Springer publishing). This a very detailed resource on explanations of each layer with helpful examples

https://ggplot2.tidyverse.org/reference/
This website presents ggplot2 functions classified by type of layer (aesthetics, geoms,…) and in which you can click on each function name to be automatically redirected to its usage

https://www.rdocumentation.org/packages/ggplot2/versions/3.3.3 The first page of this documentation page presents a “cheatsheet” that is a printable document summarizing functions in ggplot2. At the bottom of the page you will also find a list of clickable ggplot2 functions to access their usage including examples. Many links to start working with ggplot2 are also given here

## Summary of Week 3
2 comments
Seeing how to handle and analyse your data using the command-line in Week 1, and how to write and run bash scripts to automate that handling during Week 2, in Week 3 you saw how you can analyse and visually represent your data using R and RStudio.

We first reviewed together the key concepts and basics of using the R language. You had an opportunity to learn how to manipulate lists, vectors or dataframes, and how to read and interrogate data from files to create simple plots under R. We then saw how RStudio can simplify these processes, and you learned how to work under RStudio to easily create simple or more complex plots using available packages specifically designed for that purpose.

Please leave your comments in the comments section below - in particular, did you learn anything that will be applicable to your own practice or study?

We hope you enjoyed this course!
Please attempt the End of the Course Test in the next step. Getting 70% or more on the test and reviewing (marking as complete) at least 90% of steps in this course, make you eligible for the Certificate of Achievement which you can download at the end of the course.

Good luck with your next step!

