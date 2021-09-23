# Introduction to Linux for biologists

## Capturing user input


## Who am I?
We are going to create a script called whoami.sh
In that script we want to perform the following actions:
1. Prompt the user to provide their name
2. Store the user input as variable
3. Return that variable in a greeting back to thye user

* First,include the shebang at the top of our script. This tells us that we want to invoke the Bash Interpreter
* Next, we send a request for the user to enter their name in the terminal using echo.
* In our script, we will then use the read command to capture the user input.
  For this, we use the command read, followed by the variable name where we want to store the user’s response.
  It’s always worth noting that, when we declare a variable– in this case, name– for the first time, we don’t need to use the dollar symbol as a prefix
* Finally, we’ll return the user input back to the user on the command line using echo. 
  The user’s response was stored in a variable called name, so I’ll need to use this as part of the command to return that input. Below, there are some good practises to follow for reading variables. First, always prefix the variable name with the dollar symbol, enclose the variable name in curly braces to avoid unwanted consequences, and use double quotes for the same reason.

whoami.sh
```
#!/usr/bin/env bash
echo "Please enter your name."
read name
echo "Hello ${name}, it’s nice to meet you!"

```


What happens if you want to ask the user for multiple inputs? In this case, we need to split that input and store it as multiple variables
We’re going to ask the user to enter two of their favourite foods. When the user enters their two foods, which will be on a single line separated by a space, these gets stored as two variables, food1 and food2. And then within the script, we ask it to return both of those values back to the user using echo.

Lets create a script where we provide our two fav foos as apple and pears separeted by a comma

We’ll then change the permissions on our file to make it executable, and then we’ll run our script, which will prompt us to enter our two favourite foods. You’ll see here that we’re going to enter two food separated by a space. In this case, apples and pears. We’re then going to return this back to the user. So here, you’ll see that your favourite foods are apples and pears. Remember that we split our input by that space, so apples were stored as our first variable and pears were stored as the second. So when we return it as part of that echoed string, it’s then split by an and.

In this example, the user gave the same number of responses as we had variables to store them. But what would happen if the user had given three foods instead of just two? Say, for example, apples, pears, and oranges. Now, when we split their response by the space character, what would happen is apples would be stored to the variable food1, but both pears and oranges with the space between them would be stored as a string to food2. This means that, when we provide more options or responses than we have designated variables, any extra text will remain as part of that final variable.
As we’ve covered, the read command interactively captures user input, and so far, we’ve been using it in its simplest form. There are several options which are available for the read command to extend its functionality. To see them, we can type help read, which will give us the manual for the read command. 
There are two options I would like to introduce– minus p, which allows a message to be displayed before the input is read, and minus s, which will hide the user input as they’re typing it. And this is useful for secrets like passwords. To demonstrate this, let’s create a new script called login.sh, in which we’ll ask the user to provide their username and password.



food.sh

```

#!/usr/bin/env bash
echo "what are your two favourite foods?"
read food1 food2
echo "Your favourite foods are: ${food1} and ${food2}"

```


While the user is entering their password, we will be hiding the text they’re typing but storing their response as a variable in our script, which we can then either return back to the user or use later on in our script if we wanted to. Let’s take a closer look at the minus p option first, which is used in both of our read commands. When we were prompting for user input previously, the prompt and the user’s input were on two separate lines. When we use the minus p option, this allows the user input to be entered on the same line as they see the message or prompt. 
We can also look at the minus s option.
Here, we’re going to use this when we collect the password. The minus s option hides the user input as it’s being entered, so when we prompt asking them for their password, as the user is typing, the text that they are inputting won’t be shown in the terminal. But the text that they’ve entered isn’t lost. It will be stored in a variable in the background as part of our script. 
And we will then return it back to the userlogin.sh.
We’ll then update our permissions to make that script executable, and then we’ll execute our script, login.sh. First, you’ll see the prompt for your username. And you’ll be entering your response on the same line as that message. Next, you’ll see it prompt for your password. It looks like you're doing nothing here, but while typing input is hidden from the terminal.
This was stored as a variable, and then we’ve returned that back to you as an example in that message

login.sh
```
#!/usr/bin/env bash
read -p "Enter your username: " username
read -sp "Tell me your password: " password
echo -e "\nHi ${username}, your password is ${password}"

```
The ‘-e‘ option in Linux acts as interpretation of escaped characters that are backslashed.




* To get the use of echo options

   'man echo' in the terminal
   https://www.gnu.org/software/bash/manual/html_node/Escape-Character.html
   
   
 ### Passing Commandline Arguments to Bash Scripts
 
A command line argument is a parameter that we can supply to our Bash script at execution. They allow a user to dynamically affect the actions your script will    perform or the output it will generate.

To pass an argument to your Bash script, your just need to write it after the name of your script:

./fruit.sh my_argument
In our Bash script, there are several reserved/pre-defined variables which we can use to recall the user-defined parameters. The first argument is stored in $1, the second in $2, the third in $3…and so on. We cannot use $0 as that references your Bash script itself.

Let’s see how this works using an example script:

#!/usr/bin/env bash
 
echo "The first fruit is: $1"
echo "The second fruit is: $2"
echo "The third fruit is: $3"
If we run our script and don’t give an argument, we will see no output for our pre-defined variables:

./fruit.sh
The first fruit is:
The second fruit is:
The third fruit is:
Alternatively, if we provide three fruits, our script detects these and will return those values back to use via the pre-defined variables.

./fruit.sh apple pear orange
The first fruit is: apple
The second fruit is: pear
The third fruit is: orange
Sometimes, we may want to access all of the arguments. We can do this using $@.

Let’s update our example script to return all of the fruits provided as arguments to the script as well:

#!/usr/bin/env bash
 
echo "The first fruit is: $1"
echo "The second fruit is: $2"
echo "The third fruit is: $3"
echo "All fruits are: $@"
When we run our script, you can see that we now have an extra output which lists all of the fruits we gave to our script on the command line.

./fruit.sh apple pear orange
The first fruit is: apple
The second fruit is: pear
The third fruit is: orange
All fruits are: apple pear orange

## Conditional Expressions

### If and If Else Statements

we looked at how to assign values to variables. In this section, we’ll be looking at how to compare them to other variables or values and perform simple checks, like whether a file exists, using both conditional expressions and conditional statements


IF statements
Conditional statements come in many forms. The most basic form essentially says: IF our conditions are met, THEN execute the following code.

We can write our if statements in several ways:

if [[ condition ]]
then
    	command
fi
So, if the expression within the square brackets returns true, then the code found between then and fi will be executed. If the expression returns false, the script will ignore (i.e. not execute) any code which lies between then and fi. Notice that we end our if statement with fi which is if spelt backwards.

An alternative format that you might come across uses a ‘;’ to allow then to be on the same line as the conditional expression:

if [[ condition ]] ; then
    	command
fi
You’ll notice that in these examples, we’ve used spacing to indicate the code which will run if the conditional expression returned true. This is known as indenting and, although there are no requirements for it in Bash, it is a good coding practice to follow for clean, readable code.

A simple example you can try on the command line yourself is:

if [[ 1 == 1 ]] ; then
    	echo hi
fi
Here we are saying, if 1 is equal to 1 (1 == 1) then return ‘hi’ back to us. As this conditional expression is always true (as 1 is 1) then, it will always return ‘hi’.

* Important to note that white space should be added before and after the == *

### IF statements with AND/OR logic
We can use two (or more) conditional expressions with our if statement using the AND and/or OR conditions.

For example, let’s say that we have a file and we want to check that it exists, is readable and that it isn’t empty. If our file meets these criteria, we want to print “File is good” and if not, print “File is bad”.

First, let’s write our script:

#!/usr/bin/env bash
 
# Set the path for our file
file="file.txt"
 
# Check whether file exists, is readable and has data
if [[ -e ${file} ]] && [[ -r ${file} ]] && [[ -s ${file} ]]
then
    	# Execute this code if file meets those conditions
    	echo "File is good"
else
    	# Execute this code if file does not meet those conditions
    	echo "File is bad"
fi
Now, let’s run our script knowing that our file can’t meet the criteria as it doesn’t exist:

./script.sh
File is bad
Next, let’s create an empty file and try running our script again:

touch file.txt
./script.sh
Again, our script returns “File is bad” as it hasn’t met all of our conditions. Finally, let’s add some data to our file and try again:

echo "hi" > file.txt
./script.sh
Bingo! We have met all of the conditions and now our “File is good”.

IF..ELSE statements
We can extend our conditional statement to have another clause by using an if..else statement. Here we are saying, IF our conditions are met, THEN execute the following commands. However, ELSE IF these conditions are not met, execute a different set of commands.

The syntax for this looks like:

if [[ condition ]]
then
    	command1
else
    	command2
fi
IF..ELIF..ELSE statements
Sometimes, you may need to test more than one statement. The syntax for this looks like:

if [[ condition1 ]]
then
    	command1
elif [[ condition2 ]]
then
    	command2
else
    	command3
fi
It’s worth mentioning that you can have more than one elif clause in your if..elif..else statement. But, as we will soon discuss, there are more efficient ways of building that type of conditional statement.

## Conditional Expressions

Modern Bash syntax for conditional expressions encases our comparative expression inside double square brackets ([[ and ]]). In the examples below, we’ll show the syntax of different types of conditional expression.

It’s worth noting that you won’t be able to run these examples on the command line directly as we’re just showing you the syntax. Don’t worry, you’ll be able to build on these in the task at the end of the section!

The syntax for this is:

[[ option arg1 ]]
or

[[ arg1 operator arg2 ]]
A conditional expression returns a Boolean value i.e. true or false. If the condition is met, it will return true and if not, false.

It’s worth noting that the spacing is important. Here are some examples of valid and invalid conditional expression syntax.

Valid:

[[ -f ${file} ]]
Invalid:

[[ -e file]]

[[-e file]]

[[-efile]]
File and variable operators
When we process files in our Bash scripts, it is often useful to check that they exist or whether they’re empty before the rest of our script proceeds. File operators allow us to perform checks on files and give us the opportunity to handle errors gracefully.

Below are some of the most commonly used file operators.

Returns true if the file exists:

[[ -e ${file} ]]
Returns true if the file exists and is a directory:

[[ -d ${directory} ]]
Returns true if the file exists and is a regular file:

[[ -f ${file} ]]
Returns true if the file exists and is readable:

[[ -r ${file} ]]
Returns true if the file exists and has a file size > 0:

[[ -s ${file} ]]
We can also use conditional expressions to perform sanity checks on our variables.

For example, checking whether a value has been assigned to a particular variable (e.g. var):

[[ -v ${var} ]]
Or, to check that the variable length is greater than 0:

[[ -n ${string} ]]
Or, that the length of the variable is 0:

[[ -z ${string} ]]
String comparisons
Strings, as sequences of characters, can be compared. There are two string conditional expressions you need to be aware of:

Is equal to ==
Is not equal to !=
This condition will return true if string1 and string2 are identical:

[[ ${string1} == ${string2} ]]
This condition will return true if string1 and string2 are different from one another:

[[ ${string1} != ${string2} ]]
Arithmetic comparisons
In Bash, we don’t use the same syntax for string and arithmetic comparisons. There are six main arithmetic expressions:

Is equal to -eq
Is not equal to -ne
Is less than -lt
Is less than or equal to -le
Is greater than -gt
Is greater than or equal to -ge
This condition will return true if arg1 is equal to arg2:

[[ ${arg1} -eq ${arg2} ]]
Alternatively, this condition will return true if arg1 is not equal to arg2:

[[ ${arg1} -ne ${arg2} ]]
This condition will return true if arg1 is less than arg2:

[[ ${arg1} -lt ${arg2} ]]
Meanwhile, this condition will return true if arg1 less than or equal to arg2:

[[ ${arg1} -le ${arg2} ]]
This condition will return true if arg1 is greater than arg2:

[[ ${arg1} -gt ${arg2} ]]
And, finally, this condition will return true if arg1 is greater than or equal to arg2:

[[ ${arg1} -ge ${arg2} ]]
Performing multiple comparisons
Using Bash syntax, we can also combine comparisons using the && and the || operators which represent AND and OR respectively.

The following expression would only return true in the event that var1 is equal to var2 AND var3 is equal to var4:

[[ ${var1} == ${var2} ]] && [[ ${var3} == ${var4} ]]
Alternatively, the following expression would only return true when either var1 is equal to var2 OR var3 is equal to var4:

[[ ${var1} == ${var2} ]] || [[ ${var3} == ${var4} ]]
Your task:
Using the information in this and the previous section, write a Bash script called temperature.sh that:

reads in a command line argument into a variable called temperature
has a variable called min_temperature and give it a variable of 10
has a variable called max_temperature and give it a variable of 30
returns “Too hot!” if temperature is greater than max_temperature
returns “Too cold!” if temperature is less than min_temperature
returns “Just right!” if temperature is greater than the min_temperature and less than max_temperature
Test your Bash script with the following command:

./temperature.sh 22
What does your script output?

#!/usr/bin/env bash

temperature="$1"

min_temperature="10"

max_temperature="30"

if [[ ${temperature} -gt ${max_temperature} ]]

then

echo "Too hot!"

elif [[ ${temperature} -lt ${min_temperature} ]]

then

echo "Too cold!"

elif [[ ${temperature} -gt ${min_temperature} ]] && [[ ${temperature} -lt ${max_temperature} ]]

then

echo "Just right"

fi

./temperature.sh 22

Output : Just right

## Switch Case Statements
* Sometimes of conditional logic may be too complex for if statements. In these situations, the issue isn’t that it’s impossible to write the logic as nested if statements, but that to do so could result in confusing code.*

In these situations, we can use case statements to check each of our conditions in turn and process commands based on those conditions.

The case syntax looks like this:

case $string in
    	pattern_1)
      	command
      	;;
    	pattern_2)
      	alternate command
      	;;
    	*)
      	default command
      	;;
esac
Let’s break this down. First, we start with case followed by the variable or expression we want to test and then in. Next, we have our case patterns against which we want to check our variable or expression. We use the ) symbol to signify the end of each pattern. After each pattern, you can then specify one or more commands you want to execute in the event that the pattern matches the expression or variable, terminating each clause with ;;. As our last switch, it is common practice to have a default condition which is defined by having * as the pattern. Finally, we signify the end of our case statement by closing it with esac (case typed backwards!).

Here is a simple example:

#!/usr/bin/env bash
 
fruit="pineapple"
case $fruit in
    	apple)
      	echo "Your apple will cost 35p"
      	;;
    	pear)
      	echo "Your pear will cost 41p"
      	;;
    	peach)
      	echo "Your peach will cost 50p"
      	;;
    	pineapple)
      	echo "Your pineapple will cost 75p"
      	;;
    	*)
      	echo "Unknown fruit"
      	;;
esac
First, we set our variable, fruit, to have the value “pineapple”. We then compare this against several conditions looking to see whether the value of our variable matches the pattern provided. In the event that none of the patterns match our fruit, we have a default response “Unknown fruit”.

As one of the patterns is indeed “pineapple” we meet that condition and return:

Your pineapple will cost 75p
Your task
Create a Bash script called farm.sh that uses a case statement to perform the following functions:

Stores a command line argument in a variable called animal
Use a case switch statement which has the following conditions and responses
When the user enters cow, return “Here, moo”
When the user enters sheep, return “There a baa”
When the user enters duck, return “Everywhere a quack”
Otherwise, return “Old MacDonald had a farm”


```
nano farm.sh

#!/usr/bin/env bash
animal="duck"
case $animal in
    	cow)
      	echo "Here, moo"
      	;;
    	sheep)
      	echo “There a baa”
      	;;
    	duck)
      	echo "Everywhere a quack"
      	;;
    	*)
      	echo "Old MacDonald had a farm"
      	;;
esac
```

## For Loops

What is a loop? A loop is a construct which allows you to repeatedly execute the same commands. We will be discussing three types of loops: for loops, while loops and until loops.

Let’s start by looking at for loops. The basic syntax for a for loop is:

for variable in ${list}
do
    	# Execute some commands
done
We can define a list like so:

my_list="item1 item2 item3"
As a simple example, let’s create a list of fruits and use a for loop to return each item from our list:

fruits="apples pears oranges"
for fruit in ${fruits}
do
    	echo ${fruit}
done
This will output:

apples
pears
oranges
We can also use a for loop to iterate over a series of numbers. In this example, we’ll process the numbers 1 – 3 using a sequence expression. Here, you’ll see the range is specified by a beginning number (1) and an ending number (3) separated by ‘..’. This indicates that we want the sequence of numbers from the beginning to the ending number inclusive i.e. 1, 2 and 3.

for n in {1..3}
do
    	echo ${n}
done
This will output:

1
2
3
A common use of for loops is to iterate over the contents of a directory. Here is an example of how to list all files in the current directory:

for file in *
do
    	echo ${file}
done
Here we use ‘*’ as a wildcard to ask for all files and directories. We could extend this to look text files:

for file in *.txt
do
    	echo ${file}
done
This would return only those files that have a .txt file extension.

Finally, we’ll introduce you to the for loop syntax that uses three expressions: an initial value, a terminal value and an increment/decrement. Notice here that the increment uses the ‘++’ notation which simply means add 1.

for (( i=1; i<=3; i++ ))
do
    	echo $i
done
This example returns the same output as our earlier number series.

1
2
3
Your task
Create a for loop which iterates from 1 to 5 in increments of 1. If the value is 2 return “fizz” otherwise, return “buzz”.
```
for (( i=1; i<=5; i++ ))

do

if [[ ${i} -eq "2" ]]

then

echo "fizz"

else

echo "buzz"

fi

done
```
## While Loop and Until Loop
Both for loops and while loops are very similar. Typically, we use for loops where we know exactly how many iterations we need – i.e. they have a definitive start and end point. On the other hand, while loops are used where we don’t know the limitations on tasks such as read in a file or asking a user for input. They just keep iterating as long as the specified condition has been met.

The basic syntax for a while loop looks like this:

while [condition]
do
    	# Commands to run
done
First, let’s look at how not to do a while loop:

i=1
while [[ $i -eq 1 ]]
do
    	echo "hi"
done
This is what’s known as an infinite loop because the condition will always return true – i.e. nothing is changing. In this example, “hi” will just keep being printed to the terminal until we force it to stop using Ctrl+C on our keyboard.

So, that was how to use while loops in the wrong way. But, what do they look like when they are being used properly:

i=1
while [[ $i -le 3 ]]
do
   echo "$i"
   (( i++ ))
done
What we’re doing here is setting our variable to have an initial value of 1. When the loop begins, our variable is 1 (i.e. less than 3) and so the condition returns true. That means that we’re going to execute the code body, returning our variable value, 1 to the terminal. Next, we increment our variable value from 1 to 2 using the ++ notation. This continues while our variable has a value less than or equal to 3.

The result:

1
2
3
Another common use for while loops is reading in the contents of a file. Here is an example:

while read data
do
   echo "${data}"
done < infile.txt
This is what is known as a while loop. What do we mean by this? In this example, the while loop will only keep iterating while there are lines to be read from the given input file. Here, infile.txt is the name of the file that we are going to be looping over. The read command will process the file, line by line, into the data variable. Once it reaches the end of the file, the while loop will be terminated.

Until loop
We just looked at an example of a while loop. Now, we’re going to look at run-until loops. The main difference is that while loops are designed to run while a condition is satisfied and then terminate once that condition returns false. On the other hand, until loops are designed to run while the condition returns false and only terminate when the condition returns true.

The structure of until loops is almost identical to that of a while loop:

until [condition]
do
    	# Commands to run
done
For example, this loop would run until the variable is greater than 3:

i=1
until [[ $i -gt 3 ]]
do
    	echo $i
    	((i++))
done
This would output:

1
2
3


## Bash Functions

When you’re writing Bash scripts, you’ll often find that there are repetitive tasks. Instead of copying and pasting the same code to multiple places in your scripts, try using a function.

Functions are a great way of producing reusable code! They are essentially a set of commands that can be called as many times as you need in your script. What’s even better is that functions are not unique to Bash, they’re a core component of many other programming languages too.

Bash function syntax is pretty straightforward. We start off by defining the function name, followed by parentheses. The commands that we want to execute are found between the curly brackets and are known as the body of the function.

function my_function() {
    	#some code
}
There is an alternative syntax where you don’t have to prefix that first line with function:

my_function() {
    	#some code
}
However, it is much easier to pick out our functions if we use the previous syntax. It’s also a good idea to make sure that the names of your functions are relative and descriptive so that you can quickly see what they’re going to do.

When we define a function, we are not executing it. Let’s use a simple toy example to demonstrate where we are using a function to return “Hello world” back to the terminal. We’ll call our function say_hello. You can see that we don’t execute the code in the function body until we specifically call (or execute) the function with say_hello.

#!/usr/bin/env bash
 
# Define a function to print "hello world"
function say_hello() {
    	echo "Hello world"
}
 
# Execute the say_hello function
say_hello
This would output:

Hello world
We can adapt out function to take arguments using reserved variables. To access the first argument given to the function, we use the variable $1. Let’s tweak our script to use an argument, our name, that is provided to our say_hello function.

#!/usr/bin/env bash
 
# Define a function to print "hello world"
function say_hello() {
    	echo "Hello $1"
}
 
# Execute the say_hello function
say_hello "Victoria"
This would output:

Hello Victoria
Functions are one of the best ways to produce scalable and readable code. One general rule of thumb is not to make your functions too big. You can call a function within a function, so, break each function down into small, clear tasks.

Your task
Create a function called file_exists taking the first argument (a filename) which it uses to see if the file exists. If it doesn’t, return “File does not exist: “, followed by the filename.

Note: you can use the “!” notation when you want to check a negative.

If file exists:

if [[ -e $1 ]]

If file does not exist:

if [[ ! -e $1 ]]

Please try to answer the questions yourself first and then compare the results with other learners. Once you’ve tried the exercise, you can find solutions in the download area.

```
file="no_file.txt"
function file_exists() {
 if [[ ! -e $1 ]]
 then
 echo "File does not exist: $1"
 fi
}
file_exists "${file}"
#Notice that we are checking to see that a file doesn't exist, not that it does.
if [[ ! -e $1 ]]
```
#!/usr/bin/env bash

# Define a function to check file

function filecheck() {

if [[ -e $1 ]]

then

echo "$1"

elif [[ ! -e $1 ]]

then

echo "File does not exist: $1"

fi

}

# Execute the filecheck function

filecheck no-file.

## Track the Progress of Your Script and Redirect Script Outputs and Errors

Tracking the progress of your script
Now, let us imagine you have a long and complex Bash script. You execute your script, it’s started running and you’ve gone off to make a cup of tea. Ten minutes later, you come back to check on its progress but, how do you know what’s going on and where you’ve gotten up to in your script?

There are many different ways in which we can track the progress of our scripts. The simplest is to break your script down into sections and output a progress statement when you start and/or finish each section.

For example, let’s set our name as a variable and count the number of characters it contains.

#!/usr/bin/env bash
 
# Set your name as a variable
name="Victoria"
 
echo "Counting number of characters in name"
printf -- "${name}" | wc -m
As expected, we have our progress statement and the number of characters in our name:

Counting number of characters in name
            8
Now, while this may seem excessive given the simple example, it’s clear that once we start to build up our scripts, adding progress statements will be invaluable. Particularly when were discussing loops this week, where it’s possible for your scripts to get stuck in an infinite loop, failing to exit. In those situations, progress statements are absolutely essential for debugging!

Discuss with your fellow learners:
Can you see a situation where you would need to track the progress of your script?

Redirecting script outputs and errors
Despite your hardest efforts, sometimes your Bash scripts will do unexpected things. This is when we need to debug. If you have a long Bash script, it can be tricky to work out where things went wrong.

To help with debugging, we can output progress statements at key points in our code e.g. “Reading in file: x”. However, these can easily fill up your terminal and become difficult to follow. A simple solution is to write these progress statements to one or more log files.

Redirecting the output of scripts and commands to files
Simply put, redirection is the mechanism by which we can send the output of a command or script to another place. When we want to capture the output from a command or script, we usually choose to redirect those outputs into a file.

To redirect the outputs of a script, we can use the > symbol:

script.sh > output.txt
Redirection using the > symbol works not only for scripts, but any Bash command:

echo "hello world" > hello.txt
cat hello.txt
hello world
Linux streams and file descriptors
Before we take an in depth look at how we redirect our outputs and errors to log files, we first need a crash course in Linux streams and file descriptors. These streams are handled like files – i.e. you can read from them and you can write to them.

There are three streams you should be aware of:

stdin (standard input)
stdout (standard output)
stderr (standard error)
This sounds much more complicated than it really is. In a nutshell, stdout refers to the output from a command and stderr refers to the errors a command generates. The final stream, stdin refers to command line inputs which we’ll cover later in the week.

Next, we need to understand file descriptors. A file descriptor is just a (positive) integer that represents an open file. Each of our Linux streams (i.e. stdin, stdout and stderr) has been allocated a unique number in order to identify them.

All you need to remember is which of the ids below corresponds to each of the streams:

0 => stdin

1 => stdout

2 => stderr

I/O redirection
To start understanding how these streams work, let’s look at redirecting the output from a script into a single file.

Example script:

#!/usr/bin/env bash
 
# A script that tries to change directory
 
echo "Changing to a directory that doesn't exist"
cd foo
As you can see, our script returns the printed progress statement and an error that tells us that the directory we’re trying to migrate to doesn’t exist on our filesystem.

./script.sh
Changing to a directory that doesn't exist
script.sh: line 6: cd: foo: No such file or directory
These two messages are being delivered to the terminal by two different Linux streams. The first message, our progress statement, is delivered via stdout. Meanwhile, the error message is delivered via stderr.

Now, let’s see what happens when we try to redirect the outputs from that script into a file called output.txt:

./script.sh > output.txt
./script.sh: line 6: cd: foo: No such file or directory
OK, so, we can see that the stdout has been redirected to our output file but, the error is still being displayed.

cat output.txt
Changing to a directory that doesn't exist
Why is this? Well, when we use > to redirect to a file, by default, the system will only redirect the stdout.

But, what about our errors being delivered via stderr, how can we capture those?

To simplify things, let’s first look at how to redirect stdout and stderr to two different files. We’ll use the > symbol with our file descriptors (1 for stdout and 2 for stderr) to redirect our outputs to output.txt and our errors to error.txt respectively.

./script.sh 1>output.txt 2>error.txt
This command returns nothing back to our terminal. Using the cat command, we can see that, as expected, our outputs and errors have been written to output.txt and error.txt respectively.

Our stdout (progress statement returned using echo): cat output.txt Changing to a directory that doesn’t exist

And our stderr (errors):

cat error.txt
./script.sh: line 6: cd: foo: No such file or directory
In order to redirect the stdout and the stderr to the same place, we need to use a new term: 2>&1. When we use this, we redirect using the same syntax as before, but add 2>&1 to the end of our command.

This is how it works in practice:

./script.sh > combined_output.txt 2>&1
Now, if we look at our combined output file, we can see that we’ve captured both the stdout and the stderr.

cat combined_output.txt
Changing to a directory that doesn't exist
./script.sh: line 6: cd: foo: No such file or directory

Discuss with your fellow learners:

Can you see a situation where you would need to track the progress of your script?

Ans:

For example, in transcriptomics, you want to output information for a gene of interest. You could want to find the gene name "PBANKA_etc" (this is Plasmodium berghei gene annotation), and then find if it's up or down regulated, by how much, etc. I would want to go until the end of the file (so an until loop) but I would want to have something in there to say how many times that gene is represented and then give me some output for each time it appears. The gene could be represented once or it could be represented hundreds of times in multiple RNA-seq datasets. If you have many samples, you would have to run the script for each file outputted by the sequencer. I'm probably not explaining it well, but processing NGS data can take ages and I would want a progress report!

## Writing Robust Bash Scripts
Sometimes, despite having the very best intentions, subtle issues can creep into your script causing it to fail with unintended consequences. Fortunately, there are commands available to help with minimising these issues. One of these is the set command.

Let’s take a look at how the set command can help us write robust and secure Bash scripts.

First, how does the set command work? Using the set command allows us to customise the environment in which our scripts are run.

The general syntax for the set command is:

set [options]
There are more than a dozen options available for the set command. To view them, you can run the following command:

help set
In this article, we’ll be focusing on the most commonly used options.

Using set -e to catch errors
Sometimes, commands within your script may fail but, the downstream commands will continue to run. This can be extremely frustrating if you don’t see the error and assume that, as the script completed, everything has worked as expected.

Here’s an example. First, we will try to change into a directory called foo and then list the contents of that foo directory. The key here is that the foo directory doesn’t actually exist so, we can’t get its contents.

#!/usr/bin/env bash
 
cd foo
ls
What happens when we run our script?

script.sh: line 3: cd: foo: No such file or directory
File1 File2
Notice that our script generated an error when the system couldn’t find our foo directory. But, because there wasn’t an exit code, the remaining commands in the script also ran. Unfortunately, this listed the contents of our current working directory and not the foo directory as intended. Imagine if this was part of a long series of output commands and we missed the error….we may accidentally assume that our script ran correctly!

Fortunately, the set -e command comes to our rescue by ensuring that the script will fail whenever an error occurs, no matter the exit code. Try adding set -e to the top of your script:

#!/usr/bin/env bash
 
set -e
 
cd foo
ls
Bingo! This time, we can see that the script terminates as soon as it reaches the first error.

script.sh: line 5: cd:foo: No such file or directory
Using set -u to catch variables that don’t exist
By default, when executing a script, Bash will just ignore variables which don’t exist. In most cases, you won’t want this behaviour as it can have unexpected consequences!

In this example, we will first try to output a variable, $foo, which doesn’t exist and then try to output a simple string, bar.

#!/usr/bin/env bash
 
echo $foo
echo bar
When we run this script, we get the following output:

bar
Notice that the system outputs a blank line for echo $foo. This is because Bash is ignoring $foo as it doesn’t exist.

If we want the script to exit with an error instead of continuing on silently, we can add the set -u command at the top of our script.

#!/usr/bin/env bash
 
set -u
 
echo $foo
echo bar
This will result in our script exiting with the following error:

script.sh: line 6:foo: unbound variable
Notice, our script terminates before running the second echo command.

Displaying executed commands while script is running with set -x
Another default Bash behaviour is to only display results once a script has finished. This can be especially frustrating when you need to debug scripts that take a long time to run.

Let’s take an example script that outputs two simple strings, foo and bar.

#!/usr/bin/env bash
 
echo foo
echo bar
The output from this script would be:

foo
bar
Now, what if we want to know which command is producing each of the results? To find this out, we can use the set -x command which outputs the executed command before printing the command result.

#!/usr/bin/env bash
 
set -x
 
echo foo
echo bar
Running this script would give the following output:

+ echo foo
foo
+ echo bar
bar
As you can see, before executing each of the echo commands, the script first prints the command to the terminal, using a + to indicate that the output is a command. This can be especially handy when you want to debug long scripts.

Combining set options in a single command
Most of the time, you will want to use all of these options together. Instead of writing the commands out, one command per line, we can combine the options into a single command:

set -eux
Using the set command is essential to building robust Bash scripts. Not only is it part of good scripting practices but, will also save you a lot of time and frustration!


## Good Practices

For the most part, it’s easy to write Bash scripts. What’s really tricky is writing good Bash scripts.

What do we mean by this? In short, we mean writing “clean code”. Clean code is:

Easy for someone to pick up and understand
Reusable
Scalable
Expects the unexpected
There will always be situations where you write a script for a one-time only task. In these situations, people tend to cut corners and become more flexible with making a script scalable or reusable. It’s tempting, it is, we’ve all been there. But…inevitably, you’ll almost always need to do the same or similar task unexpectedly in the future. It takes time to write good code in all situations, but, I promise, it’s always worth it in the long run!

Here we’ve put together a simple list of 12 good practices which you can try to follow. This is by no means an exhaustive list and I encourage you to read more widely and continually assess and develop your scripting practices.

Plan ahead
Most of your scripting headaches will be solved by planning ahead. Think about what you want your script to do and expect the unexpected. This means not only thinking about the “happy path” where everything proceeds as you expect it to, but also the exceptions. For example, what should it do if the file its processing is empty or doesn’t exist?

2. Build your script in small steps

Depending on the complexity of the task you’re trying to accomplish, your script may be very small or somewhat larger. Now, even for the most experienced of script writers, there can be a typo or other error in their first attempt. To avoid this, no matter the scale of your script, build it up in small stages and test as you go.

3. Scale up slowly

We don’t just need to consider building the process up step by step, but also the size of the data the script is handling. A simple rule of thumb is to get your script working for a single task first. Then, build up slowly until you reach the full scale of your dataset i.e. handle one file, then 10, then 100… This allows you to predict the resources you need to process the full dataset and get a rough idea of how long this process will take.

4. Comment, comment, comment

OK, I could start by saying that you can never have too many comments in a script. But, that’s not true, at some point you won’t be able to see the wood for the trees and your script will become unnecessarily bloated. There is no hard and fast rule on how many comments you should have in your scripts, it just comes down to common sense. Simply, make sure that you have enough information so that if someone, usually you, picks your script up in 6 months that they know what it does and have a rough idea of how it does it. This isn’t just true for shell scripting, but for almost every other type of programmatic scripting you may encounter.

5. Don’t prolong the life of the script unnecessarily

This is called script longevity and links back to expecting the unexpected. Earlier in the week we looked at using the set command to handle failures and errors. A script should never fall over quietly. Trigger an exit signal when the unexpected happens. This means that a script should exit on the first error and not blithely continue running unnecessarily or fall over without you realising it.

6. Keep on top of variable management

There are a lot of things to consider for variable management. First up is syntax. I like to use upper case for environment variables and lower case for local variables.

MY_ENVIRONMENT_VARIABLE=1
my_local_variable=2
Depending on your preference you can use underscores or camel case. I prefer underscores as it keeps things consistent between variable types. But that’s up to you.

my_underscore_example=1
myCamelCaseExample=2
Variable names should be meaningful – i.e. you should know straight away what they are referring to. And, as mentioned earlier in the week, use double quotes and curly braces to avoid issues with whitespace and wildcards in the variable value.

7. Prevent code bloat by using functions

Quite often you will want to repeat the same process multiple times within a script. We could just copy and paste the code, amending it to our needs. This is inefficient and will make for a longer script. Instead, we can wrap the code in a function so that the code we want to run is only located once in our script and is referenced using a function call anywhere it is needed. Like variables, functions should be meaningfully named. They should be small, only performing a limited, clear task.

8. Don’t duplicate scripts

It can be very tempting to put paths to data directly into your scripts. Then, when you come to use the script on the new dataset, copy the script, save it under a new name and update the dataset location to point at the new data. Please don’t do this. It’s worth investing the time in making your scripts reusable and scalable by using arguments and avoiding hard coded paths.

cartoon depicting one person at the computer and another looking at that Pearson's screen over her/his shoulder saying oh my god, with a tip: never look at somebody else's computer

Source: https://imgs.xkcd.com/comics/documents.png
9. Keep debugging simple

Scripts will fail, it’s inevitable. We’ve already mentioned expecting the unexpected, but what should you do when the unexpected actually happens? How do you know at what point your script failed? This is where logging comes in – print everything your script is doing back to the system. That way, when it fails, you know at which point it stopped working and will have a much easier time debugging the code.

10. Clean up after yourself

This is simple. If you generate intermediate or temporary files, make sure that you remove them when you’ve finished with them. This should be built into the script itself and not done an afterthought once it’s run.

11. Make your code easy to read

Digestible code is always easier to work with and maintain. To quote Martin Fowler: “Any fool can write code that a computer can understand. Good programmers write code that humans can understand.”. To help, you can use linters to look over your scripts and give feedback on their readability/formatting. One of the widely used online tools for Bash script linting is https://www.shellcheck.net/.

12. Don’t walk away from new scripts

It’s oh so tempting to hit “go”, set your nice, shiny new script running…then go off and make a cup of tea. Don’t. Sit back down and make sure it runs OK for the first couple of times. Why the first couple and not just the first time I hear you ask? Because, it will inevitably be the third or fourth time you run the script that it will fail in a spectacularly dramatic fashion. I like a cup of tea as much as the next person but, please, make it before you start running your script!


