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
