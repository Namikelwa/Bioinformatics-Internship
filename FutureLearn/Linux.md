# Capturing user input


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

* echo options
man echo
https://www.gnu.org/software/bash/manual/html_node/Escape-Character.html
https://www.gnu.org/software/bash/manual/html_node/Escape-Character.htm
