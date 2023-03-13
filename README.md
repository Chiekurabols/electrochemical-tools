# electrochemistry_tools

## Introduction
This repository contains classes, collection of scripts and jupyter notebook to facilitate performing electrochemical calculations and analyse the output. 

## Prerequisites
If using mac, start with installing [Homebrew](https://brew.sh).
Then install GNU sed:
```
$ brew install gnu-sed
```
And then install git:
```
$ brew install git
```

## Getting started
Begin with cloning this git repository to your machine with the following command:

```
$ git clone https://git.uni-due.de/adk103u/electrochemistry_tools.git ./
```

Your login details should be your unnikennung and the respective password. You can ommin `./` if you want the subdirectory "electrochemistry_tools" to be created. 

Also, for the example scripts to work, you need to add "/path/to/electrochemistry_tools/internal_classes" subdirectory to `PYTHONPATH`, which can be done, for example, with the following command:

```
$ sed -i '$a export PYTHONPATH=/path/to/electrochemistry_tools/internal_classes:$PYTHONPATH' ~/.zshrc
source ~/.zshrc
```

Do not forget to change `/path/to` to the actual path on your machine.

Now you should be ready to use the provided tools!
