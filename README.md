# README #

This README briefly describes the workflow for the QUant code review

### What is this repository for? ###

This repo is to enable collaboration on code review and optimization of the QUant package.

### How do I get set up? ###

* To run the code, you start with `simulate_uncertainty.m`
* You will want to update the pathnames and filenames at the start of this m file to reflect where the data files are that you want to process
* Ensure that the folders `general` and `tools` as well as `Discharge_sm.m`, `MC_Data.m`, `OriginData_sm.m` and `clsQComp.m` are in the same folder as `simulate_uncertainty.m`
* The code will prompt the user to run `extrap` – but I don’t have the Curve Fitting Toolbox (yet) which is needed to run 'extrap' within matlab. I have this line commented out, and instead simply load a saved 'extrap' summary from a previous run by Stephanie. (as well, you can manually enter the 'extrapDev % deviation from the default extrapolation' if you don’t have the curve fitting toolbox (which is what I did on line 409).
* You can choose to run multiple transects or just one (modify the for loop on line 449)
* The code is set to determine the uncertainty when all parameters are varied and then to compute the uncertainty for when each parameter is varied. You can easily modify this if  you want to look at one parameter only, or only all the parameters (lines 773-780).
* You can vary the number of Monte Carlo realizations on line 1095 (we used 1000 for the results in the paper – but now that I think about it, this seems low….I’m not sure whether Stephanie did some sort of convergence test or not with a higher number of realizations.)
* Before the simulations are run, a figure pops up to display the parameter settings (mean and standard deviations for each parameter, and the input from the mmt file). The user can modify these values (or not) and then hit continue. We would probably want to turn this feature off if we want to do some batch runs.

### Contribution guidelines ###

* There is a folder called `code-review-notes` that contains notes on the various aspects of the review process. Please check these notes before making changes to the code
* Changes to the code need to be orchestrated via [pull requests](https://www.atlassian.com/git/tutorials/making-a-pull-request/)

Before you submit your pull request consider the following guidelines:

* Make your changes in a new git branch:

     ```
     git checkout -b my-fix-branch master
     ```

* Commit your changes using a descriptive commit message 
     ```
     git commit -a
     ```
  Note: the optional commit `-a` command line option will automatically "add" and "rm" edited files.

* Push your branch to Bitbucket:

    ```
    git push origin my-fix-branch
    ```

* In Bitbucket, create a new [pull request](https://bitbucket.org/frank-engel/quant/pull-request/new)
* If we suggest changes then:
  * Make the required updates.
  * Rebase your branch and force push to your Bitbucket repository (this will update your Pull Request):

    ```
    git rebase master -i
    git push origin my-fix-branch -f
    ```
That's it! Thank you for your contribution!

### Who do I talk to? ###

* Elizabeth Jamieson (elizabeth.jamieson@canada.ca) is the repo admin. Stephanie Moore (stephanie.moore2@canada.ca) is the first author of this code.
