# Contributing to MIRCO

Thank you for your willingness to contribute to MIRCO.  The procedure to do so is the following:

## Create a GitHub Issue

Navigate to MIRCO' [GitHub Issues page](https://github.com/imcs-compsim/MIRCO/issues) and create a new issue.  The issue can be used for any number of things &mdash; reporting a bug, suggesting an enhancement, posing a question, etc.  On the new issue creation page, you'll notice the *Description* field will be pre-populated with some text.  Follow the instructions in that template to give us as much information as you can such that we can tackle the issue as soon as is practicable.

## Work an Issue

When work is ready to commence on an issue, the workflow to use is the following:

### Fork MIRCO

* If you have not already done so, create a fork of MIRCO on GitHub under your username.
* Clone your fork of MIRCO with `git clone git@github.com:<username>/MIRCO`.
* Each time you clone your fork, `git remote add upstream git@github.com:imcs-compsim/MIRCO` to add the original MIRCO repository as the `upstream` remote.

### Create a Feature Branch

Create a local branch with branch name `<brancName>` off of `master` on which to make your changes.  `<branchName>` can be whatever you like, though we have some recommendations:

* Include the issue number in it in some way, for instance, `123-<restOfBranchName>`, or `<restOfBranchName>-123`.
* Make the branch name descriptive; that is, avoid `fixSomeStuff`, `performanceTweaks`, and generic names along those lines.

### Make Your Changes

Do whatever work is necessary to address the issue you're tackling, breaking your work into logical, compilable commits.  Feel free to commit small chunks early and often in your local repository and then use `git rebase -i` to reorganize your commits before sharing.  Make sure the commit messages you will be sharing reference the appropriate GitHub issue numbers.

### Create a Pull Request

When your changes are ready to be integrated into MIRCO's `master` branch:

* Push your local feature branch up to your fork with `git push -u origin <branchName>`.
* Navigate to your fork of MIRCO on GitHub and create a new pull request:
  * Be sure you choose:
    * base fork:  `imcs-compsim/MIRCO`
    * base:  `master`
    * head fork:  `<username>/MIRCO`
    * compare:  `<branchName>`
  * On the new pull request creation page, you'll notice the *Description* field will be pre-populated with some text.  Follow the instructions in that template to give us as much information as you can such that we can review and approve the issue as soon as is practicable.

### Feedback

At this point you'll enter into a stage where you and various MIRCO developers will iterate back and forth until your changes are in an acceptable state and can be merged in.  If you need to make changes to your pull request, make additional commits on your `<branchName>` branch and push them up to your fork.  Make sure you don't delete your remote feature branch or your fork of MIRCO before your pull request has been merged.
