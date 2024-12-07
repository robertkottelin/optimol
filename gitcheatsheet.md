# GIT CHEATSHEET

## Initialize and Clone
git clone https://github.com/robertkottelin/optimol.git
cd optimol

## Branch Management
Check out the main branch (assuming you're not already on it):
git checkout main
Create and switch to a new development branch named dev:
git checkout -b dev
List branches to see what's available:
git branch

## Basic Operations
Check the status of your repository:
git status
Add files to staging area (specific file or all files):
git add <filename>
git add .  # Adds all changes in the current directory and subdirectories
Commit changes with a message:
git commit -m "Your commit message here"

## Synchronization with Remote
Fetch the latest changes from the remote repository:
git fetch origin
Pull changes from main into your current branch (assuming you want to update dev with changes from main):
git pull origin main
Push changes from dev to the remote repository:
git push -u origin dev  # -u sets upstream for the first time

## Merging
Switch back to main to merge:
git checkout main
Pull the latest changes into main to ensure it's up-to-date before merging:
git pull origin main
Merge dev into main:
git merge dev
If you want to merge without a merge commit (assuming a fast-forward is possible), use:
git merge --ff-only dev
Resolve conflicts manually if any occur during the merge, then:
git add <conflicted_files>
git commit -m "Merge dev into main resolving conflicts"
Push the merge to the remote main:
git push origin main

## Post-Merge
Optionally, delete the dev branch if you're done with it:
git branch -d dev  # Use -D to force delete if not merged

## Additional Tips
Undo a commit (soft reset, keeping changes in your working directory):
git reset --soft HEAD~1
See commit history:
git log
Switch between branches without losing changes:
git stash
git checkout <branch>
git stash pop