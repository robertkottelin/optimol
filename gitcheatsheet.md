
# GIT CHEATSHEET

## Initialize and Clone
git clone https://github.com/robertkottelin/optimol.git
cd optimol

---
## Branch Management
- Check out the main branch (assuming you're not already on it):
  git checkout main
- Create and switch to a new development branch named `dev`:
  git checkout -b dev
- List branches to see what's available:
  git branch
---

## Basic Operations
- Check the status of your repository:
  git status
- Add files to staging area (specific file or all files):
  git add <filename>
  git add .  # Adds all changes in the current directory and subdirectories
- Commit changes with a message:
  git commit -m "Your commit message here"
---

## Synchronization with Remote
- Fetch the latest changes from the remote repository:
  git fetch origin
- Pull changes from `main` into your current branch (e.g., updating `dev` with changes from `main`):
  git pull origin main
- Push changes from `dev` to the remote repository:
  git push -u origin dev  # -u sets upstream for the first time
---

## Merging
- Switch back to `main` to merge:
  git checkout main
- Pull the latest changes into `main` to ensure it's up-to-date before merging:
  git pull origin main
- Merge `dev` into `main`:
  git merge dev
- If you want to merge without a merge commit (assuming a fast-forward is possible), use:
  git merge --ff-only dev
- Resolve conflicts manually if any occur during the merge, then:
  git add <conflicted_files>
  git commit -m "Merge dev into main resolving conflicts"
- Push the merge to the remote `main`:
  git push origin main
---

## Rebasing
- Rebase `dev` branch onto the latest `main` branch:
  git checkout dev
  git pull origin main --rebase
- If conflicts occur during the rebase:
  1. Resolve conflicts manually in the affected files.
  2. Add the resolved files:
     git add <conflicted_file>
  3. Continue the rebase:
     git rebase --continue
- Abort a rebase (if needed):
  git rebase --abort
---

## Force Push and Backup
- Create a backup of your branch before rebasing or force pushing:
  git branch backup-<branch_name>
- After rebasing, push changes with force to overwrite remote history:
  git push origin dev --force-with-lease  # Safer option to avoid overwriting unexpected changes
  **Note**: Replace `--force-with-lease` with `--force` if you are confident about overwriting.
---

## Post-Merge and Cleanup
- Optionally, delete the `dev` branch if you're done with it:
  git branch -d dev  # Use -D to force delete if not merged
---

## Reset and Undo
- Undo the last commit (soft reset, keeping changes in your working directory):
  git reset --soft HEAD~1
- Undo the last commit and discard changes:
  git reset --hard HEAD~1
- Revert a specific commit while preserving history:
  git revert <commit_hash>
---

## Backup Strategy for Safety
- Create a new branch as a backup:
  git branch backup-main
- Push the backup branch to the remote repository:
  git push origin backup-main
- Restore from backup if needed:
  git checkout backup-main
  git checkout -b main  # Or any desired branch name
---

## Additional Tips
- See detailed commit history:
  git log --oneline --graph --all
- Temporarily save uncommitted changes (stash them):
  git stash
  Retrieve them later:
  git stash pop
---