## Collaborating with Git: Git Remotes, git push, and git pull

### Creating a Shared Central Repository with Github

### Authenticating with Git Remotes

### Connecting with Git Remotes: gir remote

Let us configure our local repository to use the GitHub repository we have just created as a remote repository.
```{}
git remote add orgin git@github.com:usename/zmays-snps.git
```

By convention, *origin* is the name of your primary remote repository. We can push commits and fectch commits from it.

### Pushing Commits to a Remote Repository with git push

subcommand: git push <remote-name> <branch>
  
### Pulling Commits to a Remote Repository with git push

```{}
git pull origin master
````

### Working with your collaborators: Pushing and Pulling

Because new commits builf on top of the commit history, it is helpful to do the following to avoid problems:
  
- When pulling in changes, it helps to have your project's changes committed.

- Pull often.

### Merge Conflicts

strategy to solve them:
1. Use git status to find the conflicting file(s).
2. Open and edit those files manually to a version that fixes the conflict
3. Use git add to tell Git that you've resolved the conflict in a particular file
4. Once all conflicts are resolved, use git status to check that all changes are staged. Then, commit the resolved versions of the coflicting file(s).
### More Github Workflows: Forking and Pull Requests

By forking another person's Github repository, you are copying their repository to your own GitHub account. You can then clone your forked version and continue development in you own repository.
If you decide that you have made changes you want to share with the main repository, you can request that your commits are pulled using a pull request.

