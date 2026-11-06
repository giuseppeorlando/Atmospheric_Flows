## Git Branching Structure and Workflow

This section describes the Git branching and tag structure used when developing. In the following we differentiate between **branches in the upstream repository** and **branches in developers' forks**. Forks are used to develop features, bugfixes etc. They are created from the upstream's `develop` branch by forking to a local repo.

### Branching Structure

#### Branches in the Upstream Repository
These branches are centralized and exist in the upstream repository. The branches are never deleted.

**`master` branch**:
   - Represents the stable codebase.
   - Contains tagged releases (e.g. `v0.1`, `v0.2`).

**`develop` branch**:
   - Serves as the main branch for development integration.
   - Receives completed features and bug fixes from developers' forks.

**`legacy-support/vX.Y` branches** (used only when necessary):
   - Created for maintaining older versions (e.g. `v0.1`) when newer releases (e.g. `v0.2`) are already available.
   
#### Branches in Developers' Forks
These branches are created and reside in developers' personal forks of the upstream repository. They are temporary and merged back into the upstream repository when work is completed.

**`feature/<name>` branches**:
   - Created for new feature development.
   - Based on the `develop` branch.
   - Merged back into the upstream's `develop` branch upon completion.

**`release/vX.Y.Z` branches**:
   - Created to prepare a new release.
   - Based on the `develop` branch.
   - Merged back into the upstream's `master` and `develop` branch upon completion.

### Workflow

If you prefer to branch out from your local fork's `develop`, remember to synchronize your fork's `develop` with `upstream/develop` before starting any of the below workflows.

To be able to work with branches in the upstream repository from within your local fork, you can add the upstream as a remote:
```bash
git remote add upstream <upstream_repository_URL>
```
This is e.g. relevant, when tagging a new release on the upstream's `master` branch (see below), or if you prefer to branch out from the upstream's `develop` branch. In the latter case, please replace `origin/develop` with `upstream/develop` in below commands.

#### Developing New Features
  1. In your fork, create a new `feature/<name>` branch based on `develop`:
     ```bash
     git checkout -b feature/<name> origin/develop
     ```
  2. Implement and commit your changes.
  3. Push the branch to your fork and create a pull request to merge it into the upstream's `develop` branch. Follow the instructions in the PR.

#### Creating a New Release
1. Ensure your local fork is up to date with the upstream repository.
2. In your fork, create a new `release/vX.Y.Z` branch based on `develop`:
    ```bash
    git checkout -b release/vX.Y.Z origin/develop
    ```
3. Update the `CHANGELOG.md` with a new **[vX.Y.Z] YYYY.mm.dd** section containing all changes that are released.
4. Commit the changes:
    ```bash
    git add CHANGELOG.md 
    git commit -m "Prepare release vX.Y.Z"
    ```
5. Push the branch to your fork:
    ```bash
    git push origin release/vX.Y.Z
    ```
6. Create two new pull requests to merge the `release/vX.Y.Z` branch into the upstream's `master` and `develop` branch, respectively. Follow the instructions in the PRs.

7. Once the above pull requests have been approved and merged, create a new tag on the `master` branch, while still in your own fork:
    ```bash
    git checkout master
    git pull upstream master
    git tag -a vX.Y.Z -m "Release vX.Y.Z"
    git push upstream vX.Y.Z
    ```
8. Create a new release in GitHub
    1.  The release shall be based on the just created tag.
    2.  Click the 'Releases' section to the right on the repository's homepage near the bottom.
    3.  Click 'Draft a new release'.
    4.  Choose the relevant tag and set target branch to 'master'.
    5.  Set the name of the tag as the title of the release.
    6.  Copy and paste the changelog items relevant to this release into the description of the release.
    7.  Set the release as the latest release.
    8.  Publish the release.



### Summary of Branch Responsibilities

| **Branch Type**              | **Location**               | **Purpose**                                                                                       |
|------------------------------|----------------------------|---------------------------------------------------------------------------------------------------|
| `master`                     | Upstream                   | Holds stable releases, tagged with version numbers.                                     |
| `develop`                    | Upstream                   | Main branch for integrating features and bug fixes.                                               |
| `legacy-support/vX.Y`        | Upstream                   | For maintaining older versions and applying critical fixes to legacy releases.                          |
| `feature/<name>`             | Developers' Fork           | For developing new features. Based on `develop`.                                                  |
| `release/vX.Y.Z`          | Developers' Fork           | To prepare a new release. Based on `develop`.|
