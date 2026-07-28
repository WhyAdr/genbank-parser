# Git Workflow & Conflict Resolution Rules

1. **Verify Merge Conflict Resolution Before Staging**:
   During rebase or merge conflicts, NEVER run `git add` on conflicted files without explicitly removing all conflict markers (`<<<<<<<`, `=======`, `>>>>>>>`) and merging the content. Always verify with `Select-String -Pattern "<<<<<<<"` before staging and completing the rebase.

2. **Atomic Restructuring Commits**:
   When moving files or restructuring directories, stage the file additions and original file deletions together in a single commit. Avoid intermediate commits that introduce duplicate files.
