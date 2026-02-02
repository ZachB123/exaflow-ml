# Contributing Guide

Do not commit directly to main. All changes must go through pull requests.

## Workflow

### Starting new work
```bash
git checkout main
git pull origin main
git checkout -b prefix/your-branch-name
```

Branch prefixes:
- `feature/` - new features
- `bug/` - bug fixes
- `test/` - adding tests
- `tweak/` - minor changes

### Saving your work
```bash
git add .
git commit -m "Description of changes"
git push -u origin feature/your-feature-name
```

### Creating a pull request
1. Go to https://github.gatech.edu/VIP-CFD/real-time-inferencing
2. Click "Compare & pull request"
3. Add description and request a reviewer
4. Click "Create pull request"

### After approval
Click "Merge pull request" then "Delete branch"

```