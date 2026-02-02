# Contributing Guide

Do not commit directly to main. All changes must go through pull requests.

## Workflow

### Start new work
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

### Save your work
```bash
git add .
git commit -m "Description of changes"
git push -u origin feature/your-feature-name
```

### Create a pull request
```bash
1. Go to the main Github page
2. Create a pull request, add a description, and message Zach or whoever has context about to code to review.
```