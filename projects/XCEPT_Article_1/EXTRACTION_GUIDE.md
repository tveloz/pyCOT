# XCEPT Article 1 - Extraction Guide

This guide explains how to extract this project into a separate Git repository.

## Overview

The XCEPT Article 1 project is a conflict resolution research project that uses the pyCOT library. To make it a standalone repository, follow these steps.

## Step 1: Create the New Repository

1. Create a new repository on GitHub (e.g., `XCEPT-Conflict-Resolution`)
2. Clone it locally:
   ```bash
   git clone https://github.com/YOUR_USERNAME/XCEPT-Conflict-Resolution.git
   cd XCEPT-Conflict-Resolution
   ```

## Step 2: Copy Project Files

Copy the contents of this project folder to the new repository:

```bash
# From pyCOT root directory
cp -r projects/XCEPT_Article_1/* /path/to/XCEPT-Conflict-Resolution/
```

## Step 3: Create pyproject.toml

Create a `pyproject.toml` in the new repository:

```toml
[tool.poetry]
name = "xcept-conflict-resolution"
version = "0.1.0"
description = "Conflict Resolution Research using Chemical Organization Theory"
authors = ["Your Name <your.email@example.com>"]
readme = "README.md"

[tool.poetry.dependencies]
python = "^3.11"
pyCOT = { git = "https://github.com/YOUR_USERNAME/pyCOT.git", branch = "main" }
numpy = "^1.25.0"
matplotlib = "3.8"
pandas = "^2.1"
scipy = "^1.14.0"

[build-system]
requires = ["poetry-core"]
build-backend = "poetry.core.masonry.api"
```

## Step 4: Update Imports (if needed)

The scripts already use absolute imports from pyCOT. If you installed pyCOT via poetry, no changes are needed.

For local development alongside pyCOT:
```bash
cd XCEPT-Conflict-Resolution
poetry add --editable ../pyCOT
```

## Step 5: Create README

Create a `README.md` explaining the project.

## Step 6: Commit and Push

```bash
git add .
git commit -m "Initial commit: XCEPT Conflict Resolution project"
git push origin main
```

## Granting Access to Collaborators

### On GitHub:
1. Go to your repository: `https://github.com/YOUR_USERNAME/XCEPT-Conflict-Resolution`
2. Click **Settings** > **Collaborators and Teams**
3. Click **Add people** or **Add teams**
4. Enter the collaborator's GitHub username or email
5. Choose a role:
   - **Read**: Can view and clone
   - **Triage**: Can manage issues
   - **Write**: Can push changes
   - **Maintain**: Can manage settings (except sensitive)
   - **Admin**: Full access

### Via Git:
Collaborators with write access can:
```bash
git clone https://github.com/YOUR_USERNAME/XCEPT-Conflict-Resolution.git
cd XCEPT-Conflict-Resolution
# Make changes...
git push origin main
```

## Keeping pyCOT Updated

When pyCOT is updated:
```bash
cd XCEPT-Conflict-Resolution
poetry update pyCOT
```

Or to use a specific version/tag:
```toml
pyCOT = { git = "https://github.com/YOUR_USERNAME/pyCOT.git", tag = "v1.0.0" }
```

## Project Structure

```
XCEPT-Conflict-Resolution/
├── scripts/                    # Analysis scripts
│   ├── game_two_budget.py
│   ├── script_calibrated_conflict.py
│   ├── script_Intervention_Cost.py
│   ├── script_multiscale_conflict.py
│   ├── script_print_conflict_organizations.py
│   ├── script_semantic_partition.py
│   ├── script_transition_resolution.py
│   ├── script_two_budget_adaptive.py
│   ├── script_two_budget_game.py
│   ├── test_seasonal_climate.py
│   └── diagnostic_*.py
├── outputs/                    # Generated outputs
│   ├── visualizations/
│   ├── intervention_multiscale/
│   ├── resilience_trust_multiscale/
│   └── *.csv, *.png
├── reports/                    # Analysis reports
├── pyproject.toml
└── README.md
```

## Conflict Theory Network Data

The conflict theory networks are in the main pyCOT repository at:
```
pyCOT/data/conflict_theory/
```

You may want to copy these to your XCEPT repository or reference them from pyCOT.
