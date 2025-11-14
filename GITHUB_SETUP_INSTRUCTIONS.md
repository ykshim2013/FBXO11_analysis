# GitHub Repository Setup Instructions

Your FBXO11 variant analysis repository has been initialized locally. Follow these steps to create the GitHub repository and push your code.

## Repository Details

- **Local path:** `/Users/ykshim2025/Desktop/Code2025/FBXO11_analysis`
- **Target GitHub URL:** `https://github.com/ykshim2013/FBXO11_analysis`
- **Initial commit:** ✅ Complete (20 files, 2074 lines)

---

## Step 1: Create GitHub Repository

### Option A: Via GitHub Website (Recommended)

1. Go to [https://github.com/new](https://github.com/new)
2. Fill in repository details:
   - **Repository name:** `FBXO11_analysis`
   - **Description:** Computational framework for analyzing structural destabilization in FBXO11 missense variants
   - **Visibility:** Public (for manuscript data sharing)
   - **Initialize:** ⚠️ **DO NOT** check "Add README" or any other files (repository already initialized locally)
3. Click **Create repository**

### Option B: Via GitHub CLI (Alternative)

If you have GitHub CLI installed:

```bash
cd /Users/ykshim2025/Desktop/Code2025/FBXO11_analysis
gh repo create ykshim2013/FBXO11_analysis --public --description "Computational framework for FBXO11 variant stability analysis" --source=. --remote=origin --push
```

If using GitHub CLI, **skip to Step 3** (repository automatically pushed).

---

## Step 2: Connect Local Repository to GitHub

After creating the empty GitHub repository, connect your local repository:

```bash
cd /Users/ykshim2025/Desktop/Code2025/FBXO11_analysis

# Add GitHub remote
git remote add origin https://github.com/ykshim2013/FBXO11_analysis.git

# Verify remote
git remote -v
# Expected output:
# origin  https://github.com/ykshim2013/FBXO11_analysis.git (fetch)
# origin  https://github.com/ykshim2013/FBXO11_analysis.git (push)
```

---

## Step 3: Push Repository to GitHub

### Configure Git User (If Not Already Set)

```bash
# Set your name and email for this repository
git config user.name "Youngkyu Shim"
git config user.email "ykshim2013@gmail.com"

# Or set globally for all repositories
git config --global user.name "Youngkyu Shim"
git config --global user.email "ykshim2013@gmail.com"
```

### Push to GitHub

```bash
cd /Users/ykshim2025/Desktop/Code2025/FBXO11_analysis

# Push main branch to GitHub
git push -u origin main
```

**Expected output:**
```
Enumerating objects: 29, done.
Counting objects: 100% (29/29), done.
Delta compression using up to 8 threads
Compressing objects: 100% (24/24), done.
Writing objects: 100% (29/29), XX.XX KiB | X.XX MiB/s, done.
Total 29 (delta 2), reused 0 (delta 0)
To https://github.com/ykshim2013/FBXO11_analysis.git
 * [new branch]      main -> main
Branch 'main' set up to track remote branch 'main' from 'origin'.
```

### Authentication

If prompted for authentication:

**Option 1: Personal Access Token (Recommended)**
1. Go to GitHub → Settings → Developer settings → Personal access tokens → Tokens (classic)
2. Click "Generate new token (classic)"
3. Select scopes: `repo` (full control of private repositories)
4. Copy the token
5. Use token as password when prompted during `git push`

**Option 2: GitHub CLI**
```bash
gh auth login
```
Follow interactive prompts to authenticate.

---

## Step 4: Verify Repository Upload

1. Visit [https://github.com/ykshim2013/FBXO11_analysis](https://github.com/ykshim2013/FBXO11_analysis)
2. Check that all files are present:
   - [x] README.md with project overview
   - [x] LICENSE (MIT)
   - [x] CITATION.cff
   - [x] requirements.txt
   - [x] data/ directory structure
   - [x] scripts/ directory with analysis examples
   - [x] docs/ with methods and data dictionary
3. Verify badges display correctly (DOI will show "pending" until publication)
4. Test that relative links work (click links in README.md)

---

## Step 5: Update Manuscript with DOI (After Publication)

Once your manuscript is accepted and you create a Zenodo DOI:

1. Edit README.md and CITATION.cff with actual DOI
2. Commit changes:
   ```bash
   git add README.md CITATION.cff
   git commit -m "Update DOI after manuscript publication"
   git push
   ```

3. Create GitHub release:
   ```bash
   # Tag release matching manuscript publication
   git tag -a v1.0.0 -m "Release v1.0.0: Manuscript publication"
   git push origin v1.0.0
   ```

---

## Step 6: Optional Enhancements

### Add GitHub Actions for Continuous Integration

Create `.github/workflows/test.yml` to automatically run analysis scripts:

```yaml
name: Analysis Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      - run: pip install -r requirements.txt
      - run: python scripts/statistical_analysis/pathogenic_vs_benign.py --test
```

### Enable GitHub Pages for Documentation

1. Go to repository Settings → Pages
2. Source: Deploy from branch `main`, folder `/docs`
3. Your methods documentation will be available at: `https://ykshim2013.github.io/FBXO11_analysis/`

### Add Repository Topics

Add topics to improve discoverability:
- `fbxo11`
- `variant-analysis`
- `protein-stability`
- `computational-biology`
- `molecular-dynamics`
- `foldx`
- `rosetta`
- `gromacs`
- `alphafold`
- `neurodevelopmental-disorder`

Settings → Topics → Add topics

---

## Repository Structure Verification

After pushing, your repository should have this structure:

```
FBXO11_analysis/
├── .gitignore                         ✅ Created
├── LICENSE                            ✅ Created (MIT)
├── CITATION.cff                       ✅ Created
├── README.md                          ✅ Created (389 lines)
├── requirements.txt                   ✅ Created
│
├── data/                              ✅ Created
│   ├── foldx/
│   │   ├── README.md                  ✅ Created
│   │   ├── aggregated_results.csv     ✅ Placeholder
│   │   └── by_variant_summary.csv     ✅ Placeholder
│   ├── rosetta/
│   │   ├── README.md                  ✅ Created
│   │   ├── aggregated_results.csv     ✅ Placeholder
│   │   └── by_variant_summary.csv     ✅ Placeholder
│   ├── md_trajectories/
│   │   ├── README.md                  ✅ Created
│   │   ├── download_instructions.md   ✅ Created
│   │   ├── trajectory_summary.csv     ✅ Placeholder
│   │   └── md_analysis_results.csv    ✅ Placeholder
│   └── alphafold3_structures/
│       └── README.md                  ✅ Created
│
├── scripts/                           ✅ Created
│   ├── statistical_analysis/
│   │   └── pathogenic_vs_benign.py    ✅ Created
│   └── md_simulations/
│       └── 03_analyze_trajectories.py ✅ Created
│
├── analysis_outputs/                  ✅ Created (empty directories)
│   ├── figures/
│   ├── tables/
│   └── statistics/
│
└── docs/                              ✅ Created
    ├── METHODS.md                     ✅ Created (comprehensive)
    └── DATA_DICTIONARY.md             ✅ Created
```

---

## Troubleshooting

### Issue: "remote: Repository not found"
**Solution:** Verify repository name and that it exists on GitHub. Check that username is `ykshim2013`, not `ykshim2025`.

### Issue: "Authentication failed"
**Solution:** Use Personal Access Token instead of password. Generate token at GitHub → Settings → Developer settings → Personal access tokens.

### Issue: "Updates were rejected because the remote contains work"
**Solution:** This shouldn't happen if you created an empty repository. If it does:
```bash
git pull origin main --allow-unrelated-histories
git push origin main
```

### Issue: Large files rejected (>100 MB)
**Solution:** Placeholder CSV files are empty, so this shouldn't occur. If you add real data later:
```bash
# Use Git LFS for large files
git lfs install
git lfs track "data/**/*.xtc"
git add .gitattributes
git commit -m "Add Git LFS tracking for trajectory files"
```

---

## Next Steps

Once repository is pushed:

1. **Add actual data:** Replace placeholder CSV files with real computational results
2. **Test scripts:** Verify analysis scripts work with real data
3. **Update manuscript:** Ensure Data Availability Statement URL is correct
4. **Create Zenodo archive:** Link GitHub repository to Zenodo for DOI
5. **Share with collaborators:** Provide repository URL for review

---

## Questions?

- **Git issues:** Check [GitHub documentation](https://docs.github.com/)
- **Repository structure:** See `docs/METHODS.md` and `docs/DATA_DICTIONARY.md`
- **Data sharing:** See `data/md_trajectories/download_instructions.md`

---

**Repository created:** 2025-11-14
**Initial commit hash:** a029301
**Status:** ✅ Ready for GitHub push
