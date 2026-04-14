# Flow–Recruitment FPCA Results — GitHub Pages Site

A static website that displays per-site flow–recruitment functional regression
results with a dropdown menu, plus a cross-site synthesis tab.

## Directory layout

```
your-repo/
├── docs/                          ← GitHub Pages serves from here
│   ├── index.html                 ← Per-site dropdown viewer
│   ├── cross_site.html            ← Cross-site synthesis tab
│   ├── style.css
│   ├── app.js
│   ├── sites_index.json           ← Auto-generated list of available sites
│   └── assets/
│       ├── sites/                 ← One subfolder per site
│       │   ├── BigHole.Melrose/
│       │   │   ├── metadata.json
│       │   │   ├── page1_panels.png
│       │   │   ├── page2_summary.png
│       │   │   └── page3_all_fpcs.png
│       │   └── Ruby.Vigilante/
│       │       └── ...
│       └── cross_site/            ← Cross-site analysis figures
│           ├── beta_correlation_matrix.png
│           ├── site_dendrogram.png
│           └── ...
├── WorkingFPCA.R
├── CrossSite_Hydrograph_Analysis.R
└── ...
```

## How to populate the site

### Step 1 — Per-site results

Run `WorkingFPCA.R` against each site as you already do. The updated script
automatically writes per-site PNGs and a `metadata.json` file to
`docs/assets/sites/<SITE>/` and rebuilds `docs/sites_index.json` after the
loop finishes. Required R packages: `qpdf`, `jsonlite` (install once with
`install.packages(c("qpdf", "jsonlite"))`).

### Step 2 — Cross-site figures

Run `CrossSite_Hydrograph_Analysis.R`. The script writes its figures to
`output/plots/cross_site/` by default — copy them into `docs/assets/cross_site/`
so the website can find them. (Or modify the `out_dir` variable in the
cross-site script to point directly at the docs folder.)

```bash
mkdir -p docs/assets/cross_site
cp output/plots/cross_site/*.png docs/assets/cross_site/
```

## How to deploy on GitHub Pages

1. **Create the repository.** On GitHub, create a new repo (public or private).
   On your machine:
   ```bash
   git init
   git add docs/ WorkingFPCA.R CrossSite_Hydrograph_Analysis.R README.md
   git commit -m "Initial commit: FPCA results viewer"
   git branch -M main
   git remote add origin https://github.com/<YOUR_USERNAME>/<REPO_NAME>.git
   git push -u origin main
   ```

2. **Enable GitHub Pages.** In the repo on GitHub:
   - Click **Settings** → **Pages** (left sidebar)
   - Under **Build and deployment**, set:
     - **Source:** Deploy from a branch
     - **Branch:** `main`, folder `/docs`
   - Click **Save**

3. **Wait ~1 minute**, then visit
   `https://<YOUR_USERNAME>.github.io/<REPO_NAME>/`
   You should see the dropdown viewer.

## Updating after re-running analyses

Each time you re-run `WorkingFPCA.R` or the cross-site script:

```bash
git add docs/
git commit -m "Update site results"
git push
```

GitHub Pages will rebuild within a minute.

## Notes

- The dropdown remembers the last site you viewed (via `localStorage`) and
  also reads the URL hash, so links like
  `https://.../#BigHole.Melrose` work for direct sharing.
- All figures are PNGs at 150 DPI, sized for fast loading on the web. If you
  need higher-resolution versions for publication, re-run the analyses with
  `dpi = 300` in the relevant `ggsave()` / `png()` calls.
- The multi-page PDF (`output/plots/LL_<SITE>_reduced_model_results.pdf`) is
  still generated as before, but each page is now sized to its own content
  instead of being padded to a uniform height. This eliminates the gaps you
  were seeing.

## Privacy / sharing

If your data are sensitive, make the GitHub repo **private** and use
GitHub Pages with the [private Pages feature](https://docs.github.com/en/pages/getting-started-with-github-pages/changing-the-visibility-of-your-github-pages-site)
(requires a GitHub Enterprise or Pro plan). Otherwise, anyone with the URL
can view the public Pages site.
