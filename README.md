# Daisyworld (Shinylive)

An interactive Daisyworld model that runs **entirely in the browser** via
Shinylive (Shiny for R compiled to WebAssembly with webR). There is no server
to run or pay for: every visitor executes the app locally in their own browser,
so it scales to a whole class at no cost.

## Live app

Once deployed (see below), the app is at:

```
https://<your-github-username>.github.io/<your-repo-name>/
```

## Repository layout

```
.
├─ app/
│  └─ app.R                     # the Shiny app (single file)
├─ .github/
│  └─ workflows/
│     └─ deploy.yml             # builds + publishes to GitHub Pages on each push
└─ README.md
```

## How deployment works

On every push to `main`, the workflow in `.github/workflows/deploy.yml`:
1. sets up R on a GitHub build server,
2. installs the `shinylive` package,
3. runs `shinylive::export("app", "_site")` to compile the app and download the
   WebAssembly builds of its R package dependencies, and
4. publishes the resulting static site to GitHub Pages.

The app depends only on packages with webR builds: shiny, bslib, ggplot2,
dplyr, tidyr, tibble, readr, scales, DT.

## Editing the app

Change `app/app.R`, commit, and push to `main`. The site rebuilds automatically
in about 2–4 minutes. Watch progress under the repo's **Actions** tab.

## Notes for users

- **First load is slow** (roughly 10–30 seconds): the browser downloads webR
  and the package builds, then caches them. Later visits are fast.
- The CSV download buttons work without question. The three PNG plot-download
  buttons rely on webR's graphics device; test them once after your first
  deploy (see the deployment guide).
