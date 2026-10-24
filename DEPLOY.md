# Tablet and Render deployment

This is the same HTML frontend and Python numerical engine as the local version. No tablet-side Python is required when hosted.

## Render
1. Create a NEW repository with the contents of this folder (do not overwrite the FMM Workbench).
2. In Render, select **New > Blueprint**, connect the repository, and apply `render.yaml`.
3. After deployment, open the service URL on your tablet. The `/health` endpoint should return `status: ok`.
4. Run a small least-squares test before using stochastic optimizers.

**No deployment has been performed by the assistant.** The Blueprint requests a paid Starter instance, not a free instance; inspect Render's current pricing before confirming. Long fits execute synchronously and may time out. No authentication or job queue is included; use private access or add authentication and rate limiting before public exposure. Do not upload sensitive datasets.

## Local mode
Windows: `START_WINDOWS.bat`; Mac/Linux: `bash START_MAC_LINUX.sh`. Open the URL printed by the launcher on the same computer. To open from another device on your private Wi-Fi, you must configure network binding, firewall and authentication; the default localhost binding intentionally blocks other devices.

## Scope
The four built-in R6G cases, joint B-spline fitting and five optimizers remain available. Swanepoel extraction, Lorentz/KK fitting, experimental upload, rigorous uncertainty analysis, and experimental design remain future work. Plotly charts require an internet connection to the Plotly CDN.
