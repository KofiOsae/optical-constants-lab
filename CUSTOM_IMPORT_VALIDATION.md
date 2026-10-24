# Custom-spectrum import validation
- Implemented: CSV, TSV, TXT and JSON column-table ingestion; explicit/automatic wavelength unit selection; explicit R/T percent or fractional units; numerical validation, sorting, duplicate rejection and energy plausibility check; custom spline inverse API and browser plots.
- Verified: two parser smoke tests executed successfully (CSV µm + percent; JSON nm + fractions).
- Not verified: Flask API and browser end-to-end tests because Flask is absent in this execution environment and package installation cannot access the network.
- Not deployed: this source archive must be published to the standalone lab's repository and Render service.
- Scientific limitation: one spectrum may not uniquely constrain n, k and d; unknown thickness must be initialized and can only vary ±15% in this implementation.
