ATP Analysis (Berthold tube luminometer)


A small, reusable toolkit for quantifying ATP from Berthold tube luminometer runs. It:

	•	Fits the ATP standard curve first (on raw standards),
 
	•	Models Tris baseline drift over the run (single or two-segment),
 
	•	Treats Tris as a time-varying intercept (so we don’t double-penalize samples),
 
	•	Corrects samples by subtracting the predicted Tris at their timestamps,
 
	•	Converts corrected signal to concentration using the standard-curve slope,
 
	•	Optionally normalizes to original sample mass/volume,
 
	•	Auto-creates output folders and exports plots + CSVs,
 
	•	Provides an example notebook for reproducible runs.

What the code does (pipeline)

	1.	Integrate each sample’s RLU trace over Time (trapezoid).
 
	2.	Fit standards first (raw integrals → slope & intercept stored).
 
	3.	Fit Tris drift (single or two-segment line) vs seconds since start.
 
	4.	Correct samples by subtracting predicted Tris at each timestamp; clip negatives to 0.
 
	5.	Quantify with slope only: Concentration_ng_per_mL = Corrected / slope.
 
	6.	(Optional) Back-calculate to Original_Sample_Concentration_ng_per_<unit> using sample volumes/weights.
 
	7.	(Optional) Blanks threshold: mark samples above mean_blank + k·SD_blank.
 
	8.	(Optional) Merge metadata and plot against any column.
 
	9.	Save plots and CSVs; directories are created automatically.

Packages needed:
numpy 
matplotlib
pandas
scipy 
```
pip install numpy
pip install matplotlib
pip install pandas
pip install scipy
```
Or use a conda environment, if you prefer. 

As of 10/28/24, we have not been able to program the instrument software to output the time of measurement, so the time each tube goes into the instrument **must** be recorded. 

After exporting the data, the easiest thing to do is copy the data into a fresh csv file and change the top row to the sample names. The standards can be just the numbers and name the Tris measurments "Tris": 
<img width="784" alt="Screenshot 2024-09-26 at 8 51 08 AM" src="https://github.com/user-attachments/assets/1963d83f-f545-4f4e-9174-6f032fda9309">

Please also ensure that your sample names have a letter or are a string so that the code can easily parse them out. 

The other data to input is the time of measurement csv file, which has the sample names on the top row and the recorded time of measurment in the bottom row: 
<img width="651" alt="Screenshot 2024-09-26 at 8 50 16 AM" src="https://github.com/user-attachments/assets/f71e6c28-450e-4df4-a03b-fed245855158">

Update June 18 2025: Included a "atp_time_logger" app that opens a GUI allowing you to just click to record time of measurement given a python list of samples and then save to CSV

First, the code reads in the data and integrates the RLU values using the trapezoid approximation. Next, it combines those values with the time and converts the time into a python-readable format. 
The next function takes 4 arguments: the input dataframe, a True/False indicator whether to seperate out the Tris values, and the x bounds. The first xbound is the upper bound of the first line and second is the lower bound of the second line. 
This is because the firefly powder degrades over the course of the experiment, so it is important to track that degradation. See the difference in the slopes of RLU vs time of the Tris measurements below when you don't seperate them out: 

![482161df-b9bb-49ce-9f41-7ed0f8b7a30a](https://github.com/user-attachments/assets/1b39257e-2810-44b9-b93b-4b4996a6c76c)
![ab3db168-69da-44db-8aea-2292e77bf454](https://github.com/user-attachments/assets/e6af0e05-c0e2-4306-8133-c03bfdd18647)

Next it takes the linear equations for these, calculates the predicted Tris luminescence for the time of the sample measurement (based on the slope using the xbounds) and subtracts that from the integrated luminescence. It then sets negative values to 0. 

The next function (extract_concentration) parses samples and blanks from the standards and then puts the numeric value of the standard concentration into a new column. 

The final function (calculate_sample_concentration) generates a standard curve from the standards (while averaging triplicate measurements) and then uses that linear function to calculate the sample concentration to calculate the sample concentration in ng/mL

Example standard curve: 

![b0c76198-7a78-4d1f-9678-d273f135e208](https://github.com/user-attachments/assets/70dff8dc-c51d-4bae-acfa-5cae3e23dc63)


Please provide in the function the sample extract volume (e.g., 5 mL) and the amount of sample added to the extract. If you do not provide a sample volume (e.g., weight of sample, volume of sample, number of cells, etc), it will just provide the total ATP in the extract. 

This code will output the csv file with the calculated ATP values, the standard curve figure, and the Tris vs Time figure. 

Quickstart: 

```
from atp_analysis import ATPAnalyzer
import pandas as pd

base = "example_run/"  # change to your folder with CSVs

analyzer = (
    ATPAnalyzer(data_csv=base + "measurements.csv", time_csv=base + "timestamps.csv")
      .integrate()
      .fit_standard_curve()                            # Fit standards FIRST (raw integrals)
      .fit_tris_drift(separate=True, split_seconds=1000)  # Or separate=False for a single line
      .apply_corrections_and_quantify(
          extract_vol=4.0,                              # mL of Tris (or diluent) used
          sample_vol=pd.read_csv(base + "sample_volumes.csv"),  # None | float | DataFrame
          sample_unit="g"                               # labels output as ng per g (or "mL")
      )
      .compute_blank_threshold(
          blank_names=['Blank','Blank.1','Blank 2'],    # put your actual blank names here
          k=3.0                                         # threshold = mean_blank + k*SD_blank
      )
      .merge_metadata(
          meta=base + "metadata.tsv",                   # or .csv
          right_key="#SampleID"                         # column in metadata matching Base_Sample
      )
)

# Plots (folders auto-created)
analyzer.plot_standard_curve(save_path=base + "outputs/ATP_Standard_Curve.png", through_origin=True)
analyzer.plot_tris_drift(save_path=base + "outputs/Tris_drift.png")  # shows 1 or 2 lines as fitted
analyzer.plot_vs_metadata(x="Age", y="avg_concentration", save_path=base + "outputs/ATP_vs_age.png")

# CSVs (folder auto-created)
analyzer.save_outputs(prefix=base + "outputs/atp_")
```
Outputs

	•	*_integrals.csv — integral and timestamp per column.
 
	•	*_samples_wide.csv — aliquot-level concentrations (extract and/or original-sample units).
 
	•	*_grouped.csv — mean ± SD per Base_Sample (+ n).
 
	•	*_merged_meta.csv — grouped joined with metadata (if provided).
 
	•	ATP_Standard_Curve.png — standards mean±SD + fit (toggle through-origin for visualization).
 
	•	Tris_drift.png — Tris points + fitted drift line(s).
 
	•	ATP_vs_<metadata>.png — concentration vs chosen metadata field.
