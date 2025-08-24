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

The following contains a python script to analyze data output from the Berthold tube luminometer instrument in the lab. 
As of 10/28/24, we have not been able to program the software to output the time of measurement, so the time each tube goes into the instrument **must** be recorded. 

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
