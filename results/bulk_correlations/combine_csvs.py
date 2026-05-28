import os
import glob
import pandas as pd

# 1. Look for all CSV files in the current directory
csv_files = glob.glob("*filtered.csv")

# Quick check: did we find anything?
if not csv_files:
    print("No CSV files found in this folder!")
else:
    # 2. Create an Excel writer object
    output_excel = "combined_output_filtered.xlsx"
    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:

        for file in csv_files:
            # Extract just the filename without the '.csv' extension for the sheet name
            sheet_name = os.path.splitext(os.path.basename(file))[0]

            # Excel sheet names have a strict 31-character limit
            sheet_name = sheet_name[:31]

            print(f"Processing: {file} -> Sheet: {sheet_name}")

            # 3. Read the CSV and write it to its own Excel sheet
            df = pd.read_csv(file)
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"\nSuccess! All files combined into '{output_excel}'")